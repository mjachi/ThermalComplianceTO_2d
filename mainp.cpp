#include "mfem.hpp"
#include <iostream>
#include <fstream>
#include "main.hxx"

using namespace std;
using namespace mfem;


/**
 * @brief Creates a mesh on a unit square with the middle third cut out.
 *      
 * @param ref_levels Number of refinement levels
 */
Mesh UnitSquare_Geo(int ref_levels) {

  // Create 1x1 square of 9 quadrilaterals
  // auto mesh = Mesh::MakeCartesian2D(3, 3, mfem::Element::Type::QUADRILATERAL, true, 1.0, 1.0);
  auto mesh = Mesh::MakeCartesian2D(3, 3, mfem::Element::Type::QUADRILATERAL);

  // Set boundary conditions
  int i;
  for (i = 0; i < mesh.GetNBE(); i++) {
    Element *be = mesh.GetBdrElement(i);
    Array<int> vertices;
    be->GetVertices(vertices);

    double *coords1 = mesh.GetVertex(vertices[0]);
    double *coords2 = mesh.GetVertex(vertices[1]);

    Vector center(2);
    center(0) = 0.5 * (coords1[0] + coords2[0]);
    center(1) = 0.5 * (coords1[1] + coords2[1]);

    if ((center(0) == 0.5) && (center(1) == 1)) {
      be->SetAttribute(1);
    } else {
      be->SetAttribute(2);
    }
  }
  mesh.SetAttributes();

  // Refine mesh

  for (i = 0; i < ref_levels; i++) {
    mesh.UniformRefinement();
  }

  return mesh;
}


/**
 * @brief Nonlinear projection of ψ onto the subspace
 *        ∫_Ω sigmoid(ψ) dx = θ vol(Ω) as follows.
 *
 *        1. Compute the root of the R → R function
 *            f(c) = ∫_Ω sigmoid(ψ + c) dx - θ vol(Ω)
 *        2. Set ψ ← ψ + c.
 *
 * @param psi a GridFunction to be updated
 * @param target_volume θ vol(Ω)
 * @param tol Newton iteration tolerance
 * @param max_its Newton maximum iteration number
 * @return double Final volume, ∫_Ω sigmoid(ψ)
 */
double projit(ParGridFunction &psi, double target_volume, double tol=1e-12,
              int max_its=10)
{
    MappedGridFunctionCoefficient sigmoid_psi(&psi, sigmoid);
    MappedGridFunctionCoefficient der_sigmoid_psi(&psi, der_sigmoid);

    ParLinearForm int_sigmoid_psi(psi.ParFESpace());
    int_sigmoid_psi.AddDomainIntegrator(new DomainLFIntegrator(sigmoid_psi));
    ParLinearForm int_der_sigmoid_psi(psi.ParFESpace());
    int_der_sigmoid_psi.AddDomainIntegrator(new DomainLFIntegrator(
                                                der_sigmoid_psi));
    bool done = false;
    for (int k=0; k<max_its; k++) // Newton iteration
    {
        int_sigmoid_psi.Assemble(); // Recompute f(c) with updated ψ
        double f = int_sigmoid_psi.Sum();
        MPI_Allreduce(MPI_IN_PLACE, &f, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        f -= target_volume;

        int_der_sigmoid_psi.Assemble(); // Recompute df(c) with updated ψ
        double df = int_der_sigmoid_psi.Sum();
        MPI_Allreduce(MPI_IN_PLACE, &df, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        const double dc = -f/df;
        psi += dc;
        if (abs(dc) < tol) { done = true; break; }
    }
    if (!done)
    {
        mfem_warning("Projection reached maximum iteration without converging. Result may not be accurate.");
    }
    int_sigmoid_psi.Assemble();
    double material_volume = int_sigmoid_psi.Sum();
    MPI_Allreduce(MPI_IN_PLACE, &material_volume, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);
    return material_volume;
}


/**
 * @brief Enforce boundedness, -max_val ≤ psi ≤ max_val
 *
 * @param psi a GridFunction to be bounded (in place)
 * @param max_val upper and lower bound
 */
inline void clip(ParGridFunction &psi, const double max_val)
{
   for (auto &val : psi) { val = min(max_val, max(-max_val, val)); }
}


int main(int argc, char *argv[])
{
    // 0. Initialize MPI and HYPRE.
    Mpi::Init();
    int num_procs = Mpi::WorldSize();
    int myid = Mpi::WorldRank();
    Hypre::Init();

    // 1. Parse command-line options.
    int ref_levels = 4;
    int order = 2;
    bool visualization = true;
    double alpha = 1.0;
    double epsilon = 0.01;
    double mass_fraction = 0.5;
    int max_it = 1e3;
    double tol = 1e-4;
    double rho_min = 1e-6;
    double lambda = 1.0;
    double mu = 1.0;

    OptionsParser args(argc, argv);
    args.AddOption(&ref_levels, "-r", "--refine",
                    "Number of times to refine the mesh uniformly.");
    args.AddOption(&order, "-o", "--order",
                    "Order (degree) of the finite elements.");
    args.AddOption(&alpha, "-alpha", "--alpha-step-length",
                    "Step length for gradient descent.");
    args.AddOption(&epsilon, "-epsilon", "--epsilon-thickness",
                    "epsilon phase field thickness");
    args.AddOption(&max_it, "-mi", "--max-it",
                    "Maximum number of gradient descent iterations.");
    args.AddOption(&tol, "-tol", "--tol",
                    "Exit tolerance for ρ ");
    args.AddOption(&mass_fraction, "-mf", "--mass-fraction",
                    "Mass fraction for diffusion coefficient.");
    args.AddOption(&lambda, "-lambda", "--lambda",
                    "Lame constant λ");
    args.AddOption(&mu, "-mu", "--mu",
                    "Lame constant μ");
    args.AddOption(&rho_min, "-rmin", "--psi-min",
                    "Minimum of density coefficient.");
    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                    "--no-visualization",
                    "Enable or disable GLVis visualization.");

    args.Parse();
    if (!args.Good())
    {
        if (myid == 0)
        {
            args.PrintUsage(cout);
        }
        MPI_Finalize();
        return 1;
    }
    if (myid == 0)
    {
        mfem::out << num_procs << " number of process created.\n";
        args.PrintOptions(cout);
    }

    auto mesh = UnitSquare_Geo(ref_levels);
    int dim = mesh.Dimension();

    ParMesh pmesh(MPI_COMM_WORLD, mesh);
    mesh.Clear();

    // 4. Define the necessary finite element spaces on the mesh.
    H1_FECollection state_fec(order, dim); // space for u
    H1_FECollection filter_fec(order, dim); // space for ρ̃
    L2_FECollection control_fec(order-1, dim,
                                BasisType::GaussLobatto); // space for ψ
    ParFiniteElementSpace state_fes(&pmesh, &state_fec);
    ParFiniteElementSpace filter_fes(&pmesh, &filter_fec);
    ParFiniteElementSpace control_fes(&pmesh, &control_fec);

    HYPRE_BigInt state_size = state_fes.GlobalTrueVSize();
    HYPRE_BigInt control_size = control_fes.GlobalTrueVSize();
    HYPRE_BigInt filter_size = filter_fes.GlobalTrueVSize();
    if (myid==0)
    {
        cout << "Number of state unknowns: " << state_size << endl;
        cout << "Number of filter unknowns: " << filter_size << endl;
        cout << "Number of control unknowns: " << control_size << endl;
    }

    // 5. Set the initial guess for ρ.
    ParGridFunction u(&state_fes);
    ParGridFunction psi(&control_fes);
    ParGridFunction psi_old(&control_fes);
    ParGridFunction rho_filter(&filter_fes);
    u = 0.0;
    rho_filter = mass_fraction;
    psi = inv_sigmoid(mass_fraction);
    psi_old = inv_sigmoid(mass_fraction);

    const double sigmoid_bound = -inv_sigmoid(rho_min);

    // ρ = sigmoid(ψ)
    MappedGridFunctionCoefficient rho(&psi, sigmoid);
    // Interpolation of ρ = sigmoid(ψ) in control fes
    ParGridFunction rho_gf(&control_fes);
    // ρ - ρ_old = sigmoid(ψ) - sigmoid(ψ_old)
    DiffMappedGridFunctionCoefficient succ_diff_rho(&psi, &psi_old, sigmoid);

    // 6. Set-up the physics solver.
    int maxat = pmesh.bdr_attributes.Max();
    Array<int> ess_bdr(maxat);
    Array<int> nbc_bdr(maxat);
    ess_bdr = 0; nbc_bdr = 0;
    ess_bdr[0] = 1; nbc_bdr[1] = 1;
    ConstantCoefficient one_cp(1.0);
    ConstantCoefficient one_cp_2(1.0);
    

    DiffusionSolver *DiffSolver = new DiffusionSolver();
    DiffSolver->SetMesh(&pmesh);
    DiffSolver->SetOrder(state_fec.GetOrder());
    DiffSolver->SetRHSCoefficient(&one_cp);
    DiffSolver->SetupFEM();
    DiffSolver->SetEssentialBoundary(ess_bdr);

    // 7. Set-up the filter solver.
    ConstantCoefficient one(1.0);
    ConstantCoefficient eps2_cf(epsilon*epsilon);
    DiffusionSolver * FilterSolver = new DiffusionSolver();
    FilterSolver->SetMesh(&pmesh);
    FilterSolver->SetOrder(filter_fec.GetOrder());
    FilterSolver->SetDiffusionCoefficient(&eps2_cf);
    FilterSolver->SetMassCoefficient(&one);
    Array<int> ess_bdr_filter;
    if (pmesh.bdr_attributes.Size())
    {
        ess_bdr_filter.SetSize(pmesh.bdr_attributes.Max());
        ess_bdr_filter = 0;
    }
    FilterSolver->SetEssentialBoundary(ess_bdr_filter);
    FilterSolver->SetupFEM();

    ParBilinearForm mass(&control_fes);
    mass.AddDomainIntegrator(new InverseIntegrator(new MassIntegrator(one)));
    mass.Assemble();
    HypreParMatrix M;
    Array<int> empty;
    mass.FormSystemMatrix(empty,M);

    // 8. Define the Lagrange multiplier and gradient functions
    ParGridFunction grad(&control_fes);
    ParGridFunction w_filter(&filter_fes);

    // 9. Define some tools for later
    ConstantCoefficient zero(0.0);
    ParGridFunction onegf(&control_fes);
    onegf = 1.0;
    ParGridFunction zerogf(&control_fes);
    zerogf = 0.0;
    ParLinearForm vol_form(&control_fes);
    vol_form.AddDomainIntegrator(new DomainLFIntegrator(one));
    vol_form.Assemble();
    double domain_volume = vol_form(onegf);
    const double target_volume = domain_volume * mass_fraction;

    // 10. Connect to GLVis. Prepare for VisIt output.
    char vishost[] = "localhost";
    int  visport   = 19916;
    socketstream sout_u,sout_r,sout_rho;
    if (visualization)
    {
        sout_u.open(vishost, visport);
        sout_rho.open(vishost, visport);
        sout_r.open(vishost, visport);
        sout_u.precision(8);
        sout_rho.precision(8);
        sout_r.precision(8);
    }

    rho_gf.ProjectCoefficient(rho);
    mfem::ParaViewDataCollection paraview_dc("Elastic_compliance", &pmesh);
    paraview_dc.SetPrefixPath("ParaView");
    paraview_dc.SetLevelsOfDetail(order);
    paraview_dc.SetCycle(0);
    paraview_dc.SetDataFormat(VTKFormat::BINARY);
    paraview_dc.SetHighOrderOutput(true);
    paraview_dc.SetTime(0.0);
    paraview_dc.RegisterField("displacement",&u);
    paraview_dc.RegisterField("density",&rho_gf);
    paraview_dc.RegisterField("filtered_density",&rho_filter);

    // 11. Iterate
    int step = 0;
    //double c0 = 0.0;
    for (int k = 1; k < max_it; k++)
    {
        if (k > 1) { alpha *= ((double) k) / ((double) k-1); }
        step++;

        if (myid == 0)
        {
            cout << "\nStep = " << k << endl;
        }

        // Step 1 - Filter solve
        // Solve (ϵ^2 ∇ ρ̃, ∇ v ) + (ρ̃,v) = (ρ,v)
        FilterSolver->SetRHSCoefficient(&rho);
        FilterSolver->Solve();
        rho_filter = *FilterSolver->GetFEMSolution();

        // Step 2 - State solve
        // Solve (λ(ρ̃) ∇⋅u, ∇⋅v) + (2 μ(ρ̃) ε(u), ε(v)) = (f,v)
        SIMPInterpolationCoefficient SIMP_cf(&rho_filter,rho_min, 1.0);
        DiffSolver->SetDiffusionCoefficient(&SIMP_cf);       
        DiffSolver->SetupFEM();
        DiffSolver->Solve();
        u = *DiffSolver->GetFEMSolution();

        // Step 3 - Adjoint filter solve
        // Solve (ϵ² ∇ w̃, ∇ v) + (w̃ ,v) = (-r'(ρ̃) ( λ(ρ̃) |∇⋅u|² + 2 μ(ρ̃) |ε(u)|²),v)
        DirichletEnergyCoefficient d_cf(&u, &rho_filter, rho_min);
        FilterSolver->SetRHSCoefficient(&d_cf);
        FilterSolver->Solve();
        w_filter = *FilterSolver->GetFEMSolution();

        // Step 4 - Compute gradient
        // Solve G = M⁻¹w̃
        GridFunctionCoefficient w_cf(&w_filter);
        ParLinearForm w_rhs(&control_fes);
        w_rhs.AddDomainIntegrator(new DomainLFIntegrator(w_cf));
        w_rhs.Assemble();
        M.Mult(w_rhs,grad);

        // Step 5 - Update design variable ψ ← clip(projit(ψ - αG))
        psi.Add(-alpha, grad);
        const double material_volume = projit(psi, target_volume);
        clip(psi, sigmoid_bound);

        // Compute ||ρ - ρ_old|| in control fes.
        double norm_reduced_gradient = zerogf.ComputeL2Error(succ_diff_rho)/alpha;
        psi_old = psi;

        double compliance = (*(DiffSolver->GetLinearForm()))(u);
        MPI_Allreduce(MPI_IN_PLACE,&compliance,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        if (myid == 0)
        {
            mfem::out << "norm of reduced gradient = " << norm_reduced_gradient << endl;
            mfem::out << "compliance = " << compliance << endl;
            mfem::out << "mass_fraction = " << material_volume / domain_volume << endl;
        }

        if (visualization)
        {
            sout_u << "parallel " << num_procs << " " << myid << "\n";
            sout_u << "solution\n" << pmesh << u
                << "window_title 'Displacement u'" << flush;

            rho_gf.ProjectCoefficient(rho);
            sout_rho << "parallel " << num_procs << " " << myid << "\n";
            sout_rho << "solution\n" << pmesh << rho_gf
                    << "window_title 'Control variable ρ'" << flush;

            ParGridFunction r_gf(&filter_fes);
            r_gf.ProjectCoefficient(SIMP_cf);
            sout_r << "parallel " << num_procs << " " << myid << "\n";
            sout_r << "solution\n" << pmesh << r_gf
                << "window_title 'Design density r(ρ̃)'" << flush;

            paraview_dc.SetCycle(k);
            paraview_dc.SetTime((double)k);
            paraview_dc.Save();
        }

        if (norm_reduced_gradient < tol)
        {
            break;
        }
    }

    delete DiffSolver;
    delete FilterSolver;

    return 0;
}