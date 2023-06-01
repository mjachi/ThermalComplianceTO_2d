#include "main.hpp"
#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;
using namespace mfem;
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
double projit(GridFunction &psi, double target_volume, double tol = 1e-12,
              int max_its = 10) {
  MappedGridFunctionCoefficient sigmoid_psi(&psi, sigmoid);
  MappedGridFunctionCoefficient der_sigmoid_psi(&psi, der_sigmoid);

  LinearForm int_sigmoid_psi(psi.FESpace());
  int_sigmoid_psi.AddDomainIntegrator(new DomainLFIntegrator(sigmoid_psi));
  LinearForm int_der_sigmoid_psi(psi.FESpace());
  int_der_sigmoid_psi.AddDomainIntegrator(
      new DomainLFIntegrator(der_sigmoid_psi));
  bool done = false;
  for (int k = 0; k < max_its; k++) // Newton iteration
  {
    int_sigmoid_psi.Assemble(); // Recompute f(c) with updated ψ
    const double f = int_sigmoid_psi.Sum() - target_volume;

    int_der_sigmoid_psi.Assemble(); // Recompute df(c) with updated ψ
    const double df = int_der_sigmoid_psi.Sum();

    const double dc = -f / df;
    psi += dc;
    if (abs(dc) < tol) {
      done = true;
      break;
    }
  }
  if (!done) {
    mfem_warning("Projection reached maximum iteration without converging. "
                 "Result may not be accurate.");
  }
  int_sigmoid_psi.Assemble();
  return int_sigmoid_psi.Sum();
}

/**
 * @brief Enforce boundedness, -max_val ≤ psi ≤ max_val
 *
 * @param psi a GridFunction to be bounded (in place)
 * @param max_val upper and lower bound
 */
inline void clip(GridFunction &psi, const double max_val) {
  for (auto &val : psi) {
    val = min(max_val, max(-max_val, val));
  }
}

/**
 * @brief Driving function.
*/
int main(int argc, char *argv[]) {

  // 1. Parse command-line options.
  
  // PARAMETERS
  int ref_levels = 4;          // # of mesh refinement levels
  int order = 2;               // Order of FE
  double alpha = 0.01;          // Step size
  double epsilon = 0.01;       // Density filter radius
  double mass_fraction = 0.5;  // Mass fraction in Omega
  int max_it = 1e3;            // Maximum number of iterations for Newton solver in projection
  double tol = 1e-4;           // Convergence tolerance for change in density field
  double rho_min = 1e-3;       // Minimum of density coefficient (for modified SIMP law)

  // OPTIONS
  bool visualization = true;   // Outputs to GLVis and
  int problem_type = 1;        // Integer for problem type (affects domain and BCs)
                               //  1 - Classical unit square problem (one sink, even heating)
                               //  2 - Unit square with two sinks and even heating
                               //  3 - L-shape with even heating
                               //  4 - L-shape with inhomogeneous Neumann BCs on right-most edge
                               //  5 - "Chainlink" heat-sink problem


  // Command-line options
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
  args.AddOption(&tol, "-tol", "--tol", "Exit tolerance for ρ ");
  args.AddOption(&mass_fraction, "-mf", "--mass-fraction",
                 "Mass fraction for diffusion coefficient.");
  args.AddOption(&rho_min, "-rmin", "--psi-min",
                 "Minimum of density coefficient.");
  args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                 "--no-visualization",
                 "Enable or disable GLVis visualization.");
  args.AddOption(&problem_type, "-t", "--type", "Name problem type");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(mfem::out);
    return 1;
  }
  if (problem_type > 5 || problem_type < 1) {
    mfem::out << "Problem type integer (\"-t\" or \"--type\" argument) must be 1-5" << std::endl;
    return 1;
  }


  args.PrintOptions(mfem::out);

  Mesh mesh;
  switch (problem_type) {
    case 1:
      mesh = UnitSquare_Geo(ref_levels);
      break;
    case 2:
      mesh = UnitSquare_Geo_2Sinks(ref_levels);
      break;
    case 3:
      mesh = LShape_Homogeneous(ref_levels);
      break;
    case 4:
      mesh = LShape_Inhomo(ref_levels);
      break;
    case 5:
      mesh = HeatSink(ref_levels);
      break;
  }
  int dim = mesh.Dimension();

  // 4. Define the necessary finite element spaces on the mesh.
  H1_FECollection state_fec(order, dim);  // space for u
  H1_FECollection filter_fec(order, dim); // space for ρ̃
  L2_FECollection control_fec(order - 1, dim,
                              BasisType::GaussLobatto); // space for ψ
  FiniteElementSpace state_fes(&mesh, &state_fec);
  FiniteElementSpace filter_fes(&mesh, &filter_fec);
  FiniteElementSpace control_fes(&mesh, &control_fec);

  int state_size = state_fes.GetTrueVSize();
  int control_size = control_fes.GetTrueVSize();
  int filter_size = filter_fes.GetTrueVSize();
  mfem::out << "Number of state unknowns: " << state_size << std::endl;
  mfem::out << "Number of filter unknowns: " << filter_size << std::endl;
  mfem::out << "Number of control unknowns: " << control_size << std::endl;

  // 5. Set the initial guess for ρ.
  GridFunction u(&state_fes);
  GridFunction psi(&control_fes);
  GridFunction psi_old(&control_fes);
  GridFunction rho_filter(&filter_fes);
  u = 0.0;
  rho_filter = mass_fraction;
  psi = inv_sigmoid(mass_fraction);
  psi_old = inv_sigmoid(mass_fraction);

  const double sigmoid_bound = -inv_sigmoid(rho_min);

  // ρ = sigmoid(ψ)
  MappedGridFunctionCoefficient rho(&psi, sigmoid);
  // Interpolation of ρ = sigmoid(ψ) in control fes
  GridFunction rho_gf(&control_fes);
  // ρ - ρ_old = sigmoid(ψ) - sigmoid(ψ_old)
  DiffMappedGridFunctionCoefficient succ_diff_rho(&psi, &psi_old, sigmoid);

  // 6. Set-up the physics solver.
  int maxat = mesh.bdr_attributes.Max();
  Array<int> ess_bdr(maxat);
  Array<int> nbc_bdr(maxat);
  ess_bdr = 0; nbc_bdr = 0;

  std::vector<int> homog_nbc { 1, 2, 3 };
  auto nbc_beg = homog_nbc.begin(); auto nbc_end = homog_nbc.end();
  if (std::count(nbc_beg, nbc_end, problem_type)) {
    ess_bdr[0] = 1; nbc_bdr[1] = 1;
  } else {
    ess_bdr[0] = 1; nbc_bdr[2] = 1;
  }

  ConstantCoefficient one_cp(1.0);
  ConstantCoefficient zero_cp(0.0);
  ConstantCoefficient small_cp(0.125);
  ConstantCoefficient one_cp2(1.0);

  DiffusionSolver *DiffSolver = new DiffusionSolver();
  DiffSolver->SetMesh(&mesh);
  DiffSolver->SetOrder(state_fec.GetOrder());
  DiffSolver->SetEssentialBoundary(ess_bdr);
  if (std::count(nbc_beg, nbc_end, problem_type)) { // Constant heating inside
    DiffSolver->SetRHSCoefficient(&one_cp);
  } else {
    // Zero heating on interior
    DiffSolver->SetRHSCoefficient(&zero_cp);
    
    // Set non-homogeneous Neumann BC accordingly
    DiffSolver->SetNeumannBoundary(nbc_bdr);
    DiffSolver->SetNeumannData(&one_cp);
  }
  DiffSolver->SetupFEM();

  // 7. Set-up the filter solver.
  ConstantCoefficient one(1.0);
  ConstantCoefficient eps2_cf(epsilon * epsilon);
  DiffusionSolver *FilterSolver = new DiffusionSolver();
  FilterSolver->SetMesh(&mesh);
  FilterSolver->SetOrder(filter_fec.GetOrder());
  FilterSolver->SetDiffusionCoefficient(&eps2_cf);
  FilterSolver->SetMassCoefficient(&one);
  Array<int> ess_bdr_filter;
  if (mesh.bdr_attributes.Size()) {
    ess_bdr_filter.SetSize(mesh.bdr_attributes.Max());
    ess_bdr_filter = 0;
  }
  FilterSolver->SetEssentialBoundary(ess_bdr_filter);
  FilterSolver->SetupFEM();

  BilinearForm mass(&control_fes);
  mass.AddDomainIntegrator(new InverseIntegrator(new MassIntegrator(one)));
  mass.Assemble();
  SparseMatrix M;
  Array<int> empty;
  mass.FormSystemMatrix(empty, M);

  // 8. Define the Lagrange multiplier and gradient functions
  GridFunction grad(&control_fes);
  GridFunction w_filter(&filter_fes);

  // 9. Define some tools for later
  ConstantCoefficient zero(0.0);
  GridFunction onegf(&control_fes);
  onegf = 1.0;
  GridFunction zerogf(&control_fes);
  zerogf = 0.0;
  LinearForm vol_form(&control_fes);
  vol_form.AddDomainIntegrator(new DomainLFIntegrator(one));
  vol_form.Assemble();
  double domain_volume = vol_form(onegf);
  const double target_volume = domain_volume * mass_fraction;

  // 10. Connect to GLVis. Prepare for VisIt output.
  char vishost[] = "localhost";
  int visport = 19916;
  socketstream sout_u, sout_r, sout_rho;
  if (visualization) {
    sout_u.open(vishost, visport);
    sout_rho.open(vishost, visport);
    sout_r.open(vishost, visport);
    sout_u.precision(8);
    sout_rho.precision(8);
    sout_r.precision(8);
  }

  rho_gf.ProjectCoefficient(rho);
  mfem::ParaViewDataCollection paraview_dc("Thermal_compliance", &mesh);
  paraview_dc.SetPrefixPath("ParaView");
  paraview_dc.SetLevelsOfDetail(order);
  paraview_dc.SetCycle(0);
  paraview_dc.SetDataFormat(VTKFormat::BINARY);
  paraview_dc.SetHighOrderOutput(true);
  paraview_dc.SetTime(0.0);
  paraview_dc.RegisterField("displacement", &u);
  paraview_dc.RegisterField("density", &rho_gf);
  paraview_dc.RegisterField("filtered_density", &rho_filter);

  // 1 - Parallelize 
  // 2 - interesting domains
  //   - "Chainlink with hole in center"
  //   - "L-shaped", etc ...
  //   - non-uniform heating, ...


  // 11. Iterate
  int step = 0;
  // double c0 = 0.0;
  for (int k = 1; k < max_it; k++) {
    if (k > 1) {
      alpha *= ((double)k) / ((double)k - 1);
    }
    step++;

    mfem::out << "\nStep = " << k << std::endl;

    // Step 1 - Filter solve
    // OLD: Solve (ϵ^2 ∇ ρ̃, ∇ v ) + (ρ̃,v) = (ρ,v)
    // NEW: Same.
    FilterSolver->SetRHSCoefficient(&rho);
    FilterSolver->Solve();
    rho_filter = *FilterSolver->GetFEMSolution();

    // Step 2 - State solve
    // Solve (r(ρ̃) ∇ u, ∇ v) = (f,v)
    SIMPInterpolationCoefficient SIMP_cf(&rho_filter, rho_min, 1.0);
    DiffSolver->SetDiffusionCoefficient(&SIMP_cf);
    DiffSolver->SetupFEM();
    DiffSolver->Solve();
    u = *DiffSolver->GetFEMSolution();

    // Step 3 - Adjoint filter solve
    // Solve (ϵ² ∇ w̃, ∇ v) + (w̃, v) = (-r'(ρ̃) |∇u|²,v)
    DirichletEnergyCoefficient d_cf(&u, &rho_filter, rho_min);
    FilterSolver->SetRHSCoefficient(&d_cf);
    FilterSolver->Solve();
    w_filter = *FilterSolver->GetFEMSolution();

    // Step 4 - Compute gradient
    // Solve G = M⁻¹w̃
    GridFunctionCoefficient w_cf(&w_filter);
    LinearForm w_rhs(&control_fes);
    w_rhs.AddDomainIntegrator(new DomainLFIntegrator(w_cf));
    w_rhs.Assemble();
    M.Mult(w_rhs, grad);

    // Step 5 - Update design variable ψ ← clip(projit(ψ - αG))
    psi.Add(-alpha, grad);
    const double material_volume = projit(psi, target_volume);
    clip(psi, sigmoid_bound);

    // Compute ||ρ - ρ_old|| in control fes.
    double norm_reduced_gradient = zerogf.ComputeL2Error(succ_diff_rho) / alpha;
    psi_old = psi;

    double compliance = (*(DiffSolver->GetLinearForm()))(u);
    mfem::out << "norm of reduced gradient = " << norm_reduced_gradient
              << std::endl;
    mfem::out << "compliance = " << compliance << std::endl;
    mfem::out << "mass_fraction = " << material_volume / domain_volume
              << std::endl;
    mfem::out << "u L2 norm: " << u.Norml2() << std::endl;
    

    if (visualization) {
      sout_u << "solution\n"
             << mesh << u << "window_title 'Displacement u'" << flush;

      rho_gf.ProjectCoefficient(rho);
      sout_rho << "solution\n"
               << mesh << rho_gf << "window_title 'Control variable ρ'"
               << flush;

      GridFunction r_gf(&filter_fes);
      r_gf.ProjectCoefficient(SIMP_cf);
      sout_r << "solution\n"
             << mesh << r_gf << "window_title 'Design density r(ρ̃)'" << flush;

      paraview_dc.SetCycle(k);
      paraview_dc.SetTime((double)k);
      paraview_dc.Save();
    }

    if (norm_reduced_gradient < tol) {
      break;
    }
  }

  delete DiffSolver;
  delete FilterSolver;

  return 0;
}
