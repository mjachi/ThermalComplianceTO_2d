//                  MFEM Example 35 - Serial/Parallel Shared Code
//
//

#include "mfem.hpp"
#include <fstream>
#include <functional>
#include <iostream>



namespace mfem {
///////////////////////////////////////////////////////////////////////////////////////////////////////
// Meshing
///////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Creates a mesh on a unit square with the middle third cut out.
 *      
 * @param ref_levels Number of refinement levels
 */
Mesh UnitSquare_Geo(int ref_levels) {

  // Create 1x1 square of 9 quadrilaterals
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
 * @brief Creates a mesh on a unit square with the middle third cut out.
 *      
 * @param ref_levels Number of refinement levels
 */
Mesh UnitSquare_Geo_2Sinks(int ref_levels) {

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

    if (((center(0) == 0.5) && (center(1) == 1)) || (center(0) == 0.5 && center(1) == 0)) {
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
 * @brief Creates a mesh on an L-shaped domain with even heating.
 *      
 * @param ref_levels Number of refinement levels
 */
Mesh LShape_Homogeneous(int ref_levels) {

  mfem::Mesh mesh("../data/LShaped1.mesh", 1, 1);
  // Refine mesh

  int i;
  for (i = 0; i < ref_levels; i++) {
    mesh.UniformRefinement();
  }

  return mesh;
}

/**
 * @brief Creates a mesh on an L-shaped domain, which is functionally
 *  a unit square with the top quarter removed. Along the top-most
 *  edge, heat is allowed to escape; along the right-most edge
 *  inhomogeneous Neumann BCs are set to act as a source. There
 *  is no internal heating.
 *      
 * @param ref_levels Number of refinement levels
 */
Mesh LShape_Inhomo(int ref_levels) {
  Mesh mesh("../data/LShaped2.mesh", 1, 1);

  int i;
  for (i = 0; i < ref_levels; i++) {
    mesh.UniformRefinement();
  }

  return mesh;
}

/**
 * @brief Creates a mesh on a ring-like domain with a square
 *  cutout in the middle. Along this cutout, inhomogeneous
 *  Neumann BCs are set to act as a source. At the edge of the two
 *  "arms", heat is allowed to escape. There is no internal heating.
 *      
 * @param ref_levels Number of refinement levels
 */
Mesh HeatSink(int ref_levels) {
  Mesh mesh("../data/HeatSink.mesh", 1, 1);

  int i;
  for (i = 0; i < ref_levels; i++) {
    mesh.UniformRefinement();
  }

  return mesh;
}

// Inverse sigmoid function
double inv_sigmoid(double x) {
  double tol = 1e-12;
  x = std::min(std::max(tol, x), 1.0 - tol);
  return std::log(x / (1.0 - x));
}

// Sigmoid function
double sigmoid(double x) {
  if (x >= 0) {
    return 1.0 / (1.0 + std::exp(-x));
  } else {
    return std::exp(x) / (1.0 + std::exp(x));
  }
}

// Derivative of sigmoid function
double der_sigmoid(double x) {
  double tmp = sigmoid(-x);
  return tmp - std::pow(tmp, 2);
}

/**
 * @brief Returns f(u(x)) where u is a scalar GridFunction and f:R → R
 *
 */
class MappedGridFunctionCoefficient : public GridFunctionCoefficient {
protected:
  std::function<double(const double)> fun; // f:R → R
public:
  MappedGridFunctionCoefficient()
      : GridFunctionCoefficient(), fun([](double x) { return x; }) {}
  MappedGridFunctionCoefficient(const GridFunction *gf,
                                std::function<double(const double)> fun_,
                                int comp = 1)
      : GridFunctionCoefficient(gf, comp), fun(fun_) {}

  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip) {
    return fun(GridFunctionCoefficient::Eval(T, ip));
  }
  void SetFunction(std::function<double(const double)> fun_) { fun = fun_; }
};

/**
 * @brief Returns f(u(x)) - f(v(x)) where u, v are scalar GridFunctions and f:R
 * → R
 *
 */
class DiffMappedGridFunctionCoefficient : public GridFunctionCoefficient {
protected:
  std::function<double(const double)> fun; // f:R → R
  const GridFunction *OtherGridF;
  GridFunctionCoefficient OtherGridF_cf;

public:
  DiffMappedGridFunctionCoefficient()
      : GridFunctionCoefficient(), OtherGridF(nullptr), OtherGridF_cf(),
        fun([](double x) { return x; }) {}
  DiffMappedGridFunctionCoefficient(const GridFunction *gf,
                                    const GridFunction *other_gf,
                                    std::function<double(const double)> fun_,
                                    int comp = 1)
      : GridFunctionCoefficient(gf, comp), fun(fun_), OtherGridF(other_gf),
        OtherGridF_cf(OtherGridF) {}

  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip) {
    const double value1 = fun(GridFunctionCoefficient::Eval(T, ip));
    const double value2 = fun(OtherGridF_cf.Eval(T, ip));
    return value1 - value2;
  }
  void SetFunction(std::function<double(const double)> fun_) { fun = fun_; }
};

//  Solid isotropic material penalization (SIMP) coefficient
class SIMPInterpolationCoefficient : public Coefficient {
protected:
  GridFunction *rho_filter; // grid function
  double min_val;
  double max_val;
  double exponent;

public:
  SIMPInterpolationCoefficient(GridFunction *rho_filter_,
                               double min_val_ = 1e-6, double max_val_ = 1.0,
                               double exponent_ = 3)
      : rho_filter(rho_filter_), min_val(min_val_), max_val(max_val_),
        exponent(exponent_) {}

  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip) {
    double val = rho_filter->GetValue(T, ip);
    double coeff = min_val + pow(val, exponent) * (max_val - min_val);
    return coeff;
  }
};

// Dirichlet energy density coefficient
class DirichletEnergyCoefficient : public Coefficient {
protected:
  GridFunction *u = nullptr;          // displacement
  GridFunction *rho_filter = nullptr; // filter density
  Vector grad;                   // auxiliary matrix, used in Eval
  double exponent;
  double rho_min;

public:
  DirichletEnergyCoefficient(GridFunction *u_, GridFunction *rho_filter_,
                                double rho_min_ = 1e-6, double exponent_ = 3.0)
      : u(u_), rho_filter(rho_filter_), exponent(exponent_), rho_min(rho_min_) {
    MFEM_ASSERT(rho_min_ >= 0.0, "rho_min must be >= 0");
    MFEM_ASSERT(rho_min_ < 1.0, "rho_min must be > 1");
    MFEM_ASSERT(u, "displacement field is not set");
    MFEM_ASSERT(rho_filter, "density field is not set");
  }

  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip) {
    u->GetGradient(T, grad);
    double density = grad.Norml2();
    density = abs(density)*abs(density);
    double val = rho_filter->GetValue(T, ip);

    return (-exponent * pow(val, exponent - 1.0) * (1 - rho_min) * density);
  }
};

//  Class for solving Poisson's equation:
//
//       - ∇ ⋅(κ ∇ u) = f  in Ω
//
class DiffusionSolver {
private:
  Mesh *mesh = nullptr;
  int order = 1;
  // diffusion coefficient
  Coefficient *diffcf = nullptr;
  // mass coefficient
  Coefficient *masscf = nullptr;
  Coefficient *rhscf = nullptr;
  Coefficient *essbdr_cf = nullptr;
  Coefficient *neumann_cf = nullptr;
  VectorCoefficient *gradient_cf = nullptr;

  // FEM solver
  int dim;
  FiniteElementCollection *fec = nullptr;
  FiniteElementSpace *fes = nullptr;
  Array<int> ess_bdr;
  Array<int> neumann_bdr;
  GridFunction *u = nullptr;
  LinearForm *b = nullptr;
  bool parallel;
#ifdef MFEM_USE_MPI
  ParMesh *pmesh = nullptr;
  ParFiniteElementSpace *pfes = nullptr;
#endif

public:
  DiffusionSolver() {}
  DiffusionSolver(Mesh *mesh_, int order_, Coefficient *diffcf_,
                  Coefficient *cf_);

  void SetMesh(Mesh *mesh_) {
    mesh = mesh_;
    parallel = false;
#ifdef MFEM_USE_MPI
    pmesh = dynamic_cast<ParMesh *>(mesh);
    if (pmesh) {parallel = true;}
#endif


  }
  void SetOrder(int order_) { order = order_; }
  void SetDiffusionCoefficient(Coefficient *diffcf_) { diffcf = diffcf_; }
  void SetMassCoefficient(Coefficient *masscf_) { masscf = masscf_; }
  void SetRHSCoefficient(Coefficient *rhscf_) { rhscf = rhscf_; }
  void SetEssentialBoundary(const Array<int> &ess_bdr_) { ess_bdr = ess_bdr_; };
  void SetNeumannBoundary(const Array<int> &neumann_bdr_) {
    neumann_bdr = neumann_bdr_;
  };
  void SetNeumannData(Coefficient *neumann_cf_) { neumann_cf = neumann_cf_; }
  void SetEssBdrData(Coefficient *essbdr_cf_) { essbdr_cf = essbdr_cf_; }
  void SetGradientData(VectorCoefficient *gradient_cf_) {
    gradient_cf = gradient_cf_;
  }

  void ResetFEM();
  void SetupFEM();

  void Solve();
  GridFunction *GetFEMSolution();
  LinearForm *GetLinearForm() { return b; }
#ifdef MFEM_USE_MPI
   ParGridFunction * GetParFEMSolution();
   ParLinearForm * GetParLinearForm()
   {
      if (parallel)
      {
         return dynamic_cast<ParLinearForm *>(b);
      }
      else
      {
         MFEM_ABORT("Wrong code path. Call GetLinearForm");
         return nullptr;
      }
   }
#endif

  ~DiffusionSolver();
};

// -----------------------------------------------------------------------
// --------------------      Poisson solver     --------------------------
// -----------------------------------------------------------------------

DiffusionSolver::DiffusionSolver(Mesh *mesh_, int order_, Coefficient *diffcf_,
                                 Coefficient *rhscf_)
    : mesh(mesh_), order(order_), diffcf(diffcf_), rhscf(rhscf_) {
#ifdef MFEM_USE_MPI
   pmesh = dynamic_cast<ParMesh *>(mesh);
   if (pmesh) { parallel = true; }
#endif

  SetupFEM();
}

void DiffusionSolver::SetupFEM() {
  dim = mesh->Dimension();
  fec = new H1_FECollection(order, dim);

#ifdef MFEM_USE_MPI
   if (parallel)
   {
      pfes = new ParFiniteElementSpace(pmesh, fec);
      u = new ParGridFunction(pfes);
      b = new ParLinearForm(pfes);
   }
   else
   {
      fes = new FiniteElementSpace(mesh, fec);
      u = new GridFunction(fes);
      b = new LinearForm(fes);
   }
#else
   fes = new FiniteElementSpace(mesh, fec);
   u = new GridFunction(fes);
   b = new LinearForm(fes);
#endif
  *u = 0.0;

  if (!ess_bdr.Size()) {
    if (mesh->bdr_attributes.Size()) {
      ess_bdr.SetSize(mesh->bdr_attributes.Max());
      ess_bdr = 1;
    }
  }
}

void DiffusionSolver::Solve() {
  OperatorPtr A;
  Vector B, X;
  Array<int> ess_tdof_list;


#ifdef MFEM_USE_MPI
   if (parallel)
   {
      pfes->GetEssentialTrueDofs(ess_bdr,ess_tdof_list);
   }
   else
   {
      fes->GetEssentialTrueDofs(ess_bdr,ess_tdof_list);
   }
#else
   fes->GetEssentialTrueDofs(ess_bdr,ess_tdof_list);
#endif
  *u = 0.0;
  if (b) {
    delete b;
#ifdef MFEM_USE_MPI
      if (parallel)
      {
         b = new ParLinearForm(pfes);
      }
      else
      {
         b = new LinearForm(fes);
      }
#else
      b = new LinearForm(fes);
#endif
  }
  if (rhscf) {
    b->AddDomainIntegrator(new DomainLFIntegrator(*rhscf));
  }
  if (neumann_cf) {
    MFEM_VERIFY(neumann_bdr.Size(), "neumann_bdr attributes not provided");
    b->AddBoundaryIntegrator(new BoundaryLFIntegrator(*neumann_cf),
                             neumann_bdr);
  } else if (gradient_cf) {
    MFEM_VERIFY(neumann_bdr.Size(), "neumann_bdr attributes not provided");
    b->AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(*gradient_cf),
                             neumann_bdr);
  }

  b->Assemble();

  BilinearForm *a = nullptr;

#ifdef MFEM_USE_MPI
   if (parallel)
   {
      a = new ParBilinearForm(pfes);
   }
   else
   {
      a = new BilinearForm(fes);
   }
#else
   a = new BilinearForm(fes);
#endif
  a->AddDomainIntegrator(new DiffusionIntegrator(*diffcf));
  if (masscf) {
    a->AddDomainIntegrator(new MassIntegrator(*masscf));
  }
  a->Assemble();
  if (essbdr_cf) {
    u->ProjectBdrCoefficient(*essbdr_cf, ess_bdr);
  }
  a->FormLinearSystem(ess_tdof_list, *u, *b, A, X, B);

  CGSolver *cg = nullptr;
  Solver *M = nullptr;
#ifdef MFEM_USE_MPI
   if (parallel)
   {
      M = new HypreBoomerAMG;
      dynamic_cast<HypreBoomerAMG*>(M)->SetPrintLevel(0);
      cg = new CGSolver(pmesh->GetComm());
   }
   else
   {
      M = new GSSmoother((SparseMatrix&)(*A));
      cg = new CGSolver;
   }
#else
   M = new GSSmoother((SparseMatrix&)(*A));
   cg = new CGSolver;
#endif
  cg->SetRelTol(1e-12);
  cg->SetMaxIter(10000);
  cg->SetPrintLevel(0);
  cg->SetPreconditioner(*M);
  cg->SetOperator(*A);
  cg->Mult(B, X);
  delete M;
  delete cg;
  a->RecoverFEMSolution(X, *b, *u);
  delete a;
}

GridFunction *DiffusionSolver::GetFEMSolution() { return u; }

#ifdef MFEM_USE_MPI
ParGridFunction * DiffusionSolver::GetParFEMSolution()
{
   if (parallel)
   {
      return dynamic_cast<ParGridFunction*>(u);
   }
   else
   {
      MFEM_ABORT("Wrong code path. Call GetFEMSolution");
      return nullptr;
   }
}
#endif

void DiffusionSolver::ResetFEM() {
  delete u;
  u = nullptr;
  delete fes;
  fes = nullptr;
  delete fec;
  fec = nullptr;
  delete b;
}

DiffusionSolver::~DiffusionSolver() { ResetFEM(); }


} // namespace mfem