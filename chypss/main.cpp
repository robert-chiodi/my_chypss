//                       MFEM Example 16 - Parallel Version
//
// Compile with: make ex16p
//
// Sample runs:  mpirun -np 4 ex16p
//               mpirun -np 4 ex16p -m ../data/inline-tri.mesh
//               mpirun -np 4 ex16p -m ../data/disc-nurbs.mesh -tf 2
//               mpirun -np 4 ex16p -s 1 -a 0.0 -k 1.0
//               mpirun -np 4 ex16p -s 2 -a 1.0 -k 0.0
//               mpirun -np 8 ex16p -s 3 -a 0.5 -k 0.5 -o 4
//               mpirun -np 4 ex16p -s 14 -dt 1.0e-4 -tf 4.0e-2 -vs 40
//               mpirun -np 16 ex16p -m ../data/fichera-q2.mesh
//               mpirun -np 16 ex16p -m ../data/fichera-mixed.mesh
//               mpirun -np 16 ex16p -m ../data/escher-p2.mesh
//               mpirun -np 8 ex16p -m ../data/beam-tet.mesh -tf 10 -dt 0.1
//               mpirun -np 4 ex16p -m ../data/amr-quad.mesh -o 4 -rs 0 -rp 0
//               mpirun -np 4 ex16p -m ../data/amr-hex.mesh -o 2 -rs 0 -rp 0
//
// Description:  This example solves a time dependent nonlinear heat equation
//               problem of the form du/dt = C(u), with a non-linear diffusion
//               operator C(u) = \nabla \cdot (\kappa + \alpha u) \nabla u.
//
//               The example demonstrates the use of nonlinear operators (the
//               class ConductionOperator defining C(u)), as well as their
//               implicit time integration. Note that implementing the method
//               ConductionOperator::ImplicitSolve is the only requirement for
//               high-order implicit (SDIRK) time integration. Optional saving
//               with ADIOS2 (adios2.readthedocs.io) is also illustrated.
//
//               We recommend viewing examples 2, 9 and 10 before viewing this
//               example.

#include <cfloat>
#include <fstream>
#include <iostream>
#include "mfem.hpp"

using namespace mfem;

static ODESolver* ODESolverSelect(const int a_ode_solver_type) {
  switch (a_ode_solver_type) {
    // Implicit L-stable methods
    case 1:
      return new BackwardEulerSolver();
      break;
    case 2:
      return new SDIRK23Solver(2);
      break;
    case 3:
      return new SDIRK33Solver();
      break;
    // Explicit methods
    case 11:
      return new ForwardEulerSolver();
      break;
    case 12:
      return new RK2Solver(0.5);
      break;  // midpoint method
    case 13:
      return new RK3SSPSolver();
      break;
    case 14:
      return new RK4Solver();
      break;
    case 15:
      return new GeneralizedAlphaSolver(0.5);
      break;
    // Implicit A-stable methods (not L-stable)
    case 22:
      return new ImplicitMidpointSolver();
      break;
    case 23:
      return new SDIRK23Solver();
      break;
    case 24:
      return new SDIRK34Solver();
      break;
    default:
      std::cout << "Unknown ODE solver type: " << a_ode_solver_type << '\n';
      return nullptr;
  }
}

/** After spatial discretization, the conduction model can be written as:
 *
 *     du/dt = M^{-1}(-Ku)
 *
 *  where u is the vector representing the temperature, M is the mass matrix,
 *  and K is the diffusion operator with diffusivity depending on u:
 *  (\kappa + \alpha u).
 *
 *  Class ConductionOperator represents the right-hand side of the above ODE.
 */
class ConductionOperatorBase : public TimeDependentOperator {
 protected:
  ParFiniteElementSpace& fespace;
  Array<int> ess_tdof_list;  // this list remains empty for pure Neumann b.c.

  ParBilinearForm* M;
  ParBilinearForm* K;

  double current_dt;

  CGSolver M_solver;  // Krylov solver for inverting the mass matrix M

  double alpha, kappa;

  mutable Vector z;  // auxiliary vector
 public:
  ConductionOperatorBase(ParFiniteElementSpace& f, double alpha, double kappa,
                         const Vector& u);

  virtual void Mult(const Vector& u, Vector& du_dt) const = 0;
  virtual void ImplicitSolve(const double dt, const Vector& u, Vector& k) = 0;
  virtual void SetParameters(const Vector& u) = 0;

  virtual ~ConductionOperatorBase();
};

class ConductionOperator : public ConductionOperatorBase {
 protected:
  HypreParMatrix Mmat;
  HypreParMatrix Kmat;
  HypreParMatrix* T;     // T = M + dt K
  HypreSmoother M_prec;  // Preconditioner for the mass matrix M
  CGSolver T_solver;     // Implicit solver for T = M + dt K
  HypreSmoother T_prec;  // Preconditioner for the implicit solver
 public:
  ConductionOperator(ParFiniteElementSpace& f, double alpha, double kappa,
                     const Vector& u);

  virtual void Mult(const Vector& u, Vector& du_dt) const override final;
  /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
      This is the only requirement for high-order SDIRK implicit integration.*/
  virtual void ImplicitSolve(const double dt, const Vector& u,
                             Vector& k) override final;

  /// Update the diffusion BilinearForm K using the given true-dof vector `u`.
  virtual void SetParameters(const Vector& u) override final;

  virtual ~ConductionOperator() override final;
};

/* Same purpose as ConductionOperator above but rewritten
   to utilize Partial Assembly. Purpose of this is to test scaling
   difference when using Partial Assembly.
   NOTE: It is my understanding that only hex/quad meshes can
         utilize Partial Assembly due to their tensor-product nature.
 */
class ConductionOperatorPartial : public ConductionOperatorBase {
 protected:
  ParBilinearForm* PForm;
  CGSolver P_solver;  // Krylov solver for inverting the mass matrix M
  OperatorPtr P_op, M_op, K_op;
  OperatorJacobiSmoother* jacobi_smooth_p;
  OperatorJacobiSmoother* jacobi_smooth_m;
  ParGridFunction u_alpha_gf, dt_k_gf;

 public:
  ConductionOperatorPartial(ParFiniteElementSpace& f, double alpha,
                            double kappa, const Vector& u);

  virtual void Mult(const Vector& u, Vector& du_dt) const override final;
  /** Solve the Backward-Euler equation: k = f(u + dt*k, t), for the unknown k.
      This is the only requirement for high-order SDIRK implicit integration.*/
  virtual void ImplicitSolve(const double dt, const Vector& u,
                             Vector& k) override final;

  /// Update the diffusion BilinearForm K using the given true-dof vector `u`.
  virtual void SetParameters(const Vector& u) override final;

  virtual ~ConductionOperatorPartial() override final;
};

double InitialTemperature(const Vector& x);

int main(int argc, char* argv[]) {
  // 1. Initialize MPI.
  int num_procs, myid;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // 2. Parse command-line options.
  const char* mesh_file = "../data/star.mesh";
  int step_max = 5;
  int total_ref_levels = 0;
  int order = 2;
  int ode_solver_type = 3;
  double t_final = 0.5;
  double dt = 1.0e-5;
  double alpha = 1.0e-2;
  double kappa = 0.5;
  bool visualization = false;
  bool visit = false;
  int vis_steps = 5;
  bool adios2 = false;
  bool use_partial_assembly = false;
  bool comp = false;

  int precision = 8;
  std::cout.precision(precision);

  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&total_ref_levels, "-rt", "--refine-total",
                 "Number of times to refine the mesh uniformly.");
  args.AddOption(&order, "-o", "--order",
                 "Order (degree) of the finite elements.");
  args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                 "ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3,\n\t"
                 "\t   11 - Forward Euler, 12 - RK2, 13 - RK3 SSP, 14 - RK4.");
  args.AddOption(&step_max, "-sm", "--step-max",
                 "Number of steps to compute. Stops at t_final or step_max");
  args.AddOption(&t_final, "-tf", "--t-final", "Final time; start time is 0.");
  args.AddOption(&dt, "-dt", "--time-step", "Time step.");
  args.AddOption(&alpha, "-a", "--alpha", "Alpha coefficient.");
  args.AddOption(&kappa, "-k", "--kappa", "Kappa coefficient offset.");
  args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                 "--no-visualization",
                 "Enable or disable GLVis visualization.");
  args.AddOption(&visit, "-visit", "--visit-datafiles", "-no-visit",
                 "--no-visit-datafiles",
                 "Save data files for VisIt (visit.l1lnl.gov) visualization.");
  args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                 "Visualize every n-th timestep.");
  args.AddOption(&adios2, "-adios2", "--adios2-streams", "-no-adios2",
                 "--no-adios2-streams", "Save data using adios2 streams.");
  args.AddOption(&use_partial_assembly, "-use-pa", "--use-partial-assembly",
                 "-no-pa", "--no_partial_assembly",
                 "Construct operator using Partial Assembly from MFEM.");
  args.AddOption(
      &comp, "-comp", "--compare-solutions", "-no-comp",
      "--no-compare-solutions",
      "Compare traditional approach (from example) to more modified solution.");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cout);
    MPI_Finalize();
    return 1;
  }

  if (myid == 0 && visit) {
    args.PrintOptions(std::cout);
  }

  // 4. Define the ODE solver used for time integration. Several implicit
  //    singly diagonal implicit Runge-Kutta (SDIRK) methods, as well as
  //    explicit Runge-Kutta methods are available.
  ODESolver* ode_solver = ODESolverSelect(ode_solver_type);
  if (ode_solver == nullptr) {
    return 3;
  }

  // Create and refine mesh
  Mesh mesh(mesh_file, 1, 1);
  int dim = mesh.Dimension();
  auto refinements_to_go = total_ref_levels;
  while (refinements_to_go != 0 && mesh.GetNE() < num_procs * 100) {
    --refinements_to_go;
    mesh.UniformRefinement();
  }
  ParMesh pmesh(MPI_COMM_WORLD, mesh);
  mesh.Clear();
  for (int lev = 0; lev < refinements_to_go; lev++) {
    pmesh.UniformRefinement();
  }

  StopWatch timer_init;
  timer_init.Start();

  // 7. Define the vector finite element space representing the current and the
  //    initial temperature, u_ref.
  H1_FECollection fe_coll(order, dim);
  ParFiniteElementSpace fespace(&pmesh, &fe_coll);

  // 8. Set the initial conditions for u. All boundaries are considered
  //    natural.
  FunctionCoefficient u_0(InitialTemperature);

  Vector u;
  ParGridFunction u_gf(&fespace);
  u_gf.ProjectCoefficient(u_0);
  u_gf.GetTrueDofs(u);
  Vector u_comp = u;

  // 9. Initialize the conduction operator and the VisIt visualization.
  ConductionOperatorBase* oper = nullptr;
  if (use_partial_assembly) {
    oper = new ConductionOperatorPartial(fespace, alpha, kappa, u);
  } else {
    oper = new ConductionOperator(fespace, alpha, kappa, u);
  }

  ConductionOperator* comp_oper = nullptr;
  if (comp) {
    MFEM_VERIFY(use_partial_assembly, "");
    comp_oper = new ConductionOperator(fespace, alpha, kappa, u_comp);
  }

  VisItDataCollection* visit_dc = nullptr;
  if (visit) {
    // Setup mesh/variable for export
    u_gf.SetFromTrueDofs(u);
    // Visit specific database for export
    visit_dc = new VisItDataCollection("Example16-Parallel", &pmesh);
    visit_dc->RegisterField("temperature", &u_gf);
    visit_dc->SetCycle(0);
    visit_dc->SetTime(0.0);
    visit_dc->Save();
  }

  // 10. Perform time-integration (looping over the time iterations, ti, with a
  //     time-step dt).
  ODESolver* comp_solver = nullptr;
  if (comp) {
    comp_solver = ODESolverSelect(ode_solver_type);
    comp_solver->Init(*comp_oper);
  }
  ode_solver->Init(*oper);

  double t = 0.0;

  bool last_step = false;

  MPI_Barrier(MPI_COMM_WORLD);
  timer_init.Stop();

  StopWatch inner_timer;
  inner_timer.Start();
  double old_t, old_dt;
  for (int ti = 1; !last_step; ti++) {
    if (t + dt >= t_final - dt / 2) {
      last_step = true;
    }
    if (ti == step_max) {
      last_step = true;
    }

    if (comp) {
      old_t = t;
      old_dt = dt;
      u_comp = u;
    }
    ode_solver->Step(u, t, dt);

    if (comp) {
      comp_solver->Step(u_comp, old_t, old_dt);
      double L2_comp_local = 0.0;
      double Linf_comp_local = -DBL_MAX;
      for (int i = 0; i < u_comp.Size(); ++i) {
        L2_comp_local += u_comp(i) * u_comp(i);
        Linf_comp_local = Linf_comp_local > std::fabs(u_comp(i))
                              ? Linf_comp_local
                              : std::fabs(u_comp(i));
      }
      double L2_partial_local = 0.0;
      double Linf_partial_local = -DBL_MAX;
      for (int i = 0; i < u_comp.Size(); ++i) {
        L2_partial_local += u(i) * u(i);
        Linf_partial_local = Linf_partial_local > std::fabs(u(i))
                                 ? Linf_partial_local
                                 : std::fabs(u(i));
      }
      double L2_comp, L2_partial;
      double Linf_comp, Linf_partial;
      MPI_Allreduce(&L2_comp_local, &L2_comp, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);
      MPI_Allreduce(&L2_partial_local, &L2_partial, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);
      MPI_Allreduce(&Linf_comp_local, &Linf_comp, 1, MPI_DOUBLE, MPI_MAX,
                    MPI_COMM_WORLD);
      MPI_Allreduce(&Linf_partial_local, &Linf_partial, 1, MPI_DOUBLE, MPI_MAX,
                    MPI_COMM_WORLD);
      if (myid == 0) {
        std::cout << "Example L2: " << std::sqrt(L2_comp) << "  Partial L2 "
                  << std::sqrt(L2_partial) << '\n';
        std::cout << "Example Linf " << Linf_comp << "  Partial Linf "
                  << Linf_partial << '\n'
                  << std::endl;
      }
    }

    if (visit) {
      if (last_step || (ti % vis_steps) == 0) {
        if (myid == 0) {
          std::cout << "step " << ti << ", t = " << t << std::endl;
        }

        u_gf.SetFromTrueDofs(u);

        visit_dc->SetCycle(ti);
        visit_dc->SetTime(t);
        visit_dc->Save();
      }
    }
    oper->SetParameters(u);
    if (comp) {
      comp_oper->SetParameters(u_comp);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  inner_timer.Stop();

  StopWatch end_timer;
  end_timer.Start();

  // Delete heap allocated objects
  delete ode_solver;
  delete visit_dc;
  delete oper;
  MPI_Barrier(MPI_COMM_WORLD);
  end_timer.Stop();

  int fe_size = fespace.GlobalTrueVSize();
  int local_number_of_elements = pmesh.GetNE();
  int global_number_of_elements;
  MPI_Allreduce(&local_number_of_elements, &global_number_of_elements, 1,
                MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (myid == 0) {
    std::cout << "FE Space Size: " << fe_size << '\n';
    std::cout << "Number of elements: " << global_number_of_elements << '\n';
    std::cout << "Initialization time: "
              << timer_init.RealTime() + end_timer.RealTime() << '\n';
    std::cout << "Inner iterations " << inner_timer.RealTime() << '\n';
    std::cout << "Timed at resolution: " << timer_init.Resolution()
              << std::endl;
  }

  MPI_Finalize();

  return 0;
}

ConductionOperatorBase::ConductionOperatorBase(ParFiniteElementSpace& f,
                                               double al, double kap,
                                               const Vector& u)
    : TimeDependentOperator(f.GetTrueVSize(), 0.0),
      fespace(f),
      M(NULL),
      K(NULL),
      current_dt(0.0),
      M_solver(f.GetComm()),
      z(height) {}

ConductionOperatorBase::~ConductionOperatorBase() {
  delete M;
  delete K;
}

ConductionOperator::ConductionOperator(ParFiniteElementSpace& f, double al,
                                       double kap, const Vector& u)
    : ConductionOperatorBase(f, al, kap, u), T(NULL), T_solver(f.GetComm()) {
  const double rel_tol = 1.0e-10;
  M = new ParBilinearForm(&fespace);
  M->AddDomainIntegrator(new MassIntegrator());
  M->Assemble(0);  // keep sparsity pattern of M and K the same
  M->FormSystemMatrix(ess_tdof_list, Mmat);

  M_solver.iterative_mode = false;
  M_solver.SetRelTol(rel_tol);
  M_solver.SetAbsTol(0.0);
  M_solver.SetMaxIter(100);
  M_solver.SetPrintLevel(0);
  M_prec.SetType(HypreSmoother::Jacobi);
  M_solver.SetPreconditioner(M_prec);
  M_solver.SetOperator(Mmat);

  alpha = al;
  kappa = kap;

  T_solver.iterative_mode = false;
  T_solver.SetRelTol(rel_tol);
  T_solver.SetAbsTol(0.0);
  T_solver.SetMaxIter(100);
  T_solver.SetPrintLevel(0);
  T_solver.SetPreconditioner(T_prec);

  SetParameters(u);
}

void ConductionOperator::Mult(const Vector& u, Vector& du_dt) const {
  // Compute:
  //    du_dt = M^{-1}*-K(u)
  // for du_dt
  Kmat.Mult(u, z);
  z.Neg();  // z = -z
  M_solver.Mult(z, du_dt);
}

void ConductionOperator::ImplicitSolve(const double dt, const Vector& u,
                                       Vector& du_dt) {
  // Solve the equation:
  //    du_dt = M^{-1}*[-K(u + dt*du_dt)]
  // for du_dt
  if (!T) {
    T = Add(1.0, Mmat, dt, Kmat);
    current_dt = dt;
    T_solver.SetOperator(*T);
  }
  MFEM_VERIFY(dt == current_dt, "");  // SDIRK methods use the same dt
  Kmat.Mult(u, z);
  z.Neg();
  T_solver.Mult(z, du_dt);
}

void ConductionOperator::SetParameters(const Vector& u) {
  ParGridFunction u_alpha_gf(&fespace);
  u_alpha_gf.SetFromTrueDofs(u);
  for (int i = 0; i < u_alpha_gf.Size(); i++) {
    u_alpha_gf(i) = kappa + alpha * u_alpha_gf(i);
  }

  delete K;
  K = new ParBilinearForm(&fespace);

  GridFunctionCoefficient u_coeff(&u_alpha_gf);

  K->AddDomainIntegrator(new DiffusionIntegrator(u_coeff));
  K->Assemble(0);  // keep sparsity pattern of M and K the same
  K->FormSystemMatrix(ess_tdof_list, Kmat);
  delete T;
  T = NULL;  // re-compute T on the next ImplicitSolve
}

ConductionOperator::~ConductionOperator(void) { delete T; }

ConductionOperatorPartial::ConductionOperatorPartial(ParFiniteElementSpace& f,
                                                     double al, double kap,
                                                     const Vector& u)
    : ConductionOperatorBase(f, al, kap, u),
      PForm(nullptr),
      P_solver(f.GetComm()),
      jacobi_smooth_p(nullptr),
      jacobi_smooth_m(nullptr),
      u_alpha_gf(&fespace),
      dt_k_gf(&fespace) {
  const double rel_tol = 1.0e-10;
  M = new ParBilinearForm(&fespace);
  M->SetAssemblyLevel(AssemblyLevel::PARTIAL);
  M->AddDomainIntegrator(new MassIntegrator());
  M->Assemble(0);  // keep sparsity pattern of M and K the same
  M->FormSystemMatrix(ess_tdof_list, M_op);

  M_solver.iterative_mode = false;
  M_solver.SetRelTol(rel_tol);
  M_solver.SetAbsTol(0.0);
  M_solver.SetMaxIter(100);
  M_solver.SetPrintLevel(0);
  jacobi_smooth_m = new OperatorJacobiSmoother(*M, ess_tdof_list);
  M_solver.SetPreconditioner(*jacobi_smooth_m);
  M_solver.SetOperator(*M_op);

  P_solver.iterative_mode = false;
  P_solver.SetRelTol(rel_tol);
  P_solver.SetAbsTol(0.0);
  P_solver.SetMaxIter(100);
  P_solver.SetPrintLevel(0);

  alpha = al;
  kappa = kap;
  SetParameters(u);
}

void ConductionOperatorPartial::Mult(const Vector& u, Vector& du_dt) const {
  // Compute:
  //    du_dt = M^{-1}*-K(u)
  // for du_dt
  K_op->Mult(u, z);
  z.Neg();  // z = -z
  M_solver.Mult(z, du_dt);
}

void ConductionOperatorPartial::ImplicitSolve(const double dt, const Vector& u,
                                              Vector& du_dt) {
  // Solve the equation:
  //    du_dt = M^{-1}*[-K(u + dt*du_dt)]
  // for du_dt
  if (PForm == nullptr) {
    PForm = new ParBilinearForm(&fespace);
    PForm->SetAssemblyLevel(AssemblyLevel::PARTIAL);
    PForm->AddDomainIntegrator(new MassIntegrator());
    dt_k_gf.Distribute(u_alpha_gf);
    for (int i = 0; i < dt_k_gf.Size(); i++) {
      dt_k_gf(i) *= dt;
    }
    GridFunctionCoefficient coef_dt_kappa(&dt_k_gf);
    PForm->AddDomainIntegrator(new DiffusionIntegrator(coef_dt_kappa));
    PForm->Assemble(0);
    PForm->FormSystemMatrix(ess_tdof_list, P_op);
    if (UsesTensorBasis(fespace)) {
      delete jacobi_smooth_p;
      jacobi_smooth_p = new OperatorJacobiSmoother(*PForm, ess_tdof_list);
      P_solver.SetPreconditioner(*jacobi_smooth_p);
    }
    P_solver.SetOperator(*P_op);
    current_dt = dt;
  }
  MFEM_VERIFY(dt == current_dt, "");  // SDIRK methods use the same dt
  K_op->Mult(u, z);
  z.Neg();
  P_solver.Mult(z, du_dt);
}

void ConductionOperatorPartial::SetParameters(const Vector& u) {
  u_alpha_gf.SetFromTrueDofs(u);
  for (int i = 0; i < u_alpha_gf.Size(); i++) {
    u_alpha_gf(i) = kappa + alpha * u_alpha_gf(i);
  }

  delete K;
  K = new ParBilinearForm(&fespace);
  K->SetAssemblyLevel(AssemblyLevel::PARTIAL);
  GridFunctionCoefficient grid_coeff(&u_alpha_gf);
  K->AddDomainIntegrator(new DiffusionIntegrator(grid_coeff));
  K->Assemble(0);  // keep sparsity pattern of M and K the same
  K->FormSystemMatrix(ess_tdof_list, K_op);

  delete PForm;
  PForm =
      nullptr;  // Delete here to be reset on first entrance to SolveImplicit
}

ConductionOperatorPartial::~ConductionOperatorPartial(void) {
  delete PForm;
  delete jacobi_smooth_p;
  delete jacobi_smooth_m;
}

double InitialTemperature(const Vector& x) {
  Vector x_center = x;
  x_center -= 0.5;
  if (x_center.Norml2() < 0.25) {
    return 2.0;
  } else {
    return 0.0;
  }
}
