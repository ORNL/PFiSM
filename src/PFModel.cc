#include "PFModel.h"
#include "tools.h"

#include "samrai_internal/CellPoissonFACSolver.h"

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/math/HierarchyDataOpsReal.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/math/PatchCellDataOpsReal.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/VariableDatabase.h"

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

using namespace std;

extern "C" {
void SAMRAI_F77_FUNC(comprhs3d, COMPRHS3D)(const int&, const int&, const int&,
                                           const int&, const int&, const int&,
                                           const double*, const double*,
                                           const int&, const double&,
                                           const double*, const double&,
                                           const double&, double*);
void SAMRAI_F77_FUNC(comprhsphase3d, COMPRHSPHASE3D)(
    const int&, const int&, const int&, const int&, const int&, const int&,
    const double*, const int&, const double*, const int&, const double*,
    const double&, const double&, const double&, const double&, const double&,
    const double*);
void SAMRAI_F77_FUNC(compcforphase,
                     COMPCFORPHASE)(const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const double*, const int&, const double&,
                                    const double&, const double&, const double*,
                                    const int&);
}

PFModel::PFModel(const string& object_name, const tbox::Dimension& dim,
                 std::shared_ptr<CellPoissonFACSolver> fac_solver_temperature,
                 std::shared_ptr<CellPoissonFACSolver> fac_solver_phase,
                 std::shared_ptr<tbox::Database> input_db,
                 std::shared_ptr<geom::CartesianGridGeometry> grid_geom)
    : xfer::RefinePatchStrategy(),
      xfer::CoarsenPatchStrategy(),
      d_object_name(object_name),
      d_dim(dim),
      // we use 3 fields of depth 1
      d_temperature_var(new pdat::CellVariable<double>(dim, "temperature", 1)),
      d_phase_var(new pdat::CellVariable<double>(dim, "phase", 1)),
      d_cfield_phase_var(new pdat::CellVariable<double>(dim, "cfield", 1)),
      d_vol_var(new pdat::CellVariable<double>(dim, "vol", 1)),
      d_temperature_component(0),
      d_phase_component(1),
      d_FAC_solver_temperature(fac_solver_temperature),
      d_FAC_solver_phase(fac_solver_phase),
      d_grid_geometry(grid_geom),
      d_number_rhs_eval(0),
      d_number_precond_setup(0),
      d_number_precond_solve(0),
      d_temperature_bc_helper(new solv::CartesianRobinBcHelper(dim, "TempBC")),
      d_phase_bc_helper(new solv::CartesianRobinBcHelper(dim, "PhaseBC"))
{
   // set up variables and contexts
   // 'scratch' variables have 1 layer of ghost cells
   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   d_cur_cxt = variable_db->getContext("CURRENT");
   d_scr_cxt = variable_db->getContext("SCRATCH");

   d_temperature_cur_id =
       variable_db->registerVariableAndContext(d_temperature_var, d_cur_cxt,
                                               hier::IntVector(d_dim, 0));
   d_temperature_scr_id =
       variable_db->registerVariableAndContext(d_temperature_var, d_scr_cxt,
                                               hier::IntVector(d_dim, 1));

   d_phase_cur_id =
       variable_db->registerVariableAndContext(d_phase_var, d_cur_cxt,
                                               hier::IntVector(d_dim, 0));
   d_phase_scr_id =
       variable_db->registerVariableAndContext(d_phase_var, d_scr_cxt,
                                               hier::IntVector(d_dim, 1));

   d_cfield_phase_id =
       variable_db->registerVariableAndContext(d_cfield_phase_var, d_cur_cxt,
                                               hier::IntVector(d_dim, 0));

   d_vol_id =
       variable_db->registerVariableAndContext(d_vol_var, d_cur_cxt,
                                               hier::IntVector(d_dim, 0));

   d_current_time = 0.;

   d_print_solver_info = false;

   // Initialize object with data read from given input/restart databases.
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart) {
      getFromRestart();
   }
   getFromInput(input_db, is_from_restart);

   std::shared_ptr<tbox::Database> bc_db(
       input_db->getDatabase("BoundaryConditions"));

   d_temperature_bc_coeffs =
       new solv::LocationIndexRobinBcCoefs(d_dim, "TemperatureBCcoeffs",
                                           bc_db->getDatabase("Temperature"));

   d_temperature_bc_helper->setTargetDataId(d_temperature_scr_id);
   d_temperature_bc_helper->setCoefImplementation(d_temperature_bc_coeffs);

   d_phase_bc_coeffs =
       new solv::LocationIndexRobinBcCoefs(d_dim, "PhaseBCcoeffs",
                                           bc_db->getDatabase("Phase"));

   d_phase_bc_helper->setTargetDataId(d_phase_scr_id);
   d_phase_bc_helper->setCoefImplementation(d_phase_bc_coeffs);

   // Boundary conditions for FAC solvers should be homogeneous
   // since solver computes corrections to current guess
   d_temperature_bc_corr_coeffs =
       new solv::LocationIndexRobinBcCoefs(d_dim, "BCcorrcoeffs",
                                           bc_db->getDatabase("Temperature"));
   for (int i = 0; i < d_dim.getValue() * 2; i++) {
      double a, b, g;
      d_temperature_bc_corr_coeffs->getCoefficients(i, a, b, g);
      g = 0.;
      d_temperature_bc_corr_coeffs->setRawCoefficients(i, a, b, g);
   }

   d_FAC_solver_temperature->setBcObject(d_temperature_bc_corr_coeffs);

   d_phase_bc_corr_coeffs =
       new solv::LocationIndexRobinBcCoefs(d_dim, "PhaseBCcorrcoeffs",
                                           bc_db->getDatabase("Phase"));
   for (int i = 0; i < d_dim.getValue() * 2; i++) {
      double a, b, g;
      d_phase_bc_corr_coeffs->getCoefficients(i, a, b, g);
      g = 0.;
      d_phase_bc_corr_coeffs->setRawCoefficients(i, a, b, g);
   }

   d_FAC_solver_phase->setBcObject(d_phase_bc_corr_coeffs);

   tbox::TimerManager* tman = tbox::TimerManager::getManager();
   t_rhs_timer = tman->getTimer("PFiSM::rhs");
   t_precondset_timer = tman->getTimer("PFiSM::precondset");
   t_precondsolve_timer = tman->getTimer("PFiSM::precondsolve");
   t_factemperature_timer = tman->getTimer("PFiSM::factemperature");
   t_facphase_timer = tman->getTimer("PFiSM::facphase");
   t_factempinit_timer = tman->getTimer("PFiSM::factempinit");
   t_facphaseinit_timer = tman->getTimer("PFiSM::facphaseinit");
}

PFModel::~PFModel()
{
   std::shared_ptr<solv::SAMRAIVectorReal<double> > samvect =
       solv::Sundials_SAMRAIVector::getSAMRAIVector(d_solution_vector);
   solv::Sundials_SAMRAIVector::destroySundialsVector(d_solution_vector);

   samvect->freeVectorComponents();
   samvect.reset();
}

/*************************************************************************
 * Methods inherited from SAMRAI::mesh::StandardTagAndInitStrategy.
 ************************************************************************/
void PFModel::initializeLevelData(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const int level_number, const double time, const bool can_be_refined,
    const bool initial_time, const std::shared_ptr<hier::PatchLevel>& old_level,
    const bool allocate_data)
{
   NULL_USE(hierarchy);
   NULL_USE(level_number);
   NULL_USE(time);
   NULL_USE(can_be_refined);
   NULL_USE(initial_time);
   NULL_USE(time);
   NULL_USE(old_level);
   NULL_USE(allocate_data);

   // Empty because initialization is done by setInitialConditions().
   // Data that is not managed inside the SAMRAI CVODESolver class
   // and that must be set on the level should be initialized here
}

void PFModel::resetHierarchyConfiguration(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const int coarsest_level, const int finest_level)
{
   NULL_USE(hierarchy);
   NULL_USE(coarsest_level);
   NULL_USE(finest_level);
}

void PFModel::applyGradientDetector(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const int level_number, const double time, const int tag_index,
    const bool initial_time, const bool uses_richardson_extrapolation_too)
{
   NULL_USE(time);
   NULL_USE(initial_time);
   NULL_USE(uses_richardson_extrapolation_too);

   std::shared_ptr<hier::PatchLevel> level(
       hierarchy->getPatchLevel(level_number));

   for (hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {
      const std::shared_ptr<hier::Patch>& patch = *p;

      std::shared_ptr<pdat::CellData<int> > tag_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
              patch->getPatchData(tag_index)));
      assert(tag_data);

      // dumb implementation:
      // that tags all cells to refine everywhere
      tag_data->fillAll(1);
   }
}

/*************************************************************************
 * Methods inherited from RefinePatchStrategy.
 ***********************************************************************/
void PFModel::setPhysicalBoundaryConditions(
    hier::Patch& patch, const double time,
    const hier::IntVector& ghost_width_to_fill)
{
   d_temperature_bc_helper->setPhysicalBoundaryConditions(patch, time,
                                                          ghost_width_to_fill);

   d_phase_bc_helper->setPhysicalBoundaryConditions(patch, time,
                                                    ghost_width_to_fill);
}

void PFModel::preprocessRefine(hier::Patch& fine, const hier::Patch& coarse,
                               const hier::Box& fine_box,
                               const hier::IntVector& ratio)
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_box);
   NULL_USE(ratio);
}

void PFModel::postprocessRefine(hier::Patch& fine, const hier::Patch& coarse,
                                const hier::Box& fine_box,
                                const hier::IntVector& ratio)
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_box);
   NULL_USE(ratio);
}

/*************************************************************************
 * Methods inherited from CoarsenPatchStrategy.
 ************************************************************************/
void PFModel::preprocessCoarsen(hier::Patch& coarse, const hier::Patch& fine,
                                const hier::Box& coarse_box,
                                const hier::IntVector& ratio)
{
   NULL_USE(coarse);
   NULL_USE(fine);
   NULL_USE(coarse_box);
   NULL_USE(ratio);
}

void PFModel::postprocessCoarsen(hier::Patch& coarse, const hier::Patch& fine,
                                 const hier::Box& coarse_box,
                                 const hier::IntVector& ratio)
{
   NULL_USE(coarse);
   NULL_USE(fine);
   NULL_USE(coarse_box);
   NULL_USE(ratio);
}

/*************************************************************************
 * Methods inherited from CVODEAbstractFunction
 ************************************************************************/
int PFModel::evaluateRHSFunction(double time, solv::SundialsAbstractVector* y,
                                 solv::SundialsAbstractVector* y_dot)
{
   t_rhs_timer->start();

   // Convert Sundials vectors to SAMRAI vectors
   std::shared_ptr<solv::SAMRAIVectorReal<double> > y_samvect(
       solv::Sundials_SAMRAIVector::getSAMRAIVector(y));
   std::shared_ptr<solv::SAMRAIVectorReal<double> > y_dot_samvect(
       solv::Sundials_SAMRAIVector::getSAMRAIVector(y_dot));

   std::shared_ptr<hier::PatchHierarchy> hierarchy(
       y_samvect->getPatchHierarchy());

   // fill ghost values for temperature and phase variables
   std::shared_ptr<xfer::RefineAlgorithm> bdry_fill_alg(
       new xfer::RefineAlgorithm());
   std::shared_ptr<hier::RefineOperator> refine_op(
       d_grid_geometry->lookupRefineOperator(d_temperature_var,
                                             "CONSERVATIVE_LINEAR_REFINE"));
   bdry_fill_alg->registerRefine(d_temperature_scr_id,  // dest
                                 y_samvect->getComponentDescriptorIndex(
                                     d_temperature_component),  // src
                                 d_temperature_scr_id,          // scratch
                                 refine_op);

   std::shared_ptr<hier::RefineOperator> phase_refine_op(
       d_grid_geometry->lookupRefineOperator(d_phase_var,
                                             "CONSERVATIVE_LINEAR_REFINE"));
   bdry_fill_alg->registerRefine(d_phase_scr_id,  // dest
                                 y_samvect->getComponentDescriptorIndex(
                                     d_phase_component),  // src
                                 d_phase_scr_id,          // scratch
                                 phase_refine_op);

   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
      if (!level->checkAllocated(d_temperature_scr_id)) {
         level->allocatePatchData(d_temperature_scr_id);
         level->allocatePatchData(d_phase_scr_id);
      }

      // a pointer to "this" tells the refine schedule to invoke
      // the setPhysicalBCs defined in this class.
      std::shared_ptr<xfer::RefineSchedule> bdry_fill_alg_schedule(
          bdry_fill_alg->createSchedule(level, ln - 1, hierarchy, this));

      bdry_fill_alg_schedule->fillData(time);
   }

   // now actually compute rhs
   int y_dot_phase_id =
       y_dot_samvect->getComponentDescriptorIndex(d_phase_component);
   int y_dot_temperature_id =
       y_dot_samvect->getComponentDescriptorIndex(d_temperature_component);

   evaluateRHSPhase(hierarchy, y_dot_phase_id);

   evaluateRHSTemperature(hierarchy, y_dot_temperature_id, y_dot_phase_id);

   // free up temporary array allocations
   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      hierarchy->getPatchLevel(ln)->deallocatePatchData(d_temperature_scr_id);
      hierarchy->getPatchLevel(ln)->deallocatePatchData(d_phase_scr_id);
   }

   d_current_time = time;
   ++d_number_rhs_eval;

   t_rhs_timer->stop();

   return 0;
}
void PFModel::evaluateRHSPhase(std::shared_ptr<hier::PatchHierarchy> hierarchy,
                               const int y_dot_phase_id)
{
   // Compute rhs for phase first, since it is used in rhs for temperature
   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

      for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end();
           ++ip) {
         const std::shared_ptr<hier::Patch>& patch = *ip;

         std::shared_ptr<pdat::CellData<double> > phase(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_phase_scr_id)));
         std::shared_ptr<pdat::CellData<double> > temperature(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_temperature_scr_id)));
         std::shared_ptr<pdat::CellData<double> > rhs(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(y_dot_phase_id)));
         assert(phase);
         assert(temperature);
         assert(rhs);

         const hier::Index ifirst(patch->getBox().lower());
         const hier::Index ilast(patch->getBox().upper());

         const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
             SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                                    hier::PatchGeometry>(
                 patch->getPatchGeometry()));
         assert(patch_geom);
         const double* dx = patch_geom->getDx();

         SAMRAI_F77_FUNC(comprhsphase3d, COMPRHSPHASE3D)
         (ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
          phase->getPointer(), phase->getGhostCellWidth()[0],
          temperature->getPointer(), temperature->getGhostCellWidth()[0], dx,
          d_mobility, d_well_height, d_epsilon, d_latent_heat, d_Tmelting,
          rhs->getPointer());

      }  // loop over patches
   }     // loop over levels
}

void PFModel::evaluateRHSTemperature(
    std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int y_dot_temperature_id, const int y_dot_phase_id)
{
   // now compute rhs for temperature
   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

      for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end();
           ++ip) {
         const std::shared_ptr<hier::Patch>& patch = *ip;

         std::shared_ptr<pdat::CellData<double> > y(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_temperature_scr_id)));
         std::shared_ptr<pdat::CellData<double> > rhs(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(y_dot_temperature_id)));
         std::shared_ptr<pdat::CellData<double> > phi_dot(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(y_dot_phase_id)));

         assert(y);
         assert(rhs);
         assert(phi_dot);

         const hier::Index ifirst(patch->getBox().lower());
         const hier::Index ilast(patch->getBox().upper());

         const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
             SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                                    hier::PatchGeometry>(
                 patch->getPatchGeometry()));
         assert(patch_geom);
         const double* dx = patch_geom->getDx();

         SAMRAI_F77_FUNC(comprhs3d, COMPRHS3D)
         (ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2), dx,
          y->getPointer(), y->getGhostCellWidth()[0], d_temperature_diffusion,
          phi_dot->getPointer(), d_latent_heat, d_cp, rhs->getPointer());

      }  // loop over patches
   }     // loop over levels
}

// Initialize FAC solvers.
// Should be called after regridding or when matrix coefficients are changed
void PFModel::initializeSolvers(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
   const int coarsest_solve_ln = 0;
   const int finest_solve_ln = hierarchy->getFinestLevelNumber();

   t_factempinit_timer->start();

   d_FAC_solver_temperature->deallocateSolverState();

   d_FAC_solver_temperature->initializeSolverState(d_temperature_scr_id,
                                                   d_temperature_cur_id,
                                                   hierarchy, coarsest_solve_ln,
                                                   finest_solve_ln);

   t_factempinit_timer->stop();

   t_facphaseinit_timer->start();

   d_FAC_solver_phase->deallocateSolverState();

   d_FAC_solver_phase->initializeSolverState(d_phase_scr_id, d_phase_cur_id,
                                             hierarchy, coarsest_solve_ln,
                                             finest_solve_ln);

   t_facphaseinit_timer->stop();
}

/*****************************************************************
 * Set up FAC preconditioner for Jacobian system.
 *****************************************************************/
int PFModel::CVSpgmrPrecondSet(
    double t,
    solv::SundialsAbstractVector* y,  // current value of variable vector,
                                      // the predicted value of y(t)
    solv::SundialsAbstractVector* fy, int jok, int* jcurPtr, double gamma)
{
   NULL_USE(t);
   NULL_USE(fy);
   NULL_USE(jok);
   NULL_USE(jcurPtr);

   t_precondset_timer->start();

   tbox::plog << "CVSpgmrPrecondSet..." << endl;

   std::shared_ptr<solv::SAMRAIVectorReal<double> > y_samvect(
       solv::Sundials_SAMRAIVector::getSAMRAIVector(y));

   std::shared_ptr<hier::PatchHierarchy> hierarchy(
       y_samvect->getPatchHierarchy());

   int y_indx = y_samvect->getComponentDescriptorIndex(d_phase_component);
   PrecondSetPhase(hierarchy, y_indx, gamma);

   PrecondSetTemperature(gamma);

   initializeSolvers(hierarchy);

   ++d_number_precond_setup;

   t_precondset_timer->stop();

   // assume success and return 0
   return 0;
}

void PFModel::PrecondSetPhase(std::shared_ptr<hier::PatchHierarchy> hierarchy,
                              int y_indx, const double gamma)
{
   // Construct coarsen algorithm to fill interiors on coarser levels
   // with solution on finer level.
   xfer::CoarsenAlgorithm fill_phase_interior_on_coarser(d_dim);
   std::shared_ptr<hier::CoarsenOperator> coarsen_op(
       d_grid_geometry->lookupCoarsenOperator(d_phase_var,
                                              "CONSERVATIVE_COARSEN"));

   fill_phase_interior_on_coarser.registerCoarsen(y_indx, y_indx, coarsen_op);

   for (int amr_level = hierarchy->getFinestLevelNumber(); amr_level >= 0;
        --amr_level) {
      std::shared_ptr<hier::PatchLevel> level(
          hierarchy->getPatchLevel(amr_level));

      // Construct a coarsen schedule for all levels larger than coarsest,
      // and fill interiors of solution vector on coarser levels using fine
      // data.
      if (amr_level > 0) {
         std::shared_ptr<hier::PatchLevel> coarser_level(
             hierarchy->getPatchLevel(amr_level - 1));

         std::shared_ptr<xfer::CoarsenSchedule>
             fill_phase_interior_on_coarser_sched(
                 fill_phase_interior_on_coarser.createSchedule(coarser_level,
                                                               level));

         fill_phase_interior_on_coarser_sched->coarsenData();
      }

   }  // level loop

   // setup phase FAC solver
   setCforPhase(hierarchy, y_indx, gamma);

   d_FAC_solver_phase->setCPatchDataId(d_cfield_phase_id);
   d_FAC_solver_phase->setDConstant(-d_mobility * d_epsilon * d_epsilon);
}

void PFModel::PrecondSetTemperature(const double gamma)
{
   // setup temperature FAC solver
   d_FAC_solver_temperature->setCConstant(1.0 / gamma);
   d_FAC_solver_temperature->setDConstant(d_temperature_diffusion);
}

/*************************************************************************
 * Apply preconditioner.
 * r:  right-hand-side
 * z: solution
 * Assumes preconditioner has been setup already.
 * Return 0 if preconditioner fails; 1 otherwise.
 *************************************************************************/
int PFModel::CVSpgmrPrecondSolve(double t, solv::SundialsAbstractVector* y,
                                 solv::SundialsAbstractVector* fy,
                                 solv::SundialsAbstractVector* r,
                                 solv::SundialsAbstractVector* z, double gamma,
                                 double delta, int lr)
{
   NULL_USE(t);
   NULL_USE(y);
   NULL_USE(fy);
   NULL_USE(delta);
   NULL_USE(lr);

   t_precondsolve_timer->start();

   // plog<<"CVSpgmrPrecondSolve..."<<endl;

   // Convert passed-in CVODE vectors into SAMRAI vectors
   std::shared_ptr<solv::SAMRAIVectorReal<double> > r_samvect(
       solv::Sundials_SAMRAIVector::getSAMRAIVector(r));
   std::shared_ptr<solv::SAMRAIVectorReal<double> > z_samvect(
       solv::Sundials_SAMRAIVector::getSAMRAIVector(z));

   std::shared_ptr<hier::PatchHierarchy> hierarchy(
       r_samvect->getPatchHierarchy());

   int r0_indx =
       r_samvect->getComponentDescriptorIndex(d_temperature_component);
   int z0_indx =
       z_samvect->getComponentDescriptorIndex(d_temperature_component);

   bool converge0 = PrecondSolveTemperature(hierarchy, r0_indx, z0_indx, gamma);

   int r1_indx = r_samvect->getComponentDescriptorIndex(d_phase_component);
   int z1_indx = z_samvect->getComponentDescriptorIndex(d_phase_component);

   bool converge1 = PrecondSolvePhase(hierarchy, r1_indx, z1_indx, gamma);

   bool converge = (converge0 && converge1);

   if (d_print_solver_info) {
      double avg_convergence, final_convergence;
      d_FAC_solver_temperature->getConvergenceFactors(avg_convergence,
                                                      final_convergence);
      tbox::pout << "TEMPERATURE:" << '\n';
      tbox::pout << "   \t\t\tFinal Residual Norm: "
                 << d_FAC_solver_temperature->getResidualNorm() << '\n';
      tbox::pout << "   \t\t\tFinal Convergence Error: " << final_convergence
                 << endl;
      tbox::pout << "   \t\t\tFinal Convergence Rate: " << avg_convergence
                 << endl;
      tbox::pout << "PHASE:" << '\n';
      d_FAC_solver_phase->getConvergenceFactors(avg_convergence,
                                                final_convergence);
      tbox::pout << "   \t\t\tFinal Residual Norm: "
                 << d_FAC_solver_temperature->getResidualNorm() << endl;
      tbox::pout << "   \t\t\tFinal Convergence Error: " << final_convergence
                 << endl;
      tbox::pout << "   \t\t\tFinal Convergence Rate: " << avg_convergence
                 << endl;
   }

   ++d_number_precond_solve;

   int ret_val = (converge == true) ? 0 : 1;

   t_precondsolve_timer->stop();

   return ret_val;
}

bool PFModel::PrecondSolveTemperature(
    std::shared_ptr<hier::PatchHierarchy> hierarchy, int r0_indx, int z0_indx,
    double gamma)
{
   // We need to supply to the FAC solver a "version" of the z vector
   // that contains ghost cells.  The operations below allocate
   // on the patches a scratch context of the solution vector z and
   // fill it with 0 (initial guess for computed correction)
   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

      if (!level->checkAllocated(d_temperature_scr_id)) {
         level->allocatePatchData(d_temperature_scr_id);
      }

      for (hier::PatchLevel::iterator p(level->begin()); p != level->end();
           ++p) {

         const std::shared_ptr<hier::Patch>& patch = *p;

         // Scale RHS by 1/gamma
         math::PatchCellDataOpsReal<double> math_ops;
         std::shared_ptr<pdat::CellData<double> > r_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(r0_indx)));
         assert(r_data);
         math_ops.scale(r_data, 1.0 / gamma, r_data, r_data->getBox());

         // Set initial guess to 0
         std::shared_ptr<pdat::CellData<double> > z_scr_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_temperature_scr_id)));
         assert(z_scr_data);
         z_scr_data->fillAll(0.);
      }
   }

   t_factemperature_timer->start();

   bool converge0 =
       d_FAC_solver_temperature->solveSystem(d_temperature_scr_id, r0_indx);

   t_factemperature_timer->stop();

   // computed solutions are stored in scratch data space
   // copy them into zvector
   math::HierarchyCellDataOpsReal<double> cell_ops(hierarchy);
   cell_ops.copyData(z0_indx, d_temperature_scr_id, false);

   return converge0;
}

bool PFModel::PrecondSolvePhase(std::shared_ptr<hier::PatchHierarchy> hierarchy,
                                int r1_indx, int z1_indx, double gamma)
{
   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

      if (!level->checkAllocated(d_phase_scr_id)) {
         level->allocatePatchData(d_phase_scr_id);
      }

      for (hier::PatchLevel::iterator p(level->begin()); p != level->end();
           ++p) {

         const std::shared_ptr<hier::Patch>& patch = *p;

         // Scale RHS by 1/gamma
         math::PatchCellDataOpsReal<double> math_ops;
         std::shared_ptr<pdat::CellData<double> > r_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(r1_indx)));
         assert(r_data);
         math_ops.scale(r_data, 1.0 / gamma, r_data, r_data->getBox());

         //  Set initial guess to 0
         std::shared_ptr<pdat::CellData<double> > z_scr_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_phase_scr_id)));
         assert(z_scr_data);
         z_scr_data->fillAll(0.);
      }
   }

   t_facphase_timer->start();

   bool converge1 = d_FAC_solver_phase->solveSystem(d_phase_scr_id, r1_indx);

   t_facphase_timer->stop();

   // computed solutions are stored in scratch data space
   // copy them into zvector
   math::HierarchyCellDataOpsReal<double> cell_ops(hierarchy);
   cell_ops.copyData(z1_indx, d_phase_scr_id, false);

   return converge1;
}

/*************************************************************************
 * Methods specific to PFModel class.
 ************************************************************************/
void PFModel::setupSolutionVector(
    std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   // create SAMRAIVector with temperature and phase field
   std::shared_ptr<solv::SAMRAIVectorReal<double> > samvect(
       new solv::SAMRAIVectorReal<double>("solution", hierarchy, 0,
                                          hierarchy->getFinestLevelNumber()));

   samvect->addComponent(d_temperature_var, d_temperature_cur_id);
   samvect->addComponent(d_phase_var, d_phase_cur_id);

   samvect->allocateVectorData();

   d_solution_vector =
       solv::Sundials_SAMRAIVector::createSundialsVector(samvect);

   // Allocate memory for preconditioner variables.
   const int nlevels = hierarchy->getNumberOfLevels();
   for (int ln = 0; ln < nlevels; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
      assert(level);
      level->allocatePatchData(d_cfield_phase_id);
   }
}

/*************************************************************************
 * Set initial conditions for temperature and phase fields in CVODE solver
 *************************************************************************/
void PFModel::setInitialConditions()
{
   std::shared_ptr<solv::SAMRAIVectorReal<double> > init_samvect(
       solv::Sundials_SAMRAIVector::getSAMRAIVector(d_solution_vector));

   std::shared_ptr<hier::PatchHierarchy> hierarchy(
       init_samvect->getPatchHierarchy());

   const double zmin = d_grid_geometry->getXLower()[2];
   const double zmax = d_grid_geometry->getXUpper()[2];
   const double lz = zmax - zmin;
   tbox::pout << "Length in z-direction: " << lz << std::endl;

   for (int ln = 0; ln < hierarchy->getNumberOfLevels(); ++ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

      for (hier::PatchLevel::iterator p(level->begin()); p != level->end();
           ++p) {
         const std::shared_ptr<hier::Patch>& patch = *p;

         std::shared_ptr<geom::CartesianPatchGeometry> pg(
             SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                                    hier::PatchGeometry>(
                 patch->getPatchGeometry()));
         assert(pg);

         const double* h = pg->getDx();

         // Set initial conditions for temperature
         std::shared_ptr<pdat::CellData<double> > temperature_init(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 init_samvect->getComponentPatchData(d_temperature_component,
                                                     *patch)));
         assert(temperature_init);

         temperature_init->fillAll(d_temperature_init);

         // Set initial conditions for phase variable:
         // 1 for low z, 0 for high z
         std::shared_ptr<pdat::CellData<double> > phase_init(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 init_samvect->getComponentPatchData(d_phase_component,
                                                     *patch)));
         assert(phase_init);

         const hier::Box patch_box = patch->getBox();
         pdat::CellIterator ic(pdat::CellGeometry::begin(patch_box));
         pdat::CellIterator icend(pdat::CellGeometry::end(patch_box));

         for (; ic != icend; ++ic) {
            hier::IntVector icell = *ic;
            const double zval = (h[2] * (0.5 + icell[2]));

            if (zval > d_init_solid_fraction * lz)
               (*phase_init)(*ic) = 0.;
            else
               (*phase_init)(*ic) = 1.;
         }
      }
   }

   for (int ln = 0; ln < hierarchy->getNumberOfLevels(); ++ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
      level->allocatePatchData(d_vol_id);
   }

   computeVectorWeights(hierarchy, d_vol_id, 0,
                        hierarchy->getFinestLevelNumber());
}

/*************************************************************************
 * Set C in operator -div*D*grad*phi+C*phi
 *************************************************************************/
void PFModel::setCforPhase(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, int phase_id,
    const double gamma)
{
   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

      for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end();
           ++ip) {
         const std::shared_ptr<hier::Patch>& patch = *ip;

         std::shared_ptr<pdat::CellData<double> > phi(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(phase_id)));
         std::shared_ptr<pdat::CellData<double> > cfield(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_cfield_phase_id)));

         assert(phi);
         assert(cfield);

         const hier::Index ifirst(patch->getBox().lower());
         const hier::Index ilast(patch->getBox().upper());

         SAMRAI_F77_FUNC(compcforphase, COMPCFORPHASE)
         (ifirst(0), ilast(0), ifirst(1), ilast(1), ifirst(2), ilast(2),
          phi->getPointer(), phi->getGhostCellWidth()[0], d_mobility,
          d_well_height, gamma, cfield->getPointer(),
          cfield->getGhostCellWidth()[0]);
      }  // loop over patches
   }     // loop over levels
}

/*************************************************************************
 * Compute solid fraction
 *************************************************************************/
double PFModel::computeSolidFraction(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
   assert(d_vol_id >= 0);

   const double* lo = d_grid_geometry->getXLower();
   const double* up = d_grid_geometry->getXUpper();
   const double volume = (up[0] - lo[0]) * (up[1] - lo[1]) * (up[2] - lo[2]);

   math::HierarchyCellDataOpsReal<double> cmat(hierarchy);

   return cmat.integral(d_phase_cur_id, d_vol_id) / volume;
}


void PFModel::printCounters(const double final_time)
{
   tbox::plog << "\n\nEnd Timesteps - final time = " << final_time
              << "\n\tTotal number of RHS evaluations = " << d_number_rhs_eval
              << "\n\tTotal number of precond setups = "
              << d_number_precond_setup
              << "\n\tTotal number of precond solves = "
              << d_number_precond_solve << endl;
}

/*************************************************************************
 * Read data from input database.
 *************************************************************************/
void PFModel::getFromInput(std::shared_ptr<tbox::Database> input_db,
                           bool is_from_restart)
{
   NULL_USE(is_from_restart);

   std::shared_ptr<tbox::Database> temperature_db(
       input_db->getDatabase("Temperature"));

   std::shared_ptr<tbox::Database> pfm_db(input_db->getDatabase("PFM"));

   // read initial conditions
   d_init_solid_fraction = pfm_db->getDouble("initial_solid_fraction");
   d_temperature_init = temperature_db->getDouble("initial_value");

   // read PFM parameters
   d_epsilon = pfm_db->getDouble("epsilon");
   d_mobility = pfm_db->getDouble("mobility");
   d_well_height = pfm_db->getDouble("well_height");

   d_Tmelting = temperature_db->getDouble("melting");
   d_latent_heat = temperature_db->getDouble("latent_heat");
   d_cp = temperature_db->getDouble("cp");

   d_temperature_diffusion = temperature_db->getDouble("diffusion_coeff");

   // conversion from [J/mol] to [pJ/(mu m)^3]
   double molar_volume = input_db->getDouble("molar_volume");
   d_cp *= (1.e-6 / molar_volume);
   d_latent_heat = (1.e-6 / molar_volume);

   d_print_solver_info =
       input_db->getBoolWithDefault("print_solver_info", d_print_solver_info);
}

/*************************************************************************
 * Write data to  restart database.
 *************************************************************************/
void PFModel::putToRestart(
    const std::shared_ptr<tbox::Database>& restart_db) const
{
   NULL_USE(restart_db);

   assert(restart_db);
}

/*************************************************************************
 * Read data from restart database.
 *************************************************************************/
void PFModel::getFromRestart()
{
   std::shared_ptr<tbox::Database> root_db(
       tbox::RestartManager::getManager()->getRootDatabase());

   if (!root_db->isDatabase(d_object_name)) {
      TBOX_ERROR("Restart database corresponding to "
                 << d_object_name << " not found in the restart file.");
   }
   std::shared_ptr<tbox::Database> db(root_db->getDatabase(d_object_name));
}

/*************************************************************************
 * Register data to be plotted with VisIt
 *************************************************************************/
void PFModel::registerVisItDataWriter(
    std::shared_ptr<appu::VisItDataWriter> viz_writer)
{
   assert(viz_writer);
   d_visit_writer = viz_writer;

   if (d_visit_writer) {
      d_visit_writer->registerPlotQuantity("temperature", "SCALAR",
                                           d_temperature_cur_id, 0);
      d_visit_writer->registerPlotQuantity("phase", "SCALAR", d_phase_cur_id,
                                           0);
   }
}

void PFModel::printClassData(ostream& os) const
{
   fflush(stdout);
   os << "ptr PFModel = " << (PFModel*)this << endl;
   os << "d_object_name = " << d_object_name << endl;
   os << endl;
}
