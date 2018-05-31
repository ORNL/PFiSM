/*************************************************************************
 *
 * Adapted from SAMRAI/source/test/sundials
 *
 ************************************************************************/

#include "CVODEModel.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/math/HierarchyDataOpsReal.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/math/PatchCellDataOpsReal.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/VariableDatabase.h"

using namespace std;

extern "C" {
void SAMRAI_F77_FUNC(comprhs3d, COMPRHS3D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const double *,
   const double *,
   const double *, const double *, const double *,
   double *);
}

/*************************************************************************
 *
 * Constructor and Destructor for CVODEModel class.
 *
 ************************************************************************/
CVODEModel::CVODEModel(
   const string& object_name,
   const Dimension& dim,
   std::shared_ptr<CellPoissonFACSolver> fac_solver,
   std::shared_ptr<Database> input_db,
   std::shared_ptr<CartesianGridGeometry> grid_geom):
   RefinePatchStrategy(),
   CoarsenPatchStrategy(),
   d_object_name(object_name),
   d_dim(dim),
   d_temperature_var(new CellVariable<double>(dim, "temperature", 1)),
   d_FAC_solver(fac_solver),
   d_grid_geometry(grid_geom)
{
   /*
    * set up variables and contexts
    */
   VariableDatabase* variable_db = VariableDatabase::getDatabase();

   d_cur_cxt = variable_db->getContext("CURRENT");
   d_scr_cxt = variable_db->getContext("SCRATCH");

   d_temperature_cur_id = variable_db->registerVariableAndContext(d_temperature_var,
         d_cur_cxt,
         IntVector(d_dim, 0));
   d_temperature_scr_id = variable_db->registerVariableAndContext(d_temperature_var,
         d_scr_cxt,
         IntVector(d_dim, 1));

   d_diff_var.reset(new SideVariable<double>(d_dim, "diffusion",
         hier::IntVector::getOne(d_dim), 1));

   d_diff_id = variable_db->registerVariableAndContext(d_diff_var,
         d_cur_cxt,
         IntVector(d_dim, 0));

   /*
    * Set default values for preconditioner.
    */
   d_current_temperature_time = 0.;

   /*
    * Print solver data.
    */
   d_print_solver_info = false;

   /*
    * Counters.
    */
   d_number_rhs_eval = 0;
   d_number_precond_setup = 0;
   d_number_precond_solve = 0;

   /*
    * Initialize object with data read from given input/restart databases.
    */
   bool is_from_restart = RestartManager::getManager()->isFromRestart();
   if (is_from_restart) {
      getFromRestart();
   }
   getFromInput(input_db, is_from_restart);

   d_temperature_bc_helper =
      new solv::CartesianRobinBcHelper(d_dim,"BChelper");
   d_temperature_bc_coeffs =
      new solv::LocationIndexRobinBcCoefs(d_dim,"BCcoeffs",
         input_db->getDatabase( "BoundaryConditions" )); 

   d_temperature_bc_helper->setTargetDataId( d_temperature_scr_id );
   d_temperature_bc_helper->setCoefImplementation( d_temperature_bc_coeffs );

   /*
    * Boundary conditions for FAC solvers should be homogeneous
    * since solver computes corrections to current guess
    */
   d_temperature_bc_corr_coeffs =
      new solv::LocationIndexRobinBcCoefs(d_dim,"BCcorrcoeffs",
         input_db->getDatabase( "BoundaryConditions" ));
   for(int i=0;i<d_dim.getValue()*2;i++){
      double a,b,g;
      d_temperature_bc_corr_coeffs->getCoefficients(i,a,b,g);
      g=0.;
      d_temperature_bc_corr_coeffs->setRawCoefficients(i,a,b,g);
   }

   d_FAC_solver->setBcObject(d_temperature_bc_corr_coeffs);
}

CVODEModel::~CVODEModel()
{
   std::shared_ptr<SAMRAIVectorReal<double> > temperature_samvect =
      Sundials_SAMRAIVector::getSAMRAIVector(d_solution_vector);
   Sundials_SAMRAIVector::destroySundialsVector(d_solution_vector);

   temperature_samvect->freeVectorComponents();
   temperature_samvect.reset();
}

/*************************************************************************
 *
 * Methods inherited from mesh::StandardTagAndInitStrategy.
 *
 ************************************************************************/
void CVODEModel::initializeLevelData(
   const std::shared_ptr<PatchHierarchy>& hierarchy,
   const int level_number,
   const double time,
   const bool can_be_refined,
   const bool initial_time,
   const std::shared_ptr<PatchLevel>& old_level,
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

   // This method is empty because initialization is taken care of
   // by the setInitialConditions() method below.  If there is any
   // data that is not managed inside the SAMRAI CVODESolver class
   // but that must be set on the level, do it here.
}

void CVODEModel::resetHierarchyConfiguration(
   const std::shared_ptr<PatchHierarchy>& hierarchy,
   const int coarsest_level,
   const int finest_level)
{
   NULL_USE(hierarchy);
   NULL_USE(coarsest_level);
   NULL_USE(finest_level);

   // This method is empty because this example does not exercise the
   // situation when the grid changes, so it effectively is never called.
   // This is a subject for future work...
}

/*************************************************************************
 *
 * Cell tagging and patch level data initialization routines declared
 * in the GradientDetectorStrategy interface.  They are used to
 * construct the hierarchy initially.
 *
 *************************************************************************/
void CVODEModel::applyGradientDetector(
   const std::shared_ptr<PatchHierarchy>& hierarchy,
   const int level_number,
   const double time,
   const int tag_index,
   const bool initial_time,
   const bool uses_richardson_extrapolation_too)
{
   NULL_USE(time);
   NULL_USE(initial_time);
   NULL_USE(uses_richardson_extrapolation_too);

   std::shared_ptr<PatchLevel> level(
      hierarchy->getPatchLevel(level_number));

   for (PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {
      const std::shared_ptr<Patch>& patch = *p;

      std::shared_ptr<CellData<int> > tag_data(
         SAMRAI_SHARED_PTR_CAST<CellData<int>, PatchData>(
            patch->getPatchData(tag_index)));
      TBOX_ASSERT(tag_data);

      // dumb implementation that tags all cells.
      tag_data->fillAll(1);
   }
}

/*************************************************************************
 *
 * Methods inherited from RefinePatchStrategy.
 *
 ***********************************************************************/
void CVODEModel::setPhysicalBoundaryConditions(
   Patch& patch,
   const double time,
   const IntVector& ghost_width_to_fill)
{
   d_temperature_bc_helper->setPhysicalBoundaryConditions(
      patch,
      time,
      ghost_width_to_fill);

   //plog << "----Boundary Conditions "  << endl;
   //std::shared_ptr<CellData<double> > temperature_data(
   //   SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
   //      patch.getPatchData(d_temperature_scr_id)));
   //temperature_data->print(temperature_data->getGhostBox());
}

void CVODEModel::preprocessRefine(
   Patch& fine,
   const Patch& coarse,
   const Box& fine_box,
   const IntVector& ratio)
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_box);
   NULL_USE(ratio);
}

void CVODEModel::postprocessRefine(
   Patch& fine,
   const Patch& coarse,
   const Box& fine_box,
   const IntVector& ratio)
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_box);
   NULL_USE(ratio);
}

/*************************************************************************
 *
 * Methods inherited from CoarsenPatchStrategy.
 *
 ************************************************************************/
void CVODEModel::preprocessCoarsen(
   Patch& coarse,
   const Patch& fine,
   const Box& coarse_box,
   const IntVector& ratio)
{
   NULL_USE(coarse);
   NULL_USE(fine);
   NULL_USE(coarse_box);
   NULL_USE(ratio);
}

void CVODEModel::postprocessCoarsen(
   Patch& coarse,
   const Patch& fine,
   const Box& coarse_box,
   const IntVector& ratio)
{
   NULL_USE(coarse);
   NULL_USE(fine);
   NULL_USE(coarse_box);
   NULL_USE(ratio);
}

/*************************************************************************
 *
 * Methods inherited from CVODEAbstractFunction
 *
 ************************************************************************/
int CVODEModel::evaluateRHSFunction(
   double time,
   SundialsAbstractVector* y,
   SundialsAbstractVector* y_dot)
{
   /*
    * Convert Sundials vectors to SAMRAI vectors
    */
   std::shared_ptr<SAMRAIVectorReal<double> > y_samvect(
      Sundials_SAMRAIVector::getSAMRAIVector(y));
   std::shared_ptr<SAMRAIVectorReal<double> > y_dot_samvect(
      Sundials_SAMRAIVector::getSAMRAIVector(y_dot));

   std::shared_ptr<PatchHierarchy> hierarchy(y_samvect->getPatchHierarchy());

   if (d_print_solver_info) {
      pout << "\t\tEval RHS: "
           << "\n   \t\t\ttime = " << time
           << "\n   \t\t\ty_maxnorm = " << y_samvect->maxNorm()
           << endl;
   }

   /*
    * Allocate scratch space and fill ghost cells in the solution vector
    * 1) Create a refine algorithm
    * 2) Register with the algorithm the current & scratch space, along
    *    with a refine operator.
    * 3) Use the refine algorithm to construct a refine schedule
    * 4) Use the refine schedule to fill data on fine level.
    */
   std::shared_ptr<RefineAlgorithm> bdry_fill_alg(
      new RefineAlgorithm());
   std::shared_ptr<RefineOperator> refine_op(
      d_grid_geometry->lookupRefineOperator(d_temperature_var,
                                            "CONSERVATIVE_LINEAR_REFINE"));
   bdry_fill_alg->registerRefine(d_temperature_scr_id,  // dest
      y_samvect->getComponentDescriptorIndex(0), //src
      d_temperature_scr_id,                            // scratch
      refine_op);

   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      std::shared_ptr<PatchLevel> level(hierarchy->getPatchLevel(ln));
      if (!level->checkAllocated(d_temperature_scr_id)) {
         level->allocatePatchData(d_temperature_scr_id);
      }

      // Note:  a pointer to "this" tells the refine schedule to invoke
      // the setPhysicalBCs defined in this class.
      std::shared_ptr<RefineSchedule> bdry_fill_alg_schedule(
         bdry_fill_alg->createSchedule(level,
            ln - 1,
            hierarchy,
            this));

      bdry_fill_alg_schedule->fillData(time);
   }

   /*
    * Step through the levels and compute rhs
    */
   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      std::shared_ptr<PatchLevel> level(hierarchy->getPatchLevel(ln));

      for (PatchLevel::iterator ip(level->begin()); ip != level->end(); ++ip) {
         const std::shared_ptr<Patch>& patch = *ip;

         std::shared_ptr<CellData<double> > y(
            SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
               patch->getPatchData(d_temperature_scr_id)));
         std::shared_ptr<SideData<double> > diff(
            SAMRAI_SHARED_PTR_CAST<SideData<double>, PatchData>(
               patch->getPatchData(d_diff_id)));
         std::shared_ptr<CellData<double> > rhs(
            SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
               patch->getPatchData(
                  y_dot_samvect->getComponentDescriptorIndex(0))));
         TBOX_ASSERT(y);
         TBOX_ASSERT(diff);
         TBOX_ASSERT(rhs);

         const Index ifirst(patch->getBox().lower());
         const Index ilast(patch->getBox().upper());

         const std::shared_ptr<CartesianPatchGeometry> patch_geom(
            SAMRAI_SHARED_PTR_CAST<CartesianPatchGeometry, PatchGeometry>(
               patch->getPatchGeometry()));
         TBOX_ASSERT(patch_geom);
         const double* dx = patch_geom->getDx();

         IntVector ghost_cells(y->getGhostCellWidth());

         SAMRAI_F77_FUNC(comprhs3d, COMPRHS3D) (
            ifirst(0), ilast(0),
            ifirst(1), ilast(1),
            ifirst(2), ilast(2),
            ghost_cells(0), ghost_cells(1), ghost_cells(2),
            dx,
            y->getPointer(),
            diff->getPointer(0),
            diff->getPointer(1),
            diff->getPointer(2),
            rhs->getPointer());

      } // loop over patches
   } // loop over levels

   /*
    * Deallocate scratch space.
    */
   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      hierarchy->getPatchLevel(ln)->deallocatePatchData(d_temperature_scr_id);
   }

   /*
    * record current time and increment counter for number of RHS
    * evaluations.
    */
   d_current_temperature_time = time;
   ++d_number_rhs_eval;

   return 0;
}

/*****************************************************************
 *
 * Set up FAC preconditioner for Jacobian system.  Here we
 * use the FAC hierarchy solver in SAMRAI which automatically sets
 * up the composite grid system and uses hypre as a solver on each
 * level.
 *
 *****************************************************************/
int CVODEModel::CVSpgmrPrecondSet(
   double t,
   SundialsAbstractVector* y, // current value of variable vector,
                              // the predicted value of y(t)
   SundialsAbstractVector* fy,
   int jok,
   int* jcurPtr,
   double gamma,
   SundialsAbstractVector* vtemp1,
   SundialsAbstractVector* vtemp2,
   SundialsAbstractVector* vtemp3)
{
   NULL_USE(fy);
   NULL_USE(jok);
   NULL_USE(jcurPtr);
   NULL_USE(vtemp1);
   NULL_USE(vtemp2);
   NULL_USE(vtemp3);

   //plog<<"CVSpgmrPrecondSet..."<<endl;

   /*
    * Convert passed-in CVODE vectors into SAMRAI vectors
    */
   std::shared_ptr<SAMRAIVectorReal<double> > y_samvect(
      Sundials_SAMRAIVector::getSAMRAIVector(y));

   int y_indx = y_samvect->getComponentDescriptorIndex(0);

   /*
    * Construct refine algorithm to fill boundaries of solution vector
    */
   //RefineAlgorithm fill_temperature_vector_bounds;
   //std::shared_ptr<RefineOperator> refine_op(d_grid_geometry->
   //                                            lookupRefineOperator(d_temperature_var,
   //                                               "CONSERVATIVE_LINEAR_REFINE"));
   //fill_temperature_vector_bounds.registerRefine(d_temperature_scr_id,
   //   y_samvect->getComponentDescriptorIndex(0),
   //   d_temperature_scr_id,
   //   refine_op);

   /*
    * Construct coarsen algorithm to fill interiors on coarser levels
    * with solution on finer level.
    */
   CoarsenAlgorithm fill_temperature_interior_on_coarser(d_dim);
   std::shared_ptr<CoarsenOperator> coarsen_op(
      d_grid_geometry->lookupCoarsenOperator(d_temperature_var,
                                             "CONSERVATIVE_COARSEN"));

   fill_temperature_interior_on_coarser.registerCoarsen(y_indx,
      y_indx,
      coarsen_op);

   /*
    * Step through levels - largest to smallest
    */
   std::shared_ptr<PatchHierarchy> hierarchy(
      y_samvect->getPatchHierarchy());

   for (int amr_level = hierarchy->getFinestLevelNumber();
        amr_level >= 0;
        --amr_level) {
      std::shared_ptr<PatchLevel> level(
         hierarchy->getPatchLevel(amr_level));

      //std::shared_ptr<RefineSchedule> fill_temperature_vector_bounds_sched =
      //   fill_temperature_vector_bounds.createSchedule(level,
      //      amr_level - 1,
      //      hierarchy,
      //      this);

      //if (!level->checkAllocated(d_temperature_scr_id)) {
      //   level->allocatePatchData(d_temperature_scr_id);
      //}
      //fill_temperature_vector_bounds_sched->fillData(t);

      /*
       * Construct a coarsen schedule for all levels larger than coarsest,
       * and fill interiors of solution vector on coarser levels using fine
       * data.
       */
      if (amr_level > 0) {
         std::shared_ptr<PatchLevel> coarser_level(
            hierarchy->getPatchLevel(amr_level - 1));

         std::shared_ptr<CoarsenSchedule> fill_temperature_interior_on_coarser_sched(
            fill_temperature_interior_on_coarser.createSchedule(coarser_level,
               level));

         fill_temperature_interior_on_coarser_sched->coarsenData();
      }

      for (PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {
         const std::shared_ptr<Patch>& patch = *p;

         const Index ifirst(patch->getBox().lower());
         const Index ilast(patch->getBox().upper());

         std::shared_ptr<SideData<double> > diffusion(
            SAMRAI_SHARED_PTR_CAST<SideData<double>, PatchData>(
               patch->getPatchData(d_diff_id)));
         TBOX_ASSERT(diffusion);

         diffusion->fillAll(d_temperature_diffusion);

      } // patch loop

      //level->deallocatePatchData(d_temperature_scr_id);

   } // level loop

   d_FAC_solver->setCConstant(1.0 / gamma);
   d_FAC_solver->setDPatchDataId(d_diff_id);

   /*
    * increment counter for number of precond setup calls
    */
   ++d_number_precond_setup;

   /*
    * We return 0 or 1 here - 0 if it passes, 1 if it fails.  For now,
    * just be optimistic and return 0. Eventually we should add some
    * assertion handling above to set what this value should be.
    */
   return 0;
}

/*************************************************************************
 *
 * Apply preconditioner where right-hand-side is "r" and "z" is the
 * solution.   This routine assumes that the preconditioner setup call
 * has already been invoked.  Return 0 if preconditioner fails;
 * return 1 otherwise.
 *
 *************************************************************************/
int CVODEModel::CVSpgmrPrecondSolve(
   double t,
   SundialsAbstractVector* y,
   SundialsAbstractVector* fy,
   SundialsAbstractVector* r,
   SundialsAbstractVector* z,
   double gamma,
   double delta,
   int lr,
   SundialsAbstractVector* vtemp)
{
   NULL_USE(y);
   NULL_USE(fy);
   NULL_USE(vtemp);
   NULL_USE(delta);
   NULL_USE(lr);

   //plog<<"CVSpgmrPrecondSolve..."<<endl;

   /*
    * Convert passed-in CVODE vectors into SAMRAI vectors
    */
   std::shared_ptr<SAMRAIVectorReal<double> > r_samvect(
      Sundials_SAMRAIVector::getSAMRAIVector(r));
   std::shared_ptr<SAMRAIVectorReal<double> > z_samvect(
      Sundials_SAMRAIVector::getSAMRAIVector(z));

   std::shared_ptr<PatchHierarchy> hierarchy(
      r_samvect->getPatchHierarchy());

   int r_indx = r_samvect->getComponentDescriptorIndex(0);
   int z_indx = z_samvect->getComponentDescriptorIndex(0);
   /*
    * We need to supply to the FAC solver a "version" of the z vector
    * that contains ghost cells.  The operations below allocate
    * on the patches a scratch context of the solution vector z and
    * fill it with z vector data
    *
    * Set initial guess for z (if applicable) and copy z data into the
    * solution scratch context.
    */
   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      std::shared_ptr<PatchLevel> level(hierarchy->getPatchLevel(ln));

      if (!level->checkAllocated(d_temperature_scr_id)) {
         level->allocatePatchData(d_temperature_scr_id);
      }

      for (PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {

         const std::shared_ptr<Patch>& patch = *p;

         std::shared_ptr<CellData<double> > z_data(
            SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
               patch->getPatchData(z_indx)));
         TBOX_ASSERT(z_data);

         /*
          * Set initial guess for z here.
          */
         z_data->fillAll(0.);

         /*
          * Scale RHS by 1/gamma
          */
         PatchCellDataOpsReal<double> math_ops;
         std::shared_ptr<CellData<double> > r_data(
            SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
               patch->getPatchData(r_indx)));
         TBOX_ASSERT(r_data);
         math_ops.scale(r_data, 1.0 / gamma, r_data, r_data->getBox());

         /*
          * Copy interior data from z vector to temperature_scratch
          */
         std::shared_ptr<CellData<double> > z_scr_data(
            SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
               patch->getPatchData(d_temperature_scr_id)));
         TBOX_ASSERT(z_scr_data);
         z_scr_data->fillAll(0.);
      }

   }

   /******************************************************************
   *
   * Apply the FAC solver.  It solves the system Az=r with the
   * format "solveSystem(z, r)". A was constructed in the precondSetup()
   * method.
   *
   ******************************************************************/

   if (d_print_solver_info) {
      pout << "\t\tBefore FAC Solve (Az=r): "
           << "\n   \t\t\tz_l2norm = " << z_samvect->L2Norm()
           << "\n   \t\t\tz_maxnorm = " << z_samvect->maxNorm()
           << "\n   \t\t\tr_l2norm = " << r_samvect->L2Norm()
           << "\n   \t\t\tr_maxnorm = " << r_samvect->maxNorm()
           << endl;
   }
   /*
    * Set paramemters in the FAC solver.  It solves the system Az=r.
    * Here we supply the max norm of r in order to scale the
    * residual (i.e. residual = Az - r) to properly scale the convergence
    * error.
    */

   const int coarsest_solve_ln = 0;
   const int finest_solve_ln = hierarchy->getFinestLevelNumber();
   bool converge = d_FAC_solver->solveSystem(d_temperature_scr_id,
         r_indx,
         hierarchy,
         coarsest_solve_ln,
         finest_solve_ln);

   if (d_print_solver_info) {
      double avg_convergence, final_convergence;
      d_FAC_solver->getConvergenceFactors(avg_convergence, final_convergence);
      pout << "   \t\t\tFinal Residual Norm: "
           << d_FAC_solver->getResidualNorm() << endl;
      pout << "   \t\t\tFinal Convergence Error: "
           << final_convergence << endl;
      pout << "   \t\t\tFinal Convergence Rate: "
           << avg_convergence << endl;
   }

   /******************************************************************
   *
   * The FAC solver has computed a solution to z but it is stored
   * in the temperature_scratch data space.  Copy it from temperature_scratch back
   * into the z vector, including ghost values
   *
   ******************************************************************/
   math::HierarchyCellDataOpsReal<double> cell_ops( hierarchy );
   cell_ops.copyData( z_indx, d_temperature_scr_id, false);

   if (d_print_solver_info) {
      double avg_convergence, final_convergence;
      d_FAC_solver->getConvergenceFactors(avg_convergence, final_convergence);
      pout << "\t\tAfter FAC Solve (Az=r): "
           << "\n   \t\t\tz_l2norm = " << z_samvect->L2Norm()
           << "\n   \t\t\tz_maxnorm = " << z_samvect->maxNorm()
           << "\n   \t\t\tResidual Norm: " << d_FAC_solver->getResidualNorm()
           << "\n   \t\t\tConvergence Error: " << final_convergence
           << endl;
   }

   int ret_val = 0;
   if (converge != true) {
      ret_val = 1;
   }

   /*
    * Increment counter for number of precond solves
    */
   ++d_number_precond_solve;

   return ret_val;
}

/*************************************************************************
 *
 * Methods specific to CVODEModel class.
 *
 ************************************************************************/
void CVODEModel::setupSolutionVector(
   std::shared_ptr<PatchHierarchy> hierarchy)
{
   /* create SAMRAIVector */
   std::shared_ptr<SAMRAIVectorReal<double> > temperature_samvect(
      new SAMRAIVectorReal<double>(
         "solution",
         hierarchy,
         0,
         hierarchy->getFinestLevelNumber()));
   temperature_samvect->addComponent(d_temperature_var, d_temperature_cur_id);

   /* allocate memory for vectors. */
   temperature_samvect->allocateVectorData();

   /* create SundialsAbstractVector */
   d_solution_vector =
      Sundials_SAMRAIVector::createSundialsVector(temperature_samvect);

   /*
    * Allocate memory for preconditioner variables.
    */

   const int nlevels = hierarchy->getNumberOfLevels();

   for (int ln = 0; ln < nlevels; ++ln) {
      std::shared_ptr<PatchLevel> level(hierarchy->getPatchLevel(ln));
      TBOX_ASSERT(level);
      level->allocatePatchData(d_diff_id);
   }

}

SundialsAbstractVector *
CVODEModel::getSolutionVector(
   void)
{
   return d_solution_vector;
}

/*************************************************************************
 *
 * Set initial conditions for CVODE solver
 *
 *************************************************************************/
void CVODEModel::setInitialConditions(
   SundialsAbstractVector* temperature_init)
{
   std::shared_ptr<SAMRAIVectorReal<double> > temperature_init_samvect(
      Sundials_SAMRAIVector::getSAMRAIVector(temperature_init));

   std::shared_ptr<PatchHierarchy> hierarchy(
      temperature_init_samvect->getPatchHierarchy());

   for (int ln = 0; ln < hierarchy->getNumberOfLevels(); ++ln) {
      std::shared_ptr<PatchLevel> level(hierarchy->getPatchLevel(ln));

      for (int cn = 0; cn < temperature_init_samvect->getNumberOfComponents(); ++cn) {
         for (PatchLevel::iterator p(level->begin()); p != level->end(); ++p) {
            const std::shared_ptr<Patch>& patch = *p;

            std::shared_ptr<geom::CartesianPatchGeometry> pg(
               SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
               patch->getPatchGeometry()));
            TBOX_ASSERT(pg);

            const double* h = pg->getDx();

            /*
             * Set initial conditions for y
             */
            std::shared_ptr<CellData<double> > y_init(
               SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
                  temperature_init_samvect->getComponentPatchData(cn, *patch)));
            TBOX_ASSERT(y_init);

            const hier::Box patch_box = patch->getBox();
            pdat::CellIterator ic(pdat::CellGeometry::begin(patch_box));
            pdat::CellIterator icend(pdat::CellGeometry::end(patch_box));

            for ( ; ic != icend; ++ic) {

               hier::IntVector icell = *ic;
               const double yval=h[1]*(0.5+icell[1]);
               const double val=2.*yval+sin(2.*M_PI*yval)+1.;

               (*y_init)(*ic)=val;
            }

            /*
             * Set initial diffusion coeff values.
             */
            std::shared_ptr<SideData<double> > diffusion(
               SAMRAI_SHARED_PTR_CAST<SideData<double>, PatchData>(
                  patch->getPatchData(d_diff_id)));
            TBOX_ASSERT(diffusion);

            diffusion->fillAll(d_temperature_diffusion);
         }
      }
   }
}

/*************************************************************************
 *
 * Print program counters.  Currently, the array holds the
 * following entries:
 *    1) number of RHS evaluations
 *    2) number of precond setup calls
 *    3) number of precond solve calls
 *
 *************************************************************************/
void CVODEModel::printCounters(const double final_time)
{
   std::vector<int> counters(3);
   counters[0] = d_number_rhs_eval;
   counters[1] = d_number_precond_setup;
   counters[2] = d_number_precond_solve;

   tbox::plog << "\n\nEnd Timesteps - final time = " << final_time
              << "\n\tTotal number of RHS evaluations = " << counters[0]
              << "\n\tTotal number of precond setups = " << counters[1]
              << "\n\tTotal number of precond solves = " << counters[2]
              << endl;
}

/*
 *************************************************************************
 *
 * Get data from input database.
 *
 *************************************************************************
 */
void
CVODEModel::getFromInput(
   std::shared_ptr<Database> input_db,
   bool is_from_restart)
{
   NULL_USE(is_from_restart);

   std::shared_ptr<Database> temperature_db(
      input_db->getDatabase("Temperature"));

   d_temperature_diffusion = temperature_db->getDouble("diffusion_value");

   d_print_solver_info =
      input_db->getBoolWithDefault("print_solver_info", d_print_solver_info);
}

/*
 *************************************************************************
 *
 * Write data to  restart database.
 *
 *************************************************************************
 */
void CVODEModel::putToRestart(
   const std::shared_ptr<Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   restart_db->putDouble("d_temperature_diffusion", d_temperature_diffusion);
}

/*************************************************************************
 *
 * Read data from restart database.
 *
 *************************************************************************/
void CVODEModel::getFromRestart()
{
   std::shared_ptr<Database> root_db(
      RestartManager::getManager()->getRootDatabase());

   if (!root_db->isDatabase(d_object_name)) {
      TBOX_ERROR("Restart database corresponding to "
         << d_object_name << " not found in the restart file.");
   }
   std::shared_ptr<Database> db(root_db->getDatabase(d_object_name));

   d_temperature_diffusion = db->getDouble("d_temperature_diffusion");
}

/*************************************************************************
 *
 * Register VisIt data writer to write data to plot files that may
 * be postprocessed by the VisIt tool.
 *
 *************************************************************************/
void CVODEModel::registerVisItDataWriter(
   std::shared_ptr<appu::VisItDataWriter> viz_writer)
{
   TBOX_ASSERT(viz_writer);
   d_visit_writer = viz_writer;

   if (d_visit_writer) {
      d_visit_writer->
         registerPlotQuantity("temperature", "SCALAR",
         d_temperature_cur_id, 0);
   }
}

/*************************************************************************
 *
 * Prints class data - writes out info in class if assertion is thrown
 *
 *************************************************************************/
void CVODEModel::printClassData(
   ostream& os) const
{
   fflush(stdout);

   os << "ptr CVODEModel = " << (CVODEModel *)this << endl;

   os << "d_object_name = " << d_object_name << endl;

   os << "d_temperature_cur_id = " << d_temperature_cur_id << endl;
   os << "d_temperature_scr_id = " << d_temperature_scr_id << endl;

   os << "d_temperature_diffusion = " << d_temperature_diffusion << endl;
   os << endl;
}

void CVODEModel::setPrintSolverInfo(
   const bool info)
{
   d_print_solver_info = info;
}
