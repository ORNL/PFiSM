/*************************************************************************
 *
 * Adapted from SAMRAI/source/test/sundials
 *
 ************************************************************************/

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Weffc++"

#include "SAMRAI/SAMRAI_config.h"

/*
 * Header files for SAMRAI classes
 */
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/xfer/CoarsenPatchStrategy.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/solv/CellPoissonFACSolver.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"
#include "SAMRAI/solv/CartesianRobinBcHelper.h"

/*
 * Header files for CVODE wrapper classes
 */
#include "SAMRAI/solv/SundialsAbstractVector.h"
#include "SAMRAI/solv/Sundials_SAMRAIVector.h"
#include "SAMRAI/solv/CVODEAbstractFunctions.h"

#pragma GCC diagnostic pop

#include <vector>
#include <iostream>

using namespace SAMRAI;

class PFModel:
   public mesh::StandardTagAndInitStrategy,
   public xfer::RefinePatchStrategy,
   public xfer::CoarsenPatchStrategy,
   public solv::CVODEAbstractFunctions
{
public:

   PFModel(
      const std::string& object_name,
      const tbox::Dimension& dim,
      std::shared_ptr<solv::CellPoissonFACSolver> fac_solver_temperature,
      std::shared_ptr<solv::CellPoissonFACSolver> fac_solver_phase,
      std::shared_ptr<tbox::Database> input_db,
      std::shared_ptr<geom::CartesianGridGeometry> grid_geom);

   ~PFModel();

/*************************************************************************
 *
 * Methods inherited from StandardTagAndInitStrategy.
 *
 ************************************************************************/

   /**
    * Initialize data on a new level after it is inserted into an AMR patch
    * hierarchy by the gridding algorithm.  The level number indicates
    * that of the new level.
    *
    * Generally, when data is set, it is interpolated from coarser levels
    * in the hierarchy.  If the old level pointer in the argument list is
    * non-null, then data is copied from the old level to the new level
    * on regions of intersection between those levels before interpolation
    * occurs.   In this case, the level number must match that of the old
    * level.  The specific operations that occur when initializing level
    * data are determined by the particular solution methods in use; i.e.,
    * in the subclass of this abstract base class.
    *
    * The boolean argument initial_time indicates whether the level is
    * being introduced for the first time (i.e., at initialization time),
    * or after some regrid process during the calculation beyond the initial
    * hierarchy construction.  This information is provided since the
    * initialization of the data may be different in each of those
    * circumstances.  The can_be_refined boolean argument indicates whether
    * the level is the finest allowable level in the hierarchy.
    */
   void initializeLevelData(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double time,
      const bool can_be_refined,
      const bool initial_time,
      const std::shared_ptr<hier::PatchLevel>& old_level =
         std::shared_ptr<hier::PatchLevel>(),
      const bool allocate_data = true);

   /**
    * After hierarchy levels have changed and data has been initialized on
    * the new levels, this routine can be used to reset any information
    * needed by the solution method that is particular to the hierarchy
    * configuration.  For example, the solution procedure may cache
    * communication schedules to amortize the cost of data movement on the
    * AMR patch hierarchy.  This function will be called by the gridding
    * algorithm after the initialization occurs so that the algorithm-specific
    * subclass can reset such things.  Also, if the solution method must
    * make the solution consistent across multiple levels after the hierarchy
    * is changed, this process may be invoked by this routine.  Of course the
    * details of these processes are determined by the particular solution
    * methods in use.
    *
    * The level number arguments indicate the coarsest and finest levels
    * in the current hierarchy configuration that have changed.  It should
    * be assumed that all intermediate levels have changed as well.
    */
   void resetHierarchyConfiguration(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int coarsest_level,
      const int finest_level);

   /**
    * Set tags to the specified tag value where refinement of the given
    * level should occur using the user-supplied gradient detector.  The
    * value "tag_index" is the index of the cell-centered integer tag
    * array on each patch in the hierarchy.  The boolean argument indicates
    * whether cells are being tagged on the level for the first time;
    * i.e., when the hierarchy is initially constructed.  If it is false,
    * it should be assumed that cells are being tagged at some later time
    * after the patch hierarchy was initially constructed.  This information
    * is provided since the application of the error estimator may be
    * different in each of those circumstances.
    */
   void applyGradientDetector(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double time,
      const int tag_index,
      const bool initial_time,
      const bool uses_richardson_extrapolation_too);

   void setPrintSolverInfo(const bool info){
      d_print_solver_info = info; 
  }

/*************************************************************************
 *
 * Methods inherited from RefinePatchStrategy.
 *
 ************************************************************************/

   /**
    * Set the data at patch boundaries corresponding to the physical domain
    * boundary.  The specific boundary conditions are determined by the user.
    */
   void setPhysicalBoundaryConditions(
      hier::Patch& patch,
      const double time,
      const hier::IntVector& ghost_width_to_fill);

   /**
    * Perform user-defined refining operations.  This member function
    * is called before the other refining operators.  The preprocess
    * function should refine data from the scratch components of the
    * coarse patch into the scratch components of the fine patch on the
    * specified fine box region.  This version of the preprocess function
    * operates on a a single box at a time.  The user must define this
    * routine in the subclass.
    */
   void preprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio);

   /**
    * Perform user-defined refining operations.  This member function
    * is called after the other refining operators.  The postprocess
    * function should refine data from the scratch components of the
    * coarse patch into the scratch components of the fine patch on the
    * specified fine box region.  This version of the postprocess function
    * operates on a a single box at a time.  The user must define this
    * routine in the subclass.
    */
   void postprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio);

   /**
    * Return maximum stencil width needed for user-defined
    * data interpolation operations.  Default is to return
    * zero, assuming no user-defined operations provided.
    */
   hier::IntVector getRefineOpStencilWidth(
      const tbox::Dimension& dim) const
   {
      return hier::IntVector(dim, 0);
   }

/*************************************************************************
 *
 * Methods inherited from CoarsenPatchStrategy.
 *
 ************************************************************************/

   /**
    * Perform user-defined coarsening operations.  This member function
    * is called before the other coarsening operators.  The preprocess
    * function should copy data from the source components of the fine
    * patch into the source components of the destination patch on the
    * specified coarse box region.
    */
   void preprocessCoarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio);

   /**
    * Perform user-defined coarsening operations.  This member function
    * is called after the other coarsening operators.  The postprocess
    * function should copy data from the source components of the fine
    * patch into the source components of the destination patch on the
    * specified coarse box region.
    */
   void postprocessCoarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio);

   hier::IntVector getCoarsenOpStencilWidth(
      const tbox::Dimension& dim) const
   {
      return hier::IntVector(dim, 0);
   }

   const tbox::Dimension& getDim() const
   {
      return d_dim;
   }

/*************************************************************************
 *
 * Methods inherited from CVODEAbstractFunctions
 *
 ************************************************************************/

   /**
    * User-supplied right-hand side function evaluation.
    *
    * The function arguments are:
    *
    * - \b t        (INPUT) {current value of the independent variable}
    * - \b y        (INPUT) {current value of dependent variable vector}
    * - \b y_dot   (OUTPUT){current value of the derivative of y}
    *
    * IMPORTANT: This function must not modify the vector y.
    */
   int evaluateRHSFunction(
      double time,
      solv::SundialsAbstractVector* y,
      solv::SundialsAbstractVector* y_dot);

   int CVSpgmrPrecondSet(
      double t,
      solv::SundialsAbstractVector* y,
      solv::SundialsAbstractVector* fy,
      int jok,
      int* jcurPtr,
      double gamma,
      solv::SundialsAbstractVector* vtemp1,
      solv::SundialsAbstractVector* vtemp2,
      solv::SundialsAbstractVector* vtemp3);

   int CVSpgmrPrecondSolve(
      double t,
      solv::SundialsAbstractVector* y,
      solv::SundialsAbstractVector* fy,
      solv::SundialsAbstractVector* r,
      solv::SundialsAbstractVector* z,
      double gamma,
      double delta,
      int lr,
      solv::SundialsAbstractVector* vtemp);

/*************************************************************************
 * Methods specific to PFModel class.
 ************************************************************************/

   void setupSolutionVector(std::shared_ptr<hier::PatchHierarchy> hierarchy);

   solv::SundialsAbstractVector* getSolutionVector()
   {
      return d_solution_vector;
   }

   /**
    * Set initial conditions for problem.
    */
   void setInitialConditions();

   /**
    * Print program counters.
    */
   void printCounters(const double);

   /**
    * Writes state of PFModel object to the specified restart database.
    * This routine is a concrete implementation of the function
    * declared in the tbox::Serializable abstract base class.
    */
   void putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const;

   /**
    * Register a VisIt data writer so this class will write
    * plot files that may be postprocessed with VisIt
    */
   void registerVisItDataWriter(
      std::shared_ptr<appu::VisItDataWriter> viz_writer);

   // Prints all class data members, if assertion is thrown.
   void printClassData(std::ostream& os) const;

private:
   PFModel(const PFModel&);
   PFModel& operator=(const PFModel&);

   /*
    * These private member functions read data from input and restart.
    * When beginning a run from a restart file, all data members are read
    * from the restart file.  If the boolean flag is true when reading
    * from input, some restart values may be overridden by those in the
    * input file.
    */
   void getFromInput(
      std::shared_ptr<tbox::Database> input_db,
      bool is_from_restart);

   void getFromRestart();

   /*
    * Object name used for error/warning reporting and as a label
    * for restart database entries.
    */
   std::string d_object_name;

   const tbox::Dimension d_dim;

   solv::SundialsAbstractVector* d_solution_vector;

   // Variables
   std::shared_ptr<pdat::CellVariable<double> > d_temperature_var;
   std::shared_ptr<pdat::CellVariable<double> > d_phase_var;
   std::shared_ptr<pdat::CellVariable<double> > d_cfield_phase_var;

   // Variable Contexts
   std::shared_ptr<hier::VariableContext> d_cur_cxt;
   std::shared_ptr<hier::VariableContext> d_scr_cxt;

   // hier::Patch Data ids
   int d_temperature_cur_id;
   int d_temperature_scr_id;
   int d_phase_cur_id;
   int d_phase_scr_id;
   int d_cfield_phase_id;

   // component indexes for SAMRAI vectors
   const int d_temperature_component;
   const int d_phase_component;

   // FAC solvers
   std::shared_ptr<solv::CellPoissonFACSolver> d_FAC_solver_temperature;
   std::shared_ptr<solv::CellPoissonFACSolver> d_FAC_solver_phase;

   bool d_level_solver_allocated;

   double d_current_time;

   // Print CVODE solver information
   bool d_print_solver_info;

   // Grid geometry
   std::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;

   std::shared_ptr<appu::VisItDataWriter> d_visit_writer;

   // Temperature equation parameters
   double d_temperature_diffusion;
   double d_temperature_init;
   double d_Tmelting;
   double d_latent_heat;
   double d_cp;

   // phase-field parameters
   double d_epsilon;
   double d_mobility;
   double d_well_height;
   double d_init_solid_fraction;

   /* Program counters
    *   1 - number of RHS evaluations
    *   2 - number of precond setups
    *   3 - number of precond solves
    */
   int d_number_rhs_eval;
   int d_number_precond_setup;
   int d_number_precond_solve;

   // Utilities to setup physical boundary conditions
   solv::CartesianRobinBcHelper* d_temperature_bc_helper;
   solv::LocationIndexRobinBcCoefs* d_temperature_bc_coeffs;

   solv::CartesianRobinBcHelper* d_phase_bc_helper;
   solv::LocationIndexRobinBcCoefs* d_phase_bc_coeffs;

   // Utilities to setup physical boundary conditions in FAC solver
   solv::LocationIndexRobinBcCoefs* d_temperature_bc_corr_coeffs;
   solv::LocationIndexRobinBcCoefs* d_phase_bc_corr_coeffs;

   void setCforPhase(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      std::shared_ptr<solv::SAMRAIVectorReal<double> > y_samvect,
      const double gamma);

};
