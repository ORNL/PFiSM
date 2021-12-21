#ifndef PFiSM_PFModelARK_H
#define PFiSM_PFModelARK_H

#include "Model.h"

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Weffc++"
#endif

#include "SAMRAI/SAMRAI_config.h"

// Header files for SAMRAI classes
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
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"
#include "SAMRAI/solv/CartesianRobinBcHelper.h"

// Header files for ARKODE wrapper classes
#include "SAMRAI/solv/SundialsAbstractVector.h"
#include "SAMRAI/solv/Sundials_SAMRAIVector.h"
#include "ARKODEAbstractFunctions.h"

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#include <vector>
#include <iostream>

class CellPoissonFACSolver;

using namespace SAMRAI;

class PFModelARK : public mesh::StandardTagAndInitStrategy,
                   public xfer::RefinePatchStrategy,
                   public xfer::CoarsenPatchStrategy,
                   public ARKODEAbstractFunctions,
                   public Model
{
 public:
   PFModelARK(const std::string& object_name, const tbox::Dimension& dim,
              bool evolve_temperature,
              std::shared_ptr<CellPoissonFACSolver> fac_solver_temperature,
              std::shared_ptr<CellPoissonFACSolver> fac_solver_phase,
              std::shared_ptr<tbox::Database> input_db,
              std::shared_ptr<geom::CartesianGridGeometry> grid_geom);

   virtual ~PFModelARK();

   /*************************************************************************
    * Methods inherited from StandardTagAndInitStrategy.
    ************************************************************************/
   void initializeLevelData(
       const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
       const int level_number, const double time, const bool can_be_refined,
       const bool initial_time,
       const std::shared_ptr<hier::PatchLevel>& old_level =
           std::shared_ptr<hier::PatchLevel>(),
       const bool allocate_data = true);

   void resetHierarchyConfiguration(
       const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
       const int coarsest_level, const int finest_level);

   void applyGradientDetector(
       const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
       const int level_number, const double time, const int tag_index,
       const bool initial_time, const bool uses_richardson_extrapolation_too);

   /*************************************************************************
    * Methods inherited from RefinePatchStrategy.
    ************************************************************************/
   void setPhysicalBoundaryConditions(
       hier::Patch& patch, const double time,
       const hier::IntVector& ghost_width_to_fill);

   void preprocessRefine(hier::Patch& fine, const hier::Patch& coarse,
                         const hier::Box& fine_box,
                         const hier::IntVector& ratio);

   void postprocessRefine(hier::Patch& fine, const hier::Patch& coarse,
                          const hier::Box& fine_box,
                          const hier::IntVector& ratio);

   hier::IntVector getRefineOpStencilWidth(const tbox::Dimension& dim) const
   {
      return hier::IntVector(dim, 0);
   }

   /*************************************************************************
    * Methods inherited from CoarsenPatchStrategy.
    ************************************************************************/
   void preprocessCoarsen(hier::Patch& coarse, const hier::Patch& fine,
                          const hier::Box& coarse_box,
                          const hier::IntVector& ratio);

   void postprocessCoarsen(hier::Patch& coarse, const hier::Patch& fine,
                           const hier::Box& coarse_box,
                           const hier::IntVector& ratio);

   hier::IntVector getCoarsenOpStencilWidth(const tbox::Dimension& dim) const
   {
      return hier::IntVector(dim, 0);
   }

   /*************************************************************************
    * Methods inherited from ARKODEAbstractFunctions
    ************************************************************************/
   int evaluateRHSFunction(double time, solv::SundialsAbstractVector* y,
                           solv::SundialsAbstractVector* y_dot);

   int evaluateRHSFunctionImp(double time, solv::SundialsAbstractVector* y,
                              solv::SundialsAbstractVector* y_dot);

   int evaluateRHSFunctionExp(double time, solv::SundialsAbstractVector* y,
                              solv::SundialsAbstractVector* y_dot);

   int evaluateJTimesRHSFunction(double t, solv::SundialsAbstractVector* y,
                                 solv::SundialsAbstractVector* y_dot)
   {
      return evaluateRHSFunction(t, y, y_dot);
   }

   int CVSpgmrPrecondSet(double t, solv::SundialsAbstractVector* y,
                         solv::SundialsAbstractVector* fy, int jok,
                         int* jcurPtr, double gamma);

   int CVSpgmrPrecondSolve(double t, solv::SundialsAbstractVector* y,
                           solv::SundialsAbstractVector* fy,
                           solv::SundialsAbstractVector* r,
                           solv::SundialsAbstractVector* z, double gamma,
                           double delta, int lr);

   int applyProjection(double time, solv::SundialsAbstractVector* y,
                       solv::SundialsAbstractVector* corr, double epsProj,
                       solv::SundialsAbstractVector* err)
   {
      (void)time;

      // Zero all components of the correction
      corr->setToScalar(0.);
      return 0;
   }

   /*************************************************************************
    * Methods specific to PFModelARK class.
    ************************************************************************/
   void setPrintSolverInfo(const bool info) { d_print_solver_info = info; }

   void setupSolutionVector(std::shared_ptr<hier::PatchHierarchy> hierarchy);

   solv::SundialsAbstractVector* getSolutionVector()
   {
      return d_solution_vector;
   }

   /**
    * Set initial conditions
    */
   void setInitialConditions();

   double computeSolidFraction(
       const std::shared_ptr<hier::PatchHierarchy>& hierarchy);

   /**
    * Print program counters.
    */
   void printCounters(const double);

   /**
    * Writes state of PFModelARK object to the specified restart database.
    * This routine is a concrete implementation of the function
    * declared in the tbox::Serializable abstract base class.
    */
   void putToRestart(const std::shared_ptr<tbox::Database>& restart_db) const;

   /**
    * Register a VisIt data writer so this class will write
    * plot files that may be postprocessed with VisIt
    */
   void registerVisItDataWriter(
       std::shared_ptr<appu::VisItDataWriter> viz_writer);

   // Prints all class data members, if assertion is thrown.
   void printClassData(std::ostream& os) const;

 private:
   PFModelARK(const PFModelARK&);
   PFModelARK& operator=(const PFModelARK&);

   void getFromInput(std::shared_ptr<tbox::Database> input_db,
                     bool is_from_restart);

   void getFromRestart();

   // Object name used for error/warning reporting and as a label
   // for restart database entries.
   std::string d_object_name;

   const tbox::Dimension d_dim;

   solv::SundialsAbstractVector* d_solution_vector;

   // Variables
   std::shared_ptr<pdat::CellVariable<double> > d_temperature_var;
   std::shared_ptr<pdat::CellVariable<double> > d_phase_var;
   std::shared_ptr<pdat::CellVariable<double> > d_cfield_phase_var;
   std::shared_ptr<pdat::CellVariable<double> > d_vol_var;

   // Variable Contexts
   std::shared_ptr<hier::VariableContext> d_cur_cxt;
   std::shared_ptr<hier::VariableContext> d_scr_cxt;

   // hier::Patch Data ids
   int d_temperature_cur_id;
   int d_temperature_scr_id;
   int d_phase_cur_id;
   int d_phase_scr_id;
   int d_cfield_phase_id;
   int d_vol_id;

   // component indexes for SAMRAI vectors
   const int d_temperature_component;
   const int d_phase_component;

   const bool d_evolve_temperature;

   // FAC solvers
   std::shared_ptr<CellPoissonFACSolver> d_FAC_solver_temperature;
   std::shared_ptr<CellPoissonFACSolver> d_FAC_solver_phase;

   bool d_level_solver_allocated;

   double d_current_time;

   // Print ARKODE solver information
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

   // Moving frame
   double d_frame_velocity;

   // Program counters
   int d_number_rhs_eval;
   int d_number_precond_setup;
   int d_number_precond_solve;

   // Timers
   std::shared_ptr<tbox::Timer> t_rhs_timer;
   std::shared_ptr<tbox::Timer> t_precondset_timer;
   std::shared_ptr<tbox::Timer> t_precondsolve_timer;
   std::shared_ptr<tbox::Timer> t_factemperature_timer;
   std::shared_ptr<tbox::Timer> t_factempinit_timer;
   std::shared_ptr<tbox::Timer> t_facphase_timer;
   std::shared_ptr<tbox::Timer> t_facphaseinit_timer;

   // Utilities to setup physical boundary conditions
   solv::CartesianRobinBcHelper* d_temperature_bc_helper;
   solv::LocationIndexRobinBcCoefs* d_temperature_bc_coeffs;

   solv::CartesianRobinBcHelper* d_phase_bc_helper;
   solv::LocationIndexRobinBcCoefs* d_phase_bc_coeffs;

   // Utilities to setup physical boundary conditions in FAC solver
   solv::LocationIndexRobinBcCoefs* d_temperature_bc_corr_coeffs;
   solv::LocationIndexRobinBcCoefs* d_phase_bc_corr_coeffs;

   void setCforPhase(const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
                     const int phase_id, const double gamma);

   void initializeSolvers(
       const std::shared_ptr<hier::PatchHierarchy>& hierarchy);

   void evaluateRHSPhase(std::shared_ptr<hier::PatchHierarchy> hierarchy,
                         const int y_dot_phase_id);
   void evaluateRHSTemperature(std::shared_ptr<hier::PatchHierarchy> hierarchy,
                               const int y_dot_temperature_id,
                               const int y_dot_phase_id);


   bool PrecondSolveTemperature(std::shared_ptr<hier::PatchHierarchy> hierarchy,
                                int r0_indx, int z0_indx, double gamma);
   bool PrecondSolvePhase(std::shared_ptr<hier::PatchHierarchy> hierarchy,
                          int r1_indx, int z1_indx, double gamma);

   void PrecondSetTemperature(const double gamma);
   void PrecondSetPhase(std::shared_ptr<hier::PatchHierarchy> hierarchy,
                        int y_indx, const double gamma);
};

#endif
