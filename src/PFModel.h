#ifndef PFiSM_PFModel_H
#define PFiSM_PFModel_H

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Weffc++"
#endif

#include "SAMRAI/SAMRAI_config.h"

// Header files for SAMRAI classes
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/Box.h"
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

// Header files for CVODE wrapper classes
#include "SAMRAI/solv/SundialsAbstractVector.h"

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#include <vector>
#include <iostream>

class CellPoissonFACSolver;

using namespace SAMRAI;

class PFModel : public mesh::StandardTagAndInitStrategy,
                public xfer::RefinePatchStrategy,
                public xfer::CoarsenPatchStrategy
{
 public:
   PFModel(const std::string& object_name, const tbox::Dimension& dim,
           bool evolve_temperature,
           std::shared_ptr<CellPoissonFACSolver> fac_solver_temperature,
           std::shared_ptr<CellPoissonFACSolver> fac_solver_phase,
           std::shared_ptr<tbox::Database> input_db,
           std::shared_ptr<geom::CartesianGridGeometry> grid_geom);

   virtual ~PFModel();

   ///
   /// Methods inherited from StandardTagAndInitStrategy.
   ///
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

   ///
   /// Methods inherited from RefinePatchStrategy.
   ///
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

   ///
   /// Methods inherited from CoarsenPatchStrategy.
   ///
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

   ///
   /// Methods specific to PFModel class.
   ///
   void printConvergenceFactors();

   void setPrintSolverInfo(const bool info) { d_print_solver_info = info; }

   void setupSolutionVector(std::shared_ptr<hier::PatchHierarchy> hierarchy);

   solv::SundialsAbstractVector* getSolutionVector()
   {
      return d_solution_vector;
   }

   ///
   /// Set initial conditions
   ///
   void setInitialConditions();

   ///
   /// evaluate solid fraction
   ///
   double computeSolidFraction(
       const std::shared_ptr<hier::PatchHierarchy>& hierarchy);

   ///
   /// Print program counters.
   ///
   void printCounters(const double);

   ///
   /// Writes state of PFModel object to the specified restart database.
   /// This routine is a concrete implementation of the function
   /// declared in the tbox::Serializable abstract base class.
   ///
   void putToRestart(const std::shared_ptr<tbox::Database>& restart_db) const;

   ///
   /// Register a VisIt data writer so this class will write
   /// plot files that may be postprocessed with VisIt
   ///
   void registerVisItDataWriter(
       std::shared_ptr<appu::VisItDataWriter> viz_writer);

   bool evolve_temperature() const { return d_evolve_temperature; }

 protected:
   ///
   /// Internal functions to evaluate RHS and preconditioner
   ///
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

   void fillScratch(const double time,
                    std::shared_ptr<hier::PatchHierarchy> hierarchy,
                    const int temperature_src_id, const int phase_src_id);

   void allocateScratchData(std::shared_ptr<hier::PatchHierarchy> hierarchy);
   void deallocateScratchData(std::shared_ptr<hier::PatchHierarchy> hierarchy);

   void initializeSolvers(
       const std::shared_ptr<hier::PatchHierarchy>& hierarchy);

 private:
   PFModel(const PFModel&);
   PFModel& operator=(const PFModel&);

   ///
   /// Read model parameters from input file
   ///
   void getFromInput(std::shared_ptr<tbox::Database> input_db,
                     bool is_from_restart);

   void getFromRestart();

   ///
   /// Object name used for error/warning reporting and as a label
   /// for restart database entries.
   ///
   std::string d_object_name;

   ///
   /// problem dimension
   ///
   const tbox::Dimension d_dim;

   solv::SundialsAbstractVector* d_solution_vector;

   ///
   /// Variables
   ///
   std::shared_ptr<pdat::CellVariable<double> > d_temperature_var;
   std::shared_ptr<pdat::CellVariable<double> > d_phase_var;
   std::shared_ptr<pdat::CellVariable<double> > d_cfield_phase_var;
   std::shared_ptr<pdat::CellVariable<double> > d_vol_var;

   ///
   /// Variable Contexts
   ///
   std::shared_ptr<hier::VariableContext> d_cur_cxt;

   ///
   /// scratch context
   ///
   std::shared_ptr<hier::VariableContext> d_scr_cxt;

   ///
   /// hier::Patch Data ids
   ///
   int d_temperature_cur_id;
   int d_temperature_scr_id;
   int d_phase_cur_id;
   int d_phase_scr_id;
   int d_cfield_phase_id;
   int d_vol_id;

   ///
   /// component indexes for SAMRAI vectors
   ///
   const int d_temperature_component;
   const int d_phase_component;

   const bool d_evolve_temperature;

   ///
   /// FAC solvers to be used in preconditioner
   ///
   std::shared_ptr<CellPoissonFACSolver> d_FAC_solver_temperature;
   std::shared_ptr<CellPoissonFACSolver> d_FAC_solver_phase;

   ///
   /// Flag to print FAC solver convergence information
   ///
   bool d_print_solver_info;

   ///
   /// Grid geometry
   ///
   std::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;

   std::shared_ptr<appu::VisItDataWriter> d_visit_writer;

   ///
   /// Temperature equation parameters
   ///
   double d_temperature_diffusion;
   double d_temperature_init;
   double d_Tmelting;
   double d_latent_heat;
   double d_cp;

   ///
   /// phase-field parameters
   ///
   double d_epsilon;
   double d_mobility;
   double d_well_height;
   double d_init_solid_fraction;

   ///
   /// Timers
   ///
   std::shared_ptr<tbox::Timer> t_factemperature_timer;
   std::shared_ptr<tbox::Timer> t_factempinit_timer;
   std::shared_ptr<tbox::Timer> t_facphase_timer;
   std::shared_ptr<tbox::Timer> t_facphaseinit_timer;

   ///
   /// Utilities to setup physical boundary conditions
   ///
   solv::CartesianRobinBcHelper* d_temperature_bc_helper;
   solv::LocationIndexRobinBcCoefs* d_temperature_bc_coeffs;

   solv::CartesianRobinBcHelper* d_phase_bc_helper;
   solv::LocationIndexRobinBcCoefs* d_phase_bc_coeffs;

   ///
   /// set C field in elliptic preconditioner for phase variable
   ///
   void setCforPhase(const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
                     const int phase_id, const double gamma);
};

#endif
