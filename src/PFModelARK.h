#ifndef PFiSM_PFModelARK_H
#define PFiSM_PFModelARK_H

#include "PFModel.h"

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Weffc++"
#endif

#include "SAMRAI/SAMRAI_config.h"

// Header files for SAMRAI classes
#include "SAMRAI/tbox/Database.h"

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

class PFModelARK : public PFModel, public ARKODEAbstractFunctions
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
    * Methods inherited from ARKODEAbstractFunctions
    ************************************************************************/
   int evaluateRHSFunction(double time, solv::SundialsAbstractVector* y,
                           solv::SundialsAbstractVector* y_dot);

   int evaluateRHSFunctionImp(double time, solv::SundialsAbstractVector* y,
                              solv::SundialsAbstractVector* y_dot);

   int evaluateRHSFunctionExp(double time, solv::SundialsAbstractVector* y,
                              solv::SundialsAbstractVector* y_dot);

   int ARKSpgmrPrecondSet(double t, solv::SundialsAbstractVector* y,
                          solv::SundialsAbstractVector* fy, int jok,
                          int* jcurPtr, double gamma);

   int ARKSpgmrPrecondSolve(double t, solv::SundialsAbstractVector* y,
                            solv::SundialsAbstractVector* fy,
                            solv::SundialsAbstractVector* r,
                            solv::SundialsAbstractVector* z, double gamma,
                            double delta, int lr);

   ///
   /// Print number of RHS evaluations, precond solves, ...
   ///
   void printCounters(const double);

 private:
   PFModelARK(const PFModelARK&);
   PFModelARK& operator=(const PFModelARK&);

   bool d_have_solver_temperature;

   // component indexes for SAMRAI vectors
   const int d_temperature_component;
   const int d_phase_component;

   // Moving frame
   double d_frame_velocity;

   // Program counters
   int d_number_rhs_eval;
   int d_number_rhs_exp_eval;
   int d_number_rhs_imp_eval;
   int d_number_precond_setup;
   int d_number_precond_solve;

   // Timers
   std::shared_ptr<tbox::Timer> t_rhs_timer;
   std::shared_ptr<tbox::Timer> t_precondset_timer;
   std::shared_ptr<tbox::Timer> t_precondsolve_timer;
};

#endif
