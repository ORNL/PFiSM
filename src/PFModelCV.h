#ifndef PFiSM_PFModelCV_H
#define PFiSM_PFModeliCV_H

#include "PFModel.h"

// Header files for CVODE wrapper classes
#include "SAMRAI/solv/SundialsAbstractVector.h"
#include "SAMRAI/solv/Sundials_SAMRAIVector.h"
#include "SAMRAI/solv/CVODEAbstractFunctions.h"

class CellPoissonFACSolver;

using namespace SAMRAI;


class PFModelCV : public PFModel, public solv::CVODEAbstractFunctions
{
 public:
   PFModelCV(const std::string& object_name, const tbox::Dimension& dim,
             bool evolve_temperature,
             std::shared_ptr<CellPoissonFACSolver> fac_solver_temperature,
             std::shared_ptr<CellPoissonFACSolver> fac_solver_phase,
             std::shared_ptr<tbox::Database> input_db,
             std::shared_ptr<geom::CartesianGridGeometry> grid_geom);

   virtual ~PFModelCV();

   ///
   /// Methods inherited from CVODEAbstractFunctions
   ///
   int evaluateRHSFunction(double time, solv::SundialsAbstractVector* y,
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

   ///
   /// Print number of RHS evaluations, precond solves, ...
   ///
   void printCounters(const double final_time);

 private:
   PFModelCV(const PFModelCV&);
   PFModelCV& operator=(const PFModelCV&);

   bool d_have_solver_temperature;

   ///
   /// component indexes for SAMRAI vectors
   ///
   const int d_temperature_component;
   const int d_phase_component;

   // Moving frame
   double d_frame_velocity;

   ///
   /// Program counters
   ///
   int d_number_rhs_eval;
   int d_number_precond_setup;
   int d_number_precond_solve;

   ///
   /// Timers
   ///
   std::shared_ptr<tbox::Timer> t_rhs_timer;
   std::shared_ptr<tbox::Timer> t_precondset_timer;
   std::shared_ptr<tbox::Timer> t_precondsolve_timer;
};

#endif
