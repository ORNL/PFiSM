#include "PFModelARK.h"
#include "tools.h"
#include "FortranInterface.h"

#include "samrai_internal/CellPoissonFACSolver.h"

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

#include "SAMRAI/solv/SAMRAIVectorReal.h"

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

using namespace std;

PFModelARK::PFModelARK(
    const string& object_name, const tbox::Dimension& dim,
    bool evolve_temperature,
    std::shared_ptr<CellPoissonFACSolver> fac_solver_temperature,
    std::shared_ptr<CellPoissonFACSolver> fac_solver_phase,
    std::shared_ptr<tbox::Database> input_db,
    std::shared_ptr<geom::CartesianGridGeometry> grid_geom)
    : PFModel(object_name, dim, evolve_temperature, fac_solver_temperature,
              fac_solver_phase, input_db, grid_geom),
      d_temperature_component(0),
      d_phase_component(1),
      d_number_rhs_eval(0),
      d_number_precond_setup(0),
      d_number_precond_solve(0)
{
   d_have_solver_temperature = fac_solver_temperature ? true : false;

   d_frame_velocity = input_db->getDouble("frame_velocity");

   tbox::TimerManager* tman = tbox::TimerManager::getManager();
   t_rhs_timer = tman->getTimer("PFiSM::rhs");
   t_precondset_timer = tman->getTimer("PFiSM::precondset");
   t_precondsolve_timer = tman->getTimer("PFiSM::precondsolve");
}

PFModelARK::~PFModelARK() {}

/*************************************************************************
 * Methods inherited from ARKODEAbstractFunction for fully explicit
 * or implicit case (complete r.h.s.)
 ************************************************************************/
int PFModelARK::evaluateRHSFunction(double time,
                                    solv::SundialsAbstractVector* y,
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

   // fill scratch with data in y_samvect
   int y_phase_id = y_samvect->getComponentDescriptorIndex(d_phase_component);
   int y_temperature_id =
       y_samvect->getComponentDescriptorIndex(d_temperature_component);
   fillScratch(time, hierarchy, y_temperature_id, y_phase_id);

   // now actually compute rhs
   int y_dot_phase_id =
       y_dot_samvect->getComponentDescriptorIndex(d_phase_component);

   evaluateRHSPhaseWithVelocity(hierarchy, y_dot_phase_id, d_frame_velocity);

   if (evolve_temperature()) {
      int y_dot_temperature_id =
          y_dot_samvect->getComponentDescriptorIndex(d_temperature_component);

      evaluateRHSTemperatureWithVelocity(hierarchy, y_dot_temperature_id,
                                         y_dot_phase_id, d_frame_velocity);
   }

   ++d_number_rhs_eval;

   t_rhs_timer->stop();

   return 0;
}

/*************************************************************************
 * Implicit:  Methods inherited from ARKODEAbstractFunction
 *************************************************************************/
int PFModelARK::evaluateRHSFunctionImp(double time,
                                       solv::SundialsAbstractVector* y,
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

   // fill scratch with data in y_samvect
   int y_phase_id = y_samvect->getComponentDescriptorIndex(d_phase_component);
   int y_temperature_id =
       y_samvect->getComponentDescriptorIndex(d_temperature_component);
   fillScratch(time, hierarchy, y_temperature_id, y_phase_id);

   int y_dot_phase_id =
       y_dot_samvect->getComponentDescriptorIndex(d_phase_component);

   if (evolve_temperature()) {
      // compute rhs for phase to be used in temperature equation
      evaluateRHSPhaseWithVelocity(hierarchy, y_dot_phase_id, d_frame_velocity);

      int y_dot_temperature_id =
          y_dot_samvect->getComponentDescriptorIndex(d_temperature_component);

      evaluateRHSTemperatureDiffusion(hierarchy, y_dot_temperature_id);

      evaluateRHSTemperature(hierarchy, y_dot_temperature_id, y_dot_phase_id);
   }

   // recompute rhs for phase without velocity term
   evaluateRHSPhaseDiffusion(hierarchy, y_dot_phase_id);
   evaluateRHSPhase(hierarchy, y_dot_phase_id);

   ++d_number_rhs_imp_eval;

   t_rhs_timer->stop();

   return 0;
}

/*************************************************************************
 * Explicit: Methods inherited from ARKODEAbstractFunction
 *************************************************************************/
int PFModelARK::evaluateRHSFunctionExp(double time,
                                       solv::SundialsAbstractVector* y,
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

   // fill scratch with data in y_samvect
   int y_phase_id = y_samvect->getComponentDescriptorIndex(d_phase_component);
   int y_temperature_id =
       y_samvect->getComponentDescriptorIndex(d_temperature_component);
   fillScratch(time, hierarchy, y_temperature_id, y_phase_id);

   int y_dot_phase_id =
       y_dot_samvect->getComponentDescriptorIndex(d_phase_component);
   int y_dot_temperature_id =
       y_dot_samvect->getComponentDescriptorIndex(d_temperature_component);

   evaluateRHSmovingFrame(hierarchy, y_dot_temperature_id, y_dot_phase_id,
                          d_frame_velocity);

   ++d_number_rhs_exp_eval;

   t_rhs_timer->stop();

   return 0;
}

/*****************************************************************
 * Set up FAC preconditioner for Jacobian system.
 *****************************************************************/
int PFModelARK::ARKSpgmrPrecondSet(
    double t,
    solv::SundialsAbstractVector* y,  // current value of variable vector,
                                      // the predicted value of y(t)
    solv::SundialsAbstractVector* fy, booleantype jok, int* jcurPtr,
    double gamma)
{
   NULL_USE(t);
   NULL_USE(fy);
   NULL_USE(jok);
   NULL_USE(jcurPtr);

   t_precondset_timer->start();

   // tbox::plog << "ARKSpgmrPrecondSet..." << endl;

   std::shared_ptr<solv::SAMRAIVectorReal<double> > y_samvect(
       solv::Sundials_SAMRAIVector::getSAMRAIVector(y));

   std::shared_ptr<hier::PatchHierarchy> hierarchy(
       y_samvect->getPatchHierarchy());

   // preconditioner for phase depends on values of phase field
   int y_indx = y_samvect->getComponentDescriptorIndex(d_phase_component);
   PrecondSetPhase(hierarchy, y_indx, gamma);

   if (d_have_solver_temperature) {
      // preconditioner for temperature does not depend on temperature field
      PrecondSetTemperature(gamma);
   }

   initializeSolvers(hierarchy);
   ++d_number_precond_setup;

   t_precondset_timer->stop();

   // assume success and return 0
   return 0;
}

/*************************************************************************
 * Apply preconditioner.
 * r:  right-hand-side
 * z: solution
 * Assumes preconditioner has been setup already.
 * Return 0 if preconditioner fails; 1 otherwise.
 *************************************************************************/
int PFModelARK::ARKSpgmrPrecondSolve(double t, solv::SundialsAbstractVector* y,
                                     solv::SundialsAbstractVector* fy,
                                     solv::SundialsAbstractVector* r,
                                     solv::SundialsAbstractVector* z,
                                     double gamma, double delta, int lr)
{
   NULL_USE(t);
   NULL_USE(y);
   NULL_USE(fy);
   NULL_USE(delta);
   NULL_USE(lr);

   t_precondsolve_timer->start();

   // tbox::plog<<"ARKSpgmrPrecondSolve..."<<std::endl;

   // Convert passed-in ARKODE vectors into SAMRAI vectors
   std::shared_ptr<solv::SAMRAIVectorReal<double> > r_samvect(
       solv::Sundials_SAMRAIVector::getSAMRAIVector(r));
   std::shared_ptr<solv::SAMRAIVectorReal<double> > z_samvect(
       solv::Sundials_SAMRAIVector::getSAMRAIVector(z));

   std::shared_ptr<hier::PatchHierarchy> hierarchy(
       r_samvect->getPatchHierarchy());

   bool converge0 = true;
   if (d_have_solver_temperature) {
      int r0_indx =
          r_samvect->getComponentDescriptorIndex(d_temperature_component);
      int z0_indx =
          z_samvect->getComponentDescriptorIndex(d_temperature_component);

      converge0 = PrecondSolveTemperature(hierarchy, r0_indx, z0_indx, gamma);
   }

   int r1_indx = r_samvect->getComponentDescriptorIndex(d_phase_component);
   int z1_indx = z_samvect->getComponentDescriptorIndex(d_phase_component);

   bool converge1 = PrecondSolvePhase(hierarchy, r1_indx, z1_indx, gamma);

   bool converge = (converge0 && converge1);

   printConvergenceFactors();

   ++d_number_precond_solve;

   int ret_val = (converge == true) ? 0 : 1;

   t_precondsolve_timer->stop();

   return ret_val;
}

void PFModelARK::printCounters(const double final_time)
{
   tbox::plog
       << "\n\nEnd Timesteps - final time = " << final_time
       << "\n\tTotal number of RHS evaluations = " << d_number_rhs_eval
       << "\n\tTotal number of explicit RHS    = " << d_number_rhs_exp_eval
       << "\n\tTotal number of implicit RHS    = " << d_number_rhs_imp_eval
       << "\n\tTotal number of precond setups = " << d_number_precond_setup
       << "\n\tTotal number of precond solves = " << d_number_precond_solve
       << endl;
}
