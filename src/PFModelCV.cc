#include "PFModelCV.h"

#include "samrai_internal/CellPoissonFACSolver.h"

#include "SAMRAI/hier/RefineOperator.h"

PFModelCV::PFModelCV(
    const std::string& object_name, const tbox::Dimension& dim,
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

   tbox::TimerManager* tman = tbox::TimerManager::getManager();
   t_rhs_timer = tman->getTimer("PFiSM::rhs");
   t_precondset_timer = tman->getTimer("PFiSM::precondset");
   t_precondsolve_timer = tman->getTimer("PFiSM::precondsolve");
}

PFModelCV::~PFModelCV() {}
/*************************************************************************
 * Methods inherited from CVODEAbstractFunction
 ************************************************************************/
int PFModelCV::evaluateRHSFunction(double time, solv::SundialsAbstractVector* y,
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

   evaluateRHSPhase(hierarchy, y_dot_phase_id);

   if (evolve_temperature()) {
      int y_dot_temperature_id =
          y_dot_samvect->getComponentDescriptorIndex(d_temperature_component);

      evaluateRHSTemperature(hierarchy, y_dot_temperature_id, y_dot_phase_id);
   }

   ++d_number_rhs_eval;

   t_rhs_timer->stop();

   return 0;
}

/*****************************************************************
 * Set up FAC preconditioner for Jacobian system.
 *****************************************************************/
int PFModelCV::CVSpgmrPrecondSet(
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

   tbox::plog << "CVSpgmrPrecondSet..." << std::endl;

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
int PFModelCV::CVSpgmrPrecondSolve(double t, solv::SundialsAbstractVector* y,
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

   // tbox::plog<<"CVSpgmrPrecondSolve..."<<std::endl;

   // Convert passed-in CVODE vectors into SAMRAI vectors
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

void PFModelCV::printCounters(const double final_time)
{
   tbox::plog << "\n\nEnd Timesteps - final time = " << final_time
              << "\n\tTotal number of RHS evaluations = " << d_number_rhs_eval
              << "\n\tTotal number of precond setups = "
              << d_number_precond_setup
              << "\n\tTotal number of precond solves = "
              << d_number_precond_solve << std::endl;
}
