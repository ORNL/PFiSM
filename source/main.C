#include "_hypre_utilities.h"

#include "SAMRAI/SAMRAI_config.h"

#include "PFModel.h"
#include "PfmFACSolver.h"

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Weffc++"
#endif

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"

#include "SAMRAI/tbox/BalancedDepthFirstTree.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"

#include "SAMRAI/solv/SundialsAbstractVector.h"
#include "SAMRAI/solv/CVODESolver.h"

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

// CVODE includes
#ifndef included_cvspgmr_h
#define included_cvspgmr_h
#include "cvode/cvode_spgmr.h"
#endif

#ifndef _MSC_VER
#include <unistd.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>

using namespace SAMRAI;

int main(
   int argc,
   char* argv[])
{
   // Initialize tbox::MPI and SAMRAI
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   HYPRE_Init(argc, argv);

   {
      if (argc != 2) {
         tbox::pout << "USAGE:  " << argv[0] << " <input filename> " << std::endl;
         exit(-1);
      }

      std::string input_filename = argv[1];

      std::string run_name =
         input_filename.substr( 0, input_filename.rfind( "." ) );

      tbox::pout<<"Run name: "<<run_name<<std::endl;
      std::string log_file_name = run_name + ".log";
      tbox::PIO::logOnlyNodeZero( log_file_name );

      // Parse all data in input file
      tbox::pout<<"Parse input data..."<<std::endl;
      std::shared_ptr<tbox::InputDatabase> input_db(
         new tbox::InputDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(
         input_filename, input_db);

      // Retreive "Main" section of input db.
      std::shared_ptr<tbox::Database> main_db(input_db->getDatabase("Main"));

      const tbox::Dimension dim(static_cast<unsigned short>(3));

      int max_order = main_db->getInteger("max_order");
      int max_steps = main_db->getInteger("max_steps");
      double init_time = main_db->getDouble("init_time");
      int init_cycle = main_db->getInteger("init_cycle");
      double print_interval = main_db->getDouble("print_interval");
      bool visit_flag = main_db->getBoolWithDefault("visit_output", false);

      double relative_tolerance = main_db->getDouble("relative_tolerance");
      double absolute_tolerance = main_db->getDouble("absolute_tolerance");
      bool uses_preconditioning =
         main_db->getBoolWithDefault("uses_preconditioning", false);
      bool solution_logging =
         main_db->getBoolWithDefault("solution_logging", false);

      tbox::pout<<"Build Geometry, Hierarchy,..."<<std::endl;

      // Cartesian mesh geometry management
      std::shared_ptr<geom::CartesianGridGeometry> geometry(
         new geom::CartesianGridGeometry(
            dim,
            "Geometry",
            input_db->getDatabase("Geometry")));

      // maintains the patch levels that define the AMR hierarchy
      std::shared_ptr<hier::PatchHierarchy> hierarchy(
         new hier::PatchHierarchy(
            "Hierarchy",
            geometry,
            input_db->getDatabase("PatchHierarchy")));

      // create two FAC solvers for the two blocks of the preconditioner
      PfmFACSolver phase_fac_solver("PhaseFACsolver", dim, input_db);
      std::shared_ptr<CellPoissonFACSolver> fac_solver_temperature(
         phase_fac_solver.getCellPoissonFACSolver());

      PfmFACSolver temperature_fac_solver("TemperatureFACsolver",dim,input_db);
      std::shared_ptr<CellPoissonFACSolver> fac_solver_phase(
         temperature_fac_solver.getCellPoissonFACSolver());

      // construct main object
      std::shared_ptr<PFModel> pf_model(
         new PFModel(
            "PFModel",
            dim,
            fac_solver_temperature,
            fac_solver_phase,
            input_db->getDatabase("PFModel"),
            geometry));

      // defines an implementation for level initialization and cell
      // tagging routines needed by the GriddingAlgorithm class
      std::shared_ptr<mesh::StandardTagAndInitialize> error_est(
         new mesh::StandardTagAndInitialize(
            "StandardTagAndInitialize",
            pf_model.get(),
            input_db->getDatabase("StandardTagAndInitialize")));

      std::shared_ptr<mesh::BergerRigoutsos> box_generator(
         new mesh::BergerRigoutsos(dim,
            input_db->getDatabase("BergerRigoutsos")));

      std::shared_ptr<mesh::TreeLoadBalancer> load_balancer(
         new mesh::TreeLoadBalancer(
            dim,
            "LoadBalancer",
            input_db->getDatabase("LoadBalancer")));
      load_balancer->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());

      // Create gridding algorithm objects that will handle construction of
      // of the patch levels in the hierarchy.
      std::shared_ptr<mesh::GriddingAlgorithm> gridding_algorithm(
         new mesh::GriddingAlgorithm(
            hierarchy,
            "GriddingAlgorithm",
            input_db->getDatabase("GriddingAlgorithm"),
            error_est,
            box_generator,
            load_balancer));

      input_db->printClassData(tbox::plog);

      // create coarse level
      gridding_algorithm->makeCoarsestLevel(init_time);

      // creates finer levela
      std::vector<int> tag_buffer_array(hierarchy->getMaxNumberOfLevels());
      for (int il = 0; il < hierarchy->getMaxNumberOfLevels(); ++il) {
         tag_buffer_array[il] = 1;
      }
      
      bool done = false;
      bool initial_cycle = true;
      for (int ln = 0; hierarchy->levelCanBeRefined(ln) && !done;
           ++ln) {
         gridding_algorithm->makeFinerLevel(
            tag_buffer_array[ln],
            initial_cycle,
            init_cycle,
            init_time);
         done = !(hierarchy->finerLevelExists(ln));
      }

      // Setup timer manager for profiling code.
      tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
      std::shared_ptr<tbox::Timer> t_cvode_solve(
         tbox::TimerManager::getManager()->
         getTimer("PFiSM::time_integrator"));

      // Set up Visualization plot file writer
      int visit_number_procs_per_file=1;
      const std::string visit_dump_dirname 
         = "v."+input_filename.substr( 0, input_filename.rfind( "." ) );
      std::shared_ptr<appu::VisItDataWriter> visit_data_writer(
         new appu::VisItDataWriter(
            dim,
            "PFiSM VisIt Writer",
            visit_dump_dirname,
            visit_number_procs_per_file));
      pf_model->registerVisItDataWriter(visit_data_writer);

      // Setup solution vector, and initialize it.
      pf_model->setupSolutionVector(hierarchy);

      pf_model->setInitialConditions();

      /***********************************************************************
      * Setup CVODESolver object.
      ***********************************************************************/
      solv::CVODESolver* cvode_solver =
         new solv::CVODESolver("cvode_solver",
            pf_model.get(),
            uses_preconditioning);

      cvode_solver->setIterationType( CV_NEWTON );
      cvode_solver->setRelativeTolerance(relative_tolerance);
      cvode_solver->setAbsoluteTolerance(absolute_tolerance);
      cvode_solver->setMaximumNumberOfInternalSteps(max_steps);
      cvode_solver->setSteppingMethod(CV_ONE_STEP);
      cvode_solver->setMaximumLinearMultistepMethodOrder(max_order);
      if (uses_preconditioning) {
         cvode_solver->setPreconditioningType(PREC_LEFT);
      }

      cvode_solver->setInitialValueOfIndependentVariable(init_time);
      solv::SundialsAbstractVector* solution_vector =
         pf_model->getSolutionVector();
      cvode_solver->setInitialConditionVector(solution_vector);
      cvode_solver->initialize(solution_vector);

      /**********************************************************************
      * Start time-stepping.
      ***********************************************************************/

      std::vector<double> time(max_steps);
      std::vector<double> maxnorm(max_steps);

      double final_time = init_time;
      double print_time=0.;
      for (int interval = 1; interval <= max_steps; ++interval) {

         //tbox::plog << "interval = "<<interval<<std::endl;

         final_time += print_interval;
         cvode_solver->setFinalValueOfIndependentVariable(final_time, false);

         /*
          * Perform CVODE solve to the requested interval time.
          */
         t_cvode_solve->start();
         int ret = cvode_solver->solve();
         t_cvode_solve->stop();
         if( ret!=0 )tbox::plog << "return code = " << ret << std::endl;

         double actual_time =
            cvode_solver->getActualFinalValueOfIndependentVariable();
         double dt =
            cvode_solver->getStepSizeForLastInternalStep();
         tbox::pout << "# step = "<<interval
                    <<", time = "<<actual_time
                    <<", dt = "<<dt<<std::endl;

         /*
          * Print statistics
          */
         std::shared_ptr<solv::SAMRAIVectorReal<double> > y_result(
            solv::Sundials_SAMRAIVector::getSAMRAIVector(solution_vector));
         std::shared_ptr<hier::PatchHierarchy> result_hierarchy(
            y_result->getPatchHierarchy());

         time[interval - 1] = actual_time;
         maxnorm[interval - 1] = y_result->maxNorm();

         if (solution_logging) {
            tbox::plog << "CVODE stastistics:"<<std::endl;
            cvode_solver->printStatistics(tbox::pout);
         }

         if( actual_time > print_time && visit_flag ){
            visit_data_writer->writePlotData(
               result_hierarchy,
               interval,
               actual_time);
            print_time+=print_interval;
         }

      } // end of timestep loop

      /*************************************************************************
       * Write summary information
       ************************************************************************/
      if (solution_logging) {
         pf_model->printCounters(final_time);
      }
      if (solution_logging) {
         tbox::pout << "\n\nTimestep Summary of solution vector y()\n"
                    << "  time                   \t"
                    << "  Max Norm  \n";

         for (int interval = 0; interval < max_steps; ++interval) {
            tbox::pout.precision(18);
            tbox::pout << "  " << time[interval] << "  \t";
            tbox::pout.precision(6);
            tbox::pout << "  " << maxnorm[interval] << std::endl;
         }
      }

      tbox::TimerManager::getManager()->print(tbox::pout);

      tbox::pout<<"Delete solver..."<<std::endl;
      delete cvode_solver;

      tbox::pout<<"Clean up..."<<std::endl;
      pf_model.reset();
      gridding_algorithm.reset();
      error_est.reset();
      load_balancer.reset();
      box_generator.reset();
      hierarchy.reset();
      geometry.reset();
      visit_data_writer.reset();

      tbox::pout<<"Close scope..."<<std::endl;
   }

   tbox::pout<<"Finalize hypre..."<<std::endl;
   HYPRE_Finalize();

   // Shutdown SAMRAI and tbox::MPI.
   tbox::pout<<"Finalize SAMRAI..."<<std::endl;
   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return 0;
}
