#include "TestHyprePoisson.h"

#include <string>
#include <iostream>
#include <mpi.h>

int main(int argc, char *argv[])
{
   // Initialize MPI
   MPI_Init(&argc, &argv);

   int ret = 0;

   {
      // create main object
      TestHyprePoisson test(MPI_COMM_WORLD);

      /*
       * Process command line arguments.  For each run, the input
       * filename must be specified.  Usage is:
       *
       *    executable <input file name>
       */
      if (argc != 2) {
         std::cerr << "USAGE:  " << argv[0] << " <input filename> "
                   << std::endl;
      }

      std::string input_filename = argv[1];

      ret = test.run(input_filename);
   }

   MPI_Finalize();

   return ret;
}
