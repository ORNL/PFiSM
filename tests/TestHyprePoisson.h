#ifndef included_TestHyprePoisson
#define included_TestHyprePoisson

#include <mpi.h>
#include <string>

class TestHyprePoisson
{
public:
   TestHyprePoisson(MPI_Comm comm);

   ~TestHyprePoisson();

   int run(const std::string input_filename);
};

#endif
