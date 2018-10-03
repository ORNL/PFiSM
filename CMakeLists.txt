cmake_minimum_required(VERSION 3.1)
project (PFiSM)

#set(MPI_CXX_COMPILER CMAKE_CXX_COMPILER)
find_package(MPI)
if(MPI_FOUND)
  message(STATUS
    "MPI was found:\n"
    "  MPI_CXX_COMPILER: ${MPI_CXX_COMPILER}\n"
    "  MPIEXEC_NUMPROC_FLAG: ${MPIEXEC_NUMPROC_FLAG}\n"
    "  MPIEXEC_MAX_NUMPROCS: ${MPIEXEC_MAX_NUMPROCS}\n"
    "  MPIEXEC: ${MPIEXEC}\n"
 )
endif()

message(STATUS "MPI_EXEC: ${MPI_EXEC}\n")

enable_testing()

set (CMAKE_CXX_STANDARD 11)

enable_language(Fortran C CXX)

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -std=c++11")

add_subdirectory(source)

add_subdirectory(tests)



