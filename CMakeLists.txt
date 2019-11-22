cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project (PFiSM LANGUAGES C CXX Fortran CUDA)

#set(MPI_CXX_COMPILER CMAKE_CXX_COMPILER)
find_package(MPI REQUIRED)
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

set(HDF5_USE_STATIC_LIBRARIES ON)
find_package(HDF5 REQUIRED)
message(STATUS ${HDF5_VERSION} )
include_directories(${HDF5_INCLUDE_DIR})
set(PROJECT_LINK_LIBS_HDF5 ${HDF5_LIBRARIES})

link_directories( ${SAMRAI_DIR}/lib )
link_directories( ${HYPRE_DIR}/lib )
link_directories( ${SUNDIALS_DIR}/lib )

include_directories( ${SAMRAI_DIR}/include )
include_directories( ${HYPRE_DIR}/include )
include_directories( ${SUNDIALS_DIR}/include )

enable_testing()

set (CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -std=c++11")

add_subdirectory(source)

#add_subdirectory(tests)

