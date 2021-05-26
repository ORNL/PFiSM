cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project (PFiSM LANGUAGES C CXX Fortran)

set(PFISM_CUDA_ARCHITECTURES "" CACHE STRING "target gpu architectures")
if(${PFISM_USE_GPU})
  enable_language(CUDA)
endif(${PFISM_USE_GPU})

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

set(HDF5_USE_STATIC_LIBRARIES ON)
find_package(HDF5 REQUIRED)
message(STATUS ${HDF5_VERSION} )
include_directories(${HDF5_INCLUDE_DIR})
set(PROJECT_LINK_LIBS_HDF5 ${HDF5_LIBRARIES})

link_directories( ${SAMRAI_DIR}/lib )
link_directories( ${HYPRE_DIR}/lib )
link_directories( ${SUNDIALS_DIR}/lib64 )

include_directories( ${SAMRAI_DIR}/include )
include_directories( ${HYPRE_DIR}/include )
include_directories( ${SUNDIALS_DIR}/include )

set( PROJECT_LINK_LIBS1 libSAMRAI_tbox.a )
set( PROJECT_LINK_LIBS2 libSAMRAI_math.a )
set( PROJECT_LINK_LIBS3 libSAMRAI_hier.a )
set( PROJECT_LINK_LIBS4 libSAMRAI_pdat.a )
set( PROJECT_LINK_LIBS5 libSAMRAI_xfer.a )
set( PROJECT_LINK_LIBS6 libSAMRAI_geom.a )
set( PROJECT_LINK_LIBS7 libSAMRAI_mesh.a )
set( PROJECT_LINK_LIBS8 libSAMRAI_solv.a )
set( PROJECT_LINK_LIBS9 libSAMRAI_appu.a )
set( PROJECT_LINK_LIBS_SAMRAI ${PROJECT_LINK_LIBS9}
                              ${PROJECT_LINK_LIBS8}
                              ${PROJECT_LINK_LIBS7}
                              ${PROJECT_LINK_LIBS6}
                              ${PROJECT_LINK_LIBS5}
                              ${PROJECT_LINK_LIBS4}
                              ${PROJECT_LINK_LIBS3}
                              ${PROJECT_LINK_LIBS2}
                              ${PROJECT_LINK_LIBS1} )

set( PROJECT_LINK_LIBS10 libHYPRE.so )
set( PROJECT_LINK_LIBS11 libsundials_cvode.a )

enable_testing()

set (CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -std=c++11")

add_subdirectory(source)

add_subdirectory(tests)

