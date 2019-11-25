set( TEST_DIR ${CMAKE_SOURCE_DIR}/tests )

if(${PFISM_USE_GPU})
  set_source_files_properties(TestHyprePoisson.C PROPERTIES LANGUAGE CUDA)
  set_source_files_properties(HyprePoisson.C PROPERTIES LANGUAGE CUDA)
endif(${PFISM_USE_GPU})

add_executable(testMPI testMPI.C)
add_executable(testHyprePoisson mainTestHyprePoisson.C)

add_library(testHyprePoissonLib STATIC
   TestHyprePoisson.C
   HyprePoisson.C
   hyprepoisson3d.f)

add_test(NAME testMPI
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testMPI)
add_test(NAME testHyprePoisson
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testHyprePoisson 3d.input
         WORKING_DIRECTORY ${TEST_DIR})

target_link_libraries(testMPI ${PROJECT_LINK_LIBS1})
target_link_libraries(testHyprePoisson testHyprePoissonLib
                                       ${PROJECT_LINK_LIBS_SAMRAI}
                                       ${PROJECT_LINK_LIBS10}
                                       ${PROJECT_LINK_LIBS_HDF5})

target_include_directories(testHyprePoissonLib PUBLIC
                           ${CMAKE_SOURCE_DIR}/source)

target_include_directories(testHyprePoisson PUBLIC ${CMAKE_SOURCE_DIR}/source)

target_compile_features(testHyprePoissonLib PRIVATE cxx_std_11)
target_compile_features(testHyprePoisson PRIVATE cxx_std_11)

if(${PFISM_USE_GPU})
  set_target_properties(testHyprePoisson PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif(${PFISM_USE_GPU})
