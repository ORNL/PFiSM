set( TEST_DIR ${CMAKE_SOURCE_DIR}/tests )

if(${PFISM_USE_GPU})
  set_source_files_properties(TestHyprePoisson.C PROPERTIES LANGUAGE CUDA)
  set_source_files_properties(HyprePoisson.C PROPERTIES LANGUAGE CUDA)
  set_source_files_properties(${CMAKE_SOURCE_DIR}/src/samrai_internal/CellPoissonHypreSolver.C
                              PROPERTIES LANGUAGE CUDA)
  set_source_files_properties(${CMAKE_SOURCE_DIR}/src/PFiSM_Hypre.C
                              PROPERTIES LANGUAGE CUDA)
endif(${PFISM_USE_GPU})

add_executable(testMPI testMPI.C)
add_executable(testHyprePoisson mainTestHyprePoisson.cc)

add_library(testHyprePoissonLib STATIC
   TestHyprePoisson.cc
   ${CMAKE_SOURCE_DIR}/src/samrai_internal/CellPoissonHypreSolver.C
   ${CMAKE_SOURCE_DIR}/src/PFiSM_Hypre.cc
   HyprePoisson.cc
   hyprepoisson3d.f)

#unit tests
add_test(NAME testMPI
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testMPI)
add_test(NAME testHyprePoisson
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testHyprePoisson 3d.input
         WORKING_DIRECTORY ${TEST_DIR})

#regression tests
add_test(NAME PlanarFront
         COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/PlanarFront/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/main
         ${CMAKE_CURRENT_SOURCE_DIR}/PlanarFront/pfm3d.input
         )
add_test(NAME PlanarFrontNoP
         COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/PlanarFront/testNoP.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/main
         ${CMAKE_CURRENT_SOURCE_DIR}/PlanarFront/pfm3dNoP.input
         )
add_test(NAME PlanarFrontUniformT
         COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/PlanarFrontUniformT/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/main
         ${CMAKE_CURRENT_SOURCE_DIR}/PlanarFrontUniformT/3d.input
         )
add_test(NAME PlanarFrontIMEX
         COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/PlanarFront/testIMEX.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/main
         ${CMAKE_CURRENT_SOURCE_DIR}/PlanarFront/pfm3dIMEX.input
         )
add_test(NAME PlanarFrontEX
         COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/PlanarFront/testEX.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/main
         ${CMAKE_CURRENT_SOURCE_DIR}/PlanarFront/pfm3dExplicit.input
         )
add_test(NAME PlanarFrontARKimplicit
         COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/PlanarFront/testARKimplicit.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/main
         ${CMAKE_CURRENT_SOURCE_DIR}/PlanarFront/pfm3dARKimplicit.input
         )

target_link_libraries(testMPI ${PROJECT_LINK_LIBS1})
target_link_libraries(testHyprePoisson testHyprePoissonLib
                                       ${PROJECT_LINK_LIBS_SAMRAI}
                                       ${PROJECT_LINK_LIBS10}
                                       ${PROJECT_LINK_LIBS_HDF5})

target_include_directories(testHyprePoissonLib PUBLIC
                           ${CMAKE_SOURCE_DIR}/src)

target_include_directories(testHyprePoisson PUBLIC ${CMAKE_SOURCE_DIR}/src)

target_compile_features(testHyprePoissonLib PRIVATE cxx_std_11)
target_compile_features(testHyprePoisson PRIVATE cxx_std_11)

if(${PFISM_USE_GPU})
  set_target_properties(testHyprePoisson PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  set_property(TARGET testHyprePoisson PROPERTY CUDA_ARCHITECTURES
               ${PFISM_CUDA_ARCHITECTURES})
  set_property(TARGET testHyprePoissonLib PROPERTY CUDA_ARCHITECTURES
               ${PFISM_CUDA_ARCHITECTURES})
endif(${PFISM_USE_GPU})
