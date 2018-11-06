set( TEST_DIR ${CMAKE_SOURCE_DIR}/tests )

add_executable(testMPI ${CMAKE_SOURCE_DIR}/tests/testMPI.C)
add_executable(testHyprePoisson ${CMAKE_SOURCE_DIR}/tests/testHyprePoisson.C
                                ${CMAKE_SOURCE_DIR}/tests/HyprePoisson.C
                                ${CMAKE_SOURCE_DIR}/tests/hyprepoisson3d.f)

add_test(NAME testMPI
         COMMAND ${MPI_EXEC} ${MPIEXEC_NUMPROC_FLAG} 
                 ${MPIEXEC_MAX_NUMPROCS} 
                 ${CMAKE_CURRENT_BINARY_DIR}/testMPI)
add_test(NAME testHyprePoisson
         COMMAND ${MPI_EXEC} ${MPIEXEC_NUMPROC_FLAG} 
                 ${MPIEXEC_MAX_NUMPROCS} 
                 ${CMAKE_CURRENT_BINARY_DIR}/testHyprePoisson 3d.input
         WORKING_DIRECTORY ${TEST_DIR})

target_link_libraries(testMPI ${PROJECT_LINK_LIBS1})
target_link_libraries(testHyprePoisson ${PROJECT_LINK_LIBS_SAMRAI}
                                       ${PROJECT_LINK_LIBS10}
                                       ${PROJECT_LINK_LIBS_HDF5})

