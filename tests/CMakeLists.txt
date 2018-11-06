set( TEST_DIR ${CMAKE_SOURCE_DIR}/tests )

add_executable(testMPI ${CMAKE_SOURCE_DIR}/tests/testMPI.C)

add_test(NAME testMPI
         COMMAND ${MPI_EXEC} ${MPIEXEC_NUMPROC_FLAG} 
                 ${MPIEXEC_MAX_NUMPROCS} 
                 ${CMAKE_CURRENT_BINARY_DIR}/testMPI)

target_link_libraries(testMPI ${PROJECT_LINK_LIBS1})

