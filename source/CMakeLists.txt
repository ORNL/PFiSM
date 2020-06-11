set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)


# we need the M4 macro processor
find_program( M4_EXECUTABLE m4 DOC "The M4 macro processor" )
if( NOT M4_EXECUTABLE )
   message( SEND_ERROR "Failed to find the M4 macro processor." )
endif( NOT M4_EXECUTABLE )

# script used to convert m4 files
MACRO(m4tof filename)
  IF(NOT EXISTS ${CMAKE_SOURCE_DIR}/source/fortran/${filename}.f)
    ADD_CUSTOM_COMMAND(
      OUTPUT ${CMAKE_SOURCE_DIR}/source/fortran/${filename}.f
      COMMAND m4 
      ARGS -DPDAT_FORTDIR=${SAMRAI_DIR}/include/SAMRAI/pdat/fortran 
          ${CMAKE_SOURCE_DIR}/source/fortran/${filename}.m4 >
              ${CMAKE_SOURCE_DIR}/source/fortran/${filename}.f
      VERBATIM
      )
  ENDIF()
ENDMACRO(m4tof)

m4tof(rhs3d)
m4tof(computeCforPhase)

add_library(pfismlib STATIC
   PFModel.C
   PfmFACSolver.C
   PFiSM_Hypre.C
   CellPoissonFACOps.C
   CellPoissonHypreSolver.C
   CellPoissonFACSolver.C
   tools.C
   fortran/rhs3d.f
   fortran/computeCforPhase.f
   )

if(${PFISM_USE_GPU})
  set_source_files_properties(PFiSM_Hypre.C PROPERTIES LANGUAGE CUDA)
  set_source_files_properties(CellPoissonHypreSolver.C PROPERTIES LANGUAGE CUDA)
endif(${PFISM_USE_GPU})

target_compile_features(pfismlib PRIVATE cxx_std_11)

add_executable(main main.C)

if(${PFISM_USE_GPU})
  set_target_properties(main PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif(${PFISM_USE_GPU})

target_link_libraries(main PRIVATE pfismlib
                           ${PROJECT_LINK_LIBS_SAMRAI}
                           ${PROJECT_LINK_LIBS10}
                           ${PROJECT_LINK_LIBS11}
                           ${PROJECT_LINK_LIBS_HDF5}
                      )

