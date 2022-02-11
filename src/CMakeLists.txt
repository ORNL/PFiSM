set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)


# we need the M4 macro processor
find_program( M4_EXECUTABLE m4 DOC "The M4 macro processor" )
if( NOT M4_EXECUTABLE )
   message( SEND_ERROR "Failed to find the M4 macro processor." )
endif( NOT M4_EXECUTABLE )

# script used to convert m4 files
MACRO(m4tof filename)
  IF(NOT EXISTS ${CMAKE_SOURCE_DIR}/src/fortran/${filename}.f)
    ADD_CUSTOM_COMMAND(
      OUTPUT ${CMAKE_SOURCE_DIR}/src/fortran/${filename}.f
      COMMAND m4 
      ARGS -DPDAT_FORTDIR=${CMAKE_SOURCE_DIR}/src/fortran
          ${CMAKE_SOURCE_DIR}/src/fortran/${filename}.m4 >
              ${CMAKE_SOURCE_DIR}/src/fortran/${filename}.f
      VERBATIM
      )
  ENDIF()
ENDMACRO(m4tof)

m4tof(rhs3d)
m4tof(computeCforPhase)

add_library(pfismlib STATIC
   arkode/ARKODESolver.C
   PFModel.cc
   PFModelCV.cc
   PFModelARK.cc
   PfmFACSolver.cc
   PFiSM_Hypre.cc
   samrai_internal/CellPoissonFACOps.C
   samrai_internal/CellPoissonHypreSolver.C
   samrai_internal/CellPoissonFACSolver.C
   tools.cc
   fortran/rhs3d.f
   fortran/computeCforPhase.f
   )

if(${PFISM_USE_GPU})
  set_source_files_properties(PFiSM_Hypre.C PROPERTIES LANGUAGE CUDA)
endif(${PFISM_USE_GPU})

target_compile_features(pfismlib PRIVATE cxx_std_11)

add_executable(main main.cc)

if(${PFISM_USE_GPU})
  set_target_properties(main PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  set_property(TARGET main PROPERTY CUDA_ARCHITECTURES
               ${PFISM_CUDA_ARCHITECTURES})
  set_property(TARGET pfismlib PROPERTY CUDA_ARCHITECTURES
               ${PFISM_CUDA_ARCHITECTURES})
endif(${PFISM_USE_GPU})

target_link_libraries(main PRIVATE pfismlib
                           ${PROJECT_LINK_LIBS_SAMRAI}
                           ${PROJECT_LINK_LIBS10}
                           ${PROJECT_LINK_LIBS11}
			   ${PROJECT_LINK_LIBS12}
                           ${PROJECT_LINK_LIBS_HDF5}
                      )

