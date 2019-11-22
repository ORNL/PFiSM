set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

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
   fortran/rhs3d.f
   fortran/computeCforPhase.f
   )

set_source_files_properties(PfmFACSolver.C PROPERTIES LANGUAGE CUDA)
set_source_files_properties(PFiSM_Hypre.C PROPERTIES LANGUAGE CUDA)
set_source_files_properties(CellPoissonFACOps.C PROPERTIES LANGUAGE CUDA)
set_source_files_properties(CellPoissonHypreSolver.C PROPERTIES LANGUAGE CUDA)
set_source_files_properties(CellPoissonFACSolver.C PROPERTIES LANGUAGE CUDA)

target_compile_features(pfismlib PUBLIC cxx_std_11)

add_executable(main main.C)

set_target_properties(main PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(main PRIVATE pfismlib
                           ${PROJECT_LINK_LIBS_SAMRAI}
                           ${PROJECT_LINK_LIBS10}
                           ${PROJECT_LINK_LIBS11}
                           ${PROJECT_LINK_LIBS_HDF5}
                      )

