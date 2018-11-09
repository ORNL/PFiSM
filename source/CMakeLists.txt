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

set(SOURCES ${CMAKE_SOURCE_DIR}/source/main.C
            ${CMAKE_SOURCE_DIR}/source/PFModel.C
            ${CMAKE_SOURCE_DIR}/source/PfmFACSolver.C
            ${CMAKE_SOURCE_DIR}/source/fortran/rhs3d.f
            ${CMAKE_SOURCE_DIR}/source/fortran/computeCforPhase.f )

add_executable(main ${SOURCES} )

add_library(PFiSM_Hypre SHARED ${CMAKE_SOURCE_DIR}/source/PFiSM_Hypre.C)

target_link_libraries(main ${PROJECT_LINK_LIBS_SAMRAI}
                           ${PROJECT_LINK_LIBS10}
                           ${PROJECT_LINK_LIBS11}
                           ${PROJECT_LINK_LIBS_HDF5}
                           ${PFiSM_Hypre}
                      )

