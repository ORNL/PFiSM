set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

link_directories( ${SAMRAI_DIR}/lib )
link_directories( ${HYPRE_DIR}/lib )
link_directories( ${SUNDIALS_DIR}/lib )

set( PROJECT_LINK_LIBS1 libSAMRAI_tbox.a )
set( PROJECT_LINK_LIBS2 libSAMRAI_math.a )
set( PROJECT_LINK_LIBS3 libSAMRAI_hier.a )
set( PROJECT_LINK_LIBS4 libSAMRAI_pdat.a )
set( PROJECT_LINK_LIBS5 libSAMRAI_xfer.a )
set( PROJECT_LINK_LIBS6 libSAMRAI_geom.a )
set( PROJECT_LINK_LIBS7 libSAMRAI_mesh.a )
set( PROJECT_LINK_LIBS8 libSAMRAI_solv.a )
set( PROJECT_LINK_LIBS9 libSAMRAI_appu.a )

set( PROJECT_LINK_LIBS10 libHYPRE.so )

set( PROJECT_LINK_LIBS11 libsundials_cvode.a )

include_directories( ${SAMRAI_DIR}/include )
include_directories( ${HYPRE_DIR}/include )
include_directories( ${SUNDIALS_DIR}/include )

set(HDF5_USE_STATIC_LIBRARIES ON)
find_package(HDF5 REQUIRED)
message(STATUS ${HDF5_VERSION} )
include_directories(${HDF5_INCLUDE_DIR})
set(_hdf5_libs hdf5 hdf5_cpp)

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
          ${CMAKE_SOURCE_DIR}/fortran/${filename}.m4 > ${CMAKE_SOURCE_DIR}/fortran/${filename}.f
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

target_link_libraries(main ${PROJECT_LINK_LIBS9} 
                           ${PROJECT_LINK_LIBS8} 
                           ${PROJECT_LINK_LIBS7} 
                           ${PROJECT_LINK_LIBS6} 
                           ${PROJECT_LINK_LIBS5} 
                           ${PROJECT_LINK_LIBS4}
                           ${PROJECT_LINK_LIBS3}
                           ${PROJECT_LINK_LIBS2}
                           ${PROJECT_LINK_LIBS1}
                           ${PROJECT_LINK_LIBS10}
                           ${PROJECT_LINK_LIBS11}
                           ${_hdf5_libs}
                      )

