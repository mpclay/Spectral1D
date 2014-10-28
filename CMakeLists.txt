CMAKE_MINIMUM_REQUIRED(VERSION 2.8.9)

PROJECT(Spectral)
ENABLE_LANGUAGE(Fortran)

# Main source directory.
ADD_SUBDIRECTORY(src)

# Since FFTW doesn't have a FindFFTW.cmake, we include it in the library. The
# user must have the environment variable FFTW_ROOT set properly.
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

# Find the FFTW library, which is required.
SET(FFTW_USE_STATIC_LIBS TRUE)
SET(FFTW_FORTRAN_ENABLED TRUE)
SET(FFTW_USE_MPI FALSE)
FIND_PACKAGE(FFTW REQUIRED)
INCLUDE_DIRECTORIES(${FFTW_INCLUDES})

# The main program.
OPTION(BUILD_ADVECTION "Build the scalar advection program." OFF)
IF(BUILD_ADVECTION)
   ADD_EXECUTABLE(advection
                  ${ADVECTION_SOURCE_FILES}
                  ${CORE_SOURCE_FILES})
   TARGET_LINK_LIBRARIES(advection ${FFTW_LIBRARIES})
ENDIF(BUILD_ADVECTION)

# The differentiation program.
OPTION(BUILD_DIFFERENTIATE "Build the differentiation program." OFF)
IF(BUILD_DIFFERENTIATE)
   ADD_EXECUTABLE(differentiate
                  ${DIFFERENTIATE_SOURCE_FILES}
                  ${CORE_SOURCE_FILES})
   TARGET_LINK_LIBRARIES(differentiate ${FFTW_LIBRARIES})
ENDIF(BUILD_DIFFERENTIATE)

# The tranform program.
OPTION(BUILD_TRANSFORM "Build the transform program." OFF)
IF(BUILD_TRANSFORM)
   ADD_EXECUTABLE(transform
                  ${TRANSFORM_SOURCE_FILES}
                  ${CORE_SOURCE_FILES})
   TARGET_LINK_LIBRARIES(transform ${FFTW_LIBRARIES})
ENDIF(BUILD_TRANSFORM)

