# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# Copyright (C) 2014 by Matthew Clay <mpclay@gmail.com>

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
   INSTALL(TARGETS advection DESTINATION ${CMAKE_INSTALL_PREFIX})
ENDIF(BUILD_ADVECTION)

# Build the FD advection program.
OPTION(BUILD_ADVECTION_FD "Build the scalar advection FD program." OFF)
IF(BUILD_ADVECTION_FD)
   ADD_EXECUTABLE(advectionFD
                  ${ADVECTION_FD_SOURCE_FILES}
                  ${CORE_SOURCE_FILES})
   # In theory this link isn't needed, as DFTs are not being performed. I don't,
   # however, want to make another CORE_SOURCE_FILES with the spectral module
   # stripped out just for this one executable.
   TARGET_LINK_LIBRARIES(advectionFD ${FFTW_LIBRARIES})
   INSTALL(TARGETS advectionFD DESTINATION ${CMAKE_INSTALL_PREFIX})
ENDIF(BUILD_ADVECTION_FD)

# The differentiation program.
OPTION(BUILD_DIFFERENTIATE "Build the differentiation program." OFF)
IF(BUILD_DIFFERENTIATE)
   ADD_EXECUTABLE(differentiate
                  ${DIFFERENTIATE_SOURCE_FILES}
                  ${CORE_SOURCE_FILES})
   TARGET_LINK_LIBRARIES(differentiate ${FFTW_LIBRARIES})
   INSTALL(TARGETS differentiate DESTINATION ${CMAKE_INSTALL_PREFIX})
ENDIF(BUILD_DIFFERENTIATE)

# The tranform program.
OPTION(BUILD_TRANSFORM "Build the transform program." OFF)
IF(BUILD_TRANSFORM)
   ADD_EXECUTABLE(transform
                  ${TRANSFORM_SOURCE_FILES}
                  ${CORE_SOURCE_FILES})
   TARGET_LINK_LIBRARIES(transform ${FFTW_LIBRARIES})
   INSTALL(TARGETS transform DESTINATION ${CMAKE_INSTALL_PREFIX})
ENDIF(BUILD_TRANSFORM)

# The aliasing program.
OPTION(BUILD_ALIASING "Build the aliasing program." OFF)
IF(BUILD_ALIASING)
   ADD_EXECUTABLE(aliasing
                  ${ALIASING_SOURCE_FILES}
                  ${CORE_SOURCE_FILES})
   TARGET_LINK_LIBRARIES(aliasing ${FFTW_LIBRARIES})
   INSTALL(TARGETS aliasing DESTINATION ${CMAKE_INSTALL_PREFIX})
ENDIF(BUILD_ALIASING)

# The Parseval program.
OPTION(BUILD_PARSEVAL "Build the Parseval program." OFF)
IF(BUILD_PARSEVAL)
   ADD_EXECUTABLE(parseval
                  ${PARSEVAL_SOURCE_FILES}
                  ${CORE_SOURCE_FILES})
   TARGET_LINK_LIBRARIES(parseval ${FFTW_LIBRARIES})
   INSTALL(TARGETS parseval DESTINATION ${CMAKE_INSTALL_PREFIX})
ENDIF(BUILD_PARSEVAL)

