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
#
# Program for scalar advection.
SET(ADVECTION_SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/Advection.F90
    CACHE INTERNAL "Advection source files.")

# Program for scalar advection using FD.
SET(ADVECTION_FD_SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/AdvectionFD.F90
    CACHE INTERNAL "Advection FD source files.")

# Program to differentiate a signal.
SET(DIFFERENTIATE_SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/Differentiate.F90
    CACHE INTERNAL "Differentiation source files.")

# Program to tranform a signal to and from spectral space.
SET(TRANSFORM_SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/Transform.F90
    CACHE INTERNAL "Transform a signal to/from spectral space.")

# Program to investigate aliasing arising from nonlinear products.
SET(ALIASING_SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/Aliasing.F90
    CACHE INTERNAL "Study aliasing arising from nonlinear products.")

# Program to test Parseval's identity.
SET(PARSEVAL_SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/Parseval.F90
    CACHE INTERNAL "Check Parseval's identity.")

