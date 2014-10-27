# Core source files that are program independent
SET(CORE_SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/Parameters.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/IO.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/Spectral.F90
    CACHE INTERNAL "Core source files.")

# Program for scalar advection.
SET(ADVECTION_SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/Advection.F90
    CACHE INTERNAL "Advection source files.")

# Program to differentiate a signal.
SET(DIFFERENTIATE_SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/Differentiate.F90
    CACHE INTERNAL "Differentiation source files.")

