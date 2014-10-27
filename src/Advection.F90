!> @file Advection.F90
!> @author Matthew Clay
!> @brief Execution code for the scalar advection equation.
PROGRAM Advection_p

   ! Required modules.
   USE Parameters_m
   USE Spectral_m

   IMPLICIT NONE

   ! User inputs for the simulation.
   !
   !> Number of grid points. Must be an even number.
   INTEGER(KIND=IWP),PARAMETER :: n = 64_IWP

   ! Initialize the spectral module.
   CALL SpectralSetup(n)

   ! Clean up the code.
   CALL SpectralFinalize()

END PROGRAM Advection_p

