!> @file Spectral.F90
!> @author Matthew Clay
!> @brief Module to compute the derivative of a signal using spectral methods.
MODULE Spectral_m

   ! Required modules.
   USE ISO_C_BINDING
   USE Parameters_m,ONLY: IWP, RWP

   IMPLICIT NONE

   ! FFTW procedure definitions.
   INCLUDE 'fftw3.f03'

END MODULE Spectral_m

