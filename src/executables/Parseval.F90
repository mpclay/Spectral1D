! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! Copyright (C) 2014 by Matthew Clay <mpclay@gmail.com>
!
!> @file Parseval.F90
!> @author Matthew Clay
!> @brief Program to test Parseval's identity.
!!
!! In this program we verify Perseval's identity for a signal, which states that
!! the norm of a signal in physical space can be calculated from its Fourier
!! coefficients in spectral space.
PROGRAM Parseval_p

   ! Required modules.
   USE Parameters_m,ONLY: IWP, RWP, PI
   USE Spectral_m,ONLY: SpectralSetup, Spectrum, SpectralFinalize

   IMPLICIT NONE

   ! User inputs for the simulation.
   !
   !> Number of grid points. Must be an even number.
   INTEGER(KIND=IWP),PARAMETER :: n = 256_IWP
   !> Wavenumber for the signal: sin(k1*x).
   INTEGER(KIND=IWP),PARAMETER :: k1 = 64_IWP

   ! Variables to run the code.
   !
   !> Grid in physical space.
   REAL(KIND=RWP),DIMENSION(n) :: x
   !> Local grid spacing.
   REAL(KIND=RWP) :: dx
   !> Signal in physical space.
   REAL(KIND=RWP),DIMENSION(n) :: u
   !> Wavenumbers resolved by the grid.
   INTEGER(KIND=IWP),DIMENSION(n/2+1) :: kVec
   !> Magnitude of the Fourier coefficients of u.
   REAL(KIND=RWP),DIMENSION(n/2+1) :: sVec
   !> Square magnitude of the signal from physical space calculations.
   REAL(KIND=RWP) :: physSum
   !> Square magnitude of the signal from spectral space calculations.
   REAL(KIND=RWP) :: specSum

   ! Extraneous variables.
   !
   !> Looping index.
   INTEGER(KIND=IWP) :: i

   ! Print some information to the user.
   WRITE(*,100) '------------------------------------------------'
   WRITE(*,100) 'Parseval: A 1D Code to Verify Parsevals Identity'
   WRITE(*,100) '------------------------------------------------'
   WRITE(*,100) ''
   WRITE(*,100) 'The program will evaluate the square magnitude'
   WRITE(*,100) 'of a signal in both physical and spectral space'
   WRITE(*,100) 'to verify Parsevals identity. Results from both'
   WRITE(*,100) 'physical and spectral space will be reported.'
   WRITE(*,100) ''
   WRITE(*,100) 'Parameters'
   WRITE(*,100) '----------'
   WRITE(*,150) 'Num. of points:', n
   WRITE(*,150) 'Wavenumber:', k1
   100 FORMAT (A)
   150 FORMAT (A,T17,I3.3)

   ! Variable initialization.
   x(:) = 0.0_RWP
   u(:) = 0.0_RWP

   ! Form the computational grid.
   DO i = 1, n
      x(i) = REAL(i, RWP)*2.0_RWP*PI/REAL(n, RWP)
   END DO
   dx = 2.0_RWP*PI/REAL(n, RWP)

   ! Initialize the signal.
   DO i = 1, n
      u(i) = 10.0_RWP*SIN(REAL(k1, RWP)*x(i))
   END DO

   ! Initialize the spectral module.
   CALL SpectralSetup(n)

   ! Compute the spectrum for the signal.
   CALL Spectrum(n, u, kVec, sVec)

   ! Calculate the physical space square magnitude of the signal.
   physSum = 0.0_RWP
   DO i = 1, n
      physSum = physSum + u(i)*u(i)
   END DO
   physSum = physSum/REAL(n, RWP)

   ! Calculate the Fourier space square magnitude of the signal.
   specSum = 0.0_RWP
   DO i = 1, n/2+1
      specSum = specSum + sVec(i)
   END DO

   ! Write out the information to the user.
   WRITE(*,100) ''
   WRITE(*,100) 'Results'
   WRITE(*,100) '-------'
   WRITE(*,200) 'Physical space sum:', physSum
   WRITE(*,200) 'Spectral space sum:', specSum
   200 FORMAT (A,T20,ES15.8)

   ! Clean up the code.
   CALL SpectralFinalize()

END PROGRAM Parseval_p

