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
!> @file Aliasing.F90
!> @author Matthew Clay
!> @brief Program to investigate aliasing arising from nonlinear products.
!!
!! In this program we explore how nonlinear products lead to aliasing in
!! spectral codes. The problem setup is as follows:
!!
!!    1. Establish a grid with 16 grid points, which can resolve sinusoids up
!!       to the 8th wavenumber.
!!    2. Multiply two signals, each with wavenumbers of 6, which yields a wave
!!       with wavenumber 12. The resulting wave cannot be resolved on the grid,
!!       and should be aliased back to the -4th wavenumber.
!!    3. Calculate the spectrum for each of these signals. We expect that to see
!!       spikes at 6 for the two initial waves, and spike at 4 for the aliased
!!       wave.
PROGRAM Aliasing_p

   ! Required modules.
   USE Parameters_m,ONLY: IWP, RWP, PI
   USE IO_m,ONLY: FILE_NAME_LENGTH, FileName, WriteData
   USE Spectral_m,ONLY: SpectralSetup, Transform, Spectrum, SpectralFinalize

   IMPLICIT NONE

   ! User inputs for the simulation.
   !
   !> Number of grid points. Must be an even number.
   INTEGER(KIND=IWP),PARAMETER :: n = 16_IWP
   !> Wavenumber for the first wave: sin(k1*x).
   INTEGER(KIND=IWP),PARAMETER :: k1 = 6_IWP
   !> Wavenumber for the second wave: cos(k2*x).
   INTEGER(KIND=IWP),PARAMETER :: k2 = 6_IWP

   ! Variables to run the code.
   !
   !> Grid in physical space.
   REAL(KIND=RWP),DIMENSION(n) :: x
   !> Local grid spacing.
   REAL(KIND=RWP) :: dx
   !> First signal in physical space.
   REAL(KIND=RWP),DIMENSION(n) :: u1
   !> Second signal in physical space.
   REAL(KIND=RWP),DIMENSION(n) :: u2
   !> Multiplication of u1 and u2 in physical space.
   REAL(KIND=RWP),DIMENSION(n) :: v
   !> Output when v was transformed to and from spectral space.
   REAL(KIND=RWP),DIMENSION(n) :: vOut
   !> Wavenumbers resolved by the grid.
   INTEGER(KIND=IWP),DIMENSION(n/2+1) :: kVec
   !> Spectral output from the spectrum module.
   REAL(KIND=RWP),DIMENSION(n/2+1) :: sVec
   !> File name variable.
   CHARACTER(LEN=FILE_NAME_LENGTH) :: fname

   ! Extraneous variables.
   !
   !> Looping index.
   INTEGER(KIND=IWP) :: i

   ! Variable initialization.
   x(:) = 0.0_RWP
   u1(:) = 0.0_RWP
   u2(:) = 0.0_RWP
   v(:) = 0.0_RWP
   vOut(:) = 0.0_RWP

   ! Form the computational grid.
   DO i = 1, n
      x(i) = REAL(i, RWP)*2.0_RWP*PI/REAL(n, RWP)
   END DO
   dx = 2.0_RWP*PI/REAL(n, RWP)

   ! Initialize the signals and write them to a file.
   DO i = 1, n
      u1(i) = SIN(REAL(k1, RWP)*x(i))
      u2(i) = COS(REAL(k2, RWP)*x(i))
      v(i) = u1(i)*u2(i)
   END DO
   CALL FileName('u1', 'dat', fname)
   CALL WriteData(n, x, u1, fname)
   CALL FileName('u2', 'dat', fname)
   CALL WriteData(n, x, u2, fname)
   CALL FileName('v', 'dat', fname)
   CALL WriteData(n, x, v, fname)

   ! Initialize the spectral module.
   CALL SpectralSetup(n)

   ! Compute the spectrum for each of the signals and write them to file.
   !
   ! u1
   CALL Spectrum(n, u1, kVec, sVec)
   CALL FileName('u1Spectrum', 'dat', fname)
   CALL WriteData(n/2+1, kVec, sVec, fname)
   !
   ! u2
   CALL Spectrum(n, u2, kVec, sVec)
   CALL FileName('u2Spectrum', 'dat', fname)
   CALL WriteData(n/2+1, kVec, sVec, fname)
   !
   ! v
   CALL Spectrum(n, v, kVec, sVec)
   CALL FileName('vSpectrum', 'dat', fname)
   CALL WriteData(n/2+1, kVec, sVec, fname)

   ! Clean up the code.
   CALL SpectralFinalize()

END PROGRAM Aliasing_p

