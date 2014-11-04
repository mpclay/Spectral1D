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
!> @file Transform.F90
!> @author Matthew Clay
!> @brief Program to transform a signal to and from spectral space.
!!
!! In this program we want to investigate what happens to a signal when it is
!! transformed to and from spectral space.
PROGRAM Transform_p

   ! Required modules.
   USE Parameters_m,ONLY: IWP, RWP, PI
   USE IO_m,ONLY: FILE_NAME_LENGTH, FileName, WriteData
   USE Spectral_m,ONLY: SpectralSetup, Transform, SpectralFinalize

   IMPLICIT NONE

   ! User inputs for the simulation.
   !
   !> Number of grid points. Must be an even number.
   INTEGER(KIND=IWP),PARAMETER :: n = 256_IWP

   ! Variables to run the code.
   !
   !> Grid in physical space.
   REAL(KIND=RWP),DIMENSION(n) :: x
   !> Local grid spacing.
   REAL(KIND=RWP) :: dx
   !> Signal to be transformed.
   REAL(KIND=RWP),DIMENSION(n) :: u
   !> Output when u was transformed to and from spectral space.
   REAL(KIND=RWP),DIMENSION(n) :: out
   !> LInfinity error between the transformed and exact signal.
   REAL(KIND=RWP) :: lInf
   !> L2 error between the transformed and exact derivative.
   REAL(KIND=RWP) :: l2
   !> File name variable.
   CHARACTER(LEN=FILE_NAME_LENGTH) :: fname

   ! Extraneous variables.
   !
   !> Difference between the exact and numerical derivatives.
   REAL(KIND=RWP) :: diff
   !> Looping index.
   INTEGER(KIND=IWP) :: i

   ! Variable initialization.
   x(:) = 0.0_RWP
   u(:) = 0.0_RWP
   out(:) = 0.0_RWP

   ! Form the computational grid.
   DO i = 1, n
      x(i) = REAL(i, RWP)*2.0_RWP*PI/REAL(n, RWP)
   END DO
   dx = 2.0_RWP*PI/REAL(n, RWP)

   ! Initialize the signal and write it to a file.
   DO i = 1, n
      u(i) = EXP(SIN(9.0_RWP*x(i))*COS(5.0_RWP*x(i))**2)**4
   END DO
   CALL FileName('Input', 'dat', fname)
   CALL WriteData(n, x, u, fname)

   ! Initialize the spectral module.
   CALL SpectralSetup(n)

   ! Transform the signal to and from spectral space.
   CALL Transform(n, u, out)

   ! Write out the transformed signal to file.
   CALL FileName('Output', 'dat', fname)
   CALL WriteData(n, x, out, fname)

   ! Calculate the lInf and l2 errors between the input and output signal.
   lInf = 0.0_RWP
   l2 = 0.0_RWP
   DO i = 1, n
      diff = u(i) - out(i)
      l2 = l2 + diff**2
      IF (ABS(diff) > lInf) THEN
         lInf = ABS(diff)
      END IF
   END DO
   l2 = SQRT(l2*dx)
   WRITE(*,100) 'L2 Error:', l2
   WRITE(*,100) 'LInf Error:', lInf
   100 FORMAT (A,T16,ES15.8)

   ! Clean up the code.
   CALL SpectralFinalize()

END PROGRAM Transform_p

