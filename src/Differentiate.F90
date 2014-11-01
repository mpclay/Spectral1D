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
!> @file Differentiate.F90
!> @author Matthew Clay
!> @brief Program to differentiate a signal to investigate the accuracy of using
!! spectral methods for differentiation.
PROGRAM Differentiate_p

   ! Required modules.
   USE Parameters_m,ONLY: IWP, RWP, PI
   USE IO_m,ONLY: FILE_NAME_LENGTH, FileName, WriteData
   USE Spectral_m,ONLY: SpectralSetup, Differentiate, SpectralFinalize

   IMPLICIT NONE

   ! User inputs for the simulation.
   !
   !> Number of grid points. Must be an even number.
   INTEGER(KIND=IWP),PARAMETER :: n = 4_IWP

   ! Variables to run the code.
   !
   !> Grid in physical space.
   REAL(KIND=RWP),DIMENSION(n) :: x
   !> Local grid spacing.
   REAL(KIND=RWP) :: dx
   !> Signal to be differentiated.
   REAL(KIND=RWP),DIMENSION(n) :: u
   !> Derivative of the signal.
   REAL(KIND=RWP),DIMENSION(n) :: dudx
   !> Exact derivative of the signal.
   REAL(KIND=RWP),DIMENSION(n) :: dudxExact
   !> LInfinity error between the numerical and exact derivative.
   REAL(KIND=RWP) :: lInf
   !> L2 error between the numerical and exact derivative.
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
   dudx(:) = 0.0_RWP

   ! Form the computational grid.
   DO i = 1, n
      x(i) = REAL(i, RWP)*2.0_RWP*PI/REAL(n, RWP)
   END DO
   dx = 2.0_RWP*PI/REAL(n, RWP)

   ! Initialize the signal and write it to a file.
   DO i = 1, n
      u(i) = EXP(SIN(x(i)))
   END DO
   CALL FileName('Signal', n, 'dat', fname)
   CALL WriteData(n, x, u, fname)

   ! Compute the exact derivative and write it to file.
   DO i = 1, n
      dudxExact(i) = EXP(SIN(x(i)))*COS(x(i))
   END DO
   CALL FileName('ExactDerivative', n, 'dat', fname)
   CALL WriteData(n, x, dudxExact, fname)
   
   ! Initialize the spectral module.
   CALL SpectralSetup(n)

   ! Differentiate the signal.
   CALL Differentiate(n, u, dudx)

   ! Write out the differentiated signal to file.
   CALL FileName('SpectralDerivative', n, 'dat', fname) 
   CALL WriteData(n, x, dudx, fname)

   ! Calculate the lInf and l2 errors in the numerical derivative approximation.
   lInf = 0.0_RWP
   l2 = 0.0_RWP
   DO i = 1, n
      diff = dudx(i) - dudxExact(i)
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

END PROGRAM Differentiate_p

