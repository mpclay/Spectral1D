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
!!
!! The code will also output a 2D order finite difference approximation to the
!! derivative for comparison.
!!
!! In the future this program should accept a signal as an input parameter,
!! which will make it more flexible than hard-coding the signals in this file.
PROGRAM Differentiate_p

   ! Required modules.
   USE Parameters_m,ONLY: IWP, RWP, PI
   USE IO_m,ONLY: FILE_NAME_LENGTH, FileName, WriteData
   USE Spectral_m,ONLY: SpectralSetup, Differentiate, SpectralFinalize

   IMPLICIT NONE

   ! The types of signals available for differentiation.
   !
   !> The function EXP(SIN(X)).
   INTEGER(KIND=IWP),PARAMETER :: EXP_SIN_X = 0_IWP
   !> A square wave.
   INTEGER(KIND=IWP),PARAMETER :: SQUARE_WAVE = 1_IWP

   ! User inputs for the simulation.
   !
   !> Which wave to differentiate.
   INTEGER(KIND=IWP),PARAMETER :: wave = EXP_SIN_X
   !> Number of grid points. Must be an even number.
   INTEGER(KIND=IWP),PARAMETER :: n = 256_IWP

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
   !> Derivative of the signal using FD.
   REAL(KIND=RWP),DIMENSION(n) :: dudxFD
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

   ! Print some information to the user.
   WRITE(*,100) '---------------------------------------------------------'
   WRITE(*,100) 'Differentiate: A 1D Code to Investigate Spectral Accuracy'
   WRITE(*,100) '---------------------------------------------------------'
   WRITE(*,100) ''
   WRITE(*,100) 'The program differentiates a signal using spectral means,'
   WRITE(*,100) 'and outputs the result for the user. If an analytical'
   WRITE(*,100) 'result for the derivative is available, the error between'
   WRITE(*,100) 'the numerical and exact derivative will be reported.'
   WRITE(*,100) ''
   WRITE(*,100) 'Parameters'
   WRITE(*,100) '----------'
   WRITE(*,150) 'Num. of points:', n
   IF (wave == EXP_SIN_X) THEN
      WRITE(*,200) 'Waveform:', 'EXP(SIN(X))'
   ELSE IF (wave == SQUARE_WAVE) THEN
      WRITE(*,200) 'Waveform:', 'square wave'
   END IF
   100 FORMAT (A)
   150 FORMAT (A,T17,I3.3)
   200 FORMAT (A,T17,A)

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
   SELECT CASE (wave)
      CASE (EXP_SIN_X)
         DO i = 1, n
            u(i) = EXP(SIN(x(i)))
         END DO
      CASE (SQUARE_WAVE)
         DO i = 1, n
            IF (x(i) < 2.0_RWP*PI/3.0_RWP) THEN
               u(i) = 0.0_RWP
            ELSE IF (x(i) <= 4.0_RWP*PI/3.0_RWP) THEN
               u(i) = 1.0_RWP
            ELSE
               u(i) = 0.0_RWP
            END IF
         END DO
      CASE DEFAULT
         WRITE(*,10) 'Invalid IC option. Defaulting to EXP(SIN(x)).'
         10 FORMAT (A)
   END SELECT
   CALL FileName('Signal', n, 'dat', fname)
   CALL WriteData(n, x, u, fname)

   ! Compute the exact derivative and write it to file.
   IF (wave == EXP_SIN_X) THEN
      DO i = 1, n
         dudxExact(i) = EXP(SIN(x(i)))*COS(x(i))
      END DO
      CALL FileName('ExactDerivative', n, 'dat', fname)
      CALL WriteData(n, x, dudxExact, fname)
   END IF
   
   ! Initialize the spectral module.
   CALL SpectralSetup(n)

   ! Differentiate the signal.
   CALL Differentiate(n, u, dudx)

   ! Compute a 2D order FD approximation to the derivative for comparison.
   DO i = 2, n-1
      dudxFD(i) = (u(i+1) - u(i-1))/(2.0_RWP*dx)
   END DO
   dudxFD(1) = (u(2) - u(n))/(2.0_RWP*dx)
   dudxFD(n) = (u(1) - u(n-1))/(2.0_RWP*dx)

   ! Write out the differentiated signal to file.
   CALL FileName('SpectralDerivative', n, 'dat', fname) 
   CALL WriteData(n, x, dudx, fname)

   ! Write out the FD approximation to file.
   CALL FileName('FDDerivative', n, 'dat', fname)
   CALL WriteData(n, x, dudxFD, fname)

   ! Calculate the lInf and l2 errors in the numerical derivative approximation.
   IF (wave == EXP_SIN_X) THEN
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
      WRITE(*,100) ''
      WRITE(*,100) 'Errors'
      WRITE(*,100) '------'
      WRITE(*,250) 'L2 Error:', l2
      WRITE(*,250) 'LInf Error:', lInf
      250 FORMAT (A,T16,ES15.8)
   END IF

   ! Clean up the code.
   CALL SpectralFinalize()

END PROGRAM Differentiate_p

