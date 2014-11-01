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
!> @file Advection.F90
!> @author Matthew Clay
!> @brief Execution code for the scalar advection equation.
!!
!! In this code we solve the scalar advection equation using spectral method for
!! spatial differentiation. The variables are evolved in physical space, but
!! differentiation is carried out in spectral space. In the spectral module, the
!! FFT is performed using the FFTW library. Once differentiation is performed in
!! spectral space, the inverse transform is used to get the derivative in
!! physical space.
PROGRAM Advection_p

   ! Required modules.
   USE Parameters_m,ONLY: IWP, RWP, PI
   USE IO_m,ONLY: FILE_NAME_LENGTH, FileName, WriteData
   USE Spectral_m,ONLY: SpectralSetup, SpectralFinalize

   IMPLICIT NONE

   ! The different initial conditions that can be used.
   !
   !> Sinusoudal.
   INTEGER(KIND=IWP),PARAMETER :: SINUSOIDAL = 0_IWP
   !> Exponential function and sin/cos combinations.
   INTEGER(KIND=IWP),PARAMETER :: EXPONENTIAL = 1_IWP
   !
   ! The different wave speeds that can be used.
   !
   ! NOTE: only constant wave speed is currently implemented.
   !
   !> Constant wavespeed of 2*PI.
   INTEGER(KIND=IWP),PARAMETER :: CONSTANT = 0_IWP
   !> Variable wavespeed from "Spectral Methods in MATLAB" page 24.
   INTEGER(KIND=IWP),PARAMETER :: VARIABLE = 1_IWP

   ! User inputs for the simulation.
   !
   !> Number of grid points. Must be an even number.
   INTEGER(KIND=IWP),PARAMETER :: n = 32_IWP
   !> Time step for temporal integration.
   REAL(KIND=RWP),PARAMETER :: dt = 0.00001_RWP
   !> Number of steps to take.
   INTEGER(KIND=IWP),PARAMETER :: nEnd = 100000_IWP
   !> Type of initial conditions.
   INTEGER(KIND=IWP),PARAMETER :: ics = EXPONENTIAL
   !> Type of wavespeed for the system.
   INTEGER(KIND=IWP),PARAMETER :: aType = CONSTANT
   !> Number of steps to write a file.
   INTEGER(KIND=IWP),PARAMETER :: writePeriod = 10000_IWP
   !> Number of steps to print information to the user.
   INTEGER(KIND=IWP),PARAMETER :: printPeriod = 10000_IWP
   !> Whether or not to report error analysis at the end of the simulation.
   LOGICAL,PARAMETER :: reportErrors = .TRUE.

   ! Variables required for the simulation.
   !
   !> The computational grid.
   REAL(KIND=RWP),DIMENSION(n) :: x
   !> The local grid spacing.
   REAL(KIND=RWP) :: dx
   !> Solution at the start and end of the time step.
   REAL(KIND=RWP),DIMENSION(n) :: u0
   !> Solution for the intermediate RK steps.
   REAL(KIND=RWP),DIMENSION(n) :: ui
   !> Array for the RHS of the ODE system.
   REAL(KIND=RWP),DIMENSION(n) :: Lu
   !> Array to store the initial condition for error checking.
   REAL(KIND=RWP),DIMENSION(n) :: ic
   !> The current simulation time.
   REAL(KIND=RWP) :: t
   !> The current simulation step.
   INTEGER(KIND=IWP) :: nadv

   ! Extraneous variables.
   !
   !> Variable for output file names.
   CHARACTER(LEN=FILE_NAME_LENGTH) :: fname
   !> Looping index.
   INTEGER(KIND=IWP) :: i
   !> Variables for the order of convergence checking.
   REAL(KIND=RWP) :: l2, lInf, diff

   ! Initialize the spectral module.
   CALL SpectralSetup(n)

   ! Fill in the computational grid.
   DO i = 1, n
      x(i) = REAL(i, RWP)*2.0_RWP*PI/REAL(n, RWP)
   END DO
   dx = 2.0_RWP*PI/REAL(n, RWP)

   ! Form the initial conditions.
   SELECT CASE (ics)
      CASE (SINUSOIDAL)
         DO i = 1, n
            u0(i) = SIN(2.0_RWP*PI*REAL(i, RWP)/REAL(n, RWP))
         END DO
      CASE (EXPONENTIAL)
         DO i = 1, n
            u0(i) = EXP(EXP(SIN(x(i))**2*COS(x(i))))*EXP(COS(x(i))**2*SIN(x(i)))
         END DO
      CASE DEFAULT
         WRITE(*,100) 'Invalid IC option. Defaulting to sinusoidal.'
         100 FORMAT (A)
         u0(i) = SIN(2.0_RWP*PI*REAL(i, RWP)/REAL(n, RWP))
   END SELECT
   !
   ! Store the initial condition for error checking.
   ic(:) = u0(:)

   ! Print some information to the user.
   !
   ! ... to be added.

   ! Enter the main time loop, in which the TVD RK3 scheme of Shu is used.
   nadv = 0_IWP
   t = 0.0_RWP
   tloop: DO WHILE (nadv < nEnd)
      ! Check if we need to write data to file.
      IF (MOD(nadv, writePeriod) == 0_IWP) THEN
         CALL FileName('Solution', nadv, 'dat', fname)
         CALL WriteData(n, x, u0, fname)
      END IF

      ! Check if we need to print information to the screen.
      IF (MOD(nadv, printPeriod) == 0_IWP) THEN
         WRITE(*,200) 'Simulation step number: ', nadv, &
                      '; Simulation time: ', t, &
                      '; Max: ', MAXVAL(u0), &
                      '; Min: ', MINVAL(u0)
         200 FORMAT (A,I10.10,A,ES15.8,A,ES15.8,A,ES15.8)
      END IF

      ! Calculate the RHS with the data at the start of the time step.
      CALL RHS(n, x, u0, aType, Lu)
      !
      ! Determine the first intermediate stage in the RK scheme.
      ui(:) = u0(:) + dt*Lu(:)

      ! Calculate the RHS using the first intermediate stage.
      CALL RHS(n, x, ui, aType, Lu)
      !
      ! Determine the second intermediate stage in the RK scheme.
      ui(:) = 0.75_RWP*u0(:) + 0.25_RWP*ui(:) + dt*0.25_RWP*Lu(:)

      ! Calculate the RHS using the second intermediate stage.
      CALL RHS(n, x, ui, aType, Lu)
      !
      ! Determine u at the next time step.
      u0(:) = u0(:)/3.0_RWP + &
              2.0_RWP*ui(:)/3.0_RWP + &
              dt*2.0_RWP*Lu(:)/3.0_RWP

      ! Update the simulation time and step counters.
      nadv = nadv + 1_IWP
      t = t + dt
   END DO tloop
   !
   ! Write out the final state to the user.
   WRITE(*,200) 'Simulation step number: ', nadv, &
                '; Simulation time: ', t, &
                '; Max: ', MAXVAL(u0), &
                '; Min: ', MINVAL(u0)
   
   ! Calculate the errors.
   IF (reportErrors) THEN
      l2 = 0.0_RWP
      lInf = 0.0_RWP
      DO i = 1, n
         diff = u0(i) - ic(i)
         l2 = l2 + diff**2
         IF (ABS(diff) > lInf) THEN
            lInf = ABS(diff)
         END IF
      END DO
      l2 = SQRT(l2*dx)
      WRITE(*,100) ''
      WRITE(*,110) 'L2 Error:', l2
      WRITE(*,110) 'LInf Error:', lInf
      110 FORMAT (A,T16,ES15.8)
   END IF

   ! Write out the final solution.
   CALL FileName('Solution', nadv, 'dat', fname)
   CALL WriteData(n, x, u0, fname)

   ! Clean up the code.
   CALL SpectralFinalize()

CONTAINS

   !> Subroutine to calculate the RHS of the ODE system.
   !!
   !! The RHS of the scalar advection equation is -c(x)*dudx, where c(x) is the
   !! wave speed. We evaluate the derivative dudx using the spectral module, in
   !! which the derivative is evaluated in spectral space.
   !!
   !> @param[in] n Number of grid points.
   !> @param[in] x Physical grid.
   !> @param[in] u Solution at the current RK stage.
   !> @param[in] aType The type of wave speed to use.
   !> @param[in] Lu The RHS of the ODE system.
   SUBROUTINE RHS(n, x, u, aType, Lu)
      ! Required modules.
      USE Spectral_m,ONLY: Differentiate
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWP),INTENT(IN) :: n, aType
      REAL(KIND=RWP),DIMENSION(n),INTENT(IN) :: x, u
      REAL(KIND=RWP),DIMENSION(n),INTENT(OUT) :: Lu
      ! Local variables.
      ! Derivative of the incoming signal.
      REAL(KIND=RWP),DIMENSION(n) :: dudx
      ! Looping index.
      INTEGER(KIND=IWP) :: i

      ! Calculate the derivative of the incoming signal.
      CALL Differentiate(n, u, dudx)

      ! Form the RHS of the ODE system.
      SELECT CASE (aType)
         CASE (CONSTANT)
            Lu(:) = -2.0_RWP*PI*dudx(:)
         CASE (VARIABLE)
            DO i = 1, n
               Lu(i) = -dudx(i)
            END DO
         CASE DEFAULT
            WRITE(*,100) 'Invalid wave speed type. Halting.'
            100 FORMAT (A)
            STOP
      END SELECT
   END SUBROUTINE RHS

END PROGRAM Advection_p

