!> @file Spectral.F90
!> @author Matthew Clay
!> @brief Module to compute the derivative of a signal using spectral methods.
!!
!! 1. We set the module up once for a code.
!! 2. Purpose is to calculate derivative of a physical signal.
!! 3. We assume the grid is 2*pi in length.
!! 4. In-place DFTs.
!! 5. Discuss inefficiencies due to array copying, but that we don't care.
MODULE Spectral_m

   ! Required modules.
   USE ISO_C_BINDING
   USE Parameters_m,ONLY: IWP, RWP

   IMPLICIT NONE

   ! FFTW procedure definitions.
   INCLUDE 'fftw3.f03'

   ! Important data from the calling code.
   !
   !> Number of grid points for local usage.
   INTEGER(KIND=IWP),PRIVATE :: n
   !> Number of grid points for usage with FFTW routines.
   INTEGER(KIND=C_INT),PRIVATE :: nFFTW
   !> Largest wavenumber index in spectral space.
   INTEGER(KIND=IWP),PRIVATE :: wavMax

   ! Variables for working with the FFTW library.
   !
   !> FFT plan in the FFTW library for computing the forward (R to C) FFT.
   TYPE(C_PTR),PRIVATE :: r2cPlan
   !> FFT plan in the FFTW library for computing the inverse (C to R) FFT.
   TYPE(C_PTR),PRIVATE :: c2rPlan
   !> Size of the FFT arrays for real data.
   INTEGER(KIND=IWP),PRIVATE :: rSize
   !> Size of the FFT arrays for complex data.
   INTEGER(KIND=IWP),PRIVATE :: cSize
   !> Main data array for working with the FFTW library.
   TYPE(C_PTR),PRIVATE :: data
   !> Real cast of the data array for working in physical space.
   REAL(C_DOUBLE),DIMENSION(:),POINTER,PRIVATE :: rData
   !> Complex cast of the data array for working in spectral space.
   COMPLEX(C_DOUBLE_COMPLEX),DIMENSION(:),POINTER,PRIVATE :: cData
   !> Complex imaginary unit.
   COMPLEX(C_DOUBLE_COMPLEX),PARAMETER,PRIVATE :: i1 = (0.0_RWP, 1.0_RWP)

   ! Module procedures.
   PUBLIC :: SpectralSetup, Differentiate, SpectralFinalize

CONTAINS

   !> Initialize the spectral module.
   !!
   !! In this procedure, we set up all of the working arrays and plans for
   !! interfacing with the FFTW library.
   !!
   !> @param[in] n_ Number of grid points in the simulation.
   SUBROUTINE SpectralSetup(n_)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWP),INTENT(IN) :: n_

      ! Variable initialization.
      n = n_
      nFFTW = INT(n, C_INT)
      wavMax = n/2_IWP

      ! The working arrays will be sized to store the complex data for the
      ! transform, since in-place DFTs will be performed.
      rSize = 2_IWP*(n/2_IWP + 1_IWP)
      cSize = n/2_IWP + 1_IWP
      !
      ! Allocate memory for the transforms via the FFTW library.
      data = FFTW_ALLOC_COMPLEX(INT(cSize, C_SIZE_T))
      CALL C_F_POINTER(data, rData, [rSize])
      CALL C_F_POINTER(data, cData, [cSize])
      cData(:) = (0.0_RWP, 0.0_RWP)

      ! Make the forward and reverse plans in the FFTW library.
      r2cPlan = FFTW_PLAN_DFT_R2C_1D(nFFTW, rData, cData, FFTW_ESTIMATE)
      c2rPlan = FFTW_PLAN_DFT_C2R_1D(nFFTW, cData, rData, FFTW_ESTIMATE)
   END SUBROUTINE SpectralSetup

   !> Differentiate a signal in physical space using spectral methods.
   !!
   !! Add a discussion about this.
   !!
   !> @param[in] n Size of the array.
   !> @param[in] u Signal in physical space at the grid points.
   !> @param[out] dudx Differentiated signal in physical space.
   SUBROUTINE Differentiate(n, u, dudx)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWP),INTENT(IN) :: n
      REAL(KIND=RWP),DIMENSION(n),INTENT(IN) :: u
      REAL(KIND=RWP),DIMENSION(n),INTENT(OUT) :: dudx
      ! Local variables.
      ! Index to access the complex data according to wavenumber.
      INTEGER(KIND=IWP) :: kInd
      ! Looping indices.
      INTEGER(KIND=IWP) :: i, k

      ! Zero out the FFT working arrays.
      cData(:) = (0.0_RWP, 0.0_RWP)

      ! Fill the forward FFT working array.
      DO i = 1, n
         rData(i) = u(i)
      END DO

      ! Perform the forward FFT.
      CALL FFTW_EXECUTE_DFT_R2C(r2cPlan, rData, cData)
      !
      ! Normalize the DFT.
      cData = cData/REAL(n, RWP)

      ! Loop over all wavenumbers and differentiate the signal.
      DO k = 0, wavMax
         kInd = k + 1_IWP
         cData(kInd) = i1*REAL(k, RWP)*cData(kInd)
      END DO
      !
      ! Avoid sawtooth errors by explicitly setting the N/2 mode to zero.
      cData(wavMax+1) = (0.0_RWP, 0.0_RWP)

      ! Now that differentiation has been performed in spectral space, inverse
      ! the transform to get the derivative in physical space.
      CALL FFTW_EXECUTE_DFT_C2R(c2rPlan, cData, rData)

      ! Fill in the output array for the calling code.
      DO i = 1, n
         dudx(i) = rData(i)
      END DO
   END SUBROUTINE Differentiate

   !> Procedure to finalize the spectral module.
   SUBROUTINE SpectralFinalize()
      IMPLICIT NONE

      ! Free the working arrays for the FFT.
      CALL FFTW_FREE(data)
      rData => NULL()
      cData => NULL()

      ! Clean up the FFTW plans.
      CALL FFTW_DESTROY_PLAN(r2cPlan)
      CALL FFTW_DESTROY_PLAN(c2rPlan)
   END SUBROUTINE SpectralFinalize

END MODULE Spectral_m

