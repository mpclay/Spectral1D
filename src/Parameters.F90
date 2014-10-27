!> @file Parameters.F90
!> @author Matthew Clay
!> @brief Parameters to define precision, etc.
MODULE Parameters_m

   ! Required modules.
   USE ISO_FORTRAN_ENV,ONLY: INT32, REAL64

   !> Working integer precision.
   INTEGER,PARAMETER,PUBLIC :: IWP = INT32
   !> Working real precision.
   INTEGER,PARAMETER,PUBLIC :: RWP = REAL64
   !> Pi.
   REAL(KIND=RWP),PARAMETER,PUBLIC :: PI = ACOS(-1.0_RWP)

END MODULE Parameters_m

