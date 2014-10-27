!> @file IO.F90
!> @author Matthew Clay
!> @brief Routines for file i/o.
MODULE IO_m

   ! Required modules.
   USE Parameters_m,ONLY: IWP, RWP

   IMPLICIT NONE

   !> File name length for output files.
   INTEGER(KIND=IWP),PARAMETER,PUBLIC :: FILE_NAME_LENGTH = 256_IWP

   !> Interface for the file name function.
   INTERFACE FileName
      MODULE PROCEDURE FileNameWithoutNum
      MODULE PROCEDURE FileNameWithNum
   END INTERFACE FileName

   ! Module procedures.
   PUBLIC :: WriteData
   PRIVATE :: FileNameWithoutNum, FileNameWithNum

CONTAINS

   !> Routine to write the solution to file.
   !!
   !> @param[in] n Number of grid points in the domain.
   !> @param[in] x Locations of grid points.
   !> @param[in] u Solution to be written.
   !> @param[in] fname Name of the data file to be written.
   SUBROUTINE WriteData(n, x, u, fname)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWP),INTENT(IN) :: n
      REAL(KIND=RWP),DIMENSION(1:n),INTENT(IN) :: x, u
      CHARACTER(LEN=*),INTENT(IN) :: fname
      ! Local variables.
      ! Looping index.
      INTEGER(KIND=IWP) :: i

      ! Loop over the grid points and record the data to file.
      OPEN(UNIT=100,FILE=fname,FORM="FORMATTED",ACTION="WRITE",STATUS="REPLACE")
      DO i = 1, n
         WRITE(100,20) x(i), u(i)
         20 FORMAT (ES15.8,2X,ES15.8)
      END DO
      CLOSE(UNIT=100)
   END SUBROUTINE WriteData

   !> Routine to form a file name with a root name and file suffix.
   !!
   !> @param[in] root Root name for the file.
   !> @param[in] suffix Suffix for the file.
   !> @param[out] fname File name: root.suffix
   SUBROUTINE FileNameWithoutNum(root, suffix, fname)
      IMPLICIT NONE
      ! Calling arguments.
      CHARACTER(LEN=*),INTENT(IN) :: root, suffix
      CHARACTER(LEN=FILE_NAME_LENGTH),INTENT(OUT) :: fname

      ! Form the output file name.
      WRITE(fname,10) TRIM(root), '.', TRIM(suffix)
      10 FORMAT (A,A,A)
   END SUBROUTINE FileNameWithoutNum

   !> Routine to form a file name with a root name, number, and file suffix.
   !!
   !> @param[in] root Root name for the file.
   !> @param[in] num Number of the file.
   !> @param[in] suffix Suffix for the file.
   !> @param[out] fname File name: root_num.suffix
   SUBROUTINE FileNameWithNum(root, num, suffix, fname)
      IMPLICIT NONE
      ! Calling arguments.
      CHARACTER(LEN=*),INTENT(IN) :: root, suffix
      INTEGER(KIND=IWP),INTENT(IN) :: num
      CHARACTER(LEN=FILE_NAME_LENGTH),INTENT(OUT) :: fname

      ! Form the output file name.
      WRITE(fname,10) TRIM(root), '_', num, '.', TRIM(suffix)
      10 FORMAT (A,A,I4.4,A,A)
   END SUBROUTINE FileNameWithNum

END MODULE IO_m

