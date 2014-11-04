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

   !> Interface for writing out the data.
   INTERFACE WriteData
      MODULE PROCEDURE WriteDataIntReal
      MODULE PROCEDURE WriteDataRealReal
   END INTERFACE WriteData

   ! Module procedures.
   PRIVATE :: FileNameWithoutNum, FileNameWithNum
   PRIVATE :: WriteDataIntReal, WriteDataRealReal

CONTAINS

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
      10 FORMAT (A,A,I8.8,A,A)
   END SUBROUTINE FileNameWithNum

   !> Routine to write integer and real data arrays to file.
   !!
   !> @param[in] n Size of the arrays.
   !> @param[in] i1 First array containing integer data.
   !> @param[in] x1 Second array containing real data.
   !> @param[in] fname Name of the data file to be written.
   SUBROUTINE WriteDataIntReal(n, i1, x1, fname)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWP),INTENT(IN) :: n
      INTEGER(KIND=IWP),DIMENSION(1:n),INTENT(IN) :: i1
      REAL(KIND=RWP),DIMENSION(1:n),INTENT(IN) :: x1
      CHARACTER(LEN=*),INTENT(IN) :: fname
      ! Local variables.
      ! Looping index.
      INTEGER(KIND=IWP) :: i

      ! Loop over the grid points and record the data to file.
      OPEN(UNIT=100,FILE=fname,FORM="FORMATTED",ACTION="WRITE",STATUS="REPLACE")
      DO i = 1, n
         WRITE(100,20) i1(i), x1(i)
         20 FORMAT (I10.10,2X,ES15.8)
      END DO
      CLOSE(UNIT=100)
   END SUBROUTINE WriteDataIntReal

   !> Routine to write two real arrays to file.
   !!
   !> @param[in] n Size of the arrays.
   !> @param[in] x1 First array of real data.
   !> @param[in] x2 Second array of real data.
   !> @param[in] fname Name of the data file to be written.
   SUBROUTINE WriteDataRealReal(n, x1, x2, fname)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWP),INTENT(IN) :: n
      REAL(KIND=RWP),DIMENSION(1:n),INTENT(IN) :: x1, x2
      CHARACTER(LEN=*),INTENT(IN) :: fname
      ! Local variables.
      ! Looping index.
      INTEGER(KIND=IWP) :: i

      ! Loop over the grid points and record the data to file.
      OPEN(UNIT=100,FILE=fname,FORM="FORMATTED",ACTION="WRITE",STATUS="REPLACE")
      DO i = 1, n
         WRITE(100,20) x1(i), x2(i)
         20 FORMAT (ES15.8,2X,ES15.8)
      END DO
      CLOSE(UNIT=100)
   END SUBROUTINE WriteDataRealReal

END MODULE IO_m

