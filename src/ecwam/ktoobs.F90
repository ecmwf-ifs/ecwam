! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE KTOOBS(IU06)

! ----------------------------------------------------------------------

!**** *KTOOBS* - ROUTINE TO SET KTOIS and KTOOBSTRUCT FOR THE DIFFERENT IPROPAGS

! ----------------------------------------------------------------------

USE YOWABORT , ONLY : WAM_ABORT
USE YOWUBUF  , ONLY : NPROPAGS, NANG_OBS, KTOIS, KTOOBSTRUCT

IMPLICIT NONE

INTEGER(KIND=JWIM), INTENT(IN) :: IU06

INTEGER(KIND=JWIM) :: IP, K

! ----------------------------------------------------------------------

DO IP = 0, NPROPAGS
  SELECT CASE(IP)

  CASE(0)

    DO K = 1, NANG_OBS
      SELECT CASE(K)
      CASE(1)
        KTOIS(K) = 1
        KTOOBSTRUCT(K) = 1
      CASE(2)
        KTOIS(K) = 0
        KTOOBSTRUCT(K) = 0
      CASE(3)
        KTOIS(K) = 1
        KTOOBSTRUCT(K) = 2
      CASE(4)
        KTOIS(K) = 0
        KTOOBSTRUCT(K) = 0
      CASE(5)
        KTOIS(K) = 2
        KTOOBSTRUCT(K) = 1
      CASE(6)
        KTOIS(K) = 0
        KTOOBSTRUCT(K) = 0
      CASE(7)
        KTOIS(K) = 2
        KTOOBSTRUCT(K) = 2
      CASE(8)
        KTOIS(K) = 0
        KTOOBSTRUCT(K) = 0
      CASE DEFAULT
        WRITE(IU06,*)'***ERROR IN KTOOBS: KTOIS and KTOOBSTRUCT not defined for that K !!! '
        CALL FLUSH(IU06)
        CALL WAM_ABORT(" KTOIS and KTOOBSTRUCT not defined for that K in ",__FILENAME__,__LINE__)
      END SELECT
    ENDDO

  CASE(1)

    DO K = 1, NANG_OBS
      SELECT CASE(K)
      CASE(1)
        KTOIS(K) = 1
        KTOOBSTRUCT(K) = 1
      CASE(2)
        KTOIS(K) = 1
        KTOOBSTRUCT(K) = 4
      CASE(3)
        KTOIS(K) = 1
        KTOOBSTRUCT(K) = 2
      CASE(4)
        KTOIS(K) = 2
        KTOOBSTRUCT(K) = 3
      CASE(5)
        KTOIS(K) = 2
        KTOOBSTRUCT(K) = 1
      CASE(6)
        KTOIS(K) = 2
        KTOOBSTRUCT(K) = 4
      CASE(7)
        KTOIS(K) = 2
        KTOOBSTRUCT(K) = 2
      CASE(8)
        KTOIS(K) = 1
        KTOOBSTRUCT(K) = 3
      CASE DEFAULT
        WRITE(IU06,*)'***ERROR IN KTOOBS: KTOIS and KTOOBSTRUCT not defined for that K !!! '
        CALL FLUSH(IU06)
        CALL WAM_ABORT(" KTOIS and KTOOBSTRUCT not defined for that K in ",__FILENAME__,__LINE__)
      END SELECT
    ENDDO

  CASE(2)

    DO K = 1, NANG_OBS
      SELECT CASE(K)
      CASE(1)
        KTOIS(K) = 1
        KTOOBSTRUCT(K) = 1
      CASE(2)
        KTOIS(K) = 1
        KTOOBSTRUCT(K) = 5
      CASE(3)
        KTOIS(K) = 1
        KTOOBSTRUCT(K) = 2
      CASE(4)
        KTOIS(K) = 2
        KTOOBSTRUCT(K) = 5
      CASE(5)
        KTOIS(K) = 2
        KTOOBSTRUCT(K) = 1
      CASE(6)
        KTOIS(K) = 3
        KTOOBSTRUCT(K) = 5
      CASE(7)
        KTOIS(K) = 2
        KTOOBSTRUCT(K) = 2
      CASE(8)
        KTOIS(K) = 4
        KTOOBSTRUCT(K) = 5
      CASE DEFAULT
        WRITE(IU06,*)'***ERROR IN KTOOBS: KTOIS and KTOOBSTRUCT not defined for that K !!! '
        CALL FLUSH(IU06)
        CALL WAM_ABORT(" KTOIS and KTOOBSTRUCT not defined for that K in ",__FILENAME__,__LINE__)
      END SELECT
    ENDDO


  CASE DEFAULT
    WRITE(IU06,*)'***ERROR IN KTOOBS: KTOIS and KTOOBSTRUCT not defined !!! '
    CALL FLUSH(IU06)
    CALL WAM_ABORT(" KTOIS and KTOOBSTRUCT not defined in ",__FILENAME__,__LINE__)
  END SELECT
ENDDO

END SUBROUTINE KTOOBS
