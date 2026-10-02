! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWNEMOFLDS

      USE PARKIND_WAVE, ONLY : JWIM, JWRB
      USE YOWDRVTYPE  , ONLY : WAVE2OCEAN, OCEAN2WAVE

      IMPLICIT NONE

!*     ** *NEMOFLDS* NEMO FIELDS FOR COUPLED RUNS

      TYPE(WAVE2OCEAN) :: WAM2NEMO
      TYPE(OCEAN2WAVE) :: NEMO2WAM

      LOGICAL :: LNEMOCITHICK, LNEMOICEREST

!     LAW RESTART FILE WAITING FOR ITS NEMO -> WAM FIELDS (17-21). SAVSTRESS
!     WRITES THE FILE AND KEEPS A COPY HERE; THE NEXT RECVNEMOFIELDS REWRITES
!     IT WITH THE FIELDS THAT RECEIVE BRINGS, WHICH ARE THE ONES A RESTARTED
!     RUN NEEDS FOR ITS FIRST FORCING UPDATE (PER-TASK FILES ONLY).
      LOGICAL :: LLAWPEND = .FALSE.
      INTEGER(KIND=JWIM) :: IJINFLAW, IJSUPLAW
      CHARACTER(LEN=296) :: CFILELAW
      CHARACTER(LEN=14) :: CDLAWHDR(4)
      REAL(KIND=JWRB), ALLOCATABLE :: RFIELDLAW(:,:)

!--------------------------------------------------------------------

!*    VARIABLE     TYPE      PURPOSE
!     --------     ----      -------
!     LNEMOCITHICK LOGICAL   SET TO TRUE IF SEA ICE THICKNESS IS 
!                            AVAILABLE FROM NEMO (E.G. LIM2 ACTIVE).
!     LNEMOICEREST LOGICAL   SET TO TRUE IF SEA ICE IS NOT RESCALED BY ICE COVER
!---------------------------------------------------------------------
      END MODULE YOWNEMOFLDS
