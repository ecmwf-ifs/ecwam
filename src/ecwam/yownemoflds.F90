! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWNEMOFLDS

      USE PARKIND_WAVE, ONLY : JWIM
      USE YOWDRVTYPE  , ONLY : WAVE2OCEAN, OCEAN2WAVE

      IMPLICIT NONE

!*     ** *NEMOFLDS* NEMO FIELDS FOR COUPLED RUNS

      TYPE(WAVE2OCEAN) :: WAM2NEMO
      TYPE(OCEAN2WAVE) :: NEMO2WAM

      LOGICAL :: LNEMOCITHICK, LNEMOICEREST

!     OPTIONAL NEMO -> WAM RECORDS (KIND JWRO) AFTER THE NREAL RECORDS OF THE
!     LAW RESTART FILE: SST, ICE COVER, ICE THICKNESS, U AND V CURRENTS.
      INTEGER(KIND=JWIM), PARAMETER :: NNEMOREAL=5
      LOGICAL :: LREST_NEMO_ICE = .FALSE.

!     LAW FILE WRITTEN BY SAVSTRESS AND WAITING FOR ITS NEMO RECORDS, WHICH
!     RECVNEMOFIELDS APPENDS AFTER THE NEXT RECEIVE (PER-TASK FILES ONLY).
      LOGICAL :: LLAWPEND = .FALSE.
      INTEGER(KIND=JWIM) :: IJINFLAW, IJSUPLAW
      CHARACTER(LEN=296) :: CFILELAW

!--------------------------------------------------------------------

!*    VARIABLE     TYPE      PURPOSE
!     --------     ----      -------
!     LNEMOCITHICK LOGICAL   SET TO TRUE IF SEA ICE THICKNESS IS 
!                            AVAILABLE FROM NEMO (E.G. LIM2 ACTIVE).
!     LNEMOICEREST LOGICAL   SET TO TRUE IF SEA ICE IS NOT RESCALED BY ICE COVER
!     NNEMOREAL    INTEGER   NUMBER OF OPTIONAL NEMO RECORDS IN THE LAW FILE.
!     LREST_NEMO_ICE LOGICAL TRUE IF THE LAW FILE READ AT RESTART HAD THEM.
!     LLAWPEND     LOGICAL   TRUE IF CFILELAW IS WAITING FOR THEM.
!     IJINFLAW, IJSUPLAW     POINT RANGE OF THE RECORDS IN CFILELAW.
!     CFILELAW     CHARACTER NAME OF THE LAW FILE WAITING FOR THEM.
!---------------------------------------------------------------------
      END MODULE YOWNEMOFLDS
