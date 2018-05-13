      SUBROUTINE OUTWINT

! ----------------------------------------------------------------------

!**** *OUTWINT* -

!*    PURPOSE.
!     --------


!**   INTERFACE.
!     ----------

!     METHOD.
!     -------

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : JPPFLAG ,ITOBOUT, NIPRMOUT , NINFOBOUT,     &
     &                      INFOBOUT,BOUT
      USE YOWGRID  , ONLY : IJSLOC   ,IJLLOC
      USE YOWTEST  , ONLY : IU06
      USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('OUTWINT',0,ZHOOK_HANDLE)

!!!! if I/O server pass BOUT with INFOBOUT to it

!!!  BOUT(IJSLOC:IJLLOC,NIPRMOUT) contains the local contributions of the NIPRMOUT integrated parameters
!!!  that will need to be output
!!!! NIPRMOUT is not always the same as the model could output less at step 0
!!!  These parameters are defined with INFOBOUT(NIPRMOUT,NINFOBOUT), where
!    INFOBOUT(:,1)  : GRIB TABLE NUMBER.
!    INFOBOUT(:,2)  : GRIB PARAMETER IDENTIFIER.
!    INFOBOUT(:,3)  : GRIB REFERENCE LEVEL IN FULL METER.



!!!  else not not the i/o server

       CALL OUTINT

!!!  endif


      IF (LHOOK) CALL DR_HOOK('OUTWINT',1,ZHOOK_HANDLE)

      END SUBROUTINE OUTWINT
