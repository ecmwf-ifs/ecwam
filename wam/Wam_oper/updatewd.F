      SUBROUTINE UPDATEWD (FL1, SL, IJS, IJL, IG, THWNEW, USNEW,
     &                     USOLD, TAUW, Z0NEW, U10NEW, ROAIRN, ZIDLNEW)

! ----------------------------------------------------------------------

!**** *UPDATEWD* - COMPUTATION OF WINDSPEED UPDATE AFTER ASSIMILATION

!     J. BIDLOT             ECMWF  APRIL 2001

!*    PURPOSE.
!     --------

!       UPDATE THE WIND SPEED AND WAVE STRESS USING NEW USTAR AND
!       SPRECTRUM.

!**   INTERFACE.
!     ----------

!       *CALL* *UPDATEWD (FL1, SL,IJS, IJL, IG, THWNEW, USNEW,
!                         TAUW, Z0NEW, U10NEW,ROAIRN, ZIDLNEW)*
!          *FL1* - SPECTRUM.
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY.
!          *IJS* - INDEX OF FIRST GRIDPOINT.
!          *IJL* - INDEX OF LAST GRIDPOINT.
!          *IG*  - ACTUAL BLOCK NUMBER.
!          *THWNEW*    WIND DIRECTION.
!          *USNEW*     NEW FRICTION VELOCITY IN M/S.
!          *TAUW*      WAVE STRESS IN (M/S)**2
!          *Z0NEW*     ROUGHNESS LENGTH IN M.
!          *U10NEW*    10m WIND SPEED.
!          *ROAIRN*    AIR DENSITY IN KG/M3.
!          *ZIDLNEW*   Zi/L (Zi: INVERSION HEIGHT,
!                             L: MONIN-OBUKHOV LENGTH).

!     METHOD.
!     -------

!     FIND THE NEW CHARNOCK PARAMETER BASED ON THE NEW USTAR AND TAUW
!     AND THUS UPDATE THE ROUGHNESS LENGTH. WITH THIS NEW ROUGHNESS
!     LENGTH AND THE NEW USTAR FIND THE NEW U10 USING THE LOGARITHMIC
!     WIND PROFILE.

!     EXTERNALS.
!     -----------

!       NONE.

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------

      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO
      USE YOWPCONS , ONLY : G
      USE YOWTABL  , ONLY : UMAX
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWIND  , ONLY : WSPMIN

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: IJ, IJS, IJL, IG, ILEV

      REAL,DIMENSION(NINF:NSUP,NBLO) :: USOLD, TAUW
      REAL,DIMENSION(NINF:NSUP) :: THWNEW, USNEW, Z0NEW, U10NEW,
     &                             ROAIRN,ZIDLNEW
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL1,SL

! ----------------------------------------------------------------------

      DO IJ=IJS,IJL
        U10NEW(IJ) = (USNEW(IJ)/USOLD(IJ,IG))*U10NEW(IJ)
        U10NEW(IJ) = MIN(U10NEW(IJ),UMAX)
        U10NEW(IJ) = MAX(U10NEW(IJ),WSPMIN)
      ENDDO

      ILEV=1

      CALL AIRSEA (U10NEW(IJS), TAUW(IJS,IG), USNEW(IJS), Z0NEW(IJS),
     &             IJS, IJL, ILEV)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. UPDATEWD: AIRSEA CALLED '
        CALL FLUSH (IU06)
      ENDIF

      CALL SINPUT (FL1, SL, IJS, IJL, THWNEW, USNEW, Z0NEW,
     &             ROAIRN, ZIDLNEW, SL)

!!!   SL is used twice which means that according to SINPUT it will
!!!   be at the end the wind input source function which is the part
!!!   needed by STRESSO 
!!!   Note that it is only possible because in SINPUT FL is determined
!!!   before SL.
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. UPDATEWD: SINPUT CALLED'
        CALL FLUSH (IU06)
      ENDIF

      CALL STRESSO (FL1,IJS,IJL,IG,THWNEW,USNEW,Z0NEW,
     &              ROAIRN, TAUW,SL)

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. UPDATEWD: STRESSO CALLED'
        CALL FLUSH (IU06)
      ENDIF

      RETURN
      END SUBROUTINE UPDATEWD
