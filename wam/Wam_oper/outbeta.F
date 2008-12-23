      SUBROUTINE OUTBETA (IJS, IJL, IG, US, Z0,
     &                    RCHAR, BETA)

! ----------------------------------------------------------------------

!**** *OUTBETA* - DETERMINES THE CHARNOCK PARAMETER FOR THE IFS. 

!     P.JANSSEN      KNMI/ECMWF  JANUARY 1992
!     J.BIDLOT       ECMWF       FEBRUARY 1996  MESSAGE PASSING
!     J.BIDLOT       ECMWF       AUGUST 2008  REMOVE MESSAGE PASSING 

!*    PURPOSE.
!     --------

!       COMPUTES THE CHARNOCK PARAMETER FOR ATMOSPHERIC MODEL.

!**   INTERFACE.
!     ----------

!       *CALL* *OUTBETA (IJS, IJL, IG, IGL,US,Z0,
!                        RCHAR, BETA)*
!         *IJS*    - INDEX OF FIRST GRIDPOINT.
!         *IJL*    - INDEX OF LAST GRIDPOINT.
!         *IG*     - BLOCK NUMBER.
!         *US*     - FRICTION VELOCITY IN M/S.
!         *Z0*     - ROUGHNESS LENGTH IN M.
!         *RCHAR*  - DEFAULT VALUES FOR CHARNOCK
!         *BETA*   - CHARNOCK FIELD (BLOCK ARRAY)
!         

!     EXTERNALS.
!     ----------

!       NONE.

!     METHOD.
!     -------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NBLO
      USE YOWPCONS , ONLY : G        ,EPSUS
      USE YOWICE   , ONLY : LICERUN  ,ICEMASK
      USE YOWMPP   , ONLY : NINF     ,NSUP

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IJS, IJL, IG
      INTEGER :: IJ
      REAL :: RCHAR
      REAL,DIMENSION(NINF:NSUP,NBLO) :: US, Z0
      REAL :: BETA(IJS:IJL)
      REAL ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('OUTBETA',0,ZHOOK_HANDLE)

!*    COMPUTE CHARNOCK 'CONSTANT' BETA.
!     ---------------------------------

      DO IJ = IJS,IJL
        IF (LICERUN.AND.ICEMASK(IJ,IG).EQ.0) THEN
          BETA(IJ) = RCHAR
        ELSE
          BETA(IJ) = G*Z0(IJ,IG)/MAX(US(IJ,IG)**2,EPSUS)
        ENDIF
      ENDDO

      IF (LHOOK) CALL DR_HOOK('OUTBETA',1,ZHOOK_HANDLE)

! ----------------------------------------------------------------------
      RETURN
      END SUBROUTINE OUTBETA
