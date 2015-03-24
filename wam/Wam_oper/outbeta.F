      SUBROUTINE OUTBETA (IJS, IJL, IG, U10, US, Z0, BETA)

! ----------------------------------------------------------------------

!**** *OUTBETA* - DETERMINES THE CHARNOCK PARAMETER FOR THE IFS. 

!     P.JANSSEN      KNMI/ECMWF  JANUARY 1992
!     J.BIDLOT       ECMWF       FEBRUARY 1996  MESSAGE PASSING
!     J.BIDLOT       ECMWF       AUGUST 2008  REMOVE MESSAGE PASSING 
!     J.BIDLOT       ECMWF       MARCH 2014  USE MODEL CHARNOCK FOR ALL
!                                SEA POINTS (INCLUDING SEA ICE POINTS)

!*    PURPOSE.
!     --------

!       COMPUTES THE CHARNOCK PARAMETER FOR ATMOSPHERIC MODEL.

!**   INTERFACE.
!     ----------

!       *CALL* *OUTBETA (IJS, IJL, IG, IGL,US,Z0, BETA)
!         *IJS*    - INDEX OF FIRST GRIDPOINT.
!         *IJL*    - INDEX OF LAST GRIDPOINT.
!         *IG*     - BLOCK NUMBER.
!         *U10*    - 10M WIND SPEED IN M/S.
!         *US*     - FRICTION VELOCITY IN M/S.
!         *Z0*     - ROUGHNESS LENGTH IN M.
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
      USE YOWMPP   , ONLY : NINF     ,NSUP

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IJ, IJS, IJL
      INTEGER :: IG 

      REAL :: ZHOOK_HANDLE
      REAL :: BETA(IJS:IJL)
      REAL, DIMENSION(NINF:NSUP,NBLO) :: U10, US, Z0

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTBETA',0,ZHOOK_HANDLE)
#endif

!*    COMPUTE CHARNOCK 'CONSTANT' BETA.
!     ---------------------------------

      DO IJ = IJS,IJL
        BETA(IJ) = G*Z0(IJ,IG)/MAX(US(IJ,IG)**2,EPSUS)
      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTBETA',1,ZHOOK_HANDLE)
#endif

! ----------------------------------------------------------------------
      RETURN
      END SUBROUTINE OUTBETA
