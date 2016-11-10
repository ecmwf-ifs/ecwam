      SUBROUTINE OUTBETA (IJS, IJL, US, Z0, BETA)

! ----------------------------------------------------------------------

!**** *OUTBETA* - DETERMINES THE CHARNOCK PARAMETER

!     P.JANSSEN      KNMI/ECMWF  JANUARY 1992
!     J.BIDLOT       ECMWF       FEBRUARY 1996  MESSAGE PASSING
!     J.BIDLOT       ECMWF       AUGUST 2008  REMOVE MESSAGE PASSING 
!     J.BIDLOT       ECMWF       MARCH 2014  USE MODEL CHARNOCK FOR ALL
!                                SEA POINTS (INCLUDING SEA ICE POINTS)

!*    PURPOSE.
!     --------

!       COMPUTES THE CHARNOCK PARAMETER.

!**   INTERFACE.
!     ----------

!       *CALL* *OUTBETA (IJS, IJL, US, Z0, BETA)
!         *IJS*    - INDEX OF FIRST GRIDPOINT.
!         *IJL*    - INDEX OF LAST GRIDPOINT.
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

      USE YOWPCONS , ONLY : G        ,EPSUS
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL
      REAL, DIMENSION(IJS:IJL), INTENT(IN) :: US, Z0
      REAL, DIMENSION(IJS:IJL), INTENT(OUT) :: BETA

      INTEGER :: IJ

      REAL :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTBETA',0,ZHOOK_HANDLE)
#endif

!*    COMPUTE CHARNOCK 'CONSTANT' BETA.
!     ---------------------------------

      DO IJ = IJS,IJL
        BETA(IJ) = G*Z0(IJ)/MAX(US(IJ)**2,EPSUS)
      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTBETA',1,ZHOOK_HANDLE)
#endif

! ----------------------------------------------------------------------
      RETURN
      END SUBROUTINE OUTBETA
