      SUBROUTINE WVFRICVELO(IJS, IJL, U10, BETA, G, XKAPPA, ZREF, &
     &                      USTAR, Z0W) 

! ----------------------------------------------------------------------

!**** *WVFRICVELO* - COMPUTATION OF FRICTION VELOCITY OVER SEA


!*    PURPOSE.
!     ---------

!       TO COMPUTE THE FRICTION VELOCITY OVER SEA

!**   INTERFACE.
!     ----------

!       *CALL* *WVFRICVELO(IU06,ITEST)*
!          *IU06*  -  LOGICAL UNIT FOR PRINTER OUTPUT UNIT.

!     METHOD.
!     -------

!       A STEADY STATE NEUTRAL WIND PROFILE IS ASSUMED.
!       U10 = u* LOG(ZREF/Z0)
!       THE FRICTION VELOCITY IS COMPUTED USING THE CHARNOCK RELATION
!       FOR THE ROUGHNESS LENGTH Z0:

!       Z0 = BETA* U*^2 /g


!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!     USE HERSBACH 2011 FOR CD(U10) (SEE ALSO EDSON et al. 2013)


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), INTENT(IN) :: G, XKAPPA, ZREF
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: U10, BETA
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) ::  USTAR, Z0W

      INTEGER(KIND=JWIM), PARAMETER :: NITER=10

!     USE HERSBACH 2011 FOR CD(U10) (SEE ALSO EDSON et al. 2013)
      REAL(KIND=JWRB), PARAMETER :: C1 = 1.03E-3
      REAL(KIND=JWRB), PARAMETER :: C2 = 0.04E-3
      REAL(KIND=JWRB), PARAMETER :: P1 = 1.48
      REAL(KIND=JWRB), PARAMETER :: P2 = -0.21

      INTEGER(KIND=JWIM) :: IJ, ITER

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: C_D, UST, USTFG, Z0, F, DELFM1, FOLD
      REAL(KIND=JWRB) :: BOG, XKU10, XLOGZ, XLOGINV

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WVFRICVELO',0,ZHOOK_HANDLE)

      XLOGZ=LOG(ZREF)
      DO IJ=IJS,IJL
        BOG=BETA(IJ)/G
        XKU10=XKAPPA*U10(IJ)
        C_D = (C1 + C2*U10(IJ)**P1)*U10(IJ)**P2

        USTFG=SQRT(C_D)*U10(IJ)
        UST=USTFG
        FOLD=UST

        DO ITER=1,NITER
          Z0     = BOG*UST**2
          XLOGINV= 1._JWRB/(XLOGZ-LOG(Z0))
          F      = UST-XKU10*XLOGINV
!         protection in case it fails to converge
          IF(ABS(F) == ABS(FOLD) ) THEN
            EXIT
          ELSE IF (ABS(F) > ABS(FOLD) ) THEN
            UST=USTFG
            EXIT
          ENDIF
          DELFM1 = UST/(UST-2._JWRB*XKU10*XLOGINV**2)
          UST = UST-F*DELFM1
          FOLD = F
        ENDDO
        USTAR(IJ) = UST
        Z0W(IJ) = BOG*UST**2
      ENDDO

      IF (LHOOK) CALL DR_HOOK('WVFRICVELO',1,ZHOOK_HANDLE)

      END SUBROUTINE WVFRICVELO
