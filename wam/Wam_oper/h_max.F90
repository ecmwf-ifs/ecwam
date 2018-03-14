      SUBROUTINE H_MAX(C3,C4,NSLC,IJS,IJL,AA,BB,HMAXN,SIG_HM)
!
!***  DETERMINE EXPECTED MAXIMUM WAVE HEIGHT, NORMALISED WITH
!     SIGNIFICANT WAVE HEIGHT.
!
!     PURPOSE:
!     -------
!
!     DETERMINE EXPECTED MAXIMUM ENVELOPE WAVE HEIGHT WHICH IS OBTAINED
!     FROM EXPECTED MAXIMUM WAVE ENERGY
!                 /
!         <E_MAX>=| DE E p_MAX(E)            
!                 /
!     WHERE P_MAX(E) IS DERIVED IN JANSSEN (2015) AND BASED ON WORK DONE BY
!     MONTINA ET AL (2009). MAXIMUM WAVE HEIGHT THEN FOLLOWS FROM
!     <H_MAX>=SQRT(<E_MAX>/2. 
!
!     VARIABLE       TYPE         PURPOSE
!     --------       ----         -------
!
!     C_3            REAL         SKEWNESS
!     C_4            REAL         KURTOSIS
!     NSLC           INTEGER      NUMBER OF SIGNIFICANT LEVEL CROSSINGS
!     H_C            REAL         FIRST GUESS OF EXPECTED MAX WAVE HEIGHT
!     IJS            INTEGER      FIRST INDEX
!     IJL            INTEGER      LAST INDEX
!     AA             REAL         FIRST PARAMETER OF RESIDORI PDF
!     BB             REAL         =2(AA+1)
!
!     OUTPUT:
!     ------
!
!     HMAXN          REAL         MAXIMUM WAVE HEIGHT NORMALIZED WITH
!                                 SIGNIFICANT WAVE HEIGHT H_S. 
!     SIG_HM         REAL         WIDTH OF HMAX DISTRIBUTION 
!
!     AUTHOR:
!     ------
!
!     P.A.E.M. JANSSEN, NOVEMBER 2017
!
!----------------------------------------------------------------------
!
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS , ONLY : PI
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

!----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL), INTENT(IN):: NSLC
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: C3, C4

      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: AA, BB, HMAXN, SIG_HM

      INTEGER(KIND=JWIM), PARAMETER :: NITER = 5
      INTEGER(KIND=JWIM) :: IJ, I

      REAL(KIND=JWRB), PARAMETER :: TWOSQRT6=2._JWRB*SQRT(6._JWRB)
      REAL(KIND=JWRB), PARAMETER :: GAM = 0.5772_JWRB
      REAL(KIND=JWRB), PARAMETER :: EB = 10._JWRB

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: TWOG1, G2, AE, BE, F, Z0, XNLOG
      REAL(KIND=JWRB), DIMENSION(IJS:IJL):: E, H_C, BBM1

!----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('H_MAX',0,ZHOOK_HANDLE)
#endif

!*    1. DETERMINE EXPECTED MAXIMUM WAVE HEIGHT.
!     -----------------------------------------
! 
      TWOG1 = -2._JWRB*GAM
      G2 = GAM**2+PI**2/6._JWRB
      AE = 0.5_JWRB*EB*(EB-2._JWRB)
      BE = 0.5_JWRB*EB*(EB**2-6._JWRB*EB+6._JWRB)

      DO IJ=IJS,IJL
        H_C(IJ) = 2._JWRB
        E(IJ)   = 2._JWRB*(IJ)**2
      ENDDO

      DO IJ=IJS,IJL
        IF (NSLC(IJ).GT.0) THEN
          F  = LOG(MAX(1._JWRB+C4(IJ)*AE+C3(IJ)**2*BE,0.1))

          AA(IJ) = ((EB-F)**2-2._JWRB*EB)/(2._JWRB*F)
          BB(IJ) = 2._JWRB*(1._JWRB+AA(IJ))
          BBM1(IJ) = 1._JWRB/BB(IJ)

          XNLOG = LOG(REAL(NSLC(IJ)))
          DO I=1,NITER
            Z0 = XNLOG+0.5_JWRB*LOG(E(IJ))
            E(IJ) = (G2-TWOG1*(AA(IJ)+Z0)+(2._JWRB*AA(IJ)+Z0)*Z0)*BBM1(IJ)
          ENDDO
          HMAXN(IJ) = SQRT(0.5_JWRB*E(IJ))
          SIG_HM(IJ) = PI*HMAXN(IJ)/(TWOSQRT6*(Z0+0.5_JWRB*GAM))
        ELSE
          AA(IJ) = 0._JWRB
          BB(IJ) = 2._JWRB*(1._JWRB+AA(IJ))
          BBM1(IJ) = 1._JWRB/BB(IJ)
          HMAXN(IJ) = 0._JWRB
          SIG_HM(IJ)= 0._JWRB
        ENDIF
      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('H_MAX',1,ZHOOK_HANDLE)
#endif
      END SUBROUTINE H_MAX
