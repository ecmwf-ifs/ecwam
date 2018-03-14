      SUBROUTINE PEAK_FREQ(F3,IJS,IJL,FP)
 
!***  *PEAK*   DETERMINES THE PEAK FREQUENCY OF THE SPECTRUM
 
!     PETER JANSSEN
 
!     PURPOSE.
!     --------
 
!              DETERMINATION OF PEAK PARAMETERS
 
!     INTERFACE.
!     ----------
!              *CALL*  *PEAK_FR(F1D,IJS,IJL,FP)*
 
!                       INPUT:
!                            *F3*   - FREQUENCY SPECTRUM
!                            *IJS*   - FIRST GRIDPOINT              
!                            *IJL*   - LAST GRIDPOINT  
 
!                       OUTPUT: 
!                            *FP*     - PEAK FREQUENCY
 
!     METHOD.
!     -------
!              NONE
 
!     EXTERNALS.
!     ----------
!              NONE
 
!-----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWPCONS , ONLY : G        ,PI       ,ZPI      
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,TH        ,   &
     &                 COSTH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
 
! ----------------------------------------------------------------------
 
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), INTENT(IN) :: F3(IJS:IJL,NANG,NFRE)
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: FP

      INTEGER(KIND=JWIM) :: IJ, M, K
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL) ::  MMAX

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: A,B,C,XP1,XM1,F1P1,F1M1,F10
      REAL(KIND=JWRB) :: F1D(IJS:IJL,NFRE)
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: XMAX

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('PEAK_FREQ',0,ZHOOK_HANDLE)
#endif

!***  1. DETERMINE ONE-DIMENSIONAL FREQUENCY SPECTRUM.
!     ------------------------------------------------
 
      DO M=1,NFRE
        K=1
        DO IJ=IJS,IJL
          F1D(IJ,M) = F3(IJ,K,M)*DELTH
        ENDDO
        DO K=2,NANG
          DO IJ=IJS,IJL
            F1D(IJ,M) = F1D(IJ,M)+F3(IJ,K,M)*DELTH
          ENDDO
        ENDDO
      ENDDO
 
!***  2. DETERMINE F_PEAK USING QUADRATIC FIT.
!     ---------------------------------------
 
!     MAX OF 1-D SPECTRUM
 
      DO IJ=IJS,IJL
        XMAX(IJ) = 0._JWRB
        MMAX(IJ) = NFRE-1
      ENDDO

      DO M=2,NFRE-1
        DO IJ=IJS,IJL
          IF (F1D(IJ,M).GT.XMAX(IJ)) THEN
            MMAX(IJ) = M
            XMAX(IJ) = F1D(IJ,M)
          ENDIF
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        FP(IJ) = FR(MMAX(IJ))
      ENDDO

!     DETERMINE QUADRATIC FIT.
      DO IJ=IJS,IJL
        IF (XMAX(IJ).GT.0._JWRB) THEN
          XP1 = FR(MMAX(IJ)+1)-FR(MMAX(IJ))
          XM1 = FR(MMAX(IJ)-1)-FR(MMAX(IJ))
          F10  = F1D(IJ,MMAX(IJ))
          F1P1 = F1D(IJ,MMAX(IJ)+1)-F10
          F1M1 = F1D(IJ,MMAX(IJ)-1)-F10

!!!          A = F10
          B = 1._JWRB/(XM1-XP1)*(XM1*F1P1/XP1-XP1*F1M1/XM1)
          C = 1._JWRB/(XM1-XP1)*(F1M1/XM1-F1P1/XP1)
        
          IF (C.LT.0._JWRB) THEN
            FP(IJ) = MIN(1._JWRB,FR(MMAX(IJ))-B/(2._JWRB*C))
!!!            XMAX(IJ) = A-B**2/(4._JWRB*C)
          ENDIF
        ENDIF
      ENDDO 
 
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('PEAK_FREQ',1,ZHOOK_HANDLE)
#endif
      END SUBROUTINE PEAK_FREQ
