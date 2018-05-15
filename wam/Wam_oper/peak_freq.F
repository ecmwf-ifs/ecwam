      SUBROUTINE PEAK_FREQ(F1D,F1A,IJS,IJL,                             &
     &              XMAX,FP,SIG_OM,YMAX,SIG_TH)
!
!***  *PEAK*   DETERMINES PARAMETERS AT PEAK OF SPECTRUM.
!
!     PETER JANSSEN
!
!     PURPOSE.
!     --------
!
!              DETERMINATION OF PEAK PARAMETERS.
!
!     INTERFACE.
!     ----------
!              *CALL*  *PEAK_FR(F1D,F1A,IJS,IJL,
!                               XMAX,FP,SIG_OM,YMAX,SIG_TH)*
!
!
!                       INPUT:
!                            *F1D*   - FREQUENCY SPECTRUM
!                            *F1A*   - DIRECTIONAL SPECTRUM
!                            *IJS*   - FIRST GRIDPOINT              
!                            *IJL*   - LAST GRIDPOINT  
!
!                       OUTPUT: 
!                            *XMAX*   - MAX FREQ SPECTRUM
!                            *FP*     - PEAK FREQUENCY
!                            *SIG_OM* - RELATIVE WIDTH IN FREQUENCY 
!                            *YMAX*   - MAX DIRECTIONAL SPECTRUM
!                            *SIG_TH* - RELATIVE WIDTH IN DIRECTION
!
!     METHOD.
!     -------
!              NONE
!
!     EXTERNALS.
!     ----------
!              NONE
!
!-----------------------------------------------------------------------
!
!     MODULES:
!     --------

      USE YOWPCONS , ONLY : G        ,PI       ,ZPI      
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,TH        ,
     &                 COSTH
      USE YOWPARAM , ONLY : NANG     ,NFRE
       
! ----------------------------------------------------------------------
 
      IMPLICIT NONE

      INTEGER :: IJ, M, K, IJS, IJL, KPLUS, KMIN

      INTEGER, DIMENSION(IJS:IJL) ::  MMAX, KMAX

      REAL :: ZEPSILON
      REAL :: A,B,C,XP1,XM1,F1P1,F1M1,F10,R1,SIG_TH2
      REAL :: F1D(IJS:IJL,NFRE), F1A(IJS:IJL,NANG) 
      REAL, DIMENSION(IJS:IJL) ::  SUM0,SUM1,SUM2,SUM3,SUM4,
     &                             XMAX,YMAX,XKP,FP,SIG_OM,SIG_TH,
     &                             SIG_TH1,THMAX

!
!***  0. SECURITY FIRST

      ZEPSILON=100.*EPSILON(ZEPSILON)

!
!***  1. DETERMINE F_PEAK USING QUADRATIC FIT.
!     ---------------------------------------

!     MAX OF 1-D SPECTRUM

      DO IJ=IJS,IJL
         XMAX(IJ) = 0.
         MMAX(IJ) = 1
      ENDDO

      DO M=1,NFRE
         DO IJ=IJS,IJL
            IF (F1D(IJ,M).GT.XMAX(IJ)) THEN
               MMAX(IJ) = M
               XMAX(IJ) = F1D(IJ,M)
            ENDIF
         ENDDO
      ENDDO

      DO IJ=IJS,IJL
         FP(IJ) = FR(MMAX(IJ))
         SIG_OM(IJ)=1.0      
      ENDDO

      DO IJ=IJS,IJL 
         IF (XMAX(IJ).GT.0..AND.MMAX(IJ).LT.NFRE.AND.MMAX(IJ).GT.1) THEN
            XP1 = FR(MMAX(IJ)+1)-FR(MMAX(IJ))
            XM1 = FR(MMAX(IJ)-1)-FR(MMAX(IJ))
            F10  = F1D(IJ,MMAX(IJ)) 
            F1P1 = F1D(IJ,MMAX(IJ)+1)-F10
            F1M1 = F1D(IJ,MMAX(IJ)-1)-F10

            A = F10
            B = 1./(XM1-XP1)*(XM1*F1P1/XP1-XP1*F1M1/XM1)
            C = 1./(XM1-XP1)*(F1M1/XM1-F1P1/XP1)

            IF (C.LT.0.) THEN
               FP(IJ) = FR(MMAX(IJ))-B/(2.*C)
               XMAX(IJ) = A-B**2/(4.*C)
               SIG_OM(IJ) = SQRT((B**2-4.*A*C)/MAX(8.*C**2,ZEPSILON))
     &                      / FP(IJ)
            ENDIF

         ENDIF
      ENDDO 

!***  2. DETERMINE YMAX, THMAX AND SIG_TH  USING QUADRATIC FIT.
!     --------------------------------------------------------

      DO IJ=IJS,IJL
         YMAX(IJ) = 0.
         KMAX(IJ) = 1
      ENDDO

      DO K=1,NANG
         DO IJ=IJS,IJL
            IF (F1A(IJ,K).GT.YMAX(IJ)) THEN
               KMAX(IJ) = K
               YMAX(IJ) = F1A(IJ,K)
            ENDIF
         ENDDO
      ENDDO 

      DO IJ=IJS,IJL
         THMAX(IJ) = TH(KMAX(IJ)) 
         SIG_TH1(IJ) = 1. 
      ENDDO 

      DO IJ=IJS,IJL
         IF (YMAX(IJ).GT.0.)THEN
            XP1 = DELTH
            XM1 = -DELTH
            KPLUS = KMAX(IJ)+1
            IF (KPLUS.GT.NANG) KPLUS = KPLUS-NANG
            KMIN = KMAX(IJ)-1
            IF (KMIN.LT.1) KMIN = KMIN+NANG
            F10  = F1A(IJ,KMAX(IJ)) 
            F1P1 = F1A(IJ,KPLUS)-F10
            F1M1 = F1A(IJ,KMIN)-F10

            A = F10
            B = 1./(XM1-XP1)*(XM1*F1P1/XP1-XP1*F1M1/XM1)
            C = 1./(XM1-XP1)*(F1M1/XM1-F1P1/XP1)
         
            IF (C.LT.0.) THEN
               THMAX(IJ) = THMAX(IJ)-B/(2.*C)
               YMAX(IJ) = A-B**2/(4.*C)
               SIG_TH1(IJ) = SQRT((B**2-4.*A*C)/MAX(8.*C**2,ZEPSILON))
            ENDIF
         ENDIF 
      ENDDO

      DO IJ=IJS,IJL
         SUM1(IJ) = 0.
         SUM2(IJ) = 0.
      ENDDO

      DO K=1,NANG
         DO IJ=IJS,IJL
            IF (COS(TH(K)-THMAX(IJ)).GE.0.5) THEN
              SUM1(IJ) = SUM1(IJ) +F1A(IJ,K)*DELTH
              SUM2(IJ) = SUM2(IJ) +COS(TH(K)-THMAX(IJ))*
     &                   F1A(IJ,K)*DELTH
            ENDIF
         ENDDO
      ENDDO
      DO IJ=IJS,IJL
         IF (SUM1(IJ).GT.0.) THEN
            R1 = SUM2(IJ)/SUM1(IJ)
            SIG_TH2 = SQRT(2.*MAX(0., 1.-R1))
         ELSE
            SIG_TH2 = 1.
         ENDIF
         SIG_TH(IJ) = MIN(SIG_TH1(IJ),SIG_TH2)
      ENDDO

      RETURN

      END SUBROUTINE PEAK_FREQ
