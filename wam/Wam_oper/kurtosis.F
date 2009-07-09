      SUBROUTINE KURTOSIS(F3, IJS, IJL, IG, C4, BF2, QP, FP, HMAX, TMAX)
 
!***  *KURTOSIS*   DETERMINES KURTOSIS, BENJAMIN-FEIR INDEX
!                  GODA'S PEAKEDNESS PARAMETER AND NORMALIZED
!                  MAXIMUM WAVE HEIGHT.
 
!     PETER JANSSEN       JULY 2007.
 
!     PURPOSE.
!     --------
 
!           DETERMINATION OF KURTOSIS, B-F INDEX , GODA QP AND
!                              HMAX AND TMAX
 
!     INTERFACE.
!     ----------
!           *CALL*  *KURTOSIS((F3, IJS, IJL, IG, C4, BF, QP, HMAX, TMAX)
 
!                      INPUT:
!                           *F3*    - 2-DIMENSIONAL SPECTRUM
!                           *IJS*   - STARTING INDEX
!                           *IJL*   - LAST INDEX
!                           *IG*    - BLOCK NUMBER
!                      OUTPUT: 
!                           *C4*    - KURTOSIS
!                           *BF*    - BENJAMIN-FEIR INDEX
!                           *QP*    - GODA'S QUALITY FACTOR
!                           *FP*    - PEAK FREQUENCY
!                           *HMAX*  - MAXIMUM WAVE HEIGHT
!                           *TMAX*  - MAXIMUM WAVE PERIOD
!
!     METHOD.
!     -------
!            USED NUMERICAL SIMULATIONS OF THE NONLINEAR SCHROEDINGER
!            EQUATION TO OBTAIN THE DEPENDENCE OF KURTOSIS ON THE
!            SQUARE OF THE BENJAMIN-FEIR INDEX, BF2, AND THE WIDTH OF
!            THE DIRECTIONAL DISTRIBUTION, SIG_TH. FITTING THE NUMERICAL 
!            RESULTS ONE FINDS (ONORATO, MORI AND JANSSEN, 2007)
!                 
!              C4 = A_MORI/SIG_TH * BF2 * PI/(3*SQRT(3))
!
!            WHERE A_MORI=0.031. NOTE THAT FOR NARROW SPECTRA, I.E.
!            SIG_TH => A_MORI, THE ONE-DIMENSIONAL RESULT, PREVIOUSLY
!            IMPLEMENTED, IS OBTAINED.
!
!            IN ADDITION, IT IS WELL-KNOWN THAT WHEN A WAVE GROUP
!            APPROACHES SHALLOW WATER THE NONLINEAR FOCUSSING 
!            EFFECT DISAPPEARS FOR K*D=1.363 AND WHEN K*D BECOMES
!            LESS THAN THIS CRITICAL VALUE EVEN DEFOCUSSING OCCURS.
!            AS A CONSEQUENCE, FOR K*D>1.363 C4 > 0 WHILE IN THE
!            OPPOSITE CASE C4 < 0. HERE, FOLLOWING THE WORK
!            OF JANSSEN & ONORATO (2007) WE HAVE PARAMETRIZED THIS SHALLOW
!            WATER EFFECT BY MEANS OF THE FUNCTION TRANSF2.
!
!
!     EXTERNALS.
!     ----------
!             AKI,TRANSF2,H_MAX,PEAK_FREQ
!
!     REFERENCES:
!     ---------
 
!     NONLINEAR FOUR WAVE INTERACTIONS AND FREAK WAVES
!     PETER A.E.M. JANSSEN, JPO, 863-884, APRIL 2003

!     ON KURTOSIS AND OCCURRENCE PROBABILITY OF FREAK WAVES
!     NOBUHITO MORI AND PETER A.E.M. JANSSEN, JPO, 1471-1483, JULY 2006 

!     THE INTERMEDIATE WATER DEPTH LIMIT OF THE ZAKHAROV EQUATION AND 
!     CONSEQUENCES FOR WAVE PREDICTION
!     PETER A.E.M. JANSSEN AND MIGUEL ONORATO, ACCEPTED FOR PUBLICATION
!     IN JPO, 2007

!     EFFECTS OF DIRECTIONALITY ON THE NONLINEAR FOCUSSING IN RANDOM SEAS. 
!     MIGUEL ONORATO, NOBUHITO MORI AND PETER A.E.M. JANSSEN
!     IN PREPARATION, 2007
!
!-----------------------------------------------------------------------
!     MODULES:
!     --------

      USE YOWPCONS , ONLY : G        ,PI       ,ZPI      ,ZPISQRT   ,
     &             BF2MAX  ,BF2MIN   ,C4MAX    ,C4MIN
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSHAL  , ONLY : DEPTH    ,TFAK     ,INDEP
       
! ----------------------------------------------------------------------
 
      IMPLICIT NONE

      INTEGER :: IJ, M, K, IJS, IJL, IG

      INTEGER, DIMENSION(IJS:IJL) ::  MMAX, KMAX, NW

      REAL :: ZSUM3, ZALPHA, ZBETA, ZA_MORI, ZJ , ZEPS, ZNU, 
     &        HS, SIG_OM1, DURATION, AKI, TRANSF2
      REAL :: F3(NINF-1:NSUP,NANG,NFRE)
      REAL :: F1D(IJS:IJL,NFRE), F1A(IJS:IJL,NANG) 
      REAL, DIMENSION(IJS:IJL) ::  C4,BF2,QP,HMAX,TMAX,ETA2,SUM0,SUM1,
     &                             SUM2,SUM4,XMAX,YMAX,XKP,FP,SIG_OM,
     &                             SIG_TH,EPS 
      REAL, DIMENSION(NFRE) :: DF,FAC1,FAC2,FAC3,FAC4
 
      ZSUM3  = PI/(3.*SQRT(3.))
      ZALPHA = 0.4
      ZA_MORI = 0.062
 
!***  1. DETERMINE ONE-D FREQUENCY AND ANGULAR SPECTRUM.
!     -------------------------------------------------

      DO M=1,NFRE
         DF(M) = DFIM(M)/DELTH
         FAC1(M) = FR(M)*DF(M)
         FAC2(M) = FR(M)**2*DF(M)
         FAC4(M) = 2.*FAC1(M)
      ENDDO

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

      DO K=1,NANG
         M=1
         DO IJ=IJS,IJL
            F1A(IJ,K) = F3(IJ,K,M)*DF(M)
         ENDDO
         DO M=2,NFRE
            DO IJ=IJS,IJL
               F1A(IJ,K) = F1A(IJ,K)+F3(IJ,K,M)*DF(M)
            ENDDO
         ENDDO
      ENDDO
!
!***  2. DETERMINE PEAK PARAMETERS FP,SIG_OM AND SIG_TH.
!     -------------------------------------------------
!
      CALL PEAK_FREQ(F1D,F1A,IJS,IJL,IG,XMAX,FP,SIG_OM,YMAX,SIG_TH)

 
!***  3. DETERMINE BF2, QP AND C4.
!    ----------------------------

      DO IJ=IJS,IJL
         ETA2(IJ)  = 0.
         SUM0(IJ)  = 0.
         SUM1(IJ)  = 0.
         SUM2(IJ)  = 0.
         SUM4(IJ)  = 0.
      ENDDO


      DO M=1,NFRE
         DO IJ=IJS,IJL
            ZBETA = ZALPHA*XMAX(IJ)
            IF (F1D(IJ,M).GT.ZBETA) THEN
               ETA2(IJ) = ETA2(IJ)+F1D(IJ,M)*DF(M)
               SUM4(IJ) = SUM4(IJ)+F1D(IJ,M)**2*FAC4(M)
            ENDIF 
            SUM0(IJ) = SUM0(IJ)+F1D(IJ,M)*DF(M)
            SUM1(IJ) = SUM1(IJ)+FAC1(M)*F1D(IJ,M)
            SUM2(IJ) = SUM2(IJ)+FAC2(M)*F1D(IJ,M)
         ENDDO
      ENDDO


      DO IJ=IJS,IJL
         IF (ETA2(IJ).GT.0.) THEN
            QP(IJ)     = SUM4(IJ)/ETA2(IJ)**2
            SIG_OM1    = 1./(QP(IJ)*ZPISQRT)
            SIG_OM(IJ) = MIN(SIG_OM(IJ),SIG_OM1)
            XKP(IJ)    = AKI(ZPI*FP(IJ),DEPTH(IJ,IG))
            EPS(IJ)    = XKP(IJ)*SQRT(SUM0(IJ))

            BF2(IJ) = 2.*EPS(IJ)**2/SIG_OM(IJ)**2*
     &                TRANSF2(XKP(IJ),DEPTH(IJ,IG))
            BF2(IJ) = MAX(MIN(BF2(IJ),BF2MAX),BF2MIN)
         ELSE
            BF2(IJ)    = 0.
            EPS(IJ)    = 0.
            SIG_OM(IJ) = 0.
            QP(IJ)     = 0.
         ENDIF 
      ENDDO
      DO IJ=IJS,IJL
         IF (ETA2(IJ).GT.0.) THEN
            SIG_TH(IJ) = MAX(SIG_TH(IJ),ZA_MORI)
            ZJ = ZA_MORI/SIG_TH(IJ)*ZSUM3
            C4(IJ) = ZJ*BF2(IJ)+8.*EPS(IJ)**2
            C4(IJ) = MAX(MIN(C4(IJ),C4MAX),C4MIN)
         ELSE
            C4(IJ)  = 0.
         ENDIF
      ENDDO  

! 
!***  5. DETERMINE HMAX AND TMAX. 
!       (EXPECTED IN RECORD OF LENGTH GIVEN BY DURATION)
!     -------------------------------------------------
!
      DURATION = 10800.

      DO IJ=IJS,IJL
         IF (FP(IJ).EQ.0.) THEN
            NW(IJ) = NINT(0.1*DURATION)
         ELSE
            NW(IJ) = NINT(DURATION*FP(IJ))
         ENDIF
      ENDDO

      CALL H_MAX(C4,NW,IJS,IJL,HMAX)

      DO IJ=IJS,IJL
         IF (SUM1(IJ).GT.0.) THEN
            ZNU = SQRT(SUM0(IJ)*SUM2(IJ)/SUM1(IJ)**2 - 1.)
            ZEPS = ZNU/(SQRT(2.)*HMAX(IJ))
            TMAX(IJ) = 1.+0.5*ZEPS**2+0.75*ZEPS**4
            TMAX(IJ) = TMAX(IJ)*SUM0(IJ)/SUM1(IJ)
         ELSE
            TMAX(IJ) = 0.
         ENDIF
      ENDDO

      DO IJ=IJS,IJL
         IF (SUM0(IJ).GT.0.) THEN
            HS = 4.*SQRT(SUM0(IJ))
            HMAX(IJ) = HMAX(IJ)*HS
         ELSE
            HMAX(IJ) = 0.
         ENDIF
      ENDDO 
            
      RETURN

      END SUBROUTINE KURTOSIS
