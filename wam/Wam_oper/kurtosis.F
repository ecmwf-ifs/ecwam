      SUBROUTINE KURTOSIS(F3, IJS, IJL, C4, BF, QP)
 
!***  *KURTOSIS*   DETERMINES KURTOSIS, BENJAMIN-FEIR INDEX
!                  AND GODA'S PEAKEDNESS PARAMETER.
 
!     PETER JANSSEN       JANUARY 2003.
 
!     PURPOSE.
!     --------
 
!             DETERMINATION OF KURTOSIS, B-F INDEX AND GODA QF
 
!     INTERFACE.
!     ----------
!             *CALL*  *KURTOSIS((F3, IJS, IJL, C4, BF, QP)
 
!                      INPUT:
!                           *F3*    - 2-DIMENSIONAL SPECTRUM
!                           *IJS*   - STARTING INDEX
!                           *IJL*   - LAST INDEX
!                      OUTPUT: 
!                           *C4*    - KURTOSIS
!                           *BF*    - BENJAMIN-FEIR INDEX
!                           *QP*    - GODA'S QUALITY FACTOR
!
!     METHOD.
!     -------
!            STARTING POINT IS EQ.(29) FROM THE REFERENCE BELOW.
!            THIS EQUATION IS REWRITTEN FOR THE ANGULAR FREQUENCY-
!            DIRECTION SPECTRUM AND THE INTEGRATION OVER K4 IS PERFORMED
!            NEXT WE IMPOSE THE CONDITION FOR NARROW SPECTRA
!            IN FREQUENCY AND DIRECTION, WITH THE RESULT
!                 
!              C4 = (4*G*K_0^3)/(OM_0*M_0^2)INT F1D_1*F1D_2*F1D_3/DELOM
!
!             WHERE DELOM = OM_1+OM_2-OM_3-OM_4 AND 
!
!                  OM_4 = SQRT(|OM_1^2+OM_2^2-OM_3^2|)
!
!            FURTHERMORE, F1D IS THE ONE-DIMENSIONAL FREQUENCY SPECTRUM.
!            THE INTEGRAL IS A PRINCIPLE VALUE INTEGRAL AND IS OVER
!            THE FREQUENCIES OM_1, OM_2 AND OM_3. I HAVE PROGRAMMED
!            THE ABOVE EXPRESSION, BUT I FOUND THAT RESULTS DEPEND IN
!            A SENSITIVE MANNER ON RESOLUTION IN FREQUENCY SPACE.
!            FOR THE PRESENT OPERATIONAL RESOLUTION CONSIDERABLE
!            REDUCTIONS IN THE KURTOSIS ARE FOUND.
!
!            I THEREFORE PROCEED AS FOLLOWS. I ASSUME THAT THE SPECTRAL
!            PEAK HAS A GAUSSIAN SHAPE WITH STEEPNESS AND WIDTH GIVEN
!            BY THE ACTUAL SPECTRUM, AND A SYMMETRICAL, ACCURATE 
!            INTEGRATION OF THE PRINCIPLE VALUE INTEGRAL IS PERFORMED.
!            THIS NEEDS TO BE DONE ONLY ONCE, RESULTING IN A VERY
!            FAST ALGORITHM.
!            THE KURTOSIS NOW BECOMES

!                C4 = 2*BFI**2*J
!
!             WHERE 
 
!             J = INT DX1DX2DX3/(ZPI)^1.5*EXP(-0.5*(X1^2+X2^2+X3^2)/DEL
 
!             AND
 
!             DEL = (X3-X1)*(X3-X2)
 
!             THIS INTEGRAL CAN BE EVALUATED EXACTLY WITH THE RESULT:
 
!                J = PI/(6*SQRT(3))  
!
!     EXTERNALS.
!     ----------
!             NONE
!
!     REFERENCE:
!     ---------
 
!     NONLINEAR FOUR WAVE INTERACTIONS AND FREAK WAVES
!     PETER A.E.M. JANSSEN, JPO, 863-884, APRIL 2003
 
!-----------------------------------------------------------------------
!     MODULES:
!     --------

      USE YOWPCONS , ONLY : G        ,PI       ,ZPI      
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH   
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE

       
! ----------------------------------------------------------------------
 
      IMPLICIT NONE

      INTEGER :: IJ, M, K, IJS, IJL

      INTEGER, DIMENSION(IJS:IJL) ::  MMAX

      REAL :: ZSUM3, S_SIG, ZPISQRT, ZPI2OG
      REAL :: F3(NINF-1:NSUP,NANG,NFRE)
      REAL :: F1D(IJS:IJL,NFRE)
      REAL, DIMENSION(IJS:IJL) ::  C4,BF,QP,ETA2,SUM4,XKM,XMAX
      REAL :: FAC4(NFRE)
 
      ZSUM3  = ZPI/(6.*SQRT(3.))
      ZPISQRT = SQRT(ZPI)
      ZPI2OG = ZPI**2/G
 
!***  1. DETERMINE ONE-D SPECTRUM AND PEAK WAVE NUMBER.
!     ------------------------------------------------

      DO M=1,NFRE
         DO IJ=IJS,IJL   
            F1D(IJ,M) = 0.
         ENDDO   
      ENDDO
      DO M=1,NFRE
         DO K=1,NANG
            DO IJ=IJS,IJL
               F1D(IJ,M) = F1D(IJ,M)+F3(IJ,K,M)*DELTH
            ENDDO
         ENDDO
      ENDDO 
 
      DO IJ=IJS,IJL
        MMAX(IJ) = 1
        XMAX(IJ) = F1D(IJ,1)
      ENDDO
      DO M=2,NFRE
         DO IJ=IJS,IJL
            IF (F1D(IJ,M).GT.XMAX(IJ)) THEN
              MMAX(IJ) = M
              XMAX(IJ) = F1D(IJ,M)
            ENDIF
         ENDDO
      ENDDO   
      
      DO IJ=IJS,IJL
         XKM(IJ) = ZPI2OG*FR(MMAX(IJ))**2
      ENDDO   
 
!***  2. DETERMINE BF, QP AND C4.
!    ----------------------------
 
            
      DO M=1,NFRE
         FAC4(M) = 2.*FR(M)*DFIM(M)/DELTH
      ENDDO

      DO IJ=IJS,IJL
         ETA2(IJ)  = 0.
         SUM4(IJ)  = 0.
      ENDDO

      DO M=1,NFRE
         DO IJ=IJS,IJL
            ETA2(IJ) = ETA2(IJ)+F1D(IJ,M)*DFIM(M)/DELTH
            SUM4(IJ) = SUM4(IJ)+F1D(IJ,M)**2*FAC4(M)
         ENDDO
      ENDDO   
 
      DO IJ=IJS,IJL
         IF (ETA2(IJ).GT.0.) THEN
            QP(IJ) = SUM4(IJ)/ETA2(IJ)**2
            S_SIG = XKM(IJ)*SQRT(ETA2(IJ))
            BF(IJ) = ZPISQRT*S_SIG*QP(IJ)
            C4(IJ) = BF(IJ)**2*ZSUM3
         ELSE
            C4(IJ) = 0.
            BF(IJ) = 0.
            QP(IJ) = 0.
         ENDIF
      ENDDO
            
      RETURN

      END SUBROUTINE KURTOSIS
