!
!--------------------------------------------------------------------
!
      SUBROUTINE SECSPOM(F1,F3,IJS,IJL,NFRE,NANG,NMAX,NDEPTH,DEPTHA,
     V                   DEPTHD,OMSTART,FRAC,MR,DFDTH,OMEGA,DEPTH,
     V                   AKMEAN,TA,TB,TC_QL,TT_4M,TT_4P,IM_P,IM_M)
!
!--------------------------------------------------------------------
!
!*****SECSPOM** COMPUTES SECOND ORDER SPECTRUM IN FREQUENCY SPACE.
!
!     P.JANSSEN JULY 2008
!
!     PURPOSE
!     -------
!             DETERMINES SECOND-ORDER SPECTRUM, BASED ON JANSSEN (2008)
!             THERE ARE THREE CORRECTIONS:
!                   1) GENERATION OF SECOND-HARMONICS
!                   2) QUASI-LINEAR EFFECT
!                   3) SHIFT OF SPECTRUM BECAUSE OF STOKES FREQUENCY
!                      CORRECTION.
!
!     INTERFACE
!     ---------
!             *CALL* *SECSPOM(F1,F3,IJS,IJL,NFRE,NANG,NMAX,NDEPTH,
!                             DEPTHA,DEPTHD,FRAC,MR,DFDTH,OMEGA,
!                             DEPTH,AKMEAN,TA,TB,TC_QL,
!                             TT_4M,TT_4P,IM_P,IM_M)*
!
!
!     PARAMETER   TYPE      PURPOSE.
!     ---------   ----      -------
!
!       F1        REAL      2D FREE WAVE SPECTRUM (INPUT)
!       F3        REAL      BOUND WAVES SPECTRUM (OUTPUT)
!       IJS       INTEGER   FIRST GRID POINT
!       IJL       INTEGER   LAST GRID POINT
!       NFRE      INTEGER   NUMBER OF FREQUENCIES
!       NANG      INTEGER   NUMBER OF DIRECTIONS 
!       NMAX      INTEGER   MAXIMUM INDEX CORRESPONDS TO TWICE THE CUT-OFF
!                           FREQUENCY 
!       NDEPTH    INTEGER   NUMBER OF ENTRIES IN DEPTH TABLE
!       OMSTART   REAL      STARTING FREQUENCY OF ORIGINAL FREQUENCY GRID
!       FRAC      REAL      FRACTIONAL INCREASE IN FREQUENCY SPACE
!       MR        INTEGER   THINNING FACTOR IN FREQUENCY SPACE
!       DFDTH     REAL      PRODUCT OF INCREMENT IN FREQUENCY AND DIRECTION
!       OMEGA     REAL      ANGULAR FREQUENCY ARRAY
!       DEPTH     REAL      DEPTH ARRAY
!       AKMEAN    REAL      MEAN WAVENUMBER
!       TA        REAL      TABLE FOR MINUS INTERACTIONS
!       TB        REAL      TABLE FOR PLUS INTERACTIONS
!       TC_QL     REAL      TABLE FOR QUASI-LINEAR INTERACTIONS
!       TT_4M     REAL      TABLE FOR STOKES FREQUENCY CORRECTION
!       TT_4P     REAL      TABLE FOR STOKES FREQUENCY CORRECTION
!       IM_P      INTEGER   TABLE FOR WAVENUMBER M2 PLUS
!       IM_M      INTEGER   TABLE FOR WAVENUMBER M2 MIN
!
!
!
!     METHOD
!     ------
!             EVALUATE SECOND ORDER SPECTRUM IN FREQUENCY BASED ON
!             KRASITSKII'S CANONICAL TRANSFORMATION.
!
!     EXTERNALS
!     ---------
!             NONE
!
!     REFERENCES
!     ----------
!             V.E. ZAKHAROV, HAMILTONIAN APPROACH (1968) 
!             M.A. SROKOSZ, J.G.R.,91,995-1006 (1986)
!             P.A.E.M. JANSSEN, JFM (2009)
!
!
!--------------------------------------------------------------------
!
!
!     
      USE YOWPARAM, ONLY: CLDOMAIN
      USE YOWPCONS, ONLY: G, PI, ZPI
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------
   
      IMPLICIT NONE
   
      INTEGER IJ,IJS,IJL,NFRE,NANG,NDEPTH,M,K,M1,K1,M2,M2_M,M2_P,K2,
     V        MP,MM,L,MR,NMAX,ID
      INTEGER IM_P(NFRE,NFRE),IM_M(NFRE,NFRE),JD(IJS:IJL),
     V        INDEP(IJS:IJL)

      REAL OM0,OM0H,TH0,XK0,OM1,TH1,XK1,OM2,XK2,OM0P,XK0P,OM0M,XK0M,
     V     OMSTART,FRAC,XMR,XM2,XINCR1,XINCR2,XINCR3,XINCR4,FAC1,FAC2,
     V     FAC3,T_4M,T_4P,F2K(IJS:IJL),F2KP(IJS:IJL),F2KM(IJS:IJL),
     V     F2K1,F2K2,DELM1,XD,X_MIN,DEPTHA,DEPTHD
      REAL OMEGA(NFRE), DFDTH(NFRE)
      REAL TA(NDEPTH,NANG,NFRE,NFRE),TB(NDEPTH,NANG,NFRE,NFRE),
     V     TC_QL(NDEPTH,NANG,NFRE,NFRE),TT_4M(NDEPTH,NANG,NFRE,NFRE),
     V     TT_4P(NDEPTH,NANG,NFRE,NFRE)
      REAL F1(IJS:IJL,NANG,NFRE),F2(IJS:IJL,NANG,NMAX),
     V     F3(IJS:IJL,NANG,NFRE),DEPTH(IJS:IJL),AKMEAN(IJS:IJL),
     V     PSUM(IJS:IJL)

      REAL ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('SECSPOM',0,ZHOOK_HANDLE)
!
!***  1. COMPUTATION OF TAIL OF THE SPECTRUM AND INDEX JD
!     ---------------------------------------------------
!
!
      DO M=1,NFRE
         DO K=1,NANG
            DO IJ=IJS,IJL
               F3(IJ,K,M) = 0.
            ENDDO
         ENDDO
      ENDDO

      X_MIN = 1.0
      DO IJ=IJS,IJL
         XD = MAX(X_MIN/AKMEAN(IJ),DEPTH(IJ))
         XD = LOG(XD/DEPTHA)/LOG(DEPTHD)+1.
         ID = NINT(XD)
         ID = MAX(ID,1)
         JD(IJ) = MIN(ID,NDEPTH)
      ENDDO

      DO M=1,NFRE
         DO K=1,NANG
            DO IJ=IJS,IJL
               F2(IJ,K,M) = F1(IJ,K,M)
            ENDDO
         ENDDO
      ENDDO

      DO M=NFRE+1,NMAX
         OM0 = OMSTART*(1.+FRAC)**(MR*M-1)
         DO K=1,NANG
            DO IJ=IJS,IJL
               F2(IJ,K,M) = OMEGA(NFRE)**5*F1(IJ,K,NFRE)/OM0**5
            ENDDO
         ENDDO
      ENDDO
!     
!
!
!
!***  2. COMPUTATION OF THE 2nd ORDER FREQUENCY SPECTRUM.
!     ---------------------------------------------------
!
!
      DO M=1,NFRE
         OM0 = OMEGA(M)
         OM0H = OM0/2.                  
         MP   = MIN(M+1,NFRE)
         OM0P = OMEGA(MP)
         MM   = MAX(M-1,1) 
         OM0M = OMEGA(MM)
         DELM1 = 1./(OM0P-OM0M)
         DO K=1,NANG
            K2 = K
            DO IJ = IJS,IJL
               F2K(IJ)  = F2(IJ,K,M)
               F2KP(IJ) = F2(IJ,K,MP)
               F2KM(IJ) = F2(IJ,K,MM)
               PSUM(IJ) = 0.
            ENDDO               
            DO M1=1,NFRE
               OM1 = OMEGA(M1)
               M2_M = IM_M(M1,M)
               M2_P = IM_P(M1,M)
               DO K1=1,NANG
                  L = K-K1
                  IF (L.GT.NANG) L=L-NANG
                  IF (L.LT.1) L=L+NANG
    
                  DO IJ=IJS,IJL

                     F2K1 = F2(IJ,K1,M1)
!
!*                   2.1 OM0-OM1 CASE: SECOND HARMONICS
!                    ---------------------------------- 
!           
!                    OM2 = OM0-OM1
             
                     F2K2 = F2(IJ,K2,M2_M)
                     FAC1 = TA(JD(IJ),L,M1,M)
                     FAC2 = F2K1*F2K2+F2(IJ,K2,M1)*F2(IJ,K1,M2_M)

                     XINCR1 = FAC1*FAC2
                     PSUM(IJ)= PSUM(IJ)+XINCR1

!
!*                   2.2 OM1+OM0 CASE: INFRA-GRAVITY WAVES
!                    ------------------------------------- 
!           
!                    OM2 = OM1+OM0 

                     F2K2 = F2(IJ,K2,M2_P)
                     FAC3 = 2.*TB(JD(IJ),L,M1,M)
                     XINCR2 = FAC3*F2K2
!
!*                   2.3 QUASI-LINEAR EFFECT
!                    -----------------------
!           
!  
                     XINCR3 = TC_QL(JD(IJ),L,M1,M)*F2K(IJ)
!  
!*                   2.4 STOKES-FREQUENCY CORRECTION
!                    ------------------------------- 
!
!      
                     T_4M = TT_4M(JD(IJ),L,M1,M)
                     T_4P = TT_4P(JD(IJ),L,M1,M)
                     XINCR4 = -(F2KP(IJ)*T_4P-F2KM(IJ)*T_4M)*DELM1

                     PSUM(IJ) = PSUM(IJ)+F2K1*(XINCR2+XINCR3+XINCR4)
                     
                  ENDDO
               ENDDO
            ENDDO
            DO IJ=IJS,IJL
               F3(IJ,K,M) = F3(IJ,K,M)+PSUM(IJ)
            ENDDO
         ENDDO
      ENDDO 
! 
!
!--------------------------------------------------------------------
!      
      IF (LHOOK) CALL DR_HOOK('SECSPOM',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE SECSPOM
