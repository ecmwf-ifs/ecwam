!
!--------------------------------------------------------------------
!
      SUBROUTINE SECSPOM(F1,F3,IJS,IJL,NFRE,NANG,NMAX,NDEPTH,DEPTHA,
     &                   DEPTHD,OMSTART,FRAC,MR,DFDTH,OMEGA,DEPTH,
     &                   AKMEAN,TA,TB,TC_QL,TT_4M,TT_4P,IM_P,IM_M)
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

      USE YOWPARAM, ONLY: CLDOMAIN
      USE YOWPCONS, ONLY: G, PI, ZPI
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------
   
      IMPLICIT NONE
   
      INTEGER, INTENT(IN) :: IJS, IJL, NFRE, NANG, NMAX, NDEPTH, MR
      INTEGER, DIMENSION(NFRE,NFRE), INTENT(IN) :: IM_P, IM_M 

      INTEGER :: IJ, M, K, M1, K1, M2_M, M2_P, K2, MP, MM, L, ID
      INTEGER, DIMENSION(NANG,NANG) :: LL
      INTEGER, DIMENSION(IJS:IJL) :: JD

      REAL, INTENT(IN) :: DEPTHA, DEPTHD, OMSTART, FRAC
      REAL, DIMENSION(NFRE), INTENT(IN) :: OMEGA, DFDTH
      REAL, DIMENSION(IJS:IJL), INTENT(IN) :: DEPTH, AKMEAN
      REAL, DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F1
      REAL, DIMENSION(IJS:IJL,NANG,NFRE), INTENT(OUT) :: F3 
      REAL, DIMENSION(NDEPTH,NANG,NFRE,NFRE), INTENT(IN) ::
     &                                  TA,TB,TC_QL,TT_4M,TT_4P

      REAL :: OM0, OM0P, OM0M, OM0H, OM1
      REAL :: T_4M, T_4P, DELM1, XD, X_MIN, OMRT, XLOGD, OMG5
      REAL :: ZHOOK_HANDLE
      REAL, DIMENSION(NMAX) :: OMEGA_EXT 
      REAL, DIMENSION(IJS:IJL) :: XINCR, DF2KP, DF2KM, PSUM
      REAL, DIMENSION(IJS:IJL,NANG,NMAX) :: F2

      LOGICAL :: LLSAMEDEPTH

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SECSPOM',0,ZHOOK_HANDLE)
#endif
!
!***  1. COMPUTATION OF TAIL OF THE SPECTRUM AND INDEX JD
!     ---------------------------------------------------
!
!
      X_MIN = 1.0
      XLOGD = LOG(DEPTHD)
      DO IJ=IJS,IJL
         XD = MAX(X_MIN/AKMEAN(IJ),DEPTH(IJ))
         XD = LOG(XD/DEPTHA)/XLOGD+1.
         ID = NINT(XD)
         ID = MAX(ID,1)
         JD(IJ) = MIN(ID,NDEPTH)
      ENDDO
      LLSAMEDEPTH=.TRUE.
      ID=JD(IJS)
      DO IJ=IJS+1,IJL
         LLSAMEDEPTH=(LLSAMEDEPTH.AND.(JD(IJ-1).EQ.JD(IJ)))
         IF(.NOT.LLSAMEDEPTH) EXIT
      ENDDO

      DO M=1,NFRE
         OMEGA_EXT(M)=OMEGA(M)
         DO K=1,NANG
            DO IJ=IJS,IJL
               F2(IJ,K,M) = F1(IJ,K,M)
            ENDDO
         ENDDO
      ENDDO

      OMG5 = OMEGA(NFRE)**5
      DO M=NFRE+1,NMAX
         OM0 = OMSTART*(1.+FRAC)**(MR*M-1)
         OMEGA_EXT(M)=OM0
         OMRT = OMG5/OM0**5
         DO K=1,NANG
            DO IJ=IJS,IJL
               F2(IJ,K,M) = OMRT*F1(IJ,K,NFRE)
            ENDDO
         ENDDO
      ENDDO

      DO K=1,NANG
        DO K1=1,NANG
          L = K-K1
          IF (L.GT.NANG) L=L-NANG
          IF (L.LT.1) L=L+NANG
          LL(K1,K)=L
        ENDDO
      ENDDO

!
!***  2. COMPUTATION OF THE 2nd ORDER FREQUENCY SPECTRUM.
!     ---------------------------------------------------
 
      DO M=1,NFRE
         OM0 = OMEGA(M)
         OM0H = 0.5*OM0      
!!!         MP   = MIN(M+1,NFRE)
         MP   = MIN(M+1,NMAX)
         OM0P = OMEGA_EXT(MP)
         MM   = MAX(M-1,1) 
         OM0M = OMEGA(MM)
         DELM1 = 1./(OM0P-OM0M)
         DO K=1,NANG
            K2 = K
            DO IJ = IJS,IJL
               DF2KP(IJ) = F2(IJ,K,MP)*DELM1
               DF2KM(IJ) = F2(IJ,K,MM)*DELM1
               PSUM(IJ) = 0.
            ENDDO               
            DO M1=1,NFRE
               OM1 = OMEGA(M1)
               M2_M = IM_M(M1,M)
               M2_P = IM_P(M1,M)
               DO K1=1,NANG
                  L = LL(K1,K) 

                  IF(LLSAMEDEPTH) THEN 
!*                  2.1 OM0-OM1 CASE: SECOND HARMONICS
!                   ---------------------------------- 
!                   OM2 = OM0-OM1
                    IF (ABS(OM1).LT.OM0H) THEN
!DIR$ PREFERVECTOR
                      DO IJ=IJS,IJL
                        PSUM(IJ)= PSUM(IJ)+TA(ID,L,M1,M)*
     &                                  ( F2(IJ,K1,M1)*F2(IJ,K2,M2_M)+
     &                                    F2(IJ,K2,M1)*F2(IJ,K1,M2_M) )
                      ENDDO
                    ENDIF

!*                  2.2 OM1+OM0 CASE: INFRA-GRAVITY WAVES
!                   ------------------------------------- 
!                   OM2 = OM1+OM0 
!DIR$ PREFERVECTOR
                    DO IJ=IJS,IJL
                      XINCR(IJ) = 2.*TB(ID,L,M1,M)*F2(IJ,K2,M2_P)
                    ENDDO
!
!*                  2.3 QUASI-LINEAR EFFECT
!                   -----------------------
!DIR$ PREFERVECTOR
                    DO IJ=IJS,IJL
                     XINCR(IJ)=XINCR(IJ)+TC_QL(ID,L,M1,M)*F2(IJ,K,M)
                    ENDDO
!  
!*                   2.4 STOKES-FREQUENCY CORRECTION
!                    ------------------------------- 
!DIR$ PREFERVECTOR
                    DO IJ=IJS,IJL
                       T_4M = TT_4M(ID,L,M1,M)
                       T_4P = TT_4P(ID,L,M1,M)
                       XINCR(IJ) = XINCR(IJ)
     &                          -(DF2KP(IJ)*T_4P-DF2KM(IJ)*T_4M)
                    ENDDO

                  ELSE

!*                  2.1 OM0-OM1 CASE: SECOND HARMONICS
!                   ---------------------------------- 
!                   OM2 = OM0-OM1
                    IF (ABS(OM1).LT.OM0H) THEN
!DIR$ PREFERVECTOR
                      DO IJ=IJS,IJL
                        PSUM(IJ)= PSUM(IJ)+TA(JD(IJ),L,M1,M)*
     &                                  ( F2(IJ,K1,M1)*F2(IJ,K2,M2_M)+
     &                                    F2(IJ,K2,M1)*F2(IJ,K1,M2_M) )
                      ENDDO
                    ENDIF

!*                  2.2 OM1+OM0 CASE: INFRA-GRAVITY WAVES
!                   ------------------------------------- 
!                   OM2 = OM1+OM0 
!DIR$ PREFERVECTOR
                    DO IJ=IJS,IJL
                      XINCR(IJ) = 2.*TB(JD(IJ),L,M1,M)*F2(IJ,K2,M2_P)
                    ENDDO
!
!*                  2.3 QUASI-LINEAR EFFECT
!                   -----------------------
!DIR$ PREFERVECTOR
                    DO IJ=IJS,IJL
                     XINCR(IJ)=XINCR(IJ)+TC_QL(JD(IJ),L,M1,M)*F2(IJ,K,M)
                    ENDDO
!  
!*                   2.4 STOKES-FREQUENCY CORRECTION
!                    ------------------------------- 
!DIR$ PREFERVECTOR
                    DO IJ=IJS,IJL
                       T_4M = TT_4M(JD(IJ),L,M1,M)
                       T_4P = TT_4P(JD(IJ),L,M1,M)
                       XINCR(IJ) = XINCR(IJ)
     &                          -(DF2KP(IJ)*T_4P-DF2KM(IJ)*T_4M)
                    ENDDO
                  ENDIF

                  DO IJ=IJS,IJL
                     PSUM(IJ) = PSUM(IJ)+F2(IJ,K1,M1)*XINCR(IJ)
                  ENDDO

               ENDDO
            ENDDO ! M1
            DO IJ=IJS,IJL
               F3(IJ,K,M) = PSUM(IJ)
            ENDDO
         ENDDO !K
      ENDDO !M 
 
!--------------------------------------------------------------------
       
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SECSPOM',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE SECSPOM
