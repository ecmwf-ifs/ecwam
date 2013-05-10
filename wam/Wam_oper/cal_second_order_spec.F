!
!
!-----------------------------------------------------------------------
!
      SUBROUTINE CAL_SECOND_ORDER_SPEC(F1,F3,IJS,IJL,DEPTH,INDEP,SIG)
!
!***  *CAL_SEC_ORDER_SPEC*   DETERMINES SECOND_ORDER SPECTRUM
!
!     PETER JANSSEN JULY 2010
!
!     PURPOSE. 
!     --------
!
!              DETERMINATION OF SECOND-ORDER SPECTRUM
!
!     INTERFACE.
!     ----------
!              *CALL*  *CAL_SEC_ORDER_SPEC(F1,F3,IJS,IJL,DEPTH,INDEP,SIG)*
!
!                       INPUT:
!                            *F1*    - 2-D FREE WAVE SPECTRUM
!                            *IJS*   - FIRST GRIDPOINT              
!                            *IJL*   - LAST GRIDPOINT 
!                            *DEPTH* - DEPTH ARRAY (FROM IJS:IJL)
!                            *INDEP* - INDEX DEPTH ARRAY (FROM IJS:IJL)
!                            *SIG*   - DIRECTION OF MAPPING:
!                                      FORWARD: SIG = +1
!                                      INVERSE: SIG = -1
!
!                       OUTPUT: 
!                            *F1*    - 2-D SPECTRUM INCLUDING SECOND-ORDER
!                                      CORRECTION
!
!     METHOD.
!     -------
!              IS DESCRIBED IN JANSSEN (2009), JFM, 637, 1-44.
!
!     EXTERNALS.
!     ----------
!              SECSPOM
!
!-----------------------------------------------------------------------
!
      USE YOWFRED, ONLY: FR, DFIM, DELTH,TH,FRATIO
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    
      USE YOWPARAM, ONLY: NANG,NFRE,CLDOMAIN
      USE YOWPCONS, ONLY: G, PI, ZPI
      USE YOWSHAL, ONLY: NDEPTH,DEPTHA,DEPTHD
      USE YOWTABL, ONLY : MR, XMR, MA, XMA, NFREH, NANGH, NMAX, 
     &                    OMEGA, DFDTH, THH, DELTHH, IM_P, IM_M, 
     &                    TA, TB, TC_QL, TT_4M, TT_4P
      USE YOWTEST  , ONLY : IU06

      IMPLICIT NONE

      INTEGER IJ,IJS,IJL,M,K,K0,M0,MP,KP,MM,KM,KL,KLL,ML
      INTEGER INDEP(IJS:IJL)
      REAL,ALLOCATABLE ::  PF1(:,:,:),PF3(:,:,:)
      REAL, DIMENSION(IJS:IJL) :: DEPTH,F1MEAN,AKMEAN,XKMEAN 

      REAL SIG,FRAC,CO1,DEL,DELF,D1,D2,D3,D4,C1,
     &     C2,XM,XK,OMSTART,AREA,SUM,SUM1,SUM3,GAM_B_J,ZFAC

      REAL EMAXL(IJS:IJL)
      REAL F1(IJS:IJL,NANG,NFRE),F3(IJS:IJL,NANG,NFRE)
      REAL ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('CAL_SECOND_ORDER_SPEC',0,ZHOOK_HANDLE)
!
!***  1. DETERMINE SECOND ORDER CORRECTION TO THE SPECTRUM
!     ----------------------------------------------------
!
      GAM_B_J = 0.6
      ZFAC = GAM_B_J**2/16.
      FRAC = FRATIO-1.
      OMSTART = ZPI*FR(1)

      CALL FKMEAN(F1,IJS,IJL,EMEAN(IJS),FMEAN(IJS),F1MEAN,AKMEAN,
     &            XKMEAN)
!
!***  1.1 INTERPOLATION OR NOT.
!     ------------------------
!
      IF (MR.EQ.1 .AND. MA.EQ.1) THEN
!
!***     1.11 NO INTERPOLATION.
!        ----------------------
! 
         CALL SECSPOM(F1,F3,IJS,IJL,NFRE,NANG,NMAX,NDEPTH,DEPTHA,
     &                DEPTHD,OMSTART,FRAC,MR,DFDTH,OMEGA,DEPTH,
     &                AKMEAN,TA,TB,TC_QL,TT_4M,TT_4P,IM_P,IM_M)
         DO M=1,NFRE
            DO K=1,NANG
               DO IJ=IJS,IJL
                  DELF = F3(IJ,K,M)
                  F1(IJ,K,M)=MAX(MIN(0.000001,F1(IJ,K,M)),
     &                               F1(IJ,K,M)+SIG*DELF)
               ENDDO
            ENDDO
         ENDDO

      ELSE

!
!***     1.12 ENERGY CONSERVING INTERPOLATION SCHEME
!        -------------------------------------------
!
         ALLOCATE(PF1(IJS:IJL,NANGH,NFREH))
         ALLOCATE(PF3(IJS:IJL,NANGH,NFREH))
         PF1 = 0.
         DO M=1,NFREH
            DO K=1,NANGH
               DO IJ=IJS,IJL
                  M0 = MR*M
                  K0 = MA*K+1
                  IF (K0.GT.NANG) K0 = K0-NANG
                  IF (K0.LT.1) K0 = K0+NANG

                  PF1(IJ,K,M)=PF1(IJ,K,M)+F1(IJ,K0,M0)
               ENDDO
            ENDDO
         ENDDO

!
!***     1.13 DETERMINE SECOND-ORDER SPEC
!        --------------------------------
!

         CALL SECSPOM(PF1,PF3,IJS,IJL,NFREH,NANGH,NMAX,NDEPTH,DEPTHA,
     &                DEPTHD,OMSTART,FRAC,MR,DFDTH,OMEGA,DEPTH,
     &                AKMEAN,TA,TB,TC_QL,TT_4M,TT_4P,IM_P,IM_M)
!
!***     2.24 INTERPOLATE TOWARDS HIGH-RES GRID
!        --------------------------------------
!
         DO IJ=IJS,IJL
           IF (EMEAN(IJ).LE.ZFAC*DEPTH(IJ)**2) THEN 
             EMAXL(IJ)=1.
           ELSE
             EMAXL(IJ)=0.
           ENDIF
         ENDDO

         DO M=1,NFRE
            DO K=1,NANG
               DO IJ=IJS,IJL
               
                  XM = REAL(M/MR)
                  XK = REAL((K-1)/MA)

                  M0 = MAX(1,INT(XM))
                  K0 = INT(XK)

                  D1 = REAL(M)/REAL(MR)-XM
                  D2 = 1.-D1
                  D3 = REAL(K-1)/REAL(MA)-XK
                  D4 = 1.-D3

                  IF (K0.LT.1) K0 = K0+NANGH
                  MP = MIN(NFREH,M0+1)
                  KP = K0+1
                  IF (KP.GT.NANGH) KP = KP-NANGH

                  C1 = PF3(IJ,K0,M0)*D4+PF3(IJ,KP,M0)*D3
                  C2 = PF3(IJ,KP,MP)*D3+PF3(IJ,K0,MP)*D4

                  DELF = C1*D2+C2*D1

                  F1(IJ,K,M)=MAX(MIN(0.000001,F1(IJ,K,M)),
     &                               F1(IJ,K,M)+EMAXL(IJ)*SIG*DELF)

               ENDDO
            ENDDO
         ENDDO
      ENDIF
      IF (ALLOCATED(PF1)) DEALLOCATE(PF1)
      IF (ALLOCATED(PF3)) DEALLOCATE(PF3)
!      
      IF (LHOOK) CALL DR_HOOK('CAL_SECOND_ORDER_SPEC',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE CAL_SECOND_ORDER_SPEC
