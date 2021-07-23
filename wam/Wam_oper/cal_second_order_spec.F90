SUBROUTINE CAL_SECOND_ORDER_SPEC(F1,KIJS,KIJL,DEPTH,SIG)
 
!***  *CAL_SEC_ORDER_SPEC*   DETERMINES SECOND_ORDER SPECTRUM
 
!     PETER JANSSEN JULY 2010
 
!     PURPOSE. 
!     --------
 
!              DETERMINATION OF SECOND-ORDER SPECTRUM
 
!     INTERFACE.
!     ----------
!              *CALL*  *CAL_SEC_ORDER_SPEC(F1,KIJS,KIJL,DEPTH,SIG)*
 
!                       INPUT:
!                            *F1*    - 2-D FREE WAVE SPECTRUM (at input)
!                            *KIJS*   - FIRST GRIDPOINT              
!                            *KIJL*   - LAST GRIDPOINT 
!                            *DEPTH* - DEPTH ARRAY (FROM KIJS:KIJL)
!                            *SIG*   - DIRECTION OF MAPPING:
!                                      FORWARD: SIG = +1
!                                      INVERSE: SIG = -1
 
!                       OUTPUT: 
!                            *F1*    - 2-D SPECTRUM INCLUDING SECOND-ORDER
!                                      CORRECTION (at output)
 
!     METHOD.
!     -------
!              IS DESCRIBED IN JANSSEN (2009), JFM, 637, 1-44.
 
!     EXTERNALS.
!     ----------
!              SECSPOM
 
!-----------------------------------------------------------------------
 
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED , ONLY : FR, DFIM, DELTH,TH,FRATIO
      USE YOWPARAM, ONLY : NANG, NFRE, CLDOMAIN
      USE YOWPCONS, ONLY : G, PI, ZPI
      USE YOWSHAL , ONLY : NDEPTH, DEPTHA, DEPTHD
      USE YOWTABL , ONLY : MR, XMR, MA, XMA, NFREH, NANGH, NMAX,         &
     &                     OMEGA, DFDTH, THH, DELTHH, IM_P, IM_M,        &
     &                     TA, TB, TC_QL, TT_4M, TT_4P
      USE YOWTEST , ONLY : IU06
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK

!-----------------------------------------------------------------------

      IMPLICIT NONE
#include "fkmean.intfb.h"
#include "secspom.intfb.h"

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT) :: F1
      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH
      REAL(KIND=JWRB), INTENT(IN) :: SIG

      INTEGER(KIND=JWIM) :: IJ,M,K,K0,M0,MP,KP,MM,KM,KL,KLL,ML

      REAL(KIND=JWRB) :: FRAC,CO1,DEL,DELF,D1,D2,D3,D4,C1
      REAL(KIND=JWRB) :: C2,XM,XK,OMSTART,AREA,SUM,SUM1,SUM3,GAM_B_J,ZFAC
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: EMEANALL, FMEANALL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: F1MEAN, AKMEAN, XKMEAN, EMAXL  
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE) :: F3
      REAL(KIND=JWRB), ALLOCATABLE ::  PF1(:,:,:),PF3(:,:,:)

!-----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('CAL_SECOND_ORDER_SPEC',0,ZHOOK_HANDLE)
 
!***  1. DETERMINE SECOND ORDER CORRECTION TO THE SPECTRUM
!     ----------------------------------------------------
 
      GAM_B_J = 0.6_JWRB
      ZFAC = GAM_B_J**2/16.0_JWRB
      FRAC = FRATIO-1.0_JWRB
      OMSTART = ZPI*FR(1)

      CALL FKMEAN(F1,KIJS,KIJL,KIJS,KIJL,EMEANALL,FMEANALL,F1MEAN,AKMEAN,XKMEAN)
 
!***  1.1 INTERPOLATION OR NOT.
!     ------------------------
 
      IF (MR.EQ.1 .AND. MA.EQ.1) THEN
 
!***     1.11 NO INTERPOLATION.
!        ----------------------
  
         CALL SECSPOM(F1(KIJS:KIJL,:,:),F3,KIJS,KIJL,NFRE,NANG,NMAX,NDEPTH,DEPTHA,       &
     &                DEPTHD,OMSTART,FRAC,MR,DFDTH,OMEGA,DEPTH,         &
     &                AKMEAN,TA,TB,TC_QL,TT_4M,TT_4P,IM_P,IM_M)
         DO M=1,NFRE
            DO K=1,NANG
               DO IJ=KIJS,KIJL
                  DELF = F3(IJ,K,M)
                  F1(IJ,K,M)=MAX(MIN(0.000001_JWRB,F1(IJ,K,M)),F1(IJ,K,M)+SIG*DELF)
               ENDDO
            ENDDO
         ENDDO

      ELSE

 
!***     1.12 ENERGY CONSERVING INTERPOLATION SCHEME
!        -------------------------------------------
 
         ALLOCATE(PF1(KIJS:KIJL,NANGH,NFREH))
         ALLOCATE(PF3(KIJS:KIJL,NANGH,NFREH))
         PF1(:,:,:) = 0._JWRB
         DO M=1,NFREH
            M0 = MR*M
            DO K=1,NANGH
               K0 = MA*K+1
               IF (K0.GT.NANG) K0 = K0-NANG
               IF (K0.LT.1) K0 = K0+NANG
               DO IJ=KIJS,KIJL
                  PF1(IJ,K,M)=PF1(IJ,K,M)+F1(IJ,K0,M0)
               ENDDO
            ENDDO
         ENDDO

!***     1.13 DETERMINE SECOND-ORDER SPEC
!        --------------------------------

         CALL SECSPOM(PF1,PF3,KIJS,KIJL,NFREH,NANGH,NMAX,NDEPTH,DEPTHA,   &
     &                DEPTHD,OMSTART,FRAC,MR,DFDTH,OMEGA,DEPTH,         &
     &                AKMEAN,TA,TB,TC_QL,TT_4M,TT_4P,IM_P,IM_M)
 
!***     2.24 INTERPOLATE TOWARDS HIGH-RES GRID
!        --------------------------------------
 
         DO IJ=KIJS,KIJL
           IF (EMEANALL(IJ).LE.ZFAC*DEPTH(IJ)**2) THEN 
             EMAXL(IJ)=1._JWRB
           ELSE
             EMAXL(IJ)=0._JWRB
           ENDIF
         ENDDO

         DO M=1,NFRE
            XM = REAL(M/MR)
!!!            D1 = REAL(M)/REAL(MR)-XM
            M0 = INT(XM)
            IF(M0.LT.1) THEN
              M0 = 1
              MP = M0+1
              D1 = 1._JWRB
            ELSE IF(M0.LT.NFREH) THEN
              MP = M0+1
              D1 = (FR(M)-FR(MR*M0))/(FR(MR*MP)-FR(MR*M0)) 
            ELSE
              M0 = NFREH
              MP = M0
              D1 = 0._JWRB
            ENDIF
            D2 = 1.0_JWRB-D1

            DO K=1,NANG
               XK = REAL((K-1)/MA)
               K0 = INT(XK)
               D3 = REAL(K-1)/REAL(MA)-XK
               D4 = 1.-D3
               IF (K0.LT.1) K0 = K0+NANGH
               KP = K0+1
               IF (KP.GT.NANGH) KP = KP-NANGH

               DO IJ=KIJS,KIJL
                  C1 = PF3(IJ,K0,M0)*D4+PF3(IJ,KP,M0)*D3
                  C2 = PF3(IJ,KP,MP)*D3+PF3(IJ,K0,MP)*D4
                  DELF = C1*D2+C2*D1

                  F1(IJ,K,M)=MAX(MIN(0.000001_JWRB,F1(IJ,K,M)),         &
     &                               F1(IJ,K,M)+EMAXL(IJ)*SIG*DELF)

               ENDDO
            ENDDO
         ENDDO
      ENDIF
      IF (ALLOCATED(PF1)) DEALLOCATE(PF1)
      IF (ALLOCATED(PF3)) DEALLOCATE(PF3)

      IF (LHOOK) CALL DR_HOOK('CAL_SECOND_ORDER_SPEC',1,ZHOOK_HANDLE)

END SUBROUTINE CAL_SECOND_ORDER_SPEC
