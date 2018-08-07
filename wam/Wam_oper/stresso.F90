      SUBROUTINE STRESSO (F, SPOS, IJS, IJL,                            &
     &                    MIJ, RHOWGDFTH,                               &
     &                    THWNEW, USNEW, Z0NEW, ROAIRN,                 &
     &                    TAUW, PHIWA)

! ----------------------------------------------------------------------

!**** *STRESSO* - COMPUTATION OF WAVE STRESS.

!     H. GUNTHER      GKSS/ECMWF NOVEMBER   1989 CODE MOVED FROM SINPUT.
!     P.A.E.M. JANSSEN     KNMI  AUGUST     1990
!     J. BIDLOT            ECMWF FEBRUARY   1996-97
!     S. ABDALLA           ECMWF OCTOBER    2001 INTRODUCTION OF VARIABLE
!                                                AIR DENSITY
!     J. BIDLOT            ECMWF            2007  ADD MIJ
!     P.A.E.M. JANSSEN     ECMWF            2011  ADD FLUX CALULATIONS

!*    PURPOSE.
!     --------

!       COMPUTE NORMALIZED WAVE STRESS FROM INPUT SOURCE FUNCTION

!**   INTERFACE.
!     ----------

!        *CALL* *STRESSO (F, SPOS, IJS, IJL,
!    &                    MIJ, RHOWGDFTH,
!    &                    THWNEW, USNEW, Z0NEW, ROAIRN,
!    &                    TAUW, PHIWA)*
!         *F*           - WAVE SPECTRUM.
!         *SPOS*        - POSITIVE WIND INPUT SOURCE FUNCTION ARRAY.
!         *IJS*         - INDEX OF FIRST GRIDPOINT.
!         *IJL*         - INDEX OF LAST GRIDPOINT.
!         *MIJ*         - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!         *RHOWGDFTH    - WATER DENSITY * G * DF * DTHETA
!                         FOR TRAPEZOIDAL INTEGRATION BETWEEN FR(1) and FR(MIJ) 
!                         !!!!!!!!  RHOWGDFTH=0 FOR FR > FR(MIJ)
!         *THWNEW*      - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                         NOTATION (POINTING ANGLE OF WIND VECTOR,
!                         CLOCKWISE FROM NORTH).
!         *USNEW*       - NEW FRICTION VELOCITY IN M/S.
!         *Z0NEW*       - ROUGHNESS LENGTH IN M.
!         *ROAIRN*      - AIR DENSITY IN KG/M**3.
!         *TAUW*        - KINEMATIC WAVE STRESS IN (M/S)**2
!         *PHIWA*       - ENERGY FLUX FROM WIND INTO WAVES INTEGRATED
!                         OVER THE FULL FREQUENCY RANGE.

!     METHOD.
!     -------

!       THE INPUT SOURCE FUNCTION IS INTEGRATED OVER FREQUENCY
!       AND DIRECTIONS.
!       BECAUSE ARRAY *SPOS* IS USED, ONLY THE INPUT SOURCE
!       HAS TO BE STORED IN *SPOS* (CALL FIRST SINPUT, THEN
!       STRESSO, AND THEN THE REST OF THE SOURCE FUNCTIONS)

!     EXTERNALS.
!     -----------

!       NONE.

!     REFERENCE.
!     ----------
!       R SNYDER ET AL,1981.
!       G. KOMEN, S. HASSELMANN AND K. HASSELMANN, JPO, 1984.
!       P. JANSSEN, JPO, 1985

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : ALPHA    ,TAUWSHELTER, ITSHELT,LWVFLX_SNL
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,TH       ,    &
     &            COSTH    ,SINTH    ,FR5      ,FRATIO
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI
      USE YOWSHAL  , ONLY : CINV     ,INDEP
      USE YOWTABL  , ONLY : EPS1
      USE YOWSTAT  , ONLY : IPHYS

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "tau_phi_hf.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL), INTENT(IN) :: MIJ

      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F, SPOS
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE), INTENT(IN) :: RHOWGDFTH
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: THWNEW, USNEW, Z0NEW, ROAIRN
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: TAUW, PHIWA


      INTEGER(KIND=JWIM) :: IJ, M, K, I, J, II

      REAL(KIND=JWRB) :: PHIHF_REDUC
      REAL(KIND=JWRB) :: CONST
      REAL(KIND=JWRB) :: XI, XJ, DELI1, DELI2, DELJ1, DELJ2, XK, DELK1, DELK2
      REAL(KIND=JWRB) :: PHI2
      REAL(KIND=JWRB) :: ABS_TAUWSHELTER, GM1
      REAL(KIND=JWRB) :: COSW, FCOSW2
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: TAUHF, TEMP1, CONST1, XSTRESS, YSTRESS
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: PHIHF, TEMP2, CONST2
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: TAU1, PHI1, XLEVTAIL 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: CMRHOWGDFTH
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: US2, TAUX, TAUY, TAUPX, TAUPY
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: USDIRP, UST, UST2 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: SUMT, SUMX, SUMY

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('STRESSO',0,ZHOOK_HANDLE)

!*    1. PRECOMPUTE FREQUENCY SCALING.
!        -----------------------------

      GM1 = 1.0_JWRB/G
      CONST = DELTH*(ZPI)**4*GM1

!     REDUCTION OF FOR HIGH FREQUENCY ENERGY FLUX INTO OCEAN CONTRIBUTING TO UPPER OCEAN MIXING
      IF(.NOT.LWVFLX_SNL) THEN
        PHIHF_REDUC=1.0_JWRB
      ELSE
        PHIHF_REDUC=1.0_JWRB
      ENDIF

      ABS_TAUWSHELTER=ABS(TAUWSHELTER)

      DO IJ=IJS,IJL
        CONST1(IJ)  = CONST*FR5(MIJ(IJ))*GM1
        CONST2(IJ)  = ROAIRN(IJ)*CONST*FR5(MIJ(IJ))
      ENDDO

!*    2. COMPUTE WAVE STRESS OF ACTUAL BLOCK.
!        ------------------------------------

!*    2.2 INTEGRATE INPUT SOURCE FUNCTION OVER FREQUENCY AND DIRECTIONS.
!         --------------------------------------------------------------

      DO IJ=IJS,IJL
        PHIWA(IJ)   = 0.0_JWRB
        XSTRESS(IJ) = 0.0_JWRB
        YSTRESS(IJ) = 0.0_JWRB
      ENDDO
      DO M=1,MAXVAL(MIJ(:))
!     THE INTEGRATION ONLY UP TO FR=MIJ SINCE RHOWGDFTH=0 FOR FR>MIJ
        K=1
        DO IJ=IJS,IJL
          SUMT(IJ) = SPOS(IJ,K,M)
          SUMX(IJ) = SPOS(IJ,K,M)*SINTH(K)
          SUMY(IJ) = SPOS(IJ,K,M)*COSTH(K)
        ENDDO
        DO K=2,NANG
          DO IJ=IJS,IJL
            SUMT(IJ) = SUMT(IJ) + SPOS(IJ,K,M)
            SUMX(IJ) = SUMX(IJ) + SPOS(IJ,K,M)*SINTH(K)
            SUMY(IJ) = SUMY(IJ) + SPOS(IJ,K,M)*COSTH(K)
          ENDDO
        ENDDO
        DO IJ=IJS,IJL
          PHIWA(IJ)   =  PHIWA(IJ) + SUMT(IJ)*RHOWGDFTH(IJ,M)
          CMRHOWGDFTH(IJ) = CINV(INDEP(IJ),M)*RHOWGDFTH(IJ,M)
          XSTRESS(IJ) = XSTRESS(IJ) + SUMX(IJ)*CMRHOWGDFTH(IJ)
          YSTRESS(IJ) = YSTRESS(IJ) + SUMY(IJ)*CMRHOWGDFTH(IJ)
        ENDDO
      ENDDO

!     TAUW is the kinematic wave stress !
      DO IJ=IJS,IJL
        XSTRESS(IJ) = XSTRESS(IJ)/MAX(ROAIRN(IJ),1.0_JWRB)
        YSTRESS(IJ) = YSTRESS(IJ)/MAX(ROAIRN(IJ),1.0_JWRB)
      ENDDO

!*    2.3 CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
!     ----------------------------------------------------

      DO IJ=IJS,IJL
        US2(IJ)=USNEW(IJ)**2
      ENDDO

      IF (IPHYS.EQ.0 .OR. ITSHELT.EQ.0) THEN
        DO IJ=IJS,IJL
          USDIRP(IJ)=THWNEW(IJ)
          UST(IJ)=USNEW(IJ)
          UST2(IJ)=US2(IJ)
        ENDDO
      ELSE
        DO IJ=IJS,IJL
          TAUX(IJ)=US2(IJ)*SIN(THWNEW(IJ))
          TAUY(IJ)=US2(IJ)*COS(THWNEW(IJ))
          TAUPX(IJ)=TAUX(IJ)-ABS_TAUWSHELTER*XSTRESS(IJ)
          TAUPY(IJ)=TAUY(IJ)-ABS_TAUWSHELTER*YSTRESS(IJ)
          USDIRP(IJ)=ATAN2(TAUPX(IJ),TAUPY(IJ))
          UST2(IJ)=SQRT(TAUPX(IJ)**2+TAUPY(IJ)**2)
          UST(IJ)=SQRT(UST2(IJ))
        ENDDO
      ENDIF

      K=1
      DO IJ=IJS,IJL
        COSW     = MAX(COS(TH(K)-THWNEW(IJ)),0.0_JWRB)
        FCOSW2   = F(IJ,K,MIJ(IJ))*COSW**2
        TEMP1(IJ) = FCOSW2*COSW
        TEMP2(IJ) = FCOSW2 
      ENDDO

      DO K=2,NANG
        DO IJ=IJS,IJL
          COSW     = MAX(COS(TH(K)-THWNEW(IJ)),0.0_JWRB)
          FCOSW2   = F(IJ,K,MIJ(IJ))*COSW**2
          TEMP1(IJ) = TEMP1(IJ) + FCOSW2*COSW
          TEMP2(IJ) = TEMP2(IJ) + FCOSW2 
        ENDDO
      ENDDO

      IF (ITSHELT.EQ.0) THEN
        DO IJ=IJS,IJL
          XLEVTAIL(IJ)=0.0_JWRB
        ENDDO
      ELSE
        DO IJ=IJS,IJL
          XLEVTAIL(IJ)=TAUWSHELTER*CONST1(IJ)*TEMP1(IJ)
        ENDDO
      ENDIF

      CALL TAU_PHI_HF(IJS, IJL, MIJ, UST, Z0NEW, XLEVTAIL, TAU1, PHI1)

      DO IJ=IJS,IJL
        TAUHF(IJ) = CONST1(IJ)*TEMP1(IJ)*TAU1(IJ)
        PHIHF(IJ) = CONST2(IJ)*TEMP2(IJ)*PHI1(IJ)
      ENDDO

      DO IJ=IJS,IJL
        PHIWA(IJ)   = PHIWA(IJ)   + PHIHF_REDUC*PHIHF(IJ)
        XSTRESS(IJ) = XSTRESS(IJ) + TAUHF(IJ)*SIN(USDIRP(IJ))
        YSTRESS(IJ) = YSTRESS(IJ) + TAUHF(IJ)*COS(USDIRP(IJ))
        TAUW(IJ) = SQRT(XSTRESS(IJ)**2+YSTRESS(IJ)**2)

        TAUW(IJ) = MIN(TAUW(IJ),US2(IJ)-EPS1)
        TAUW(IJ) = MAX(TAUW(IJ),0.0_JWRB)
      ENDDO

      IF (LHOOK) CALL DR_HOOK('STRESSO',1,ZHOOK_HANDLE)

      END SUBROUTINE STRESSO
