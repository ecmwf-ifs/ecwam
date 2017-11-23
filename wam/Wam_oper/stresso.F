      SUBROUTINE STRESSO (F, SL, IJS, IJL,
     &                    MIJ, RHOWGDFTH,
     &                    THWNEW, USNEW, Z0NEW, ROAIRN,
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

!        *CALL* *STRESSO (F, SL, IJS, IJL,
!    &                    MIJ, RHOWGDFTH,
!    &                    THWNEW, USNEW, Z0NEW, ROAIRN,
!    &                    TAUW, PHIWA)*
!         *F*           - WAVE SPECTRUM.
!         *SL*          - TOTAL SOURCE FUNCTION ARRAY.
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
!       BECAUSE ARRAY *SL* IS USED, ONLY THE INPUT SOURCE
!       HAS TO BE STORED IN *SL* (CALL FIRST SINPUT, THEN
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

      USE YOWCOUP  , ONLY : ALPHA    ,TAUWSHELTER, ITSHELT
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,TH       ,
     &            COSTH    ,SINTH    ,FR5      ,FRATIO
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER  ,YINVEPS
      USE YOWSHAL  , ONLY : CINV     ,INDEP
      USE YOWTABL  , ONLY : EPS1
      USE YOWSTAT  , ONLY : IPHYS

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "tau_phi_hf.intfb.h"

      INTEGER, INTENT(IN) :: IJS, IJL
      INTEGER, DIMENSION(IJS:IJL), INTENT(IN) :: MIJ

      REAL,DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F, SL
      REAL,DIMENSION(IJS:IJL,NFRE), INTENT(IN) :: RHOWGDFTH
      REAL,DIMENSION(IJS:IJL), INTENT(IN) :: THWNEW, USNEW, Z0NEW,ROAIRN
      REAL,DIMENSION(IJS:IJL), INTENT(OUT) :: TAUW, PHIWA


      INTEGER :: IJ, M, K, I, J, II

      REAL :: CONST, CNST
      REAL :: XI, XJ, DELI1, DELI2, DELJ1, DELJ2, XK, DELK1, DELK2
      REAL :: PHI2
      REAL :: ABS_TAUWSHELTER, GM1
      REAL :: COSW, FCOSW2
      REAL :: ZHOOK_HANDLE
      REAL, DIMENSION(IJS:IJL) :: TAUHF, TEMP1, CONST1, XSTRESS, YSTRESS
      REAL, DIMENSION(IJS:IJL) :: PHIHF, TEMP2, CONST2
      REAL, DIMENSION(IJS:IJL) :: TAU1, PHI1, XLEVTAIL 
      REAL, DIMENSION(IJS:IJL) :: CMRHOWGDFTH
      REAL, DIMENSION(IJS:IJL) :: US2, TAUX, TAUY, TAUPX, TAUPY
      REAL, DIMENSION(IJS:IJL) :: USDIRP, UST, UST2 
      REAL, DIMENSION(IJS:IJL) :: SUMT, SUMX, SUMY

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('STRESSO',0,ZHOOK_HANDLE)
#endif

!*    1. PRECOMPUTE FREQUENCY SCALING.
!        -----------------------------

      CONST = DELTH*(ZPI)**4/G

      ABS_TAUWSHELTER=ABS(TAUWSHELTER)

      GM1 = 1./G

      DO IJ=IJS,IJL
        CONST1(IJ)  = CONST*FR5(MIJ(IJ))*GM1
        CONST2(IJ)  = ROAIRN(IJ)*CONST*FR5(MIJ(IJ))
      ENDDO

!*    2. COMPUTE WAVE STRESS OF ACTUAL BLOCK.
!        ------------------------------------

!*    2.2 INTEGRATE INPUT SOURCE FUNCTION OVER FREQUENCY AND DIRECTIONS.
!         --------------------------------------------------------------

      DO IJ=IJS,IJL
        PHIWA(IJ)   = 0.
        XSTRESS(IJ) = 0.
        YSTRESS(IJ) = 0.
      ENDDO
      DO M=1,MAXVAL(MIJ(:))
!     THE INTEGRATION ONLY UP TO FR=MIJ SINCE RHOWGDFTH=0 FOR FR>MIJ
        K=1
        DO IJ=IJS,IJL
          SUMT(IJ)   = SL(IJ,K,M)
!         until we have figure out how to deal with negative input
          CNST = MAX(SL(IJ,K,M),0.)
          SUMX(IJ) = CNST*SINTH(K)
          SUMY(IJ) = CNST*COSTH(K)
        ENDDO
        DO K=2,NANG
          DO IJ=IJS,IJL
            SUMT(IJ) = SUMT(IJ) + SL(IJ,K,M)
!           until we have figure out how to deal with negative input
            CNST = MAX(SL(IJ,K,M),0.)
            SUMX(IJ) = SUMX(IJ) + CNST*SINTH(K)
            SUMY(IJ) = SUMY(IJ) + CNST*COSTH(K)
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
        XSTRESS(IJ) = XSTRESS(IJ)/MAX(ROAIRN(IJ),1.)
        YSTRESS(IJ) = YSTRESS(IJ)/MAX(ROAIRN(IJ),1.)
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
        COSW     = MAX(COS(TH(K)-THWNEW(IJ)),0.)
        FCOSW2   = F(IJ,K,MIJ(IJ))*COSW**2
        TEMP1(IJ) = FCOSW2*COSW
        TEMP2(IJ) = FCOSW2 
      ENDDO

      DO K=2,NANG
        DO IJ=IJS,IJL
          COSW     = MAX(COS(TH(K)-THWNEW(IJ)),0.)
          FCOSW2   = F(IJ,K,MIJ(IJ))*COSW**2
          TEMP1(IJ) = TEMP1(IJ) + FCOSW2*COSW
          TEMP2(IJ) = TEMP2(IJ) + FCOSW2 
        ENDDO
      ENDDO

      IF (ITSHELT.EQ.0) THEN
        DO IJ=IJS,IJL
          XLEVTAIL(IJ)=0.
        ENDDO
      ELSE
        DO IJ=IJS,IJL
          XLEVTAIL(IJ)=TAUWSHELTER*CONST1(IJ)*TEMP1(IJ)
        ENDDO
      ENDIF

      CALL TAU_PHI_HF(IJS, IJL, MIJ, UST, Z0NEW, XLEVTAIL,
     &                TAU1, PHI1)
      DO IJ=IJS,IJL
        TAUHF(IJ) = CONST1(IJ)*TEMP1(IJ)*UST2(IJ)*TAU1(IJ)
        PHIHF(IJ) = CONST2(IJ)*TEMP2(IJ)*UST2(IJ)*PHI1(IJ)
      ENDDO

      DO IJ=IJS,IJL
        PHIWA(IJ)   = PHIWA(IJ)   + PHIHF(IJ)
        XSTRESS(IJ) = XSTRESS(IJ) + TAUHF(IJ)*SIN(USDIRP(IJ))
        YSTRESS(IJ) = YSTRESS(IJ) + TAUHF(IJ)*COS(USDIRP(IJ))
        TAUW(IJ) = SQRT(XSTRESS(IJ)**2+YSTRESS(IJ)**2)

        TAUW(IJ) = MIN(TAUW(IJ),US2(IJ)-EPS1)
        TAUW(IJ) = MAX(TAUW(IJ),0.)
      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('STRESSO',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE STRESSO
