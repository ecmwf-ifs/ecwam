      SUBROUTINE STRESSO (F, SL, SPOS, IJS, IJL,                        &
     &                    MIJ, RHOWGDFTH,                               &
     &                    THWNEW, USNEW, Z0NEW, ROAIRN,                 &
     &                    TAUW, TAUWDIR, PHIWA, LLPHIWA)

! ----------------------------------------------------------------------

!**** *STRESSO* - COMPUTATION OF WAVE STRESS.

!     H. GUNTHER      GKSS/ECMWF NOVEMBER   1989 CODE MOVED FROM SINPUT.
!     P.A.E.M. JANSSEN     KNMI  AUGUST     1990
!     J. BIDLOT            ECMWF FEBRUARY   1996-97
!     S. ABDALLA           ECMWF OCTOBER    2001 INTRODUCTION OF VARIABLE
!                                                AIR DENSITY
!     J. BIDLOT            ECMWF            2007  ADD MIJ
!     P.A.E.M. JANSSEN     ECMWF            2011  ADD FLUX CALULATIONS

!*    PURPOSE.jjjj
!     --------

!       COMPUTE NORMALIZED WAVE STRESS FROM INPUT SOURCE FUNCTION
 
!**   INTERFACE.
!     ----------

!        *CALL* *STRESSO (F, SPOS, IJS, IJL,
!    &                    MIJ, RHOWGDFTH,
!    &                    THWNEW, USNEW, Z0NEW, ROAIRN,
!    &                    TAUW, TAUWDIR, PHIWA)*
!         *F*           - WAVE SPECTRUM.
!         *SL*          - WIND INPUT SOURCE FUNCTION ARRAY (positive and negative contributions).
!         *SPOS*        - POSITIVE WIND INPUT SOURCE FUNCTION ARRAY.
!         *IJS*         - INDEX OF FIRST GRIDPOINT.
!         *IJL*         - INDEX OF LAST GRIDPOINT.
!         *MIJ*         - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE for fluxes calculation.
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
!         *TAUWDIR*     - KINEMATIC WAVE STRESS DIRECTION
!         *PHIWA*       - ENERGY FLUX FROM WIND INTO WAVES INTEGRATED
!                         OVER THE FULL FREQUENCY RANGE.
!         *LLPHIWA*     - TRUE IF PHIWA NEEDS TO BE COMPUTED

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

      USE YOWCOUP  , ONLY : LLGCBZ0
      USE YOWFRED  , ONLY : FR       ,RHOWG_DFIM ,DELTH    ,TH       ,    &
     &            COSTH    ,SINTH    ,FR5
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,GM1      ,ZPI        ,ZPI4GM2
      USE YOWPHYS  , ONLY : TAUWSHELTER
      USE YOWSHAL  , ONLY : CINV     ,INDEP
      USE YOWTABL  , ONLY : EPS1
      USE YOWSTAT  , ONLY : IPHYS
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "tau_phi_hf.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL), INTENT(IN) :: MIJ

      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F, SL, SPOS
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE), INTENT(IN) :: RHOWGDFTH
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: THWNEW, USNEW, Z0NEW, ROAIRN
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: TAUW, TAUWDIR, PHIWA
      LOGICAL, INTENT(IN) :: LLPHIWA


      INTEGER(KIND=JWIM) :: IJ, M, K, I, J, II

      REAL(KIND=JWRB) :: TAUTOUS2
      REAL(KIND=JWRB) :: C2 
      REAL(KIND=JWRB) :: XI, XJ, DELI1, DELI2, DELJ1, DELJ2, XK, DELK1, DELK2
      REAL(KIND=JWRB) :: PHI2
      REAL(KIND=JWRB) :: COSW, FCOSW2
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: TAUHF, F1DCOS3, CONST1, XSTRESS, YSTRESS
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: PHIHF, F1DCOS2, CONST2
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: TAU1, PHI1, XLEVTAIL 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: CMRHOWGDFTH
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: US2, TAUX, TAUY, TAUPX, TAUPY
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: USDIRP, UST
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: SUMT, SUMX, SUMY

      LOGICAL :: LTAUWSHELTER

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('STRESSO',0,ZHOOK_HANDLE)

      DO IJ=IJS,IJL
        CONST1(IJ) = ZPI4GM2*FR5(MIJ(IJ))
      ENDDO

      DO IJ=IJS,IJL
        PHIWA(IJ)   = 0.0_JWRB
        XSTRESS(IJ) = 0.0_JWRB
        YSTRESS(IJ) = 0.0_JWRB
      ENDDO

!*    CONTRIBUTION TO THE WAVE STRESS FROM THE NEGATIVE PART OF THE WIND INPUT
!     ------------------------------------------------------------------------

      IF ( LLPHIWA ) THEN
!     full energy flux due to negative Sinput (SL-SPOS)
!     we assume that above NFRE, the contibutions can be neglected
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=IJS,IJL
              PHIWA(IJ) = PHIWA(IJ) + (SL(IJ,K,M)-SPOS(IJ,K,M))*RHOWG_DFIM(M)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!*    CALCULATE LOW-FREQUENCY CONTRIBUTION TO STRESS AND ENERGY FLUX (positive sinput).
!     ---------------------------------------------------------------------------------
      DO M=1,MAXVAL(MIJ(:))
!     The integration only up to FR=MIJ since RHOWGDFTH=0 for FR>MIJ
        K=1
        DO IJ=IJS,IJL
          SUMX(IJ) = SPOS(IJ,K,M)*SINTH(K)
          SUMY(IJ) = SPOS(IJ,K,M)*COSTH(K)
        ENDDO
        DO K=2,NANG
          DO IJ=IJS,IJL
            SUMX(IJ) = SUMX(IJ) + SPOS(IJ,K,M)*SINTH(K)
            SUMY(IJ) = SUMY(IJ) + SPOS(IJ,K,M)*COSTH(K)
          ENDDO
        ENDDO
        DO IJ=IJS,IJL
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

      IF ( LLPHIWA ) THEN
        DO M=1,MAXVAL(MIJ(:))
!       THE INTEGRATION ONLY UP TO FR=MIJ SINCE RHOWGDFTH=0 FOR FR>MIJ
          K=1
          DO IJ=IJS,IJL
            SUMT(IJ) = SPOS(IJ,K,M)
          ENDDO
          DO K=2,NANG
            DO IJ=IJS,IJL
              SUMT(IJ) = SUMT(IJ) + SPOS(IJ,K,M)
            ENDDO
          ENDDO
          DO IJ=IJS,IJL
            PHIWA(IJ) = PHIWA(IJ) + SUMT(IJ)*RHOWGDFTH(IJ,M)
          ENDDO
        ENDDO
      ENDIF

!*    CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS and energy flux (positive sinput).
!     ----------------------------------------------------------------------------------

      DO IJ=IJS,IJL
        US2(IJ)=USNEW(IJ)**2
      ENDDO

      IF ( IPHYS.EQ.0 .OR. TAUWSHELTER == 0.0_JWRB) THEN
        LTAUWSHELTER = .FALSE.
        DO IJ=IJS,IJL
          USDIRP(IJ)=THWNEW(IJ)
          UST(IJ)=USNEW(IJ)
        ENDDO
      ELSE
        LTAUWSHELTER = .TRUE.
        DO IJ=IJS,IJL
          TAUX(IJ)=US2(IJ)*SIN(THWNEW(IJ))
          TAUY(IJ)=US2(IJ)*COS(THWNEW(IJ))
          TAUPX(IJ)=TAUX(IJ)-TAUWSHELTER*XSTRESS(IJ)
          TAUPY(IJ)=TAUY(IJ)-TAUWSHELTER*YSTRESS(IJ)
          USDIRP(IJ)=ATAN2(TAUPX(IJ),TAUPY(IJ))
          UST(IJ)=(TAUPX(IJ)**2+TAUPY(IJ)**2)**0.25_JWRB
        ENDDO
      ENDIF

      K=1
      DO IJ=IJS,IJL
        COSW     = MAX(COS(TH(K)-THWNEW(IJ)),0.0_JWRB)
        FCOSW2   = DELTH*F(IJ,K,MIJ(IJ))*COSW**2
        F1DCOS3(IJ) = FCOSW2*COSW
        F1DCOS2(IJ) = FCOSW2
      ENDDO

      DO K=2,NANG
        DO IJ=IJS,IJL
          COSW     = MAX(COS(TH(K)-THWNEW(IJ)),0.0_JWRB)
          FCOSW2   = DELTH*F(IJ,K,MIJ(IJ))*COSW**2
          F1DCOS3(IJ) = F1DCOS3(IJ) + FCOSW2*COSW
          F1DCOS2(IJ) = F1DCOS2(IJ) + FCOSW2 
        ENDDO
      ENDDO

      IF( LTAUWSHELTER ) THEN
        DO IJ=IJS,IJL
          XLEVTAIL(IJ)=TAUWSHELTER*CONST1(IJ)*F1DCOS3(IJ)
        ENDDO
      ELSE
        XLEVTAIL(:) = 0.0_JWRB
      ENDIF
  
      CALL TAU_PHI_HF(IJS, IJL, LTAUWSHELTER, MIJ, F1DCOS2, USNEW, Z0NEW, &
     &                XLEVTAIL, UST, TAU1, PHI1, LLPHIWA)

      DO IJ=IJS,IJL
        TAUHF(IJ) = CONST1(IJ)*F1DCOS3(IJ)*TAU1(IJ)
        XSTRESS(IJ) = XSTRESS(IJ) + TAUHF(IJ)*SIN(USDIRP(IJ))
        YSTRESS(IJ) = YSTRESS(IJ) + TAUHF(IJ)*COS(USDIRP(IJ))
        TAUW(IJ) = SQRT(XSTRESS(IJ)**2+YSTRESS(IJ)**2)
        TAUW(IJ) = MAX(TAUW(IJ),0.0_JWRB)
        TAUWDIR(IJ) = ATAN2(XSTRESS(IJ),YSTRESS(IJ))
      ENDDO

!      IF (.NOT. LLGCBZ0) THEN
      TAUTOUS2 = 1.0_JWRB/(1.0_JWRB+EPS1)
      DO IJ=IJS,IJL
        TAUW(IJ) = MIN(TAUW(IJ),US2(IJ)*TAUTOUS2)
      ENDDO
!      ENDIF

      IF ( LLPHIWA ) THEN
        C2 = (ZPI)**4*GM1
        DO IJ=IJS,IJL
          CONST2(IJ) = ROAIRN(IJ)*C2*FR5(MIJ(IJ))
          PHIWA(IJ) = PHIWA(IJ) + CONST2(IJ)*F1DCOS2(IJ)*PHI1(IJ)
        ENDDO
      ENDIF

      IF (LHOOK) CALL DR_HOOK('STRESSO',1,ZHOOK_HANDLE)

      END SUBROUTINE STRESSO
