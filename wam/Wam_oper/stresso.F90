      SUBROUTINE STRESSO (IJS, IJL, KIJS, KIJL, GFL, SL, SPOS,          &
     &                    CINV,                                         &
     &                    THWNEW, USNEW, Z0NEW, ROAIRN, RNFAC,          &
     &                    TAUW, TAUWDIR, PHIWA, LLPHIWA)

! ----------------------------------------------------------------------

!**** *STRESSO* - COMPUTATION OF WAVE STRESS.

!     H. GUNTHER      GKSS/ECMWF NOVEMBER   1989 CODE MOVED FROM SINPUT.
!     P.A.E.M. JANSSEN     KNMI  AUGUST     1990
!     J. BIDLOT            ECMWF FEBRUARY   1996-97
!     S. ABDALLA           ECMWF OCTOBER    2001 INTRODUCTION OF VARIABLE
!                                                AIR DENSITY
!     P.A.E.M. JANSSEN     ECMWF            2011  ADD FLUX CALULATIONS

!*    PURPOSE.
!     --------

!       COMPUTE NORMALIZED WAVE STRESS FROM INPUT SOURCE FUNCTION

!**   INTERFACE.
!     ----------

!        *CALL* *STRESSO (IJS, IJL, KIJS, KIJL, GFL, SL, SPOS,
!    &                    CINV,
!    &                    THWNEW, USNEW, Z0NEW, ROAIRN, RNFAC,
!    &                    TAUW, TAUWDIR, PHIWA)*
!         *IJS:IJL*     - 1st DIMENSION of GFL
!         *KIJS*        - INDEX OF FIRST GRIDPOINT.
!         *KIJL*        - INDEX OF LAST GRIDPOINT.
!         *GFL*         - WAVE SPECTRUM.
!         *SL*          - WIND INPUT SOURCE FUNCTION ARRAY (positive and negative contributions).
!         *SPOS*        - POSITIVE WIND INPUT SOURCE FUNCTION ARRAY.
!         *CINV*        - INVERSE PHASE VELOCITY.
!         *THWNEW*      - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                         NOTATION (POINTING ANGLE OF WIND VECTOR,
!                         CLOCKWISE FROM NORTH).
!         *USNEW*       - NEW FRICTION VELOCITY IN M/S.
!         *Z0NEW*       - ROUGHNESS LENGTH IN M.
!         *ROAIRN*      - AIR DENSITY IN KG/M**3.
!         *RNFAC*       - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
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
      USE YOWPHYS  , ONLY : TAUWSHELTER
      USE YOWTABL  , ONLY : EPS1
      USE YOWSTAT  , ONLY : IPHYS

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "tau_phi_hf.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: GFL

      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE), INTENT(IN) :: CINV

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: SL, SPOS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: THWNEW, USNEW, Z0NEW, ROAIRN, RNFAC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: TAUW, TAUWDIR, PHIWA
      LOGICAL, INTENT(IN) :: LLPHIWA


      INTEGER(KIND=JWIM) :: IJ, M, K, I, J, II

      REAL(KIND=JWRB) :: TAUTOUS2
      REAL(KIND=JWRB) :: COSW, FCOSW2
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: XSTRESS, YSTRESS
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TAUHF, PHIHF
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: CMRHOWGDFTH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: US2, TAUX, TAUY, TAUPX, TAUPY
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: USDIRP, UST
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: SUMT, SUMX, SUMY

      LOGICAL :: LTAUWSHELTER

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('STRESSO',0,ZHOOK_HANDLE)

      DO IJ=KIJS,KIJL
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
            DO IJ=KIJS,KIJL
              PHIWA(IJ) = PHIWA(IJ) + (SL(IJ,K,M)-SPOS(IJ,K,M))*RHOWG_DFIM(M)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!*    CALCULATE LOW-FREQUENCY CONTRIBUTION TO STRESS AND ENERGY FLUX (positive sinput).
!     ---------------------------------------------------------------------------------
      DO M=1,NFRE
        K=1
        DO IJ=KIJS,KIJL
          SUMX(IJ) = SPOS(IJ,K,M)*SINTH(K)
          SUMY(IJ) = SPOS(IJ,K,M)*COSTH(K)
        ENDDO
        DO K=2,NANG
          DO IJ=KIJS,KIJL
            SUMX(IJ) = SUMX(IJ) + SPOS(IJ,K,M)*SINTH(K)
            SUMY(IJ) = SUMY(IJ) + SPOS(IJ,K,M)*COSTH(K)
          ENDDO
        ENDDO
        DO IJ=KIJS,KIJL
          CMRHOWGDFTH(IJ) = CINV(IJ,M)*RHOWG_DFIM(M)
          XSTRESS(IJ) = XSTRESS(IJ) + SUMX(IJ)*CMRHOWGDFTH(IJ)
          YSTRESS(IJ) = YSTRESS(IJ) + SUMY(IJ)*CMRHOWGDFTH(IJ)
        ENDDO
      ENDDO

!     TAUW is the kinematic wave stress !
      DO IJ=KIJS,KIJL
        XSTRESS(IJ) = XSTRESS(IJ)/MAX(ROAIRN(IJ),1.0_JWRB)
        YSTRESS(IJ) = YSTRESS(IJ)/MAX(ROAIRN(IJ),1.0_JWRB)
      ENDDO

      IF ( LLPHIWA ) THEN
        DO M=1,NFRE
          K=1
          DO IJ=KIJS,KIJL
            SUMT(IJ) = SPOS(IJ,K,M)
          ENDDO
          DO K=2,NANG
            DO IJ=KIJS,KIJL
              SUMT(IJ) = SUMT(IJ) + SPOS(IJ,K,M)
            ENDDO
          ENDDO
          DO IJ=KIJS,KIJL
            PHIWA(IJ) = PHIWA(IJ) + SUMT(IJ)*RHOWG_DFIM(M)
          ENDDO
        ENDDO
      ENDIF

!*    CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS and energy flux (positive sinput).
!     ----------------------------------------------------------------------------------

      DO IJ=KIJS,KIJL
        US2(IJ)=USNEW(IJ)**2
      ENDDO

      IF ( IPHYS.EQ.0 .OR. TAUWSHELTER == 0.0_JWRB) THEN
        LTAUWSHELTER = .FALSE.
        DO IJ=KIJS,KIJL
          USDIRP(IJ)=THWNEW(IJ)
          UST(IJ)=USNEW(IJ)
        ENDDO
      ELSE
        LTAUWSHELTER = .TRUE.
        DO IJ=KIJS,KIJL
          TAUX(IJ)=US2(IJ)*SIN(THWNEW(IJ))
          TAUY(IJ)=US2(IJ)*COS(THWNEW(IJ))
          TAUPX(IJ)=TAUX(IJ)-TAUWSHELTER*XSTRESS(IJ)
          TAUPY(IJ)=TAUY(IJ)-TAUWSHELTER*YSTRESS(IJ)
          USDIRP(IJ)=ATAN2(TAUPX(IJ),TAUPY(IJ))
          UST(IJ)=(TAUPX(IJ)**2+TAUPY(IJ)**2)**0.25_JWRB
        ENDDO
      ENDIF

      CALL TAU_PHI_HF(IJS, IJL, KIJS, KIJL, LTAUWSHELTER, USNEW, Z0NEW, &
     &                GFL, THWNEW, ROAIRN, RNFAC,                       &
     &                UST, TAUHF, PHIHF, LLPHIWA)

      DO IJ=KIJS,KIJL
        XSTRESS(IJ) = XSTRESS(IJ) + TAUHF(IJ)*SIN(USDIRP(IJ))
        YSTRESS(IJ) = YSTRESS(IJ) + TAUHF(IJ)*COS(USDIRP(IJ))
        TAUW(IJ) = SQRT(XSTRESS(IJ)**2+YSTRESS(IJ)**2)
        TAUW(IJ) = MAX(TAUW(IJ),0.0_JWRB)
        TAUWDIR(IJ) = ATAN2(XSTRESS(IJ),YSTRESS(IJ))
      ENDDO
      IF ( .NOT. LLGCBZ0) THEN
        TAUTOUS2 = 1.0_JWRB/(1.0_JWRB+EPS1)
        DO IJ=KIJS,KIJL
          TAUW(IJ) = MIN(TAUW(IJ),US2(IJ)*TAUTOUS2)
        ENDDO
     ENDIF

      IF ( LLPHIWA ) THEN
        DO IJ=KIJS,KIJL
          PHIWA(IJ) = PHIWA(IJ) + PHIHF(IJ)
        ENDDO
      ENDIF

      IF (LHOOK) CALL DR_HOOK('STRESSO',1,ZHOOK_HANDLE)

      END SUBROUTINE STRESSO
