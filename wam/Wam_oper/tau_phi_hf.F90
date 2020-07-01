      SUBROUTINE TAU_PHI_HF(IJS, IJL, LTAUWSHELTER, MIJ, USTAR, Z0,    &
     &                      XLEVTAIL, UST, TAUHF, PHIHF)

! ----------------------------------------------------------------------

!**** *TAU_PHI_HF* - COMPUTATION OF HIGH-FREQUENCY STRESS.
!                                   HIGH-FREQUENCY ENERGY FLUX.

!     PETER A.E.M. JANSSEN    KNMI      OCTOBER 90
!     JEAN BIDLOT  ECMWF  JANUARY 2017

!*    PURPOSE.
!     ---------

!       COMPUTE HIGH-FREQUENCY WAVE STRESS AND ENERGY FLUX

!**   INTERFACE.
!     ----------

!       *CALL* *TAU_PHI_HF(IJS, IJL, LTAUWSHELTER, MIJ, USTAR, UST, Z0, XLEVTAIL,
!                          TAUHF, PHIHF)
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *LTAUWSHELTER* - if true then XLEVTAIL are non zeros.
!          *MIJ* - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!          *USTAR* FRICTION VELOCITY
!          *UST*   REDUCED FRICTION VELOCITY DUE TO SHELTERING
!          *Z0*    ROUGHNESS LENGTH 
!          *XLEVTAIL* TAIL LEVEL
!          *TAUHF* HIGH-FREQUENCY STRESS
!          *PHIHF* HIGH-FREQUENCY ENERGY FLUX INTO OCEAN



!     METHOD.
!     -------

!       IT NEEDS A CALL TO INIT_X0TAUHF TO INITIALISE
!       SEE REFERENCE FOR WAVE STRESS CALCULATION.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : X0TAUHF, JTOT_TAUHF, WTAUHF, LLGCBZ0
      USE YOWFRED  , ONLY : ZPIFR
      USE YOWPCONS , ONLY : G
      USE YOWPHYS  , ONLY : ZALP   , XKAPPA
      USE YOMHOOK  , ONLY : LHOOK  , DR_HOOK

      USE YOWTEST  , ONLY : IU06
! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "omegagc.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      LOGICAL, INTENT(IN) :: LTAUWSHELTER
      INTEGER(KIND=JWIM), INTENT(IN) :: MIJ(IJS:IJL)
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: USTAR, Z0, XLEVTAIL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(INOUT) :: UST
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: TAUHF, PHIHF

      INTEGER(KIND=JWIM) :: NS
      INTEGER(KIND=JWIM) :: J, IJ

      REAL(KIND=JWRB), PARAMETER :: ZSUPMAX = 0.0_JWRB  !  LOG(1.)
      REAL(KIND=JWRB) :: XKS, OMS
      REAL(KIND=JWRB) :: OMEGA, OMEGACC
      REAL(KIND=JWRB) :: X0G
      REAL(KIND=JWRB) :: YC, Y, CM1, ZX, ZARG, ZLOG, ZBETA
      REAL(KIND=JWRB) :: FNC, FNC2, GM1
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: SQRTZ0OG, ZSUP, ZINF, DELZ
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: TAUW, XLOGGZ0, SQRTGZ0

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('TAU_PHI_HF',0,ZHOOK_HANDLE)

      GM1= 1.0_JWRB/G

!     See INIT_X0TAUHF
      X0G=X0TAUHF*G

!*    COMPUTE THE INTEGRALS 
!     ---------------------

      SQRTZ0OG(:)  = SQRT(Z0(:)*GM1)

      IF(LLGCBZ0) THEN
        DO IJ=IJS,IJL
          CALL OMEGAGC(USTAR(IJ), NS, XKS, OMS)
          ZSUP(IJ) = MIN(LOG(OMS*SQRTZ0OG(IJ)),ZSUPMAX)
        ENDDO
      ELSE
        ZSUP(:) = ZSUPMAX
      ENDIF

      DO IJ=IJS,IJL
        TAUW(IJ) = UST(IJ)**2
        XLOGGZ0(IJ) = LOG(G*Z0(IJ))
        OMEGACC = MAX(ZPIFR(MIJ(IJ)), X0G/UST(IJ))

        SQRTGZ0(IJ) = 1.0_JWRB/SQRTZ0OG(IJ)
        YC = OMEGACC*SQRTZ0OG(IJ)
        ZINF(IJ) = LOG(YC)
        DELZ(IJ) = MAX((ZSUP(IJ)-ZINF(IJ))/REAL(JTOT_TAUHF-1,JWRB),0.0_JWRB)

        TAUHF(IJ) = 0.0_JWRB
        PHIHF(IJ) = 0.0_JWRB
      ENDDO

     ! Intergrals are integrated following a change of variable : Z=LOG(Y)
      IF( LTAUWSHELTER ) THEN
        DO IJ=IJS,IJL
          DO J=1,JTOT_TAUHF
            Y         = EXP(ZINF(IJ)+REAL(J-1,JWRB)*DELZ(IJ))
            OMEGA     = Y*SQRTGZ0(IJ)
            CM1       = OMEGA*GM1
            ZX        = UST(IJ)*CM1 +ZALP
            ZARG      = XKAPPA/ZX
            ZLOG      = XLOGGZ0(IJ)+2.0_JWRB*LOG(CM1)+ZARG 
            ZLOG      = MIN(ZLOG,0.0_JWRB)
            ZBETA     = EXP(ZLOG)*ZLOG**4
            FNC2      = ZBETA*TAUW(IJ)*WTAUHF(J)*DELZ(IJ)
            TAUW(IJ)  = MAX(TAUW(IJ)-XLEVTAIL(IJ)*FNC2,0.0_JWRB)
            UST(IJ)   = SQRT(TAUW(IJ))
            TAUHF(IJ) = TAUHF(IJ) + FNC2 
            PHIHF(IJ) = PHIHF(IJ) + FNC2/Y
          ENDDO
          PHIHF(IJ) = SQRTZ0OG(IJ)*PHIHF(IJ)
        ENDDO
      ELSE
        DO IJ=IJS,IJL
          DO J=1,JTOT_TAUHF
            Y         = EXP(ZINF(IJ)+REAL(J-1,JWRB)*DELZ(IJ))
            OMEGA     = Y*SQRTGZ0(IJ)
            CM1       = OMEGA*GM1
            ZX        = UST(IJ)*CM1 +ZALP
            ZARG      = XKAPPA/ZX
            ZLOG      = XLOGGZ0(IJ)+2.0_JWRB*LOG(CM1)+ZARG 
            ZLOG      = MIN(ZLOG,0.0_JWRB)
            ZBETA     = EXP(ZLOG)*ZLOG**4
            FNC2      = ZBETA*WTAUHF(J)
            TAUHF(IJ) = TAUHF(IJ) + FNC2 
            PHIHF(IJ) = PHIHF(IJ) + FNC2/Y
          ENDDO
          TAUHF(IJ) = TAUW(IJ)*TAUHF(IJ)*DELZ(IJ)
          PHIHF(IJ) = SQRTZ0OG(IJ)*TAUW(IJ)*PHIHF(IJ)*DELZ(IJ)
        ENDDO
      ENDIF

      IF (LHOOK) CALL DR_HOOK('TAU_PHI_HF',1,ZHOOK_HANDLE)

      END SUBROUTINE TAU_PHI_HF
