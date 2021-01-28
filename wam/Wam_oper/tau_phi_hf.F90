       SUBROUTINE TAU_PHI_HF(IJS, IJL, LTAUWSHELTER, USTAR, Z0,         &
     &                      F, THWNEW, ROAIRN, RNFAC,                   &
     &                      UST, TAUHF, PHIHF, LLPHIHF)

! ----------------------------------------------------------------------

!**** *TAU_PHI_HF* - COMPUTATION OF HIGH-FREQUENCY STRESS.
!                                   HIGH-FREQUENCY ENERGY FLUX.

!     PETER A.E.M. JANSSEN    KNMI      OCTOBER 90
!     JEAN BIDLOT  ECMWF  JANUARY 2017

!*    PURPOSE.
!     ---------

!       COMPUTE HIGH-FREQUENCY WAVE STRESS AND ENERGY FLUX

!**   INTERFACE.
!     ---------

!       *CALL* *TAU_PHI_HF(IJS, IJL, LTAUWSHELTER, USTAR, UST, Z0,
!                          F, THWNEW, ROAIRN, RNFAC, &
!                          UST, TAUHF, PHIHF, LLPHIHF)
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *LTAUWSHELTER* - if true then TAUWSHELTER 
!          *F*           - WAVE SPECTRUM.
!          *THWNEW*      - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!          *ROAIRN*      - AIR DENSITY IN KG/M**3.
!          *RNFAC* - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
!          *USTAR* FRICTION VELOCITY
!          *UST*   REDUCED FRICTION VELOCITY DUE TO SHELTERING
!          *Z0*    ROUGHNESS LENGTH 
!          *TAUHF* HIGH-FREQUENCY STRESS
!          *PHIHF* HIGH-FREQUENCY ENERGY FLUX INTO OCEAN
!          *LLPHIHF* TRUE IF PHIHF NEEDS TO COMPUTED



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

      USE YOWCOUP  , ONLY : X0TAUHF, JTOT_TAUHF, WTAUHF, LLGCBZ0, LLNORMAGAM 
      USE YOWFRED  , ONLY : ZPIFR  , FR5,   TH    ,DELTH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G      , GM1       ,ZPI    , ZPI4GM1,  ZPI4GM2
      USE YOWPHYS  , ONLY : ZALP   , XKAPPA    ,TAUWSHELTER, GAMNCONST, RN1_RN
      USE YOWPHYS  , ONLY : DELTA_THETA_RN
      USE YOMHOOK  , ONLY : LHOOK  , DR_HOOK
!debile debug
      USE YOWSTAT  , ONLY : cdtpro 

      USE YOWTEST  , ONLY : IU06
! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "omegagc.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      LOGICAL, INTENT(IN) :: LTAUWSHELTER
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: USTAR, Z0
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: THWNEW, ROAIRN, RNFAC
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F

      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(INOUT) :: UST
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: TAUHF, PHIHF
      LOGICAL, INTENT(IN) :: LLPHIHF

      INTEGER(KIND=JWIM) :: NS
      INTEGER(KIND=JWIM) :: J, IJ, K

      REAL(KIND=JWRB), PARAMETER :: ZSUPMAX = 0.0_JWRB  !  LOG(1.)
      REAL(KIND=JWRB) :: XKS, OMS
      REAL(KIND=JWRB) :: OMEGA, OMEGACC
      REAL(KIND=JWRB) :: X0G
      REAL(KIND=JWRB) :: YC, Y, CM1, ZX, ZARG, ZLOG, ZBETA
      REAL(KIND=JWRB) :: FNC, FNC2
      REAL(KIND=JWRB) :: GAMNORMA ! RENORMALISATION FACTOR OF THE GROWTH RATE
      REAL(KIND=JWRB) :: ZN, GAMCFR5
      REAL(KIND=JWRB) :: COSW, FCOSW2 
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
!debile debug
      REAL(KIND=JWRB) :: gam, eps

      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: SQRTZ0OG, ZSUP, ZINF, DELZ
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: TAUL, XLOGGZ0, SQRTGZ0
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: USTPH
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: CONST, CONSTTAU, CONSTPHI
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: F1DCOS2, F1DCOS3 

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('TAU_PHI_HF',0,ZHOOK_HANDLE)

!     See INIT_X0TAUHF
      X0G=X0TAUHF*G

      IF (LLPHIHF) USTPH(:) = UST(:)

!*    COMPUTE THE INTEGRALS 
!     ---------------------

      DO IJ=IJS,IJL
        XLOGGZ0(IJ) = LOG(G*Z0(IJ))
        OMEGACC = MAX(ZPIFR(NFRE), X0G/UST(IJ))
        SQRTZ0OG(IJ)  = SQRT(Z0(IJ)*GM1)
        SQRTGZ0(IJ) = 1.0_JWRB/SQRTZ0OG(IJ)
        YC = OMEGACC*SQRTZ0OG(IJ)
        ZINF(IJ) = LOG(YC)
      ENDDO

      DO IJ=IJS,IJL
        CONSTTAU(IJ) = ZPI4GM2*FR5(NFRE)
      ENDDO

      K=1
      DO IJ=IJS,IJL
        COSW     = MAX(COS(TH(K)-THWNEW(IJ)),0.0_JWRB)
        FCOSW2   = F(IJ,K,NFRE)*COSW**2
        F1DCOS3(IJ) = FCOSW2*COSW
        F1DCOS2(IJ) = FCOSW2
      ENDDO
      DO K=2,NANG
        DO IJ=IJS,IJL
          COSW     = MAX(COS(TH(K)-THWNEW(IJ)),0.0_JWRB)
          FCOSW2   = F(IJ,K,NFRE)*COSW**2
          F1DCOS3(IJ) = F1DCOS3(IJ) + FCOSW2*COSW
          F1DCOS2(IJ) = F1DCOS2(IJ) + FCOSW2 
        ENDDO
      ENDDO
      DO IJ=IJS,IJL
        F1DCOS3(IJ) = DELTH*F1DCOS3(IJ)
        F1DCOS2(IJ) = DELTH*F1DCOS2(IJ)
      ENDDO

      IF(LLNORMAGAM) THEN
        GAMCFR5 = GAMNCONST*FR5(NFRE)
        DO IJ=IJS,IJL
          CONST(IJ) = GAMCFR5*RNFAC(IJ)*F1DCOS2(IJ)*SQRTGZ0(IJ)
        ENDDO
      ELSE
!!!!        CONST(:) = 0.0_JWRB
!!debile debug
        GAMCFR5 = GAMNCONST*FR5(NFRE)
        DO IJ=IJS,IJL
          CONST(IJ) = GAMCFR5*RNFAC(IJ)*F1DCOS2(IJ)*SQRTGZ0(IJ)
        ENDDO

      ENDIF


!     TAUHF :
      IF(LLGCBZ0) THEN
        DO IJ=IJS,IJL
          CALL OMEGAGC(USTAR(IJ), NS, XKS, OMS)
          ZSUP(IJ) = MIN(LOG(OMS*SQRTZ0OG(IJ)),ZSUPMAX)
        ENDDO
      ELSE
        ZSUP(:) = ZSUPMAX
      ENDIF

      DO IJ=IJS,IJL
        TAUL(IJ) = UST(IJ)**2
        DELZ(IJ) = MAX((ZSUP(IJ)-ZINF(IJ))/REAL(JTOT_TAUHF-1,JWRB),0.0_JWRB)
        TAUHF(IJ) = 0.0_JWRB
      ENDDO

     ! Intergrals are integrated following a change of variable : Z=LOG(Y)
      IF( LTAUWSHELTER ) THEN
        DO IJ=IJS,IJL
          DO J=1,JTOT_TAUHF
            Y         = EXP(ZINF(IJ)+REAL(J-1,JWRB)*DELZ(IJ))
            OMEGA     = Y*SQRTGZ0(IJ)
            CM1       = OMEGA*GM1
            ZX        = UST(IJ)*CM1 + ZALP
            ZARG      = XKAPPA/ZX
            ZLOG      = XLOGGZ0(IJ)+2.0_JWRB*LOG(CM1)+ZARG 
            ZLOG      = MIN(ZLOG,0.0_JWRB)
            ZBETA     = EXP(ZLOG)*ZLOG**4
            ZN        = CONST(IJ)*ZBETA*UST(IJ)*Y

!!debile debug
       eps = ROAIRN(IJ)/1025._JWRB
       gam = ZN*EPS*XKAPPA*UST(IJ)*2._JWRB*ZPI*G*OMEGA**2/(DELTA_THETA_RN*RNFAC(IJ)*F1DCOS2(IJ)*ZPIFR(NFRE)**5)
       write(iu06,'(a9,1x,a12,1x,2(f14.8,1x))') 'debile_hf',cdtpro,OMEGA**2/G, gam/OMEGA
      IF(.not. LLNORMAGAM) zn=0.0_JWRB



            GAMNORMA  = (1.0_JWRB + RN1_RN*ZN)/(1.0_JWRB + ZN)
            FNC2      = F1DCOS3(IJ)*CONSTTAU(IJ)* ZBETA*TAUL(IJ)*WTAUHF(J)*DELZ(IJ) * GAMNORMA
            TAUL(IJ)  = MAX(TAUL(IJ)-TAUWSHELTER*FNC2,0.0_JWRB)

            UST(IJ)   = SQRT(TAUL(IJ))
            TAUHF(IJ) = TAUHF(IJ) + FNC2
          ENDDO
        ENDDO
      ELSE
        DO IJ=IJS,IJL
          DO J=1,JTOT_TAUHF
            Y         = EXP(ZINF(IJ)+REAL(J-1,JWRB)*DELZ(IJ))
            OMEGA     = Y*SQRTGZ0(IJ)
            CM1       = OMEGA*GM1
            ZX        = UST(IJ)*CM1 + ZALP
            ZARG      = XKAPPA/ZX
            ZLOG      = XLOGGZ0(IJ)+2.0_JWRB*LOG(CM1)+ZARG 
            ZLOG      = MIN(ZLOG,0.0_JWRB)
            ZBETA     = EXP(ZLOG)*ZLOG**4
            FNC2      = ZBETA*WTAUHF(J)
            ZN        = CONST(IJ)*ZBETA*UST(IJ)*Y
            GAMNORMA  = (1.0_JWRB + RN1_RN*ZN)/(1.0_JWRB + ZN)
            TAUHF(IJ) = TAUHF(IJ) + FNC2 * GAMNORMA
          ENDDO
          TAUHF(IJ) = F1DCOS3(IJ)*CONSTTAU(IJ) * TAUL(IJ)*TAUHF(IJ)*DELZ(IJ)
        ENDDO
      ENDIF


      PHIHF(:) = 0.0_JWRB
      IF (LLPHIHF) THEN
!     PHIHF:
!     We are neglecting the gravity-capillary contribution 
!     Recompute DELZ over the full interval
      DO IJ=IJS,IJL
        TAUL(IJ) = USTPH(IJ)**2
        ZSUP(IJ) = ZSUPMAX
        DELZ(IJ) = MAX((ZSUP(IJ)-ZINF(IJ))/REAL(JTOT_TAUHF-1,JWRB),0.0_JWRB)
      ENDDO

      DO IJ=IJS,IJL
        CONSTPHI(IJ) = ROAIRN(IJ)*ZPI4GM1*FR5(NFRE)
      ENDDO

     ! Intergrals are integrated following a change of variable : Z=LOG(Y)
      IF( LTAUWSHELTER ) THEN
        DO IJ=IJS,IJL
          DO J=1,JTOT_TAUHF
            Y         = EXP(ZINF(IJ)+REAL(J-1,JWRB)*DELZ(IJ))
            OMEGA     = Y*SQRTGZ0(IJ)
            CM1       = OMEGA*GM1
            ZX        = USTPH(IJ)*CM1 + ZALP
            ZARG      = XKAPPA/ZX
            ZLOG      = XLOGGZ0(IJ)+2.0_JWRB*LOG(CM1)+ZARG 
            ZLOG      = MIN(ZLOG,0.0_JWRB)
            ZBETA     = EXP(ZLOG)*ZLOG**4
            ZN        = CONST(IJ)*ZBETA*USTPH(IJ)*Y
            GAMNORMA  = (1.0_JWRB + RN1_RN*ZN)/(1.0_JWRB + ZN)
            FNC2      = ZBETA*TAUL(IJ)*WTAUHF(J)*DELZ(IJ) * GAMNORMA
            TAUL(IJ)  = MAX(TAUL(IJ)-TAUWSHELTER*F1DCOS3(IJ)*CONSTTAU(IJ)*FNC2,0.0_JWRB)
            USTPH(IJ)   = SQRT(TAUL(IJ))
            PHIHF(IJ) = PHIHF(IJ) + FNC2/Y
          ENDDO
          PHIHF(IJ) = F1DCOS2(IJ)*CONSTPHI(IJ) * SQRTZ0OG(IJ)*PHIHF(IJ)
        ENDDO
      ELSE
        DO IJ=IJS,IJL
          DO J=1,JTOT_TAUHF
            Y         = EXP(ZINF(IJ)+REAL(J-1,JWRB)*DELZ(IJ))
            OMEGA     = Y*SQRTGZ0(IJ)
            CM1       = OMEGA*GM1
            ZX        = USTPH(IJ)*CM1 + ZALP
            ZARG      = XKAPPA/ZX
            ZLOG      = XLOGGZ0(IJ)+2.0_JWRB*LOG(CM1)+ZARG 
            ZLOG      = MIN(ZLOG,0.0_JWRB)
            ZBETA     = EXP(ZLOG)*ZLOG**4
            ZN        = CONST(IJ)*ZBETA*USTPH(IJ)*Y
            GAMNORMA  = (1.0_JWRB + RN1_RN*ZN)/(1.0_JWRB + ZN)
            FNC2      = ZBETA*WTAUHF(J) * GAMNORMA
            PHIHF(IJ) = PHIHF(IJ) + FNC2/Y
          ENDDO
          PHIHF(IJ) = F1DCOS2(IJ)*CONSTPHI(IJ) * SQRTZ0OG(IJ)*TAUL(IJ)*PHIHF(IJ)*DELZ(IJ)
        ENDDO
      ENDIF
      ENDIF

      IF (LHOOK) CALL DR_HOOK('TAU_PHI_HF',1,ZHOOK_HANDLE)

      END SUBROUTINE TAU_PHI_HF
