       SUBROUTINE TAU_PHI_HF(KIJS, KIJL, LTAUWSHELTER, UFRIC, Z0M, &
     &                       FL1, WDWAVE, AIRD, RNFAC,             &
     &                       UST, TAUHF, PHIHF, LLPHIHF)

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

!       *CALL* *TAU_PHI_HF(KIJS, KIJL, LTAUWSHELTER, UFRIC, UST, Z0M,
!                          FL1, WDWAVE, AIRD, RNFAC, &
!                          UST, TAUHF, PHIHF, LLPHIHF)
!          *KIJS*         - INDEX OF FIRST GRIDPOINT
!          *KIJL*         - INDEX OF LAST GRIDPOINT
!          *LTAUWSHELTER* - if true then TAUWSHELTER 
!          *FL1*          - WAVE SPECTRUM.
!          *WDWAVE*       - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!          *AIRD*       - AIR DENSITY IN KG/M**3.
!          *RNFAC*        - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
!          *UFRIC*        - FRICTION VELOCITY
!          *UST*          - REDUCED FRICTION VELOCITY DUE TO SHELTERING
!          *Z0M*          - ROUGHNESS LENGTH 
!          *TAUHF*        - HIGH-FREQUENCY STRESS
!          *PHIHF*        - HIGH-FREQUENCY ENERGY FLUX INTO OCEAN
!          *LLPHIHF*      - TRUE IF PHIHF NEEDS TO COMPUTED


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
      USE YOWPHYS  , ONLY : ZALP   , XKAPPA    ,TAUWSHELTER, GAMNCONST
      USE YOWTEST  , ONLY : IU06

      USE YOMHOOK  , ONLY : LHOOK  , DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "omegagc.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      LOGICAL, INTENT(IN) :: LTAUWSHELTER
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: UFRIC, Z0M
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: WDWAVE, AIRD, RNFAC

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: UST
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: TAUHF, PHIHF
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
      REAL(KIND=JWRB) :: ZN1, ZN2, GAMCFR5
      REAL(KIND=JWRB) :: COSW, FCOSW2 
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: SQRTZ0OG, ZSUP, ZINF, DELZ
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TAUL, XLOGGZ0, SQRTGZ0
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: USTPH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: CONST1, CONST2, CONSTTAU, CONSTPHI
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: F1DCOS2, F1DCOS3 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: F1D, F1DSIN2 

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('TAU_PHI_HF',0,ZHOOK_HANDLE)

!     See INIT_X0TAUHF
      X0G=X0TAUHF*G

      IF (LLPHIHF) USTPH(:) = UST(:)

!*    COMPUTE THE INTEGRALS 
!     ---------------------

      DO IJ=KIJS,KIJL
        XLOGGZ0(IJ) = LOG(G*Z0M(IJ))
        OMEGACC = MAX(ZPIFR(NFRE), X0G/UST(IJ))
        SQRTZ0OG(IJ)  = SQRT(Z0M(IJ)*GM1)
        SQRTGZ0(IJ) = 1.0_JWRB/SQRTZ0OG(IJ)
        YC = OMEGACC*SQRTZ0OG(IJ)
        ZINF(IJ) = LOG(YC)
      ENDDO

      DO IJ=KIJS,KIJL
        CONSTTAU(IJ) = ZPI4GM2*FR5(NFRE)
      ENDDO

      K=1
      DO IJ=KIJS,KIJL
        COSW     = MAX(COS(TH(K)-WDWAVE(IJ)),0.0_JWRB)
        FCOSW2   = FL1(IJ,K,NFRE)*COSW**2
        F1DCOS3(IJ) = FCOSW2*COSW
        F1DCOS2(IJ) = FCOSW2
        F1DSIN2(IJ) = FL1(IJ,K,NFRE)*SIN(TH(K)-WDWAVE(IJ))**2
        F1D(IJ) = FL1(IJ,K,NFRE)
      ENDDO
      DO K=2,NANG
        DO IJ=KIJS,KIJL
          COSW     = MAX(COS(TH(K)-WDWAVE(IJ)),0.0_JWRB)
          FCOSW2   = FL1(IJ,K,NFRE)*COSW**2
          F1DCOS3(IJ) = F1DCOS3(IJ) + FCOSW2*COSW
          F1DCOS2(IJ) = F1DCOS2(IJ) + FCOSW2 
          F1DSIN2(IJ) = F1DSIN2(IJ) + FL1(IJ,K,NFRE)*SIN(TH(K)-WDWAVE(IJ))**2
          F1D(IJ) = F1D(IJ) + FL1(IJ,K,NFRE)
        ENDDO
      ENDDO
      DO IJ=KIJS,KIJL
        F1DCOS3(IJ) = DELTH*F1DCOS3(IJ)
        F1DCOS2(IJ) = DELTH*F1DCOS2(IJ)
        F1DSIN2(IJ) = DELTH*F1DSIN2(IJ)
        F1D(IJ) = DELTH*F1D(IJ)
      ENDDO

      IF (LLNORMAGAM) THEN
        GAMCFR5 = GAMNCONST*FR5(NFRE)
        DO IJ=KIJS,KIJL
          CONST1(IJ) = GAMCFR5*RNFAC(IJ)*F1DSIN2(IJ)*SQRTGZ0(IJ)
          CONST2(IJ) = GAMCFR5*RNFAC(IJ)*F1D(IJ)*SQRTGZ0(IJ)
        ENDDO
      ELSE
        CONST1(:) = 0.0_JWRB
        CONST2(:) = 0.0_JWRB
      ENDIF


!     TAUHF :
      IF (LLGCBZ0) THEN
        DO IJ=KIJS,KIJL
          CALL OMEGAGC(UFRIC(IJ), NS, XKS, OMS)
          ZSUP(IJ) = MIN(LOG(OMS*SQRTZ0OG(IJ)),ZSUPMAX)
        ENDDO
      ELSE
        ZSUP(:) = ZSUPMAX
      ENDIF

      DO IJ=KIJS,KIJL
        TAUL(IJ) = UST(IJ)**2
        DELZ(IJ) = MAX((ZSUP(IJ)-ZINF(IJ))/REAL(JTOT_TAUHF-1,JWRB),0.0_JWRB)
        TAUHF(IJ) = 0.0_JWRB
      ENDDO

     ! Intergrals are integrated following a change of variable : Z=LOG(Y)
      IF ( LTAUWSHELTER ) THEN
        DO IJ=KIJS,KIJL
          DO J=1,JTOT_TAUHF
            Y         = EXP(ZINF(IJ)+REAL(J-1,JWRB)*DELZ(IJ))
            OMEGA     = Y*SQRTGZ0(IJ)
            CM1       = OMEGA*GM1
            ZX        = UST(IJ)*CM1 + ZALP
            ZARG      = XKAPPA/ZX
            ZLOG      = XLOGGZ0(IJ)+2.0_JWRB*LOG(CM1)+ZARG 
            ZLOG      = MIN(ZLOG,0.0_JWRB)
            ZBETA     = EXP(ZLOG)*ZLOG**4
            ZN1       = CONST1(IJ)*ZBETA*UST(IJ)*Y
            ZN2       = CONST2(IJ)*ZBETA*UST(IJ)*Y
            GAMNORMA  = (1.0_JWRB + ZN1)/(1.0_JWRB + ZN2)
            FNC2      = F1DCOS3(IJ)*CONSTTAU(IJ)* ZBETA*TAUL(IJ)*WTAUHF(J)*DELZ(IJ) * GAMNORMA
            TAUL(IJ)  = MAX(TAUL(IJ)-TAUWSHELTER*FNC2,0.0_JWRB)

            UST(IJ)   = SQRT(TAUL(IJ))
            TAUHF(IJ) = TAUHF(IJ) + FNC2
          ENDDO
        ENDDO
      ELSE
        DO IJ=KIJS,KIJL
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
            ZN1       = CONST1(IJ)*ZBETA*UST(IJ)*Y
            ZN2       = CONST2(IJ)*ZBETA*UST(IJ)*Y
            GAMNORMA  = (1.0_JWRB + ZN1)/(1.0_JWRB + ZN2)
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
      DO IJ=KIJS,KIJL
        TAUL(IJ) = USTPH(IJ)**2
        ZSUP(IJ) = ZSUPMAX
        DELZ(IJ) = MAX((ZSUP(IJ)-ZINF(IJ))/REAL(JTOT_TAUHF-1,JWRB),0.0_JWRB)
      ENDDO

      DO IJ=KIJS,KIJL
        CONSTPHI(IJ) = AIRD(IJ)*ZPI4GM1*FR5(NFRE)
      ENDDO

     ! Intergrals are integrated following a change of variable : Z=LOG(Y)
      IF ( LTAUWSHELTER ) THEN
        DO IJ=KIJS,KIJL
          DO J=1,JTOT_TAUHF
            Y         = EXP(ZINF(IJ)+REAL(J-1,JWRB)*DELZ(IJ))
            OMEGA     = Y*SQRTGZ0(IJ)
            CM1       = OMEGA*GM1
            ZX        = USTPH(IJ)*CM1 + ZALP
            ZARG      = XKAPPA/ZX
            ZLOG      = XLOGGZ0(IJ)+2.0_JWRB*LOG(CM1)+ZARG 
            ZLOG      = MIN(ZLOG,0.0_JWRB)
            ZBETA     = EXP(ZLOG)*ZLOG**4
            ZN1       = CONST1(IJ)*ZBETA*UST(IJ)*Y
            ZN2       = CONST2(IJ)*ZBETA*UST(IJ)*Y
            GAMNORMA  = (1.0_JWRB + ZN1)/(1.0_JWRB + ZN2)
            FNC2      = ZBETA*TAUL(IJ)*WTAUHF(J)*DELZ(IJ) * GAMNORMA
            TAUL(IJ)  = MAX(TAUL(IJ)-TAUWSHELTER*F1DCOS3(IJ)*CONSTTAU(IJ)*FNC2,0.0_JWRB)
            USTPH(IJ)   = SQRT(TAUL(IJ))
            PHIHF(IJ) = PHIHF(IJ) + FNC2/Y
          ENDDO
          PHIHF(IJ) = F1DCOS2(IJ)*CONSTPHI(IJ) * SQRTZ0OG(IJ)*PHIHF(IJ)
        ENDDO
      ELSE
        DO IJ=KIJS,KIJL
          DO J=1,JTOT_TAUHF
            Y         = EXP(ZINF(IJ)+REAL(J-1,JWRB)*DELZ(IJ))
            OMEGA     = Y*SQRTGZ0(IJ)
            CM1       = OMEGA*GM1
            ZX        = USTPH(IJ)*CM1 + ZALP
            ZARG      = XKAPPA/ZX
            ZLOG      = XLOGGZ0(IJ)+2.0_JWRB*LOG(CM1)+ZARG 
            ZLOG      = MIN(ZLOG,0.0_JWRB)
            ZBETA     = EXP(ZLOG)*ZLOG**4
            ZN1       = CONST1(IJ)*ZBETA*UST(IJ)*Y
            ZN2       = CONST2(IJ)*ZBETA*UST(IJ)*Y
            GAMNORMA  = (1.0_JWRB + ZN1)/(1.0_JWRB + ZN2)
            FNC2      = ZBETA*WTAUHF(J) * GAMNORMA
            PHIHF(IJ) = PHIHF(IJ) + FNC2/Y
          ENDDO
          PHIHF(IJ) = F1DCOS2(IJ)*CONSTPHI(IJ) * SQRTZ0OG(IJ)*TAUL(IJ)*PHIHF(IJ)*DELZ(IJ)
        ENDDO
      ENDIF
      ENDIF

      IF (LHOOK) CALL DR_HOOK('TAU_PHI_HF',1,ZHOOK_HANDLE)

      END SUBROUTINE TAU_PHI_HF
