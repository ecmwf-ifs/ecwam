      SUBROUTINE TAU_PHI_HF(IJS, IJL, MIJ, USTAR, Z0, XLEVTAIL,
     &                      TAUHF, PHIHF)

! ----------------------------------------------------------------------

!**** *TAU_PHI_HF* - COMPUTATION OF HIGH-FREQUENCY STRESS.
!                                   HIGH-FREQUENCY ENERGY FLUX.

!     PETER A.E.M. JANSSEN    KNMI      OCTOBER 90
!     JEAN BIDLOT  ECMWF  JANUARY 2017

!*    PURPOSE.
!     ---------

!       COMPUTE HIGH-FREQUENCY WAVE STRESS AND ENERGY FLUX
!       FOR BOTH ECMWF PHYSICS AND METREO FRANCE PHYSICS.

!**   INTERFACE.
!     ----------

!       *CALL* *TAU_PHI_HF(IJS, IJL, MIJ, USTAR, Z0, XLEVTAIL,
!                          TAUHF, PHIHF)
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *MIJ* - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!          *USTAR* FRICTION VELOCITY
!          *Z0*    ROUGHNESS LENGTH 
!          *XLEVTAIL* TAIL LEVEL (METEO FRANCE PHYSICS)
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

      USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,ALPHA    ,XKAPPA,
     &               X0TAUHF, JTOT_TAUHF, WTAUHF
      USE YOWFRED  , ONLY : FR
      USE YOWPCONS , ONLY : G        ,ZPI
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL
      INTEGER, INTENT(IN) :: MIJ(IJS:IJL)
      REAL,DIMENSION(IJS:IJL), INTENT(IN) :: USTAR, Z0, XLEVTAIL
      REAL,DIMENSION(IJS:IJL), INTENT(OUT) :: TAUHF, PHIHF

      INTEGER :: J, IJ

      REAL, PARAMETER :: ZSUP = 0.  !  LOG(1.)
      REAL :: OMEGA, OMEGAC, OMEGACC
      REAL :: X0G, UST, UST0, TAUW, TAUW0, DELZTAUW0M1
      REAL :: YC, Y, CM1, ZX, ZARG, ZLOG, ZBETA
      REAL :: DELZ, ZINF
      REAL :: FNC, FNC2, SQRTZ0OG, SQRTGZ0, GM1, GZ0, XLOGGZ0
      REAL :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('TAU_PHI_HF',0,ZHOOK_HANDLE)
#endif

      GM1= 1./G

!     See INIT_X0TAUHF
      X0G=X0TAUHF*G

!*    COMPUTE THE INTEGRALS 
!     ---------------------

      DO IJ=IJS,IJL
        OMEGAC    = ZPI*FR(MIJ(IJ))
        UST0      = USTAR(IJ)
        TAUW0     = UST0**2
        GZ0       = G*Z0(IJ)
        XLOGGZ0   = LOG(GZ0)
        OMEGACC   = MAX(OMEGAC,X0G/UST0)
        SQRTZ0OG  = SQRT(Z0(IJ)*GM1)
        SQRTGZ0   = 1./SQRTZ0OG
        YC        = OMEGACC*SQRTZ0OG
        ZINF      = LOG(YC)
        DELZ      = MAX((ZSUP-ZINF)/REAL(JTOT_TAUHF-1),0.)

        TAUHF(IJ)= 0.
        PHIHF(IJ)= 0.

        TAUW     = TAUW0
        DELZTAUW0M1  = DELZ/TAUW0
        UST      = UST0
        ! Intergrals are integrated following a change of variable : Z=LOG(Y)
        DO J=1,JTOT_TAUHF
          Y         = EXP(ZINF+REAL(J-1)*DELZ)
          OMEGA     = Y*SQRTGZ0
          CM1       = OMEGA*GM1
          ZX        = UST*CM1 +ZALP
          ZARG      = XKAPPA/ZX
          ZLOG      = XLOGGZ0+2.*LOG(CM1)+ZARG 
          ZLOG      = MIN(ZLOG,0.)
          ZBETA     = EXP(ZLOG)*ZLOG**4
          FNC2      = WTAUHF(J)*ZBETA*TAUW*DELZTAUW0M1
          TAUW      = MAX(TAUW-XLEVTAIL(IJ)*FNC2*TAUW0,0.)
          UST       = SQRT(TAUW)
          TAUHF(IJ) = TAUHF(IJ) + FNC2
          PHIHF(IJ) = PHIHF(IJ) + FNC2/Y
        ENDDO
        PHIHF(IJ) = SQRTZ0OG*PHIHF(IJ)

      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('TAU_PHI_HF',1,ZHOOK_HANDLE)
#endif
      END SUBROUTINE TAU_PHI_HF
