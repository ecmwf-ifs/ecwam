      SUBROUTINE INIT_X0TAUHF

! ----------------------------------------------------------------------

!**** *INIT_X0TAUHF* -

!*    PURPOSE.
!     ---------

!     INITIALISATION FOR TAU_PHI_HF


!**   INTERFACE.
!     ----------

!       *CALL* *INIT_X0TAUHF

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,ALPHA    ,XKAPPA,       &
     &             X0TAUHF, JTOT_TAUHF, WTAUHF
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: J

      REAL(KIND=JWRB) :: CONST1, X0, FF, F, DF
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('INIT_X0TAUHF',0,ZHOOK_HANDLE)


!*    1. PRELIMINARY CALCULATIONS.
!        -------------------------

!     find lowest limit for integration X0 *(G/USTAR)
!     ALPHA*X0**2*EXP(XKAPPA/(X0+ZALP))=1
      X0=0.005_JWRB
      DO J=1,30
        FF=EXP(XKAPPA/(X0+ZALP))-1.0_JWRB
        F=ALPHA*X0**2*FF-1.0_JWRB
        IF (F.EQ.0.0_JWRB) EXIT
        DF=ALPHA*FF*(2.0_JWRB*X0-XKAPPA*(X0/(X0+ZALP))**2)
        X0=X0-F/DF
      ENDDO
      X0TAUHF=X0

      CONST1 = (BETAMAX/XKAPPA**2)/3.0_JWRB

      ! Simpson Integration weights (JTOT_TAUHF must be odd) !
      WTAUHF(1)=CONST1
      DO J=2,JTOT_TAUHF-1,2
        WTAUHF(J)=4.0_JWRB*CONST1
        WTAUHF(J+1)=2.0_JWRB*CONST1
      ENDDO
      WTAUHF(JTOT_TAUHF)=CONST1


      IF (LHOOK) CALL DR_HOOK('INIT_X0TAUHF',1,ZHOOK_HANDLE)

      END SUBROUTINE INIT_X0TAUHF
