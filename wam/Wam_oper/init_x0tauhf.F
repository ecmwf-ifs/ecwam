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

      USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,ALPHA    ,XKAPPA,
     &             X0TAUHF, JTOT_TAUHF, WTAUHF
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: J

      REAL :: CONST1, X0, FF, F, DF
      REAL :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('INIT_X0TAUHF',0,ZHOOK_HANDLE)
#endif

!*    1. PRELIMINARY CALCULATIONS.
!        -------------------------

!     find lowest limit for integration X0 *(G/USTAR)
!     ALPHA*X0**2*EXP(XKAPPA/(X0+ZALP))=1
      X0=0.005
      DO J=1,30
        FF=EXP(XKAPPA/(X0+ZALP))-1.
        F=ALPHA*X0**2*FF-1.
        IF(F.EQ.0.0) EXIT
        DF=ALPHA*FF*(2*X0-XKAPPA*(X0/(X0+ZALP))**2)
        X0=X0-F/DF
      ENDDO

      CONST1 = (BETAMAX/XKAPPA**2)/3.

      ! Simpson Integration weights (JTOT_TAUHF must be odd) !
      WTAUHF(1)=CONST1
      DO J=2,JTOT_TAUHF-1,2
        WTAUHF(J)=4.*CONST1
        WTAUHF(J+1)=2.*CONST1
      ENDDO
      WTAUHF(JTOT_TAUHF)=CONST1

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('INIT_X0TAUHF',1,ZHOOK_HANDLE)
#endif
      END SUBROUTINE INIT_X0TAUHF
