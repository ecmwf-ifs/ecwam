      SUBROUTINE TABU_SWELLFT

!**** *TABU_SWELLFT* - FRICTION COEFFICIENTS IN OSCILLATORY BOUNDARY LAYERS

!     FABRICE ARDHUIN  IFREMER  2013

!*    PURPOSE.
!     --------
!     TO ESTIMATE FRICTION COEFFICIENTS IN OSCILLATORY BOUNDARY LAYERS

!**   INTERFACE.
!     ----------

!       *CALL*TABU_SWELLFT*

!     METHOD.
!     -------
!       TABULATION ON KELVIN FUNCTIONS.

!     EXTERNALS.
!     -----------

!     KERKEI  (zeroth order Kelvin function Ker and Kei)

! ----------------------------------------------------------------------
!
      USE YOWTABL  , ONLY : IAB      ,SWELLFT

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, PARAMETER      :: NITER=100
      REAL,    PARAMETER      :: ABMIN = 0.3 
      REAL,    PARAMETER      :: ABMAX = 8., KAPPA = 0.40
!     VARIABLE.   TYPE.     PURPOSE.
!      *NITER*     INTEGER   NUMBER OF ITERATIONS TO OBTAIN TOTAL STRESS
! ----------------------------------------------------------------------
      INTEGER :: I,ITER
      REAL :: DELAB
      REAL :: KER, KEI
      REAL :: ABR,ABRLOG,L10,FACT,FSUBW,FSUBWMEMO,DZETA0,DZETA0MEMO
!AR:DZETA0 was not set
      DZETA0 = 0
!
      DELAB   = (ABMAX-ABMIN)/REAL(IAB)
      L10=ALOG(10.)
      DO I=1,IAB
         ABRLOG=ABMIN+REAL(I)*DELAB
         ABR=EXP(ABRLOG*L10)
         FACT=1/ABR/(21.2*KAPPA)
         FSUBW=0.05
         DO ITER=1,NITER
            FSUBWMEMO=FSUBW
            DZETA0MEMO=DZETA0
            DZETA0=FACT*FSUBW**(-.5)
            CALL KERKEI(2.*SQRT(DZETA0),KER,KEI)
            FSUBW=.08/(KER**2+KEI**2)
            FSUBW=.5*(FSUBWMEMO+FSUBW)
            DZETA0=.5*(DZETA0MEMO+DZETA0)
         ENDDO   
         SWELLFT(I)=FSUBW
      ENDDO

      END SUBROUTINE TABU_SWELLFT
