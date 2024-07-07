! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE OUTCH (KIJS, KIJL, IMETHOD,       & 
 &                ZUST, PCHAR , PCHARHQ, CD, & 
 &                CH)

! ----------------------------------------------------------------------

!**** *OUTCH* - DETERMINES THE NEUTRAL HEAT EXCHANGE COEFFICIENT.

! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWABORT , ONLY : WAM_ABORT
USE YOWPCONS , ONLY : G , EPSUS
USE YOWPHYS  , ONLY : XKAPPA, XNLEV, RNU, RNUM

USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
INTEGER(KIND=JWIM), INTENT(IN) :: IMETHOD                      ! = 0    no sea state effect
                                                               ! = 1    based on Janssen( Tech memo 239, 1997)
                                                               ! = 2    based on Janssen and Bidlot, 2018
                                                               !        Progress in Operational Wave Forecasting, Procedia IUTAM
                                                               !        Volume 26, 2018, Pages 14-29.

REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: ZUST      ! FRICTION VELOCITY
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: PCHAR     ! CHARNOCK
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: PCHARHQ   ! EQUIVALENT CHARNOCK FIELD FOR HEAT AND MOISTURE
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: CD        ! NEUTRAL DRAG COEFFICIENT AT HEIGHT XNLEV
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT):: CH        ! NEUTRAL HEAT EXCHANGE COEFFICIENT AT HEIGHT XNLEV


INTEGER(KIND=JWIM) :: IJ

REAL(KIND=JWRB) :: RNUH
REAL(KIND=JWRB) :: Z0M, Z0W, Z0H, Z0WHQ
REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: ZUST2, ZUSTM1, Z0T, SQRTCDKAPPA

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ROUGNESS LENGTH FUNCTIONS
!     -------- ------ ---------


!     JEAN BIDLOT    E.C.M.W.F.      sea-state dependent latent and sensible transfer coefficients
!                                    after Janssen( Tech memo 239, 1997) 
!     ------------------------------------------------------------------

! AERODYNAMIC ROUGHNESS LENGTH OVER SEAS AND OCEANS:
! --------------------------------------------------
REAL(KIND=JWRB) :: PCHARNOCK ! CHARNOCK PARAMETER
REAL(KIND=JWRB) :: PUST2     ! FRICTION VELOCITY SQUARED  

REAL(KIND=JWRB) :: PZ0SEA    ! ROUGHNESS LENGTH FOR MOMEMTUM OVER SEA DUE TO WAVES
PZ0SEA(G, PCHARNOCK, PUST2) = (PCHARNOCK/G)*PUST2


! ROUGHNESS LENGTH FOR HEAT AND HUMIDITY OVER SEAS AND OCEANS:
! ------------------------------------------------------------
REAL(KIND=JWRB), PARAMETER :: Z0HQMIN=0.0000001_JWRB
REAL(KIND=JWRB) :: PZ0       ! ROUGHNESS LENGTH FOR MOMEMTUM OVER SEA 
REAL(KIND=JWRB) :: PZ0HQ     ! ROUGHNESS LENGTH FOR MOMEMTUM OVER SEA WITHOUT THE OCEAN WAVES CONTRIBUTION
REAL(KIND=JWRB) :: PZN       ! ROUGHNESS LENGTH FOR HEAT OR HUMIDITY OVER SEA 
REAL(KIND=JWRB) :: PZP
REAL(KIND=JWRB) :: PZM

REAL(KIND=JWRB) :: PZPLUS    ! LARGEST ROOT OF Z**2 + (PZN+2*PZ0HQ)*Z + PZN*PZ0HQ = 0 
PZPLUS(PZ0HQ, PZN) = -(PZ0HQ+0.5_JWRB*PZN)+SQRT(PZ0HQ**2+0.25_JWRB*PZN**2) 

REAL(KIND=JWRB) :: PZMINS    ! SMALLEST ROOT OF Z**2 + (PZN+2*PZ0HQ)*Z + PZN*PZ0HQ = 0 
PZMINS(PZ0HQ, PZN) = -(PZ0HQ+0.5_JWRB*PZN)-SQRT(PZ0HQ**2+0.25_JWRB*PZN**2) 

REAL(KIND=JWRB) :: PZZ       ! USEFUL FUNCTION
PZZ(PZ0HQ, PZP, PZM) = ABS(PZP)**((PZ0HQ+PZP)/(PZP-PZM))

REAL(KIND=JWRB) :: PZNSEA    ! ROUGHNESS LENGTH FOR HEAT OR HUMIDITY OVER SEA AND OCEANS 
PZNSEA(PZ0HQ, PZN, PZ0, PUST2) = MAX( PZZ(PZ0HQ, PZPLUS(PZ0HQ,PZN), PZMINS(PZ0HQ,PZN)) * &
                               &      PZZ(PZ0HQ, PZMINS(PZ0HQ,PZN), PZPLUS(PZ0HQ,PZN)),  &
                               &      Z0HQMIN)

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('OUTCH',0,ZHOOK_HANDLE)

RNUH = 0.40_JWRB * RNU ! REDUCED KINEMATIC AIR VISCOSITY FOR HEAT (see IFS documentation 3.2.4) 

DO IJ = KIJS,KIJL
  ZUSTM1(IJ)=1.0_JWRB/MAX(ZUST(IJ),EPSUS)
  ZUST2(IJ)=ZUST(IJ)**2
  SQRTCDKAPPA(IJ)=SQRT(CD(IJ))*XKAPPA 
ENDDO

SELECT CASE (IMETHOD)
CASE(0)
  DO IJ = KIJS,KIJL
    Z0T(IJ) = RNUH*ZUSTM1(IJ)
  ENDDO

CASE(1)
  DO IJ = KIJS,KIJL
    Z0M=RNUM*ZUSTM1(IJ)
    Z0W=PZ0SEA(G,PCHAR(IJ),ZUST2(IJ))
    Z0H=RNUH*ZUSTM1(IJ)
    Z0WHQ=PZ0SEA(G,PCHARHQ(IJ),ZUST2(IJ))
    Z0T(IJ)=PZNSEA(Z0WHQ,Z0H,Z0W,ZUST2(IJ))
  ENDDO

CASE(2)
  DO IJ = KIJS,KIJL
    Z0M=RNUM*ZUSTM1(IJ)
    Z0W=PZ0SEA(G,PCHAR(IJ),ZUST2(IJ))
    Z0H=RNUH*ZUSTM1(IJ)
    Z0T(IJ)=SQRT(Z0H*(Z0M+Z0W))
  ENDDO

CASE DEFAULT
  CALL WAM_ABORT('Invalid value of IMETHOD variable.', &
  &              __FILENAME__, __LINE__)

END SELECT

DO IJ = KIJS,KIJL
  CH(IJ)= SQRTCDKAPPA(IJ) / LOG(1.0_JWRB + XNLEV/Z0T(IJ))
ENDDO

IF (LHOOK) CALL DR_HOOK('OUTCH',1,ZHOOK_HANDLE)
! ----------------------------------------------------------------------
END SUBROUTINE OUTCH
