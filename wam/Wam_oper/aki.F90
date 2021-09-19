REAL(KIND=JWRB) FUNCTION AKI(OM,BETA)
! ----------------------------------------------------------------------

!**** *AKI* - FUNCTION TO COMPUTE WAVE NUMBER.

!     G. KOMEN, P. JANSSEN   KNMI        01/06/1986

!*    PURPOSE.
!     -------

!       *AKI* COMPUTES THE WAVE NUMBER AS FUNCTION OF
!             CIRCULAR FREQUENCY AND WATER DEPTH.

!**   INTERFACE.
!     ----------

!       *FUNCTION* *AKI (OM, BETA)*
!          *OM*      - CIRCULAR FREQUENCY.
!          *BETA*    - WATER DEPTH.

!     METHOD.
!     -------

!       NEWTONS METHOD TO SOLVE THE DISPERSION RELATION IN SHALLOW
!       WATER.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS , ONLY : G     ,DKMAX

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(IN) :: OM, BETA

!*    *PARAMETER*  RELATIVE ERROR LIMIT OF NEWTON'S METHOD.
      REAL(KIND=JWRB), PARAMETER :: EBS = 0.0001_JWRB

      REAL(KIND=JWRB) :: AKM1, AKM2, AO, AKP, BO, TH, STH  
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('AKI',0,ZHOOK_HANDLE)

!*    1. START VALUE:  MAXIMUM FROM DEEP  AND EXTREM SHALLOW WATER
!                      WAVE NUMBER.
!        ---------------------------------------------------------

      AKM1=OM**2/(4.0_JWRB*G)
      AKM2=OM/(2.0_JWRB*SQRT(G*BETA))
      AO=MAX(AKM1,AKM2)

! ----------------------------------------------------------------------

!*    2. ITERATION LOOP.
!        ---------------

 2000 CONTINUE
      AKP = AO
      BO = BETA*AO
      IF (BO.GT.DKMAX) THEN
        AKI = OM**2/G
      ELSE
        TH = G*AO*TANH(BO)
        STH = SQRT(TH)
        AO = AO+(OM-STH)*STH*2.0_JWRB/(TH/AO+G*BO/COSH(BO)**2)
        IF (ABS(AKP-AO) > EBS*AO) GO TO 2000
        AKI = AO
      ENDIF

      IF (LHOOK) CALL DR_HOOK('AKI',1,ZHOOK_HANDLE)

END FUNCTION AKI
