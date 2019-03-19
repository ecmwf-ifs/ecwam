      SUBROUTINE Z0WAVE (IJS, IJL, US, TAUW, UTOP, Z0)

! ----------------------------------------------------------------------

!**** *Z0WAVE* - DETERMINE THE SEA STATE DEPENDENT ROUGHNESS LENGTH.

!*    PURPOSE.
!     --------

!       COMPUTE ROUGHNESS LENGTH. 

!**   INTERFACE.
!     ----------

!       *CALL* *Z0WAVE (IJS, IJL, US, TAUW, UTOP, Z0)
!          *IJS*  - INDEX OF FIRST GRIDPOINT.
!          *IJL*  - INDEX OF LAST GRIDPOINT.
!          *US*   - OUTPUT BLOCK OF SURFACE STRESSES.
!          *TAUW* - INPUT BLOCK OF WAVE STRESSES.
!          *UTOP* - WIND SPEED.
!          *ZO*   - OUTPUT BLOCK OF ROUGHNESS LENGTH.

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ---------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : ALPHA, LLCAPCHNK
      USE YOWPCONS , ONLY : G
      USE YOWTABL  , ONLY : EPS1
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB),DIMENSION(IJS:IJL),INTENT(IN)  ::  US, TAUW, UTOP
      REAL(KIND=JWRB),DIMENSION(IJS:IJL),INTENT(OUT) ::  Z0

      INTEGER(KIND=JWIM) :: IJ
      REAL(KIND=JWRB) :: UST2, UST3, ARG, GM1
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: ALPHAOG

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('Z0WAVE',0,ZHOOK_HANDLE)

      GM1= 1.0_JWRB/G
      IF(LLCAPCHNK) THEN
        DO IJ=IJS,IJL
          ALPHAOG(IJ)= CHNKMIN(UTOP(IJ))*GM1
        ENDDO
      ELSE
        DO IJ=IJS,IJL
          ALPHAOG(IJ)= ALPHA*GM1
        ENDDO
      ENDIF

      DO IJ=IJS,IJL
        UST2 = US(IJ)**2
        UST3 = US(IJ)**3
        ARG = MAX(UST2-TAUW(IJ),EPS1)
        Z0(IJ) = ALPHAOG(IJ)*UST3/SQRT(ARG)
      ENDDO

      IF (LHOOK) CALL DR_HOOK('Z0WAVE',1,ZHOOK_HANDLE)

      END SUBROUTINE Z0WAVE
