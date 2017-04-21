      SUBROUTINE Z0WAVE (IJS, IJL, US, TAUW, Z0)

! ----------------------------------------------------------------------

!**** *Z0WAVE* - DETERMINE THE SEA STATE DEPENDENT ROUGHNESS LENGTH.

!*    PURPOSE.
!     --------

!       COMPUTE ROUGHNESS LENGTH. 

!**   INTERFACE.
!     ----------

!       *CALL* *Z0WAVE (IJS, IJL, US, TAUW, Z0)
!          *IJS*  - INDEX OF FIRST GRIDPOINT.
!          *IJL*  - INDEX OF LAST GRIDPOINT.
!          *US*   - OUTPUT BLOCK OF SURFACE STRESSES.
!          *TAUW* - INPUT BLOCK OF WAVE STRESSES.
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

      USE YOWCOUP  , ONLY : ALPHA
      USE YOWPCONS , ONLY : G
      USE YOWTABL  , ONLY : EPS1
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL
      REAL,DIMENSION(IJS:IJL),INTENT(IN)  ::  US, TAUW
      REAL,DIMENSION(IJS:IJL),INTENT(OUT) ::  Z0

      INTEGER :: IJ
      REAL :: UST2, UST3, ARG, ALPHAOG
      REAL :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('Z0WAVE',0,ZHOOK_HANDLE)
#endif

      ALPHAOG=ALPHA/G

      DO IJ=IJS,IJL
        UST2 = US(IJ)**2
        UST3 = US(IJ)**3
        ARG = MAX(UST2-TAUW(IJ),EPS1)
        Z0(IJ) = ALPHAOG*UST3/SQRT(ARG)
      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('Z0WAVE',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE Z0WAVE
