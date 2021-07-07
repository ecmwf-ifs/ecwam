      SUBROUTINE WAMCUR (U, V, IJS, IJL)

! ----------------------------------------------------------------------

!**** *WAMCUR* - TRANSFORMS INPUT CURRENTS TO BLOCKED WAM POINTS.

!     J. BIDLOR  AUGUST 2008 :: REDEFINE WAMCUR TO BLOCk TRANSFORM.

!*    PURPOSE.
!     --------

!       CONVERTS INPUT CURRENT FIELDS TO WAM BLOCKS FOR ALL
!       POINTS IN THE GRID ON A PE !!! INCLUDING OVER THE LOCALLY
!       DEFINED HALO !!!!.

!**   INTERFACE.
!     ----------

!       *CALL WAMCUR (U, V, IJS, IJL)*
!          *U*   - INTERPOLATED U CURRENT AT ALL POINTS AND BLOCKS.
!          *V*   - INTERPOLATED V CURRENT AT ALL POINTS.
!          *IJS*    - INDEX OF FIRST GRIDPOINT
!          *IJL*    - INDEX OF LAST GRIDPOINT

!     METHOD.
!     -------


!     EXTERNALS.
!     ----------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCURR  , ONLY : CURRENT_MAX
      USE YOWMAP   , ONLY : IFROMIJ  ,JFROMIJ
      USE YOWPARAM , ONLY : NGY
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWIND  , ONLY : FIELDG
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), DIMENSION (IJS:IJL), INTENT(OUT) :: U, V

      INTEGER (KIND=JWIM):: IJ, IX, JY
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

!*    1. TRANSFORM GRIDDED WIND INPUT INTO BLOCK
!        ----------------------------------------

      IF (LHOOK) CALL DR_HOOK('WAMCUR',0,ZHOOK_HANDLE)

      DO IJ = IJS,IJL
        IX = IFROMIJ(IJ)
        JY = JFROMIJ(IJ)
        U(IJ) = FIELDG(IX,JY)%UCUR
        U(IJ) = SIGN(MIN(ABS(U(IJ)),CURRENT_MAX),U(IJ))
        V(IJ) = FIELDG(IX,JY)%VCUR
        V(IJ) = SIGN(MIN(ABS(V(IJ)),CURRENT_MAX),V(IJ))
      ENDDO

      IF (LHOOK) CALL DR_HOOK('WAMCUR',1,ZHOOK_HANDLE)

      END SUBROUTINE WAMCUR
