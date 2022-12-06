      SUBROUTINE WAMCUR (NXS, NXE, NYS, NYE, FIELDG,    &
     &                   KIJS, KIJL, IFROMIJ, JFROMIJ,  &
     &                   U, V)

! ----------------------------------------------------------------------

!**** *WAMCUR* - TRANSFORMS INPUT CURRENTS TO BLOCKED WAM POINTS.

!     J. BIDLOT  AUGUST 2008 :: REDEFINE WAMCUR TO BLOCK TRANSFORM.

!*    PURPOSE.
!     --------

!       CONVERTS INPUT CURRENT FIELDS TO WAM BLOCKS FOR ALL
!       POINTS IN THE GRID ON A PE !!! INCLUDING OVER THE LOCALLY
!       DEFINED HALO !!!!.

!**   INTERFACE.
!     ----------

!       *CALL WAMCUR (NXS, NXE, NYS, NYE, FIELDG,
!                     KIJS, KIJL, IFROMIJ, JFROMIJ, U, V)*
!          *NXS:NXE*  FIRST DIMENSION OF FIELDG
!          *NYS:NYE*  SECOND DIMENSION OF FIELDG
!          *FIELDG* - INPUT FORCING FIELDS ON THE WAVE MODEL GRID
!          *KIJS:KIJL*  - DIMENSION OF PASSED ARRAYS
!          *IFROMIJ* - POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!          *JFROMIJ* - POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!          *U*       - INTERPOLATED U CURRENT AT ALL POINTS AND BLOCKS.
!          *V*       - INTERPOLATED V CURRENT AT ALL POINTS.

!     METHOD.
!     -------


!     EXTERNALS.
!     ----------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : FORCING_FIELDS

      USE YOWCURR  , ONLY : CURRENT_MAX

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), DIMENSION(NXS:NXE, NYS:NYE), INTENT(IN) :: FIELDG
      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(IN) :: IFROMIJ  ,JFROMIJ
      REAL(KIND=JWRB), DIMENSION (KIJS:KIJL), INTENT(OUT) :: U, V


      INTEGER (KIND=JWIM):: IJ, IX, JY
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

!*    1. TRANSFORM GRIDDED WIND INPUT INTO BLOCK
!        ----------------------------------------

      IF (LHOOK) CALL DR_HOOK('WAMCUR',0,ZHOOK_HANDLE)

      DO IJ = KIJS, KIJL
        IX = IFROMIJ(IJ)
        JY = JFROMIJ(IJ)
        U(IJ) = FIELDG(IX,JY)%UCUR
        U(IJ) = SIGN(MIN(ABS(U(IJ)),CURRENT_MAX),U(IJ))
        V(IJ) = FIELDG(IX,JY)%VCUR
        V(IJ) = SIGN(MIN(ABS(V(IJ)),CURRENT_MAX),V(IJ))
      ENDDO

      IF (LHOOK) CALL DR_HOOK('WAMCUR',1,ZHOOK_HANDLE)

      END SUBROUTINE WAMCUR