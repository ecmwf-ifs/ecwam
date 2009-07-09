      SUBROUTINE WAMCUR (U, V, KINF, KSUP, IJS, IJL, IG)

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

!       *CALL WAMCUR (U, V, KINF, KSUP, IJS, IJL, IG)*
!          *U*   - INTERPOLATED U CURRENT AT ALL POINTS AND BLOCKS.
!          *V*   - INTERPOLATED V CURRENT AT ALL POINTS.
!          *KINF-1* - LOWER BOUND ON FIRST DIMENSION OF U AND V 
!          *KSUP* - UPPER BOUND ON FIRST DIMENSION OF U AND V 
!          *IJS*    - INDEX OF FIRST GRIDPOINT
!          *IJL*    - INDEX OF LAST GRIDPOINT
!          *IG*     - BLOCK NUMBER

!     METHOD.
!     -------


!     EXTERNALS.
!     ----------

! ----------------------------------------------------------------------

      USE YOWMAP   , ONLY : IFROMIJ  ,JFROMIJ
      USE YOWPARAM , ONLY : NGY      ,NBLO
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWIND  , ONLY : FIELDG
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IJS, IJL, KINF, KSUP, IG 
      INTEGER :: IJ, IX, JY

      REAL, DIMENSION (KINF-1:KSUP,NBLO) :: U,V
      REAL :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

!*    1. TRANSFORM GRIDDED WIND INPUT INTO BLOCK
!        ----------------------------------------

      IF (LHOOK) CALL DR_HOOK('WAMCUR',0,ZHOOK_HANDLE)

      DO IJ = IJS,IJL
        IX = IFROMIJ(IJ,IG)
        JY = JFROMIJ(IJ,IG)
        U(IJ,IG) = FIELDG(IX,JY)%UCUR
        V(IJ,IG) = FIELDG(IX,JY)%VCUR
      ENDDO

      IF (LHOOK) CALL DR_HOOK('WAMCUR',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE WAMCUR
