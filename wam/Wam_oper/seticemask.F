      SUBROUTINE SETICEMASK (X, Y, IJS, IJL, IG, THRSCI)
! ----------------------------------------------------------------------

!**** *SETICEMASK* - IMPOSES THE ICE MASK IF PRESENT.

!     J. BIDLOT   ECMWF       MARCH 2000 

!*    PURPOSE.
!     --------

!    COPY X TO Y AND IMPOSE THE ICE MASK IF PRESENT


!**   INTERFACE.
!     ----------

!       *CALL* *SETICEMASK (X, Y, IJS, IJL, IG, THRSCI )
!         *X*      - INPUT BLOCK ARRAY.
!         *Y*      - OUTPUT BLOCK ARRAY.
!         *IJS*    - INDEX OF FIRST GRIDPOINT.
!         *IJL*    - INDEX OF LAST GRIDPOINT.
!         *IG*     - BLOCK NUMBER.
!         *THRSCI* - SEA ICE THRESHOLD FOR GIVEN PARAMETER.

! ----------------------------------------------------------------------

      USE YOWICE   , ONLY : LICERUN  , CICOVER
      USE YOWPCONS , ONLY : ZMISS
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IJ, IJS, IJL, IG
      REAL :: THRSCI
      REAL :: ZHOOK_HANDLE
      REAL, DIMENSION(IJS:IJL) :: X, Y

      IF (LHOOK) CALL DR_HOOK('SETICEMASK',0,ZHOOK_HANDLE)

      IG=1
      IF (LICERUN) THEN
        DO IJ = IJS,IJL
          IF (CICOVER(IJ,IG).GT.THRSCI) THEN
             Y(IJ) = ZMISS
          ELSE
             Y(IJ) = X(IJ)
          ENDIF
        ENDDO
      ELSE
        DO IJ = IJS,IJL
          Y(IJ) = X(IJ)
        ENDDO
      ENDIF

      IF (LHOOK) CALL DR_HOOK('SETICEMASK',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE SETICEMASK
