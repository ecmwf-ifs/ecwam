      SUBROUTINE SETWMASK (X, Y, IJS, IJL, THRSCI, IODP)
! ----------------------------------------------------------------------

!**** *SETWMASK* - IMPOSES MASKS
!                  THE ICE MASK IF PRESENT (see LICERUN).
!                  LAND MASK IF IODP(IJ)=0


!*    PURPOSE.
!     --------

!    COPY X TO Y AND IMPOSE THE ICE MASK AND LAND MASK IF PRESENT


!**   INTERFACE.
!     ----------

!       *CALL* *SETWMASK (X, Y, IJS, IJL, THRSCI, IODP )
!         *X*      - INPUT BLOCK ARRAY.
!         *Y*      - OUTPUT BLOCK ARRAY.
!         *IJS*    - INDEX OF FIRST GRID POINT.
!         *IJL*    - INDEX OF LAST GRID POINT.
!         *THRSCI* - SEA ICE THRESHOLD FOR GIVEN PARAMETER.
!         *IODP*   - SEA MASK (0 over LAND, 1 OVER SEA).

! ----------------------------------------------------------------------

      USE YOWICE   , ONLY : LICERUN  , CICOVER
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSTAT  , ONLY : LLSOURCE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, PARAMETER :: IG=1
      INTEGER, INTENT(IN) :: IJS, IJL
      INTEGER, DIMENSION(IJS:IJL), INTENT(IN) :: IODP 
      INTEGER :: IJ

      REAL, INTENT(IN) :: THRSCI
      REAL :: ZHOOK_HANDLE
      REAL, DIMENSION(IJS:IJL), INTENT(IN) :: X
      REAL, DIMENSION(IJS:IJL), INTENT(OUT) :: Y

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SETWMASK',0,ZHOOK_HANDLE)
#endif

      IF (LICERUN .AND. LLSOURCE) THEN
        DO IJ = IJS,IJL
          IF (CICOVER(IJ,IG).GT.THRSCI) THEN
             Y(IJ) = ZMISS
          ELSE
             Y(IJ) = X(IJ)*IODP(IJ) + (1-IODP(IJ))*ZMISS
          ENDIF
        ENDDO
      ELSE
        DO IJ = IJS,IJL
          Y(IJ) = X(IJ)*IODP(IJ) + (1-IODP(IJ))*ZMISS
        ENDDO
      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SETWMASK',1,ZHOOK_HANDLE)
#endif
      RETURN
      END SUBROUTINE SETWMASK
