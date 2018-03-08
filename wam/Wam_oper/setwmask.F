      SUBROUTINE SETWMASK (X, IJS, IJL, CICVR, THRSCI, ICM, IODP)
! ----------------------------------------------------------------------

!**** *SETWMASK* - IMPOSES MASKS
!                  THE ICE MASK IF PRESENT (see LICERUN).
!                  LAND MASK IF IODP(IJ)=0


!*    PURPOSE.
!     --------

!    IMPOSE THE ICE MASK AND LAND MASK IF PRESENT


!**   INTERFACE.
!     ----------

!       *CALL* *SETWMASK (X, IJS, IJL, CICVR, THRSCI, ICM, IODP)
!         *X*      - INPUT BLOCK ARRAY.
!         *IJS*    - INDEX OF FIRST GRID POINT.
!         *IJL*    - INDEX OF LAST GRID POINT.
!         *CICVR*  - SEA ICE COVER FIELD.
!         *THRSCI* - SEA ICE THRESHOLD FOR GIVEN PARAMETER.
!         *ICM*    - IF =1 SEA ICE MASK IS IMPOSED ON OUTPUT FIELD
!         *IODP*   - SEA MASK FIELD (0 over LAND, 1 OVER SEA).

! ----------------------------------------------------------------------

      USE YOWICE   , ONLY : LICERUN
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSTAT  , ONLY : LLSOURCE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL, ICM
      INTEGER, DIMENSION(IJS:IJL), INTENT(IN) :: IODP 
      REAL, DIMENSION(IJS:IJL), INTENT(INOUT) :: X
      REAL, INTENT(IN) :: THRSCI
      REAL, DIMENSION(IJS:IJL), INTENT(IN) :: CICVR 

      INTEGER :: IJ

      REAL :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SETWMASK',0,ZHOOK_HANDLE)
#endif

      IF (LICERUN .AND. LLSOURCE) THEN
        DO IJ = IJS,IJL
          IF (CICVR(IJ).GT.THRSCI) THEN
             X(IJ) = ZMISS
          ELSE
             X(IJ) = X(IJ)*IODP(IJ) + (1-IODP(IJ))*ZMISS
          ENDIF
        ENDDO
      ELSE
        DO IJ = IJS,IJL
          X(IJ) = X(IJ)*IODP(IJ) + (1-IODP(IJ))*ZMISS
        ENDDO
      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SETWMASK',1,ZHOOK_HANDLE)
#endif
      END SUBROUTINE SETWMASK
