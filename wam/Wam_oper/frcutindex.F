      SUBROUTINE FRCUTINDEX (IJS, IJL, FM, FMWS, CICVR, MIJ)

! ----------------------------------------------------------------------

!**** *FRCUTINDEX* - RETURNS THE LAST FREQUENCY INDEX OF
!                    PROGNOSTIC PART OF SPECTRUM.

!**   INTERFACE.
!     ----------

!       *CALL* *FRCUTINDEX (IJS, IJL, FM, FMWS, CICVR, MIJ)
!          *IJS*    - INDEX OF FIRST GRIDPOINT
!          *IJL*    - INDEX OF LAST GRIDPOINT
!          *FM*     - MEAN FREQUENCY
!          *FMWS*   - MEAN FREQUENCY OF WINDSEA
!          *CICVR*  - CICVR 
!          *MIJ*    - LAST FREQUENCY INDEX

!     METHOD.
!     -------

!*    COMPUTES LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!*    FREQUENCIES LE 2.5*MAX(FMWS,FM).

!     EXTERNALS.
!     ---------

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------

      USE YOWICE   , ONLY : CITHRSH_TAIL
      USE YOWFRED  , ONLY : FR       ,FRATIO   ,FLOGSPRDM1
      USE YOWPARAM , ONLY : NFRE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER:: IJ, IJS, IJL
      INTEGER :: MIJ(IJS:IJL)

      REAL :: TAILFACTOR, FPMH, FPM4
      REAL, DIMENSION(IJS:IJL) :: FM, FMWS, CICVR

      REAL ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('FRCUTINDEX',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

!*    COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!*    FREQUENCIES LE 2.5*MAX(FMWS,FM).
!     ------------------------------------------------------------

      TAILFACTOR=2.5
      FPMH = TAILFACTOR/FR(1)
     
      DO IJ=IJS,IJL
        IF (CICVR(IJ).LE.CITHRSH_TAIL) THEN
          FPM4 = MAX(FMWS(IJ),FM(IJ))*FPMH
          MIJ(IJ) = INT(LOG10(FPM4)*FLOGSPRDM1)+1
          MIJ(IJ) = MIN(MAX(1,MIJ(IJ)),NFRE)
        ELSE
          MIJ(IJ) = NFRE
        ENDIF
      ENDDO

      IF (LHOOK) CALL DR_HOOK('FRCUTINDEX',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE FRCUTINDEX
