      SUBROUTINE FRCUTINDEX (IJS, IJL, FM, FMWS, MIJ)

! ----------------------------------------------------------------------

!**** *FRCUTINDEX* - RETURNS THE LAST FREQUENCY INDEX OF
!                    PROGNOSTIC PART OF SPECTRUM.

!**   INTERFACE.
!     ----------

!       *CALL* *FRCUTINDEX (IJS, IJL, FM, FMWS, MIJ)
!          *IJS*    - INDEX OF FIRST GRIDPOINT
!          *IJL*    - INDEX OF LAST GRIDPOINT
!          *FM*     - MEAN FREQUENCY
!          *FMWS*   - MEAN FREQUENCY OF WINDSEA
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

      USE YOWFRED  , ONLY : FR       ,FRATIO   ,FLOGSPRDM1
      USE YOWPARAM , ONLY : NFRE

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER:: IJ, IJS, IJL
      INTEGER :: MIJ(IJS:IJL)

      REAL :: TAILFACTOR, FPMH, FPM4
      REAL, DIMENSION(IJS:IJL) :: FM, FMWS

! ----------------------------------------------------------------------

!*    COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!*    FREQUENCIES LE 2.5*MAX(FMWS,FM).
!     ------------------------------------------------------------

      TAILFACTOR=2.5
      FPMH = TAILFACTOR/FR(1)
     
      DO IJ=IJS,IJL
        FPM4 = MAX(FMWS(IJ),FM(IJ))*FPMH
        MIJ(IJ) = INT(LOG10(FPM4)*FLOGSPRDM1)+1
        MIJ(IJ) = MIN(MIJ(IJ),NFRE)
      ENDDO

      RETURN
      END SUBROUTINE FRCUTINDEX
