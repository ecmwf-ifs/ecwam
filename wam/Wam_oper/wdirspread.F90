      SUBROUTINE WDIRSPREAD (F, IJS, IJL, EMEAN, LLPEAKF, WDIRSPRD)

! ----------------------------------------------------------------------

!**** *WDIRSPREAD* - COMPUTATION OF THE DIRECTIONAL SPREAD 

!     JEAN BIDLOT    ECMWF           MARCH 2000 

!*    PURPOSE.
!     --------

!     TO COMPUTE THE DIRECTIONAL SPREAD EITHER FOR THE WHOLE FREQUENCY
!     RANGE OR BASED ON THE PEAK FREQUENCY.

!**   INTERFACE.
!     ----------

!       *CALL* *WDIRSPREAD (F, IJS, IJL, LLPEAKF, WDIRSPRD)*
!          *F*     - SPECTRUM.
!          *IJS*   - INDEX OF FIRST GRIDPOINT.
!          *IJL*   - INDEX OF LAST GRIDPOINT.
!          *EMEAN* - MEAN WAVE ENERGY FOR THE WAVE SYSTEM (INPUT).
!          *LLPEAKF* - TRUE IF THE DIRECTIONAL SPREAD IS BASED ON 
!                      THE PEAK FREQUENCY.
!          *WDIRSPRD* - DIRECTIONAL SPREAD.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!      *SCOSFL* 
!      *PEAKFRI*

!     REFERENCE.
!     ----------

!       ORIGINALLY REQUESTED BY LUIGI CAVALERI. WE WILL TRY TO PRODUCE
!       SOME FORM OF REFERENCE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,WETAIL
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : EPSMIN
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "peakfri.intfb.h"
#include "scosfl.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: EMEAN
      LOGICAL, INTENT(IN) :: LLPEAKF
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: WDIRSPRD

      INTEGER(KIND=JWIM) :: IJ, M
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL) :: IFRINDEX
      REAL(KIND=JWRB) :: COEF_FR, ONE
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: TEMP
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE) :: F1D 

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WDIRSPREAD',0,ZHOOK_HANDLE)

!*    1. INITIALIZE ARRAYS
!        -----------------

      ONE=1.0_JWRB
      COEF_FR=WETAIL*FR(NFRE)

      DO IJ = IJS,IJL
        WDIRSPRD(IJ) = 0.0_JWRB
      ENDDO

      IF(LLPEAKF) THEN
!     COMPUTATION IS BASED ON THE PEAK FREQUENCY
        CALL PEAKFRI (F, IJS, IJL, IFRINDEX, TEMP, F1D)
        CALL SCOSFL (F, IJS, IJL, IFRINDEX, WDIRSPRD)
        DO IJ = IJS,IJL
          IF(TEMP(IJ).GT.0.0_JWRB) THEN
            WDIRSPRD(IJ) = MIN(WDIRSPRD(IJ)/TEMP(IJ),ONE)
          ELSE
            WDIRSPRD(IJ) = ONE 
          ENDIF
        ENDDO

      ELSE
!     COMPUTATION IS BASED ON THE WHOLE FREQUENCY RANGE
        DO M = 1,NFRE
          IFRINDEX=M
          CALL SCOSFL (F, IJS, IJL, IFRINDEX, TEMP)
          DO IJ = IJS,IJL
            WDIRSPRD(IJ) = WDIRSPRD(IJ) + TEMP(IJ)*DFIM(M)
          ENDDO
        ENDDO
   
!       ADD TAIL CONTRIBUTION
!       note that it uses TEMP because the last call to SCOSFL is made
        DO IJ = IJS,IJL
!         the division by delth is needed since dfim contains delth
          WDIRSPRD(IJ) = WDIRSPRD(IJ)/DELTH + TEMP(IJ)*COEF_FR
        ENDDO

!       NORMALISATION
        DO IJ = IJS,IJL
          IF(EMEAN(IJ).GT.EPSMIN) THEN
            WDIRSPRD(IJ) = MIN(WDIRSPRD(IJ)/EMEAN(IJ),ONE)
          ELSE
            WDIRSPRD(IJ) = ONE 
          ENDIF
        ENDDO

      ENDIF

!     COMPUTE THE ACTUAL SPREAD

      DO IJ = IJS,IJL
        WDIRSPRD(IJ) = SQRT(2.0_JWRB*(ONE-WDIRSPRD(IJ)))
      ENDDO

      IF (LHOOK) CALL DR_HOOK('WDIRSPREAD',1,ZHOOK_HANDLE)

      END SUBROUTINE WDIRSPREAD
