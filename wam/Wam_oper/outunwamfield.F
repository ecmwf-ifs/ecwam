      SUBROUTINE OUTUNWAMFIELD (ARRAY, IJS, IJL)

! ----------------------------------------------------------------------

!**** *OUTUNWAMFIELD* - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT
!                AND MEAN WAVE NUMBER .
!                THE COMPUTATION OF THE MEAN WAVE ENERGU WAS ALSO
!                ADDED SUCH THAT A CALL TO OUTUNWAMFIELD DOES NOT NEED
!                TO BE PRECEDED BY A CALL TO SEMEAN.

!     S.D. HASSELMANN
!     MODIFIED : P.JANSSEN (INTEGRATION OF F**-4 TAIL)
!     OPTIMIZED BY : L. ZAMBRESKY AND H. GUENTHER


!*    PURPOSE.
!     --------

!       COMPUTE MEAN FREQUENCY AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *OUTUNWAMFIELD (F, IJS, IJL, EM, FM, LLAK)*
!              *F*   - SPECTRUM.
!              *IJS* - INDEX OF FIRST GRIDPOINT
!              *IJL* - INDEX OF LAST GRIDPOINT
!              *EM*  - MEAN WAVE ENERGY 
!              *FM*  - MEAN WAVE FREQUENCY 
!              *LLAK*- TRUE IF MEAN WAVE NUMBER IS COMPUTED

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWUNPOOL, ONLY : FILEDEF, XFN_TEST, DBG 


! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL    :: ARRAY(IJS:IJL)
      INTEGER :: IJ,IJS,IJL

      REAL, SAVE :: XFNOUTTIME
 
      DATA XFNOUTTIME /0./

! ----------------------------------------------------------------------

!*    1. WRITE TEST OUTPUT TO XFN 
!        ------------------------------------------------
      XFNOUTTIME = XFNOUTTIME + 1.

      WRITE(DBG%FHNDL,*) 'WRITING IN OUTUNWAMTEST', IJS, IJL
      WRITE(DBG%FHNDL,*) 'SUM ARRAY', SUM(ARRAY)

      WRITE(XFN_TEST%FHNDL) XFNOUTTIME 
      WRITE(XFN_TEST%FHNDL) (0., 0., ARRAY(IJ), IJ = IJS, IJL)

      RETURN
      END SUBROUTINE OUTUNWAMFIELD
