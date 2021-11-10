SUBROUTINE OUTWPSP (IJS, IJL, FL1, FF_NOW)
! ----------------------------------------------------------------------

!**** *OUTWPSP* - MODEL OUTPUT OF SPECTRA AT GIVEN LOCATIONS 

!*    PURPOSE.
!     --------

!       CONTROL OUTPUT OF WAVE AND WIND FIELDS.


!**   INTERFACE.
!     ----------
!      *CALL*OUTWPSP (IJS, IJL, FL1, FF_NOW) 
!      *IJS*    - INDEX OF FIRST LOCAL GRIDPOINT = IJS if not unstructured
!      *IJL*    - INDEX OF LAST LOCAL GRIDPOINT = IJL if not unstructured
!      *FL1*    - INPUT SPECTRUM.
!      *FF_NOW* - FORCING FIELDS

!     EXTERNALS.
!     ----------

!       *OUT_ONEGRDPT_SP* - ONE GRID POINT OUTPUT
!       *OUTERS*    - OUTPUT OF SATELLITE COLOCATION SPECTRA.
!   
!     METHOD.
!     -------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : FORCING_FIELDS

      USE YOWPARAM , ONLY : NANG     ,NFRE     ,CLDOMAIN
      USE YOWSTAT  , ONLY : CDATEA   ,CDTPRO   ,MARSTYPE

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "out_onegrdpt_sp.intfb.h"
#include "outers.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: FL1
      TYPE(FORCING_FIELDS), DIMENSION(IJS:IJL), INTENT(IN) :: FF_NOW


      INTEGER(KIND=JWIM) :: IJ
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('OUTWPSP',0,ZHOOK_HANDLE)

ASSOCIATE(WSWAVE => FF_NOW%WSWAVE, &
 &        WDWAVE => FF_NOW%WDWAVE, &
 &        UFRIC => FF_NOW%UFRIC)


!*    OUTPUT OF SPECTRA FOR ONE GRID POINT SIMULATION
!     -----------------------------------------------

      IF (CLDOMAIN == 's') THEN
        DO IJ=IJS, IJL
          CALL OUT_ONEGRDPT_SP(FL1(IJ:IJ,:,:),UFRIC(IJ),CDTPRO)
        ENDDO
      ENDIF

!*    OUTPUT OF SPECTRA FOR SATELLITE COLLOCATION.
!     --------------------------------------------

      IF (MARSTYPE == 'an' .OR. MARSTYPE == 'fg' .OR. MARSTYPE == '4v') THEN
        IF (CDTPRO /= CDATEA) THEN
          CALL OUTERS (FL1, IJS, IJL, CDTPRO, WSWAVE, WDWAVE, UFRIC)
        ENDIF
      ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('OUTWPSP',1,ZHOOK_HANDLE)

END SUBROUTINE OUTWPSP
