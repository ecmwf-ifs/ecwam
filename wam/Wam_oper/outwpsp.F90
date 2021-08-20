SUBROUTINE OUTWPSP (IJSLOC, IJLLOC, IJ_OFFSET, FL1, WSWAVE, WDWAVE, UFRIC)
! ----------------------------------------------------------------------

!**** *OUTWPSP* - MODEL OUTPUT OF SPECTRA AT GIVEN LOCATIONS 

!*    PURPOSE.
!     --------

!       CONTROL OUTPUT OF WAVE AND WIND FIELDS.


!**   INTERFACE.
!     ----------
!      *CALL*OUTWPSP (IJSLOC, IJLLOC, IJ_OFFSET, FL1, WSWAVE, WDWAVE, UFRIC)
!      *IJSLOC* - INDEX OF FIRST LOCAL GRIDPOINT
!      *IJLLOC* - INDEX OF LAST LOCAL GRIDPOINT
!      *IJ_OFFSET* OFFSET to point IJSLOC and IJLLOC to the global block of data
!                   only meaningful if unstructured grid
!      *FL1*    - INPUT SPECTRUM.
!      *WSWAVE* - WIND SPEED
!      *WDWAVE* - WIND DIRECTION 
!      *UFRIC*  - FRICTION VELOCITY

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

      USE YOWPARAM , ONLY : NANG     ,NFRE     ,CLDOMAIN
      USE YOWSTAT  , ONLY : CDATEA   ,CDTPRO   ,MARSTYPE

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "out_onegrdpt_sp.intfb.h"
#include "outers.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJSLOC, IJLLOC
      INTEGER(KIND=JWIM), INTENT(IN) :: IJ_OFFSET
      REAL(KIND=JWRB), DIMENSION(IJSLOC:IJLLOC,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(IJSLOC:IJLLOC), INTENT(IN) :: WSWAVE, WDWAVE, UFRIC 


      INTEGER(KIND=JWIM) :: IJ
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('OUTWPSP',0,ZHOOK_HANDLE)

!*    OUTPUT OF SPECTRA FOR ONE GRID POINT SIMULATION
!     -----------------------------------------------

      IF (CLDOMAIN == 's') THEN
        DO IJ=IJSLOC, IJLLOC
          CALL OUT_ONEGRDPT_SP(FL1(IJ:IJ,:,:),UFRIC(IJ),CDTPRO)
        ENDDO
      ENDIF

!*    OUTPUT OF SPECTRA FOR SATELLITE COLLOCATION.
!     --------------------------------------------

      IF (MARSTYPE == 'an' .OR. MARSTYPE == 'fg' .OR. MARSTYPE == '4v') THEN
        IF (CDTPRO /= CDATEA) THEN
          CALL OUTERS (FL1, IJSLOC, IJLLOC, CDTPRO, WSWAVE, WDWAVE, UFRIC)
        ENDIF
      ENDIF

      IF (LHOOK) CALL DR_HOOK('OUTWPSP',1,ZHOOK_HANDLE)

END SUBROUTINE OUTWPSP
