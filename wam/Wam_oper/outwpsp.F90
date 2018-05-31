      SUBROUTINE OUTWPSP (IJSLOC, IJLLOC, IJ_OFFSET, FL1, USTAR)
! ----------------------------------------------------------------------

!**** *OUTWPSP* - MODEL OUTPUT OF SPECTRA AT GIVEN LOCATIONS 

!*    PURPOSE.
!     --------

!       CONTROL OUTPUT OF WAVE AND WIND FIELDS.


!**   INTERFACE.
!     ----------
!      *CALL*OUTWPSP (IJSLOC, IJLLOC, IJ_OFFSET, FL1
!      *IJSLOC* - INDEX OF FIRST LOCAL GRIDPOINT
!      *IJLLOC* - INDEX OF LAST LOCAL GRIDPOINT
!      *IJ_OFFSET* OFFSET to point IJSLOC and IJLLOC to the global block of data
!                   only meaningful if unstructured grid
!      *FL1*    - INPUT SPECTRUM.
!      *USTAR*  - FRICTION VELOCITY

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
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "out_onegrdpt_sp.intfb.h"
#include "outers.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJSLOC, IJLLOC
      INTEGER(KIND=JWIM), INTENT(IN) :: IJ_OFFSET
      REAL(KIND=JWRB), DIMENSION(IJSLOC:IJLLOC,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(IJSLOC:IJLLOC), INTENT(IN) :: USTAR 

      INTEGER(KIND=JWIM) :: IJ
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTWPSP',0,ZHOOK_HANDLE)
#endif

!*    OUTPUT OF SPECTRA FOR ONE GRID POINT SIMULATION
!     -----------------------------------------------

      IF(CLDOMAIN.EQ.'s') THEN
        DO IJ=IJSLOC, IJLLOC
          CALL OUT_ONEGRDPT_SP(FL1(IJ:IJ,:,:),USTAR(IJ),CDTPRO)
        ENDDO
      ENDIF

!*    OUTPUT OF SPECTRA FOR SATELLITE COLLOCATION.
!     --------------------------------------------

      IF (MARSTYPE.EQ.'an'.OR.MARSTYPE.EQ.'fg'.OR.MARSTYPE.EQ.'4v') THEN
        IF (CDTPRO.NE.CDATEA) THEN
          CALL OUTERS (FL1, IJSLOC, IJLLOC, CDTPRO)
          IF (ITEST.GE.3) THEN
            WRITE(IU06,*) '      SUB. OUTWPSP: OUTPUT OF SPECTRA',      &
     &       ' FOR SATELLITE COLLOCATION DONE FOR ', CDTPRO
          ENDIF
        ENDIF
      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTWPSP',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE OUTWPSP
