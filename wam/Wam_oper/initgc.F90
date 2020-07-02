      SUBROUTINE INITGC

! ----------------------------------------------------------------------

!**** *INITGC* -

!*    PURPOSE.
!     ---------

!     INITIALISATION FOR THE GRAVITY-CAPILLARY PARAMTERISATION


!**   INTERFACE.
!     ----------

!       *CALL* *INITGC

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : NWAV_GC, KRATIO_GC, XKS_GC, XKL_GC,  &
     &                      XK_GC, OMEGA_GC, VG_GC, C_GC, C2OSQRTVG, ZFAK_GC
      USE YOWPCONS , ONLY : G,  SURFT

      USE YOMHOOK  ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: I

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

!     INCLUDE FUNCTIONS FROM GRAVITY-CAPILLARY DISPERSION REALTIONS
#include "gc_dispersion.h"

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('INITGC',0,ZHOOK_HANDLE)

      NWAV_GC = NINT(LOG(XKL_GC/XKS_GC)/(LOG(KRATIO_GC)))

      IF(.NOT.ALLOCATED(XK_GC)) ALLOCATE(XK_GC(NWAV_GC))
      IF(.NOT.ALLOCATED(OMEGA_GC)) ALLOCATE(OMEGA_GC(NWAV_GC))
      IF(.NOT.ALLOCATED(VG_GC)) ALLOCATE(VG_GC(NWAV_GC))
      IF(.NOT.ALLOCATED(C_GC)) ALLOCATE(C_GC(NWAV_GC))
      IF(.NOT.ALLOCATED(C2OSQRTVG)) ALLOCATE(C2OSQRTVG(NWAV_GC))
      IF(.NOT.ALLOCATED(ZFAK_GC)) ALLOCATE(ZFAK_GC(NWAV_GC))

!     COMPUTATION OF WAVENUMBER AND INCREMENT ARRAYS.
!     -------------------------------------------------
      DO I=1,NWAV_GC
        XK_GC(I) = XKS_GC*(KRATIO_GC)**(I-1)
        OMEGA_GC(I) = FOMEG_GC(XK_GC(I))
        VG_GC(I) = FVG_GC(XK_GC(I))
        C_GC(I) = FC_GC(XK_GC(I))      
        C2OSQRTVG(I) = C_GC(I)**2 / SQRT(VG_GC(I))
        ZFAK_GC(I) = OMEGA_GC(I)**2 / (G*XK_GC(I))
      ENDDO

      IF (LHOOK) CALL DR_HOOK('INITGC',1,ZHOOK_HANDLE)

      END SUBROUTINE INITGC
