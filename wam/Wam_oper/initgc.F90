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
     &                      OMEGA_GC, XK_GC, DELK_GC
      USE YOWPCONS , ONLY : G, SURFT

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

      IF(.NOT.ALLOCATED(OMEGA_GC)) ALLOCATE(OMEGA_GC(NWAV_GC))
      IF(.NOT.ALLOCATED(XK_GC)) ALLOCATE(XK_GC(NWAV_GC))
      IF(.NOT.ALLOCATED(DELK_GC)) ALLOCATE(DELK_GC(NWAV_GC))

!     COMPUTATION OF WAVENUMBER AND INCREMENT ARRAYS.
!     -------------------------------------------------
      DO I=1,NWAV_GC
        XK_GC(I) = XKS_GC*(KRATIO_GC)**(I-1)
        OMEGA_GC(I) = OMEG(XK_GC(I))
      ENDDO

      DELK_GC(1) = 0.5*(XK_GC(2)-XK_GC(1))
      DO I=2,NWAV_GC-1
        DELK_GC(N) = 0.5_JWRB*(XK_GC(I+1)-XK_GC(I-1))
      ENDDO
      DELK_GC(NWAV_GC) = 0.5_JWRB*(XK_GC(NWAV_GC)-XK_GC(NWAV_GC-1))

      IF (LHOOK) CALL DR_HOOK('INITGC',1,ZHOOK_HANDLE)

      END SUBROUTINE INITGC
