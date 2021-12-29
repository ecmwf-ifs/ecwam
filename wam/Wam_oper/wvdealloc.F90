      SUBROUTINE WVDEALLOC

! ----------------------------------------------------------------------

!**** *WVDEALLOC* - WAVE MODEL DEALLOCATION 

!     J. BIDLOT     ECMWF   JANUARY 1997 ATMOSPHERIC COUPLING

!     MODIFICATION.
!     -------------
!     S. ABDALLA    ECMWF   OCTOBER 2001 INCLUSION OF AIR DENSITY & Zi/L


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMEAN  , ONLY : INTFLDS

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WVDEALLOC',0,ZHOOK_HANDLE)

!     1.  DEALLOCATE NECESSARY ARRAYS
!         -------------------------

      IF (ALLOCATED(INTFLDS)) DEALLOCATE(INTFLDS)

      IF (LHOOK) CALL DR_HOOK('WVDEALLOC',1,ZHOOK_HANDLE)

      END SUBROUTINE WVDEALLOC
