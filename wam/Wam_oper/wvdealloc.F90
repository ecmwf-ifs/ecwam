      SUBROUTINE WVDEALLOC

! ----------------------------------------------------------------------

!**** *WVDEALLOC* - WAVE MODEL DEALLOCATION 

!     J. BIDLOT     ECMWF   JANUARY 1997 ATMOSPHERIC COUPLING

!     MODIFICATION.
!     -------------
!     S. ABDALLA    ECMWF   OCTOBER 2001 INCLUSION OF AIR DENSITY & Zi/L


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,PHIEPS   ,              &
     &                      USTOKES  ,VSTOKES  ,STRNMS   ,              &
     &                      PHIAW    ,TAUOC    ,TAUXD    ,TAUYD    ,    &
     &                      TAUOCXD  ,TAUOCYD  ,PHIOCD   ,              &
     &                      WSEMEAN  ,WSFMEAN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WVDEALLOC',0,ZHOOK_HANDLE)

!     1.  DEALLOCATE NECESSARY ARRAYS
!         -------------------------

      IF (ALLOCATED(EMEAN)) DEALLOCATE(EMEAN)
      IF (ALLOCATED(FMEAN)) DEALLOCATE(FMEAN)

      IF (ALLOCATED(USTOKES)) DEALLOCATE(USTOKES)
      IF (ALLOCATED(VSTOKES)) DEALLOCATE(VSTOKES)
      IF (ALLOCATED(STRNMS)) DEALLOCATE(STRNMS)

      IF (ALLOCATED(PHIEPS)) DEALLOCATE(PHIEPS)
      IF (ALLOCATED(PHIAW)) DEALLOCATE(PHIAW)
      IF (ALLOCATED(TAUOC)) DEALLOCATE(TAUOC)
      IF (ALLOCATED(TAUXD)) DEALLOCATE(TAUXD)
      IF (ALLOCATED(TAUYD)) DEALLOCATE(TAUYD)
      IF (ALLOCATED(TAUOCXD)) DEALLOCATE(TAUOCXD)
      IF (ALLOCATED(TAUOCYD)) DEALLOCATE(TAUOCYD)
      IF (ALLOCATED(PHIOCD)) DEALLOCATE(PHIOCD)
      IF (ALLOCATED(WSEMEAN)) DEALLOCATE(WSEMEAN)
      IF (ALLOCATED(WSFMEAN)) DEALLOCATE(WSFMEAN)

      IF (LHOOK) CALL DR_HOOK('WVDEALLOC',1,ZHOOK_HANDLE)

      END SUBROUTINE WVDEALLOC
