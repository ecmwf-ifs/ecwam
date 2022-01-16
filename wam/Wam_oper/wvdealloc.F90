      SUBROUTINE WVDEALLOC

! ----------------------------------------------------------------------

!**** *WVDEALLOC* - WAVE MODEL DEALLOCATION 

!     J. BIDLOT     ECMWF   JANUARY 1997 ATMOSPHERIC COUPLING

!     MODIFICATION.
!     -------------
!     S. ABDALLA    ECMWF   OCTOBER 2001 INCLUSION OF AIR DENSITY & Zi/L


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWNEMOCOU
      USE YOWMEAN  , ONLY : INTFLDS
      USE YOWWIND  , ONLY : FF_NEXT

      USE YOWNEMOFLDS , ONLY : WAM2NEMO, NEMO2WAM

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WVDEALLOC',0,ZHOOK_HANDLE)

!     1.  DEALLOCATE NECESSARY ARRAYS
!         -------------------------

      IF (ALLOCATED(INTFLDS)) DEALLOCATE(INTFLDS)

      IF (ALLOCATED(FF_NEXT)) DEALLOCATE(FF_NEXT)

      IF (.NOT. LWNEMOCOU) THEN
        IF (ALLOCATED(WAM2NEMO)) DEALLOCATE(WAM2NEMO) 
        IF (ALLOCATED(NEMO2WAM)) DEALLOCATE(NEMO2WAM) 
      ENDIF

      IF (LHOOK) CALL DR_HOOK('WVDEALLOC',1,ZHOOK_HANDLE)

      END SUBROUTINE WVDEALLOC
