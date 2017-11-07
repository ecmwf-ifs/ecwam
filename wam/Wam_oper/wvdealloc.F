      SUBROUTINE WVDEALLOC

! ----------------------------------------------------------------------

!**** *WVDEALLOC* - WAVE MODEL DEALLOCATION 

!     J. BIDLOT     ECMWF   JANUARY 1997 ATMOSPHERIC COUPLING

!     MODIFICATION.
!     -------------
!     S. ABDALLA    ECMWF   OCTOBER 2001 INCLUSION OF AIR DENSITY & Zi/L


! ----------------------------------------------------------------------

      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,THQ      ,PHIEPS   ,
     &                      PHIAW    ,TAUOC
      USE YOWSPEC  , ONLY : U10NEW   ,THWNEW   ,USNEW    ,Z0NEW    ,
     &            ROAIRN   ,ZIDLNEW  ,FL3
      USE YOWWIND  , ONLY : CDTNEXT  ,FF_NEXT
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('WVDEALLOC',0,ZHOOK_HANDLE)
#endif

!     1.  DEALLOCATE NECESSARY ARRAYS
!         -------------------------
!      IF(ALLOCATED(U10NEW)) DEALLOCATE(U10NEW)
!      IF(ALLOCATED(THWNEW)) DEALLOCATE(THWNEW)
!      IF(ALLOCATED(USNEW)) DEALLOCATE(USNEW)
!      IF(ALLOCATED(Z0NEW)) DEALLOCATE(Z0NEW)

!      IF(ALLOCATED(ROAIRN)) DEALLOCATE(ROAIRN)
!      IF(ALLOCATED(ZIDLNEW)) DEALLOCATE(ZIDLNEW)

      IF(ALLOCATED(FL3)) DEALLOCATE(FL3)


      IF(ALLOCATED(CDTNEXT)) DEALLOCATE(CDTNEXT)
      IF(ALLOCATED(FF_NEXT)) DEALLOCATE(FF_NEXT)

      IF(ALLOCATED(EMEAN)) DEALLOCATE(EMEAN)
      IF(ALLOCATED(FMEAN)) DEALLOCATE(FMEAN)
      IF(ALLOCATED(THQ)) DEALLOCATE(THQ)
      IF(ALLOCATED(PHIEPS)) DEALLOCATE(PHIEPS)
      IF(ALLOCATED(PHIAW)) DEALLOCATE(PHIAW)
      IF(ALLOCATED(TAUOC)) DEALLOCATE(TAUOC)

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('WVDEALLOC',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE WVDEALLOC
