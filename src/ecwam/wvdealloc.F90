! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

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
      USE YOWGRID  , ONLY : NCHNK
      USE YOWSPEC  , ONLY : FF_NOW   ,VARS_4D, MIJ
      USE YOWSHAL  , ONLY : WVPRPT
      USE YOWFRED  , ONLY : WVPRPT_LAND

      USE YOWNEMOFLDS , ONLY : WAM2NEMO, NEMO2WAM

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      INTEGER(KIND=JWIM) :: ICHNK

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WVDEALLOC',0,ZHOOK_HANDLE)

!     1.  DEALLOCATE NECESSARY ARRAYS
!         -------------------------

      IF (INTFLDS%LALLOC) THEN
        CALL INTFLDS%DEALLOC()
      ENDIF

      IF (FF_NEXT%LALLOC) THEN
         CALL FF_NEXT%DEALLOC()
      ENDIF

      IF (FF_NOW%LALLOC) THEN
         CALL FF_NOW%DEALLOC()
      ENDIF

      IF (VARS_4D%LALLOC) THEN
         CALL VARS_4D%DEALLOC()
      ENDIF

      IF (WVPRPT%LALLOC)THEN
         CALL WVPRPT%DEALLOC()
      ENDIF

      IF (WVPRPT_LAND%LALLOC)THEN
         CALL WVPRPT_LAND%DEALLOC()
      ENDIF

      IF (.NOT. LWNEMOCOU) THEN
        IF (WAM2NEMO%LALLOC) THEN
           CALL WAM2NEMO%DEALLOC()
        ENDIF
        IF (NEMO2WAM%LALLOC) THEN
           CALL NEMO2WAM%DEALLOC()
        ENDIF
      ENDIF

      IF(MIJ%LALLOC) CALL MIJ%DEALLOC()

      IF (LHOOK) CALL DR_HOOK('WVDEALLOC',1,ZHOOK_HANDLE)

      END SUBROUTINE WVDEALLOC
