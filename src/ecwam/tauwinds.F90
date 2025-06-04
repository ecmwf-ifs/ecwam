! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

      FUNCTION TAUWINDS(SDENSIG,CINV,DSII) RESULT(TAU_WINDS)

! ----------------------------------------------------------------------------
!
!  1. Purpose :
!
!      Wind stress (tau) computation from wind-momentum-input
!      function which can be obtained from wind-energy-input (Sin).
!
!                            / FRMAX
!      tau = g * rho_water * | Sin(f)/C(f) df
!                            /

!----------------------------------------------------------------------
!
!     INTERFACE VARIABLES.
!     --------------------

!     ORIGIN.
!     ----------
!     Adapted from Babanin Young Donelan & Banner (BYDB) physics 
!     as implemented as ST6 in WAVEWATCH-III 
!     WW3 module:       W3SRC6MD    
!     WW3 subroutine:   TAUWINDS
!     Implementation into ECWAM DECEMBER 2021 by J. Kousal 

! ----------------------------------------------------------------------------
!

        USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
        USE YOWPCONS , ONLY : G        ,ZPI    ,ROWATER
        USE YOMHOOK  , ONLY : LHOOK   ,DR_HOOK, JPHOOK

!----------------------------------------------------------------------

        IMPLICIT NONE

        REAL(KIND=JWRB), INTENT(IN)  :: SDENSIG(:) ! Sin(sigma) in [m2/rad-Hz]
        REAL(KIND=JWRB), INTENT(IN)  :: CINV(:)    ! inverse phase speed
        REAL(KIND=JWRB), INTENT(IN)  :: DSII(:)    ! freq. bandwidths in [radians]

        REAL(KIND=JWRB) :: TAU_WINDS
        REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------------
!

        IF (LHOOK) CALL DR_HOOK('TAUWINDS',0,ZHOOK_HANDLE)

        TAU_WINDS = G * ROWATER * SUM(SDENSIG*CINV*DSII)

        IF (LHOOK) CALL DR_HOOK('TAUWINDS',1,ZHOOK_HANDLE)

        END FUNCTION TAUWINDS
