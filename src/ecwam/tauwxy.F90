! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

      SUBROUTINE TAUWXY(SDENSX_IN, SDENSY_IN, CINV_IN, DSII_IN, TAUWX_OUT, TAUWY_OUT)

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
!     Adapted from Babanin Young Donelan & Banner (ZBRY) physics 
!     as implemented as ST6 in WAVEWATCH-III 
!     Implementation into ECWAM DECEMBER 2021 by J. Kousal 

! ----------------------------------------------------------------------------
!

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWFRED  , ONLY : NFRE_EXT
      USE YOWPCONS , ONLY : G, ROWATER
      USE YOMHOOK  , ONLY : LHOOK   ,DR_HOOK, JPHOOK

!----------------------------------------------------------------------


      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(IN)  :: SDENSX_IN(NFRE_EXT), SDENSY_IN(NFRE_EXT)
      REAL(KIND=JWRB), INTENT(IN)  :: CINV_IN(NFRE_EXT), DSII_IN(NFRE_EXT)
      REAL(KIND=JWRB), INTENT(OUT) :: TAUWX_OUT, TAUWY_OUT

      INTEGER(KIND=JWIM) :: M
      REAL(KIND=JWRB)    :: SUMX, SUMY
      REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

!----------------------------------------------------------------------
!

      IF (LHOOK) CALL DR_HOOK('TAUWXY',0,ZHOOK_HANDLE)

      SUMX = 0.0_JWRB
      SUMY = 0.0_JWRB

      DO M = 1, NFRE_EXT
         SUMX = SUMX + SDENSX_IN(M) * CINV_IN(M) * DSII_IN(M)
         SUMY = SUMY + SDENSY_IN(M) * CINV_IN(M) * DSII_IN(M)
      END DO

      TAUWX_OUT = G * ROWATER * SUMX
      TAUWY_OUT = G * ROWATER * SUMY

      IF (LHOOK) CALL DR_HOOK('TAUWXY',1,ZHOOK_HANDLE)

      END SUBROUTINE TAUWXY
