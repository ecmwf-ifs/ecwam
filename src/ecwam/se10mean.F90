! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SE10MEAN (KIJS, KIJL, FL1, E10)

! ----------------------------------------------------------------------

!**** *SE10MEAN* - COMPUTATION OF SPECTRAL VARIANCE FOR ALL WAVES
!                 WITH PERIOD LARGER THAN 10s

!*    PURPOSE.
!     --------

!       TO COMPUTE TOTAL ENERGY AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *SE10MEAN(KIJS, KIJL, FL1, E10)*
!          *KIJS*    - INDEX OF FIRST GRIDPOINT
!          *KIJL*    - INDEX OF LAST GRIDPOINT
!          *FL1      - SPECTRUM.
!          *E10*     - MEAN ENERGY

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR
      USE YOWPARAM , ONLY : NANG     ,NFRE

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "sebtmean.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJL, NANG, NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(OUT) :: E10

      REAL(KIND=JWRB) :: TEWHMIN, TEWHMAX
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SE10MEAN',0,ZHOOK_HANDLE)

      TEWHMIN = 10.0_JWRB
      TEWHMAX = 1.0_JWRB / FR(1)
      CALL SEBTMEAN (KIJS, KIJL, FL1, TEWHMIN, TEWHMAX, E10)

      IF (LHOOK) CALL DR_HOOK('SE10MEAN',1,ZHOOK_HANDLE)

      END SUBROUTINE SE10MEAN
