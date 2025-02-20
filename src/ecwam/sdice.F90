! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SDICE  (KIJS, KIJL, FL1, FLD, SL, SLICE,     &
     &                  WAVNUM, CGROUP,                       &
     &                  CICV,CITH, ALPFAC)
! ----------------------------------------------------------------------

!**** *SDICE* - CALLING OF DIFFERENT SEA ICE ATTENUATION TERMS


!     JOSH KOUSAL       ECMWF 2023


!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!       *CALL* *SDICE (KIJS, KIJL, FL1, FLD,SL,*
!                       WAVNUM, CGROUP,
!                       CICV,CITH, ALPFAC)*
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *FL1*    - SPECTRUM.
!          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*     - TOTAL SOURCE FUNCTION ARRAY
!          *SLICE*  - TOTAL SOURCE FUNCTION ARRAY, ICE
!          *INDEP*  - DEPTH INDEX
!          *WAVNUM* - WAVE NUMBER
!          *CGROUP* - GROUP SPEED
!          *CICV*   - SEA ICE COVER
!          *CITH*   - SEA ICE THICKNESS
!          *ALPFAC* - FACTOR TO REDUCE ATTENUATION IN ICE. EQUAL TO 1 IN SOLID ICE

!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!     Jie Yu * , W. Erick Rogers, David W. Wang, 2022 
!     DOI:10.1016/j.coldregions.2022.103582

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG    ,NFRE

      USE YOWICE   , ONLY : LCIWA1   ,LCIWA2    ,LCIWA3

      USE YOMHOOK  , ONLY : LHOOK   ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "sdice1.intfb.h"
#include "sdice2.intfb.h"
#include "sdice3.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: FLD, SL
      REAL(KIND=JWRB), DIMENSION(KIJS,NANG,NFRE), INTENT(OUT) ::        SLICE
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: WAVNUM, CGROUP
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: CICV
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: CITH
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: ALPFAC
      
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SDICE',0,ZHOOK_HANDLE)

!     Attenuation of waves in ice (type 1): scattering
      IF(LCIWA1) CALL SDICE1 (KIJS, KIJL, FL1, FLD, SL, SLICE, WAVNUM, CGROUP, CICV, CITH)

!     Attenuation of waves in ice (type 2): bottom friction
      IF(LCIWA2) CALL SDICE2 (KIJS, KIJL, FL1, FLD, SL, SLICE, WAVNUM, CGROUP, CICV      ) 

!     Attenuation of waves in ice (type 3): viscous friction
      IF(LCIWA3) CALL SDICE3 (KIJS, KIJL, FL1, FLD, SL, SLICE, WAVNUM, CGROUP, CICV, CITH, ALPFAC)
      
      IF (LHOOK) CALL DR_HOOK('SDICE',1,ZHOOK_HANDLE)

      END SUBROUTINE SDICE
