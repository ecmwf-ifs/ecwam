! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE FRCUTINDEX (KIJS, KIJL, FM, FMWS, UFRIC, CICOVER,     &
     &                       MIJ, RHOWGDFTH)

! ----------------------------------------------------------------------

!**** *FRCUTINDEX* - RETURNS THE LAST FREQUENCY INDEX OF
!                    PROGNOSTIC PART OF SPECTRUM.

!**   INTERFACE.
!     ----------

!       *CALL* *FRCUTINDEX (KIJS, KIJL, FM, FMWS, CICOVER, MIJ, RHOWGDFTH)
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *FM*     - MEAN FREQUENCY
!          *FMWS*   - MEAN FREQUENCY OF WINDSEA
!          *UFRIC*  - FRICTION VELOCITY IN M/S
!          *CICOVER*- CICOVER 
!          *MIJ*    - LAST FREQUENCY INDEX for imposing high frequency tail
!          *RHOWGDFTH - WATER DENSITY * G * DF * DTHETA
!                       FOR TRAPEZOIDAL INTEGRATION BETWEEN FR(1) and FR(MIJ) 
!                       !!!!!!!!  RHOWGDFTH=0 FOR FR > FR(MIJ)


!     METHOD.
!     -------

!*    COMPUTES LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!*    FREQUENCIES LE 2.5*MAX(FMWS,FM).


!!! be aware that if this is NOT used, for iphys=1, the cumulative dissipation has to be
!!! re-activated (see module yowphys) !!!


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,DFIM       ,FRATIO   ,FLOGSPRDM1, &
     &                ZPIFR,                                            &
     &                DELTH          ,RHOWG_DFIM ,FRIC
      USE YOWICE   , ONLY : CITHRSH_TAIL
      USE YOWPARAM , ONLY : NFRE
      USE YOWPCONS , ONLY : G        ,EPSMIN
      USE YOWPHYS  , ONLY : TAILFACTOR, TAILFACTOR_PM
      USE YOWSTAT  , ONLY : IPHYS

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "frcutindex_default.intfb.h"
#include "frcutindex_zbry.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM), INTENT(OUT) :: MIJ(KIJL)
      REAL(KIND=JWRB),DIMENSION(KIJL), INTENT(IN) :: FM, FMWS, UFRIC, CICOVER
      REAL(KIND=JWRB),DIMENSION(KIJL,NFRE), INTENT(OUT) :: RHOWGDFTH 


      INTEGER(KIND=JWIM) :: IJ, M

      REAL(KIND=JWRB) :: FPMH, FPPM, FM2, FPM, FPM4
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('FRCUTINDEX',0,ZHOOK_HANDLE)

      SELECT CASE (IPHYS)
      CASE(0,1)
         CALL FRCUTINDEX_DEFAULT(KIJS, KIJL, FM, FMWS, UFRIC, CICOVER,     &
         &                       MIJ)
      CASE(2)
         CALL FRCUTINDEX_ZBRY (KIJS, KIJL, FM, UFRIC, CICOVER,     &
         &                       MIJ)
      END SELECT

!     SET RHOWGDFTH
      DO IJ=KIJS,KIJL
        DO M=1,MIJ(IJ)
          RHOWGDFTH(IJ,M) = RHOWG_DFIM(M)
        ENDDO
        IF (MIJ(IJ) /= NFRE) RHOWGDFTH(IJ,MIJ(IJ))=0.5_JWRB*RHOWGDFTH(IJ,MIJ(IJ))
        DO M=MIJ(IJ)+1,NFRE
          RHOWGDFTH(IJ,M) = 0.0_JWRB
        ENDDO
      ENDDO
      
      IF (LHOOK) CALL DR_HOOK('FRCUTINDEX',1,ZHOOK_HANDLE)

      END SUBROUTINE FRCUTINDEX
