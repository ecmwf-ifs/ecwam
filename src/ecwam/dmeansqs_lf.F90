! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!


      SUBROUTINE DMEANSQS(XKMSS, KIJS, KIJL, F, WAVNUM, MSSXX, MSSYY, MSSXY)

! ----------------------------------------------------------------------

!**** *DMEANSQS* - COMPUTATION OF DIRECTIONAL MEAN SQUARE SLOPE UP TO WAVE
!                  NUMBER XKMSS.

!*    PURPOSE.
!     --------

!       COMPUTE DIRECTIONAL MEAN SQUARE SLOPE COMPONENTS AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *DMEANSQS (XKMSS, KIJS, KIJL, F, WAVNUM, MSSXX, MSSYY, MSSXY)*
!              *XKMSS*      - WAVE NUMBER CUT OFF (M-1).
!              *KIJS*       - INDEX OF FIRST GRIDPOINT.
!              *KIJL*       - INDEX OF LAST GRIDPOINT.
!              *F*          - SPECTRUM.
!              *WAVNUM*     - WAVE NUMBER (M-1), FUNCTION OF FREQUENCY.
!              *MSSXX*      - MEAN SQUARE SLOPE ALONG MODEL X-AXIS, S(TH)^2 (OUTPUT).
!              *MSSYY*      - MEAN SQUARE SLOPE ALONG MODEL Y-AXIS, C(TH)^2 (OUTPUT).
!              *MSSXY*      - CROSS TERM MEAN SQUARE SLOPE, C(TH)*S(TH) (OUTPUT).

!     METHOD.
!     -------

!     THIS ROUTINE IS A VERSION OF MEANSQS.F90 BUT WITH AN ADDITIONAL DIRECTION
!     WEIGHTING AND WITHOUT ANY CONTRIBUTIONS FROM BEYOND THE FR(NFRE)
!     THE WAVE-NUMBER CUT-OFF XKMSS DEFINES A CUT-OFF FREQUENCY FCUT. THE
!     SPECTRAL INTEGRATION IS PERFORMED UP TO THE CORRESPONDING FREQUENCY
!     INDEX. THE COMPONENTS ARE OBTAINED BY PASSING DIRECTIONAL WEIGHTS
!     TO THE INTEGRATION ROUTINE (WANG).

!     EXTERNALS.
!     ----------

!       MEANSQS_LF

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB
      USE YOWPCONS , ONLY : G, ZPI
      USE YOWFRED  , ONLY : FR, TH, FRATIO
      USE YOWPARAM , ONLY : NFRE, NANG
      USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "meansqs_lf.intfb.h"

      REAL(KIND=JWRB), INTENT(IN) :: XKMSS
      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: WAVNUM
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(IN) :: F

      REAL(KIND=JWRB), DIMENSION(KIJL),           INTENT(OUT) :: MSSXX
      REAL(KIND=JWRB), DIMENSION(KIJL),           INTENT(OUT) :: MSSYY
      REAL(KIND=JWRB), DIMENSION(KIJL),           INTENT(OUT) :: MSSXY

      INTEGER(KIND=JWIM) :: NFRE_MSS, NFRE_EFF
      REAL(KIND=JWRB) :: FCUT
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NANG) :: WANG


! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('DMEANSQS',0,ZHOOK_HANDLE)

!     CALCULATE INTEGRATION UPPER LIMIT (IN UNITS OF INDEX)
      IF (XKMSS == 0._JWRB) THEN
        
!       SET UPPER INTEGRATION LIMIT TO MODEL MAXIMUM FREQUENCY
        NFRE_EFF = NFRE
        FCUT = FR(NFRE)

      ELSE

!       CALCULATE INDEX UPTO WHICH TO CALCULATE MSS FROM SPECTRUM
        FCUT = SQRT(G*XKMSS)/ZPI
        NFRE_MSS = INT(LOG(FCUT/FR(1))/LOG(FRATIO))+1 
        NFRE_EFF = MIN(NFRE,NFRE_MSS)

      ENDIF

!     X-COMPONENT
      WANG(:) = SIN(TH(:))*SIN(TH(:))
      CALL MEANSQS_LF(NFRE_EFF, KIJS, KIJL, F, WAVNUM, MSSXX, WANG)

!     Y-COMPONENT
      WANG(:) = COS(TH(:))*COS(TH(:))
      CALL MEANSQS_LF(NFRE_EFF, KIJS, KIJL, F, WAVNUM, MSSYY, WANG)

!     CROSS TERM
      WANG(:) = COS(TH(:))*SIN(TH(:))
      CALL MEANSQS_LF(NFRE_EFF, KIJS, KIJL, F, WAVNUM, MSSXY, WANG)

      IF (LHOOK) CALL DR_HOOK('DMEANSQS',1,ZHOOK_HANDLE)

      END SUBROUTINE DMEANSQS
