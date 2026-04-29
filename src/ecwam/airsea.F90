! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE AIRSEA (KIJS, KIJL, &
&                        HALP, U10, U10DIR, TAUW, TAUWDIR, RNFAC,  &
&                        US, Z0, Z0B, CHRNCK, ICODE_WND, IUSFG)

! ----------------------------------------------------------------------

!**** *AIRSEA* - DETERMINE TOTAL STRESS IN SURFACE LAYER.

!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990
!     JEAN BIDLOT         ECMWF     FEBRUARY 1999 : TAUT is already
!                                                   SQRT(TAUT)
!     JEAN BIDLOT         ECMWF     OCTOBER 2004: QUADRATIC STEP FOR
!                                                 TAUW

!*    PURPOSE.
!     --------

!       COMPUTE TOTAL STRESS.

!**   INTERFACE.
!     ----------

!       *CALL* *AIRSEA (KIJS, KIJL,
!                       U10, U10DIR, TAUW, TAUWDIR,
!                       US, Z0, Z0B, CHRNCK, ICODE_WND, IUSFG,
!                       HALP, RNFAC)*

!          *KIJS*    - INDEX OF FIRST GRIDPOINT.
!          *KIJL*    - INDEX OF LAST GRIDPOINT.
!          *U10*     - WINDSPEED U10.
!          *U10DIR*  - WINDSPEED DIRECTION.
!          *TAUW*    - WAVE STRESS.
!          *TAUWDIR* - WAVE STRESS DIRECTION.
!          *US*      - OUTPUT OR OUTPUT BLOCK OF FRICTION VELOCITY.
!          *Z0*      - OUTPUT BLOCK OF ROUGHNESS LENGTH.
!          *Z0B*     - BACKGROUND ROUGHNESS LENGTH.
!          *CHRNCK*  - CHARNOCK COEFFICIENT
!          *ICODE_WND* SPECIFIES WHICH OF U10 OR US HAS BEEN FILED UPDATED:
!                     U10: ICODE_WND=3 --> US will be updated
!                     US:  ICODE_WND=1 OR 2 --> U10 will be updated
!          *IUSFG*   - IF = 1 THEN USE THE FRICTION VELOCITY (US) AS FIRST GUESS in TAUT_Z0
!                           0 DO NOT USE THE FIELD US 
!          *HALP*    - OPTIONAL 1/2 PHILLIPS PARAMETER (required for JAN branch).
!          *RNFAC*   - OPTIONAL WIND-DEPENDENT GROWTH RENORMALISATION FACTOR (required for JAN branch).


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM, ONLY : NANG     ,NFRE
      USE YOWPHYS,  ONLY : XKAPPA, XNLEV
      USE YOWTEST,  ONLY : IU06
      USE YOWWIND,  ONLY : WSPMIN
      USE YOWSTAT,  ONLY : IPHYS

      USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

#include "abort1.intfb.h"
#include "airsea_jan.intfb.h"
#include "airsea_zbry.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL, ICODE_WND, IUSFG
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT (IN), OPTIONAL :: HALP, RNFAC
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT (IN) :: U10DIR, TAUW, TAUWDIR
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT (INOUT) :: U10, US, CHRNCK
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT (OUT) :: Z0, Z0B 

      INTEGER(KIND=JWIM) :: IJ, I, J

      REAL(KIND=JWRB) :: XI, XJ, DELI1, DELI2, DELJ1, DELJ2, UST2, ARG, SQRTCDM1
      REAL(KIND=JWRB) :: XKAPPAD, XLOGLEV
      REAL(KIND=JWRB) :: XLEV
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK ('AIRSEA', 0, ZHOOK_HANDLE)

      SELECT CASE (IPHYS)
      CASE(0,1)
            IF ((.NOT.PRESENT(HALP)) .OR. (.NOT.PRESENT(RNFAC))) THEN
                  WRITE (IU06, * ) ' ++++++++++++++++++++++++++++++++++++++++++'
                  WRITE (IU06, * ) ' + AIRSEA : HALP/RNFAC REQUIRED FOR JAN   +'
                  WRITE (IU06, * ) ' ++++++++++++++++++++++++++++++++++++++++++'
                  CALL ABORT1
            ENDIF
            CALL AIRSEA_JAN (KIJS, KIJL, &
&                              HALP, U10, U10DIR, TAUW, TAUWDIR, RNFAC,  &
&                              US, Z0, Z0B, CHRNCK, ICODE_WND, IUSFG)
      CASE(2) 
            CALL AIRSEA_ZBRY(KIJS, KIJL, &
            &                U10, U10DIR, TAUW, TAUWDIR,  &
            &                US, Z0, Z0B, CHRNCK, ICODE_WND, IUSFG)
            
      END SELECT
      
      IF (LHOOK) CALL DR_HOOK ('AIRSEA', 1, ZHOOK_HANDLE)

      END SUBROUTINE AIRSEA
