! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE AIRSEA_ZBRY (KIJS, KIJL, &
&                        HALP, U10, U10DIR, TAUW, TAUWDIR, RNFAC,  &
&                        US, Z0, Z0B, CHRNCK, ICODE_WND, IUSFG)

! ----------------------------------------------------------------------

!**** *AIRSEA_ZBRY* - DETERMINE TOTAL STRESS IN SURFACE LAYER.

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

!       *CALL* *AIRSEA_ZBRY (KIJS, KIJL, FL1, WAVNUM,
!                       HALP, U10, U10DIR, TAUW, TAUWDIR, RNFAC,
!                       US, Z0, Z0B, CHRNCK, ICODE_WND, IUSFG)*

!          *KIJS*    - INDEX OF FIRST GRIDPOINT.
!          *KIJL*    - INDEX OF LAST GRIDPOINT.
!          *FL1*     - SPECTRA
!          *WAVNUM*  - WAVE NUMBER
!          *HALP*    - 1/2 PHILLIPS PARAMETER
!          *U10*     - WINDSPEED U10.
!          *U10DIR*  - WINDSPEED DIRECTION.
!          *TAUW*    - WAVE STRESS.
!          *TAUWDIR* - WAVE STRESS DIRECTION.
!          *RNFAC*   - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
!          *US*      - OUTPUT OR OUTPUT BLOCK OF FRICTION VELOCITY.
!          *Z0*      - OUTPUT BLOCK OF ROUGHNESS LENGTH.
!          *Z0B*     - BACKGROUND ROUGHNESS LENGTH.
!          *CHRNCK*  - CHARNOCK COEFFICIENT
!          *ICODE_WND* SPECIFIES WHICH OF U10 OR US HAS BEEN FILED UPDATED:
!                     U10: ICODE_WND=3 --> US will be updated
!                     US:  ICODE_WND=1 OR 2 --> U10 will be updated
!          *IUSFG*   - IF = 1 THEN USE THE FRICTION VELOCITY (US) AS FIRST GUESS in TAUT_Z0
!                           0 DO NOT USE THE FIELD US 


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM, ONLY : NANG     ,NFRE
      USE YOWPHYS,  ONLY : XKAPPA, XNLEV, CDFAC
      USE YOWSTAT,  ONLY : IPHYS2_AIRSEA
      USE YOWPCONS, ONLY : G
      USE YOWTEST,  ONLY : IU06
      USE YOWWIND,  ONLY : WSPMIN

      USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

#include "abort1.intfb.h"
#include "taut_z0.intfb.h"
#include "z0wave.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL, ICODE_WND, IUSFG
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT (IN) :: HALP, U10DIR, TAUW, TAUWDIR, RNFAC
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT (INOUT) :: U10, US, CHRNCK
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT (OUT) :: Z0, Z0B

      INTEGER(KIND=JWIM) :: IJ, I, J

      REAL(KIND=JWRB) :: ZNLEV
      REAL(KIND=JWRB), PARAMETER :: RKAP = 0.4_JWRB
      REAL(KIND=JWRB), PARAMETER :: ZRN=1.65E-6_JWRB  ! effective kinematic viscosity (0.11*1.5e-5)

      ! for the ietrative scheme
      INTEGER(KIND=JWIM), PARAMETER :: NITER=15

!      CD=ACD+BCD*U10
      REAL(KIND=JWRB), PARAMETER :: ACD=0.0008_JWRB
      REAL(KIND=JWRB), PARAMETER :: BCD=0.00008_JWRB

!     CD = ACDLIN + BCDLIN*SQRT(PCHAR) * U10
      REAL(KIND=JWRB), PARAMETER :: ACDLIN=0.0008_JWRB
      REAL(KIND=JWRB), PARAMETER :: BCDLIN=0.00047_JWRB
      REAL(KIND=JWRB), PARAMETER :: XEPS=0.00001_JWRB
      REAL(KIND=JWRB), PARAMETER :: USTMIN=0.000001_JWRB
      REAL(KIND=JWRB), PARAMETER :: PCHARMAX=0.1_JWRB
      REAL(KIND=JWRB), PARAMETER :: Z0FG=0.01_JWRB

      INTEGER(KIND=JWIM) :: ITER
      REAL(KIND=JWRB) :: XZNLEV, PCHAROG, XKUTOP, XOLOGZ0
      REAL(KIND=JWRB) :: CDLIN, UST, USTOLD, Z0CH, Z0VIS, F, DELF

      REAL(KIND=JWRB) :: XI, XJ, DELI1, DELI2, DELJ1, DELJ2, UST2, ARG, SQRTCDM1
      REAL(KIND=JWRB) :: XKAPPAD, XLOGLEV
      REAL(KIND=JWRB) :: XLEV, FLX4A0, CD
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK ('AIRSEA_ZBRY', 0, ZHOOK_HANDLE)

!*    2. DETERMINE TOTAL STRESS AND ROUGHNESS (if needed)
!        ----------------------------------

      IF (ICODE_WND == 3) THEN

!       Wind height
        ZNLEV    = 10._JWRB

        !$loki inline
        SELECT CASE (IPHYS2_AIRSEA)
            
        ! implementation of Hwang (2011) as in ST6
        CASE(0)
          FLX4A0 = CDFAC
          DO IJ=KIJS,KIJL
            IF (U10(IJ) .GE. 50.33_JWRB) THEN
                US(IJ) = 2.026_JWRB * SQRT(FLX4A0)
                CD     = (US(IJ)/U10(IJ))**2
            ELSE
                CD     = FLX4A0 * ( 8.058_JWRB + 0.967_JWRB*U10(IJ) - 0.016_JWRB*U10(IJ)**2 ) * 1E-4_JWRB
                US(IJ) = U10(IJ) * SQRT(CD)
            END IF
      !
            Z0(IJ)  = ZNLEV * EXP ( -0.4_JWRB / SQRT(CD) )
          ENDDO

        ! implementation of iterative scheme
        CASE(1,2)

            DO IJ=KIJS,KIJL
                  ! --------------------------------------------
                  ! Iterative method
                  
                  XKUTOP  = RKAP*U10(IJ)

                  ! Start with old charnock (and protect the scheme)
                  PCHAROG = MIN(CHRNCK(IJ),PCHARMAX)/G

                  ! Cd as a linear relation with slope function of Charnock
                  CDLIN= ACDLIN + BCDLIN*SQRT(PCHAROG*G) * U10(IJ)

                  ! first guess for u*
      !          UST = U10(IJ)*SQRT(ACD+BCD*U10(IJ))   ! Use linear approx
      !          UST = SQRT(CD)*U10(IJ)                ! Use Hersbach approx
                  UST = SQRT(CDLIN)*U10(IJ)              ! Use Hersbach approx

                  ! iterate
                  DO ITER=1,NITER
                        USTOLD   = MAX(UST,USTMIN)
                        Z0CH     = PCHAROG*UST**2
                        Z0VIS    = ZRN/UST
                        Z0(IJ)   = Z0CH+Z0VIS
                        XZNLEV    = ZNLEV/(ZNLEV+Z0(IJ))
                        XOLOGZ0  = 1.0_JWRB/LOG(1.0_JWRB+ZNLEV/Z0(IJ))
                        F        = UST-XKUTOP*XOLOGZ0
                        DELF = 1.0_JWRB-XKUTOP*XOLOGZ0**2*XZNLEV*                   &
                  &            (2.0_JWRB*Z0CH-Z0VIS)/(UST*Z0(IJ))
                        IF(DELF /= 0.0_JWRB) UST = UST-F/DELF

                        IF(ABS(UST-USTOLD)<=UST*XEPS .AND. ABS(F)<=XEPS) EXIT
                  ENDDO

                  ! Update Z0, US and then charnock
                  IF(ITER > NITER) THEN
            !           failed to iterate
                        Z0(IJ) = Z0FG
                        US(IJ) = XKUTOP/LOG(1.0+ZNLEV/Z0(IJ))
                  ELSE
                        US(IJ) = MAX(UST,USTMIN)
            !            Z0(IJ) = Z0CH ! Commented out -> Z0=Z0TOT
                  ENDIF
                  
                  
            ENDDO
        END SELECT

      ELSEIF (ICODE_WND == 1 .OR. ICODE_WND == 2) THEN

!*    3. DETERMINE ROUGHNESS LENGTH (if needed).
!        ---------------------------

        !$loki inline
        CALL Z0WAVE (KIJS, KIJL, US, TAUW, U10, Z0, Z0B, CHRNCK)

!*    3. DETERMINE U10 (if needed).
!        ---------------------------

        XKAPPAD = 1.0_JWRB / XKAPPA
        XLOGLEV = LOG (XNLEV)

        DO IJ = KIJS, KIJL
          U10 (IJ) = XKAPPAD * US (IJ) * (XLOGLEV - LOG (Z0 (IJ)))
          U10 (IJ) = MAX (U10 (IJ), WSPMIN)
        ENDDO

      ELSE
        WRITE (IU06, * ) ' ++++++++++++++++++++++++++++++++++++++++++'
        WRITE (IU06, * ) ' + AIRSEA_ZBRY : INVALID VALUE OF ICODE_WND    +'
        WRITE (IU06, * ) ' ICODE_WND = ', ICODE_WND
        WRITE (IU06, * ) ' ++++++++++++++++++++++++++++++++++++++++++'
        CALL ABORT1
      ENDIF

      IF (LHOOK) CALL DR_HOOK ('AIRSEA_ZBRY', 1, ZHOOK_HANDLE)

      END SUBROUTINE AIRSEA_ZBRY
