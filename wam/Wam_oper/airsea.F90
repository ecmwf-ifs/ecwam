      SUBROUTINE AIRSEA (FL1, U10, U10DIR, ROAIRN, TAUW, TAUWDIR, RNFAC, &
&                        US, Z0, IJS, IJL, ICODE_WND, IUSFG)

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

!       *CALL* *AIRSEA (FL1, U10, U10DIR, ROAIRN, TAUW, TAUWDIR, RNFAC,
!                       US, Z0, IJS, IJL, ICODE_WND, IUSFG)*
!          *FL1*  - SPECTRA
!          *U10*  - INPUT OR OUTPUT BLOCK OF WINDSPEED U10.
!          *U10DIR*  - INPUT OR OUTPUT BLOCK OF WINDSPEED DIRECTION.
!          *ROAIRN* - AIR DENSITY IN KG/M3
!          *TAUW* - INPUT BLOCK OF WAVE STRESS.
!          *TAUWDIR* - INPUT BLOCK OF WAVE STRESS DIRECTION.
!          *RNFAC* - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
!          *US*   - OUTPUT OR OUTPUT BLOCK OF FRICTION VELOCITY.
!          *ZO*   - OUTPUT BLOCK OF ROUGHNESS LENGTH.
!          *IJS*  - INDEX OF FIRST GRIDPOINT.
!          *IJL*  - INDEX OF LAST GRIDPOINT.
!          *ICODE_WND* SPECIFIES WHICH OF U10 OR US HAS BEEN FILED UPDATED:
!                     U10: ICODE_WND=3 --> US will be updated
!                     US:  ICODE_WND=1 OR 2 --> U10 will be updated
!          *IUSFG* - IF = 1 THEN USE THE FRICTION VELOCITY (US) AS FIRST GUESS in TAUT_Z0
!                         0 DO NOT USE THE FIELD US 



!     METHOD.
!     -------

!       CALL TAUT_Z0

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ---------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM, ONLY : NANG     ,NFRE
      USE YOWPHYS,  ONLY : XKAPPA, XNLEV
      USE YOWTEST,  ONLY : IU06
      USE YOWWIND,  ONLY : WSPMIN
      USE YOMHOOK,  ONLY : LHOOK, DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "abort1.intfb.h"
#include "taut_z0.intfb.h"
#include "z0wave.intfb.h"

      INTEGER(KIND=JWIM), INTENT (IN) :: IJS, IJL, ICODE_WND, IUSFG

      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION (IJS:IJL), INTENT (IN) :: U10DIR, ROAIRN, TAUW, TAUWDIR, RNFAC
      REAL(KIND=JWRB), DIMENSION (IJS:IJL), INTENT (INOUT) :: U10, US
      REAL(KIND=JWRB), DIMENSION (IJS:IJL), INTENT (OUT) :: Z0

      INTEGER(KIND=JWIM) :: IJ, I, J

      REAL(KIND=JWRB) :: XI, XJ, DELI1, DELI2, DELJ1, DELJ2, UST2, ARG, SQRTCDM1
      REAL(KIND=JWRB) :: XKAPPAD, XLOGLEV
      REAL(KIND=JWRB) :: XLEV
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK ('AIRSEA', 0, ZHOOK_HANDLE)

!*    2. DETERMINE TOTAL STRESS (if needed)
!        ----------------------------------

      IF (ICODE_WND == 3) THEN
        CALL TAUT_Z0 (IJS, IJL, IUSFG, FL1, U10, U10DIR, ROAIRN, TAUW, TAUWDIR, RNFAC, US, Z0)

      ELSEIF (ICODE_WND == 1 .OR. ICODE_WND == 2) THEN

!*    3. DETERMINE ROUGHNESS LENGTH (if needed).
!        ---------------------------

        CALL Z0WAVE (IJS, IJL, US, TAUW, U10, Z0)

!*    3. DETERMINE U10 (if needed).
!        ---------------------------

        XKAPPAD = 1.0_JWRB / XKAPPA
        XLOGLEV = LOG (XNLEV)

        DO IJ = IJS, IJL
          U10 (IJ) = XKAPPAD * US (IJ) * (XLOGLEV - LOG (Z0 (IJ)))
          U10 (IJ) = MAX (U10 (IJ), WSPMIN)
        ENDDO

      ELSE
        WRITE (IU06, * ) ' ++++++++++++++++++++++++++++++++++++++++++'
        WRITE (IU06, * ) ' + AIRSEA : INVALID VALUE OF ICODE_WND    +'
        WRITE (IU06, * ) ' ICODE_WND = ', ICODE_WND
        WRITE (IU06, * ) ' ++++++++++++++++++++++++++++++++++++++++++'
        CALL ABORT1
      ENDIF

      IF (LHOOK) CALL DR_HOOK ('AIRSEA', 1, ZHOOK_HANDLE)

      END SUBROUTINE AIRSEA
