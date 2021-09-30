SUBROUTINE GETFRSTWND (CDTWIS, CDTWIE,                 &
 &                     IJS, IJL,                       &
 &                     UCUR, VCUR,                     &
 &                     WSWAVE, WDWAVE, UFRIC, Z0M,     &
 &                     AIRD, WSTAR, CICOVER, CITHICK,  &
 &                     IREAD, LWCUR, LLMORE)

! ----------------------------------------------------------------------

!*    PURPOSE.
!     --------

!       GETFRSTWND PROCESS INITIAL WINDFIELDS.

!**   INTERFACE.
!     ----------  

!       *CALL* *GETFRSTWND (CDTWIS, CDTWIE,
!    &                 IJS, IJL,
!    &                 WSWAVE,WDWAVE,UFRIC,Z0M,
!    &                 AIRD, WSTAR, CICOVER, CITHICK,
!                      IREAD, LWCUR, LLMORE)
!          *CDTWIS*   - DATE OF FIRST WIND FIELD.
!          *CDTWIE*   - DATE OF LAST FIRST WIND FIELD.
!          *IJS:IJL   - ARRAYS DIMENSION
!        i *UCUR*     - U-COMPONENT OF THE SURFACE CURRENT
!          *VCUR*     - V-COMPONENT OF THE SURFACE CURRENT
!          *WSWAVE*   - WIND SPEED.
!          *WDWAVE*   - WIND DIRECTION (RADIANS).
!          *UFRIC*    - FRICTION VELOCITY.
!          *Z0M*      - ROUGHNESS LENGTH IN M.
!          *AIRD*     - AIR DENSITY IN KG/M3.
!          *WSTAR*    - CONVECTIVE VELOCITY. 
!          *CICOVER*  - SEA ICE COVER.
!          *CITHICK*  - SEA ICE THICKNESS. 
!          *IREAD*    - PROCESSOR WHICH WILL ACCESS THE FILE ON DISK
!          *LWCUR*    - LOGICAL INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS
!          *LLMORE*   - TRUE IF MORE WIND DATA NEEDED (i.e. call to *NOTIM*)


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWSTAT  , ONLY : IDELWO   ,NPROMA_WAM
      USE YOWPHYS  , ONLY : XNLEV
      USE YOWTEST  , ONLY : IU06
      USE YOWWIND  , ONLY : CDATEWL

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "cdustarz0.intfb.h"
#include "getwnd.intfb.h"
#include "incdate.intfb.h"

      CHARACTER(LEN=14), INTENT(INOUT) :: CDTWIS
      CHARACTER(LEN=14), INTENT(IN) :: CDTWIE
      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), DIMENSION (IJS:IJL), INTENT(IN) :: UCUR, VCUR
      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(INOUT) :: WSWAVE, WDWAVE, UFRIC, Z0M
      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(INOUT) :: AIRD, WSTAR, CICOVER, CITHICK
      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      LOGICAL, INTENT(IN) :: LWCUR
      LOGICAL, INTENT(OUT) :: LLMORE 


      INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA
      INTEGER(KIND=JWIM) :: IJ
      INTEGER(KIND=JWIM) :: ICODE_WND

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB),DIMENSION(IJS:IJL) :: CD

      CHARACTER(LEN=14) :: CDTWIH, ZERO

      LOGICAL :: LWNDFILE, LCLOSEWND

! ----------------------------------------------------------------------

!*    1. INITIALISE WIND REQUEST DATE AND PROCESS FIRST WINDFIELD
!*       IF COLD START.
!        --------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('GETFRSTWND',0,ZHOOK_HANDLE)

      ZERO = ' '
      CDTWIH = CDTWIS

      NPROMA=NPROMA_WAM

      LCLOSEWND=.FALSE.
      IF (LWCOU) THEN
        LWNDFILE=.FALSE.
      ELSE
        LWNDFILE=.TRUE.
      ENDIF

      CDATEWL = CDTWIS
      CALL GETWND (IJS, IJL,                             &
     &             UCUR, VCUR,                           &
     &             WSWAVE, UFRIC,                        &
     &             WDWAVE,                               &
     &             AIRD, WSTAR,                          &
     &             CICOVER, CITHICK,                     &
     &             CDATEWL, LWNDFILE, LCLOSEWND, IREAD,  &
     &             LWCUR, ICODE_WND)


        CALL GSTATS(1493,0)
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO=IJS,IJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL)
          CALL CDUSTARZ0 (KIJS, KIJL, WSWAVE(KIJS), XNLEV,            &
     &                    CD(KIJS), UFRIC(KIJS), Z0M(KIJS))
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1493,1)

      IF (CDATEWL == CDTWIE) THEN
        LLMORE = .FALSE.
      ELSE
        CALL INCDATE (CDTWIS, IDELWO)
        LLMORE = .TRUE.
      ENDIF

      IF (LHOOK) CALL DR_HOOK('GETFRSTWND',1,ZHOOK_HANDLE)

END SUBROUTINE GETFRSTWND
