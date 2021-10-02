SUBROUTINE GETFRSTWND (CDTWIS, CDTWIE,                 &
 &                     IJS, IJL, IFROMIJ, JFROMIJ,     &
 &                     WVENVI, FF_NOW,                 &
 &                     IREAD, LWCUR, NEMO2WAM,         & 
 &                     LLMORE)

! ----------------------------------------------------------------------

!*    PURPOSE.
!     --------

!       GETFRSTWND PROCESS INITIAL WINDFIELDS.

!**   INTERFACE.
!     ----------  

!       *CALL* *GETFRSTWND (CDTWIS, CDTWIE,
!                      IJS, IJL,  IFROMIJ, JFROMIJ,
!                      WVENVI, FF_NOW,
!                      IREAD, LWCUR, NEMO2WAM,
!                      LLMORE)
!          *CDTWIS*   - DATE OF FIRST WIND FIELD.
!          *CDTWIE*   - DATE OF LAST FIRST WIND FIELD.
!          *IJS:IJL   - ARRAYS DIMENSION
!          *IFROMIJ*  - POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!          *JFROMIJ*  - POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!          *WVENVI*   - WAVE ENVIRONMENT.
!          *FF_NOW*   - DATA STRUCTURE WITH THE CURRENT FORCING FIELDS
!          *IREAD*    - PROCESSOR WHICH WILL ACCESS THE FILE ON DISK
!          *LWCUR*    - LOGICAL INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS
!          *NEMO2WAM* - FIELDS FRON OCEAN MODEL to WAM
!          *LLMORE*   - TRUE IF MORE WIND DATA NEEDED (i.e. call to *NOTIM*)


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : ENVIRONMENT, FORCING_FIELDS, OCEAN2WAVE

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
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL), INTENT(IN) :: IFROMIJ  ,JFROMIJ
      TYPE(ENVIRONMENT), DIMENSION(IJS:IJL), INTENT(INOUT) :: WVENVI
      TYPE(FORCING_FIELDS), DIMENSION(IJS:IJL), INTENT(INOUT) :: FF_NOW
      TYPE(OCEAN2WAVE), DIMENSION(IJS:IJL), INTENT(INOUT) :: NEMO2WAM

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

ASSOCIATE(UCUR => WVENVI%UCUR, &
 &        VCUR => WVENVI%VCUR, &
 &        WSWAVE => FF_NOW%WSWAVE, &
 &        WDWAVE => FF_NOW%WDWAVE, &
 &        UFRIC => FF_NOW%UFRIC, &
 &        Z0M => FF_NOW%Z0M, &
 &        AIRD => FF_NOW%AIRD, &
 &        WSTAR => FF_NOW%WSTAR, &
 &        CICOVER => FF_NOW%CICOVER, &
 &        CITHICK => FF_NOW%CITHICK, &
 &        NEMOCICOVER => NEMO2WAM%NEMOCICOVER, &
 &        NEMOCITHICK => NEMO2WAM%NEMOCITHICK)


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
      CALL GETWND (IJS, IJL, IFROMIJ, JFROMIJ,            &
     &             UCUR, VCUR,                            &
     &             WSWAVE, UFRIC,                         &
     &             WDWAVE,                                &
     &             AIRD, WSTAR,                           &
     &             CICOVER, CITHICK,                      &
     &             CDATEWL, LWNDFILE, LCLOSEWND, IREAD,   &
     &             LWCUR, NEMOCICOVER, NEMOCITHICK,       & 
     &             ICODE_WND)


        CALL GSTATS(1493,0)
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO=IJS,IJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL)
          CALL CDUSTARZ0 (KIJS, KIJL, WSWAVE(KIJS:KIJL), XNLEV,            &
     &                    CD(KIJS), UFRIC(KIJS:KIJL), Z0M(KIJS:KIJL))
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1493,1)

      IF (CDATEWL == CDTWIE) THEN
        LLMORE = .FALSE.
      ELSE
        CALL INCDATE (CDTWIS, IDELWO)
        LLMORE = .TRUE.
      ENDIF

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('GETFRSTWND',1,ZHOOK_HANDLE)

END SUBROUTINE GETFRSTWND
