SUBROUTINE GETFRSTWND (CDTWIS, CDTWIE,                 &
 &                     BLK2LOC, WVENVI, FF_NOW,        &
 &                     IREAD, LWCUR, NEMO2WAM,         & 
 &                     LLMORE)

! ----------------------------------------------------------------------

!*    PURPOSE.
!     --------

!       GETFRSTWND PROCESS INITIAL WINDFIELDS.

!**   INTERFACE.
!     ----------  

!       *CALL* *GETFRSTWND (CDTWIS, CDTWIE,
!                      BLK2LOC, WVENVI, FF_NOW,
!                      IREAD, LWCUR, NEMO2WAM,
!                      LLMORE)
!          *CDTWIS*   - DATE OF FIRST WIND FIELD.
!          *CDTWIE*   - DATE OF LAST FIRST WIND FIELD.
!          *BLK2LOC*  - POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!          *WVENVI*   - WAVE ENVIRONMENT.
!          *FF_NOW*   - DATA STRUCTURE WITH THE CURRENT FORCING FIELDS
!          *IREAD*    - PROCESSOR WHICH WILL ACCESS THE FILE ON DISK
!          *LWCUR*    - LOGICAL INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS
!          *NEMO2WAM* - FIELDS FRON OCEAN MODEL to WAM
!          *LLMORE*   - TRUE IF MORE WIND DATA NEEDED (i.e. call to *NOTIM*)


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDLOC, ENVIRONMENT, FORCING_FIELDS, OCEAN2WAVE

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWSTAT  , ONLY : IDELWO
      USE YOWPHYS  , ONLY : XNLEV
      USE YOWWIND  , ONLY : CDATEWL

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "cdustarz0.intfb.h"
#include "getwnd.intfb.h"
#include "incdate.intfb.h"

      CHARACTER(LEN=14), INTENT(INOUT) :: CDTWIS
      CHARACTER(LEN=14), INTENT(IN) :: CDTWIE
      TYPE(WVGRIDLOC), DIMENSION(NPROMA_WAM, NCHNK), INTENT(IN) :: BLK2LOC
      TYPE(ENVIRONMENT), DIMENSION(NPROMA_WAM, NCHNK), INTENT(INOUT) :: WVENVI
      TYPE(FORCING_FIELDS), DIMENSION(NPROMA_WAM, NCHNK), INTENT(INOUT) :: FF_NOW
      TYPE(OCEAN2WAVE), DIMENSION(NPROMA_WAM, NCHNK), INTENT(INOUT) :: NEMO2WAM

      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      LOGICAL, INTENT(IN) :: LWCUR
      LOGICAL, INTENT(OUT) :: LLMORE 


      INTEGER(KIND=JWIM) :: ICHNK, KIJS, KIJL
      INTEGER(KIND=JWIM) :: IJ
      INTEGER(KIND=JWIM) :: ICODE_WND

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB),DIMENSION(NPROMA_WAM, NCHNK) :: CD

      CHARACTER(LEN=14) :: CDTWIH, ZERO

      LOGICAL :: LWNDFILE, LCLOSEWND

! ----------------------------------------------------------------------

!*    1. INITIALISE WIND REQUEST DATE AND PROCESS FIRST WINDFIELD
!*       IF COLD START.
!        --------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GETFRSTWND',0,ZHOOK_HANDLE)

ASSOCIATE(IFROMIJ => BLK2LOC%IFROMIJ, &
 &        JFROMIJ => BLK2LOC%JFROMIJ, &
 &        UCUR => WVENVI%UCUR, &
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

      LCLOSEWND=.FALSE.
      IF (LWCOU) THEN
        LWNDFILE=.FALSE.
      ELSE
        LWNDFILE=.TRUE.
      ENDIF

      CDATEWL = CDTWIS
      CALL GETWND (IFROMIJ, JFROMIJ,                      &
     &             UCUR, VCUR,                            &
     &             WSWAVE, UFRIC,                         &
     &             WDWAVE,                                &
     &             AIRD, WSTAR,                           &
     &             CICOVER, CITHICK,                      &
     &             CDATEWL, LWNDFILE, LCLOSEWND, IREAD,   &
     &             LWCUR, NEMOCICOVER, NEMOCITHICK,       & 
     &             ICODE_WND)


        CALL GSTATS(1444,0)
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ICHNK, KIJS, KIJL)
        DO ICHNK = 1, NCHNK
          KIJS=1
          KIJL=NPROMA_WAM
          CALL CDUSTARZ0 (KIJS, KIJL, WSWAVE(:,ICHNK), XNLEV,            &
     &                    CD(:,ICHNK), UFRIC(:,ICHNK), Z0M(:,ICHNK))
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1444,1)

      IF (CDATEWL == CDTWIE) THEN
        LLMORE = .FALSE.
      ELSE
        CALL INCDATE (CDTWIS, IDELWO)
        LLMORE = .TRUE.
      ENDIF

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('GETFRSTWND',1,ZHOOK_HANDLE)

END SUBROUTINE GETFRSTWND
