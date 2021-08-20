SUBROUTINE NOTIM (CDTWIS, CDTWIE,                       &
&                 IJS, IJL, FF_NEXT,                    &
&                 IREAD, LWCUR)

! ----------------------------------------------------------------------

!**** *NOTIM* - STEERING MODULE IF NO TIME INTERPOLATION WANTED.

!*    PURPOSE.
!     --------

!       NOTIM NO TIME INTERPOLATION: PROCESS FORCING FIELDS.

!**   INTERFACE.
!     ----------  

!       *CALL* *NOTIM (CDTWIS, CDTWIE,
!    &                 IJS, IJL,
!                      IREAD, LWCUR)
!          *CDTWIS*   - DATE OF FIRST WIND FIELD.
!          *CDTWIE*   - DATE OF LAST FIRST WIND FIELD.
!          *IJS:IJL   - ARRAYS DIMENSION
!          *IREAD*  - PROCESSOR WHICH WILL ACCESS THE FILE ON DISK
!          *LWCUR*  - LOGICAL INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : FORCING_FIELDS

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWSTAT  , ONLY : IDELPRO  ,IDELWO
      USE YOWPHYS  , ONLY : XNLEV
      USE YOWTEST  , ONLY : IU06
      USE YOWWIND  , ONLY : CDATEWL  ,CDTNEXT

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "getwnd.intfb.h"
#include "incdate.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      TYPE(FORCING_FIELDS), DIMENSION(IJS:IJL), INTENT(OUT) :: FF_NEXT
      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      CHARACTER(LEN=14), INTENT(IN) :: CDTWIS, CDTWIE
      LOGICAL, INTENT(IN) :: LWCUR


      INTEGER(KIND=JWIM) :: ICODE_WND

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      CHARACTER(LEN=14) :: CDTWIH

      LOGICAL :: LWNDFILE, LCLOSEWND

! ----------------------------------------------------------------------

!*    1. INITIALISE WIND REQUEST DATE AND PROCESS FIRST WINDFIELD
!*       IF COLD START.
!        --------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('NOTIM',0,ZHOOK_HANDLE)

ASSOCIATE(U10 => FF_NEXT%WSWAVE, &
 &        US => FF_NEXT%UFRIC, &
 &        THW => FF_NEXT%WDWAVE, &
 &        ADS => FF_NEXT%AIRD, &
 &        WSTAR => FF_NEXT%WSTAR, &
 &        CICR => FF_NEXT%CICOVER, &
 &        CITH => FF_NEXT%CITHICK )


      CDTWIH = CDTWIS

      LCLOSEWND=.FALSE.
      IF (LWCOU) THEN
        LWNDFILE=.FALSE.
      ELSE
        LWNDFILE=.TRUE.
      ENDIF

! ----------------------------------------------------------------------

!*    3. LOOP OVER OUTPUT WIND TIMES.
!        ----------------------------

!*    3.1 READ ONE WIND FIELD AND TRANSFORM TO BLOCKS.
!         --------------------------------------------

!*    3.2 SAVE BLOCKED WIND FIELD.
!         ------------------------

      IF (IDELPRO > IDELWO) THEN
!*      WIND FIELD IS NOT CONSTANT FOR ONE PROPAGATION TIME STEP:
!       THIS OPTION IS NO LONGER AVAILABLE !!!
        WRITE (IU06,*) ' '
        WRITE (IU06,*) ' **********************************'
        WRITE (IU06,*) ' *                                *'
        WRITE (IU06,*) ' *   FATAL ERROR IN SUB. NOTIM:   *'
        WRITE (IU06,*) ' *   ==========================   *'
        WRITE (IU06,*) ' *   THE OPTION IDELPRO.GT.IDELWO *'
        WRITE (IU06,*) ' *      IS LONGER AVAILABLE !!!!! *'
        WRITE (IU06,*) ' *                                *'
        WRITE (IU06,*) ' *          PROGRAM ABORTS        *'
        WRITE (IU06,*) ' *                                *'
        WRITE (IU06,*) ' **********************************'
        CALL ABORT1
      ELSE
!*      WIND FIELD IS CONSTANT FOR ONE PROPAGATION TIME STEP:
!       -----------------------------------------------------

        CDTNEXT=CDTWIH

        CALL GETWND (IJS, IJL,                              &
     &               U10, US,                               &
     &               THW,                                   &
     &               ADS, WSTAR,                            &
     &               CICR, CITH,                            &
     &               CDTWIH, LWNDFILE, LCLOSEWND, IREAD,    &
     &               LWCUR, ICODE_WND)


!*      UPDATE WIND FIELD REQUEST TIME.
        CALL INCDATE (CDTWIH,IDELWO)

!*      IF TIME LEFT BRANCH NOT ALLOWED TO LOOP ANYMORE 
        IF (CDTWIH <= CDTWIE) THEN
          WRITE (IU06,*) ' '
          WRITE (IU06,*) ' ********************************************'
          WRITE (IU06,*) ' *                                          *'
          WRITE (IU06,*) ' *       FATAL ERROR IN SUB. NOTIM:         *'
          WRITE (IU06,*) ' *       ==========================         *'
          WRITE (IU06,*) ' * ERROR WHEN WRITTING ON WIND OUTPUT FILE  *'
          WRITE (IU06,*) ' * MORE THAN ONE LOOP NOT ALLOWED ANY LONGER*'
          WRITE (IU06,*) ' * WHEN USING ARRAY INSTEAD OF WRITE TO FILE*'
          WRITE (IU06,*) ' *                                          *'
          WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS         *'
          WRITE (IU06,*) ' *                                          *'
          WRITE (IU06,*) ' ********************************************'
          CALL ABORT1
        ENDIF

      ENDIF

END ASSOCIATE

      IF (LHOOK) CALL DR_HOOK('NOTIM',1,ZHOOK_HANDLE)

END SUBROUTINE NOTIM
