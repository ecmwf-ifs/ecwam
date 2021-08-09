      SUBROUTINE NOTIM (CDTWIS, CDTWIE,                       &
     &                  IJS, IJL,                             &
     &                  U10OLD, THWOLD, USOLD, Z0OLD,         &
     &                  ROAIRO, ZIDLOLD, CICOVER, CITHICK,    &
     &                  IREAD, LWCUR)

! ----------------------------------------------------------------------

!     MODIFIED  J. BIDLOT  FEBRARY 1996  MESSAGE PASSING
!               S. ABDALLA OCTOBER 2001  INCLUSION OF AIR DENSITY & Zi/L
!               J BIDLOT AUGUST 2006 INCLUSION OF ICE MASK.
!               J BIDLOT AUGUST 2008 REMOVE INFORMATION FROM COUPLED RUNS
!                                    AND CLEAN UP OBSOLETE OPTIONS.


!**** *NOTIM* - STEERING MODULE IF NO TIME INTERPOLATION WANTED.

!*    PURPOSE.
!     --------

!       NOTIM NO TIME INTERPOLATION: PROCESS WINDFIELDS.

!**   INTERFACE.
!     ----------  

!       *CALL* *NOTIM (CDTWIS, CDTWIE,
!    &                 IJS, IJL,
!    &                 U10OLD,THWOLD,USOLD,Z0OLD,
!    &                 ROAIRO, ZIDLOLD, CICOVER, CITHICK,
!                      IREAD, LWCUR)
!          *CDTWIS* - DATE OF FIRST WIND FIELD.
!          *CDTWIE* - DATE OF LAST FIRST WIND FIELD.
!          *IJS:IJL - ARRAYS DIMENSION
!          *U10OLD* - WIND SPEED.
!          *THWOLD* - WIND DIRECTION (RADIANS).
!          *USOLD*  - FRICTION VELOCITY.
!          *Z0OLD*  - ROUGHNESS LENGTH IN M.
!          *ROAIRO* - AIR DENSITY IN KG/M3.
!          *ZIDLOLD*- Zi/L 
!                     (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!          *CICOVER*- SEA ICE COVER.
!          *CITHICK*- SEA ICE THICKNESS. 
!          *IREAD*  - PROCESSOR WHICH WILL ACCESS THE FILE ON DISK
!          *LWCUR*  - LOGICAL INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS


!     METHOD.
!     -------

!       NO TIME INTERPOLATION:
!       WINDFIELDS ARE PROCESSED EVERY IDELWI SECONDS (U,V),
!       THE WINDS INTERPOLATED IN SPACE ONLY (US,DS)

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.
!       *AIRSEA*    - TOTAL STRESS IN SURFACE LAYER.
!       *GETWND*    - READ A WINDFIELD AND COMPUTE WIND
!                     FOR ALL BLOCKS (US,DS).
!       *INCDATE*   - INCREMENT DATE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWSTAT  , ONLY : IDELPRO  ,IDELWO   ,NPROMA_WAM
      USE YOWPHYS  , ONLY : XNLEV
      USE YOWTEST  , ONLY : IU06
      USE YOWWIND  , ONLY : CDA      ,CDTNEXT  ,NSTORE   ,FF_NEXT

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "cdustarz0.intfb.h"
#include "getwnd.intfb.h"
#include "incdate.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      
      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(INOUT) ::            &
     &               U10OLD, THWOLD, USOLD, Z0OLD,                    &
     &               ROAIRO, ZIDLOLD, CICOVER, CITHICK

      CHARACTER(LEN=14), INTENT(IN) :: CDTWIS, CDTWIE

      LOGICAL, INTENT(IN) :: LWCUR


      INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA
      INTEGER(KIND=JWIM) :: IJ
      INTEGER(KIND=JWIM) :: ICODE_WND

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB),DIMENSION(IJS:IJL) :: U10, US, THW, ADS, ZIDL, CICR, CITH, CD

      CHARACTER(LEN=14) :: CDTWIH, CDT, ZERO

      LOGICAL :: LWNDFILE, LCLOSEWND

! ----------------------------------------------------------------------

!*    1. INITIALISE WIND REQUEST DATE AND PROCESS FIRST WINDFIELD
!*       IF COLD START.
!        --------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('NOTIM',0,ZHOOK_HANDLE)

      ZERO = ' '
      CDTWIH = CDTWIS

      NPROMA=NPROMA_WAM

      LCLOSEWND=.FALSE.
      IF (LWCOU) THEN
        LWNDFILE=.FALSE.
      ELSE
        LWNDFILE=.TRUE.
      ENDIF

      IF (CDA == ZERO) THEN
        CDA = CDTWIS
        CALL GETWND (IJS, IJL,                          &
     &               U10OLD, USOLD,                     &
     &               THWOLD,                            &
     &               ROAIRO, ZIDLOLD,                   &
     &               CICOVER, CITHICK,                  &
     &               CDA, LWNDFILE, LCLOSEWND, IREAD,   &
     &               LWCUR, ICODE_WND)


          CALL GSTATS(1493,0)
!$OMP     PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
          DO JKGLO=IJS,IJL,NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,IJL)
            CALL CDUSTARZ0 (KIJS, KIJL, U10OLD(KIJS), XNLEV,            &
     &                      CD(KIJS), USOLD(KIJS), Z0OLD(KIJS))
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1493,1)

        IF (CDA == CDTWIE) THEN
          IF (LHOOK) CALL DR_HOOK('NOTIM',1,ZHOOK_HANDLE)
          RETURN
        ENDIF
        CALL INCDATE (CDTWIH, IDELWO)
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
        NSTORE=1

        IF (.NOT.ALLOCATED(CDTNEXT)) ALLOCATE(CDTNEXT(NSTORE))

        IF (.NOT.ALLOCATED(FF_NEXT)) ALLOCATE(FF_NEXT(IJS:IJL,NSTORE))

        CDTNEXT(1)=CDTWIH

        CALL GETWND (IJS, IJL,                              &
     &               U10, US,                               &
     &               THW,                                   &
     &               ADS, ZIDL,                             &
     &               CICR, CITH,                            &
     &               CDTWIH, LWNDFILE, LCLOSEWND, IREAD,    &
     &               LWCUR, ICODE_WND)

        IF (ICODE_WND == 3 ) THEN
          FF_NEXT(IJS:IJL,1)%WSWAVE = U10(IJS:IJL)
        ELSE
          FF_NEXT(IJS:IJL,1)%USTAR  = US(IJS:IJL)
        ENDIF
        FF_NEXT(IJS:IJL,1)%WDWAVE = THW(IJS:IJL)
        FF_NEXT(IJS:IJL,1)%AIRD   = ADS(IJS:IJL)
        FF_NEXT(IJS:IJL,1)%ZIDL   = ZIDL(IJS:IJL)
        FF_NEXT(IJS:IJL,1)%CIFR   = CICR(IJS:IJL)
        FF_NEXT(IJS:IJL,1)%CITH   = CITH(IJS:IJL)


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

        IF (LHOOK) CALL DR_HOOK('NOTIM',1,ZHOOK_HANDLE)
        RETURN

      ENDIF

      IF (LHOOK) CALL DR_HOOK('NOTIM',1,ZHOOK_HANDLE)

      END SUBROUTINE NOTIM
