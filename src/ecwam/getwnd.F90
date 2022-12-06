SUBROUTINE GETWND (IFROMIJ, JFROMIJ,                      &
 &                 NXS, NXE, NYS, NYE, FIELDG,            &
 &                 UCUR, VCUR,                            &
 &                 U10, US,                               &
 &                 THW,                                   &
 &                 ADS, WSTAR,                            &
 &                 CICOVER, CITHICK,                      &
 &                 CDTWIS, LWNDFILE, LCLOSEWND, IREAD,    &
 &                 LWCUR, NEMOCICOVER, NEMOCITHICK,       &
 &                 ICODE_WND)

! ----------------------------------------------------------------------

!**** *GETWND* - ROUTINE TO READ AND PROCESS ONE WINDFIELD.

!*    PURPOSE.
!     --------

!        READ A WINDFIELD FROM THE WINDFILE (SEARCH FOR IT)
!        AND CALCULATES THE WIND VELOCITY  AND DIRECTION
!        FOR ALL WAM BLOCKS.
!        ALSO INPUT THE OTHER FORCING FIELD.

!**   INTERFACE.
!     ----------

!       *CALL* *GETWND (IFROMIJ, JFROMIJ,
!                       NXS, NXE, NYS, NYE, FIELDG,
!                       UCUR, VCUR,
!                       U10, THW, ADS, WSTAR, CICOVER, CITHICK,
!                       CDTWIS, LWNDFILE, LCLOSEWND,
!                       LWCUR, NEMOCICOVER, NEMOCITHICK, 
!                       ICODE_WND)*
!         *IFROMIJ*  POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!         *JFROMIJ*  POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!         *NXS:NXE*  FIRST DIMENSION OF FIELDG
!         *NYS:NYE*  SECOND DIMENSION OF FIELDG
!         *FIELDG* - INPUT FORCING FIELDS ON THE WAVE MODEL GRID
!         *UCUR*   - U-COMPONENT OF THE SURFACE CURRENT
!         *VCUR*   - V-COMPONENT OF THE SURFACE CURRENT
!         *U10*    - MAGNITUDE OF 10m WIND AT EACH POINT AND BLOCK.
!         *THW*    - DIRECTION OF 10m WIND AT EACH POINT AND BLOCK.
!         *ADS*    - AIR DENSITY AT EACH POINT AND BLOCK.
!         *WSTAR*  - CONVECTIVE VELOCITY.
!         *CICOVER*  - SEA ICE COVER.
!         *CITHICK*   - SEA ICE THICKNESS.
!         *CDTWIS* - DATE OF WIND FIELD TO BE LOOKED FOR.
!         *LWNDFILE - FLAG USED TO DETERMINE WHETHER WINDS ARE READ FROM
!                     FILE OR ARE AVAILABLE IN ARRAY FIELDG (SEE *IFSTOWAM).
!         *LCLOSEWND* IF TRUE THE INPUT FILE WILL BE CLOSED AND
!                     THE UNIT RESET
!         *IREAD*  - PROCESSOR WHICH WILL ACCESS THE FILE ON DISK
!                    (IF NEEDED)
!         *LWCUR*  -  LOGICAL INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS
!         *NEMOCICOVER  NEMO SEA ICE COVER (if used)
!         *NEMOCITHICK  NEMO SEA ICE THICKNESS (if used)
!         *ICODE_WND* SPECIFIES WHICH OF U10 OR US HAS BEEN UPDATED:
!                     U10: ICODE_WND=3
!                     US:  ICODE_WND=1 OR 2

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU, JWRO
      USE YOWDRVTYPE  , ONLY : FORCING_FIELDS

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWICE   , ONLY : IPARAMCI
      USE YOWMPP   , ONLY : IRANK
      USE YOWPARAM , ONLY : CLDOMAIN , LWDINTS
      USE YOWTEST  , ONLY : IU06
      USE YOWWIND  , ONLY : IUNITW
      USE YOWWNDG  , ONLY : ICODE    ,ICODE_CPL

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE YOWGRIB  , ONLY : IGRIB_CLOSE_FILE

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "micep.intfb.h"
#include "readwind.intfb.h"
#include "wamwnd.intfb.h"

      INTEGER(KIND=JWIM), DIMENSION(NPROMA_WAM, NCHNK), INTENT(IN) :: IFROMIJ  ,JFROMIJ
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), DIMENSION(NXS:NXE, NYS:NYE), INTENT(INOUT) :: FIELDG
      REAL(KIND=JWRB), DIMENSION (NPROMA_WAM, NCHNK), INTENT(IN) :: UCUR, VCUR
      REAL(KIND=JWRB), DIMENSION (NPROMA_WAM, NCHNK), INTENT(INOUT) :: U10, US 
      REAL(KIND=JWRB), DIMENSION (NPROMA_WAM, NCHNK), INTENT(OUT) :: THW
      REAL(KIND=JWRB), DIMENSION (NPROMA_WAM, NCHNK), INTENT(OUT) :: ADS, WSTAR
      REAL(KIND=JWRB), DIMENSION (NPROMA_WAM, NCHNK), INTENT(OUT) :: CICOVER, CITHICK
      CHARACTER(LEN=14), INTENT(IN) :: CDTWIS
      LOGICAL, INTENT(IN) :: LWNDFILE, LCLOSEWND
      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      LOGICAL, INTENT(IN) :: LWCUR
      REAL(KIND=JWRO), DIMENSION (NPROMA_WAM, NCHNK), INTENT(IN) :: NEMOCICOVER, NEMOCITHICK
      INTEGER(KIND=JWIM), INTENT(OUT) :: ICODE_WND


      INTEGER(KIND=JWIM) :: ICHNK, KIJS, KIJL 
      INTEGER(KIND=JWIM) :: I, J

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      LOGICAL :: LLNOTOPENED

      CHARACTER(LEN=14) :: CDTWIR
      CHARACTER(LEN=24) :: FILNM

      SAVE CDTWIR
      SAVE LLNOTOPENED

! ----------------------------------------------------------------------

!*    1. WIND DATA ARE READ
!        ------------------

IF (LHOOK) CALL DR_HOOK('GETWND',0,ZHOOK_HANDLE)

1000 CONTINUE

      IF (IUNITW == 0) THEN
        LLNOTOPENED=.TRUE.
      ELSE
        LLNOTOPENED=.FALSE.
      ENDIF

!     GET FORCING FIELDS FROM INPUT FILES (if needed)
!     -----------------------------------

      IF (LWNDFILE) THEN
        CALL READWIND (CDTWIR, FILNM, LLNOTOPENED, IREAD,  &
     &                 NXS, NXE, NYS, NYE, FIELDG)

        ICODE_WND = ICODE

!       CHECK WIND FIELD DATE

        IF (CDTWIR < CDTWIS) THEN
!         DATE OF INPUT FIELD IS BEFORE REQUESTED DATE
!         TRY AGAIN
          IF (LWNDFILE) THEN
            GOTO 1000
          ELSE
            WRITE (IU06,*) ' ****************************************'
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' *      FATAL ERROR SUB. GETWND          *'
            WRITE (IU06,*) ' *      =======================          *'
            WRITE (IU06,*) ' * WIND DATE IS EARLIER THAN EXPECTED    *'
            WRITE (IU06,*) ' * DECODED DATE IS  CDTWIR = ', CDTWIR
            WRITE (IU06,*) ' * DATE EXPECTED IS CDTWIS = ', CDTWIS
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS      *'
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' ****************************************'
            CALL ABORT1
          ENDIF
        ELSEIF (CDTWIR > CDTWIS) THEN

!         DATE OF INPUT FIELD IS LATER THAN REQUESTED DATE
          WRITE (IU06,*) ' ****************************************'
          WRITE (IU06,*) ' *                                      *'
          WRITE (IU06,*) ' *      FATAL ERROR SUB. GETWND         *'
          WRITE (IU06,*) ' *      =======================         *'
          WRITE (IU06,*) ' * WIND DATE IS LATER THAN EXPECTED     *'
          IF (LWNDFILE) THEN
            WRITE (IU06,*) ' * DATE READ IS     CDTWIR = ', CDTWIR
          ELSE
            WRITE (IU06,*) ' * DECODED DATE IS  CDTWIR = ', CDTWIR
          ENDIF
          WRITE (IU06,*) ' * DATE EXPECTED IS CDTWIS = ', CDTWIS
          WRITE (IU06,*) ' *                                      *'
          WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS     *'
          WRITE (IU06,*) ' *                                      *'
          WRITE (IU06,*) ' ****************************************'
          CALL ABORT1
        ENDIF

        IF (LCLOSEWND .AND. LWNDFILE .AND.                              &
     &     .NOT.(CLDOMAIN == 's' .OR. LWDINTS) ) THEN
          IF (IRANK == IREAD) THEN
            CALL IGRIB_CLOSE_FILE(IUNITW)
            LLNOTOPENED = .TRUE.
            IUNITW=0
          ENDIF
        ENDIF

      ELSE
        ICODE_WND = ICODE_CPL
      ENDIF
! ----------------------------------------------------------------------

!*    3. INTERPOLATE AND BLOCK WINDFIELD
!        -------------------------------

        CALL GSTATS(1444,0)
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ICHNK, KIJS, KIJL)
        DO ICHNK = 1, NCHNK
          KIJS=1
          KIJL=NPROMA_WAM
          CALL WAMWND (KIJS, KIJL,                          &
     &                 IFROMIJ(:,ICHNK), JFROMIJ(:,ICHNK),  &
     &                 NXS, NXE, NYS, NYE, FIELDG,          &
     &                 UCUR(:,ICHNK), VCUR(:,ICHNK),        &
     &                 U10(:,ICHNK), US(:,ICHNK),           &
     &                 THW(:,ICHNK), ADS(:,ICHNK),          &
     &                 WSTAR(:,ICHNK), CITHICK(:,ICHNK),    &
     &                 LWCUR, ICODE_WND)


          CALL MICEP(IPARAMCI, KIJS, KIJL, IFROMIJ(:,ICHNK), JFROMIJ(:,ICHNK),  &
     &               NXS, NXE, NYS, NYE, FIELDG,                                &
     &               CICOVER(:,ICHNK), CITHICK(:,ICHNK),                        &
     &               NEMOCICOVER(:,ICHNK), NEMOCITHICK(:,ICHNK))
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1444,1)

IF (LHOOK) CALL DR_HOOK('GETWND',1,ZHOOK_HANDLE)

END SUBROUTINE GETWND