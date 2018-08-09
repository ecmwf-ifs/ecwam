      SUBROUTINE GETWND (MIJS, MIJL,                                    &
     &                   U10, US,                                       &
     &                   THW,                                           &
     &                   ADS, ZIDL,                                     &
     &                   CICVR, CITH,                                   &
     &                   CDTWIS, LWNDFILE, LCLOSEWND, IREAD,            &
     &                   LWCUR, ICODE_WND)

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

!       *CALL* *GETWND (MIJS, MIJL,
!                       U10, THW, ADS, ZIDL, CICVR, CITH,
!                       CDTWIS, LWNDFILE, LCLOSEWND,
!                       LWCUR, ICODE_WND)*
!         *MIJS*    - INDEX OF FIRST GRIDPOINT
!         *MIJL*    - INDEX OF LAST GRIDPOINT
!         *U10*    - MAGNITUDE OF 10m WIND AT EACH POINT AND BLOCK.
!         *THW*    - DIRECTION OF 10m WIND AT EACH POINT AND BLOCK.
!         *ADS*    - AIR DENSITY AT EACH POINT AND BLOCK.
!         *ZIDL*   - Zi/L  AT EACH POINT AND BLOCK.
!         *CICVR*  - SEA ICE COVER.
!         *CITH*   - SEA ICE THICKNESS.
!         *CDTWIS* - DATE OF WIND FIELD TO BE LOOKED FOR.
!         *LWNDFILE - FLAG USED TO DETERMINE WHETHER WINDS ARE READ FROM
!                     FILE OR ARE AVAILABLE IN ARRAY FIELDG (SEE *IFSTOWAM).
!         *LCLOSEWND* IF TRUE THE INPUT FILE WILL BE CLOSED AND
!                     THE UNIT RESET
!         *IREAD*  - PROCESSOR WHICH WILL ACCESS THE FILE ON DISK
!                    (IF NEEDED)
!         *LWCUR*  -  LOGICAL INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS
!         *ICODE_WND* SPECIFIES WHICH OF U10 OR US HAS BEEN UPDATED:
!                     U10: ICODE_WND=3
!                     US:  ICODE_WND=1 OR 2

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.
!       *READWIND*   - READING WINDS.
!       *WAMWND*    - CALCULATE WIND IN WAM POINTS.

!     REFERENCE.
!     ----------

!       NONE.


!    MODIFIED BY:
!    ------------
!    B. HANSEN    ECMWF 1997
!                 RESTRUCTURE CALL TO READWIND.
!
!    S. ABDALLA   ECMWF OCTOBER 2001
!                 MODIFICATION THE CALL TO READWIND WAMWND; AND 
!                 INCLUSION OF AIR DENSITY AND Zi/L.
!    J. BIDLOT    ECMWF NOVEMEBR 2003
!                 INTRODUCE openMP
!    J. BIDLOT    ECMWF AUGUST 2006
!                 SEA ICE FRACTION IF COUPLED.
!    J. BIDLOT    ECMWF AUGUST 2008
!                 READWIND WAS SPLIT BETWEEN READWIND AND IFSTOWAM:
!                 IN COUPLED RUNS, READWIND IS NOT CALLED 
! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWICE   , ONLY : IPARAMCI
      USE YOWMPP   , ONLY : IRANK
      USE YOWPARAM , ONLY : CLDOMAIN , LWDINTS
      USE YOWSTAT  , ONLY : NPROMA_WAM 
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWIND  , ONLY : WSPMIN   ,IUNITW
      USE YOWWNDG  , ONLY : ICODE    ,ICODE_CPL
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWUNPOOL, ONLY : LLUNSTR
      USE UNWAM, ONLY : USE_DIRECT_WIND_FILE
      USE UNSTRUCT_WIND, ONLY : SET_WIND_UNSTRUCTURED
      USE GRIB_API_INTERFACE

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "micep.intfb.h"
#include "readwind.intfb.h"
#include "wamwnd.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      INTEGER(KIND=JWIM), INTENT(IN) :: MIJS, MIJL
      INTEGER(KIND=JWIM), INTENT(OUT) :: ICODE_WND

      REAL(KIND=JWRB), DIMENSION (MIJS:MIJL), INTENT(INOUT) :: U10, US 
      REAL(KIND=JWRB), DIMENSION (MIJS:MIJL), INTENT(OUT) :: THW
      REAL(KIND=JWRB), DIMENSION (MIJS:MIJL), INTENT(OUT) :: ADS, ZIDL
      REAL(KIND=JWRB), DIMENSION (MIJS:MIJL), INTENT(OUT) :: CICVR, CITH

      CHARACTER(LEN=14), INTENT(IN) :: CDTWIS

      LOGICAL, INTENT(IN) :: LWNDFILE, LCLOSEWND, LWCUR


      INTEGER(KIND=JWIM) :: JKGLO,KIJS,KIJL,NPROMA
      INTEGER(KIND=JWIM) :: I, J

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      LOGICAL :: LLNOTOPENED
      LOGICAL, SAVE :: LONLYONCE
      LOGICAL :: IsAssigned

      CHARACTER(LEN=14) :: CDTWIR
      CHARACTER(LEN=24) :: FILNM

      DATA LONLYONCE /.TRUE./

      SAVE CDTWIR
      SAVE LLNOTOPENED

! ----------------------------------------------------------------------

!*    1. WIND DATA ARE READ
!        ------------------

      IF (LHOOK) CALL DR_HOOK('GETWND',0,ZHOOK_HANDLE)

 1000 CONTINUE

      IF(IUNITW.EQ.0) THEN
        LLNOTOPENED=.TRUE.
      ELSE
        LLNOTOPENED=.FALSE.
      ENDIF

      NPROMA=NPROMA_WAM

!     GET FORCING FIELDS FROM INPUT FILES (if needed)
!     -----------------------------------
      IsAssigned=.FALSE.
      IF (LLUNSTR .and. USE_DIRECT_WIND_FILE) THEN
        IsAssigned=.TRUE.
        CALL SET_WIND_UNSTRUCTURED
      END IF


      IF(LWNDFILE .AND. (IsAssigned .EQV. .FALSE.)) THEN
        CALL READWIND (CDTWIR, FILNM, LLNOTOPENED, IREAD)

        ICODE_WND = ICODE

!       CHECK WIND FIELD DATE

        IF (CDTWIR.LT.CDTWIS) THEN
!         DATE OF INPUT FIELD IS BEFORE REQUESTED DATE
!         TRY AGAIN
          IF(LWNDFILE) THEN
            IF (ITEST.GT.1) THEN
              WRITE(IU06,*) ' SUB. GETWND - BEFORE REQUESTED DATE '
              WRITE(IU06,*) ' CDTWIR= ',CDTWIR
              WRITE(IU06,*) ' CDTWIS= ',CDTWIS
              WRITE(IU06,*) ' SUB. GETWND - CALLING READWIND AGAIN'
              CALL FLUSH(IU06)
            ENDIF
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
        ELSEIF (CDTWIR.GT.CDTWIS) THEN

!         DATE OF INPUT FIELD IS LATER THAN REQUESTED DATE
          WRITE (IU06,*) ' ****************************************'
          WRITE (IU06,*) ' *                                      *'
          WRITE (IU06,*) ' *      FATAL ERROR SUB. GETWND         *'
          WRITE (IU06,*) ' *      =======================         *'
          WRITE (IU06,*) ' * WIND DATE IS LATER THAN EXPECTED     *'
          IF(LWNDFILE) THEN
            WRITE (IU06,*) ' * DATE READ IS    CDTWIR = ', CDTWIR
          ELSE
            WRITE (IU06,*) ' * DECODED DATE IS CDTWIR = ', CDTWIR
          ENDIF
          WRITE (IU06,*) ' * DATE EXPECTED IS CDTWIS = ', CDTWIS
          WRITE (IU06,*) ' *                                      *'
          WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS     *'
          WRITE (IU06,*) ' *                                      *'
          WRITE (IU06,*) ' ****************************************'
          CALL ABORT1
        ENDIF

        IF (LCLOSEWND .AND. LWNDFILE .AND.                              &
     &     .NOT.(CLDOMAIN.EQ.'s' .OR. LWDINTS) ) THEN
          IF(IRANK.EQ.IREAD) THEN
            CALL IGRIB_CLOSE_FILE(IUNITW)
            LLNOTOPENED = .TRUE.
            IUNITW=0
            IF (ITEST.GT.1) WRITE(IU06,*) ' SUB. GETWND - CLOSE ', FILNM
          ENDIF
        ENDIF

      ELSE
        ICODE_WND = ICODE_CPL
      ENDIF
! ----------------------------------------------------------------------

!*    3. INTERPOLATE AND BLOCK WINDFIELD
!        -------------------------------

! Mod for OPENMP
        CALL GSTATS(1444,0)
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO=MIJS,MIJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,MIJL)
          CALL WAMWND (KIJS, KIJL,                                      &
     &                 U10(KIJS), US(KIJS),                             &
     &                 THW(KIJS), ADS(KIJS), ZIDL(KIJS), CITH(KIJS),    &
     &                 LWCUR, ICODE_WND)
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1444,1)

        IF(LONLYONCE) THEN
          WRITE (IU06,*) ' '
          WRITE (IU06,*) ' SUB. GETWND : '
          WRITE (IU06,*) ' '
          WRITE (IU06,*) ' WIND SPEEDS LOWER THAN ',WSPMIN, ' M/S'
          WRITE (IU06,*) ' WERE RESET TO  ',WSPMIN, ' M/S'
          WRITE (IU06,*) ' '
          CALL FLUSH(IU06)
          LONLYONCE=.FALSE.
        ENDIF

!       USE THE SEA ICE FRACTION TO DEFINE THE SEA ICE BOUNDARY

        CALL GSTATS(1444,0)
! Mod for OPENMP
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO=MIJS,MIJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,MIJL)
          CALL MICEP(IPARAMCI, CICVR(KIJS), CITH(KIJS), KIJS, KIJL)
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1444,1)

      IF (LHOOK) CALL DR_HOOK('GETWND',1,ZHOOK_HANDLE)

      END SUBROUTINE GETWND
