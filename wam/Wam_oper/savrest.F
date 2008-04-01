      SUBROUTINE SAVREST (IU06, IU04, IU12, IU14, IU15, IGL,
     &                    CDATE, CDATEF, CDATER )

! ----------------------------------------------------------------------

!**** *SAVREST* - SAVES RESTART FILES.

!     HEINZ GUNTHER      GKSS/ECMWF     MAY 1989
!     J. BIDLOT          ECMWF          NOVEMBER 1996 : ADDED  CDATER
!     B. HANSEN                         APRIL   1997 : RESTART FACILITY
!     J. BIDLOT          ECMWF          SEPTEMBER : USE OF PBIO 

!*    PURPOSE.
!     --------

!       COPIES RESTART FILES.

!**   INTERFACE.
!     ----------

!       *CALL* *SAVREST (IU06, IU04, IU12, IU14, IU15, IGL,
!                        CDATE, CDATEF,CDATER)*
!         *IU06*   INTEGER   OUTPUT UNIT FOR PRINTER MESSAGES.
!         *IU04*   INTEGER   OUTPUT UNIT WAMINFO FILE.
!         *IU12*   INTEGER   OUTPUT UNIT BLOCKS OF SPECTRA.
!         *IU14*   INTEGER   OUTPUT UNIT SECOND LAT OF BLOCKS.
!         *IU15*   INTEGER   OUTPUT UNIT LAST WINDFIELDS.
!         *IGL*    INTEGER   NUMBER OF BLOCKS.
!         *CDATE*  CHAR*14   DATE/TIME OF LAST FIELD IN FILES.
!         *CDATEF* CHAR*14   START DATE OF FORECAST PERIOD.
!         *CDATER* CHAR*14   OUTPUT DATE FOR BOTH RESTART FILES
!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       *GSFILE*    - GETS OR SAVES FILES.

!     REFERENCE.
!     ----------

!       NONE.

! ---------------------------------------------------------------------

      CHARACTER*14 CDATE,CDATEF,CDATER

! ---------------------------------------------------------------------

!*    1. ANALYSIS OR FORECATS FILES.
!        ---------------------------

      ICH = 0

      IF (CDATE.LE.CDATEF) THEN

!*    1.1 ANALYSIS FILES.
!         ---------------

!*    1.1.1 SPECTRA FILE.
!           -------------

        CALL GSFILE (IU06, IU12, ICH, CDATE, CDATEF, 'BLS', 'S')

!*    1.1.2 SECOND LATITUDE FILE IF MULTIBLOCK.
!           -----------------------------------

        IF (IGL.GT.1) THEN
          CALL GSFILE (IU06, IU14, ICH, CDATE, CDATEF, 'SLA', 'S')
        ENDIF

!*    1.1.3 LAST WIND FIELDS.
!           -----------------

        IF((CDATE .EQ. CDATER ).OR.(CDATER .EQ. '00000000000000')) THEN
          CALL GSFILE (IU06, IU15, ICH, CDATE, CDATEF, 'LAW', 'S')
        ENDIF

      ELSE

!*    1.1 FORECAST FILES.
!         ---------------

!*    1.2.1 SPECTRA FILE.
!           -------------

        CALL GSFILE (IU06, IU12, ICH, CDATE, CDATEF, 'BLS', 'S')

!*    1.2.2 SECOND LATITUDE FILE IF MULTIBLOCK.
!           -----------------------------------

        IF (IGL.GT.1) THEN
          CALL GSFILE (IU06, IU14, ICH, CDATE, CDATEF, 'SLA', 'S')
        ENDIF

!*    1.2.3 LAST WIND FIELDS.
!           -----------------

        IF((CDATE .EQ. CDATER ).OR.(CDATER .EQ. '00000000000000')) THEN
          CALL GSFILE (IU06, IU15, ICH, CDATE, CDATEF, 'LAW', 'S')
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE SAVREST


