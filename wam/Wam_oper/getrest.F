      SUBROUTINE GETREST (IU06, IU12, IU14, IU15, IGL,
     &                    CDATE, CDATEF)
! ----------------------------------------------------------------------

!**** *GETREST*  GETS RESTART FILES

!     HEINZ GUNTHER   GKSS/ECMWF     MAY 1989

!*    PURPOSE.
!     --------

!       FETCHES RESTART FILES.

!**   INTERFACE.
!     ----------

!       *CALL* *GETREST (IU06, IU12, IU14, IU15, IGL,
!                        CDATE, CDATEF)*

!         *IU06*    INTEGER    OUTPUT UNIT FOR PRINTER MESSAGES.
!         *IU12*    INTEGER    OUTPUT UNIT BLOCKS OF SPECTRA.
!         *IU14*    INTEGER    OUTPUT UNIT SECOND LAT OF BLOCKS.
!         *IU15*    INTEGER    OUTPUT UNIT LAST WINDFIELDS.
!         *IGL*     INTEGER    NUMBER OF BLOCKS.
!         *CDATE*   CHAR*14    DATE/TIME OF LAST FIELD IN FILES.
!         *CDATEF*  CHAR*14    START DATE OF FORECAST PERIOD.

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

      CHARACTER*14 CDATE,CDATEF

! ----------------------------------------------------------------------

!*    1. ANALYSIS OR FORECATS FILES.
!        ---------------------------


      ICH = 0

      IF (CDATE.LE.CDATEF) THEN

!*    1.1 ANALYSIS FILES.
!         ---------------

!*    1.1.1 SPECTRA FILE.
!           -------------

        CALL GSFILE (IU06, IU12, ICH, CDATE, CDATEF, 'BLS', 'G')

!*    1.1.2 SECOND LATITUDE FILE IF MULTIBLOCK VERSION.
!           -------------------------------------------

        IF (IGL.GT.1) THEN
          CALL GSFILE (IU06, IU14, ICH, CDATE, CDATEF, 'SLA', 'G')
        ENDIF

!*    1.1.3 LAST WIND FIELDS.
!           -----------------

        CALL GSFILE (IU06, IU15, ICH, CDATE, CDATEF, 'LAW', 'G')
      ELSE

!*    1.2 FORECAST FILES.
!         ---------------

!*    1.2.1 SPECTRA FILE.
!           -------------

        CALL GSFILE (IU06, IU12, ICH, CDATE, CDATEF, 'BLS', 'G')

!*    1.2.2 SECOND LATITUDE FILE IF MULTIBLOCK VERSION.
!           -------------------------------------------

        IF (IGL.GT.1) THEN
          CALL GSFILE (IU06, IU14, ICH, CDATE, CDATEF, 'SLA', 'G')
        ENDIF

!*    1.2.3 LAST WIND FIELDS.
!           -----------------

        CALL GSFILE (IU06, IU15, ICH, CDATE, CDATEF, 'LAW', 'G')
      ENDIF

      RETURN
      END SUBROUTINE GETREST


