! ----------------------------------------------------------------------

      SUBROUTINE GRSTNAME (CDATE, CDATEF, FILEID, KCPLEN, CPAD,
     &                     FILENAME)

! ----------------------------------------------------------------------

!**** *GRSTNAME*  DETERMINES FILENAME (INCLUDING PATH) FOR RESTART FILES

!     J. BIDLOT     ECMWF         NOVEMBER 1996      PARALLEL RUN 
!     B. HANSEN     ECMWF         APRIL    1996      MOVE CYMPPATH TO
!                                                    COMTEXT.

!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!       *CALL* *GRSTNAME (CDATE, CDATEF, FILEID, KCPLEN, CPAD,
!                         FILENAME)*

!         *CDATE*    CHAR*12    DATE (YYMMDDHHMM) OF FILE TO BE SAVED.
!         *CDATEF*   CHAR*12    DATE TIME GROUP OF FORECAST START.
!         *FILEID*   CHARACTER  FILE ID.
!         *KCPLEN*   INTEGER    LENGTH OF CPAD
!         *CPAD*     CHARACTER  PATH FOR OUTPUT TO DISK  
!         *FILENAME* CHARACTER  FILE NAME INCLUDING PATH

!     EXTERNALS.
!     ----------

!       *DIFDATE*   - COMPUTE TIME DIFFERENCES.

!     METHOD.
!     -------

!        FILE NAME IS BUILT FROM THE 3 CHARACTER FILEID FOLLOWED
!        BY YYMMDDHHFFFF WHERE YYMMDDHH IS THE LAST ANALYSIS TIME AND
!        FFFF IS THE FORCECAST PERIOD IN HOURS.
!        THESE NUMBERS ARE CONSTRUCTED FROM CDATE AND CDATEF.

!        THEN IF KCPLEN .NOT.EQ.0
!        THE PATH NAME IS ADDED. THE PATH IS DETERMINED CPAD
!        IN THE INPUT NAMELIST 

!     REFERENCES.
!     -----------

!       NONE.

! ----------------------------------------------------------------------

      CHARACTER FILEID*3, FILENA*17
      CHARACTER CDATE*12, CDATEF*12, CDATEH*10
      CHARACTER*(*) CPAD, FILENAME

!*    1.0  CONSTRUCT FILE NAME ALONE 
!          -------------------------

      FILENA = FILEID
      IF (CDATE.LE.CDATEF) THEN
        CDATEH = CDATE(1:10)
        ISHIFT = 0
      ELSE
        CDATEH = CDATEF(1:10)
        CALL DIFDATE (CDATEF, CDATE, ISHIFT)
        ISHIFT = ISHIFT/3600
      ENDIF
      FILENA(4:13) = CDATEH
      WRITE (FILENA(14:17),'(I4.4)') ISHIFT
      LNAME = LEN_TRIM(FILENA)

!     2.0   GET PATH AND ADD IT TO FILE NAME
!           --------------------------------

      IF(KCPLEN.EQ.0) THEN
        FILENAME=FILENA(1:LNAME)
      ELSE
        FILENAME=CPAD(1:KCPLEN)//'/'//FILENA(1:LNAME)
      ENDIF

      RETURN
      END SUBROUTINE GRSTNAME
