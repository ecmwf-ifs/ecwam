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

!         *CDATE*    CHAR*12    DATE (YYMMDDHHmm) OF FILE TO BE SAVED.
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
!        BY YYYYMMDDHHmm_ddddhhmi WHERE YYMMDDHHmm IS THE LAST ANALYSIS
!        TIME AND ddddhhmi IS THE FORCECAST RANGE 
!        IN dddd DAYS, hh HOURS AND mi MINUTES.
!        ddddhhmi=00000000  FOR ANALYSIS FILES.
!        THESE NUMBERS ARE CONSTRUCTED FROM CDATE AND CDATEF.

!        THEN IF KCPLEN .NOT.EQ.0
!        THE PATH NAME IS ADDED. THE PATH IS DETERMINED CPAD
!        IN THE INPUT NAMELIST 

!     REFERENCES.
!     -----------

!       NONE.

! ----------------------------------------------------------------------

      CHARACTER FILEID*3, FILENA*24
      CHARACTER CDATE*12, CDATEF*12, CDATEH*12
      CHARACTER*(*) CPAD, FILENAME

!*    1.0  CONSTRUCT FILE NAME ALONE 
!          -------------------------

      FILENA = FILEID
      IF (CDATE.LE.CDATEF) THEN
        CDATEH = CDATE
        ISHIFT = 0
      ELSE
        CDATEH = CDATEF
        CALL DIFDATE (CDATEF, CDATE, ISHIFT)
      ENDIF
      FILENA(4:15) = CDATEH
      FILENA(16:16) = "_"

      IDDDD=ISHIFT/86400
      IHH=(ISHIFT-IDDDD*86400)/3600
      IMIMI=(ISHIFT-IDDDD*86400-IHH*3600)/60

      WRITE (FILENA(17:24),'(I4.4,I2.2,I2.2)') IDDDD,IHH,IMIMI
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
