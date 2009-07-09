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

!         *CDATE*    CHAR*14    DATE (YYYYMMDDHHmmss) OF FILE TO BE SAVED.
!         *CDATEF*   CHAR*14    DATE TIME GROUP OF FORECAST START.
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
!        BY YYYYMMDDHHmmss_ddddddhhmmss WHERE YYYYMMDDHHmmss IS THE LAST ANALYSIS
!        TIME AND ddddddhhmmss IS THE FORCECAST RANGE 
!        IN dddddd DAYS, hh HOURS, mi MINUTES AND ss seconds.
!        ddddddhhmmss=000000000000  FOR ANALYSIS FILES.
!        THESE NUMBERS ARE CONSTRUCTED FROM CDATE AND CDATEF.

!        THEN IF KCPLEN .NOT.EQ.0
!        THE PATH NAME IS ADDED. THE PATH IS DETERMINED CPAD
!        IN THE INPUT NAMELIST 

!     REFERENCES.
!     -----------

!       NONE.

! ----------------------------------------------------------------------
      INTEGER :: IMAXYEAR
      CHARACTER FILEID*3, FILENA*30
      CHARACTER*14 :: CDATE, CDATEF, CDATEH, CDT
      CHARACTER*(*) CPAD, FILENAME

!*    1.0  CONSTRUCT FILE NAME ALONE 
!          -------------------------

      IREP=KIND(IMAXYEAR)
      NDAY=366
      IYEARSEC=NDAY*86400
      IMAXYEAR=(2**(8*IREP-1)-1)/IYEARSEC

      FILENA = FILEID
      IF (CDATE.LE.CDATEF) THEN
        CDATEH = CDATE
        ISHIFT = 0
        ISHIFTDAY=0
      ELSE
        CDATEH = CDATEF

        READ(CDATEF,'(I4, 5I2)')IYEAR1,IMON1,IDAY1,IHOUR1,IMIN1,ISEC1
        READ(CDATE ,'(I4, 5I2)')IYEAR2,IMON2,IDAY2,IHOUR2,IMIN2,ISEC2
        ISHIFTDAY=0
        CDT=CDATE
!       IN CASE THE DATE DIFFERENCE IS LARGER THAN WHAT DIFDATE CAN DO.
        IYRDIF=IYEAR2-IYEAR1
        IF(ABS(IYRDIF).GE.IMAXYEAR) THEN
          NCHUNCK=ABS(IYRDIF)/IMAXYEAR
          DO IC=1,NCHUNCK
            ISFT=-SIGN(1,IYRDIF)*IMAXYEAR*IYEARSEC
            CALL INCDATE(CDT,ISFT)
            ISHIFTDAY=ISHIFTDAY-ISFT/86400
          ENDDO
        ENDIF

        CALL DIFDATE (CDATEF, CDT, ISHIFT)
      ENDIF
      FILENA(4:17) = CDATEH
      FILENA(18:18) = "_"

      IDDDDDD=ISHIFT/86400
      IHH=(ISHIFT-IDDDDDD*86400)/3600
      IMIMI=(ISHIFT-IDDDDDD*86400-IHH*3600)/60
      ISS=ISHIFT-IDDDDDD*86400-IHH*3600-IMIMI*60

      IDDDDDD=IDDDDDD+ISHIFTDAY

      WRITE (FILENA(19:30),'(I6.6,3I2.2)')IDDDDDD,IHH,IMIMI,ISS
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
