      SUBROUTINE READSAT (IU06, IUME, CDTPRO,
     &                    CDATE, RLAT, RLON, SWH, WS, flagged, EOFD)

!--------------------------------------------------------------------

!**** *READSAT* - READ RADAR ALTIMETER DATA.

!     H. GUNTHER     GKSS/ECMWF       AUGUST 1991

!     PURPOSE.
!     --------

!       ACCESS ERS1 FLAGGED DATA FILES AND
!       READ ONE MEASUREMENT.

!**   INTERFACE.
!     ----------

!      *CALL* *READSAT (IU06, CDTPRO, CDATE, RLAT, RLON, SWH, WS, EOFD)*

!          *IU06*     - PRINTER UNIT.
!          *IUME*     - INPUT UNIT FOR MEASUREMENTS.
!          *CDTPRO*   - MODEL TIME (YYYYMMDDHHMM).
!          *CDATE*    - DATE OF MEASUREMENT (YYYYMMDDHHMMSS).
!          *RLAT*     - LATITUDE OF MEASUREMENT (DEGREE).
!          *RLON*     - LONGITUDE OF MEASUREMENT (DEGREES).
!          *SWH*      - WAVE HEIGHT.
!          *WS*       - WIND SPEED.
!          *FLAGGED*  - DATA QC-FLAG (.F.: PASSEDQC, .T.:REJECTED BY QC)
!          *EOFD*     - END OF DATA. (.TRUE. , .FALSE.)

!     METHOD.
!     -------

!         NONE.

!      EXTERNALS.
!      ----------

!         NONE.

!      MODIFICATIONS.
!      --------------

!         FEB 92 Bjoern Hansen:
!                Call to ecfile dismissed. Files are assumed to be
!                online. File name convention changed to 3 character
!                file identifier and 10 digit DATE/TIME group rather
!                than 4 character file identifier and 9 digit ops like
!                DATE/TYPE group.

!        JULY 1998 : changed to 12 digit DATE/TIME group

!-----------------------------------------------------------------------

      PARAMETER ( jpv=42, jpf=6)

      REAL pvalu(jpv)
      INTEGER iflag(jpf)

      CHARACTER* 3 yofid
      CHARACTER*12 CDTPRO, CDTGEO, CDTGE1, ZERO, CTIME
      CHARACTER*14 CDATE
      PARAMETER(ZERO=' ')
      LOGICAL EOFD, FLAGGED, LLTEST

      DATA yofid /'RFL'/
      DATA CDTGEO, CDTGE1 /ZERO, ZERO/

! ----------------------------------------------------------------------

!*    1. INITIALIZE.
!        -----------

      IF (CDTGE1.EQ.ZERO) THEN
        CDTGE1 = CDTPRO
        EOFD = .TRUE.
      ENDIF

! ----------------------------------------------------------------------

!*    2. GENERATE FIRST FILE DATES AND FETCH FILE FROM ECFILE.
!        -----------------------------------------------------

      IF (EOFD) THEN
        LLTEST=.TRUE.
        CALL CONFILE (IU06, IUME, CDTGE1, YOFID, IFAIL, LLTEST)
        IF(IFAIL.NE.0) THEN
          WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++'
          WRITE(IU06,*) ' +                                       +'
          WRITE(IU06,*) ' +      WARNING ERROR IN --READSAT--     +'
          WRITE(IU06,*) ' +      ==========================       +'
          WRITE(IU06,*) ' +                                       +'
          WRITE(IU06,*) ' +      NO ALTIMETER MEASUREMENTS        +'
          WRITE(IU06,*) ' + MODEL TIME IS  CDTPRO = ', CDTPRO
          WRITE(IU06,*) ' + FILE DATES ARE CDTGEO = ', CDTGEO
          WRITE(IU06,*) ' +                CDTGE1 = ', CDTGE1
          WRITE(IU06,*) ' + IFAIL = ', IFAIL
          WRITE(IU06,*) ' + EOFD  = ', EOFD
          WRITE(IU06,*) ' + PROCESSING WILL BE CONTINUED          +'
          WRITE(IU06,*) ' +                                       +'
          WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++'
          CALL FLUSH(IU06)
          CDTGE1 = ZERO
          RETURN
        ENDIF
        EOFD = .FALSE.
        READ (iume,END=4000) ctime
      ENDIF

! ----------------------------------------------------------------------

!*    3. READ THE MEASUREMENT.
!        ---------------------

      READ(IUME,END=4000) CDATE, pvalu, iflag
      rlat = pvalu(10)
      rlon = pvalu(11)
      ws   = pvalu(13)
      swh  = pvalu(15)
      ifl=0
      DO j=1,6
        ifl=ifl+iflag(j)
      ENDDO
      IF ( ifl .EQ. 0 ) flagged = .FALSE.
      IF (CDATE.LT.'19911111090100') SWH = SWH*1.25

      RETURN

! ----------------------------------------------------------------------

!*    4. END OF FILE, FETCH NEXT FILE AND TRY AGAIN.
!        -------------------------------------------

 4000 CONTINUE
      EOFD = .TRUE.
      CDTGE1 = ZERO
!     CALL INCDATE (CDTGE1, 6*3600)
      RETURN
      END SUBROUTINE READSAT
