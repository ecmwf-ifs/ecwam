      SUBROUTINE READSAT (IU06, IUME, CDTPRO,
     1                    CDATE, RLAT, RLON, SWH, WS, flagged, EOFD)
C
C--------------------------------------------------------------------
C
C**** *READSAT* - READ RADAR ALTIMETER DATA. (ERS_1)
C
C     H. GUNTHER     GKSS/ECMWF       AUGUST 1991
C
C     PURPOSE.
C     --------
C
C       ACCESS ERS1 FLAGGED DATA FILES IN ECFILE AND
C       READ ONE MEASUREMENT.
C
C**   INTERFACE.
C     ----------
C
C      *CALL* *READSAT (IU06, CDTPRO, CDATE, RLAT, RLON, SWH, WS, EOFD)*
C
C          *IU06*     - PRINTER UNIT.
C          *IUME*     - INPUT UNIT FOR MEASUREMENTS.
C          *CDTPRO*   - MODEL TIME (YYMMDDHHMM).
C          *CDATE*    - DATE OF MEASUREMENT (YYMMDDHHMM).
C          *RLAT*     - LATITUDE OF MEASUREMENT (DEGREE).
C          *RLON*     - LONGITUDE OF MEASUREMENT (DEGREES).
C          *SWH*      - WAVE HEIGHT.
C          *WS*       - WIND SPEED.
C          *FLAGGED*  - DATA QC-FLAG (.F.: PASSED QC, .T.: REJECTED BY QC)
C          *EOFD*     - END OF DATA. (.TRUE. , .FALSE.)
C
C     METHOD.
C     -------
C
C         NONE.
C
C      EXTERNALS.
C      ----------
C
C         NONE.
C
C      MODIFICATIONS.
C      --------------
C
C         FEB 92 Bjoern Hansen:
C                Call to ecfile dismissed. Files are assumed to be
C                online. File name convention changed to 3 character
C                file identifier and 10 digit DATE/TIME group rather
C                than 4 character file identifier and 9 digit ops like
C                DATE/TYPE group.
C
C-----------------------------------------------------------------------

      PARAMETER ( jpv=42, jpf=6)
C
      REAL pvalu(jpv)
      INTEGER iflag(jpf)
C
      CHARACTER* 3 yofid
      CHARACTER*10 CDATE, CDTPRO, CDTGEO, CDTGE1, ZERO, CTIME
      CHARACTER*14 C14
      PARAMETER(ZERO=' ')
      LOGICAL EOFD, FLAGGED

      DATA yofid /'RFL'/
      DATA CDTGEO, CDTGE1 /ZERO, ZERO/
C
C ----------------------------------------------------------------------
C
C*    1. INITIALIZE.
C        -----------
C
 1000 CONTINUE
      IF (CDTGE1.EQ.ZERO) THEN
         CDTGE1 = CDTPRO
         EOFD = .TRUE.
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    2. GENERATE FIRST FILE DATES AND FETCH FILE FROM ECFILE.
C        -----------------------------------------------------
C
 2000 CONTINUE
      IF (EOFD) THEN
         CALL confile (iu06, iume, cdtge1, yofid, ifail)
         IF(IFAIL.NE.0) THEN
            WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++'
            WRITE(IU06,*) ' +                                       +'
            WRITE(IU06,*) ' +      WARNING ERROR IN --READSAT--     +'
            WRITE(IU06,*) ' +      ==========================       +'
            WRITE(IU06,*) ' +                                       +'
            WRITE(IU06,*) ' +      NO SATELLITE MEASUREMENTS        +'
            WRITE(IU06,*) ' + MODEL TIME IS  CDTPRO = ', CDTPRO
            WRITE(IU06,*) ' + FILE DATES ARE CDTGEO = ', CDTGEO
            WRITE(IU06,*) ' +                CDTGE1 = ', CDTGE1
            WRITE(IU06,*) ' + IFAIL = ', IFAIL
            WRITE(IU06,*) ' + EOFD  = ', EOFD
            WRITE(IU06,*) ' + PROCESSING WILL BE CONTINUED          +'
            WRITE(IU06,*) ' +                                       +'
            WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++'
            CDTGE1 = ZERO
            RETURN
         ENDIF
         EOFD = .FALSE.
         READ (iume,END=4000) ctime
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    3. READ THE MEASUREMENT.
C        ---------------------
C
 3000 CONTINUE
cccc      READ(IUME,END=4000) C14, (pvalu(j),j=1,jpv), (iflag(j),j=1,jpf)
      READ(IUME,END=4000) C14, pvalu, iflag
      rlat = pvalu(10)
      rlon = pvalu(11)
      ws   = pvalu(13)
      swh  = pvalu(15)
      ifl=0
      DO 3100 j=1,6
        ifl=ifl+iflag(j)
 3100 CONTINUE
      IF ( ifl .EQ. 0 ) flagged = .FALSE.
      CDATE = C14(3:12)
      IF (CDATE.LT.'9111110901') SWH = SWH*1.25

      RETURN
C
C ----------------------------------------------------------------------
C
C*    4. END OF FILE, FETCH NEXT FILE AND TRY AGAIN.
C        -------------------------------------------
C
 4000 CONTINUE
      EOFD = .TRUE.
      CDTGE1 = ZERO
C     CALL INCDATE (CDTGE1, 6*3600)
      RETURN
      END
