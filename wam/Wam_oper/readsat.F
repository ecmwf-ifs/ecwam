      SUBROUTINE READSAT (IU06, IUME, CDTPRO, IDENTI, ISENSOR,
     &                    CDATE, RLAT, RLON, IDES_HS, SWH,
     &                    IDES_WS, WS, FLAGGED, EOFD)

!--------------------------------------------------------------------

!**** *READSAT* - READ RADAR ALTIMETER DATA.

!     H. GUNTHER     GKSS/ECMWF       AUGUST 1991

!     PURPOSE.
!     --------

!       ACCESS ERS1 FLAGGED DATA FILES AND
!       READ ONE MEASUREMENT.

!**   INTERFACE.
!     ----------

!          *IU06*     - PRINTER UNIT.
!          *IUME*     - INPUT UNIT FOR MEASUREMENTS.
!          *CDTPRO*   - MODEL TIME (YYYYMMDDHHMM).
!          *IDENTI*   - SATELLITE IDENTIFIER.
!          *ISENSOR*  - SENSOR IDENTIFIER.
!          *CDATE*    - DATE OF MEASUREMENT (YYYYMMDDHHMMSS).
!          *IDENTI*   - SATELLITE IDENTIFIER.
!          *RLAT*     - LATITUDE OF MEASUREMENT (DEGREE).
!          *RLON*     - LONGITUDE OF MEASUREMENT (DEGREES).
!          *IDES_HS*  - DATA DESCRIPTOR FOR WAVE HEIGHT.
!          *SWH*      - WAVE HEIGHT.
!          *IDES_WS*  - DATA DESCRIPTOR FOR WIND SPEED.
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

! make sure that the size of the arrays are consistent with what is used
! in the altiemter data preprocessing
#include <parameters.h>

      REAL pvalu(jpelem)
      INTEGER iflag(jpflan)

      CHARACTER* 3 yofid
      CHARACTER*12 CDTPRO, CDTGEO, CDTGE1, ZERO, CTIME
      CHARACTER*14 CDATE
      PARAMETER(ZERO=' ')
      LOGICAL EOFD, FLAGGED, LLTEST
      LOGICAL:: LLONG=.TRUE.

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

      IF (LLONG) THEN
        READ(IUME,END=4000,IOSTAT=IRS) CDATE,
     &                                (pvalu(j),j=1,jpelem),
     &                                (iflag(j),j=1,jpflan)
        IF (IRS /= 0 ) THEN
          PVALU=1.7E38
          BACKSPACE (IUME)
          READ(IUME,END=4000,IOSTAT=JRS) CDATE,
     &                                  (pvalu(j),j=1,28),
     &                                  (iflag(j),j=1,jpflan)
          LLONG = .FALSE.
          IF (JRS /= 0 ) THEN
            WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++'
            WRITE(IU06,*) ' +                               +'
            WRITE(IU06,*) ' +    ERROR IN --READSAT--       +'
            WRITE(IU06,*) ' +    ==================         +'
            WRITE(IU06,*) ' +                               +'
            WRITE(IU06,*) ' +  INPUT DATA FORMAT UNKNOWN    +'
            WRITE(IU06,*) ' +  JRS = ',JRS 
            WRITE(IU06,*) ' +  PROGRAM WILL ABORT           +'
            WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++'
            CALL ABORT1
          ENDIF
        ENDIF
      ELSE
        READ(IUME,END=4000           ) CDATE,
     &                                (pvalu(j),j=1,28),
     &                                (iflag(j),j=1,jpflan)
      ENDIF

      identi = INT(pvalu(1))
      isensor = INT(pvalu(3))
      rlat = pvalu(10)
      rlon = pvalu(11)
      ides_ws = 11012
      ws   = pvalu(13)
      ides_hs = 22021
      swh  = pvalu(15)
      ifl=0
      DO j=1,jpflan
        ifl=ifl+iflag(j)
      ENDDO
      IF ( ifl .EQ. 0 ) flagged = .FALSE.

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
