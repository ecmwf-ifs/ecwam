      SUBROUTINE READSAT (IU06, IUME, CDTPRO, IDENTI, ISENSOR,
     &                    CDATE, RLAT, RLON, IDES_HS, SWH,
     &                    IDES_WS, WS, LLEOFD, YOFID)

!--------------------------------------------------------------------

!**** *READSAT* - READ RADAR ALTIMETER DATA.

!     H. GUNTHER     GKSS/ECMWF       AUGUST 1991

!     PURPOSE.
!     --------

!       ACCESS RADAR ALTIMETER REDUCED QUALITY-CONTROL PASSED
!       DATA FILE (QCPYYYYMMDDHH0000) AND READ ONE MEASUREMENT.

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
!          *LLEOFD*     - END OF DATA. (.TRUE. , .FALSE.)
!          *YOFID*    - CHARACTER*3   FILE ID.


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

!        FEB 2003: JEAN BIDLOT
!                  ALLOW MORE THAN ONE ALTIMETER

!        JUL 2003: SALEH ABDALLA
!                  CHANGE INPUT FILE FROM FLAGGED ALTIMETER DATA (RFL)
!                  TO REDUCED QUALITY-CONTROL PASSED DATA (QCP)

!-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IU06, IUME
      INTEGER, INTENT(INOUT) :: IDENTI, ISENSOR, IDES_HS, IDES_WS

      REAL, INTENT(INOUT) :: RLAT, RLON, SWH, WS

      CHARACTER(LEN=3), INTENT(IN) :: YOFID
      CHARACTER(LEN=14), INTENT(IN) :: CDTPRO
      CHARACTER(LEN=14), INTENT(INOUT) :: CDATE

      LOGICAL, INTENT(INOUT) :: LLEOFD


      INTEGER :: IFAIL, IRS, JRS

      CHARACTER(LEN=12) :: CTIME
      CHARACTER(LEN=14) :: CDTGE1

      CHARACTER(LEN=14), PARAMETER :: ZERO='              '

      LOGICAL :: LLTEST
      LOGICAL :: LLONG=.TRUE.

      DATA CDTGE1 /ZERO/

! ----------------------------------------------------------------------

!*    1. INITIALIZE.
!        -----------

      IF (CDTGE1.EQ.ZERO) THEN
        CDTGE1 = CDTPRO
        LLEOFD = .TRUE.
      ENDIF

! ----------------------------------------------------------------------

!*    2. GENERATE FIRST FILE DATES AND FETCH FILE FROM ECFILE.
!        -----------------------------------------------------

      IF (LLEOFD) THEN
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
          WRITE(IU06,*) ' +                CDTGE1 = ', CDTGE1
          WRITE(IU06,*) ' + IFAIL = ', IFAIL
          WRITE(IU06,*) ' + LLEOFD  = ', LLEOFD
          WRITE(IU06,*) ' + PROCESSING WILL BE CONTINUED          +'
          WRITE(IU06,*) ' +                                       +'
          WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++'
          CALL FLUSH(IU06)
          CDTGE1 = ZERO
          RETURN
        ENDIF
        LLEOFD = .FALSE.
        READ (IUME,END=4000) CTIME
      ENDIF

! ----------------------------------------------------------------------

!*    3. READ THE MEASUREMENT.
!        ---------------------

3000  CONTINUE
      READ(IUME,END=4000,IOSTAT=IRS) CDATE,
     &                               IDENTI, ISENSOR,
     &                               RLAT, RLON,
     &                               IDES_WS, WS,
     &                               IDES_HS, SWH
      IF (IRS /= 0 ) THEN
        BACKSPACE (IUME)
        READ(IUME,END=4000,IOSTAT=JRS) CTIME
        IF (JRS /= 0 ) THEN
          WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++'
          WRITE(IU06,*) ' +                               +'
          WRITE(IU06,*) ' +    ERROR IN --READSAT--       +'
          WRITE(IU06,*) ' +    ==================         +'
          WRITE(IU06,*) ' +                               +'
          WRITE(IU06,*) ' +  INPUT DATA FORMAT UNKNOWN    +'
          WRITE(IU06,*) ' +  IRS = ',IRS 
          WRITE(IU06,*) ' +  PROGRAM WILL ABORT           +'
          WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++'
          CALL ABORT1
        ELSE
!         THERE IS MORE THAN ONE ALTIMETER
          GOTO 3000
        ENDIF
      ENDIF

      RETURN

! ----------------------------------------------------------------------

!*    4. END OF FILE, FETCH NEXT FILE AND TRY AGAIN.
!        -------------------------------------------

 4000 CONTINUE
      LLEOFD = .TRUE.
      CDTGE1 = ZERO

      RETURN
      END SUBROUTINE READSAT
