MODULE WAV_NETCDF_FCT
      CONTAINS
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAV_GENERIC_NETCDF_ERROR(CallFct, idx, iret)
      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal
      USE NETCDF
      implicit none
      integer, intent(in) :: iret, idx
      character(*), intent(in) :: CallFct
      character(len=500) :: CHRERR
      character(len=1500) :: CHRERR_tot
      IF (iret .NE. nf90_noerr) THEN
        CHRERR = nf90_strerror(iret)
        WRITE(CHRERR_tot,*) TRIM(CallFct), ' -', idx, '-', TRIM(CHRERR)
        WRITE(*,*) TRIM(CHRERR_tot)
        WRITE(740+MyRankGlobal,*) TRIM(CHRERR_tot)
        WRITE(740+MyRankGlobal,*) 'Debug the netcdf'
        FLUSH(740+MyRankGlobal)
        STOP 'Debug the netcdf'
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAV_GET_ETIMEDAY(eTimeDay, Increment)
      USE YOWWAMI  , ONLY : CBPLTDT
      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal
      IMPLICIT NONE
      real*8, intent(out) :: eTimeDay
      real*8, intent(in) :: Increment
      real*8 eJD
      character(len=4) eYear
      character(len=2) eMonth, eDay, eHour, eMin, eSec
      integer year, month, day, hour, min, sec
      eYear(1:1)  = CBPLTDT(1:1)
      eYear(2:2)  = CBPLTDT(2:2)
      eYear(3:3)  = CBPLTDT(3:3)
      eYear(4:4)  = CBPLTDT(4:4)
      eMonth(1:1) = CBPLTDT(5:5)
      eMonth(2:2) = CBPLTDT(6:6)
      eDay(1:1)   = CBPLTDT(7:7)
      eDay(2:2)   = CBPLTDT(8:8)
      eHour(1:1)  = CBPLTDT(9:9)
      eHour(2:2)  = CBPLTDT(10:10)
      eMin(1:1)   = CBPLTDT(11:11)
      eMin(2:2)   = CBPLTDT(12:12)
      eSec(1:1)   = CBPLTDT(13:13)
      eSec(2:2)   = CBPLTDT(14:14)
      read(eYear , '(i10)' ) year
      read(eMonth, '(i10)' ) month
      read(eDay  , '(i10)' ) day
      read(eHour , '(i10)' ) hour
      read(eMin  , '(i10)' ) min
      read(eSec  , '(i10)' ) sec
      CALL DATE_ConvertSix2mjd(year, month, day, hour, min, sec, eJD)
      eTimeDay = eJD + Increment / DBLE(86400)
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*)  'CBPLTDT=', CBPLTDT
      WRITE(740+MyRankGlobal,*)  'year=', year
      WRITE(740+MyRankGlobal,*)  'month=', month
      WRITE(740+MyRankGlobal,*)  'day=', day
      WRITE(740+MyRankGlobal,*)  'hour=', hour
      WRITE(740+MyRankGlobal,*)  'min=', min
      WRITE(740+MyRankGlobal,*)  'sec=', sec
      WRITE(740+MyRankGlobal,*)  'eJD=', eJD
      WRITE(740+MyRankGlobal,*)  'eTimeDay=', eTimeDay
      WRITE(740+MyRankGlobal,*)  'Increment=', Increment
      FLUSH(740+MyRankGlobal)
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CF_EXTRACT_TIME(eStrUnitTime, ConvertToDay, eTimeStart)
      USE PARKIND_WAVE, ONLY : JWRB, JWRU
      IMPLICIT NONE
      character(len=100), intent(in) :: eStrUnitTime
      real(KIND=JWRU), intent(out) :: ConvertToDay, eTimeStart
      character(len=100) :: Xname, Yname
      character(len=10) :: YnameYear, YnameMonth, YnameDay
      character(len=10) :: YnameHour, YnameMin, YnameSec
      character(len=50) :: YnameB, YnameD, YnameE
      character(len=50) :: YnameDate, YnameTime, YnameTimeP
      character(len=15) :: eStrTime
      integer alenB, alenC, alenD, alenE, alenTime, alenDate
      integer alen, posBlank
      integer lenHour, lenMin, lenSec, lenMonth, lenDay, posSepDateTime
      alen=LEN_TRIM(eStrUnitTime)
      posBlank=INDEX(eStrUnitTime(1:alen), ' ')
      Xname=eStrUnitTime(1:posBlank-1) ! should be days/hours/seconds
      IF (TRIM(Xname) .eq. 'days') THEN
        ConvertToDay=1
      ELSEIF (TRIM(Xname) .eq. 'hours') THEN
        ConvertToDay=1/24
      ELSEIF (TRIM(Xname) .eq. 'seconds') THEN
        ConvertToDay=1/86400
      ELSE
        Print *, 'Error in the code for conversion'
        STOP
      END IF
      !
      Yname=eStrUnitTime(posBlank+1:alen)
      alenB=LEN_TRIM(Yname)
      posBlank=INDEX(Yname(1:alenB), ' ')
      YnameB=Yname(posBlank+1:alenB) ! should be 1990-01-01 0:0:0
      !
      alenC=LEN_TRIM(YnameB)
      posSepDateTime=INDEX(YnameB(1:alenC), ' ')
      IF (posSepDateTime .gt. 0) THEN
        YnameDate=YnameB(1:posSepDateTime-1) ! should be 1990-01-01
        YnameTimeP=YnameB(posSepDateTime+1:alenC) ! should be 0:0:0
        alenC=LEN_TRIM(YnameTimeP)
        posBlank=INDEX(YnameTimeP(1:alenC), ' ')
        IF (posBlank .eq. 0) THEN
          YnameTime=YnameTimeP
        ELSE
          YnameTime=YnameTimeP(1:posBlank-1)
        END IF
      ELSE
        YnameDate=YnameB
        eStrTime(10:10)='0'
        eStrTime(11:11)='0'
        eStrTime(12:12)='0'
        eStrTime(13:13)='0'
        eStrTime(14:14)='0'
        eStrTime(15:15)='0'
      END IF
      !
      alenDate=LEN_TRIM(YnameDate)
      posBlank=INDEX(YnameDate(1:alenDate), '-')
      YnameYear=YnameDate(1:posBlank-1) ! should be 1990
      YnameD=YnameDate(posBlank+1:alenDate)
      alenD=LEN_TRIM(YnameD)
      posBlank=INDEX(YnameD(1:alenD), '-')
      YnameMonth=YnameD(1:posBlank-1) ! should be 01
      YnameDay=YnameD(posBlank+1:alenD) ! should be 01
      !
      ! year
      eStrTime( 1: 1)=YnameYear( 1: 1)
      eStrTime( 2: 2)=YnameYear( 2: 2)
      eStrTime( 3: 3)=YnameYear( 3: 3)
      eStrTime( 4: 4)=YnameYear( 4: 4)
      !
      ! month
      lenMonth=LEN_TRIM(YnameMonth)
      IF (lenMonth .eq. 2) THEN
        eStrTime( 5: 5)=YnameMonth( 1: 1)
        eStrTime( 6: 6)=YnameMonth( 2: 2)
      ELSE
        IF (lenMonth .eq. 1) THEN
          eStrTime( 5: 5)='0'
          eStrTime( 6: 6)=YnameMonth( 1: 1)
        ELSE
          Print *, 'DIE in trying to get the month'
          STOP
        END IF
      END IF
      !
      ! day
      lenDay=LEN_TRIM(YnameDay)
      IF (lenDay .eq. 2) THEN
        eStrTime( 7: 7)=YnameDay( 1: 1)
        eStrTime( 8: 8)=YnameDay( 2: 2)
      ELSE
        IF (lenDay .eq. 1) THEN
          eStrTime( 7: 7)='0'
          eStrTime( 8: 8)=YnameDay( 1: 1)
        ELSE
          Print *, 'DIE in trying to get the day'
          STOP
        END IF
      END IF
      !
      eStrTime( 9: 9)='.'
      !
      IF (posSepDateTime .gt. 0) THEN
        !
        alenTime=LEN_TRIM(YnameTime)
        posBlank=INDEX(YnameTime(1:alenTime), ':')
        YnameHour=YnameTime(1:posBlank-1) ! should be 0
        YnameE=YnameTime(posBlank+1:alenTime)
        alenE=LEN_TRIM(YnameE)
        posBlank=INDEX(YnameE(1:alenE), ':')
        YnameMin=YnameE(1:posBlank-1) ! should be 0
        YnameSec=YnameE(posBlank+1:alenE) ! should be 0
        !
        !
        ! Hour
        lenHour=LEN_TRIM(YnameHour)
        IF (lenHour .eq. 2) THEN
          eStrTime(10:10)=YnameHour( 1: 1)
          eStrTime(11:11)=YnameHour( 2: 2)
        ELSE
          IF (lenHour .eq. 1) THEN
            eStrTime(10:10)='0'
            eStrTime(11:11)=YnameHour( 1: 1)
          ELSE
            Print *, 'DIE in trying to get the hour'
            STOP
          END IF
        END IF
        !
        ! Min
        lenMin=LEN_TRIM(YnameMin)
        IF (lenMin .eq. 2) THEN
          eStrTime(12:12)=YnameMin( 1: 1)
          eStrTime(13:13)=YnameMin( 2: 2)
        ELSE
          IF (lenMin .eq. 1) THEN
            eStrTime(12:12)='0'
            eStrTime(13:13)=YnameMin( 1: 1)
          ELSE
            Print *, 'DIE in trying to get the min'
          END IF
        END IF
        !
        ! Sec
        lenSec=LEN_TRIM(YnameSec)
        IF (lenSec .eq. 2) THEN
          eStrTime(14:14)=YnameSec( 1: 1)
          eStrTime(15:15)=YnameSec( 2: 2)
        ELSE
          IF (lenSec .eq. 1) THEN
            eStrTime(14:14)='0'
            eStrTime(15:15)=YnameSec( 1: 1)
          ELSE
            Print *, 'YnameSec=', TRIM(Ynamesec)
            Print *, 'lenSec=', lenSec
            Print *, 'DIE in trying to get the sec'
            STOP
          END IF
        END IF
      END IF
      CALL CT2MJD(eStrTime, eTimeStart)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
END MODULE
