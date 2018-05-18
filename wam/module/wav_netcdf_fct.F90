MODULE WAV_NETCDF_FCT
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal
      TYPE TIMEPERIOD
        INTEGER(KIND=JWIM) :: nbTime
        REAL(KIND=JWRU), allocatable :: ListTime(:)
      END TYPE TIMEPERIOD
      CONTAINS
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE TEST_FILE_EXIST_DIE(eFile, errmsg)
      IMPLICIT NONE
      CHARACTER(LEN=*) :: eFile
      CHARACTER(LEN=*) :: errmsg
      LOGICAL :: LFLIVE
      INQUIRE( FILE = TRIM(eFile), EXIST = LFLIVE )
      IF ( .NOT. LFLIVE ) THEN
        Print *, 'Missing file =', TRIM(eFile)
        Print *, 'Reason for call =', TRIM(errmsg)
        STOP
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE FIND_MATCH_TIME(RecTime, eTime, iTime1, w1, iTime2, w2)
      IMPLICIT NONE
      TYPE(TIMEPERIOD), intent(in) :: RecTime
      REAL(KIND=JWRU), intent(in) :: eTime
      INTEGER(KIND=JWIM), intent(out) :: iTime1, iTime2
      real(KIND=JWRU), intent(out) :: w1, w2
      REAL(KIND=JWRU), parameter :: tolDay = 0.00000001
      REAL(KIND=JWRU) :: DeltaTime
      DO iTime2=2,RecTime % nbTime
        iTime1=iTime2-1
        DeltaTime=RecTime % ListTime(iTime2) - RecTime % ListTime(iTime1)
        w1=(RecTime % ListTime(iTime2) - eTime) / DeltaTime
        w2=(eTime - RecTime % ListTime(iTime1)) / DeltaTime
        !        WRITE(740+MyRankGlobal,*) 'iTime1=', iTime1, ' iTime2=', iTime2
        !        WRITE(740+MyRankGlobal,*) 'eTime=', eTime
        !        WRITE(740+MyRankGlobal,*) 'DeltaTime=', DeltaTime
        !        WRITE(740+MyRankGlobal,*) 'ListTime(iTime1)=', ListTime(iTime1)
        !        WRITE(740+MyRankGlobal,*) 'ListTime(iTime2)=', ListTime(iTime2)
        !        WRITE(740+MyRankGlobal,*) 'w1=', w1, ' w2=', w2
        IF ((w1 + tolDay .ge. 0.).and.(w2 + tolDay .ge. 0.)) THEN
          RETURN
        END IF
      END DO
      WRITE(740+MyRankGlobal,*) 'JWRU=', JWRU
      WRITE(740+MyRankGlobal,*) 'We did not find the time'
      WRITE(740+MyRankGlobal,*) 'nbTime=', RecTime % nbTime
      WRITE(740+MyRankGlobal,*) 'eTime=', eTime
      WRITE(740+MyRankGlobal,*) 'ListTime(1)     = ', RecTime % ListTime(1)
      WRITE(740+MyRankGlobal,*) 'ListTime(nbTime) = ', RecTime % ListTime(RecTime % nbTime)
      FLUSH(740+MyRankGlobal)
      STOP
      END SUBROUTINE FIND_MATCH_TIME
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_LIST_TIME(ncid, RecTime)
      USE NETCDF
      IMPLICIT NONE
      INTEGER(KIND=JWIM), intent(in) :: ncid
      TYPE(TIMEPERIOD), intent(inout) :: RecTime
      character(len=*), parameter :: CallFct = "READ_LIST_TIME"
      REAL(KIND=JWRU), allocatable :: ListTime_mjd(:)
      character (len=100) :: eStrUnitTime
      REAL(KIND=JWRU) :: ConvertToDay, eTimeStart, eTimeBnd
      REAL(KIND=JWRU) :: eTime
      INTEGER(KIND=JWIM), dimension(nf90_max_var_dims) :: dimids
      INTEGER(KIND=JWIM) :: nbtime_mjd
      INTEGER(KIND=JWIM) :: istat
      INTEGER(KIND=JWIM) :: iTime
      INTEGER(KIND=JWIM) :: varid
      !
      ! reading the time
      !
      ISTAT = nf90_inq_varid(ncid, "ocean_time", varid)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)
      !
      ISTAT = nf90_get_att(ncid, varid, "units", eStrUnitTime)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)
      CALL CF_EXTRACT_TIME(eStrUnitTime, ConvertToDay, eTimeStart)
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'eStrUnitTime = ', eStrUnitTime
      WRITE(740+MyRankGlobal,*) 'ConvertToDay = ', ConvertToDay
      WRITE(740+MyRankGlobal,*) '  eTimeStart = ', eTimeStart
#endif
      !
      ISTAT = nf90_inquire_variable(ncid, varid, dimids=dimids)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)
      !
      ISTAT = nf90_inquire_dimension(ncid, dimids(1), len=nbtime_mjd)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)
      allocate(ListTime_mjd(nbtime_mjd), stat=istat)
      !
      ISTAT = nf90_get_var(ncid, varid, ListTime_mjd)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)
      !
      RecTime % nbTime = nbtime_mjd
      allocate(RecTime % ListTime(nbtime_mjd), stat=istat)
      DO iTime=1,nbtime_mjd
         eTime = ListTime_mjd(iTime)*ConvertToDay + eTimeStart
         RecTime % ListTime(iTime) = eTime
      END DO
      deallocate(ListTime_mjd)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAV_GENERIC_NETCDF_ERROR(CallFct, idx, iret)
      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal
      USE NETCDF
      implicit none
      INTEGER(KIND=JWIM), intent(in) :: iret, idx
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
      REAL(KIND=JWRU), intent(out) :: eTimeDay
      REAL(KIND=JWRU), intent(in) :: Increment
      REAL(KIND=JWRU) :: eJD
      character(len=4) eYear
      character(len=2) eMonth, eDay, eHour, eMin, eSec
      INTEGER(KIND=JWIM) :: year, month, day, hour, min, sec
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
      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal
      IMPLICIT NONE
      character(len=100), intent(in) :: eStrUnitTime
      real(KIND=JWRU), intent(out) :: ConvertToDay, eTimeStart
      character(len=100) :: Xname, Yname
      character(len=10) :: YnameYear, YnameMonth, YnameDay
      character(len=10) :: YnameHour, YnameMin, YnameSec
      character(len=50) :: YnameB, YnameD, YnameE
      character(len=50) :: YnameDate, YnameTime, YnameTimeP
      character(len=15) :: eStrTime
      INTEGER(KIND=JWIM) :: alenB, alenC, alenD, alenE, alenTime, alenDate
      INTEGER(KIND=JWIM) :: alen, posBlank
      INTEGER(KIND=JWIM) :: lenHour, lenMin, lenSec, lenMonth, lenDay, posSepDateTime
      alen=LEN_TRIM(eStrUnitTime)
!      WRITE(740+MyRankGlobal,*) 'alen=', alen
      posBlank=INDEX(eStrUnitTime(1:alen), ' ')
!      WRITE(740+MyRankGlobal,*) 'posBlank=', posBlank
      Xname=eStrUnitTime(1:posBlank-1) ! should be days/hours/seconds
!      WRITE(740+MyRankGlobal,*) 'Xname=', Xname
      IF (TRIM(Xname) .eq. 'days') THEN
        ConvertToDay=REAL(1,JWRU)
      ELSEIF (TRIM(Xname) .eq. 'hours') THEN
        ConvertToDay=REAL(1,JWRU)/REAL(24,JWRU)
      ELSEIF (TRIM(Xname) .eq. 'seconds') THEN
        ConvertToDay=REAL(1,JWRU)/REAL(86400,JWRU)
!        WRITE(740+MyRankGlobal,*) 'assignment here for seconds'
      ELSE
        Print *, 'Error in the code for conversion'
        STOP
      END IF
!      WRITE(740+MyRankGlobal,*) 'ConvertToDay=', ConvertToDay
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
      WRITE(740+MyRankGlobal,*) 'eStrTime=', eStrTime
      CALL CT2MJD(eStrTime, eTimeStart)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
END MODULE
