!=======================================================================
      SUBROUTINE WGRIB2FDB (KUSO, KTEST, 
     &                      IGRIB_HANDLE, KLEN, KGRIB,
     &                      CDFDBSF, KFDB, LFDBOPEN,
     &                      IMDLGRBID_G,IMDLGRBID_M,
     &                      KERR)
!=======================================================================
!**
!**** NAME    *WGRIB2FDB*
!**** ----
!**
!**   PURPOSE
!**   -------
!**      WRITE FIELD TO THE FIELD DATA BASE.
!**
!**   INTERFACE
!**   ---------
!**      INPUT:*KUSO*         LOGICAL UNIT FOR STANDARD OUTPUT.
!**            *KTEST*        SWITCH DIAGNOSTICS OUTPUT ON IF KTEST GT 1.
!**            *IGRIB_HANDLE* GRIB HANDLE CONTAINING THE ENCODED DATA.
!**            *KLEN*         GRIB MESSAGE LENGTH.
!**            *KGRIB*        GRIB MESSAGE.
!**            *CDFDBSF*
!**            *KFDB*         DATA BASE REFERENCE
!**            *LFDBOPEN*     TRUE IF DATA BASE IS OPENED.
!**            *KERR*         0 IF NO ERRORS ENCOUNTERED.
!**
!**   METHOD
!**   ------
!**      FIELD DATA BASE ROUTINES ARE USED TO OPEN THE FIELD
!**      DATA BASE AND FILE, AND TO WRITE/READ THE FIELD INTO/
!**      FROM THE FIELD DATA BASE.
!**
!**
!**   AUTHOR
!**   ------
!**     JEAN BIDLOT   ECMWF   AUGUST 2009. 

!     ------------------------------------------------------------------

      USE YOWGRIBHD, ONLY : NWINOFF
      USE YOWMPP   , ONLY : NPRECI

      USE FDBSUBS_MOD
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

!     ------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: KUSO, KTEST, IGRIB_HANDLE
      INTEGER, INTENT(IN) :: KLEN 
      INTEGER, INTENT(IN) :: IMDLGRBID_M, IMDLGRBID_G 
      INTEGER, INTENT(INOUT) :: KFDB
      INTEGER, INTENT(OUT) :: KERR
      INTEGER, DIMENSION(KLEN), INTENT(INOUT) :: KGRIB
      CHARACTER(LEN=*), INTENT(INOUT) :: CDFDBSF
      LOGICAL, INTENT(INOUT) :: LFDBOPEN

      INTEGER :: ICDFDBSF
      INTEGER :: ISTREAM, IMDL, ISTEP, ILVTP, ILEN
      INTEGER :: NLOCGRB, NENSFNB, NLEG, ICENTRE
      INTEGER :: ISTAT, ISTATUS
      INTEGER(KIND=JPKSIZE_T) :: KBYTES
      INTEGER, SAVE :: IFILE_HANDLE = -999

      REAL :: ZHOOK_HANDLE

      CHARACTER(LEN=1) :: CLDOMAIN
      CHARACTER(LEN=1),save :: CLDOMAIN_PREV = " "
      CHARACTER(LEN=1) :: CLEG

      CHARACTER(LEN=2) :: CLTYPE
      CHARACTER(LEN=2),save :: CLTYPE_PREV = "  "

      CHARACTER(LEN=2) :: CLCLASS
      CHARACTER(LEN=2),save :: CLCLASS_PREV ="  "


      CHARACTER(LEN=2) :: CLDIR
      CHARACTER(LEN=2),save :: CLDIR_PREV ="  "

      CHARACTER(LEN=2) :: CLFRE
      CHARACTER(LEN=2),save :: CLFRE_PREV="  "

      CHARACTER(LEN=3) :: CLPARAM
      CHARACTER(LEN=3),save :: CLPARAM_PREV="   "

      CHARACTER(LEN=3) :: CLNUMBER

      CHARACTER(LEN=3) :: CLWINOFF
      CHARACTER(LEN=3),save :: CLWINOFF_PREV="   "

      CHARACTER(LEN=4) :: CLTIME
      CHARACTER(LEN=4),save :: CLTIME_PREV ="    "

      CHARACTER(LEN=4) :: CLEXPVER
      CHARACTER(LEN=4),save :: CLEXPVER_PREV ="    "

      CHARACTER(LEN=4) :: CLSTREAM
      CHARACTER(LEN=4),save :: CLSTREAM_PREV ="    "

      CHARACTER(LEN=4) :: CDORIGIN

      CHARACTER(LEN=5) :: CLSYSTEM

      CHARACTER(LEN=5) :: CLMETHOD

      CHARACTER(LEN=6) :: CLSTEP
      CHARACTER(LEN=6),save :: CLSTEP_PREV = "      "

      CHARACTER(LEN=8) :: CLDATE
      CHARACTER(LEN=8),save :: CLDATE_PREV = "        "

      CHARACTER(LEN=8) :: HDATE

      CHARACTER(LEN=8) :: CDATEREF

      CHARACTER(LEN=12) :: C12

      CHARACTER(LEN=1),save :: CLLEVTY_PREV=" "
      CHARACTER(LEN=2),save :: CLREPR_PREV ="  "

      CHARACTER(LEN=64) :: CLFILE

!     ------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('WGRIB2FDB',0,ZHOOK_HANDLE)
#endif
!*    1.   DERIVE FIELD DATA BASE VARIABLES FROM GRIB_HANDLE

!     ------------------------------------------------------------------

!for debugging
!!      CLFILE = 'allfields_%p.grib'
!!      CALL EXPAND_STRING(IRANK,NPROC,0,0,CLFILE,1)
!!      IF(IFILE_HANDLE == -999) THEN
!!        CALL IGRIB_OPEN_FILE(IFILE_HANDLE,TRIM(CLFILE),'w')
!!      ENDIF


!     THIS VALUE WILL BE RETURNED IF EVERYTHING GOES OK.
      KERR = 0

      ILEN = KLEN
!*    OPEN FDB (only once)
      IF(.NOT.LFDBOPEN)  THEN
        CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'stream',ISTREAM)
        CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'localDefinitionNumber',
     &                       NLOCGRB)

!       OPEN THE APPROPRIATE FDB:
        IF (ISTREAM == 1082 .OR.
     &      ISTREAM == 1222 .OR.
     &      ISTREAM == 1095 .OR.
     &      ISTREAM == 1203 .OR.
     &      ISTREAM == 1204) THEN
          ISTAT = IOPENFDBSUBS( 'seas', KFDB, 'w' )
        ELSE
          ISTAT = IOPENFDBSUBS( 'fdb', KFDB, 'w' )
        ENDIF
        IF ( ISTAT /= 0 ) THEN
          WRITE(kuso,'("Error\ /WGRIB2FDB/ iopenfdb return status:",
     &    I3)') ISTAT
          KERR = 1
          RETURN
        ELSEIF (KTEST > 1) THEN
           WRITE(kuso,'("\ /WGRIB2FDB/ iopenfdb status:", i3)') istat
        ENDIF
        LFDBOPEN = .TRUE.

!       DEFINE FDB ROOT DIR:
        ICDFDBSF = LEN_TRIM(CDFDBSF)
        IF (ICDFDBSF > 0 ) THEN
          ISTAT = ISET_FDBSUBS_ROOT(KFDB, CDFDBSF(1:ICDFDBSF) )
          IF ( ISTAT /= 0 ) THEN
            WRITE(KUSO,'("Error\ /WGRIB2FDB/ Root dir specified: ",a)')
     &      CDFDBSF(1:ICDFDBSF)
            WRITE(kuso,'("Error\ /WGRIB2FDB/ iopenfdb return status:",
     &      I3)') ISTAT
            KERR = 1
            RETURN
          ELSEIF (KTEST > 1) THEN
            WRITE(KUSO,'("\ /WGRIB2FDB/ iset_fdb_root status:", i3,
     &      " Root: ", a25)') ISTAT, CDFDBSF(1:ICDFDBSF)
          ENDIF
        ENDIF

!       ENSEMBLE NUMBER AND LEG (ONLY SET WHEN NEEDED !):
        IF(ISTREAM == 1081 .OR.
     &     ISTREAM == 1083 .OR.
     &     ISTREAM == 1084 .OR.
     &     ISTREAM == 1078 .OR.
     &     ISTREAM == 1079 .OR.
     &     ISTREAM == 1086 ) THEN

          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'perturbationNumber',
     &                         NENSFNB)
          WRITE( CLNUMBER, '(I3.3)' ) NENSFNB
          IF ( NENSFNB > 0 ) THEN
            ISTAT = ISETVALFDBSUBS( KFDB, 'number', CLNUMBER)
          ENDIF

          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'legNumber',NLEG,ISTATUS)
          IF(ISTATUS == 0 ) THEN
            IF( NLEG > 0 ) THEN
              WRITE( CLEG, '(I1.1)' ) NLEG 
              ISTAT = ISETVALFDBSUBS( KFDB, 'leg', CLEG)
            ENDIF
          ENDIF

        ELSE IF(ISTREAM == 1203 .OR. ISTREAM == 1204) THEN
          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'perturbationNumber',
     &                         NENSFNB)
          WRITE( CLNUMBER, '(I3.3)' ) NENSFNB
          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'type',C12)
          CLTYPE=C12(1:2)
          IF ( NENSFNB > 0 .OR. CLTYPE == 'an') THEN
            ISTAT = ISETVALFDBSUBS( KFDB, 'number', CLNUMBER)
          ENDIF
        ELSE IF(ISTREAM == 1082 .OR.
     &          ISTREAM == 1088 .OR.
     &          ISTREAM == 1250 .OR.
     &          ISTREAM == 1095 .OR.
     &          ISTREAM == 1222) THEN
          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'perturbationNumber',
     &                         NENSFNB)
          WRITE( CLNUMBER, '(I3.3)' ) NENSFNB
          ISTAT = ISETVALFDBSUBS( KFDB, 'number', CLNUMBER)
        ENDIF
!       SYSTEM AND METHOD:
        IF (ISTREAM == 1082 .OR.
     &      ISTREAM == 1095 .OR.
     &      ISTREAM == 1222 .OR.
     &      ISTREAM == 1203 .OR.
     &      ISTREAM == 1204) THEN
!!        if the class is not od then system is set to off
          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'class',C12)
          CLCLASS=C12(1:2)
          IF(CLCLASS /= 'od' ) THEN
            CLSYSTEM='off'
          ELSE
            CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'system',C12)
            CLSYSTEM=C12(1:5)
          ENDIF
          ISTAT = ISETVALFDBSUBS( KFDB, 'system', CLSYSTEM)

          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'method',C12)
          CLMETHOD=C12(1:5)
          ISTAT = ISETVALFDBSUBS( KFDB, 'method', CLMETHOD)
          ISTAT = ISETVALFDBSUBS( KFDB, 'level','0000')

          IF (KTEST.GT.3) WRITE(KUSO,'("\ /WGRIB2FDB/ system:", a5,
     &                    " method ", a5)') CLSYSTEM,CLMETHOD 
        ENDIF

!       CENTRE ORIGIN :
        IF (NLOCGRB == 18) THEN
          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'consensusCount',ICENTRE,
     &                         ISTATUS)
          IF(ISTATUS /= 0 ) ICENTRE=0

          IF(ICENTRE == 0 ) then
          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'origin',C12)
            CDORIGIN = C12(1:4)
            CALL WAM_U2L1CR(CDORIGIN)
            ISTAT = ISETVALFDBSUBS( KFDB, 'origin', CDORIGIN )
          ELSE
            ISTAT = ISETVALFDBSUBS( KFDB, 'origin', 'consensus' )
          ENDIF
        ELSE IF (NLOCGRB == 23) THEN
          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'centre',ICENTRE)
          IF(ICENTRE == 98) then
            CDORIGIN = 'ecmf'
          ELSE
            WRITE(kuso,'( "Error\ /WGRIB2FDB/ origin unknown ", i4)' )
     &            ICENTRE 
            KERR = 1
            RETURN
          ENDIF
          ISTAT = ISETVALFDBSUBS( KFDB, 'origin', CDORIGIN )
        ENDIF

!       HINDCAST REFERENCE DATE:

        IF ( ISTREAM == 1204 ) THEN
          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'refdate',C12)
          CDATEREF=c12(1:8)
          ISTAT = ISETVALFDBSUBS(KFDB, 'refdate', CDATEREF)
        ELSE IF(ISTREAM == 1084 .OR.
     &          ISTREAM == 1085 ) THEN
          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'referenceDate',C12)
          CDATEREF=c12(1:8)
          ISTAT = ISETVALFDBSUBS(KFDB, 'refdate', CDATEREF)
        ELSE IF (ISTREAM == 1078 .OR.
     &           ISTREAM == 1079 ) THEN
          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'dataDate',C12)
          HDATE=C12(1:8)
          ISTAT = ISETVALFDBSUBS(KFDB,'hdate', HDATE)
        ENDIF

      ENDIF

!     ------------------------------------------------------------------

!*    2. MAKE DATA BASE NAME.

!     SET:

!     STREAM:
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'stream',C12)
      CLSTREAM=C12(1:4)
      if(clstream.ne.clstream_prev) then
        ISTAT = ISETVALFDBSUBS( KFDB, 'stream', CLSTREAM )
        clstream_prev=clstream
      endif

      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'stream',ISTREAM)
!     CLASS:
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'class',C12)
      CLCLASS=C12(1:2)
      if(clclass.ne.clclass_prev) then
        ISTAT = ISETVALFDBSUBS( KFDB, 'class',  CLCLASS )
        clclass_prev=clclass
      endif

!     EXPVER:
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'expver',C12)
      CLEXPVER=C12(1:4)
      if(clexpver.ne.clexpver_prev) then
        ISTAT = ISETVALFDBSUBS( KFDB, 'expver', CLEXPVER )
        clexpver_prev=clexpver
      endif

!     DOMAIN:
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'generatingProcessIdentifier',
     &                     IMDL)
      IF(IMDL == IMDLGRBID_G) THEN
        CLDOMAIN = 'g' 
      ELSEIF(IMDL == IMDLGRBID_M) THEN
        CLDOMAIN = 'm' 
      ELSE
        WRITE(kuso,'( "Error\ /WGRIB2FDB/ domain unknown ", i4)' )IMDL
        KERR = 1
        RETURN
      ENDIF
      if(cldomain.ne.cldomain_prev) then
        ISTAT = ISETVALFDBSUBS( KFDB, 'domain', CLDOMAIN )
        cldomain_prev=cldomain
      endif

!     REPRESENTATION:
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'gridType',C12)
      if(c12(9:10).ne.clrepr_prev) then
        ISTAT = ISETVALFDBSUBS( KFDB, 'repr', c12(9:10) )
        clrepr_prev=c12(9:10)
      endif

!     TYPE:
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'type',C12)
      CLTYPE=C12(1:2)
      if(CLTYPE.ne.cltype_prev) then
        ISTAT = ISETVALFDBSUBS( KFDB, 'type', CLTYPE )
        cltype_prev=c12(1:2)
      endif

!     PARAMETER:
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'param',C12)
      CLPARAM=C12(1:3)
      if(CLPARAM.ne.clparam_prev) then
        ISTAT = ISETVALFDBSUBS( KFDB, 'param', CLPARAM )
        clparam_prev=CLPARAM
      endif

!     DATE:
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'dataDate',C12)
      CLDATE=C12(1:8)
!!!   For VAREPS hindcast we need to reset CLDATE to the reference date
!!!   if non zero
      C12='000000000000'
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'referenceDate',C12,ISTATUS)
      IF(ISTATUS == 0 .AND. C12(1:1).NE. '0' ) THEN
        CLDATE=C12(1:8)
      ENDIF
      if(CLDATE.ne.cldate_prev) then
        ISTAT = ISETVALFDBSUBS( KFDB, 'date', CLDATE )
        cldate_prev=cldate
      endif

!     TIME:
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'time',C12)
      CLTIME=C12(1:4)
      if(CLTIME.ne.cltime_prev) then
        ISTAT = ISETVALFDBSUBS( KFDB, 'time', CLTIME )
        cltime_prev=cltime
      endif
!     FC STEP:
      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'stepUnits','h')
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'endStep',ISTEP)
!     make sure the step is in hours
      ISTEP=ISTEP
      WRITE(CLSTEP, '( I6.6 )' ) ISTEP 
      if(CLSTEP.ne.clstep_prev) then
        ISTAT = ISETVALFDBSUBS( KFDB, 'step', CLSTEP )
        clstep_prev=CLSTEP
      endif

!     LEVEL TYPE:
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'levtype',ilvtp)
      IF(ILVTP == 209) THEN
        if('wv'.ne.cllevty_prev) then
          ISTAT = ISETVALFDBSUBS( KFDB, 'levty', 'wv'  )
          cllevty_prev='wv'
        endif
      ELSE IF(ILVTP == 212) THEN
        if('ws'.ne.cllevty_prev) then
          ISTAT = ISETVALFDBSUBS( KFDB, 'levty', 'ws'  )
          cllevty_prev='ws'
        endif
      ELSE
        if('s'.ne.cllevty_prev) then
          ISTAT = ISETVALFDBSUBS( KFDB, 'levty', 's'   )
          cllevty_prev='s'
        endif
      ENDIF 

!     DIRECTION AND FREQUENCY:
      IF(CLPARAM == '251' ) THEN
         CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'directionNumber',C12)
         CLDIR=C12(1:2)
         if(CLDIR.ne.cldir_prev) then
           ISTAT = ISETVALFDBSUBS( KFDB, 'direction', CLDIR )
           cldir_prev=CLDIR
         endif
         CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'frequencyNumber',C12)
         CLFRE=C12(1:2)
         if(CLFRE.ne.clfre_prev) then
           ISTAT = ISETVALFDBSUBS( KFDB, 'frequency', CLFRE )
           clfre_prev=CLFRE
         endif
      ELSE
        CLDIR="00"
        CLFRE="00"
        IF (ISTREAM /= 1082 .AND.
     &      ISTREAM /= 1095 .AND.
     &      ISTREAM /= 1203 .AND.
     &      ISTREAM /= 1204 .AND.
     &      ISTREAM /= 1222) THEN
         if(CLDIR.ne.cldir_prev) then
            ISTAT = ISETVALFDBSUBS( KFDB, 'direction', CLDIR )
            cldir_prev=CLDIR
         endif
         if(CLFRE.ne.clfre_prev) then
            ISTAT = ISETVALFDBSUBS( KFDB, 'frequency', CLFRE )
            clfre_prev=CLFRE
         endif
        ENDIF
      ENDIF

      IF (KTEST.GT.3) THEN
        WRITE(kuso,'("\ /WGRIB2FDB/ CLASS: ", a8)') CLCLASS
        WRITE(kuso,'("\ /WGRIB2FDB/ EXPVER:", a8)') CLEXPVER
        WRITE(kuso,'("\ /WGRIB2FDB/ DOMAIN:", a8)') CLDOMAIN
        WRITE(kuso,'("\ /WGRIB2FDB/ TY:    ", a8)') CLTYPE
        WRITE(kuso,'("\ /WGRIB2FDB/ PARAM: ", a8)') CLPARAM
        WRITE(kuso,'("\ /WGRIB2FDB/ DATE:  ", a8)') CLDATE
        WRITE(kuso,'("\ /WGRIB2FDB/ TIME:  ", a8)') CLTIME
        WRITE(kuso,'("\ /WGRIB2FDB/ STEP:  ", a8)') CLSTEP
        WRITE(kuso,'("\ /WGRIB2FDB/ DIR:   ", a8)') CLDIR
        WRITE(kuso,'("\ /WGRIB2FDB/ FREQ:  ", a8)') CLFRE
        CALL FLUSH (KUSO)
      ENDIF
! LONG WINDOW OFFSET
        IF (ISTREAM == 1250 .OR. ISTREAM == 1248 ) THEN
          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'type',C12)
          CLTYPE=C12(1:2)
!          IF (CLTYPE /= 'an') THEN
            CALL IGRIB_GET_VALUE(IGRIB_HANDLE,
     &      'offsetToEndOf4DvarWindow',NWINOFF)
!          ENDIF
          WRITE(CLWINOFF,'(I3.3)') NWINOFF
          if(CLWINOFF.ne.clwinoff_prev) then
            ISTAT = ISETVALFDBSUBS( KFDB, 'anoffset', CLWINOFF)
            clwinoff_prev=clwinoff
          endif
        ENDIF

!     ------------------------------------------------------------------

!*    3.  GET AND WRITE FIELD TO FDB

!     ------------------------------------------------------------------

!for debugging
!!      write(*,*) 'also writting to file the data destined to fdb !!!'
!!      CALL IGRIB_GET_MESSAGE_SIZE(IGRIB_HANDLE,KBYTES)
!!      CALL IGRIB_WRITE_BYTES(IFILE_HANDLE,KGRIB,KBYTES)


      ISTAT = IWRITEFDBSUBS( KFDB, KGRIB, ILEN ) 
      IF ( ISTAT .NE. 0 ) THEN
        WRITE(KUSO, '( "Error\ /WGRIB2FDB/ Failed to write to fdb")' )
        CALL GSTATS(1787,0)
        ISTAT = ICLOSEFDBSUBS( KFDB )
        CALL GSTATS(1787,1)
        KERR = 1
        RETURN
      ENDIF
      IF (KTEST.GT.1) WRITE(KUSO,'("\ /WGRIB2FDB/ iwritefdb ilen:", i10,
     &                      " status:", i4)') ILEN, ISTAT

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('WGRIB2FDB',1,ZHOOK_HANDLE)
#endif
      END SUBROUTINE WGRIB2FDB
