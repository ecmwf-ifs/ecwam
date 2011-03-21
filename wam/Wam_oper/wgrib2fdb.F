!=======================================================================
      SUBROUTINE WGRIB2FDB (KUSO, KTEST, IGRIB_HANDLE,
     &                      CDFDBSF, KFDB, LFDBOPEN,
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

      USE YOWGRIBHD, ONLY : IMDLGRBID_G,IMDLGRBID_M
      USE YOWMPP   , ONLY : NPRECI

      USE FDBSUBS_MOD
      USE GRIB_API_INTERFACE
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

!     ------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER ::  KUSO, KTEST,IGRIB_HANDLE, KERR, KFDB
      INTEGER :: KBYTES, KLEN
      INTEGER :: ICDFDBSF
      INTEGER :: ISTREAM, IMDL, ISTEP, ILVTP
      INTEGER :: NLOCGRB, NENSFNB, NLEG, ICENTRE
      INTEGER :: ISTAT, ISTATUS
      INTEGER, ALLOCATABLE :: KGRIB(:)

      REAL :: ZHOOK_HANDLE

      CHARACTER(LEN=*) :: CDFDBSF

      CHARACTER(LEN=1) :: CLDOMAIN
      CHARACTER(LEN=1) :: CLEG

      CHARACTER(LEN=2) :: CLTYPE
      CHARACTER(LEN=2) :: CLCLASS
      CHARACTER(LEN=2) :: CLDIR
      CHARACTER(LEN=2) :: CLFRE

      CHARACTER(LEN=3) :: CLPARAM
      CHARACTER(LEN=3) :: CLNUMBER

      CHARACTER(LEN=4) :: CLTIME
      CHARACTER(LEN=4) :: CLEXPVER
      CHARACTER(LEN=4) :: CLSTREAM
      CHARACTER(LEN=4) :: CDORIGIN

      CHARACTER(LEN=5) :: CLSYSTEM
      CHARACTER(LEN=5) :: CLMETHOD

      CHARACTER(LEN=6) :: CLSTEP

      CHARACTER(LEN=8) :: CLDATE
      CHARACTER(LEN=8) :: HDATE
      CHARACTER(LEN=8) :: CDATEREF

      CHARACTER(LEN=12) :: C12

      LOGICAL :: LFDBOPEN

!     ------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WGRIB2FDB',0,ZHOOK_HANDLE)
  
!*    1.   DERIVE FIELD DATA BASE VARIABLES FROM GRIB_HANDLE

!     ------------------------------------------------------------------

!     THIS VALUE WILL BE RETURNED IF EVERYTHING GOES OK.
      KERR = 0

!*    OPEN FDB (only once)


      IF(.NOT.LFDBOPEN)  THEN

!       OPEN THE APPROPRIATE FDB:

        CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'stream',ISTREAM)
        CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'localDefinitionNumber',
     &                       NLOCGRB)

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
            CALL U2L1CR(CDORIGIN)
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
          ISTAT = ISETVALFDB( KFDB, 'origin', CDORIGIN )
        ENDIF

!       HINDCAST REFERENCE DATE:

        IF ( ISTREAM == 1204 ) THEN
          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'refdate',C12)
          CDATEREF=c12(1:8)
          ISTAT = ISETVALFDB(KFDB, 'refdate', CDATEREF)
        ELSE IF(ISTREAM == 1084 .OR.
     &          ISTREAM == 1085 ) THEN
          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'referenceDate',C12)
          CDATEREF=c12(1:8)
          ISTAT = ISETVALFDB(KFDB, 'refdate', CDATEREF)
        ELSE IF (ISTREAM == 1078 .OR.
     &           ISTREAM == 1079 ) THEN
          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'dataDate',C12)
          HDATE=C12(1:8)
          ISTAT = ISETVALFDB(KFDB,'hdate', HDATE)
        ENDIF

      ENDIF

!     ------------------------------------------------------------------

!*    2. MAKE DATA BASE NAME.

!     SET:

!     STREAM:
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'stream',C12)
      CLSTREAM=C12(1:4)
      ISTAT = ISETVALFDBSUBS( KFDB, 'stream', CLSTREAM )

!     CLASS:
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'class',C12)
      CLCLASS=C12(1:2)
      ISTAT = ISETVALFDBSUBS( KFDB, 'class',  CLCLASS )

!     EXPVER:
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'expver',C12)
      CLEXPVER=C12(1:4)
      ISTAT = ISETVALFDBSUBS( KFDB, 'expver', CLEXPVER )

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
      ISTAT = ISETVALFDBSUBS( KFDB, 'domain', CLDOMAIN )

!     REPRESENTATION:
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'gridType',C12)
      ISTAT = ISETVALFDBSUBS( KFDB, 'repr', c12(9:10) )

!     TYPE:
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'type',C12)
      CLTYPE=C12(1:2)
      ISTAT = ISETVALFDBSUBS( KFDB, 'type', CLTYPE )

!     PARAMETER:
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'param',C12)
      CLPARAM=C12(1:3)
      ISTAT = ISETVALFDBSUBS( KFDB, 'param', CLPARAM )

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
      ISTAT = ISETVALFDBSUBS( KFDB, 'date', CLDATE )

!     TIME:
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'time',C12)
      CLTIME=C12(1:4)
      ISTAT = ISETVALFDBSUBS( KFDB, 'time', CLTIME )

!     FC STEP:
      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'stepUnits','h')
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'endStep',ISTEP)
!     make sure the step is in hours
      ISTEP=ISTEP
      WRITE(CLSTEP, '( I6.6 )' ) ISTEP 
      ISTAT = ISETVALFDBSUBS( KFDB, 'step', CLSTEP )

!     LEVEL TYPE:
      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'levtype',ilvtp)
      IF(ILVTP == 209) THEN
        ISTAT = ISETVALFDBSUBS( KFDB, 'levty', 'wv'  )
      ELSE IF(ILVTP == 212) THEN
        ISTAT = ISETVALFDBSUBS( KFDB, 'levty', 'ws'  )
      ELSE
        ISTAT = ISETVALFDBSUBS( KFDB, 'levty', 's'   )
      ENDIF 

!     DIRECTION AND FREQUENCY:
      IF(CLPARAM == '251' ) THEN
         CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'directionNumber',C12)
         CLDIR=C12(1:2)
         ISTAT = ISETVALFDBSUBS( KFDB, 'direction', CLDIR )
         CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'frequencyNumber',C12)
         CLFRE=C12(1:2)
         ISTAT = ISETVALFDBSUBS( KFDB, 'frequency', CLFRE )
      ELSE
        CLDIR="00"
        CLFRE="00"
        IF (ISTREAM /= 1082 .AND.
     &      ISTREAM /= 1095 .AND.
     &      ISTREAM /= 1203 .AND.
     &      ISTREAM /= 1204 .AND.
     &      ISTREAM /= 1222) THEN
          ISTAT = ISETVALFDBSUBS( KFDB, 'direction', CLDIR )
          ISTAT = ISETVALFDBSUBS( KFDB, 'frequency', CLFRE )
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

!     ------------------------------------------------------------------

!*    3.  GET AND WRITE FIELD TO FDB

!     ------------------------------------------------------------------

      CALL IGRIB_GET_MESSAGE_SIZE(IGRIB_HANDLE,KBYTES)
      KLEN=(KBYTES+NPRECI)/NPRECI
      ALLOCATE(KGRIB(KLEN))
      CALL IGRIB_GET_MESSAGE(IGRIB_HANDLE,KGRIB)
      ISTAT = IWRITEFDBSUBS( KFDB, KGRIB, KLEN ) 
      DEALLOCATE(KGRIB)

      IF ( ISTAT .NE. 0 ) THEN
        WRITE(KUSO, '( "Error\ /WGRIB2FDB/ Failed to write to fdb")' )
        CALL GSTATS(1787,0)
        ISTAT = ICLOSEFDBSUBS( KFDB )
        CALL GSTATS(1787,1)
        KERR = 1
        RETURN
      ENDIF
      IF (KTEST.GT.1) WRITE(KUSO,'("\ /WGRIB2FDB/ iwritefdb klen:", i10,
     &                      " status:", i4)') KLEN, ISTAT

!     ------------------------------------------------------------------

!*    4.  RETURN

      IF (LHOOK) CALL DR_HOOK('WGRIB2FDB',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE WGRIB2FDB
