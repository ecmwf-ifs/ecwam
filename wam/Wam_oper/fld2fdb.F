!=======================================================================

      SUBROUTINE FLD2FDB (kuso, ktest, ksec1, kfield, klen, cdfun, kerr,
     &                    cdclass, cdexpver, cdfdbsf, cddomain,
     &                    kdir, kfre, kfdb, lfdbopen)

!=======================================================================

!**
!**** NAME    *FLD2FDB*
!**** ----
!**
!**   PURPOSE
!**   -------
!**      WRITE/READ FIELD TO/FROM THE FIELD DATA BASE.
!**
!**   INTERFACE
!**   ---------
!**      INPUT: KUSO      - LOGICAL UNIT FOR STANDARD OUTPUT.
!**             KTEST     - SWITCH DIAGNOSTICS OUTPUT ON IF KTEST GT 1.
!**             KSEC1     - GRIB SECTION 1 OF THE KFIELD.
!**             CDFUN     - FUNCTION: 'WRITE' - WRITE KFIELD INTO FDB.
!**                                   'READ'  - READ FDB FIELD.
!**                                   'CLOSE' - CLOSE FDB.
!**             CDCLASS   - CLASS OF CURRENT RUN.
!**             CDEXPVER  - EXPERIMENT VERSION NUMBER.
!**             CDDOMAIN  - FDB DOMAIN IDENTIFIER.
!**             KDIR      - CURRENT DIRECTION IF PARAM=251, ELSE  = 0 
!**             KFRE      - CURRENT FREQUENCY IF PARAM=251, ELSE = 0 
!**
!**      INPUT FOR THE CDFUN = 'WRITE':
!**
!**              KFIELD   - INTEGER ARRAY OF GRIB CODE.
!**              KLEN     - LEN OF THE ARRAY KFIELD.
!**
!**
!**      OUTPUT FOR THE CDFUN = 'READ':
!**
!**              KFIELD   - INTEGER ARRAY OF GRIB CODE.
!**              KLEN     - LEN OF THE ARRAY KFIELD.
!**
!**      OUTPUT FOR ALL FUNCTIONS: KERR - 0; NO ERRORS ENCOUNTERED
!**                                KFDB - DATA BASE REFERENCE
!**
!**   METHOD
!**   ------
!**      FIELD DATA BASE ROUTINES ARE USED TO OPEN THE FIELD
!**      DATA BASE AND FILE, AND TO WRITE/READ THE FIELD INTO/
!**      FROM THE FIELD DATA BASE.
!**
!**              DATABASE NAME:'ODXXTYYYYYMMDDTT', WHERE
!**                             TY     - DATA TYPE ('AN','FC')
!**                             YYYYMMDD - DATE (YEAR,MONTH,DAY)
!**                             TT     - BASE TIME
!**
!**              FILE NAMES:   'REXXX0TTT00', WHERE
!**                             RE     - DATA REPRESENTATION 'LL'
!**                             TTTTTT - TIME STEP
!**
!**              FIELD NAMES:  'PNNNLLLLHHHHXX', WHERE
!**                             NNN  - PARAMTER NUMBER
!**                             LLLL - LEVEL TYPE 'SFCX'
!**                             HHHH - LEVEL HIGHT '0000'
!**
!**   EXTERNALS  $FDBLIB
!**   ---------
!**      * IOPENFDB   *  -  OPEN       DATA BASE
!**      * ICLOSEFDB  *  -  CLOSE      DATA BASE
!**      * ISETVALFDB *  -  SET PARAMETERS ( CLASS, TYPE, DATE ETC )
!**      * IREADFDB   *  -  READ FIELD FROM DATA BASE
!**      * IWRITEFDB  *  -  WRITE FIELD INTO DATA BASE
!**
!**
!**   AUTHOR
!**   ------
!**      D. JOKIC,  ECMWF,  JUN-1992.
!**
!**   MODIFICATIONS.
!**   --------------
!**      B. HANSEN,  ECMWF,  AUG-1994.
!**      CREATES NOW ALSO FDB DATA BASES FOR RD EXPERIMENTS
!**
!**      B. HANSEN,  ECMWF,  MAR-1995.
!**      PASSES A REFERENCE TO THE DATABASE BACK TO THE CALLER.
!**
!**      G. KONSTANDINIDIS, ECMWF, MAY-1996
!**        MAJOR REWRITE WHILE KEEPING THE SAME INTERFACE AND SIMILAR
!**        FUNCTIONALITY.
!**
!**      B. HANSEN,  ECMWF, NOV-1996
!**        A) TUNING (DO FUNCTION CALLS ONLY WHEN NEEDED).
!**        B) CATER FOR TWO NEW DOMAINS 'N' FOR NORTHER HEMISPHERE
!**                                     's' for Southern Hemisphere
!        J. BIDLOT ECMWF NOV 1997 : allow new stream waef and type
!                                   cf or pf from ensemble runs.
!        J. BIDLOT ECMWF MAY 2000 : allow new stream wasf and 
!                                   system and method numbers.

!        B. HANSEN,  ECMWF, FEB-1998
!          MOVE INITIALIZATION TO USERIN AT THE BOTTOM OF SECTION 1.4. 

!     ------------------------------------------------------------------

!*        1. VARIABLE DECLARATIONS

!     ------------------------------------------------------------------

!     the module is needed in order to resolve the problem with the 
!     output of 251 for the multi-analysis problem.
!     because 251 does not contain the origin in the grib header
!     it is not very clean !!!!!
      USE YOWGRIBHD , ONLY : KSEK1
      USE FDBSUBS

      INTEGER       klen
      DIMENSION ksec1(*), kfield(klen)

      INTEGER       kdir, kfre, ndir, nfre 

      CHARACTER*(*) cdfun
      CHARACTER*(*) cdfdbsf
      CHARACTER*(*) cddomain

      CHARACTER*1  cldomain
      CHARACTER*1  cleg

      CHARACTER*2  cltime
      CHARACTER*2  cltype
      CHARACTER*2  cdclass
      CHARACTER*2  clclass
      CHARACTER*2  cldir
      CHARACTER*2  clfre
      CHARACTER*2  cdum

      CHARACTER*3  clparam
      CHARACTER*3  clnumber

      CHARACTER*4  cdexpver
      CHARACTER*4  clexpver
      CHARACTER*5  clsystem
      CHARACTER*5  clmethod
      CHARACTER*6  clstep
      CHARACTER*4  clstream
      CHARACTER*4  cdorigin

      CHARACTER*8  cldate
      CHARACTER*8  hdate
      CHARACTER*8  cdateref

      LOGICAL :: LASTREAM 
      LOGICAL :: lfdbopen 

!     ------------------------------------------------------------------

!*        1.   DERIVE FIELD DATA BASE VARIABLES FROM ARGUMENT LIST.

!     ------------------------------------------------------------------

!     THIS VALUE WILL BE RETURNED IF EVERYTHING GO OK.
      kerr = 0

      CALL u2l1cr( cdfun )

      IF ( cdfun .NE. 'close' ) THEN
        idum=0
        call WSTREAM_STRG(ksec1(40),clstream,idum,idum,cdum,
     &                    idum,lastream)
        if(clstream.eq.'****') then 
          WRITE(kuso,'( "Error\ /Fld2fdb/ stream unknown ", i4)' )
     &    ksec1(40) 
          kerr = 1
          return
        endif
        IF (ksec1(39) .EQ. 9 ) THEN
          cltype='fc'
        ELSE IF  (ksec1(39) .EQ. 10 ) THEN
          cltype='cf'
        ELSE IF  (ksec1(39) .EQ. 11 ) THEN
          cltype='pf'
        ELSE IF  (ksec1(39) .EQ. 6 ) THEN
          cltype='4v'
        ELSE IF  (ksec1(39) .EQ. 2 ) THEN
          cltype='an'
        ELSE IF  (ksec1(39) .EQ. 65 ) THEN
          cltype='cv'
        ELSE IF  (ksec1(39) .EQ. 1 ) THEN
          cltype='fg'
        ELSE
          WRITE(kuso,'("Error\ /Fld2fdb/type in grib header not known",
     &    I6)') ksec1(39) 
        ENDIF

        cldomain = cddomain

        IF ( cltype .EQ. 'fg' .and. ksec1(16) .EQ. 0 ) then
          clstep = '000006'
        ELSE
!         make sure the step is in hours
          if(ksec1(15).eq.1) then
            istep=ksec1(16)
          else if(ksec1(15).eq.2) then
            istep=ksec1(16)*24
          else if(ksec1(15).eq.12) then
            istep=ksec1(16)*12
          else if(ksec1(15).eq.11) then
            istep=ksec1(16)*6
          else if(ksec1(15).eq.10) then
            istep=ksec1(16)*3
          else if(ksec1(15).eq.0) then
            if(mod(ksec1(16),60).ne.0) then
              WRITE(kuso,'("Error\ /Fld2fdb/step not multiple of hours",
     &        I6)') ksec1(16) 
              kerr = 1
              RETURN
            endif
            istep=ksec1(16)/60
          else
            WRITE(kuso,'("Error\ /Fld2fdb/doesnot know how to use step",
     &      I3)') ksec1(15) 
            kerr = 1
            RETURN
          endif
          WRITE(clstep, '( i6.6 )' ) istep 
        ENDIF

        WRITE( clparam, '(i3.3)' ) ksec1(6)
        iyyyy=(ksec1(21)-1)*100+ksec1(10)
        WRITE( cldate(01:04), '(i4.4)' ) iyyyy 
        WRITE( cldate(05:06), '(i2.2)' ) ksec1(11)
        WRITE( cldate(07:08), '(i2.2)' ) ksec1(12)
        WRITE( cltime, '(i2.2)' ) ksec1(13)

!
! For VAREPS hindcast we use a different logic 
!

        hdate=cldate
        if (ksec1(37) == 27 .or. ksec1(37) == 30) then
          NREFERENCE=ksec1(49)
          if(NREFERENCE.gt.0) then
              write( cldate, '(i8.8)' )  NREFERENCE
          endif
        elseif(ksec1(37) == 13 ) then
            if (ksec1(40) == 1078 .or.
     &          ksec1(40) == 1079 ) then
                   ndir=ksec1(46)
                   nfre=ksec1(47)
                   NREFERENCE=ksec1(52+ndir+nfre)
                   if(NREFERENCE.gt.0) then
                      write( cldate, '(i8.8)' )  NREFERENCE
                   endif
            endif
        endif

        clclass  = cdclass
        clexpver = cdexpver

        CALL u2l1cr( cltype )
        CALL u2l1cr( clclass )
        CALL u2l1cr( clexpver )
      ENDIF


!     ------------------------------------------------------------------

!*        2.   OPEN - CLOSE - CHECK ACTION REQUESTED

!     ------------------------------------------------------------------

      IF ( cdfun .EQ. 'write' .OR. cdfun .EQ. 'read' ) THEN
        IF ( .NOT. lfdbopen )  THEN
          if (ksec1(40).eq.1082 .or.
     &        ksec1(40).eq.1222 .or.
     &        ksec1(40).eq.1095 .or.
     &        ksec1(40).eq.1203 .or.
     &        ksec1(40).eq.1204) then
            istat = iopenfdbsubs( 'seas', kfdb, cdfun(1:1) )
          else
            istat = iopenfdbsubs( 'fdb', kfdb, cdfun(1:1) )
          endif
          IF ( istat .NE. 0 ) THEN
            WRITE(kuso,'("Error\ /Fld2fdb/ Action specified: ",a)')
     &      cdfun
            WRITE(kuso,'("Error\ /Fld2fdb/ iopenfdb return status:",
     &      I3)') istat
            kerr = 1
            RETURN
          ELSEIF (KTEST .GT. 1) THEN
            WRITE(kuso,'("\ /Fld2fdb/ iopenfdb status:", i3)') istat
          ENDIF

          icdfdbsf = LEN_TRIM(cdfdbsf)
          IF (icdfdbsf .GT. 0 ) THEN
            istat = iset_fdbsubs_root(kfdb,    cdfdbsf(1:icdfdbsf) )
            IF ( istat .NE. 0 ) THEN
              WRITE(kuso,'("Error\ /Fld2fdb/ Root dir specified: ",a)')
     &        cdfdbsf(1:icdfdbsf)
              WRITE(kuso,'("Error\ /Fld2fdb/ iopenfdb return status:",
     &        I3)') istat
              kerr = 1
              RETURN
            ELSEIF (KTEST .GT. 1) THEN
              WRITE(kuso,'("\ /Fld2fdb/ iset_fdb_root status:", i3,
     &        " Root: ", a25)') istat, cdfdbsf(1:icdfdbsf)
            ENDIF
          ENDIF

          lfdbopen = .TRUE.

          idum=0
          call WSTREAM_STRG(ksec1(40),clstream,idum,idum,cdum,
     &                    idum,lastream)

          if(clstream.eq.'****') then 
            WRITE(kuso,'( "Error\ /Fld2fdb/ stream unknown ", i4)' )
     &      ksec1(40) 
            kerr = 1
            return
          endif

          istat = isetvalfdbsubs( kfdb,   'STREAM', clstream )

          istat = isetvalfdbsubs( kfdb,   'REPR',    'll' )

!         ensemble number and leg (only set when needed)
          WRITE( clnumber, '(i3.3)' ) ksec1(42)

          if ( ksec1(40).eq.1081 .OR.
     &         ksec1(40).eq.1083 .OR.
     &         ksec1(40).eq.1084 .OR.
     &         ksec1(40).eq.1078 .OR.
     &         ksec1(40).eq.1079 .OR.
     &         ksec1(40).eq.1086 ) then
            if ( ksec1(42).gt.0 ) then
              istat = isetvalfdbsubs( kfdb,   'NUMBER',  clnumber)
            endif
            if (kdir.NE.0 .AND. kfre.NE.0 .AND. 
     &         (clparam .EQ. '251' .OR. ksec1(7).eq.212) ) THEN
              ndir=ksec1(46)
              nfre=ksec1(47)
              if( ksec1(57+ndir+nfre).gt.0 ) then
                WRITE( cleg, '(i1.1)' ) ksec1(57+ndir+nfre)
                istat = isetvalfdbsubs( kfdb,   'LEG',  cleg)
              endif
            else
              if( ksec1(48).gt.0 ) then
                WRITE( cleg, '(i1.1)' ) ksec1(48)
                istat = isetvalfdbsubs( kfdb,   'LEG',  cleg)
              endif
            endif
          else if(ksec1(40).eq.1203 .or. ksec1(40).eq.1204) then
            if ( ksec1(42).gt.0 .or. ksec1(39).eq.2) then
              istat = isetvalfdb( kfdb,   'NUMBER',  clnumber)
            endif
          else if(ksec1(40).eq.1082 .OR.
     &            ksec1(40).eq.1095 .OR.
     &            ksec1(40).eq.1222) then
              istat = isetvalfdbsubs( kfdb,   'NUMBER',  clnumber)
          endif

!         system and method

          if (ksec1(40).eq.1082 .or.
     &        ksec1(40).eq.1095 .or.
     &        ksec1(40).eq.1222 .or.
     &        ksec1(40).eq.1203 .or.
     &        ksec1(40).eq.1204) then
            if (kdir.NE.0 .AND. kfre.NE.0 .AND. 
     &         (clparam .EQ. '251' .OR. ksec1(7).eq.212) ) THEN
              ndir=ksec1(46)
              nfre=ksec1(47)

!!            if the class is not od then system is set to off
              if(ksec1(38).ne.1 ) then
                clsystem='off'
              else
                write( clsystem, '(i5.5)' ) ksec1(50+ndir+nfre)
              endif
              istat = isetvalfdbsubs( kfdb,   'system',  clsystem)

              write( clmethod, '(i5.5)' ) ksec1(51+ndir+nfre)
              istat = isetvalfdbsubs( kfdb,   'method',  clmethod)
            else
!!            if the class is not od then system is set to off
              if(ksec1(38).ne.1 ) then
                clsystem='off'
              else
                write( clsystem, '(i5.5)' ) ksec1(44)
              endif
              istat = isetvalfdbsubs( kfdb,   'system',  clsystem)

              write( clmethod, '(i5.5)' ) ksec1(45)
              istat = isetvalfdbsubs( kfdb,   'method',  clmethod)
              istat = isetvalfdbsubs( kfdb, 'level','0000')
            endif

            IF (KTEST.GT.3) WRITE(kuso,'("\ /Fld2fdb/ system:", a5,
     &                      " method ", a5)') clsystem,clmethod 
          endif

!         centre origin 
          if (ksec1(37) == 18) then
            if(KSEK1(46) == 0 ) then
              if(KSEK1(44) == 78) then
                cdorigin = 'edzw'
              else if(KSEK1(44) == 85) then
                cdorigin = 'lfpw'
              else if(KSEK1(44) == 7) then
                cdorigin = 'kwbc'
              else if(KSEK1(44) == 74) then
                cdorigin = 'egrr'
              else
                WRITE(kuso,'( "Error\ /Fld2fdb/ origin unknown ", i4)' )
     &                KSEK1(44)
                kerr = 1
                return
              endif
              istat = isetvalfdbsubs( kfdb,  'origin', cdorigin )
            else
              istat = isetvalfdbsubs( kfdb,  'origin', 'consensus' )
            endif
          else if(ksec1(37) == 23) then
            if(ksec1(2).eq.98) then
              cdorigin = 'ecmf'
            else
              WRITE(kuso,'( "Error\ /Fld2fdb/ origin unknown ", i4)' )
     &              ksec1(2)
              kerr = 1
              return
            endif
            istat = isetvalfdb( kfdb,  'origin', cdorigin )
          endif

!         hindcast reference date
          if(ksec1(37) == 23 ) then
            if(ksec1(40) == 1204) then
              write( cdateref, '(i8.8)' ) ksec1(49) 
              istat = isetvalfdb( kfdb, 'refdate', cdateref)
            endif
          elseif(ksec1(37) == 26 ) then
            write( cdateref, '(i8.8)' ) ksec1(44) 
            istat = isetvalfdb( kfdb, 'refdate', cdateref)
          elseif(ksec1(37) == 27 .or. ksec1(37) == 30 ) then
            if(ksec1(40) == 1078.or.ksec1(40) == 1079) then
              istat = isetvalfdb( kfdb, 'hdate', hdate)
            endif
          elseif(ksec1(37) == 26 ) then
            write( cdateref, '(i8.8)' ) ksec1(44) 
            istat = isetvalfdb( kfdb, 'refdate', cdateref)
          elseif(ksec1(37) == 13 ) then
            ndir=ksec1(46)
            nfre=ksec1(47)
            if(ksec1(40) == 1084 .or.
     &         ksec1(40) == 1085 .or.
     &         ksec1(40) == 1204 ) then
              write( cdateref, '(i8.8)' ) ksec1(52+ndir+nfre)
              istat = isetvalfdb(kfdb,'refdate',cdateref)
            elseif (ksec1(40) == 1078 .or.
     &              ksec1(40) == 1079 ) then
               istat = isetvalfdb(kfdb,'hdate', hdate)
            endif
          endif

        ENDIF
      ELSEIF ( cdfun .EQ. 'close' ) THEN
        CALL GSTATS(1787,0)
        istat = iclosefdbsubs( kfdb )
        CALL GSTATS(1787,1)
        IF ( istat .NE. 0 ) THEN
          WRITE(kuso,'( "Error\ /Fld2fdb/ Action specified: ", a)' )
     &    cdfun
          WRITE(kuso, '( "Error\ /Fld2fdb/ iclosefdb return status:",
     &    I3)') istat
          kerr = 1
          RETURN
        ELSE
          IF (KTEST.GT.1)
     &     WRITE(kuso, '( "/Fld2fdb/ Closed FDB and DB: " )' )
          kerr = 0
          lfdbopen = .FALSE.
          RETURN
        ENDIF
      ELSE
        WRITE(kuso, '( "Error\ /Fld2fdb/ Unrecognised action",
     &            " specified: ", a )' ) cdfun
        kerr = 1
        RETURN
      ENDIF

!     ------------------------------------------------------------------

!*    3. MAKE DATA BASE NAME.

!     SET : CLASS - STREAM - EXPVER - DOMAIN - TYPE - PARAMETER
!           LEVELTYPE - LEVEL - REPRESENTATION 
!           DATE - TIME - STEP

!     ------------------------------------------------------------------
      istat = isetvalfdbsubs( kfdb, 'stream', clstream )
      istat = isetvalfdbsubs( kfdb, 'class',  clclass  )
      istat = isetvalfdbsubs( kfdb, 'expver', clexpver )
      istat = isetvalfdbsubs( kfdb, 'domain', cldomain )
      istat = isetvalfdbsubs( kfdb, 'repr',    'll' )
      istat = isetvalfdbsubs( kfdb, 'type',   cltype   )
      istat = isetvalfdbsubs( kfdb, 'param',  clparam  )
      istat = isetvalfdbsubs( kfdb, 'date',   cldate   )
      istat = isetvalfdbsubs( kfdb, 'time',   cltime   )
      istat = isetvalfdbsubs( kfdb, 'step',   clstep   )
      if(ksec1(7).eq.209) then
        istat = isetvalfdbsubs( kfdb, 'levty',   'wv'     )
      else if(ksec1(7).eq.212) then
        istat = isetvalfdbsubs( kfdb, 'levty',   'ws'     )
      else
        istat = isetvalfdbsubs( kfdb, 'levty',   's'     )
      endif 

      IF (kdir.NE.0 .AND. kfre.NE.0 .AND. 
     &    (clparam .EQ. '251' .OR. ksec1(7).eq.212) ) THEN
        idir = kdir 
        ifre = kfre 
        WRITE( cldir, '(i2.2)' ) idir 
        WRITE( clfre, '(i2.2)' ) ifre 
        istat = isetvalfdbsubs( kfdb,   'direction', cldir )
        istat = isetvalfdbsubs( kfdb,   'frequency', clfre )
      ELSE
        if (ksec1(40).ne.1082 .and.
     &      ksec1(40).ne.1095 .and.
     &      ksec1(40).ne.1203 .and.
     &      ksec1(40).ne.1204 .and.
     &      ksec1(40).ne.1222) then
          cldir="00"
          clfre="00"
          istat = isetvalfdbsubs( kfdb,   'direction', cldir )
          istat = isetvalfdbsubs( kfdb,   'frequency', clfre )
        endif
      ENDIF

      IF (KTEST.GT.3) THEN
        WRITE(kuso,'("\ /Fld2fdb/ CLASS: ", a8)') clclass
        WRITE(kuso,'("\ /Fld2fdb/ EXPVER:", a8)') clexpver
        WRITE(kuso,'("\ /Fld2fdb/ DOMAIN:", a8)') cldomain
        WRITE(kuso,'("\ /Fld2fdb/ TY:    ", a8)') cltype
        WRITE(kuso,'("\ /Fld2fdb/ PARAM: ", a8)') clparam
        WRITE(kuso,'("\ /Fld2fdb/ DATE:  ", a8)') cldate
        WRITE(kuso,'("\ /Fld2fdb/ TIME:  ", a8)') cltime
        WRITE(kuso,'("\ /Fld2fdb/ STEP:  ", a8)') clstep
        WRITE(kuso,'("\ /Fld2fdb/ DIR:   ", a8)') cldir
        WRITE(kuso,'("\ /Fld2fdb/ FREQ:  ", a8)') clfre
        CALL flush (kuso)
      ENDIF

!     ------------------------------------------------------------------

!*    4.  READ OR WRITE FIELD FROM/TO FDB

!     ------------------------------------------------------------------

      IF ( cdfun .EQ. 'write' )      THEN
        istat = iwritefdbsubs( kfdb, kfield, klen ) 
        IF ( istat .NE. 0 ) THEN
          WRITE(kuso, '( "Error\ /Fld2fdb/ Failed to write to fdb")' )
          CALL GSTATS(1787,0)
          istat = iclosefdbsubs( kfdb )
          CALL GSTATS(1787,1)
          kerr = 1
          RETURN
        ENDIF
        IF (KTEST.GT.1) WRITE(kuso,'("\ /Fld2fdb/ iwritefdb klen:", i10,
     &                  " status:", i4)') klen, istat
      ELSEIF ( cdfun .EQ. 'read' )  THEN
        istat = ireadfdbsubs( kfdb, kfield, klen )
        IF ( istat .NE. 0 ) THEN
          WRITE(kuso, '( "Error\ /Fld2fdb/ Failed to read from fdb")' )
          CALL GSTATS(1787,0)
          istat = iclosefdbsubs( kfdb )
          CALL GSTATS(1787,1)
          kerr = 1
          RETURN
        ENDIF
      ELSE
!       NORMALLY YOU SHOULD NOT ARRIVE AT THIS POINT
        WRITE(kuso, '( "Error\ /Fld2fdb/ Unrecognised action",
     &  " specified: ", a )' ) cdfun
        kerr = 1
        RETURN
      ENDIF

!     ------------------------------------------------------------------

!*    5.  RETURN

!     ------------------------------------------------------------------

      RETURN
      END SUBROUTINE FLD2FDB


