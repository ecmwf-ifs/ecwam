!=======================================================================

      SUBROUTINE FLD2FDB (kuso, ktest, ksec1, kfield, klen, cdfun, kerr,
     &                    cdtype, cdclass, cdexpver, cdfdbsf, cddomain,
     &                    kdir, kfre, kfdb)

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
!**             CDTYPE    - TYPE OF FIELD
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

      DIMENSION ksec1(*), kfield(*)

      INTEGER       klen, kdir, kfre, ndir, nfre 

      CHARACTER*(*) cdfun
      CHARACTER*(*) cdfdbsf
      CHARACTER*(*) cddomain

      CHARACTER*1  cldomain

      CHARACTER*2  cltime
      CHARACTER*2  cdtype
      CHARACTER*2  cltype
      CHARACTER*2  cdclass
      CHARACTER*2  clclass
      CHARACTER*2  cldir
      CHARACTER*2  clfre

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

!     HELPS TO CALL IOPENFDB ONLY ONCE.
!     ---------------------------------

      LOGICAL     lopenfdb 
      DATA        lopenfdb /.FALSE./
      SAVE        lopenfdb


!     KEEP THE VALUES OF the FDB REFERENCE

      INTEGER IADDR1
      SAVE IADDR1 

!     ------------------------------------------------------------------

!*        1.   DERIVE FIELD DATA BASE VARIABLES FROM ARGUMENT LIST.

!     ------------------------------------------------------------------

!     THIS VALUE WILL BE RETURNED IF EVERYTHING GO OK.
      kerr = 0

      CALL u2l1cr( cdfun )

      IF ( cdfun .NE. 'close' ) THEN
        IF ( ksec1(16) .NE. 0 .and. cdtype.ne.'4v' )  THEN
          if(ksec1(43).ne.0) then
            if(ksec1(42).ne.0) then
              cltype   = 'pf'
            else
              cltype   = 'cf'
            endif
          else
            cltype   = 'fc'
          endif
        ELSE
          cltype   = cdtype
        ENDIF
        if(cdtype.eq.'sf') then
          if(ksec1(43).ne.0) then
            if(ksec1(42).ne.0) then
              cltype   = 'fc'
            else
              cltype   = 'cf'
            endif
          else
            cltype   = 'fc'
          endif
        endif

        cldomain = cddomain

        IF ( cltype .EQ. 'fg' .and. ksec1(16) .EQ. 0 ) then
          clstep = '000006'
        ELSE
          WRITE(clstep, '( i6.6 )' ) ksec1(16)
        ENDIF

        WRITE( clparam, '(i3.3)' ) ksec1(6)
        iyyyy=(ksec1(21)-1)*100+ksec1(10)
        WRITE( cldate(01:04), '(i4.4)' ) iyyyy 
        WRITE( cldate(05:06), '(i2.2)' ) ksec1(11)
        WRITE( cldate(07:08), '(i2.2)' ) ksec1(12)
        WRITE( cltime, '(i2.2)' ) ksec1(13)

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
        IF ( .NOT. lopenfdb )  THEN
          if (cdtype.eq.'sf' .or. ksec1(40).eq.1082
     &                       .or. ksec1(40).eq.1095) then
            istat = iopenfdb( 'seas', iaddr1, cdfun(1:1) )
          else
            istat = iopenfdb( 'fdb', iaddr1, cdfun(1:1) )
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
            istat = iset_fdb_root(iaddr1,    cdfdbsf(1:icdfdbsf) )
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

          lopenfdb = .TRUE.

          call WSTREAM_STRG(ksec1(40),clstream)

          if(clstream.eq.'****') then 
            WRITE(kuso,'( "Error\ /Fld2fdb/ stream unknown ", i4)' )
     &      ksec1(40) 
            kerr = 1
            return
          endif

          istat = isetvalfdb( iaddr1,   'STREAM', clstream )

          istat = isetvalfdb( iaddr1,   'REPR',    'll' )

!         ensemble number (only set when needed)

          WRITE( clnumber, '(i3.3)' ) ksec1(42)

          if (cdtype.eq.'cf' .or. cdtype.eq.'pf'
     &                       .or. ksec1(40).eq.1081 .or.
     &        cdtype.eq.'sf' .or. ksec1(40).eq.1082
     &                       .or. ksec1(40).eq.1095) then
            istat = isetvalfdb( iaddr1,   'NUMBER',  clnumber)
          endif

          IF (KTEST.GT.3) WRITE(kuso,'("\ /Fld2fdb/ STREAM:", a5,
     &                    " NUMBER ", a3)') clstream, clnumber

!         system and method

          if (cdtype.eq.'sf' .or. ksec1(40).eq.1082
     &                       .or. ksec1(40).eq.1095) then
            if (kdir.NE.0 .AND. kfre.NE.0 .AND. 
     &          clparam .EQ. '251' ) THEN
              ndir=ksec1(44)
              nfre=ksec1(45)
              write( clsystem, '(i5.5)' ) ksec1(50+ndir+nfre)
              istat = isetvalfdb( iaddr1,   'system',  clsystem)
              write( clmethod, '(i5.5)' ) ksec1(51+ndir+nfre)
              istat = isetvalfdb( iaddr1,   'method',  clmethod)
            else
              write( clsystem, '(i5.5)' ) ksec1(44)
              istat = isetvalfdb( iaddr1,   'system',  clsystem)
              write( clmethod, '(i5.5)' ) ksec1(45)
              istat = isetvalfdb( iaddr1,   'method',  clmethod)
              istat = isetvalfdb( iaddr1, 'level','0000')
            endif

            IF (KTEST.GT.3) WRITE(kuso,'("\ /Fld2fdb/ system:", a5,
     &                      " method ", a5)') clsystem,clmethod 
          endif

!         multi-analysis ensemble runs (see setwgribhd)
          if (ksec1(40) == 1083) then
            if(ksec1(46) == 0 ) then
              if(ksec1(44) == 78) then
                cdorigin = 'edzw'
              else if(ksec1(44) == 85) then
                cdorigin = 'lfpw'
              else if(ksec1(44) == 7) then
                cdorigin = 'kwbc'
              else if(ksec1(44) == 74) then
                cdorigin = 'egrr'
              else
                WRITE(kuso,'( "Error\ /Fld2fdb/ origin unknown ", i4)' )
     &                ksec1(44) 
                kerr = 1
                return
              endif
            else
              cdorigin = 'cons'
            endif
            istat = isetvalfdb( iaddr1,  'origin', cdorigin )
          endif

        ENDIF
      ELSEIF ( cdfun .EQ. 'close' ) THEN
        istat = iclosefdb( iaddr1 )
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
          lopenfdb = .FALSE.
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

      istat = isetvalfdb( iaddr1, 'class',  clclass  )
      istat = isetvalfdb( iaddr1, 'expver', clexpver )
      istat = isetvalfdb( iaddr1, 'domain', cldomain )
      istat = isetvalfdb( iaddr1, 'type',   cltype   )
      istat = isetvalfdb( iaddr1, 'param',  clparam  )
      istat = isetvalfdb( iaddr1, 'date',   cldate   )
      istat = isetvalfdb( iaddr1, 'time',   cltime   )
      istat = isetvalfdb( iaddr1, 'step',   clstep   )
      istat = isetvalfdb( iaddr1, 'levty',   's'        )
!     istat = isetvalfdb( iaddr1, 'level',   '0000'     )
!     istat = isetvalfdb( iaddr1, 'fstep',  '000000'    )

      IF (kdir.NE.0 .AND. kfre.NE.0 .AND. clparam .EQ. '251' ) THEN
        idir = kdir 
        ifre = kfre 
        WRITE( cldir, '(i2.2)' ) idir 
        WRITE( clfre, '(i2.2)' ) ifre 
        istat = isetvalfdb( iaddr1,   'direction', cldir )
        istat = isetvalfdb( iaddr1,   'frequency', clfre )
      ELSE
        if (cdtype.ne.'sf' .and. ksec1(40).ne.1082
     &                     .and. ksec1(40).ne.1095) then
          cldir="00"
          clfre="00"
          istat = isetvalfdb( iaddr1,   'direction', cldir )
          istat = isetvalfdb( iaddr1,   'frequency', clfre )
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
        istat = iwritefdb( iaddr1, kfield, klen ) 
        IF ( istat .NE. 0 ) THEN
          WRITE(kuso, '( "Error\ /Fld2fdb/ Failed to write to fdb")' )
          istat = iclosefdb( iaddr1 )
          kerr = 1
          RETURN
        ENDIF
        IF (KTEST.GT.1) WRITE(kuso,'("\ /Fld2fdb/ iwritefdb klen:", i10,
     &                  " status:", i4)') klen, istat
      ELSEIF ( cdfun .EQ. 'read' )  THEN
        istat = ireadfdb( iaddr1, kfield, klen )
        IF ( istat .NE. 0 ) THEN
          WRITE(kuso, '( "Error\ /Fld2fdb/ Failed to read from fdb")' )
          istat = iclosefdb( iaddr1 )
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

      kfdb = iaddr1
      RETURN
      END SUBROUTINE FLD2FDB


