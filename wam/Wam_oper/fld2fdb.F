C=======================================================================
C
      SUBROUTINE FLD2FDB (kuso, ksec1, ksec2, kfield, klen, cdfun, kerr,
     .                    cdtype, cdclass, cdexpver, cdfdbsf, kfdb)
C
C=======================================================================
C
C**
C**** NAME    *FLD2FDB*
C**** ----
C**
C**   PURPOSE
C**   -------
C**      WRITE/READ FIELD TO/FROM THE FIELDS DATA BASE.
C**
C**   INTERFACE
C**   ---------
C**      INPUT: KUSO      - LOGICAL UNIT FOR STANDARD OUTPUT.
C**             KSEC1     - GRIB SECTION 1 OF THE KFIELD.
C**             KSEC2     - GRIB SECTION 2 OF THE KFIELD.
C**             CDFUN     - FUNCTION: 'WRITE' - WRITE KFIELD INTO FDB.
C**                                   'READ'  - READ FDB FIELD.
C**                                   'CLOSE' - CLOSE FDB.
C**             CDTYPE    - TYPE OF FIELD ( FG, AN OR FC).
C**             CDCLASS   - CLASS OF CURRENT RUN.
C**             CDEXPVER  - EXPERIMENT VERSION NUMBER.
C**
C**      INPUT FOR THE CDFUN = 'WRITE':
C**
C**              KFIELD   - INTEGER ARRAY OF GRIB CODE.
C**              KLEN     - LEN OF THE ARRAY KFIELD.
C**
C**
C**      OUTPUT FOR THE CDFUN = 'READ':
C**
C**              KFIELD   - INTEGER ARRAY OF GRIB CODE.
C**              KLEN     - LEN OF THE ARRAY KFIELD.
C**
C**      OUTPUT FOR ALL FUNCTIONS: KERR - 0; NO ERRORS ENCOUNTERED
C**                                KFDB - DATA BASE REFERENCE
C**
C**   METHOD
C**   ------
C**      FIELDS DATA BASE ROUTINES ARE USED TO OPEN THE FIELDS
C**      DATA BASE AND FILE, AND TO WRITE/READ THE FIELD INTO/
C**      FROM THE FIELDS DATA BASE.
C**
C**              DATABASE NAME:'ODXXTYYYMMDDTT', WHERE
C**                             TY     - DATA TYPE ('AN','FC')
C**                             YYMMDD - DATE (YEAR,MONTH,DAY)
C**                             TT     - BASE TIME
C**
C**              FILE NAMES:   'REXXX0TTT00', WHERE
C**                             RE  - DATA REPRESENTATION 'LL'
C**                             TTT - TIME STEP
C**
C**              FIELD NAMES:  'PNNNLLLLHHHHXX', WHERE
C**                             NNN  - PARAMTER NUMBER
C**                             LLLL - LEVEL TYPE 'SFCX'
C**                             HHHH - LEVEL HIGHT '0000'
C**
C**   EXTERNALS  $FDBLIB
C**   ---------
C**      * IINITFDB   *  -  INITIALISE DATA BASE
C**      * IOPENFDB   *  -  OPEN       DATA BASE
C**      * ICLOSEFDB  *  -  CLOSE      DATA BASE
C**      * ISETVALFDB *  -  SET PARAMETERS ( CLASS, TYPE, DATE ETC )
C**      * IREADFDB   *  -  READ FIELD FROM DATA BASE
C**      * IWRITEFDB  *  -  WRITE FIELD INTO DATA BASE
C**
C**
C**   AUTHOR
C**   ------
C**      D. JOKIC,  ECMWF,  JUN-1992.
C**
C**   MODIFICATIONS.
C**   --------------
C**      B. HANSEN,  ECMWF,  AUG-1994.
C**      CREATES NOW ALSO FDB DATA BASES FOR RD EXPERIMENTS
C**
C**      B. HANSEN,  ECMWF,  MAR-1995.
C**      PASSES A REFERENCE TO THE DATABASE BACK TO THE CALLER.
C**
C**      G. KONSTANDINIDIS, ECMWF, MAY-1996
C**        MAJOR REWRITE WHILE KEEPING THE SAME INTERFACE AND SIMILAR
C**        FUNCTIONALITY.
C**
C**      B. HANSEN,  ECMWF, NOV-1996
C**        A) TUNING (DO FUNCTION CALLS ONLY WHEN NEEDED).
C**        B) CATER FOR TWO NEW DOMAINS 'N' FOR NORTHER HEMISPHERE
C**                                     's' for Southern Hemisphere
C        J. BIDLOT ECMWF NOV 1997 : allow new stream waef and type cf or pf
C                                   from ensemble runs.
C
C     ------------------------------------------------------------------
C
C*        1. VARIABLE DECLARATIONS
C
C     ------------------------------------------------------------------
C
      DIMENSION ksec1(*), ksec2(*), kfield(*)
C
      INTEGER       klen
C
      CHARACTER*(*) cdfun
      CHARACTER*(*) cdfdbsf
C
      CHARACTER*1  cldomain
C
      CHARACTER*2  cltime
      CHARACTER*2  cdtype
      CHARACTER*2  cltype
      CHARACTER*2  cdclass
      CHARACTER*2  clclass
C
      CHARACTER*3  clparam
      CHARACTER*3  clnumber
C
      CHARACTER*4  cdexpver
      CHARACTER*4  clexpver
      CHARACTER*4  clstep
C
      CHARACTER*6  cldate
C
      CHARACTER*1  cl_p_domain
      CHARACTER*2  cl_p_time
      CHARACTER*2  cl_p_type
      CHARACTER*2  cl_p_class
      CHARACTER*3  cl_p_param
      CHARACTER*4  cl_p_expver
      CHARACTER*4  cl_p_step
      CHARACTER*6  cl_p_date
C
      DATA cl_p_domain, cl_p_time, cl_p_type, cl_p_class, cl_p_param,
     .     cl_p_expver, cl_p_step, cl_p_date / 8*"  " /
C
C     HELPS TO CALL IINITFDB ONLY ONCE.
C     ---------------------------------
C
      LOGICAL     linitfdb 
      DATA        linitfdb /.FALSE./
      SAVE        linitfdb
C
C     HELPS TO CALL IOPENFDB ONLY ONCE.
C     ---------------------------------
C
      LOGICAL     lopenfdb 
      DATA        lopenfdb /.FALSE./
      SAVE        lopenfdb
C
C
C     KEEP THE VALUES OF the FDB REFERENCE
C
      INTEGER IADDR1
      SAVE IADDR1 
C
C     ------------------------------------------------------------------
C
C*        1.   DERIVE FIELDS DATA BASE VARIABLES FROM ARGUMENT LIST.
C
C     ------------------------------------------------------------------
C
C     THIS VALUE WILL BE RETURNED IF EVERYTHING GO OK.
      kerr = 0
C
      CALL u2l1cr( cdfun )
C
      IF ( cdfun .NE. 'close' ) THEN
        IF ( ksec1(16) .NE. 0 )  THEN
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
C
        IF ( ksec1(3).ge.204 ) THEN
          cldomain = 'm'
        ELSE
          IF ( ksec2(7) .GE. 0 ) THEN
            cldomain = 'n'
          ELSEIF ( ksec2(4) .LT. 0 ) THEN
            cldomain = 's'
          ELSE
            cldomain = 'g'
          ENDIF
        ENDIF
C
        IF ( cltype .EQ. 'fg' .and. ksec1(16) .EQ. 0 ) then
          clstep = '0006'
        ELSE
          WRITE(clstep, '( i4.4 )' ) ksec1(16)
        ENDIF
C
        WRITE( clparam, '(i3.3)' ) ksec1(6)
        WRITE( cldate(01:02), '(i2.2)' ) ksec1(10)
        WRITE( cldate(03:04), '(i2.2)' ) ksec1(11)
        WRITE( cldate(05:06), '(i2.2)' ) ksec1(12)
        WRITE( cltime, '(i2.2)' ) ksec1(13)
C
        clclass  = cdclass
        clexpver = cdexpver
C
        CALL u2l1cr( cltype )
        CALL u2l1cr( clclass )
        CALL u2l1cr( clexpver )
      ENDIF
C
C
C     ------------------------------------------------------------------
C
C*        2.   INITIALISE - OPEN - CLOSE - CHECK ACTION REQUESTED
C
C     ------------------------------------------------------------------
C
      IF ( .NOT. linitfdb )  THEN
        istat = iinitfdb()
        Linitfdb = .TRUE.
      ENDIF
C
      IF ( cdfun .EQ. 'write' .OR. cdfun .EQ. 'read' ) THEN
        IF ( .NOT. lopenfdb )  THEN
          istat = iopenfdb( 'fdb', iaddr1, cdfun(1:1) )
          IF ( istat .NE. 0 ) THEN
            WRITE(kuso,'("Error\ /Fld2fdb/ Action specified: ",a)')
     .      cdfun
            WRITE(kuso,'("Error\ /Fld2fdb/ iopenfdb return status:",
     .      I3)') istat
            kerr = 1
            RETURN
          ENDIF
C
          icdfdbsf = LEN_TRIM(cdfdbsf)
          istat = iset_fdb_root(iaddr1,    cdfdbsf(1:icdfdbsf) )
          IF ( istat .NE. 0 ) THEN
            WRITE(kuso,'("Error\ /Fld2fdb/ Root dir specified: ",a)')
     .      cdfdbsf(1:icdfdbsf)
            WRITE(kuso,'("Error\ /Fld2fdb/ iopenfdb return status:",
     .      I3)') istat
            kerr = 1
            RETURN
          ENDIF
C
          lopenfdb = .TRUE.
C
          if(cltype.eq.'cf'.or.cltype.eq.'pf') then
            istat = isetvalfdb( iaddr1,   'STREAM', 'waef' )
          else
            istat = isetvalfdb( iaddr1,   'STREAM', 'wave' )
          endif

          istat = isetvalfdb( iaddr1,   'LEVTY',   's' )
          istat = isetvalfdb( iaddr1,   'LEVEL',   '0000' )
          istat = isetvalfdb( iaddr1,   'REPR',    'll' )
C
c
c         ensemble number
c
          WRITE( clnumber, '(i3.3)' ) ksec1(42)
          istat = isetvalfdb( iaddr1,   'NUMBER',  clnumber)

        ENDIF
      ELSEIF ( cdfun .EQ. 'close' ) THEN
        istat = iclosefdb( iaddr1 )
        IF ( istat .NE. 0 ) THEN
          WRITE(kuso,'( "Error\ /Fld2fdb/ Action specified: ", a)' )
     .    cdfun
          WRITE(kuso, '( "Error\ /Fld2fdb/ iclosefdb return status:",
     .    I3)') istat
          kerr = 1
          RETURN
        ELSE
          WRITE(kuso, '( "/Fld2fdb/ Closed FDB and DB: " )' )
          kerr = 0
          lopenfdb = .FALSE.
          RETURN
        ENDIF
      ELSE
        WRITE(kuso, '( "Error\ /Fld2fdb/ Unrecognised action",
     .            " specified: ", a )' ) cdfun
        kerr = 1
        RETURN
      ENDIF
C
C
C
C     ------------------------------------------------------------------
C
C*    3. MAKE DATA BASE NAME.
C
C     SET : CLASS - STREAM - EXPVER - DOMAIN - TYPE - PARAMETER
C           LEVELTYPE - LEVEL - REPRESENTATION 
C           DATE - TIME - STEP
C
C     ------------------------------------------------------------------
C
c      WRITE(*,'(2a8)') cl_p_class, clclass
c      WRITE(*,'(2a8)') cl_p_expver, clexpver
c      WRITE(*,'(2a8)') cl_p_domain, cldomain
c      WRITE(*,'(2a8)') cl_p_type, cltype
c      WRITE(*,'(2a8)') cl_p_param, clparam
c      WRITE(*,'(2a8)') cl_p_date, cldate
c      WRITE(*,'(2a8)') cl_p_time, cltime
c      WRITE(*,'(2a8)') cl_p_step, clstep
C
      IF ( cl_p_class .NE. clclass ) THEN
        istat = isetvalfdb( iaddr1,   'CLASS',   clclass )
        cl_p_class=clclass
      ENDIF
C
      IF ( cl_p_expver .NE. clexpver ) THEN
        istat = isetvalfdb( iaddr1,   'expver',  clexpver )
        cl_p_expver=clexpver
      ENDIF
C
      IF ( cl_p_domain .NE. cldomain ) THEN
        istat = isetvalfdb( iaddr1, 'domain', cldomain )
        cl_p_domain=cldomain
      ENDIF
C
      IF ( cl_p_type .NE. cltype ) THEN
        istat = isetvalfdb( iaddr1, 'TY',     cltype )
        cl_p_type=cltype
      ENDIF
C
      IF ( cl_p_param .NE. clparam ) THEN
        istat = isetvalfdb( iaddr1, 'PARAM',  clparam )
        cl_p_param=clparam
      ENDIF
C
      IF ( cl_p_date .NE. cldate ) THEN
        istat = isetvalfdb( iaddr1, 'date',   cldate )
        cl_p_date=cldate
      ENDIF
C
      IF ( cl_p_time .NE. cltime ) THEN
        istat = isetvalfdb( iaddr1, 'time',   cltime )
        cl_p_time=cltime
      ENDIF
C
      IF ( cl_p_step .NE. clstep ) THEN
        istat = isetvalfdb(iaddr1,  'step',   clstep )
        cl_p_step=clstep
      ENDIF
C
C     ------------------------------------------------------------------
C
C*    4.  READ OR WRITE FIELD FROM/TO FDB
C
C     ------------------------------------------------------------------
C
      IF ( cdfun .EQ. 'write' )      THEN
        istat = iwritefdb( iaddr1, kfield, klen ) 
        IF ( istat .NE. 0 ) THEN
          WRITE(kuso, '( "Error\ /Fld2fdb/ Failed to write to fdb")' )
          istat = iclosefdb( iaddr1 )
          kerr = 1
          RETURN
        ENDIF
      ELSEIF ( cdfun .EQ. 'read' )  THEN
        istat = ireadfdb( iaddr1, kfield, klen )
        IF ( istat .NE. 0 ) THEN
          WRITE(kuso, '( "Error\ /Fld2fdb/ Failed to read from fdb")' )
          istat = iclosefdb( iaddr1 )
          kerr = 1
          RETURN
        ENDIF
      ELSE
C       NORMALLY YOU SHOULD NOT ARRIVE AT THIS POINT
        WRITE(kuso, '( "Error\ /Fld2fdb/ Unrecognised action",
     .  " specified: ", a )' ) cdfun
        kerr = 1
        RETURN
      ENDIF
C
C     ------------------------------------------------------------------
C
C*    5.  RETURN
C
C     ------------------------------------------------------------------
C
      kfdb = iaddr1
      RETURN
      END
