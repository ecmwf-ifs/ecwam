C=======================================================================
C
      SUBROUTINE FLD2FDB (kuso, ktest, ksec1, kfield, klen, cdfun, kerr,
     .                    cdtype, cdclass, cdexpver, cdfdbsf, cddomain,
     .                    kdir, kfre, kfdb)
C
C=======================================================================
C
C**
C**** NAME    *FLD2FDB*
C**** ----
C**
C**   PURPOSE
C**   -------
C**      WRITE/READ FIELD TO/FROM THE FIELD DATA BASE.
C**
C**   INTERFACE
C**   ---------
C**      INPUT: KUSO      - LOGICAL UNIT FOR STANDARD OUTPUT.
C**             KTEST     - SWITCH DIAGNOSTICS OUTPUT ON IF KTEST GT 1.
C**             KSEC1     - GRIB SECTION 1 OF THE KFIELD.
C**             CDFUN     - FUNCTION: 'WRITE' - WRITE KFIELD INTO FDB.
C**                                   'READ'  - READ FDB FIELD.
C**                                   'CLOSE' - CLOSE FDB.
C**             CDTYPE    - TYPE OF FIELD ( FG, AN OR FC).
C**             CDCLASS   - CLASS OF CURRENT RUN.
C**             CDEXPVER  - EXPERIMENT VERSION NUMBER.
C**             CDDOMAIN  - FDB DOMAIN IDENTIFIER.
C**             KDIR      - CURRENT DIRECTION IF PARAM=251, OTHERWISE SET TO 0 
C**             KFRE      - CURRENT FREQUENCY IF PARAM=251, OTHERWISE SET TO 0 
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
C**      FIELD DATA BASE ROUTINES ARE USED TO OPEN THE FIELD
C**      DATA BASE AND FILE, AND TO WRITE/READ THE FIELD INTO/
C**      FROM THE FIELD DATA BASE.
C**
C**              DATABASE NAME:'ODXXTYYYYYMMDDTT', WHERE
C**                             TY     - DATA TYPE ('AN','FC')
C**                             YYYYMMDD - DATE (YEAR,MONTH,DAY)
C**                             TT     - BASE TIME
C**
C**              FILE NAMES:   'REXXX0TTT00', WHERE
C**                             RE     - DATA REPRESENTATION 'LL'
C**                             TTTTTT - TIME STEP
C**
C**              FIELD NAMES:  'PNNNLLLLHHHHXX', WHERE
C**                             NNN  - PARAMTER NUMBER
C**                             LLLL - LEVEL TYPE 'SFCX'
C**                             HHHH - LEVEL HIGHT '0000'
C**
C**   EXTERNALS  $FDBLIB
C**   ---------
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
C        J. BIDLOT ECMWF NOV 1997 : allow new stream waef and type
C                                   cf or pf from ensemble runs.
C
C        B. HANSEN,  ECMWF, FEB-1998
C          MOVE INITIALIZATION TO USERIN AT THE BOTTOM OF SECTION 1.4. 
C
C     ------------------------------------------------------------------
C
C*        1. VARIABLE DECLARATIONS
C
C     ------------------------------------------------------------------
C
      DIMENSION ksec1(*), kfield(*)
C
      INTEGER       klen, kdir, kfre 
C
      CHARACTER*(*) cdfun
      CHARACTER*(*) cdfdbsf
      CHARACTER*(*) cddomain
C
      CHARACTER*1  cldomain
C
      CHARACTER*2  cltime
      CHARACTER*2  cdtype
      CHARACTER*2  cltype
      CHARACTER*2  cdclass
      CHARACTER*2  clclass
      CHARACTER*2  cldir
      CHARACTER*2  clfre
C
      CHARACTER*3  clparam
      CHARACTER*3  clnumber
C
      CHARACTER*4  cdexpver
      CHARACTER*4  clexpver
      CHARACTER*6  clstep
      CHARACTER*4  clstream
C
      CHARACTER*8  cldate
C
      CHARACTER*1  cl_p_domain
      CHARACTER*2  cl_p_time
      CHARACTER*2  cl_p_type
      CHARACTER*2  cl_p_class
      CHARACTER*3  cl_p_param
      CHARACTER*4  cl_p_expver
      CHARACTER*6  cl_p_step
      CHARACTER*8  cl_p_date
C
      DATA cl_p_domain, cl_p_time, cl_p_type, cl_p_class, cl_p_param,
     .     cl_p_expver, cl_p_step, cl_p_date / 8*"  " /
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
C*        1.   DERIVE FIELD DATA BASE VARIABLES FROM ARGUMENT LIST.
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
        cldomain = cddomain
C
        IF ( cltype .EQ. 'fg' .and. ksec1(16) .EQ. 0 ) then
          clstep = '000006'
        ELSE
          WRITE(clstep, '( i6.6 )' ) ksec1(16)
        ENDIF
C
        WRITE( clparam, '(i3.3)' ) ksec1(6)
        iyyyy=(ksec1(21)-1)*100+ksec1(10)
        WRITE( cldate(01:04), '(i4.4)' ) iyyyy 
        WRITE( cldate(05:06), '(i2.2)' ) ksec1(11)
        WRITE( cldate(07:08), '(i2.2)' ) ksec1(12)
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
C*        2.   OPEN - CLOSE - CHECK ACTION REQUESTED
C
C     ------------------------------------------------------------------
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
          ELSEIF (KTEST .GT. 1) THEN
            WRITE(kuso,'("\ /Fld2fdb/ iopenfdb status:", i3)') istat
          ENDIF
C
          icdfdbsf = LEN_TRIM(cdfdbsf)
          IF (icdfdbsf .GT. 0 ) THEN
            istat = iset_fdb_root(iaddr1,    cdfdbsf(1:icdfdbsf) )
            IF ( istat .NE. 0 ) THEN
              WRITE(kuso,'("Error\ /Fld2fdb/ Root dir specified: ",a)')
     .        cdfdbsf(1:icdfdbsf)
              WRITE(kuso,'("Error\ /Fld2fdb/ iopenfdb return status:",
     .        I3)') istat
              kerr = 1
              RETURN
            ELSEIF (KTEST .GT. 1) THEN
              WRITE(kuso,'("\ /Fld2fdb/ iset_fdb_root status:", i3,
     .        " Root: ", a25)') istat, cdfdbsf(1:icdfdbsf)
            ENDIF
          ENDIF
C
          lopenfdb = .TRUE.
C
          if(cltype.eq.'cf'.or.cltype.eq.'pf') then
            clstream = 'waef'
          else
            clstream = 'wave'
          endif
C
          istat = isetvalfdb( iaddr1,   'STREAM', clstream )
C
          istat = isetvalfdb( iaddr1,   'REPR',    'll' )
C
C         ensemble number
C
          WRITE( clnumber, '(i3.3)' ) ksec1(42)
          istat = isetvalfdb( iaddr1,   'NUMBER',  clnumber)
C
          IF (KTEST.GT.1) WRITE(kuso,'("\ /Fld2fdb/ STREAM:", a5,
     .                    " NUMBER ", a3)') clstream, clnumber

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
          IF (KTEST.GT.1)
     .     WRITE(kuso, '( "/Fld2fdb/ Closed FDB and DB: " )' )
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
      IF ( cl_p_class .NE. clclass ) THEN
        cl_p_class=clclass
      ENDIF
C
      IF ( cl_p_expver .NE. clexpver ) THEN
        cl_p_expver=clexpver
      ENDIF
C
      IF ( cl_p_domain .NE. cldomain ) THEN
        cl_p_domain=cldomain
      ENDIF
C
      IF ( cl_p_type .NE. cltype ) THEN
        cl_p_type=cltype
      ENDIF
C
      IF ( cl_p_param .NE. clparam ) THEN
        cl_p_param=clparam
      ENDIF
C
      IF ( cl_p_date .NE. cldate ) THEN
        cl_p_date=cldate
      ENDIF
C
      IF ( cl_p_time .NE. cltime ) THEN
        cl_p_time=cltime
      ENDIF
C
      IF ( cl_p_step .NE. clstep ) THEN
        cl_p_step=clstep
      ENDIF
C
      istat = isetvalfdb( iaddr1, 'class',  cl_p_class  )
      istat = isetvalfdb( iaddr1, 'expver', cl_p_expver )
      istat = isetvalfdb( iaddr1, 'domain', cl_p_domain )
      istat = isetvalfdb( iaddr1, 'type',   cl_p_type   )
      istat = isetvalfdb( iaddr1, 'param',  cl_p_param  )
      istat = isetvalfdb( iaddr1, 'date',   cl_p_date   )
      istat = isetvalfdb( iaddr1, 'time',   cl_p_time   )
      istat = isetvalfdb( iaddr1, 'step',   cl_p_step   )
      istat = isetvalfdb( iaddr1, 'levty',   's'        )
C     istat = isetvalfdb( iaddr1, 'level',   '0000'     )
C     istat = isetvalfdb( iaddr1, 'fstep',  '000000'    )
C
      IF (kdir.NE.0 .AND. kfre.NE.0 .AND. cl_p_param .EQ. '251' ) THEN
          idir = kdir 
          ifre = kfre 
          WRITE( cldir, '(i2.2)' ) idir 
          WRITE( clfre, '(i2.2)' ) ifre 
          istat = isetvalfdb( iaddr1,   'direction', cldir )
          istat = isetvalfdb( iaddr1,   'frequency', clfre )
      ELSE
          cldir="00"
          clfre="00"
          istat = isetvalfdb( iaddr1,   'direction', cldir )
          istat = isetvalfdb( iaddr1,   'frequency', clfre )
      ENDIF
C
      IF (KTEST.GT.1) THEN
        WRITE(kuso,'("\ /Fld2fdb/ CLASS: ", a8)') cl_p_class
        WRITE(kuso,'("\ /Fld2fdb/ EXPVER:", a8)') cl_p_expver
        WRITE(kuso,'("\ /Fld2fdb/ DOMAIN:", a8)') cl_p_domain
        WRITE(kuso,'("\ /Fld2fdb/ TY:    ", a8)') cl_p_type
        WRITE(kuso,'("\ /Fld2fdb/ PARAM: ", a8)') cl_p_param
        WRITE(kuso,'("\ /Fld2fdb/ DATE:  ", a8)') cl_p_date
        WRITE(kuso,'("\ /Fld2fdb/ TIME:  ", a8)') cl_p_time
        WRITE(kuso,'("\ /Fld2fdb/ STEP:  ", a8)') cl_p_step
        WRITE(kuso,'("\ /Fld2fdb/ DIR:   ", a8)') cldir
        WRITE(kuso,'("\ /Fld2fdb/ FREQ:  ", a8)') clfre
        CALL flush (kuso)
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
        IF (KTEST.GT.1) WRITE(kuso,'("\ /Fld2fdb/ iwritefdb klen:", i10,
     .                  " status:", i4)') klen, istat
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
