C ======================================================================
C
      PROGRAM CPFDB
C
C ======================================================================
C
C**** CPFDB
C****           copy data from a given fdb type to another one.
C
C        B. HANSEN   *ECMWF*   13/03/98
C
C     Interface.
C     ----------
C         INPUT.
C         ------
C           *command line*   -d ccssssdxxxxyyyymmdd -t tthh  -s step
C                            -D ccssssdxxxxyyyymmdd -T tthh  -S step
C
C                            -p parameter -l level
C
C            where ccssssdxxxxyyyymmdd may be rdwavesat7k19980224
C            and   tthh                may be an21
C                  step                may be 000001
C                       ie. cc = class
C                       ie. ssss = stream
C                       ie. d = domain
C                       ie. xxxx = experiment id
C                       ie. yyyy = year
C                       ie. mm = month
C                       ie. dd = day
C                       ie. tt = type
C                       ie. hh = time
C             -d from_directory  -t from_type [ -r from_root ]
C             -D to_directory    -t to_type   [ -R to_root   ]
C             [ -x decode ] decode = decode
C
C         OUTPUT.
C         -------
C           *iu_grs*    output unit opened for grs file.
C
C     Method.
C     -------
C       NONE.
C
C     Modification.
C     -------------
C       NONE.
C
C
C ----------------------------------------------------------------------
C
      IMPLICIT NONE
C.....!!!!!!!!!!!!!
C
C ----------------------------------------------------------------------
C
      INTEGER, PARAMETER :: JPUIN  =  5, ! input unit for user information.
     .                      JPUSO  =  6, ! output unit for printer output.
     .                      JPNIB  = 1129120
C
      INTEGER IINITFDB, ISTAT, IOPENFDB, JLEVEL, IERR, ILEVEL,
     .        ISETVALFDB, ILEN, IREAD, IREADFDB, ICLOSEFDB,
     .        ISET_FDB_ROOT, IWRITEFDB, IRANK, ISETRANKFDB
C
      INTEGER,             DIMENSION(2)     :: IUNIT, ILCR
      INTEGER,             DIMENSION(JPNIB) :: INGRIB
C
      CHARACTER (LEN=  2)               :: CLPAR
      CHARACTER (LEN=  4)               :: CLLEVEL
      CHARACTER (LEN=  4), DIMENSION(2) :: CLTYPE
      CHARACTER (LEN=  6), DIMENSION(2) :: CLSTEP
      CHARACTER (LEN= 19), DIMENSION(2) :: CLDIR, CLROOT
C
C
C ----------------------------------------------------------------------
C
C*       PRESET VARIABLES.
C        -----------------
C
      clpar = '250'
      cltype(1) = 'fc18'
      cltype(2) = 'an21'
      cldir(1)  = 'rdwavegat7k19980224'
      cldir(2)  = 'rdwavegat7k19980224'
      clstep(1) = '000003'
      clstep(2) = '000000'
      clroot(1) = '/vpp700/fdb50'
      clroot(2) = '/vpp700/fdb48'
C
C     WRITE(jpuso,'(" cpfdb: Mon Mar 16 19:45:35 GMT 1998 ")')
C
C ----------------------------------------------------------------------
C
C*       CRACK COMMAND LINE ARGUMENTS.
C        -----------------------------
      CALL command (jpuso, cldir, clstep, clroot, cllevel, cltype, ierr)
      IF (ierr.ne.0) CALL abort
C
C
C ----------------------------------------------------------------------
C
C*       LOOP ON ALL LEVELS AND COPY FIELDS.
C        -----------------------------------
C
      READ (cllevel, '(I4.4)') ILEVEL
      ilcr(1)   = LEN_TRIM(clroot(1))
      ilcr(2)   = LEN_TRIM(clroot(2))
      ISTAT = IINITFDB()
      DO JLEVEL = 1, ILEVEL
        ISTAT = IOPENFDB   ('fdb', IUNIT(1), 'r'   )
        IF (ilcr(1) .GT. 1) ISTAT = ISET_FDB_ROOT(IUNIT(1),
     .                              clroot(1)(1:ilcr(1)) )
        WRITE (cllevel, '(I4.4)') JLEVEL
        ISTAT = ISETVALFDB (IUNIT(1), 'levt',   'm'      )
        ISTAT = ISETVALFDB (IUNIT(1), 'param',  '250'    )
        ISTAT = ISETVALFDB (IUNIT(1), 'repres', 'll'     )
        ISTAT = ISETVALFDB (IUNIT(1), 'fstep',  '000000' )
        ISTAT = ISETVALFDB (IUNIT(1), 'class',  cldir(1)( 1: 2))
        ISTAT = ISETVALFDB (IUNIT(1), 'stream', cldir(1)( 3: 6))
        ISTAT = ISETVALFDB (IUNIT(1), 'domain', cldir(1)( 7: 7))
        ISTAT = ISETVALFDB (IUNIT(1), 'expver', cldir(1)( 8:11))
        ISTAT = ISETVALFDB (IUNIT(1), 'date',   cldir(1)(12:19))
        ISTAT = ISETVALFDB (IUNIT(1), 'type',   cltype(1)(1:2))
        ISTAT = ISETVALFDB (IUNIT(1), 'time',   cltype(1)(3:4))
        ISTAT = ISETVALFDB (IUNIT(1), 'step',   clstep(1) )
        ISTAT = ISETVALFDB (IUNIT(1), 'leve',   cllevel   )
        WRITE (jpuso,'("  copy ",a2, a4, a1, a4, a8, "/", a2, a2,
     .                 " step ", a4, " level ", a4,
     .                    " to ",a2, a4, a1, a4, a8, "/", a2, a2
     .                 " step ", a4, " level ", a4)')
     .  cldir(1)( 1: 2), cldir(1)( 3: 6), cldir(1)( 7: 7),
     .  cldir(1)( 8:11), cldir(1)(12:19), cltype(1)(1:2), 
     .  cltype(1)(3:4), clstep(1), cllevel,
     .  cldir(2)( 1: 2), cldir(2)( 3: 6), cldir(2)( 7: 7),
     .  cldir(2)( 8:11), cldir(2)(12:19), cltype(2)(1:2),
     .  cltype(2)(3:4), clstep(2), cllevel
C
        ILEN  = JPNIB
        IREAD = IREADFDB   (IUNIT, INGRIB, ILEN)
        WRITE (JPUSO,'(" READ: LEVEL=",A4, " LENGTH=", i9,
     .               " STATUS=", i3)') cllevel, ilen, iread
C
        IF (IREAD .NE. 0 ) CALL abort
        ISTAT = ICLOSEFDB   (IUNIT(1))
        ISTAT = IOPENFDB   ('fdb', IUNIT(2), 'w'   )
        IF (ilcr(2) .GT. 1) ISTAT = ISET_FDB_ROOT(IUNIT(2),
     .                              clroot(2)(1:ilcr(2)) )
C
        ISTAT = ISETVALFDB (IUNIT(2), 'levt',   'm'      )
        ISTAT = ISETVALFDB (IUNIT(2), 'param',  '250'    )
        ISTAT = ISETVALFDB (IUNIT(2), 'repres', 'll'     )
        ISTAT = ISETVALFDB (IUNIT(2), 'fstep',  '000000' )
        ISTAT = ISETVALFDB (IUNIT(2), 'class',  cldir(2)( 1: 2))
        ISTAT = ISETVALFDB (IUNIT(2), 'stream', cldir(2)( 3: 6))
        ISTAT = ISETVALFDB (IUNIT(2), 'domain', cldir(2)( 7: 7))
        ISTAT = ISETVALFDB (IUNIT(2), 'expver', cldir(2)( 8:11))
        ISTAT = ISETVALFDB (IUNIT(2), 'date',   cldir(2)(12:19))
        ISTAT = ISETVALFDB (IUNIT(2), 'type',   cltype(2)(1:2))
        ISTAT = ISETVALFDB (IUNIT(2), 'time',   cltype(2)(3:4))
        ISTAT = ISETVALFDB (IUNIT(2), 'step',   clstep(2) )
        ISTAT = ISETVALFDB (IUNIT(2), 'leve',   cllevel   )
        ISTAT = ISETVALFDB (IUNIT(2), 'frank',  cllevel(2:4)   )
        IRANK = ISETRANKFDB(IUNIT(2), jlevel)
C
        ISTAT = IWRITEFDB  (IUNIT(2), INGRIB, ILEN)
        IF (ISTAT .NE. 0) CALL abort
        ISTAT = ICLOSEFDB  (IUNIT(2))
C
      ENDDO
C
      WRITE(jpuso,'(" NORMAL END")')
C
      END
