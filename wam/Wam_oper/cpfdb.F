! ======================================================================

      PROGRAM CPFDB

! ======================================================================

!**** CPFDB
!****           copy data from a given fdb type to another one.

!        B. HANSEN   *ECMWF*   13/03/98

!     Interface.
!     ----------
!         INPUT.
!         ------
!           *command line*   -d ccssssdxxxxyyyymmdd -t tthh  -s step
!                            -D ccssssdxxxxyyyymmdd -T tthh  -S step

!                            -p parameter -l level

!            where ccssssdxxxxyyyymmdd may be rdwavesat7k19980224
!            and   tthh                may be an21
!                  step                may be 000001
!                       ie. cc = class
!                       ie. ssss = stream
!                       ie. d = domain
!                       ie. xxxx = experiment id
!                       ie. yyyy = year
!                       ie. mm = month
!                       ie. dd = day
!                       ie. tt = type
!                       ie. hh = time
!             -d from_directory  -t from_type [ -r from_root ]
!             -D to_directory    -t to_type   [ -R to_root   ]
!             [ -x decode ] decode = decode

!         OUTPUT.
!         -------
!           *iu_grs*    output unit opened for grs file.

!     Method.
!     -------
!       NONE.

!     Modification.
!     -------------
!       NONE.


! ----------------------------------------------------------------------

      USE FDBSUBS
      IMPLICIT NONE
!.....!!!!!!!!!!!!!

! ----------------------------------------------------------------------

      INTEGER, PARAMETER :: JPUIN =  5,! input unit for user information
     &                      JPUSO =  6,! output unit for printer output
     &                      JPNIB = 1129120

      INTEGER ISTAT, JLEVEL, IERR, ILEVEL,
     &        ILEN, IREAD,
     &        IRANK
      INTEGER jlevel_temp

      INTEGER KFIELD

      INTEGER,             DIMENSION(2)     :: IUNIT, ILCR
      INTEGER,             DIMENSION(JPNIB) :: INGRIB

      CHARACTER (LEN=  2)               :: CLPAR
      CHARACTER (LEN=  4)               :: CLLEVEL
      CHARACTER (LEN=  4), DIMENSION(2) :: CLTYPE
      CHARACTER (LEN=  6), DIMENSION(2) :: CLSTEP
      CHARACTER (LEN= 19), DIMENSION(2) :: CLDIR, CLROOT


! ----------------------------------------------------------------------

!*       PRESET VARIABLES.
!        -----------------

      clpar = '250'
      cltype(1) = 'fc18'
      cltype(2) = 'an21'
      cldir(1)  = 'rdwavegat7k19980224'
      cldir(2)  = 'rdwavegat7k19980224'
      clstep(1) = '000003'
      clstep(2) = '000000'
      clroot(1) = '/vpp700/fdb50'
      clroot(2) = '/vpp700/fdb48'

!     WRITE(jpuso,'(" cpfdb: Mon Mar 16 19:45:35 GMT 1998 ")')

! ----------------------------------------------------------------------

!*       CRACK COMMAND LINE ARGUMENTS.
!        -----------------------------
      CALL command (jpuso, cldir, clstep, clroot, cllevel, cltype, ierr)
      IF (ierr.ne.0) CALL abort


! ----------------------------------------------------------------------

!*       LOOP ON ALL LEVELS AND COPY FIELDS.
!        -----------------------------------

      READ (cllevel, '(I4.4)') ILEVEL
      ilcr(1)   = LEN_TRIM(clroot(1))
      ilcr(2)   = LEN_TRIM(clroot(2))
      ISTAT = IINITFDBSUBS()
      KFIELD=0
      DO JLEVEL = 1, ILEVEL
        ISTAT = IOPENFDBSUBS   ('fdb', IUNIT(1), 'r'   )
        IF (ilcr(1) .GT. 1) ISTAT = ISET_FDBSUBS_ROOT(IUNIT(1),
     &                              clroot(1)(1:ilcr(1)) )
        WRITE (cllevel, '(I4.4)') JLEVEL
        ISTAT = ISETVALFDBSUBS (IUNIT(1), 'levt',   'm'      )
        ISTAT = ISETVALFDBSUBS (IUNIT(1), 'param',  '250'    )
        ISTAT = ISETVALFDBSUBS (IUNIT(1), 'repres', 'll'     )
        ISTAT = ISETVALFDBSUBS (IUNIT(1), 'fstep',  '000000' )
        ISTAT = ISETVALFDBSUBS (IUNIT(1), 'class',  cldir(1)( 1: 2))
        ISTAT = ISETVALFDBSUBS (IUNIT(1), 'stream', cldir(1)( 3: 6))
        ISTAT = ISETVALFDBSUBS (IUNIT(1), 'domain', cldir(1)( 7: 7))
        ISTAT = ISETVALFDBSUBS (IUNIT(1), 'expver', cldir(1)( 8:11))
        ISTAT = ISETVALFDBSUBS (IUNIT(1), 'date',   cldir(1)(12:19))
        ISTAT = ISETVALFDBSUBS (IUNIT(1), 'type',   cltype(1)(1:2))
        ISTAT = ISETVALFDBSUBS (IUNIT(1), 'time',   cltype(1)(3:4))
        ISTAT = ISETVALFDBSUBS (IUNIT(1), 'step',   clstep(1) )
        ISTAT = ISETVALFDBSUBS (IUNIT(1), 'leve',   cllevel   )
        WRITE (jpuso,'("  copy ",a2, a4, a1, a4, a8, "/", a2, a2,
     &                 " step ", a4, " level ", a4,
     &                    " to ",a2, a4, a1, a4, a8, "/", a2, a2,
     &                 " step ", a4, " level ", a4)')
     &  cldir(1)( 1: 2), cldir(1)( 3: 6), cldir(1)( 7: 7),
     &  cldir(1)( 8:11), cldir(1)(12:19), cltype(1)(1:2), 
     &  cltype(1)(3:4), clstep(1), cllevel,
     &  cldir(2)( 1: 2), cldir(2)( 3: 6), cldir(2)( 7: 7),
     &  cldir(2)( 8:11), cldir(2)(12:19), cltype(2)(1:2),
     &  cltype(2)(3:4), clstep(2), cllevel

        ILEN  = JPNIB
        IREAD = IREADFDBSUBS   (IUNIT(1), INGRIB, ILEN)
        WRITE (JPUSO,'(" READ: LEVEL=",A4, " LENGTH=", i9,
     &               " STATUS=", i3)') cllevel, ilen, iread

        IF (IREAD .NE. 0 ) CALL abort
        ISTAT = ICLOSEFDBSUBS   (IUNIT(1))
        ISTAT = IOPENFDBSUBS   ('fdb', IUNIT(2), 'w'   )
        IF (ilcr(2) .GT. 1) ISTAT = ISET_FDBSUBS_ROOT(IUNIT(2),
     &                              clroot(2)(1:ilcr(2)) )

        ISTAT = ISETVALFDBSUBS (IUNIT(2), 'levt',   'm'      )
        ISTAT = ISETVALFDBSUBS (IUNIT(2), 'param',  '250'    )
        ISTAT = ISETVALFDBSUBS (IUNIT(2), 'repres', 'll'     )
        ISTAT = ISETVALFDBSUBS (IUNIT(2), 'fstep',  '000000' )
        ISTAT = ISETVALFDBSUBS (IUNIT(2), 'class',  cldir(2)( 1: 2))
        ISTAT = ISETVALFDBSUBS (IUNIT(2), 'stream', cldir(2)( 3: 6))
        ISTAT = ISETVALFDBSUBS (IUNIT(2), 'domain', cldir(2)( 7: 7))
        ISTAT = ISETVALFDBSUBS (IUNIT(2), 'expver', cldir(2)( 8:11))
        ISTAT = ISETVALFDBSUBS (IUNIT(2), 'date',   cldir(2)(12:19))
        ISTAT = ISETVALFDBSUBS (IUNIT(2), 'type',   cltype(2)(1:2))
        ISTAT = ISETVALFDBSUBS (IUNIT(2), 'time',   cltype(2)(3:4))
        ISTAT = ISETVALFDBSUBS (IUNIT(2), 'step',   clstep(2) )
        ISTAT = ISETVALFDBSUBS (IUNIT(2), 'leve',   cllevel   )
        ISTAT = ISETVALFDBSUBS (IUNIT(2), 'frank',  cllevel(2:4)   )
        jlevel_temp=jlevel
        IRANK = ISETRANKFDBSUBS(IUNIT(2), jlevel_temp)

        ISTAT = IWRITEFDBSUBS  (IUNIT(2), INGRIB, ILEN)
        KFIELD=KFIELD+1
        IF (ISTAT .NE. 0) CALL abort
        ISTAT = ICLOSEFDBSUBS  (IUNIT(2))

      ENDDO

      WRITE(jpuso,'(" NORMAL END")')

      END
