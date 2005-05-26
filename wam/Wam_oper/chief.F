      PROGRAM chief 

! ----------------------------------------------------------------------

!**** *CHIEF* - SUPERVISES WAVE MODEL EXECUTION.

!     LIANA ZAMBRESKY      GKSS/ECMWF  JUNE 1989
!     H. GUNTHER           ECMWF       JUNE 1990  MODIFIED FOR CYCLE_4.
!     J. BIDLOT            ECMWF       FEBRUARY 1996-97 MESSAGE PASSING
!     J. DOYLE             ECMWF       OCTOBER 1996 ATMOSPHERIC COUPLING
!     J. BIDLOT            ECMWF       APRIL 97 ADD ZDELATM TO CALL TO
!                                               WAVEMDL
!     B. HANSEN            ECMWF       APRIL 97 SIGNAL HANDLING.
!     S. ABDALLA           ECMWF       OCTOBER 2001 MODIFICATION OF THE
!                                                   CALL TO WAVEMDL

!*    PURPOSE.
!     --------

!       THIS PROGRAM SUPERVISES THE EXECUTION OF THE WAM MODEL.

!**   INTERFACE.
!     ----------

!       IN ORDER FOR THE WAM MODEL TO EXECUTE, IT NEEDS
!       FILES FROM ESSENTIALLY FIVE SOURCES.

!       1. THE UNFORMATED FILES CREATED BY THE JOB PREPROC

!       2. USER INPUT FILE

!       3. THE WIND INPUT FILE.

!       4  THE BOUNDARY VALUE INPUT FILES CREATED BY JOB BOUINT.
!          THESE FILES ARE DYNAMICALLY ASSIGNED.

!       5. THE START FILES:
!          THE RESTART FILES HAVE TO BE CREATED BY JOB
!          PRESET, IF A COLD START HAS TO BE DONE.
!          THESE FILES OR FILES FROM A PREVIOUS MODEL RUN
!          ARE AUTOMATICALLY ASSIGNED. (SEE SUB GSFILE).

!       EXPLANATIONS FOR ALL FILES ARE GIVEN IN DETAIL IN SUB INITMDL

!     LIBRARIES.
!     ----------

!         NONE.

!     METHOD.
!     -------

!       THIS VERSION OF THE WAM MODEL HAS BEEN PRODUCED
!       BY MERGING AND CORRECTLY INTERFACING WHAT USED
!       TO BE THE STAND ALONE PROGRAMS:
!               PREWIND AND THE WAM MODEL.
!       PREWIND REFORMATS WINDS INTO THE WAM MODEL BLOCKED
!       STRUCTURE.  STARTING WITH THE INITIAL SEA STATE
!       FILES, THE WAM MODEL CAN THEN INTEGRATE FORWARD
!       IN TIME, DRIVEN BY THE REFORMATTED WINDS.
!       THE SEA STATE AND RESULT FILES ARE SAVED IN REGULAR
!       INTERVALLS. THE SEA STATE FILE SERVE AS THE INITIAL
!       CONDITION FOR A RESTART.

!       EACH CALL OF THE SUB WAVEMDL INTEGRATES FORWARD IN
!       TIME BY ONE WIND INPUT TIMESTEP OR ONE PROPAGATION
!       TIMESTEP, WHAT EVER IS LONGER.
!       IN THE FIRST CALL TO WAVEMDL AN INITIALIZATION IS
!       DONE IN ADDITION.

!     EXTERNALS.
!     ----------

!       *WAVEMDL*   - SUPERVISES THE OVERALL FLOW THROUGH
!                     THE MAIN MODULES: INITMDL, PREWIND
!                     AND WAMODEL.

!     REFERENCE.
!     ----------

!       EACH MODULE IS OF ITSELF THOROUGHLY DOCUMENTED.

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,THQ      ,AKMEAN
      USE YOWMESPAS, ONLY : LMESSPASS 
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPREVIOUS,NNEXT    ,
     &            NINF     ,NSUP     ,MPMAXLENGTH
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO     ,NIBLO    ,
     &            LL1D
      USE YOWSTAT  , ONLY : CDATEE   ,CDTPRO   ,ISIGHUP  ,ISIGINT 
      USE YOWSPEC  , ONLY : NSTART   ,NEND     ,KLENTOP  ,KLENBOT  ,
     &            NFROMPE  ,NTOPE    ,NIJSTART ,IJTOPE   ,NTOPELST ,
     &            NFROMPELST,
     &            U10NEW   ,U10OLD   ,THWNEW   ,THWOLD   ,USNEW    ,
     &            USOLD    ,Z0NEW    ,Z0OLD    ,TAUW     ,FCONST   ,
     &            ROAIRN   ,ROAIRO   ,ZIDLNEW  ,ZIDLOLD  ,
     &            FL1      ,FL2      ,FL3      ,SL
      USE MPL_MODULE

! ----------------------------------------------------------------------

      IMPLICIT LOGICAL(L)

! DIMENSION DUMMY COUPLED VARIABLES

      INTEGER, PARAMETER:: NLONW=1
      INTEGER, PARAMETER :: NLATW=1
      INTEGER, PARAMETER :: NWVFIELDS=1
      INTEGER, PARAMETER :: NC=1
      INTEGER, PARAMETER :: NR=1
      INTEGER, PARAMETER :: NGPTOTG=NC*NR
      INTEGER, PARAMETER :: NFIELDS=1
! ----------------------------------------------------------------------
      REAL WVFLDG(NLONW,NLATW,NWVFIELDS), ZDELATM(NLATW)

      INTEGER ISEC1(255), ISEC2(22+NR)
      REAL FIELDS(NGPTOTG,NFIELDS)
      INTEGER MASK_IN(NGPTOTG)
      INTEGER MASK_OUT(NLONW,NLATW)
! ----------------------------------------------------------------------

      INTEGER IINITFDB_VPP
      CHARACTER *3 DBNAME
      CHARACTER*12 ZERO,CBEGDAT
      LOGICAL LFDB, LLSTOP, LLWRRE, LLRESTARTED, LLIRANK
      DATA LLSTOP, LLWRRE / 2*.FALSE. /

      time0=-wam_user_clock()
      IU06=6

!   0.0 GET MACHINE ARCHTECTURE TYPE AND FDB OPTION FLAG
!       ------------------------------------------------
!!!!! INPUT UNIT FOR THE USER FILE IS NOW 3
      IU05 = 3

      CALL MPUSERIN(IU05,LMESSPASS,LFDB,LWCOU,LL1D)

      IF(LWCOU) THEN
        IU06=6
        WRITE (IU06,*) ' *********************************'
        WRITE (IU06,*) ' *                               *'
        WRITE (IU06,*) ' * PROBLEM IN CHIEF              *'
        WRITE (IU06,*) ' * ============================= *'
        WRITE (IU06,*) ' * LWCOU=TRUE WITH NO ATM. MODEL *'
        WRITE (IU06,*) ' *            TO COUPLE TO       *'
        WRITE (IU06,*) ' *                               *'
        WRITE (IU06,*) ' *                               *'
        WRITE (IU06,*) ' *********************************'
        CALL ABORT1
      ENDIF

      IF(LMESSPASS) THEN

!     0.1 INITIALISE MESSAGE PASSING PROTOCOL 
!         -----------------------------------

        CALL MPL_INIT(KERROR=KERROR)
        IF(KERROR.LT.0) THEN 
          IU06=6
          WRITE (IU06,*) ' ******************************************'
          WRITE (IU06,*) ' *                                        *'
          WRITE (IU06,*) ' *      FATAL ERROR PROGRAM CHIEF         *'
          WRITE (IU06,*) ' *      =========================         *'
          WRITE (IU06,*) ' *                                        *'
          WRITE (IU06,*) ' *            PROBLEM WITH                *'
          WRITE (IU06,*) ' *      MESSAGE PASSING INITIALISATION    *'
          WRITE (IU06,*) ' *                                        *'
          WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS       *'
          WRITE (IU06,*) ' *                                        *'
          WRITE (IU06,*) ' ******************************************'
          CALL ABORT1
        ELSE

!     0.2 GET RANK AND TOTAL NUMBER OF PROCESSORS  
!         ---------------------------------------

          IRANK = MPL_MYRANK()
          NPROC = MPL_NPROC()
          IF (IRANK .EQ. 1) THEN
            WRITE (IU06,*) ' ******************************************'
            WRITE (IU06,*) ' *                                        *'
            WRITE (IU06,*) ' *               PROGRAM CHIEF            *'
            WRITE (IU06,*) ' *               =============            *'
            WRITE (IU06,*) ' *                                        *'
            WRITE (IU06,*) ' *                                        *'
            WRITE (IU06,*) ' *      MESSAGE PASSING INITIALISATION    *'
            WRITE (IU06,*) ' *                                        *'
            WRITE (IU06,*) ' * RUN ON ', NPROC, ' PEs    THIS: ', IRANK
            WRITE (IU06,*) ' *                                        *'
            WRITE (IU06,*) ' ******************************************'
          ENDIF
        ENDIF

      ELSE
        IRANK = 1 
        NPROC = 1 
      ENDIF

      NPREVIOUS=IRANK-1
      IF(IRANK.EQ.NPROC) THEN
        NNEXT=0
      ELSE
        NNEXT=IRANK+1
      ENDIF

!     0.3 DETERMINE GRID DOMAIN DECOMPOSITION 
!         -----------------------------------

      NPR=NPROC
      LLIRANK=.FALSE.
      CALL MPDECOMP(NPR,MAXLEN,LLIRANK)
      MPMAXLENGTH=MAXLEN

!     1.  ALLOCATE NECESSARY ARRAYS
!         -------------------------

      ALLOCATE (U10OLD(NINF:NSUP,NBLO),THWOLD(NINF:NSUP,NBLO),
     &          USOLD(NINF:NSUP,NBLO),U10NEW(NINF:NSUP),
     &          THWNEW(NINF:NSUP),USNEW(NINF:NSUP),
     &          Z0OLD(NINF:NSUP,NBLO),Z0NEW(NINF:NSUP),
     &          TAUW(NINF:NSUP,NBLO))

      ALLOCATE (ROAIRN(NINF:NSUP),ROAIRO(NINF:NSUP,NBLO),
     &          ZIDLNEW(NINF:NSUP),ZIDLOLD(NINF:NSUP,NBLO))
 
      ALLOCATE (FL1(NINF-1:NSUP,NANG,NFRE),FL3(NINF-1:NSUP,NANG,NFRE))
      IF(.NOT.LMESSPASS) ALLOCATE (FL2(NINF-1:NSUP,NANG,NFRE))
      ALLOCATE (SL(NINF-1:NSUP,NANG,NFRE))
      ALLOCATE (FCONST(NINF:NSUP,NFRE))

      ALLOCATE (EMEAN(NSTART(IRANK):NEND(IRANK)))
      ALLOCATE (FMEAN(NSTART(IRANK):NEND(IRANK)))
      ALLOCATE (THQ(NSTART(IRANK):NEND(IRANK)))
      ALLOCATE (AKMEAN(NSTART(IRANK):NEND(IRANK)))

!     0.4 INITIALIZE SIGNAL HANDLER.
!         --------------------------

      ISIGHUP = 0  !    1 /* hangup */
      ISIGINT = 0  !    2 /* interrupt (rubout) */

!!!!  the call to IFSSIG is specific to signal handling for runs at
!!!!  ECMWF, it can be commented out for other configuration.

      IF(LMESSPASS) THEN
        WRITE(IU06,*) ' INITIALIZE SIGNAL HANDLER on PE ', IRANK 
        CALL IFSSIG (ISIGHUP, ISIGINT, IRANK)
        WRITE(IU06,*) ' SIGNAL HANDLER on PE ', IRANK, ' OK'
      ENDIF  

!     ------------------------------------------------------------------

!*    2. CALLS TO WAVEMDL UNTIL MODEL DATE REACHES END DATE.
!*       EACH CALL INTEGRATES ONE WIND INPUT TIMESTEP, OR ONE
!*       PROPAGATION TIMESTEP, WHAT EVER IS LONGER.
!        ---------------------------------------------------

      ZERO   = ' '
      CDTPRO = ZERO

!* DEFINE DUMMY PARAMETERS TO FILL COUPLED ARRAYS

      IUGB=0
      IVGB=0
      ILEN=1
      CBEGDAT='999999999999'
      IDUM=0
      IDURAT=-99

 20   CONTINUE

      CALL WAVEMDL(CBEGDAT, IDUM, IDURAT, 
     &             NFIELDS, NGPTOTG, NC, NR, FCRANGE,
     &             ISEC1, ISEC2, FIELDS,
     &             NWVFIELDS, WVFLDG,
     &             NLONW, NLATW, LLSTOP, LLWRRE,
     &             LLRESTARTED, ZDELATM, KQGAUSS,
     &             MASK_IN, MASK_OUT)
      IF (LLSTOP) GOTO 30
      IF (CDTPRO.LT.CDATEE) GOTO 20

! DEALLOCATE ARRAYS

 30   CONTINUE

      DEALLOCATE (U10OLD,THWOLD,USOLD,U10NEW,THWNEW,USNEW,
     &            Z0OLD,Z0NEW,TAUW,FL1,FL3,SL,FCONST,
     &            ROAIRN,ROAIRO,ZIDLNEW,ZIDLOLD,
     &            NSTART,NEND,KLENTOP,KLENBOT,NFROMPE,NTOPE,
     &            NIJSTART,IJTOPE,NTOPELST,NFROMPELST)
      IF(ALLOCATED(FL2)) DEALLOCATE(FL2)

!    3.  TERMINATE MESSAGE PASSING PROTOCOL 
!        -----------------------------------

      time=time0+wam_user_clock()
      time=time*1E-06
      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++'
      WRITE (IU06,*) ' + TOTAL USER TIME IN SECONDS +'
      WRITE (IU06,*) ' + ', time 
      IF(LMESSPASS) THEN
        WRITE (IU06,*) ' +                            +'
        WRITE (IU06,*) ' + ON PE : ',IRANK   
      ENDIF
      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++'

      IF(LMESSPASS) THEN
        CALL MPCLOSE_UNIT
        CALL MPL_END
      ENDIF

      STOP
      END
