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

      USE YOWCOUP  , ONLY : LWCOU    ,LWFLUX 
      USE YOWCURR  , ONLY : U        ,V
      USE YOWGRIBHD, ONLY : NKSEK1
      USE YOWICE   , ONLY : ICEMASK
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,THQ      ,PHIEPS   ,
     &                      PHIAW    ,TAUOC
      USE YOWMESPAS, ONLY : LMESSPASS 
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPREVIOUS,NNEXT    ,
     &            NINF     ,NSUP     ,MPMAXLENGTH
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO     ,NIBLO    ,
     &            LL1D
      USE YOWSTAT  , ONLY : CDATEE   ,CDTPRO   ,ISIGHUP  ,ISIGINT  ,
     &            IPROPAGS ,LSUBGRID ,IREFRA 
      USE YOWSPEC  , ONLY : NSTART   ,NEND     ,KLENTOP  ,KLENBOT  ,
     &            NFROMPE  ,NTOPE    ,NIJSTART ,IJTOPE   ,NTOPELST ,
     &            NFROMPELST,
     &            U10NEW   ,U10OLD   ,THWNEW   ,THWOLD   ,USNEW    ,
     &            USOLD    ,Z0NEW    ,Z0OLD    ,TAUW     ,FCONST   ,
     &            ROAIRN   ,ROAIRO   ,ZIDLNEW  ,ZIDLOLD  ,
     &            FL1      ,FL2      ,FL3      ,SL
      USE MPL_MODULE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT LOGICAL(L)

! DIMENSION DUMMY COUPLED VARIABLES

      INTEGER, PARAMETER:: NLONW=1
      INTEGER, PARAMETER :: NLATW=1
      INTEGER, PARAMETER :: NWVFIELDS=8
      INTEGER, PARAMETER :: NC=1
      INTEGER, PARAMETER :: NR=1
      INTEGER, PARAMETER :: NGPTOTG=NC*NR
      INTEGER, PARAMETER :: NFIELDS=7
! ----------------------------------------------------------------------
      INTEGER :: NLON, NLAT
      INTEGER :: IGRIB_HANDLE_DUM
      INTEGER ISEC1(NKSEK1), ISEC2(22+NR)
      INTEGER MASK_IN(NGPTOTG)
      INTEGER MASK_OUT(NLONW,NLATW)

      REAL :: RSOUTW, RNORTW
      REAL WVFLDG(NLONW,NLATW,NWVFIELDS), ZDELATM(NLATW)
      REAL FIELDS(NGPTOTG,NFIELDS)
      REAL :: ZHOOK_HANDLE


! ----------------------------------------------------------------------

      INTEGER IINITFDB_VPP
      INTEGER :: NADV
      CHARACTER *3 DBNAME
      CHARACTER*14 ZERO,CBEGDAT
      LOGICAL LFDB, LLSTOP, LLWRRE, LLRESTARTED, LLIRANK
      LOGICAL :: LDWCOUNORMS
      LOGICAL :: FRSTIME
      LOGICAL :: LMESSPASS_NL
      LOGICAL :: LWCUR
      LOGICAL :: LWSTOKES
      DATA LLSTOP, LLWRRE / 2*.FALSE. /

      time0=-wam_user_clock()
      IU06=6
      LHOOK = .TRUE.
      IF (LHOOK) CALL DR_HOOK('CHIEF',0,ZHOOK_HANDLE)

!   0.0 SET MACHINE ARCHTECTURE TYPE
!       ----------------------------
!!!!! BECAUSE SOME MPI ENVIRONMENT REQUIRES THAT THE VERY FIRST CALL IS TO
!!!!! BE TO MPI_INIT, WE NEED TO ASSUME THAT THE CODE IS RUN UNDER MPI
!!!!! IF IT IS NOT THE CASE CHANGE LMESSPASS TO FALSE.

      LMESSPASS=.TRUE.

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

!     GET CONTROLLING FLAG FROM INPUT NAMELIST 
!     (the rest will be read in userin)
      IU05 =  I_GET_UNIT (IU06, 'wam_namelist', 's', 'f', 0)
      CALL MPUSERIN(IU05,LMESSPASS_NL,LFDB,LWCOU,LL1D,LSUBGRID,
     &              IPROPAGS, IREFRA)

      IF(LMESSPASS_NL.NEQV.LMESSPASS) THEN
        WRITE (IU06,*) ' *********************************'
        WRITE (IU06,*) ' *                               *'
        WRITE (IU06,*) ' *       PROBLEM IN CHIEF        *'
        WRITE (IU06,*) ' * ============================= *'
        WRITE (IU06,*) ' * LMESSPASS AND LMESSPASS_NL    *'
        WRITE (IU06,*) ' *     SHOULD BE THE SAME        *'
        WRITE (IU06,*) ' * LMESPASS    =',LMESSPASS
        WRITE (IU06,*) ' * LMESPASS_NL =',LMESSPASS_NL
        WRITE (IU06,*) ' * LMESPASS_NL IS HARD CODED,    *'
        WRITE (IU06,*) ' * YOU MIGHT HAVE TO CHANGE IT !!*'
        WRITE (IU06,*) ' *                               *'
        WRITE (IU06,*) ' *********************************'
        CALL ABORT1
      ENDIF

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


      NPREVIOUS=IRANK-1
      IF(IRANK.EQ.NPROC) THEN
        NNEXT=0
      ELSE
        NNEXT=IRANK+1
      ENDIF

!     0.2 GET MODEL PARAMETERS
!         --------------------

      CALL WVWAMINIT (LWCOU,IU06,NLON,NLAT,RSOUTW,RNORTW)


!     0.3 DETERMINE GRID DOMAIN DECOMPOSITION 
!         -----------------------------------

      NADV=0
      FRSTIME=.TRUE.

      CALL INIWCST

      NPR=NPROC
      LLIRANK=.FALSE.
      CALL MPDECOMP(NPR,MAXLEN,LLIRANK)
      MPMAXLENGTH=MAXLEN


!     1.  ALLOCATE NECESSARY ARRAYS
!         -------------------------

      LWFLUX=.FALSE. ! will be reset to true if ocean fluxes are output.

      ALLOCATE (U10OLD(NINF:NSUP,NBLO),THWOLD(NINF:NSUP,NBLO),
     &          USOLD(NINF:NSUP,NBLO),U10NEW(NINF:NSUP),
     &          THWNEW(NINF:NSUP),USNEW(NINF:NSUP),
     &          Z0OLD(NINF:NSUP,NBLO),Z0NEW(NINF:NSUP),
     &          TAUW(NINF:NSUP,NBLO))

      ALLOCATE (ROAIRN(NINF:NSUP),ROAIRO(NINF:NSUP,NBLO),
     &          ZIDLNEW(NINF:NSUP),ZIDLOLD(NINF:NSUP,NBLO))
 
      ALLOCATE (ICEMASK(NINF:NSUP,NBLO))

      ALLOCATE (FL1(NINF-1:NSUP,NANG,NFRE),FL3(NINF-1:NSUP,NANG,NFRE))
      IF(.NOT.LMESSPASS) ALLOCATE (FL2(NINF-1:NSUP,NANG,NFRE))
      ALLOCATE (SL(NINF-1:NSUP,NANG,NFRE))
      ALLOCATE (FCONST(NINF:NSUP,NFRE))

      ALLOCATE (EMEAN(NSTART(IRANK):NEND(IRANK)))
      ALLOCATE (FMEAN(NSTART(IRANK):NEND(IRANK)))
      ALLOCATE (THQ(NSTART(IRANK):NEND(IRANK)))

      IF(IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN
        IF(.NOT.ALLOCATED(U)) ALLOCATE(U(NINF-1:NSUP,NBLO))
        U=0.
        IF(.NOT.ALLOCATED(V)) ALLOCATE(V(NINF-1:NSUP,NBLO))
        V=0.
      ENDIF


!     0.4 INITIALIZE SIGNAL HANDLER.
!         --------------------------

      ISIGHUP = 0  !    1 /* hangup */
      ISIGINT = 0  !    2 /* interrupt (rubout) */

!!!!  the call to IFSSIG and sigmaster are specific to signal handling for runs at
!!!!  ECMWF, it can be commented out for other configuration.

      IF(LMESSPASS) THEN
        WRITE(IU06,*) ' INITIALIZE SIGNAL HANDLER on PE ', IRANK 
        CALL IFSSIG (ISIGHUP, ISIGINT, IRANK)

        IF(IRANK.EQ.1) CALL SIGMASTER()

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

      CBEGDAT='99999999999999'
      IDUM=0
      IDURAT=-99
      NFDBREF=0  ! only used in coupled model
      IGRIB_HANDLE_DUM=-99 ! only used in coupled model
      LWCUR=.FALSE. ! only used in coupled runs
      LWSTOKES=.FALSE.  ! only used in coupled runs

      RMISS=-999. ! missing data indicator

 20   CONTINUE

      CALL WAVEMDL(CBEGDAT, IDUM, IDURAT, 
     &             NFIELDS, NGPTOTG, NC, NR,
     &             IGRIB_HANDLE_DUM, ISEC1, ISEC2, RMISS, FIELDS,
     &             LWCUR, LWSTOKES,
     &             NWVFIELDS, WVFLDG,
     &             NLONW, NLATW, LLSTOP, LLWRRE,
     &             LLRESTARTED, ZDELATM, KQGAUSS,
     &             LDWCOUNORMS, MASK_IN, MASK_OUT,
     &             NFDBREF,
     &             FRSTIME, NADV)

      IF (LLSTOP) GOTO 30
      IF (CDTPRO.LT.CDATEE) GOTO 20

! DEALLOCATE ARRAYS

 30   CONTINUE

      DEALLOCATE (U10OLD,THWOLD,USOLD,U10NEW,THWNEW,USNEW,
     &            Z0OLD,Z0NEW,TAUW,ICEMASK,FL1,FL3,SL,FCONST,
     &            ROAIRN,ROAIRO,ZIDLNEW,ZIDLOLD,
     &            NSTART,NEND,KLENTOP,KLENBOT,NFROMPE,NTOPE,
     &            NIJSTART,IJTOPE,NTOPELST,NFROMPELST)
      IF(ALLOCATED(FL2)) DEALLOCATE(FL2)

      DEALLOCATE (EMEAN)
      DEALLOCATE (FMEAN)
      DEALLOCATE (THQ)
      IF (LWFLUX) THEN
         DEALLOCATE (PHIEPS)
         DEALLOCATE (PHIAW)
         DEALLOCATE (TAUOC)
      ENDIF


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

      IF (LHOOK) CALL DR_HOOK('CHIEF',1,ZHOOK_HANDLE)

      STOP
      END
