      SUBROUTINE RUNWAM(LNEMOIO)

! ----------------------------------------------------------------------

!**** *RUNWAM* - SUPERVISES WAVE MODEL EXECUTION.

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

      USE MPL_MPIF
      USE YOWCOUP  , ONLY : LWCOU    ,LWFLUX   ,LWNEMOCOU          ,
     &                      NEMOINIDATE, NEMOINITIME               ,
     &                      NEMOITINI,   NEMOITEND                 ,
     &                      NEMOTSTEP,   NEMOFRCO                  ,
     &                      NEMONSTEP,   NEMOCSTEP, NEMOWSTEP
      USE YOWCURR  , ONLY : U        ,V
      USE YOWICE   , ONLY : CICOVER  ,CITHICK  ,CIWA
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,THQ      ,PHIEPS   ,
     &                      PHIAW    ,TAUOC
      USE YOWMESPAS, ONLY : LMESSPASS 
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPREVIOUS,NNEXT    ,
     &            NINF     ,NSUP     ,MPMAXLENGTH
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO     ,NIBLO    ,
     &            LL1D
      USE YOWSTAT  , ONLY : CDATEE   ,CDTPRO   ,ISIGHUP  ,ISIGINT  ,
     &            IPROPAGS ,LSUBGRID ,IREFRA   ,IDELPRO
      USE YOWSPEC  , ONLY : NSTART   ,NEND     ,KLENTOP  ,KLENBOT  ,
     &            NFROMPE  ,NTOPE    ,NIJSTART ,IJTOPE   ,NTOPELST ,
     &            NFROMPELST,
     &            U10NEW   ,U10OLD   ,THWNEW   ,THWOLD   ,USNEW    ,
     &            USOLD    ,Z0NEW    ,Z0OLD    ,TAUW     ,
     &            ROAIRN   ,ROAIRO   ,ZIDLNEW  ,ZIDLOLD  ,
     &            FL1      ,FL3      ,SL
      USE YOWUNPOOL, ONLY : LLUNSTR
      USE MPL_MODULE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
#if defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
      USE coupling_var, only : WAV_COMM_WORLD
#endif
      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal, MyRankLocal

! ----------------------------------------------------------------------

      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: LNEMOIO

! DIMENSION DUMMY COUPLED VARIABLES
      INTEGER, PARAMETER :: NLONW=1
      INTEGER, PARAMETER :: NLATW=1
      INTEGER, PARAMETER :: NWVFIELDS=8
      INTEGER, PARAMETER :: NC=1
      INTEGER, PARAMETER :: NR=1
      INTEGER, PARAMETER :: NGPTOTG=NC*NR
      INTEGER, PARAMETER :: NFIELDS=7

      INTEGER :: NATMFLX
      INTEGER :: NLON, NLAT
      INTEGER :: IGRIB_HANDLE_DUM
      INTEGER :: NADV
      INTEGER :: KSTOP, KSTPW
      INTEGER :: IDUM
      INTEGER :: KQGAUSS
      INTEGER :: NFDBREF, NPR
      INTEGER :: MAXLEN, IU06
      INTEGER :: KERROR
      INTEGER :: I_GET_UNIT
      INTEGER :: MASK_IN(NGPTOTG)
      INTEGER :: MASK_OUT(NLONW,NLATW)

      REAL :: PSTEP
      REAL :: RSOUTW, RNORTW
      REAL :: PRPLRADI, PRPLRG
      REAL :: WVFLDG(NLONW,NLATW,NWVFIELDS), ZDELATM(NLATW)
      REAL :: FIELDS(NGPTOTG,NFIELDS)
      REAL :: ZHOOK_HANDLE
      REAL :: RMISS
      REAL :: time0, time
      REAL :: wam_user_clock

      CHARACTER(LEN=3) :: DBNAME
      CHARACTER(LEN=14) :: ZERO,CBEGDAT

      LOGICAL :: LLSTOP, LLWRRE, LLRESTARTED, LLIRANK
      LOGICAL :: LDWCOUNORMS
      LOGICAL :: FRSTIME
      LOGICAL :: LWCUR
      LOGICAL :: LWSTOKES
      LOGICAL :: LLRNL
      LOGICAL :: LDWCOU2W, LFDBIFS

      DATA LLSTOP, LLWRRE / 2*.FALSE. /

! ----------------------------------------------------------------------
#ifdef ECMWF 
      IF (LHOOK) CALL DR_HOOK('RUNWAM',0,ZHOOK_HANDLE)
#endif

#if defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
      LMPLUSERCOMM = .TRUE.
      MPLUSERCOMM = WAV_COMM_WORLD
#else
      IF (.NOT.LNEMOIO) THEN
         LMPLUSERCOMM = .TRUE.
         MPLUSERCOMM = MPI_COMM_WORLD
      ENDIF
#endif

      time0=-wam_user_clock()
      IU06=6

!     0.1 INITIALISE MESSAGE PASSING PROTOCOL 
!         -----------------------------------

#if !defined MODEL_COUPLING_ATM_WAV && !defined MODEL_COUPLING_OCN_WAV
      MyRankGlobal=0
#endif

      CALL MPL_INIT(KERROR=KERROR)
      IF(KERROR.LT.0) THEN 
        IU06=6
        WRITE (IU06,*) ' ******************************************'
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' *      FATAL ERROR PROGRAM RUNWAM         *'
        WRITE (IU06,*) ' *      =========================         *'
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' *            PROBLEM WITH                *'
        WRITE (IU06,*) ' *      MESSAGE PASSING INITIALISATION    *'
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS       *'
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' ******************************************'
        CALL ABORT1
      ENDIF

!     0.2 GET MODEL PARAMETERS
!         --------------------

      LWCOU=.FALSE.
      LDWCOU2W=.FALSE.
      LWFLUX=.FALSE. ! will be reset to true if ocean fluxes are output.
      LWCUR=.FALSE. ! only used in coupled runs
      LFDBIFS=.FALSE.


      PRPLRADI=1.0
      PRPLRG=1.0
      LLRNL=.TRUE.

      CALL WVWAMINIT (LWCOU,IU06,LLRNL,NLON,NLAT,RSOUTW,RNORTW)

      CALL WVWAMINIT1 (LWCOU,LDWCOU2W,LWFLUX,LWCUR,LFDBIFS)

!     0.3 DETERMINE GRID DOMAIN DECOMPOSITION 
!         -----------------------------------

      NADV=0
      FRSTIME=.TRUE.

      CALL INIWCST(PRPLRADI)

      CALL WVWAMDECOMP(LWCUR)


!     1.  ALLOCATE NECESSARY ARRAYS
!         -------------------------

      CALL WVALLOC(LWCUR)


!     0.4 INITIALIZE SIGNAL HANDLER.
!         --------------------------

      ISIGHUP = 0  !    1 /* hangup */
      ISIGINT = 0  !    2 /* interrupt (rubout) */

!!!!  the call to IFSSIG and sigmaster are specific to signal handling for runs at
!!!!  ECMWF, it can be commented out for other configuration.
#ifdef ECMWF 
      WRITE(IU06,*) ' INITIALIZE SIGNAL HANDLER on PE ', IRANK 
      CALL IFSSIG (ISIGHUP, ISIGINT, IRANK)

      IF(IRANK.EQ.1) CALL SIGMASTER()

      WRITE(IU06,*) ' SIGNAL HANDLER on PE ', IRANK, ' OK'
#endif


!     ------------------------------------------------------------------

!*    2. CALLS TO WAVEMDL UNTIL MODEL DATE REACHES END DATE.
!*       EACH CALL INTEGRATES ONE WIND INPUT TIMESTEP, OR ONE
!*       PROPAGATION TIMESTEP, WHAT EVER IS LONGER.
!        ---------------------------------------------------

      ZERO   = ' '
      CDTPRO = ZERO

!* DEFINE DUMMY PARAMETERS TO FILL COUPLED ARRAYS

      CBEGDAT='99999999999999'
      PSTEP=0. ! only used in coupled model
      KSTOP=0  ! only used in coupled model
      KSTPW=0  ! only used in coupled model
      IDUM=0
      NFDBREF=0  ! only used in coupled model
      IGRIB_HANDLE_DUM=-99 ! only used in coupled model
      NATMFLX=0
      IF (.NOT.LWNEMOCOU) THEN
        LWCUR=.FALSE. ! only used in coupled runs
        LWSTOKES=.FALSE.  ! only used in coupled runs
      ENDIF
      RMISS=-999. ! missing data indicator

      ! WAM-NEMO COUPLING

      IF (LWNEMOCOU) THEN
        IF(IRANK.EQ.1) WRITE (IU06,*)'CALLING NEMOGCMCOUP_INIT'
#ifdef ECMWF
        CALL NEMOGCMCOUP_INIT( MPL_COMM, NEMOINIDATE, NEMOINITIME, 
     &                         NEMOITINI, NEMOITEND, NEMOTSTEP,
     &                         .TRUE., -1, .FALSE. )
#endif
        IF(IRANK.EQ.1) THEN
          WRITE(IU06,*)'NEMO INITIAL DATE         : ',NEMOINIDATE
          WRITE(IU06,*)'NEMO INITIAL TIME         : ',NEMOINITIME
          WRITE(IU06,*)'NEMO INITIAL STEP         : ',NEMOITINI
          WRITE(IU06,*)'NEMO FINAL STEP           : ',NEMOITEND
          WRITE(IU06,*)'NEMO TIME STEP            : ',NEMOTSTEP
        ENDIF
        IF (NEMOFRCO==0) THEN
          WRITE(IU06,*) ' NEMOFRCO NEEDS TO BE SET IN WAM-NEMO MODE.'
          CALL ABORT1
        ELSE
          NEMOCSTEP=NEMOITINI
          NEMONSTEP=NINT(NEMOFRCO*IDELPRO/NEMOTSTEP)
          NEMOWSTEP=0
        ENDIF
        IF ((NEMOFRCO*IDELPRO)/=NINT(NEMOTSTEP*NEMONSTEP)) THEN
          WRITE(IU06,*)'Inconsistent NEMO and WAM coupling intervals:'
          WRITE(IU06,*)'WAM coupling interval is : ',NEMOFRCO*IDELPRO
          WRITE(IU06,*)'NEMO coupling interval is: ',NEMOTSTEP*NEMONSTEP
          CALL ABORT1
        ELSE
          IF (IRANK.EQ.1) THEN
            WRITE(IU06,*)'NEMONSTEP                 : ',NEMONSTEP
            WRITE(IU06,*)'NEMO coupling interval is : ',
     &                   NEMOTSTEP*NEMONSTEP
          ENDIF
        ENDIF
        LWSTOKES=.TRUE.
      ENDIF

 20   CONTINUE

      CALL WAVEMDL(CBEGDAT, PSTEP, KSTOP, KSTPW,
     &             NFIELDS, NGPTOTG, NC, NR,
     &             IGRIB_HANDLE_DUM, RMISS, FIELDS,
     &             NATMFLX,
     &             LWCUR, LWSTOKES,
     &             NWVFIELDS, WVFLDG,
     &             NLONW, NLATW, LLSTOP, LLWRRE,
     &             LLRESTARTED, ZDELATM, KQGAUSS,
     &             LDWCOUNORMS, MASK_IN, MASK_OUT,
     &             NFDBREF,
     &             FRSTIME, NADV, PRPLRADI, PRPLRG,IDUM)

      IF (LLSTOP) GOTO 30
      IF (CDTPRO.LT.CDATEE) GOTO 20

 30   CONTINUE

!    3.  TERMINATE MESSAGE PASSING PROTOCOL 
!        -----------------------------------

      time=time0+wam_user_clock()
      time=time*1E-06
      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++'
      WRITE (IU06,*) ' + TOTAL USER TIME IN SECONDS +'
      WRITE (IU06,*) ' + ', time 
      WRITE (IU06,*) ' +                            +'
      WRITE (IU06,*) ' + ON PE : ',IRANK   
      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++'

#ifdef WITH_NEMO
      IF (LWNEMOCOU) CALL NEMOGCMCOUP_FINAL
#endif

      CALL MPCLOSE_UNIT
      IF (.NOT.LNEMOIO) CALL MPL_END

#ifdef ECMWF 
      IF (LHOOK) CALL DR_HOOK('RUNWAM',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE RUNWAM
