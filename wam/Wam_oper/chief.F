      PROGRAM chief 
C
C ----------------------------------------------------------------------
C
C**** *CHIEF* - SUPERVISES WAVE MODEL EXECUTION.
C
C     LIANA ZAMBRESKY      GKSS/ECMWF  JUNE 1989
C     H. GUNTHER           ECMWF       JUNE 1990  MODIFIED FOR CYCLE_4.
C     J. BIDLOT            ECMWF       FEBRUARY 1996-97 MESSAGE PASSING
C     J. DOYLE             ECMWF       OCTOBER  1996 ATMOSPHERIC COUPLING
C     J. BIDLOT            ECMWF       APRIL 97 ADD ZDELATM TO CALL TO
C                                               WAVEMDL
C     B. HANSEN            ECMWF       APRIL 97 SIGNAL HANDLING.
C
C*    PURPOSE.
C     --------
C
C       THIS PROGRAM SUPERVISES THE EXECUTION OF THE WAM MODEL.
C
C**   INTERFACE.
C     ----------
C
C       IN ORDER FOR THE WAM MODEL TO EXECUTE, IT NEEDS
C       FILES FROM ESSENTIALLY FIVE SOURCES.
C
C       1. THE UNFORMATED FILES CREATED BY THE JOB PREPROC
C
C       2. USER INPUT FILE
C
C       3. THE WIND INPUT FILE.
C
C       4  THE BOUNDARY VALUE INPUT FILES CREATED BY JOB BOUINT.
C          THESE FILES ARE DYNAMICALLY ASSIGNED.
C
C       5. THE START FILES:
C          THE RESTART FILES HAVE TO BE CREATED BY JOB
C          PRESET, IF A COLD START HAS TO BE DONE.
C          THESE FILES OR FILES FROM A PREVIOUS MODEL RUN
C          ARE AUTOMATICALLY ASSIGNED. (SEE SUB GSFILE).
C
C       EXPLANATIONS FOR ALL FILES ARE GIVEN IN DETAIL IN SUB INITMDL
C
C     LIBRARIES.
C     ----------
C
C         NONE.
C
C     METHOD.
C     -------
C
C       THIS VERSION OF THE WAM MODEL HAS BEEN PRODUCED
C       BY MERGING AND CORRECTLY INTERFACING WHAT USED
C       TO BE THE STAND ALONE PROGRAMS:
C               PREWIND AND THE WAM MODEL.
C       PREWIND REFORMATS WINDS INTO THE WAM MODEL BLOCKED
C       STRUCTURE.  STARTING WITH THE INITIAL SEA STATE
C       FILES, THE WAM MODEL CAN THEN INTEGRATE FORWARD
C       IN TIME, DRIVEN BY THE REFORMATTED WINDS.
C       THE SEA STATE AND RESULT FILES ARE SAVED IN REGULAR
C       INTERVALLS. THE SEA STATE FILE SERVE AS THE INITIAL
C       CONDITION FOR A RESTART.
C
C       EACH CALL OF THE SUB WAVEMDL INTEGRATES FORWARD IN
C       TIME BY ONE WIND INPUT TIMESTEP OR ONE PROPAGATION
C       TIMESTEP, WHAT EVER IS LONGER.
C       IN THE FIRST CALL TO WAVEMDL AN INITIALIZATION IS
C       DONE IN ADDITION.
C
C     EXTERNALS.
C     ----------
C
C       *WAVEMDL*   - SUPERVISES THE OVERALL FLOW THROUGH
C                     THE MAIN MODULES: INITMDL, PREWIND
C                     AND WAMODEL.
C
C     REFERENCE.
C     ----------
C
C       EACH MODULE IS OF ITSELF THOROUGHLY DOCUMENTED.
C
C ----------------------------------------------------------------------
C
      USE WAVE_MODULE
C
      IMPLICIT LOGICAL(L)
C
#include "param.h"
C
#include "comcoup.h"
C
#include "commesspass.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
#include "comstat.h"
C
#include "txtstat.h"
C
C ----------------------------------------------------------------------
C
C DIMENSION DUMMY COUPLED VARIABLES
C
      PARAMETER(NLONW=1,NLATW=1)
      INTEGER IUGB(1),IVGB(1)
      REAL BETAG(NLONW,NLATW), ZDELATM(NLATW)
C ----------------------------------------------------------------------
C
      INTEGER IINITFDB_VPP
      CHARACTER *3 DBNAME
      CHARACTER*12 ZERO,CBEGDAT
      LOGICAL LDTRMSG, LFDB, LLSTOP, LLWRRE, LLRESTARTED
      DATA LLSTOP, LLWRRE / 2*.FALSE. /
C
      time0=-wam_user_clock()
      IU06=6
      WRITE(IU06,'("
     . WAM Cycle4 Release8 CI 4 (common with IFS CY18R6)
     .    CC=wam_sth_CY18R5_oam_for_18r6")')
C
C   0.0 GET MACHINE ARCHTECTURE TYPE AND FDB OPTION FLAG
C       ------------------------------------------------
C!!!! INPUT UNIT FOR THE USER FILE IS NOW 3
      IU05 = 3
C
      CALL MPUSERIN(IU05,LMESSPASS,LFDB,LWCOU)
C
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
C
      IF(LMESSPASS) THEN
C
C     0.01 INITIALISE FDB SERVER FOR DISTRIBUTED MEMORY ARCHITECHTURE
C          ----------------------------------------------------------
       IF(LFDB) THEN
        DBNAME='fdb'
        IFDBST = IINITFDB_VPP(DBNAME)
C
        IF(IFDBST.LT.0) THEN 
          IU06=6
          WRITE (IU06,*) ' *********************************'
          WRITE (IU06,*) ' *                               *'
          WRITE (IU06,*) ' * PROBLEM WITH IINITFDB_VPP:    *'
          WRITE (IU06,*) ' * =========================     *'
          WRITE (IU06,*) ' *                               *'
          IF(IFDBST.EQ.-1)
     $    WRITE (IU06,*) ' * USER ERROR                    *'
          IF(IFDBST.EQ.-2)
     $    WRITE (IU06,*) ' * TCP/IP ERROR                  *'
          IF(IFDBST.EQ.-3)
     $    WRITE (IU06,*) ' * MEMORY ERROR                  *'
          WRITE (IU06,*) ' *                               *'
          WRITE (IU06,*) ' *********************************'
          ENDIF
       ENDIF
C
C     0.1 INITIALISE MESSAGE PASSING PROTOCOL 
C         -----------------------------------
C
       CALL MPE_INIT(LDTRMSG,KERROR)
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
C
C     0.2 GET RANK AND TOTAL NUMBER OF PROCESSORS  
C         ---------------------------------------
C
         IRANK = MPE_MYRANK()
         NPROC = MPE_NPROC()
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
C
      ELSE
         IRANK = 1 
         NPROC = 1 
      ENDIF
C
      NPREVIOUS=IRANK-1
      IF(IRANK.EQ.NPROC) THEN
        NNEXT=0
      ELSE
        NNEXT=IRANK+1
      ENDIF
C
C     0.3 DETERMINE GRID DOMAIN DECOMPOSITION 
C         -----------------------------------
      ALLOCATE (NSTART(NPROC),NEND(NPROC),KLENBOT(NPROC),KLENTOP(NPROC))
C
      NPR=NPROC
      CALL MPDECOMP(NPR,NSTART,NEND,KLENTOP,KLENBOT,MAXLEN)
      MPMAXLENGTH=MAXLEN

C
C     1.  ALLOCATE NECESSARY ARRAYS
C         -------------------------
C
      NINF=NSTART(IRANK)-KLENBOT(IRANK)
      NSUP=NEND(IRANK)+KLENTOP(IRANK)
      IF(.NOT.LMESSPASS) NSUP=NIBLO 
C
      ALLOCATE (U10OLD(NINF:NSUP,NBLO),THWOLD(NINF:NSUP,NBLO),
     1          USOLD(NINF:NSUP,NBLO),U10NEW(NINF:NSUP),
     2          THWNEW(NINF:NSUP),USNEW(NINF:NSUP),
     3          Z0OLD(NINF:NSUP,NBLO),Z0NEW(NINF:NSUP),
     4          TAUW(NINF:NSUP,NBLO))
C
      ALLOCATE (FL1(NINF-1:NSUP,NANG,NFRE),FL3(NINF-1:NSUP,NANG,NFRE))
      IF(.NOT.LMESSPASS) ALLOCATE (FL2(NINF-1:NSUP,NANG,NFRE))
      ALLOCATE (SL(NINF-1:NSUP,NANG,NFRE))
      ALLOCATE (FCONST(NINF:NSUP,NFRE))
C
C     0.4 INITIALIZE SIGNAL HANDLER.
C         --------------------------
C
      ISIGHUP = 0  !    1 /* hangup */
      ISIGINT = 0  !    2 /* interrupt (rubout) */
C
      WRITE(IU06,*) ' INITIALIZE SIGNAL HANDLER on PE ', IRANK 
      CALL IFSSIG (ISIGHUP, ISIGINT, IRANK)
      WRITE(IU06,*) ' SIGNAL HANDLER on PE ', IRANK, ' OK'
C
C     ------------------------------------------------------------------
C
C*    2. CALLS TO WAVEMDL UNTIL MODEL DATE REACHES END DATE.
C*       EACH CALL INTEGRATES ONE WIND INPUT TIMESTEP, OR ONE
C*       PROPAGATION TIMESTEP, WHAT EVER IS LONGER.
C        ---------------------------------------------------
C
      ZERO   = ' '
      CDTPRO = ZERO
C
C* DEFINE DUMMY PARAMETERS TO FILL COUPLED ARRAYS
C
      IUGB=0
      IVGB=0
      ILEN=1
      CBEGDAT='999999999999'
      IDUM=0
      IDURAT=-99
C
 20   CONTINUE
C
      CALL WAVEMDL(CBEGDAT, IDUM, IDURAT, IUGB, IVGB, ILEN,
     &             BETAG, NLONW, NLATW, LLSTOP, LLWRRE,
     &             LLRESTARTED, ZDELATM)
      IF (LLSTOP) GOTO 30
      IF (CDTPRO.LT.CDATEE) GOTO 20
C
C DEALLOCATE ARRAYS
C
 30   CONTINUE
C
      DEALLOCATE (U10OLD,THWOLD,USOLD,U10NEW,THWNEW,USNEW,
     1            Z0OLD,Z0NEW,TAUW,FL1,FL3,SL,FCONST,
     2            NSTART,NEND,KLENTOP,KLENBOT)
      IF(ALLOCATED(FL2)) DEALLOCATE(FL2)
C
C    3.  TERMINATE MESSAGE PASSING PROTOCOL 
C        -----------------------------------
C
      time=time0+wam_user_clock()
      time=time*1E-06
      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++'
      WRITE (IU06,*) ' + TOTAL USER TIME IN SECONDS +'
      WRITE (IU06,*) ' +                            +'
      WRITE (IU06,*) ' + ', time 
      IF(LMESSPASS) THEN
      WRITE (IU06,*) ' +                            +'
      WRITE (IU06,*) ' + ON PE : ',IRANK   
      ENDIF
      WRITE (IU06,*) ' +                            +'
      WRITE (IU06,*) ' ++++++++++++++++++++++++++++++'
C
      IF(LMESSPASS) THEN
       CALL MPCLOSE_UNIT
       CALL MPE_END(KERROR)
      ENDIF
C
      STOP
      END
