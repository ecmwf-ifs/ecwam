      SUBROUTINE WAVEMODEL(IBEGDAT,IDURAT,IUGB,IVGB,ILEN,BETAG,NLONW
     1                  ,NLATW)
C
C ----------------------------------------------------------------------
C
C**** *WAVEMODEL* - SUBTROUTINE CALL OF WAVE MODEL FOR COUPLING TO IFS.
C
C     LIANA ZAMBRESKY      GKSS/ECMWF  JUNE 1989
C     H. GUNTHER           ECMWF       JUNE 1990  MODIFIED FOR CYCLE_4.
!     J. BIDLOT            ECMWF       FEBRUARY 1996 MEASSAGE PASSING
C     J. DOYLE             ECMWF       OCTOBER  1996 ATMOSPHERIC COUPLING
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
!
#include "param.h"
C
!
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

! DEFINE THE ALLOCATABLE ARRAYS
 
      INTEGER,ALLOCATABLE :: NSTART(:),NEND(:),KLENTOP(:),KLENBOT(:)
      REAL,ALLOCATABLE :: U10NEW(:), U10OLD(:,:),
     1                    THWNEW(:), THWOLD(:,:),
     2                    USNEW (:), USOLD (:,:),
     3                    Z0NEW (:), Z0OLD (:,:), TAUW(:,:)
      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: FL1, FL2, FL3

C
C ----------------------------------------------------------------------
      INTEGER IINITFDB_VPP
      CHARACTER *3 DBNAME
      CHARACTER*10 ZERO
      LOGICAL LDTRMSG,LFDB
C COUPLING ARRAYS
      REAL BETAG(NLONW,NLATW)
      INTEGER IUGB(ILEN),IVGB(ILEN)     

      time0=-user_clock()
      IU06=6
C
!   0.0 GET MACHINE ARCHTECTURE TYPE AND FDB OPTION FLAG
!       ------------------------------------------------
!!!!! INPUT UNIT FOR THE USER FILE IS NOW 3
      IU05 = 3

      CALL MPUSERIN(IU05,LMESSPASS,LFDB,LWCOU)

      IF(LWCOU) THEN
        IF(NLONW.NE.NGX.AND.NLATW.NE.NGY) THEN
          IU06=6
          WRITE (IU06,*) ' *********************************'
          WRITE (IU06,*) ' *                               *'
          WRITE (IU06,*) ' * PROBLEM IN WAVEMODEL........  *'
          WRITE (IU06,*) ' * PROBLEM WITH NLONW AND NLATW  *'
          WRITE (IU06,*) ' * NOT EQUAL TO NGX   AND NGY  : *'
          WRITE (IU06,*) ' * ============================= *'
          WRITE (IU06,*) ' *                               *'
          WRITE (IU06,*) ' * NLONW=',NLONW
          WRITE (IU06,*) ' * NLATW=',NLATW
          WRITE (IU06,*) ' * NGX=',NGX  
          WRITE (IU06,*) ' * NGY=',NGY  
          WRITE (IU06,*) ' *                               *'
          WRITE (IU06,*) ' *                               *'
          WRITE (IU06,*) ' *********************************'
          CALL ABORT1
        ENDIF
      ENDIF      

      IF(LMESSPASS) THEN

C     0.01 INITIALISE FDB SERVER FOR DISTRIBUTED MEMORY ARCHITECHTURE
C          ----------------------------------------------------------
        IF(.NOT.LWCOU) THEN
          IF(LFDB) THEN
            DBNAME='fdb'
            ISTAT = IINITFDB_VPP(DBNAME)
  
            IF(ISTAT.LT.0) THEN 
              IU06=6
              WRITE (IU06,*) ' *********************************'
              WRITE (IU06,*) ' *                               *'
              WRITE (IU06,*) ' * PROBLEM WITH IINITFDB_VPP:    *'
              WRITE (IU06,*) ' * =========================     *'
              WRITE (IU06,*) ' *                               *'
              IF(ISTAT.EQ.-1)
     $        WRITE (IU06,*) ' * USER ERROR                    *'
              IF(ISTAT.EQ.-2)
     $        WRITE (IU06,*) ' * TCP/IP ERROR                  *'
              IF(ISTAT.EQ.-3)
     $        WRITE (IU06,*) ' * MEMORY ERROR                  *'
              WRITE (IU06,*) ' *                               *'
              WRITE (IU06,*) ' *********************************'
            ENDIF
          ENDIF

!     0.1 INITIALISE MESSAGE PASSING PROTOCOL 
!         -----------------------------------

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
          ENDIF
        ENDIF

!     0.2 GET RANK AND TOTAL NUMBER OF PROCESSORS  
!         ---------------------------------------

         IRANK = MPE_MYRANK()
         NPROC = MPE_NPROC()
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
      ALLOCATE (NSTART(NPROC),NEND(NPROC),KLENBOT(NPROC),KLENTOP(NPROC))

      CALL MPDECOMP(NSTART,NEND,KLENTOP,KLENBOT)

!     1.  ALLOCATE NECESSARY ARRAYS
!         -------------------------

      NINF=NSTART(IRANK)-KLENBOT(IRANK)
      NSUP=NEND(IRANK)+KLENTOP(IRANK)
      IF(.NOT.LMESSPASS) NSUP=NIBLO 

      ALLOCATE (U10OLD(NINF:NSUP,NBLO),THWOLD(NINF:NSUP,NBLO),
     1          USOLD(NINF:NSUP,NBLO),U10NEW(NINF:NSUP),
     2          THWNEW(NINF:NSUP),USNEW(NINF:NSUP),
     3          Z0OLD(NINF:NSUP,NBLO),Z0NEW(NINF:NSUP),
     4          TAUW(NINF:NSUP,NBLO))

      ALLOCATE (FL1(NINF-1:NSUP,NANG,NFRE),FL2(NINF-1:NSUP,NANG,NFRE),
     1          FL3(NINF-1:NSUP,NANG,NFRE))

C
C*    2. CALLS TO WAVEMDL UNTIL MODEL DATE REACHES END DATE.
C*       EACH CALL INTEGRATES ONE WIND INPUT TIMESTEP, OR ONE
C*       PROPAGATION TIMESTEP, WHAT EVER IS LONGER.

C        ---------------------------------------------------
C
      ZERO   = ' '
      CDTPRO = ZERO
 20   CONTINUE
      CALL WAVEMDL(FL1,FL2,FL3,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     1             U10NEW,THWNEW,USNEW,Z0NEW,
     2             NSTART,NEND,KLENTOP,KLENBOT,
     3             IBEGDAT,IDURAT,IUGB,IVGB,ILEN,BETAG,NLONW,NLATW)
      IF (CDTPRO.LT.CDATEE) GOTO 20

! DEALLOCATE ARRAYS

      DEALLOCATE (U10OLD,THWOLD,USOLD,U10NEW,THWNEW,USNEW,
     1            Z0OLD,Z0NEW,TAUW,FL1,FL2,FL3,
     2            NSTART,NEND,KLENTOP,KLENBOT)

!    3.  TERMINATE MESSAGE PASSING PROTOCOL 
!        -----------------------------------

      time=time0+user_clock()
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

      IF(LMESSPASS.AND..NOT.LWCOU) THEN
       CALL MPCLOSE_UNIT
       CALL MPE_END(KERROR)
      ENDIF

      RETURN
      END
