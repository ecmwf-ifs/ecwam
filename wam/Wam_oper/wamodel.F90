SUBROUTINE WAMODEL (NADV, LDSTOP, LDWRRE, BLK2GLO,             &
 &                  WVENVI, WVPRPT, FF_NOW, FF_NEXT, INTFLDS,  &
 &                  WAM2NEMO, NEMO2WAM, FL1)

! ----------------------------------------------------------------------

!**** *WAMODEL* - 3-G WAM MODEL - TIME INTEGRATION OF WAVE FIELDS.
!                                 AND OUTPUTS

!     S.D. HASSELMANN  MPI       1.12.85

!     G. KOMEN         KNMI         6.86  MODIFIED FOR SHALLOW WATER
!     P. JANSSEN                          ASPECTS.

!     S.D. HASSELMANN  MPI       15.2.87  MODIFIED FOR CYBER 205.

!     P. LIONELLO      ISDGM      6.3.87  MODIFIED TO OUTPUT SWELL.

!     S.D. HASSELMANN  MPI        1.6.87  ALL VERSIONS COMBINED INTO
!                                         ONE MODEL. DEEP AND SHALLOW
!                                         WATER , CRAY AND CYBER 205
!                                         VERSION.

!     CYCLE_2 MODICIFATIONS:
!     ----------------------

!     L. ZAMBRESKY     GKSS        10.87  OPTIMIZED FOR CRAY, CYBER 205
!     H. GUNTHER

!     A. SPEIDEL       MPI          4.88  VARIABLE DIMENSIONS, INTERNAL
!                                         CHECKS (CFL-CRITERION).

!     A. SPEIDEL       MPI         11.88  CHANGES FOR CRAY-2.

!     K. HUBBERT       POL          6.89  DEPTH AND CURRENT REFRACTION.
!                                         PRECALCULATION OF TERMS IN
!                                         *PROPDOT*.
!                                         SOLVE WAVE ACTION EQUATION
!                                         FOR CURRENT REFRACTION.

!     CYCLE_3 MODICIFATIONS:
!     ----------------------

!     R. PORTZ , S.D. HASSELMANN   MPI          1990

!      - RESTRUCTURE MODEL TO CALL THE ACTUAL INTEGRATION IN TIME
!        AS A SUBROUTINE: WAMODEL. A SHELL PROGRAM "WAMSHELL" READS
!        OUTPUT FROM PREPROC AND COMPUTES THE WIND ARRAYS FOR THE
!        INTEGRATION PERIOD FROM PREWIND, WHICH HAS BEEN INCORPORATED
!        AS A SUBROUTINE.
!      - ALL INTERMEDIATE AND RESTART I/O IS DONE IN THE SUBROUTINE
!        WAMODEL AND INPREST.
!      - THE YOWMON BLOCK IN THE PREPROCESSOR AND MODEL ARE MADE
!        COMPATIBLE.
!      - THE YOWPUTATION OF SEVERAL PARAMETERS HAS BEEN TRANSFERRED
!        FROM THE MODEL TO PREPROC.
!      - DEPTH AND CURRENT REFRACTION HAS BEEN INCORPORATED INTO THE
!        MODEL.
!      - OPEN BOUNDARIES ARE INCORPORATED IN THE MODEL.
!      - SEVERAL MINOR ERRORS HAVE BEEN REMOVED.
!      - THE BUFFERED I/O FOR THE CYBER 205 HAS BEEN CHANGED INTO A
!        BINARY READ AND WRITE.

!     CYCLE_4 MODICIFATIONS:
!     ----------------------

!     L. ZAMBRESKY   GKSS/ECMWF   6.89  ECMWF SUB VERSION
!                                       BASED ON CYCLE_2.

!     H. GUNTHER     GKSS/ECMWF 10.89  ECMWF SUB VERSION REORGANIZED.
!                                      - COMMON BLOCK STRUCTURE.
!                                      - BLOCKING STRUCTURE.
!                                      - TIME COUNTING.
!                                      - GRIDDED OUTPUT FIELDS.
!                                      - HEADERS ADDED TO OUTPUT FILES.
!                                      - ERRORS IN PROPAGATION CORRECTED

!     P.A.E.M. JANSSEN KNMI      1990  COUPLED MODEL.

!     H. GUNTHER     GKSS/ECMWF  8.91  LOGARITHMIC DEPTH TABLES.
!                                      MPI CYCLE_3 AND ECMWF VERSIONS
!                                      COMBINED INTO CYCLE_4.

!     J. BIDLOT ECMWF 1996   MESSAGE PASSING
!     J. DOYLE  ECMWF OCTOBER 1996   ATMOSPHERIC COUPLING
!     J. BIDLOT ECMWF FEBRUARY 1997   MODULE
!     J. BIDLOT ECMWF MARCH 1997  ADD SAVSTRESS AND SAVSPEC 
!     B. HANSEN ECMWF MARCH 1997  SIGNAL HANDLING.
!          LDSTOP* - SET .TRUE. IF STOP SIGNAL RECEIVED.
!          LDWRRE* - SET .TRUE. IF RESTART SIGNAL RECEIVED.
!     S. ABDALLA  ECMWF  OCTOBER 2000  INCLUDE AIR DENSITY & Zi/L
!                                      & CALL OUTBETA MOVED TO WAVEMDL
!

!*    PURPOSE.
!     --------

!       COMPUTATION OF THE 2-D FREQUENCY-DIRECTION WAVE SPECTRUM AT ALL
!       GRID POINTS FOR A GIVEN INITIAL SPECTRUM AND FORCING SURFACE
!       STRESS FIELD.

!**   INTERFACE.
!     ----------

!     *CALL* *WAMODEL (NADV, LDSTOP, LDWRRE, BLK2GLO,
!    &                 WVENVI, WVPRPT, FF_NOW, FF_NEXT, INTFLDS,
!    &                 WAM2NEMO, NEMO2WAM, FL1)
!        *NADV*      NUMBER OF ADVECTION ITERATIONS
!                    PER CALL OF WAMODEL, OUTPUT PARAMETER.
!        *LDSTOP*    SET .TRUE. IF STOP SIGNAL RECEIVED.
!        *LDWRRE*    SET .TRUE. IF RESTART SIGNAL RECEIVED.
!        *BLK2GLO*   BLOCK TO GRID TRANSFORMATION
!        *WVENVI*    WAVE ENVIRONMENT FIELDS
!        *WVPRPT*    WAVE PROPERTIES FIELDS
!        *FF_NOW*    FORCING FIELDS AT CURRENT TIME.
!        *FF_NEXT*   DATA STRUCTURE WITH THE NEXT FORCING FIELDS
!        *INTFLDS*   INTEGRATED/DERIVED PARAMETERS
!        *WAM2NEMO*  WAVE FIELDS PASSED TO NEMO
!        *NEMO2WAM*  FIELDS FRON OCEAN MODEL to WAM

!     REFERENCE.
!     ----------

!       SNYDER, R.L., F.W. DOBSON, J.A. ELLIOT, AND R.B. LONG:
!          ARRAY MEASUREMENTS OF ATMOSPHERIC PRESSURE FLUCTUATIONS
!          ABOVE SURFACE GRAVITY WAVES. J.FLUID MECH. 102, 1-59 ,1981.
!       G. KOMEN, S. HASSELMANN, K. HASSELMANN:
!          ON THE EXISTENCE OF A FULLY DEVELOPED WIND SEA SPECTRUM.
!          JPO,1984.
!       S. HASSELMANN, K. HASSELMANN, J.H. ALLENDER, T.P. BARNETT:
!          IMPROVED METHODS OF COMPUTING AND PARAMETERIZING THE
!          NONLINEAR ENERGY TRANSFER IN A GRAVITY WAVE SPECTRUM.
!          JPO, 1985.
!       S. HASSELMANN, K. HASSELMANN: A GLOBAL WAVE MODEL,
!          WAM REPORT,JUNE,30/1985.
!       P. JANSSEN, G. KOMEN: A SHALLOW WATER EXTENSION OF THE
!          3-G WAM-MODEL. WAM REPORT 1985.
!       THE WAMDI GROUP: THE WAM MODEL - A THIRD GENERATION OCEAN
!          WAVE PREDICTION MODEL. JPO, VOL. 18, NO. 12, 1988.
!       P.A.E.M JANSSEN: JPO, 1989 AND 1991.
!       K. HASSELMANN: TRANSPORT EQUATION OF FINITE DEPTH SURFACE
!          WAVE SPECTRUM IN TIME DEPENDANT CURRENT AND DEPTH FIELD USING
!          NONCANONICAL SPATIAL (SPHERICAL) AND WAVE NUMBER (FRQUENCY-
!          DIRECTION) COORDINATES. WAM REPORT 1988.

! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO, ENVIRONMENT, FREQUENCY, FORCING_FIELDS,  &
     &                         INTGT_PARAM_FIELDS, WAVE2OCEAN, OCEAN2WAVE

      USE YOWCPBO  , ONLY : IBOUNC   ,NBOUNC    ,GBOUNC  , IPOGBO  ,    &
     &            CBCPREF
      USE YOWCOUP  , ONLY : LWCOU    ,                                  &
     &                      LWNEMOCOU,                                  &
     &                      NEMOWSTEP, NEMOFRCO     ,                   &
     &                      NEMOCSTEP, NEMONSTEP
      USE YOWCOUT  , ONLY : COUTT    ,COUTS    ,FFLAG20  ,GFLAG20  ,    &
     &            FFLAG    ,GFLAG    ,                                  &
     &            NGOUT    ,LLOUTERS ,                                  &
     &            NIPRMOUT ,                                            &
     &            LFDB     ,NOUTT    ,NOUTS    ,                        &
     &            CASS     ,NASS     ,LOUTINT  ,                        &
     &            LRSTPARALW, LRSTINFDAT,                               &
     &            LRSTST0  ,LWAMANOUT
      USE YOWCURR  , ONLY : CDTCUR
      USE YOWFPBO  , ONLY : IBOUNF
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRIBHD, ONLY : LGRHDIFS 
      USE YOWGRID  , ONLY : IJS      ,IJL       ,                       &
     &            IJSLOC   ,IJLLOC   ,IJGLOBAL_OFFSET
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE
      USE YOWMESPAS, ONLY : LFDBIOOUT,LGRIBOUT ,LNOCDIN  ,LWAVEWIND 
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,KTAG 
      USE YOWPARAM , ONLY : NIBLO    ,NANG     ,NFRE
      USE YOWPCONS , ONLY : ZMISS    ,DEG      ,EPSMIN
      USE YOWSTAT  , ONLY : CDATEE   ,CDATEF   ,CDTPRO   ,CDTRES   ,    &
     &            CDATER   ,CDATES   ,CDTINTT  ,IDELPRO  ,IDELT    ,    &
     &            IDELWI   ,IREST    ,IDELRES  ,IDELINT  ,              &
     &            IASSI    ,                                            &
     &            CDTBC    ,IDELBC   ,                                  &
     &            IPROPAGS ,                                            &
     &            CDATEA   ,MARSTYPE ,                                  &
     &            LLSOURCE ,                                            &
     &            LANAONLY ,LFRSTFLD ,NPROMA_WAM,IREFDATE
      USE YOWSPEC, ONLY   : NBLKS    ,NBLKE
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH    ,CWI      ,LRESTARTED
      USE YOWUNIT  , ONLY : IU02     ,IU19     ,IU20
      USE YOWWAMI  , ONLY : CBPLTDT  ,CEPLTDT  ,IANALPD  ,IFOREPD  ,    &
     &            IDELWIN  ,NFCST    ,ISTAT
      USE YOWWIND  , ONLY : CDATEWO

      USE MPL_MODULE
      USE FDBSUBS_MOD, ONLY : IFLUSHFDBSUBS
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
      
! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "mpcrtbl.intfb.h"
#include "outwint.intfb.h"
#include "outwpsp.intfb.h"
#include "abort1.intfb.h"
#include "bouinpt.intfb.h"
#include "chesig.intfb.h"
#include "difdate.intfb.h"
#include "gsfile_new.intfb.h"
#include "headbc.intfb.h"
#include "iwam_get_unit.intfb.h"
#include "incdate.intfb.h"
#include "outbc.intfb.h"
#include "outbs.intfb.h"
#include "outint.intfb.h"
#include "outspec.intfb.h"
#include "outstep0.intfb.h"
#include "savspec.intfb.h"
#include "savstress.intfb.h"
#include "unsetice.intfb.h"
#include "updnemofields.intfb.h"
#include "updnemostress.intfb.h"
#include "wamintgr.intfb.h"
#include "writsta.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: NADV
      LOGICAL, INTENT(INOUT) :: LDSTOP, LDWRRE
      TYPE(WVGRIDGLO), DIMENSION(NIBLO), INTENT(IN) :: BLK2GLO
      TYPE(ENVIRONMENT), DIMENSION(IJS:IJL), INTENT(INOUT) :: WVENVI
      TYPE(FREQUENCY), DIMENSION(IJS:IJL,NFRE), INTENT(INOUT) :: WVPRPT
      TYPE(FORCING_FIELDS), DIMENSION(IJS:IJL), INTENT(INOUT) :: FF_NOW
      TYPE(FORCING_FIELDS), DIMENSION(IJS:IJL), INTENT(IN) :: FF_NEXT
      TYPE(INTGT_PARAM_FIELDS), DIMENSION(IJS:IJL), INTENT(INOUT) :: INTFLDS
      TYPE(WAVE2OCEAN), DIMENSION(IJS:IJL), INTENT(INOUT) :: WAM2NEMO
      TYPE(OCEAN2WAVE), DIMENSION(IJS:IJL), INTENT(IN) :: NEMO2WAM
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: FL1


      INTEGER(KIND=JWIM) :: IJ, K, M, J, IRA, KADV, ICH
      INTEGER(KIND=JWIM) :: IFIL, IC, ICL, ICR, II
      INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA
      INTEGER(KIND=JWIM) :: JSTPNEMO, IDATE, ITIME
      INTEGER(KIND=JWIM) :: IU04
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL) :: MIJ

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE) :: XLLWS

      CHARACTER(LEN= 2) :: MARSTYPEBAK
      CHARACTER(LEN=14) :: CDATEWH, CZERO
      CHARACTER(LEN=14) :: CDATE, CDTPRA, CDTIMP, CDTIMPNEXT, CDTRCF

      LOGICAL :: LLFLUSH
      LOGICAL :: LSV, LRST, LOUT
      LOGICAL :: LLNONASSI

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WAMODEL',0,ZHOOK_HANDLE)

ASSOCIATE(WSWAVE => FF_NOW%WSWAVE, &
 &        WDWAVE => FF_NOW%WDWAVE, &
 &        CICOVER => FF_NOW%CICOVER, &
 &        USTOKES => INTFLDS%USTOKES, &
 &        VSTOKES => INTFLDS%VSTOKES, &
 &        STRNMS  => INTFLDS%STRNMS )


!     0.0 INITIALISATION
!         --------------

      CZERO = ' '
      LLFLUSH = .FALSE.
      KTAG=200
      LRSTST0=.FALSE.
      NPROMA=NPROMA_WAM

!     TIME FOR THE NEXT SOURCE TERM INTEGRATION
      CDTIMPNEXT=CDTPRO
      CALL INCDATE(CDTIMPNEXT,IDELT)   
!     TIME FOR WIND INPUT UPDATE (SEE NEWWIND)
      CDTIMP=CDTPRO

!     0.1 MINIMUM ENERGY
!         --------------
      IF (CDTPRO == CDATEA .AND. LLSOURCE ) THEN
!       INSURE THERE IS SOME WAVE ENERGY FOR GRID POINTS THAT HAVE BEEN
!       FREED FROM SEA ICE (ONLY DONE INITIALLY AND IF THE MODEL IS NOT
!       RESTARTED).
!       IT ALSO RESETS THE MIMIMUM ENERGY LEVEL THAT MIGHT HAVE BEEN LOST
!       WHEN GETTING THE DATA FROM GRIB.
        CALL GSTATS(1236,0)
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO=IJS,IJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL)
          CALL UNSETICE(KIJS, KIJL, WVENVI(KIJS), FF_NOW(KIJS), FL1(KIJS:KIJL,:,:))  
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1236,1)
      ENDIF


!     0.2 OUTPUT INITIAL CONDITION AND/OR FORECAST STEP 0
!         -----------------------------------------------
      IF (CDTPRO == CDATEA .OR. CDTPRO == CDATEF) THEN
         CALL OUTSTEP0 (WVENVI, WVPRPT, FF_NOW, INTFLDS,  &
 &                      WAM2NEMO, NEMO2WAM, FL1)
      ENDIF




!*    1. ADVECTION/PHYSICS TIME LOOP.
!        ----------------------------

      ADVECTION : DO KADV = 1,NADV

!*      1.1 FIX END DATE OF THIS PROPAGATION STEP AND OUTPUT TIMES.
!           -------------------------------------------------------

        CDTPRA = CDTPRO
        CALL INCDATE(CDTPRO,IDELPRO)

!       UPDATE OUTPUT TIMES.
        IF (NOUTT > 0) THEN
          CDTINTT = CZERO
          DO J=1,NOUTT
            IF (CDTPRO == COUTT(J)) THEN
              IF (FFLAG20 .OR. GFLAG20) CDTINTT = COUTT(J)
            ENDIF
          ENDDO
        ELSE
          IF ((FFLAG20.OR.GFLAG20) .AND. CDTINTT.LT.CDTPRO) CALL INCDATE (CDTINTT,IDELINT)
        ENDIF

!       UPDATE SPECTRA OUTPUT DATE
        IF (NOUTS > 0) THEN
!         reset CDATES to insure that spectra output is only controlled
!         by list COUTS
          CDATES='000000000000'
          DO J=1,NOUTS
            IF (CDTPRO == COUTS(J)) THEN
              CDTRES=CDTPRO
              CDATER=CDTRES
              CDATES=CDTRES
              EXIT
            ENDIF
          ENDDO
        ELSE
          IF (CDTRES < CDTPRO) CALL INCDATE(CDTRES,IDELRES)
        ENDIF

!NEST (not used at ECMWF)
        IF ((IBOUNC == 1 .OR. IBOUNF == 1) .AND. CDTBC < CDTPRO) CALL INCDATE(CDTBC,IDELBC)
!NEST



!            -------------------------------------------------------
!*      1.2  THIS IS THE CORE OF THE WAVE MODEL:
!*           COMPUTATION OF PROPAGATION
!*           INTEGRATION OF SOURCE TERMS OVER SUB TIME STEPS BETWEEN
!*           PROPAGATION TIME STEPS.
!            -------------------------------------------------------
!        SET TIME COUNTER.
         CDATE   = CDTPRA
         CDATEWH = CDATEWO

         DO WHILE (CDTIMPNEXT <= CDTPRO)
           CALL WAMINTGR (IJS, IJL, BLK2GLO,                          &
 &                        CDTPRA, CDATE, CDATEWH, CDTIMP, CDTIMPNEXT, &
 &                        WVENVI, WVPRPT, FF_NOW, FF_NEXT, INTFLDS,   &
 &                        WAM2NEMO, MIJ, FL1, XLLWS)
         ENDDO



#ifdef ECMWF
        IF (.NOT.LWCOU .AND. .NOT. LDSTOP) THEN
!!!!      the call to CHESIG is a signal handeling facility which is
!!!!      specific to running WAM at ECMWF, it can be ignored when
!!!!      WAM is not run at ECMWF.
            CALL CHESIG (IU06, ITEST, IRANK, NPROC, LDSTOP, LDWRRE)
        ENDIF
#endif

!       1.3 CHECK WHETHER OUTPUT(s) NEEDED
!           ------------------------------
        LRST=(LDWRRE .AND. KADV == NADV )
        IF (LRST) THEN
          WRITE(IU06,*) ' '
          WRITE(IU06,*) '  ******************************************'
          IF (LDSTOP) THEN
          WRITE(IU06,*) '  AN INTERRUPT SIGNAL HAS BEEN RECEIVED '
          ENDIF
          WRITE(IU06,*) '  THE NECESSARY BINARY RESTART FILES WILL BE'
          WRITE(IU06,*) '  GENERATED.'
          WRITE(IU06,*) '  ******************************************'
          WRITE(IU06,*) ' '
          CALL FLUSH (IU06)
        ENDIF


!NEST (not used at ECMWF)
!*      1.4.1 INPUT OF BOUNDARY VALUES.
!           -------------------------
        IF (IBOUNF == 1) CALL BOUINPT (IU02, FL1, IJS, IJL, NBLKS, NBLKE)
!*      1.4.2 OUTPUT OF BOUNDARY POINTS.
!           --------------------------
        IF (IBOUNC == 1) CALL OUTBC (FL1, IJS, IJL, IU19) 
!NEST


!*      1.5 POINT OUTPUT (not usually used at ECMWF) 
!           ----------------------------------------
        IF ( NGOUT > 0 .AND. (CDTINTT == CDTPRO .OR. LRST) ) THEN
!           OUTPUT POINT SPECTRA (not usually used at ECMWF)
            IF (LLOUTERS) CALL OUTWPSP (IJS, IJL, FL1, FF_NOW)
        ENDIF


!       1.6 COMPUTE OUTPUT PARAMETERS FIELDS AND PRINT OUT NORMS
!           ----------------------------------------------------
        IF ( (CDTINTT == CDTPRO .OR. LRST) .AND. NIPRMOUT > 0 ) THEN
            CALL OUTBS (IJS, IJL, MIJ, FL1, XLLWS,                   &
     &                  WVPRPT, WVENVI, FF_NOW, INTFLDS, NEMO2WAM)
        ENDIF


!*      1.7 ONE PROPAGATION TIMESTEP DONE
!           -----------------------------
        WRITE(IU06,*) ' !!!!!!!!!!!!!! WAVE FIELDS INTEGRATED FOR DATE : ', CDTPRO

!       IS THIS AN ANALYSIS TIME FROM THE INPUT LIST ?
        LLNONASSI=.TRUE.
        IF (IASSI == 1) THEN
          DO J=1,NASS
            IF (CDTPRO == CASS(J)) THEN
              LLNONASSI=.FALSE.
              EXIT
            ENDIF
          ENDDO
        ENDIF


!*      1.8 SAVE FIELD OF 2D SPECTRA AND/OR BINARY RESTART FILES.
!           -----------------------------------------------------
!         THIS WILL HAPPEN WHEN IT IS NOT IN DATA ASSIMILATION MODE AND
!         IT IS EITHER A DETERMINED OUTPUT TIME
!         OR THE INTERUPT SIGNAL HAS BEEN TRIGGERED and it will wait
!         until the end of the advection loop.
!         OTHERWISE THE OUTPUT WILL OCCUR IN WAMASSI.

        LSV=(CDTRES == CDTPRO .OR. CDATEE == CDTPRO .OR. CDTPRO == CDATER)

        IF (LSV .OR. LRST) THEN

          LOUT = ((IREST == 1) .AND. (CDTPRO == CDATER .OR. CDTPRO <= CDATES)) .AND. LSV .AND. LWAMANOUT

          IF ( LOUT .OR. LRST ) THEN

!           1.8.1 SAVE SPECTRUM IN GRIB
!                 ---------------------
            IF (LOUT .AND. LGRIBOUT) THEN
!             we have insured that the spectra will be written to FDB
!             even when the restart option is triggered and it is an
!             output step for the spectra.

!             IF THE OUTPUT TIME IS NOT AN ANALYSIS TIME THEN TYPE FG or 4V
!             BECOMES TYPE AN (i.e. pseudo analysis)
              MARSTYPEBAK=MARSTYPE
              IF ((MARSTYPE == 'fg' .AND. KADV < NADV) .OR.             &
     &            (MARSTYPE == '4v' .AND. LLNONASSI) ) THEN
                MARSTYPE='an'
              ENDIF

              CALL OUTSPEC(IJS, IJL, FL1, CICOVER)
              LLFLUSH = .TRUE.

              MARSTYPE=MARSTYPEBAK

              WRITE(IU06,*) ' '
              WRITE(IU06,*) '  GRIB WAVE SPECTRA DISPOSED AT........ CDTPRO  = ', CDTPRO
              WRITE(IU06,*) ' '
            ENDIF  

!           1.8.2 SAVE RESTART FILES IN PURE BINARY FORM (in needed)
!                 --------------------------------------
            IF ( .NOT.LGRIBOUT .OR. LDWRRE ) THEN
              CALL SAVSTRESS(IJS, IJL, WVENVI, FF_NOW, NBLKS, NBLKE, CDTPRO, CDATEF)
              WRITE(IU06,*) ' '
              WRITE(IU06,*) '  BINARY STRESS FILE DISPOSED AT........ CDTPRO  = ', CDTPRO
              WRITE(IU06,*) ' '

              CALL SAVSPEC(IJS, IJL, FL1, NBLKS, NBLKE, CDTPRO, CDATEF, CDATER)
              WRITE(IU06,*) '  BINARY WAVE SPECTRA DISPOSED AT........ CDTPRO  = ', CDTPRO
              WRITE(IU06,*) ' '
              CALL FLUSH(IU06)
            ENDIF


!*          1.8.3 UPDATE, WRITE AND SAVE WAMINFO FILE.
!                 -----------------------------------
            IF (LRST .AND. IRANK == 1) THEN
              ICH = 7 
              CALL DIFDATE (CDATEF,CDATEE,IFOREPD)
              IF (CDTPRO <= CDATEF) THEN
                CALL DIFDATE (CDTPRO,CDATEF,IANALPD)
                CBPLTDT = CDTPRO
                NFCST = 1
              ELSE
                NFCST = 0
                IANALPD = 0
                CBPLTDT = CDATEF
                CALL DIFDATE (CDTPRO,CDATEE,IFOREPD)
              ENDIF
              ISTAT(:) = 0
              IF (CDATE == CDATEE) ISTAT(1) = 1
              IDELWIN = IDELWI

              CEPLTDT = CDATEF

              IU04 = IWAM_GET_UNIT (IU06,CWI(1:ICPLEN+8) , 'w', 'f', 0)

              CALL WRITSTA (IU04, CDTPRO, CDATEE, IANALPD, IFOREPD,     &
     &                      IDELWIN, CDATER, CDATES, CBPLTDT, CEPLTDT,  &
     &                      IASSI, NFCST, ISTAT, CDTCUR,                &
     &                      LRSTPARALW, NPROC)
 
              CLOSE (IU04)
              WRITE(IU06,*) ' WAMINFO FILE WRITTEN FOR RESTART...',     &
     &         ' CDTPRO  = ', CDTPRO
              WRITE(IU06,*) '                                    ',     &
     &         ' CDATEF  = ', CDATEF
              WRITE(IU06,*) ' TO ', CWI(1:ICPLEN+8)
              CALL FLUSH(IU06)

              IF (LRSTINFDAT) THEN
!               WRITE AN ADDITIONAL wamfile WITH DATE/TIME INFO added to filename
                CDTRCF=CDTPRO
                CALL INCDATE(CDTRCF,IDELPRO)
                IU04 =  IWAM_GET_UNIT (IU06,CWI(1:ICPLEN+8)//'.'//      &
     &                              CDTRCF(1:8)//'_'//CDTRCF(9:14),     &
     &                              'w', 'f', 0)
                CALL WRITSTA (IU04, CDTPRO, CDATEE, IANALPD, IFOREPD,   &
     &                        IDELWIN, CDATER, CDATES, CBPLTDT, CEPLTDT,&
     &                        IASSI, NFCST, ISTAT, CDTCUR,              &
     &                        LRSTPARALW, NPROC)
                CLOSE (IU04)
              ENDIF

            ENDIF  ! waminfo

          ENDIF  ! END SAVE RESTART FIELDS

        ENDIF  ! END SAVE 1.8


!       1.9 WRITE INTEGRATED PARAMETER DATA TO FDB OR TO FILE
!           -------------------------------------------------
        IF (LWAMANOUT .AND. CDTINTT == CDTPRO .AND. NIPRMOUT > 0 ) THEN
!         IF THE OUTPUT TIME IS NOT AN ANALYSIS TIME THEN TYPE FG or 4V
!         BECOMES TYPE AN (i.e. speudo analysis)
          MARSTYPEBAK=MARSTYPE
          IF ((MARSTYPE == 'fg' .AND. KADV < NADV) .OR.                 &
     &        (MARSTYPE == '4v' .AND. LLNONASSI) ) THEN
            MARSTYPE='an'
          ENDIF

          CALL OUTWINT
          LLFLUSH = .TRUE.

          MARSTYPE=MARSTYPEBAK

          CALL GSTATS(753,0)
          CALL MPL_BARRIER(CDSTRING='WAMODEL:')
          CALL GSTATS(753,1)
        ENDIF


!*      1.10 FLUSH FDB IF IT HAS BEEN USED AND IT IS NOT AN ANALYSIS (it will be done in *wamassi*)
!            -------------------------------------------------------

        IF ( LFDB .AND. LLFLUSH .AND. (IASSI /= 1 .OR. CDTPRO > CDATEF) ) THEN
          CALL GSTATS(1976,0)
          CALL IFLUSHFDBSUBS()
          CALL GSTATS(1976,1)
          WRITE(IU06,*) ' ' 
          WRITE(IU06,*) '  FDB FLUSHED AT ',  CDTPRO, ' FROM WAMODEL. '
          CALL FLUSH (IU06)
          LLFLUSH=.FALSE.
        ENDIF


!       1.11 FOR PURE BINARY DATA (obsolete option at ECMWF !!!!!!!):
!*           OUTPUT FILES AND RECOVERY FILES ARE DISPOSED WHEN
!            TIME REACHES THE DISPOSE DATE OR WHEN THE MODEL
!            HAS BEEN SIGNALLED TO DO SO.
!            -------------------------------------------------
        IF (FFLAG20) THEN
          IF (CDATEE == CDTPRO .AND. LOUTINT .AND. IRANK == 1  ) THEN
            CALL GSFILE (IU06, IU20, 0, CDTPRO, CDATEF, 'MAP', 'S')
          ENDIF
        ENDIF

!NEST (not used at ECMWF)
!       SAVE BOUNDARY VALUE FILE.
        IF (CDTBC == CDTPRO) THEN
          IF (IBOUNC == 1 .AND. IRANK == 1 ) THEN
            DO II=1,GBOUNC
            CALL GSFILE(IU06, IU19(II), 0, CDTBC, CDTBC,                &
     &        CBCPREF(II), 'S')
            IF (CDTBC < CDATEE)                                         &
     &        CALL HEADBC (IPOGBO(II)-IPOGBO(II-1), IDELPRO,            &
     &                     TH(1), FR(1), IU19(II), IU06) 
            ENDDO
          ENDIF
        ENDIF
!NEST
 

!*      1.12 WAM-NEMO COUPLING (!!!!! WHEN NO atmospheric model !!!!!!) (currently not used at ECMWF)
!       (when coupled see cnt4 in ifs)
!       ----------------------------------------------------
        IF (LWNEMOCOU .AND. (.NOT.LWCOU)) THEN
          NEMOWSTEP=NEMOWSTEP+1

          IF (MOD(NEMOWSTEP,NEMOFRCO) == 0) THEN
            CALL UPDNEMOFIELDS
            CALL UPDNEMOSTRESS

#ifdef WITH_NEMO
            DO JSTPNEMO=NEMOCSTEP,NEMOCSTEP+NEMONSTEP-1
               ! Advance the NEMO model 1 time step
               CALL NEMOGCMCOUP_STEP( JSTPNEMO, IDATE, ITIME )
               WRITE(IU06,*)'NEMO TIME IS : ',JSTPNEMO, IDATE, ITIME
            ENDDO
#endif
            NEMOCSTEP=NEMOCSTEP+NEMONSTEP
          ENDIF
        ENDIF


!*    BRANCHING BACK TO 1.0 FOR NEXT PROPAGATION STEP.
      ENDDO ADVECTION

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('WAMODEL',1,ZHOOK_HANDLE)

END SUBROUTINE WAMODEL
