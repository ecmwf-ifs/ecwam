      SUBROUTINE WAMODEL (NADV, BETAG, NLONW, NLATW, LDSTOP, LDWRRE)

C ----------------------------------------------------------------------
C
C**** *WAMODEL* - 3-G WAM MODEL - TIME INTEGRATION OF WAVE FIELDS.
C
C     S.D. HASSELMANN  MPI       1.12.85
C
C     G. KOMEN         KNMI         6.86  MODIFIED FOR SHALLOW WATER
C     P. JANSSEN                          ASPECTS.
C
C     S.D. HASSELMANN  MPI       15.2.87  MODIFIED FOR CYBER 205.
C
C     P. LIONELLO      ISDGM      6.3.87  MODIFIED TO OUTPUT SWELL.
C
C     S.D. HASSELMANN  MPI        1.6.87  ALL VERSIONS COMBINED INTO
C                                         ONE MODEL. DEEP AND SHALLOW
C                                         WATER , CRAY AND CYBER 205
C                                         VERSION.
C
C     CYCLE_2 MODICIFATIONS:
C     ----------------------
C
C     L. ZAMBRESKY     GKSS        10.87  OPTIMIZED FOR CRAY, CYBER 205
C     H. GUNTHER
C
C     A. SPEIDEL       MPI          4.88  VARIABLE DIMENSIONS, INTERNAL
C                                         CHECKS (CFL-CRITERION).
C
C     A. SPEIDEL       MPI         11.88  CHANGES FOR CRAY-2.
C
C     K. HUBBERT       POL          6.89  DEPTH AND CURRENT REFRACTION.
C                                         PRECALCULATION OF TERMS IN
C                                         *PROPDOT*.
C                                         SOLVE WAVE ACTION EQUATION
C                                         FOR CURRENT REFRACTION.
C
C     CYCLE_3 MODICIFATIONS:
C     ----------------------
C
C     R. PORTZ , S.D. HASSELMANN   MPI          1990
C
C      - RESTRUCTURE MODEL TO CALL THE ACTUAL INTEGRATION IN TIME
C        AS A SUBROUTINE: WAMODEL. A SHELL PROGRAM "WAMSHELL" READS
C        OUTPUT FROM PREPROC AND COMPUTES THE WIND ARRAYS FOR THE
C        INTEGRATION PERIOD FROM PREWIND, WHICH HAS BEEN INCORPORATED
C        AS A SUBROUTINE.
C      - ALL INTERMEDIATE AND RESTART I/O IS DONE IN THE SUBROUTINE
C        WAMODEL AND INPREST.
C      - THE COMMON BLOCK IN THE PREPROCESSOR AND MODEL ARE MADE
C        COMPATIBLE.
C      - THE COMPUTATION OF SEVERAL PARAMETERS HAS BEEN TRANSFERRED
C        FROM THE MODEL TO PREPROC.
C      - DEPTH AND CURRENT REFRACTION HAS BEEN INCORPORATED INTO THE
C        MODEL.
C      - OPEN BOUNDARIES ARE INCORPORATED IN THE MODEL.
C      - SEVERAL MINOR ERRORS HAVE BEEN REMOVED.
C      - THE BUFFERED I/O FOR THE CYBER 205 HAS BEEN CHANGED INTO A
C        BINARY READ AND WRITE.
C
C     CYCLE_4 MODICIFATIONS:
C     ----------------------
C
C     L. ZAMBRESKY   GKSS/ECMWF   6.89  ECMWF SUB VERSION
C                                       BASED ON CYCLE_2.
C
C     H. GUNTHER     GKSS/ECMWF 10.89  ECMWF SUB VERSION REORGANIZED.
C                                      - COMMON BLOCK STRUCTURE.
C                                      - BLOCKING STRUCTURE.
C                                      - TIME COUNTING.
C                                      - GRIDDED OUTPUT FIELDS.
C                                      - HEADERS ADDED TO OUTPUT FILES.
C                                      - ERRORS IN PROPAGATION CORRECTED
C
C     P.A.E.M. JANSSEN KNMI      1990  COUPLED MODEL.
C
C     H. GUNTHER     GKSS/ECMWF  8.91  LOGARITHMIC DEPTH TABLES.
C                                      MPI CYCLE_3 AND ECMWF VERSIONS
C                                      COMBINED INTO CYCLE_4.
C
C     J. BIDLOT ECMWF 1996   MESSAGE PASSING
C     J. DOYLE  ECMWF OCTOBER 1996   ATMOSPHERIC COUPLING
C     J. BIDLOT ECMWF FEBRUARY 1997   MODULE
C     J. BIDLOT ECMWF MARCH 1997  ADD SAVSTRESS AND SAVSPEC 
C     B. HANSEN ECMWF MARCH 1997  SIGNAL HANDLING.
C          LDSTOP* - SET .TRUE. IF STOP SIGNAL RECEIVED.
C          LDWRRE* - SET .TRUE. IF RESTART SIGNAL RECEIVED.
C
CSHALLOW
C          DIFFERENCES FOR SHALLOW WATER RUNS TO DEEP WATER RUNS
C          ARE ENCLOSED IN COMMENT LINES : 'CSHALLOW'.
CSHALLOW
CNEST
C          DIFFERENCES FOR NESTED GRID RUNS TO NORMAL RUNS
C          ARE ENCLOSED IN COMMENT LINES : 'CNEST'.
CNEST
CREFRA
C          DIFFERENCES FOR REFRACTION RUNS TO NORMAL RUNS
C          ARE ENCLOSED IN COMMENT LINES : 'CREFRA'.
CREFRA
C
C*    PURPOSE.
C     --------
C
C       COMPUTATION OF THE 2-D FREQUENCY-DIRECTION WAVE SPECTRUM AT ALL
C       GRID POINTS FOR A GIVEN INITIAL SPECTRUM AND FORCING SURFACE
C       STRESS FIELD.
C
C**   INTERFACE.
C     ----------
C
C     *CALL* *WAMODEL (NADV,BETAG,NLONW,NLATW)
C        *NADV*      INTEGER   NUMBER OF ADVECTION ITERATIONS
C                              PER CALL OF WAMODEL, OUTPUT PARAMETER.
C        *BETAG*     REAL      CHARNOCK PARAMETER TO BE PASSED BACK TO
C                              THE ATMOSPHERIC MODEL
C        *NLONW*     INTEGER   FIRST DIMENSION OF BETAG
C        *NLATW*     INTEGER   SECOND DIMENSION OF BETAG
C
C     METHOD.
C     -------
C
C       GRID POINTS ARE LAT - LONG,VECTORIZATION IS ACHIEVED BY RUNNING
C       THROUGH THE GRID POINTS IN AN INNER LOOP ORGANIZED AS 1-D ARRAY
C       IN BLOCKS,-ALL COMPUTATIONS ARE CARRIED OUT FOR ONE BLOCK AT A
C       TIME (SEE "BLOCK STRUCTURE" BELOW)
C
C       ALL COMPONENTS OF THE SPECTRUM ARE COMPUTED PROGNOSTICALLY FROM
C       THE SPECTRAL TRANSPORT EQUATION UP TO A VARIABLE CUT-OFF
C       FREQUENCY = MAX(4*FPM,2.5*FMEAN),WHERE FPM IS THE
C       PIERSON MOSKOVITZ FREQUENCY AND FMEAN IS THE MEAN FREQUENCY,
C       BEYOND THE PROGNOSTIC CUTOFF A DIAGNOSTIC F**-5 TAIL IS ATTACHED
C       CONTINUOUSLY FOR EACH DIRECTION,
C
C       SOURCE FUNCTIONS ARE TAKEN FROM KOMEN ET AL(1984)
C
C       THE NONLINEAR TRANSFER IS PARAMETERIZED BY THE DISCRETE INTER-
C       ACTION APPROXIMATION OF HASSELMANN ET AL (1985B)
C
C       THE SOURCE FUNCTION AND THE ADVECTION TERM ARE INTEGRATED ON TWO
C       DIFFERENT TIME STEP LEVELS AND WITH DIFFERENT METHODS,-THE
C       ADVECTION TIME STEP IS A MULTIPLE OF THE SOURCE FUNCTION
C       TIME STEP.
C
C       THE SOURCE FUNCTIONS ARE INTEGRATED IMPLICITLY ACCORDING TO
C       HASSELMANN AND HASSELMANN (1985A),-THE RELEVANT FUNCTIONAL
C       DERIVATIVES OF THE INDIVIDUAL SOURCE FUNCTIONS REQUIRED FOR THE
C       SOLUTION OF THE IMPLICIT EQUATION ARE COMPUTED WITHIN THE SOURCE
C       FUNCTION SUBS,- THE TIME STEP IS TYPICALLY 20 MIN,
C
C       THE ADVECTION IS INTEGRATED BY A FIRST ORDER UPWIND SCHEME,ALSO
C       ACCORDING TO HASSELMANN AND HASSELMANN (1985A),-THE ADVECTIVE
C       TIMESTEP IS DEPENDENT ON THE FREQUENCY AND SPATIAL GRID IN
C       ACCORDANCE WITH CFL,
C
C       WINDS ARE READ IN EVERY WIND TIME STEP.IF THE WIND TIME STEP IS
C       GREATER THAN THE SOURCE TERM TIME STEP DELTWIND/DELTSOURCE STEPS
C       ARE INTEGRATED WITH CONSTANT WINDS,
C       WIND TIME STEP,PROPAGATION TIME STEP AND SOURCE TERM TIME STEP
C       SHOULD HAVE INTEGER RATIOS, THEY ARE GIVEN IN SECONDS AT
C       FULL MINUTES.
C
CNEST
C       ZERO ENERGY INFLUX IS ASSUMED AT COAST LINES. OPEN BOUNDARIES
C       ARE INCORPORATED IN THE MODEL, IF IT RUNS AS A NESTED GRID.
CNEST
C
C       BLOCK STRUCTURE (SEE PREPROC FOR DETAILS):
C       SEA POINTS ARE COLLECTED INTO A 1-DIMENSIONAL ARRAY.
C       BLOCKS OF MAXIMALLY NIBLO ELEMENTS.
C       SEA POINTS ARE COUNTED ALONG LINES OF LATITUDES FROM LEFT COAST
C       TO RIGHT COAST WORKING FROM SOUTH TO NORTH.
C       BLOCKS OVERLAP OVER TWO LATITUDE LINES,TO COMPUTE NORTH-SOUTH
C       ADVECTION TERMS, SEE ALSO COMMON GRIDPAR AND UBUF.
C
C       THE WIND FILES FOR THE BLOCKED WINDS CREATED BY PREWIND ARE
C       READ AND DELETED IN SUB IMPLSCH (IU17 AND IU18). THE FILE
c       NAMES ARE CREATED IN SUB CREWFN AND AN IMPLICIT OPEN IS USED.
C
C       THE FILE HANDLING SUBS OPENFIL, GSFILE AND CREWFN ARE COMPUTER
C       DEPENDENT AND MAY BE ADOPTED BY THE USER.
C       THE PROGRAM CLOSES AND DELETES ALL WORK FILES.
C
C       ALL PARAMETERS HAVE TO BE THE VALUES GIVEN AT THE END OF THE
C       PREPROC OUTPUT IN COLUMN 'REQUIRED'.
C
C     EXTERNALS.
C     ----------
C
C       *ABORT1*    - TERMINATES PROCESSING.
C       *AIRSEA*    - SURFACE LAYER STRESS.
CNEST
C       *BOUINPT*   - BOUNDARY VALUE INPUT.
CNEST
C       *CREWFN*    - CREATES A WIND FILE NAME.
CREFRA
C       *DOTDC*     - READ COMMON REFDOT.
CREFRA
C       *FEMEAN*    - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT.
C       *FILLBL*    - ADD LATITUDES TO A BLOCK.
C       *GSFILE*    - ROUTINE TO DYNAMICALLY FETCH OR DISPOSE FILES.
C       *HEADBC*    - WRITE BOUNDARY OUTPUT FILE HEADER.
C       *iecf_len*  - LENGTH OF A CHARACTER ARRAY.
C       *IMPLSCH*   - IMPLICIT SCHEME FOR INTEGRATION OF SOURCE
C                     FUNCTIONS IN TIME AND INPUT OF WINDS.
C       *INCDATE*   - UPDATE DATE TIME GROUP.
CREFRA
C       *INTPOL*    - MAP SPECTRUM FROM SIGMA TO OMEGA SPACE.
CREFRA
CNEST
C       *INTSPEC*   - INTERPOLATION OF SPECTRA.
CNEST
C       *MAKEGRID*  - MAKES GRIDDED FIELDS.
C       *MPEXCHNG*  - PERFORMED MESSAGE PASSING OF THE SPECTRUM
C       *OPENFIL*   - OPENS A FILE.
C       *OUTBC*     - OUTPUT OF BOUNDARY VALUES.
C       *OUTBS*     - CONTROLS OUTPUT FROM BLOCKS.
C       *OUTGRID*   - SAVE BLOCKED PARAMETERS INTO GRID ARRAYS.
C       *OUTINT*    - OUTPUT OF INTEGRATED PARAMETERS.
C       *OUTPP*     - WRITES GRIDDED FIELDS TO PRINTER.
C       *OUTSPP*    - OUTPUT OF SPECTRA AT SELECTED POINTS.
C       *PEAKFR*    - COMUTE PEAK FREQUENCY.
C       *PROPAG*    - PROPAGATION SCHEME.
C       *PRSPP*     - PRINT A SPECTRUM.
C       *PRSPPS*    - PRINT A SWELL SPECTRUM.
CNEST
C       *ROTSPEC*   - ROTATE A SPECTRUM.
CNEST
C       *SAVSTRESS* - DISPOSE STRESS/WIND RESTART FILES.
C       *SAVSPEC    - DISPOSE SPECTRUM RESTART FILES.
CSHALLOW
C       *SBOTTOM*   - COMPUTES BOTTOM DISSIPATION SOURCE TERM AND
C                     LINEAR CONTRIBUTION TO FUNCTIONAL MATRIX.
CSHALLOW
C       *SDISSIP*   - COMPUTATION OF DISSIPATION SOURCE FUNCTION
C                     AND LINEAR CONTRIBUTION OF DISSIPATION TO
C                     FUNCTIONAL MATRIX IN IMPLICIT SCHEME.
C       *SEMEAN*    - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.
C       *SEPWISW*   - COMPUTATION OF 2-DIMENSIONAL SWELL DISTRIBUTION
C                     TOTAL SWELL ENERGY, MEAN SWELL DIRECTION, AND
C                     MEAN SWELL FREQUENCY AT EACH GRID POINT.
CICE
C       *SETICE*    - SET SPECTRA ON ICE EDGE TO ZERO.
CICE
C       *SINPUT*    - COMPUTATION OF INPUT SOURCE FUNCTION, AND
C                     LINEAR CONTRIBUTION OF INPUT SOURCE FUNCTION
C                     TO FUNCTIONAL MATRIX IN IMPLICIT SCHEME.
C       *SNONLIN*   - COMPUTATION OF NONLINEAR TRANSFER RATE AND
C                     DIAGONAL LINEAR CONTRIBUTION OF NONLINEAR SOURCE
C                     FUNCTION TO FUNCTIONAL MATRIX.
C       *SPLITBL*   - STORING OVERLAPPING LATITUDES OF BLOCKS.
C       *STHQ*      - COMPUTATION OF MEAN WAVE DIRECTION AT EACH
C                     GRID POINT.
C       *STRESSO*   - COMPUTATION OF WAVE STRESS.
CNEST
C       *STRSPEC*   - STRETCH A SPECTRUM.
CNEST
C       *WAMNORM*   - COMPUTES A FEW NORMS OF GRIDDED FIELDS
C
C     REFERENCE.
C     ----------
C
C       SNYDER, R.L., F.W. DOBSON, J.A. ELLIOT, AND R.B. LONG:
C          ARRAY MEASUREMENTS OF ATMOSPHERIC PRESSURE FLUCTUATIONS
C          ABOVE SURFACE GRAVITY WAVES. J.FLUID MECH. 102, 1-59 ,1981.
C       G. KOMEN, S. HASSELMANN, K. HASSELMANN:
C          ON THE EXISTENCE OF A FULLY DEVELOPED WIND SEA SPECTRUM.
C          JPO,1984.
C       S. HASSELMANN, K. HASSELMANN, J.H. ALLENDER, T.P. BARNETT:
C          IMPROVED METHODS OF COMPUTING AND PARAMETERIZING THE
C          NONLINEAR ENERGY TRANSFER IN A GRAVITY WAVE SPECTRUM.
C          JPO, 1985.
C       S. HASSELMANN, K. HASSELMANN: A GLOBAL WAVE MODEL,
C          WAM REPORT,JUNE,30/1985.
C       P. JANSSEN, G. KOMEN: A SHALLOW WATER EXTENSION OF THE
C          3-G WAM-MODEL. WAM REPORT 1985.
C       THE WAMDI GROUP: THE WAM MODEL - A THIRD GENERATION OCEAN
C          WAVE PREDICTION MODEL. JPO, VOL. 18, NO. 12, 1988.
C       P.A.E.M JANSSEN: JPO, 1989 AND 1991.
C       K. HASSELMANN: TRANSPORT EQUATION OF FINITE DEPTH SURFACE
C          WAVE SPECTRUM IN TIME DPENDANT CURRENT AND DEPTH FIELD USING
C          NONCANONICAL SPACIAL (SPHERICAL) AND WAVE NUMBER (FRQUENCY-
C          DIRECTION) COORDINATES. WAM REPROT 1988.
C
C -------------------------------------------------------------------
      USE WAVE_MODULE
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "txtpara.h"
C
CNEST
#include "comcbou.h"
C
#include "txtcbou.h"
CNEST
C
#include "comcoup.h"
C
#include "txtcoup.h"
C
#include "comcout.h"
C
#include "txtcout.h"
C
CREFRA
#include "comcurr.h"
C
#include "txtcurr.h"
CREFRA
C
CNEST
#include "comfbou.h"
C
#include "txtfbou.h"
CNEST
C
#include "comfred.h"
C
#include "txtfred.h"
C
#include "comgrid.h"
C
#include "txtgrid.h"
C
CICE
#include "comice.h"
C
#include "txtice.h"
CICE
C
#include "comindn.h"
C
#include "txtindn.h"
C
#include "comintp.h"
C
#include "txtintp.h"
C
#include "comints.h"
C
#include "txtints.h"
C
#include "commap.h"
C
#include "txtmap.h"
C
#include "commean.h"
C
#include "txtmean.h"
C
#include "comover.h"
C
#include "txtover.h"
C
CREFRA
#include "comrefd.h"
C
#include "txtrefd.h"
CREFRA
C
#include "comshal.h"
C
#include "txtshal.h"
C
#include "comstat.h"
C
#include "txtstat.h"
C
#include "comswel.h"
C
#include "txtswel.h"
C
#include "comtabl.h"
C
#include "txttabl.h"
C
#include "comtest.h"
C
#include "txttest.h"
C
#include "comtext.h"
C
#include "txttext.h"
C
#include "comubuf.h"
C
#include "txtubuf.h"
C
#include "comunit.h"
C
#include "txtunit.h"
C
#include "comwind.h"
C
#include "txtwind.h"
#include "parcons.h"
C
!
#include "commesspass.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
#include "comwami.h"
C
      CHARACTER*10 CDATE, CDATEWH, CDTPRA, ZERO
      CHARACTER*10 CDTINTTBAK,CDTINTSBAK
      CHARACTER*80 FILENAME
      LOGICAL LOPENED, LDSTOP, LDWRRE, LLFLUSH
C
C ----------------------------------------------------------------------
C     define local allocatable arrays

      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: FL4
C ----------------------------------------------------------------------
C
#include "commltsk.h"
C
C ----------------------------------------------------------------------
      REAL BETAG(NLONW,NLATW)
C ----------------------------------------------------------------------
C
C*    *ADDITIONAL PARAMETERS FOR MULTITASKING
C
      INTEGER ITSKCA(9,NTASKS)
      EXTERNAL PHYS
      CHARACTER*1 CMODE
      LOGICAL LMLTSK
C
      LOGICAL NEWREAD, NEWFILE
      DATA NEWFILE / .FALSE. /
      ZERO = ' '
C
C ----------------------------------------------------------------------
C
C*    0. WRITE OUT NORMS FOR WAVE ENERGY, PEAK PERIOD AND WIND SPEED
C        -----------------------------------------------------------
      IF (ITEST .GE. 2) THEN
         WRITE(IU06,*)
         WRITE(IU06,*) '   SUB. WAMODEL: JUST ENTERED FOR CDTPRO= ',
     .   CDTPRO, ' AND CDATEA= ', CDATEA
         CALL FLUSH(IU06)
      ENDIF
      LLFLUSH = .TRUE.
C
      KTAG=200
C
      IF(CDTPRO.EQ.CDATEA) THEN
C
        CDTINTTBAK=CDTINTT
        CDTINTSBAK=CDTINTS
        CDTINTT=CDTPRO
        CDTINTS=CDTPRO
C
        DO IG=1,IGL
C
         IGH1=(IG/2)*2
C
         DO IJ=IJS(IG),IJL(IG)
            U10NEW(IJ) = U10OLD(IJ,IG)
            THWNEW(IJ) = THWOLD(IJ,IG)
            USNEW(IJ) = USOLD(IJ,IG)
         ENDDO
C
         IF (IGH1.EQ.IG) THEN
CREFRA
             IF (IREFRA.EQ.2) THEN
                CALL INTPOL (FL1, SL, IJS(IG), IJL(IG), IG)
                CALL OUTBS (SL, FL2, IJS(IG), IJL(IG), IG, IGL,
     1          IU25, IU26,THWOLD,USOLD,U10NEW,THWNEW,USNEW,TAUW,
     2          NSTART,NEND)
CREFRA
             ELSE
                CALL OUTBS (FL1, FL2, IJS(IG), IJL(IG), IG, IGL,
     1          IU25, IU26,THWOLD,USOLD,U10NEW,THWNEW,USNEW,TAUW,
     2          NSTART,NEND)
             ENDIF
          ELSE
CREFRA
             IF (IREFRA.EQ.2) THEN
                CALL INTPOL (FL1, SL, IJS(IG), IJL(IG), IG)
                CALL OUTBS (SL, FL1, IJS(IG), IJL(IG), IG, IGL,
     1          IU25, IU26,THWOLD,USOLD,U10NEW,THWNEW,USNEW,TAUW,
     2          NSTART,NEND)
CREFRA
             ELSE
                CALL OUTBS (FL1, FL3, IJS(IG), IJL(IG), IG, IGL,
     1          IU25, IU26,THWOLD,USOLD,U10NEW,THWNEW,USNEW,TAUW,
     2          NSTART,NEND)
             ENDIF
          ENDIF
C
        ENDDO
C
        IF(GFLAG(1)) THEN
          CALL WAMNORM(WHGTTG,NGX,NGY,ZMISS,EPSMIN,WHGTAVG,WHGTMIN,
     &                 WHGTMAX,LMESSPASS,IRANK,IPFGTBL(1),1,KTAG)
          KTAG=KTAG+1
          IF(IRANK.EQ.1) THEN
            WRITE(IU06,*) ' '
            WRITE(IU06,*) '   WAMNORM FOR WAVE ENERGY ON ',CDTPRO
            WRITE(IU06,*) '          AVERAGE            MIMIMUM  ',
     &                    '           MAXIMUM'
            WRITE(IU06,*) '  ',WHGTAVG,WHGTMIN,WHGTMAX
            WRITE(IU06,*) ' '
          ENDIF
        ENDIF
        IF(GFLAG(6)) THEN
          CALL WAMNORM(WPKFTG,NGX,NGY,ZMISS,0.,WPKFAVG,WPKFMIN,WPKFMAX,
     &                 LMESSPASS,IRANK,IPFGTBL(6),1,KTAG)
          KTAG=KTAG+1
          IF(IRANK.EQ.1) THEN
            WRITE(IU06,*) '   WAMNORM FOR WAVE PEAK FREQUENCY: '
            WRITE(IU06,*) '          AVERAGE            MIMIMUM  ',
     &                    '           MAXIMUM'
            WRITE(IU06,*) '  ',WPKFAVG,WPKFMIN,WPKFMAX
            WRITE(IU06,*) ' '
          ENDIF
        ENDIF
        IF(GFLAG(10)) THEN
          CALL WAMNORM(U10G,NGX,NGY,ZMISS,0.,U10AVG,U10MIN,U10MAX,
     &                 LMESSPASS,IRANK,IPFGTBL(10),1,KTAG)
          KTAG=KTAG+1
          IF(IRANK.EQ.1) THEN
            WRITE(IU06,*) '   WAMNORM FOR WAVE WIND SPEED: '
            WRITE(IU06,*) '          AVERAGE          MIMIMUM  ',
     &                    '      MAXIMUM'
            WRITE(IU06,*) '  ',U10AVG,U10MIN,U10MAX
          ENDIF
        ENDIF
C
      ENDIF
C
C                                                   
C*    0.1 SAVE INITIAL INTEGRATED FIELDS FOR FC RUN 
C         -----------------------------------------
      IF((CDTPRO.EQ.CDATEA).AND.
     &(MARSTYPE.EQ.'cf'.OR.MARSTYPE.EQ.'pf'.OR.MARSTYPE.EQ.'fc')) THEN
        CALL OUTINT (IU06, IU20, IU21, IU30, ITEST, MFDBSF)
        ISTAT = IFLUSHFDB (MFDBSF)

        CDTINTT=CDTINTTBAK
        CDTINTS=CDTINTSBAK
        IF (ITEST .GE. 2) THEN
           WRITE(IU06,*)
             WRITE(IU06,*) '   SUB. WAMODEL: INITIAL FIELDS SAVED '
             WRITE(IU06,*) '                 FOR ENSEMBLE RUN     '
             WRITE(IU06,*) '    '
             CALL FLUSH(IU06)
        ENDIF
      ENDIF
C                                                   
C*    1. ADVECTION TIME LOOP.
C        --------------------
C
      DO 1000 KADV = 1,NADV

         IF(.NOT.ALLOCATED(FL4).AND..NOT.LMESSPASS)
     1      ALLOCATE(FL4(NINF-1:NSUP,NANG,NFRE))
C
C*    1.1 FIX END DATE OF THIS PROPAGATION STEP AND OUTPUT TIMES.
C         -------------------------------------------------------
C
      CDTPRA = CDTPRO
      CALL INCDATE(CDTPRO,IDELPRO)
      IF (ITEST.GE.2) THEN
         WRITE(IU06,*) '   SUB. WAMODEL: START OF PROPAGATION '
         WRITE(IU06,*) '     START DATE IS    CDTPRA = ',CDTPRA
         WRITE(IU06,*) '     END DATE WILL BE CDTPRO = ',CDTPRO
         CALL FLUSH (IU06)
      ENDIF
C
C         UPDATE OUTPUT TIMES.
C
      IF (NOUTT.GT.0) THEN
         CDTINTT = ZERO
         CDTINTS = ZERO
         CDTSPT  = ZERO
         CDTSPS  = ZERO
         DO 1001 J=1,NOUTT
            IF (CDTPRO.EQ.COUTT(J)) THEN
               IF (PFLAG20.OR.FFLAG20.OR.GFLAG20) CDTINTT = COUTT(J)
               IF (PFLAG21.OR.FFLAG21.OR.GFLAG21) CDTINTS = COUTT(J)
               IF (PFLAG25.OR.FFLAG25.OR.GFLAG25) CDTSPT  = COUTT(J)
               IF (PFLAG26.OR.FFLAG26.OR.GFLAG26) CDTSPS  = COUTT(J)
            ENDIF
 1001    CONTINUE
      ELSE
         IF ((PFLAG20.OR.FFLAG20.OR.GFLAG20).AND.CDTINTT.LT.CDTPRO)
     1      CALL INCDATE (CDTINTT,IDELINT)
         IF ((PFLAG21.OR.FFLAG21.OR.GFLAG21).AND.CDTINTS.LT.CDTPRO)
     1      CALL INCDATE (CDTINTS,IDELINS)
         IF ((PFLAG25.OR.FFLAG25.OR.GFLAG25).AND.CDTSPT.LT.CDTPRO)
     1      CALL INCDATE (CDTSPT,IDELSPT)
         IF ((PFLAG26.OR.FFLAG26.OR.GFLAG26).AND.CDTSPS.LT.CDTPRO)
     1      CALL INCDATE (CDTSPS,IDELSPS)
      ENDIF
      IF (CDTRES.LT.CDTPRO) CALL INCDATE(CDTRES,IDELRES)
C
C*    1.2 RESET FILES.
C         ------------
      IF(.NOT.LMESSPASS) REWIND IU12

      IF (IGL.NE.1) THEN
         REWIND IU08
         REWIND IU11
         REWIND IU13
         REWIND IU14
C
C*    1.3 CHANGE WORKING FILE UNITS.
C         --------------------------
C
         IFIL=IU11
         IU11=IU12
         IU12=IFIL
         IFIL=IU13
         IU13=IU14
         IU14=IFIL
C
C*    1.4 READ FIRST BLOCK OF SPECTRA, IF MULTI BLOCK VERSION.
C         ----------------------------------------------------
C
         READ(IU11) FL1
         IF (ITEST.GE.2)
     1      WRITE(IU06,*) '   SUB. WAMODEL: START OF READING',
     2                    ' FIRST BLOCK OF SPECTRA'
      ENDIF
C
C*    1.5  LOOP FOR BLOCKS OF LATITUDES.
C          -----------------------------
C
      DO 1500 IG=1,IGL
C
C
C*    1.5.1 EVEN OR ODD BLOCK NUMBER
C          -------------------------
C
         IGH1=(IG/2)*2
C
C*    1.5.2 INPUT NEIGHBOURING GRID POINT INDICES (COMMON BLOCK UBUF).
C           ----------------------------------------------------------
C
         IF (IGL.NE.1) THEN
            READ (IU08) KLAT, KLON
CSHALLOW
C
C*    1.5.3 COMPUTE SHALLOW WATER TABLE INDICES.
C           ------------------------------------
C
            IF (ISHALLO.NE.1) THEN
               DO 1530 IJ=NINF,IJLT(IG)
                  XD = LOG(DEPTH(IJ,IG)/DEPTHA)/LOG(DEPTHD)+1.
                  ID = NINT(XD)
                  ID = MAX(ID,1)
                  INDEP(IJ) = MIN(ID,NDEPTH)
 1530          CONTINUE
            ENDIF
CSHALLOW
            IF (ITEST.GE.2) THEN
               IF (ITESTB.GE.IG) THEN
                  WRITE(IU06,*) '   SUB. WAMODEL: COMMON UBUF ',
     1                 'READ AND INITIALIZED FOR BLOCK IG = ', IG
               ENDIF
            ENDIF
C
C*    1.5.4 COUPLING WITH NEIGHBOURING BLOCKS IG +- 1 AND START
C*          INPUT OF SPECTRA FOR BLOCK IG+1.
C           ----------------------------------------------------
C
            IF(IGH1.EQ.IG) THEN
               CALL FILLBL(FL2, IG, IU13)
               IF (IG.NE.IGL)
     1            READ(IU11) FL1
            ELSE
               CALL FILLBL(FL1, IG, IU13)
               IF (IG.NE.IGL)
     1            READ(IU11) FL2
            ENDIF
            IF (ITEST.GE.2) THEN
               IF (ITESTB.GE.IG) THEN
                  WRITE(IU06,*) '   SUB. WAMODEL: LATITUDES ADDED',
     1                          ' TO BLOCK AND START OF READING',
     2                       ' FOR NEXT BLOCK'
               ENDIF
            ENDIF
         ENDIF
C
C
C*    1.5.5 SET TIME COUNTER.
C           -----------------
C
         CDATE   = CDTPRA
         CDATEWH = CDATEWO    
         NEWREAD = .FALSE.
C
 1550    CONTINUE                       
C
C
C*    1.5.5.1 READING OF NEW WINDS IF NEEDED.
C             -------------------------------
C        
         CALL NEWWIND(IJS(IG),IJL(IG),IG,IGL,CDATE,CDATEWH,
     1                NEWREAD,NEWFILE,U10OLD,THWOLD,U10NEW,THWNEW)
         IF (ITEST.GE.2) THEN
            IF (ITESTB.GE.IG) THEN
               WRITE(IU06,*) '   SUB. WAMODEL: NEWWIND CALLED' 
            ENDIF
         ENDIF
C
C*    1.5.5.2 COMPUTATION OF PROPAGATION &&& 
C*            INTEGRATION OF SOURCE TERMS OVER SUB TIME STEPS BETWEEN
C*            PROPAGATION TIME STEPS.
C             -------------------------------------------------------
C
! DETERMINE WHAT FORM OF PARALELLISM IS USED BASED ON THE TYPE
! OF MACHINE ARCHTECTURE

       IF(LMESSPASS) THEN
C
         IF(IGH1.NE.IG) THEN
! IF CDATE CORRESPONDS TO A PROPAGATION TIME , THEN CALL MPEXCHNG
! TO OBTAIN INFORMATION AT NEIGHBORING GRID POINTS LOCATED ON ADJACENT
! PROCESSES.
           IF(CDATE.EQ.CDTPRA) THEN
              CALL MPEXCHNG(FL1,NSTART,NEND,KLENBOT,KLENTOP)
              IF (ITEST.GE.2) THEN
               WRITE(IU06,*) '   SUB. WAMODEL: MPEXCHNG CALLED' 
              ENDIF
           ENDIF

           IF (CDATE.EQ.CDTPRA) THEN
              CALL PROPAGS(FL1,FL3,IJS(IG),IJL(IG),IG,SL,FCONST)
              IF (ITEST.GE.2) THEN
                WRITE(IU06,*) '   SUB. WAMODEL: PROPAGS CALLED' 
              ENDIF
           ENDIF
           CALL IMPLSCH (FL3, FL1, IJS(IG), IJL(IG), IG,
     1                   THWOLD,USOLD,TAUW,Z0OLD,
     2                   U10NEW,THWNEW,USNEW,Z0NEW,SL,FCONST)

           IF (ITEST.GE.2) THEN
             WRITE(IU06,*) '   SUB. WAMODEL: IMPLSCH CALLED' 
             CALL FLUSH (IU06)
           ENDIF
         ELSE
           WRITE (IU06,*) '******************************************'
           WRITE (IU06,*) '*                                        *'
           WRITE (IU06,*) '*    FATAL ERROR                         *'
           WRITE (IU06,*) '*    ===========                         *'
           WRITE (IU06,*) '* A MULTI BLOCK VERSION RUN WAS REQUESTED*'
           WRITE (IU06,*) '* THIS OPTION IS NOT IMPLEMENTED FOR     *'
           WRITE (IU06,*) '* MESSAGE PASSING PROTOCOL               *'
           WRITE (IU06,*) '*                                        *'
           WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.      *'
           WRITE (IU06,*) '*                                        *'
           WRITE (IU06,*) '******************************************'
           CALL ABORT1
         END IF

       ELSE
         NCHUNK = 0
         NLENGTH = (IJL(IG)-IJS(IG)+1)
         NWORK   = NLENGTH/NSIZE
         IF (MOD(NLENGTH,NSIZE).NE.0) NWORK = NWORK + 1
C
         NTSKCAP = 3
         NSTKSZ  = 1
         NSTKINC = 1
         DO 15520 JROC=1,NTASKS-1
            ITSKCA(1,JROC) = NTSKCAP
            ITSKCA(3,JROC) = JROC
            ITSKCA(4,JROC) = 1 ! !!!! JUST TO GET PAST COMP 'STACKSZW'
            ITSKCA(5,JROC) = NSTKSZ
            ITSKCA(6,JROC) = 1 ! !!!! JUST TO GET PAST COMP 'STACKINW'
            ITSKCA(7,JROC) = NSTKINC
C
            IF(IGH1.NE.IG) THEN
               CALL TSKSTART(ITSKCA(1,JROC),PHYS,FL3,FL1,FL4,IJS(IG),
     1                       IJL(IG),IG,NSIZE,NWORK,CDATE,CDTPRA)
            ELSE
               CALL TSKSTART(ITSKCA(1,JROC),PHYS,FL3,FL2,FL4,IJS(IG),
     1                       IJL(IG),IG,NSIZE,NWORK,CDATE,CDTPRA)
            END IF
15520    CONTINUE
         IF (ITEST.GE.2) CALL FLUSH(IU06)
C
         IF(IGH1.NE.IG) THEN
           CALL PHYS (FL3,FL1,FL4,IJS(IG),IJL(IG),IG,NSIZE,NWORK,
     1                CDATE, CDTPRA,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     2                U10NEW,THWNEW,USNEW,Z0NEW,SL,FCONST)
         ELSE
           CALL PHYS (FL3,FL2,FL4,IJS(IG),IJL(IG),IG,NSIZE, NWORK,
     1                CDATE, CDTPRA,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     2                U10NEW,THWNEW,USNEW,Z0NEW,SL,FCONST)
         END IF
C
         DO 15521 JROC=1,NTASKS-1
            CALL TSKWAIT(ITSKCA(1,JROC))
15521    CONTINUE

       ENDIF

C
C
C*    1.5.5.3 COPY AND CLOSE FILES IF NEEDED.
C             -------------------------------
C        
         CALL CLOSEND(IJS(IG),IJL(IG),IG,IGL,CDATE,CDATEWH,
     1                NEWREAD,NEWFILE,U10OLD,THWOLD,U10NEW,THWNEW)
         IF (ITEST.GE.2) THEN
            IF (ITESTB.GE.IG) THEN
               WRITE(IU06,*) '   SUB. WAMODEL: CLOSEND CALLED' 
            ENDIF
         ENDIF
C        
C*    1.5.5.4 UPDATE TIME;IF TIME LEFT BRANCH BACK TO 1.5.5 
C             ---------------------------------------------
C
         CALL INCDATE(CDATE,IDELT)   
         IF (CDATE.LT.CDTPRO) GO TO 1550
C
C        END OF TIME LOOP ALL TIME STEPS DONE
C
         IF (ITEST.GE.2) THEN
            IF (ITESTB.GE.IG) THEN
               WRITE(IU06,*) '   SUB. WAMODEL: SOURCE FUNCTIONS',
     1                       ' INTEGRATED:'
            ENDIF
         ENDIF
CICE
C
C*    1.5.6 SET FL3 ON ICE POINTS TO ZERO
C           -----------------------------
C
         IF (LICERUN) THEN
            IF (ITEST.GE.1) THEN
               WRITE(IU06,*) '   SUB. WAMODEL: SPECTRUM = 0 AT',
     1                       ' ICEPOINTS'
            ENDIF
            CALL SETICE(FL3,IJS(IG),IJL(IG),IG)
         ENDIF
C
CICE
CNEST
C
C*    1.5.7 INPUT OF BOUNDARY VALUES.
C           -------------------------
C
         IF (IBOUNF.EQ.1) THEN
            CALL BOUINPT (FL3, IG, IU02)
            IF (ITEST.GE.2) THEN
               IF (ITESTB.GE.IG) WRITE(IU06,*)
     1            '   SUB. WAMODEL: BOUNDARY VALUES INSERTED'
            ENDIF
         ENDIF
CNEST
C
C*    1.5.8 OUTPUT OF SECOND LATITUDES FOR NEXT TIMESTEP
C*          START OF OUTPUT OF ALL SPECTRA IN BLOCK.
C           --------------------------------------------
C
         IF (IGL.NE.1) THEN
            IF (IG.NE.1) CALL SPLITBL(FL3, IJS(IG), IJL2(IG), IU14)
            WRITE(IU12) FL3
            IF (ITEST.GE.2) THEN
               IF (ITESTB.GE.IG) WRITE(IU06,*)
     1            '   SUB. WAMODEL: OVERLAPPING LATITUDES',
     2            ' SAVED AND START OF BLOCK OUTPUT'
            ENDIF
         ENDIF
CNEST
C
C*    1.5.9 OUTPUT OF BOUNDARY POINTS.
C           --------------------------
C
         IF (IBOUNC.EQ.1) THEN
            CALL OUTBC (FL3, IJS(IG), IJL(IG), IG, IU19, THWNEW, USNEW)
            IF (ITEST.GE.2) THEN
               IF (ITESTB.GE.IG) WRITE(IU06,*)
     1            '   SUB. WAMODEL: BOUNDARY OUTP',
     2            'UT (COURSE GRID) DONE IN SUB OUTBC'
            ENDIF
         ENDIF
CNEST
C
C*    1.5.10 MODEL OUTPUT INTEGRATED DATA ARE SAVED IN COMMON.
C*           SPECTRA ARE SAVED ON DISK AND PRINTED.
C            -------------------------------------------------
C
         IF (ITEST.GE.2) THEN
           IF (ITESTB.GE.IG) THEN
             WRITE(IU06,*) '   SUB. WAMODEL: MODEL ',
     1       'OUTPUT? CDTINTT=',CDTINTT,' CDTINTS=',CDTINTS
             WRITE(IU06,*) '   SUB. WAMODEL: MODEL ',
     1       'OUTPUT? CDTSPT=',CDTSPT,' CDTSPS=',CDTSPS
             WRITE(IU06,*) '   SUB. WAMODEL: MODEL ',
     1       'OUTPUT? IGH1=',IGH1,' IG=',IG
           ENDIF
         ENDIF
C
         IF(.NOT.LWCOU .AND. .NOT. LDSTOP) THEN
           CALL CHESIG (IU06, ITEST, IRANK, NPROC, LDSTOP, LDWRRE)
         ENDIF
C
         IF (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO .OR.
     1       CDTSPT.EQ.CDTPRO .OR. CDTSPS.EQ.CDTPRO .OR.
     2       CDATEE.EQ.CDTPRO .OR. (LDWRRE .AND. KADV.EQ.NADV ))
     3     THEN
C
            IF (IGH1.EQ.IG) THEN
CREFRA
               IF (IREFRA.EQ.2) THEN
                  CALL INTPOL (FL3, SL, IJS(IG), IJL(IG), IG)
                  CALL OUTBS (SL, FL2, IJS(IG), IJL(IG), IG, IGL,
     1            IU25, IU26,THWOLD,USOLD,U10NEW,THWNEW,USNEW,TAUW,
     2            NSTART,NEND)

CREFRA
               ELSE
                  CALL OUTBS (FL3, FL2, IJS(IG), IJL(IG), IG, IGL,
     1            IU25, IU26,THWOLD,USOLD,U10NEW,THWNEW,USNEW,TAUW,
     2            NSTART,NEND)
               ENDIF
            ELSE
CREFRA
               IF (IREFRA.EQ.2) THEN
                  CALL INTPOL (FL3, SL, IJS(IG), IJL(IG), IG)
                  CALL OUTBS (SL, FL1, IJS(IG), IJL(IG), IG, IGL,
     1            IU25, IU26,THWOLD,USOLD,U10NEW,THWNEW,USNEW,TAUW,
     2            NSTART,NEND)
CREFRA
               ELSE
                  CALL OUTBS (FL3, FL1, IJS(IG), IJL(IG), IG, IGL,
     1            IU25, IU26,THWOLD,USOLD,U10NEW,THWNEW,USNEW,TAUW,
     2            NSTART,NEND)
               ENDIF
            ENDIF

            IF (ITEST.GE.2) THEN
               IF (ITESTB.GE.IG) WRITE(IU06,*)
     1             '   SUB. WAMODEL: MODEL OUTPUT DONE IN SUB OUTBS'
            ENDIF
         ENDIF
C
C
C*    1.5.11 STORE CHARNOCK CONSTANT FOR ATMOSPHERIC MODEL.
C            ----------------------------------------------
C
      IF(LWCOU
     1   .AND.
     2   (KADV.EQ.NADV)
     3   .AND.
     4   ((CDTPRO.NE.CDATEF).OR.((CDTPRO.EQ.CDATEF).AND.(IASSI.NE.1)))
     5    ) THEN
        CALL OUTBETA ( IJS(IG), IJL(IG), IG, IGL,USNEW,Z0NEW,
     1                 NSTART, NEND, BETAG, NLONW, NLATW)
      ENDIF
C
C
C*    BRANCHING BACK TO 1.5 FOR NEXT BLOCK OF LATITUDES
C
 1500 CONTINUE
C
C*    1.6 IF ONE BLOCK VERSION COPY RESULTS.
C         ----------------------------------
C
      IF (IGL.EQ.1) THEN
         DO M=1,NFRE
           DO K=1,NANG
             DO IJ=NINF-1,NSUP
               FL1(IJ,K,M) = FL3(IJ,K,M)
             ENDDO
           ENDDO
         ENDDO
         IF (ITEST.GE.2)
     1      WRITE(IU06,*) '   SUB. WAMODEL: BLOCK COPIED'
      ENDIF
C
C*    1.7 ONE PROPAGATION TIMESTEP DONE FOR ALL BLOCKS.
C         ---------------------------------------------
C
C
C         WRITE INTEGRATED DATA TO FILE AND/OR PRINTER
C         DATA WERE COLLECTED INSIDE THE BLOCK LOOP.
C
      IF (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO) THEN
C
C       WRITE OUT NORMS FOR WAVE ENERGY, PEAK FREQUENCY AND WIND SPEED 
C
        IF(GFLAG(1)) THEN
          CALL WAMNORM(WHGTTG,NGX,NGY,ZMISS,EPSMIN,WHGTAVG,WHGTMIN,
     &                 WHGTMAX,LMESSPASS,IRANK,IPFGTBL(1),1,KTAG)
          KTAG=KTAG+1
          IF(IRANK.EQ.1) THEN
            WRITE(IU06,*) ' '
            WRITE(IU06,*) '   WAMNORM FOR WAVE ENERGY ON ',CDTPRO
            WRITE(IU06,*) '          AVERAGE            MIMIMUM  ',
     &                    '           MAXIMUM'
            WRITE(IU06,*) '  ',WHGTAVG,WHGTMIN,WHGTMAX
            WRITE(IU06,*) ' '
          ENDIF
        ENDIF
        IF(GFLAG(6)) THEN
          CALL WAMNORM(WPKFTG,NGX,NGY,ZMISS,0.,WPKFAVG,WPKFMIN,WPKFMAX,
     &                 LMESSPASS,IRANK,IPFGTBL(6),1,KTAG)
          KTAG=KTAG+1
          IF(IRANK.EQ.1) THEN
            WRITE(IU06,*) '   WAMNORM FOR WAVE PEAK FREQUENCY: '
            WRITE(IU06,*) '          AVERAGE            MIMIMUM  ',
     &                    '           MAXIMUM'
            WRITE(IU06,*) '  ',WPKFAVG,WPKFMIN,WPKFMAX
            WRITE(IU06,*) ' '
          ENDIF
        ENDIF
        IF(GFLAG(10)) THEN
          CALL WAMNORM(U10G,NGX,NGY,ZMISS,0.,U10AVG,U10MIN,U10MAX,
     &                 LMESSPASS,IRANK,IPFGTBL(10),1,KTAG)
          KTAG=KTAG+1
          IF(IRANK.EQ.1) THEN
            WRITE(IU06,*) '   WAMNORM FOR WAVE WIND SPEED: '
            WRITE(IU06,*) '          AVERAGE          MIMIMUM  ',
     &                    '      MAXIMUM'
            WRITE(IU06,*) '  ',U10AVG,U10MIN,U10MAX
          ENDIF
        ENDIF
C
        CALL OUTINT (IU06, IU20, IU21, IU30, ITEST, MFDBSF)
        IF(LMESSPASS) THEN
          CALL MPE_BARRIER(IERR)
          IF(IERR.LT.0)CALL MPEI_ABORT('MPE_BARRIER ERROR in WAMODEL')
        ENDIF
      ENDIF

C
C         PRINT TIME.
C
         WRITE(IU06,112) CDTPRO
         IF (ITEST .GE. 2) CALL FLUSH (IU06)
  112 FORMAT(/,3X,'!!!!!!!!!!!!!! ',
     1            'WAVE FIELDS INTEGRATED  DATE IS: ',A10,
     2            '  !!!!!!!!!!!!!! ')
C
C*    1.8 OUTPUT FILES AND RECOVERY FILES ARE DISPOSED WHEN
C         TIME REACHES THE DISPOSE DATE OR WHEN THE MODEL
C         HAS BEEN SIGNALLED TO DO SO.
C         -------------------------------------------------
C
      IF (ITEST.GE.1) THEN
        WRITE(IU06,*) " WAMODEL: 1.8  CDTRES=", CDTRES, " CDTPRO=",
     .   CDTPRO, " LDWRRE=", LDWRRE, "KADV=", KADV, " NADV=", NADV
        WRITE(IU06,*)
     . " CDATEF=", CDATEF, "CDATEE=", CDATEE, " CDATER=", CDATER,
     .   " CDATES=", CDATES, "IREST=", IREST
      ENDIF
C
      IF (CDTRES.EQ.CDTPRO.OR.CDATEE.EQ.CDTPRO.OR.CDTPRO.EQ.CDATER.OR.
     1    (LDWRRE .AND. KADV.EQ.NADV)) THEN
C
C*    1.8.2 STORE 2-D SPECTRA AT SELECTED GRID POINTS,
C           INTEGRATED PARAMETERS OF ENTIRE GRID.
C           ------------------------------------------
C
C           SAVE INTEGRATED DATA FILE.
C
         IF (GFLAG20 .AND. .NOT. LFDB) THEN
           IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
             CALL GSFILE (IU06, IU30, 0, CDTPRO, CDATEF, 'MPP', 'S')
           ENDIF
         ENDIF
         IF (FFLAG20) THEN
           IF((LMESSPASS.AND.(IRANK.EQ.NPROC)).OR.(.NOT.LMESSPASS)) THEN
             CALL GSFILE (IU06, IU20, 0, CDTPRO, CDATEF, 'MAP', 'S')
           ENDIF
         ENDIF
C
C           SAVE INTEGRATED SWELL AND SEA DATA FILE.
C           ----------------------------------------
C
         IF (FFLAG21) THEN
            IF((LMESSPASS.AND.(IRANK.EQ.MIN(2,NPROC)))
     $         .OR.(.NOT.LMESSPASS)) THEN
            CALL GSFILE (IU06, IU21, 0, CDTPRO, CDATEF, 'SWE', 'S')
            ENDIF
         ENDIF
C
C           SAVE SPECTRA FILE.
C
         IF (FFLAG25) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
              CALL GSFILE (IU06, IU25, 0, CDTPRO, CDATEF, 'OUT', 'S')
          ENDIF
         ENDIF
C
C           SAVE SWELL SPECTRA FILE.
C
         IF (FFLAG26) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
              CALL GSFILE (IU06, IU26, 0, CDTPRO, CDATEF, 'SWS', 'S')
          ENDIF
         ENDIF
C
CNEST
C
C           SAVE BOUNDARY VALUE FILE.
C
         IF (IBOUNC.EQ.1) THEN
            CALL GSFILE (IU06, IU19, 0, CDTPRO, CDATEF, 'CBO', 'S')
            IF (CDTPRO.LT.CDATEE)
     1         CALL HEADBC (NBOUNC, IDELPRO, TH(1), FR(1), IU19, IU06)
         ENDIF
CNEST
C
C*    1.8.3 SAVE RESTART FIELDS.
C           --------------------
C

        IF (   (IASSI.NE.1 .OR. CDTPRO.GT.CDATEF)
     1         .AND.
     2         (  ( (IREST.GE.1 .AND. IREST.LE.4)
     3              .AND.
     4              (CDTPRO.EQ.CDATER .OR. CDTPRO.LE.CDATES) )
     5            .OR.
     6            (LDWRRE .AND. KADV.EQ.NADV)  )   ) THEN


C*       STRESS FIELD IS DUMPED FOR RESTART.
C        -----------------------------------
C        NOTE :
C        STRESS FIELD IS DUMPED FOR RESTART
C        ONLY IF CURRENT TIME CORRESPONDS TO CDATER
C        EXCEPT IF CDATER SET TO 0000000000

         IF((CDTPRO .EQ. CDATER ) .OR.
     .      (CDATER .EQ. '0000000000' ) .OR.
     .      (LDWRRE)) THEN

            IF ( .NOT.LGRIBIN ) THEN
              CALL SAVSTRESS(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,NSTART,NEND,
     1                       LDSTOP, LDWRRE, IREST,CDTPRO,CDATEF)
            ENDIF

         ENDIF
C
C         SAVE SPECTRUM
C         -------------

          IF ( LGRIBIN .AND. LMESSPASS) THEN
            CALL SPEC2FDB (IU06, NPROC, IRANK, LWCOU, TH, FR,
     1                     FL1, NINF, NSUP, NANG, NFRE, NGX, NGY,IRGG,
     2                     NLONRGG, KXLT, IXLG, NIBLO, NBLO, ZMISS,
     3                     AMONOP, AMOWEP, AMOEAP, AMOSOP, XDELLA,
     4                     XDELLO, NSTART, NEND, 
     5                     CDTPRO, CDATEF, CDATER, CDATES, CDATEA,
     6                     MARSTYPE, YCLASS, YEXPVER, CFDB2DSP,
     7                     NENSFNB, NTOTENS, LFDBIOOUT, ITEST, MFDBSF)
            ISTAT = IFLUSHFDB (MFDBSF)
            LLFLUSH = .FALSE.
          ELSE
            CALL SAVSPEC(FL3,NSTART,NEND,LDSTOP,LDWRRE,IREST,
     1                   CDTPRO,CDATEF,CDATER,CDATES,MARSTYPE,
     2                   YCLASS, YEXPVER, CFDB2DSP)
          ENDIF
          WRITE(IU06,*) ' '
          WRITE(IU06,*) '  WAVE SPECTRA FILES DISPOSED AT........',
     1                      ' CDTPRO  = ', CDTPRO
          WRITE(IU06,*) ' '

C
C*    1.8.4 UPDATE, WRITE AND SAVE WAMINFO FILE.
C           ------------------------------------
C
          IF (LDWRRE .AND. KADV.EQ.NADV .AND.
     1        (LMESSPASS .AND. IRANK .EQ. 1 .OR. .NOT. LMESSPASS )) THEN
            ICH = 7 
            CALL DIFDATE (CDATEF,CDATEE,IFOREPD)
            IF (CDTPRO.LE.CDATEF) THEN
               CALL DIFDATE (CDTPRO,CDATEF,IANALPD)
               CBPLTDT = CDTPRO
               NFCST = 1
            ELSE
               NFCST = 0
               IANALPD = 0
               CBPLTDT = CDATEF
               CALL DIFDATE (CDTPRO,CDATEE,IFOREPD)
            ENDIF
            ISTAT = 0
            IF (CDATE.EQ.CDATEE) ISTAT(1) = 1
            IDELWIN = IDELWI
C
            CEPLTDT = CDATEF
C
            OPEN (IU04,FILE=CWI(1:ICPLEN+8),
     1            STATUS='UNKNOWN',FORM='FORMATTED')
            CALL WRITSTA (IU04, CDTPRO, CDATEE, IANALPD, IFOREPD,
     1                    IDELWIN, CDATER, CDATES, CBPLTDT, CEPLTDT,
     2                    IASSI, NFCST, ISTAT)
 
            CLOSE (IU04)
            WRITE(IU06,*) ' WAMINFO FILE WRITTEN FOR RESTART...',
     1                      ' CDTPRO  = ', CDTPRO
            WRITE(IU06,*) '                                    ',
     1                      ' CDATEF  = ', CDATEF
            WRITE(IU06,*) ' TO ', CWI(1:ICPLEN+8)
          ENDIF
        END IF
      END IF
      IF ( LFDB .AND. LLFLUSH .AND. 
     .     (IASSI.NE.1 .OR. CDTPRO.GT.CDATEF) .AND.
     .     (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO) ) THEN
        ISTAT = IFLUSHFDB (MFDBSF)
        WRITE(IU06,*) ' DB ', MFDBSF , ' FLUSHED AT ',
     .    CDTPRO, ' FROM WAMODEL. '
        CALL FLUSH (IU06)
      ENDIF
C
C*    BRANCHING BACK TO 1.0 FOR NEXT PROPAGATION STEP.
C
 1000 CONTINUE
C
      IF(ALLOCATED(FL4)) DEALLOCATE(FL4)
C
      RETURN
      END
