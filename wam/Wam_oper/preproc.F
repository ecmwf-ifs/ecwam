      PROGRAM preproc 

C ----------------------------------------------------------------------
C
C**** *PROGRAM PREPROC* - PREPARE DATA (BUT NOT WINDS) FOR INPUT
C                         TO WAM WAVE MODELS.
C
C     SUSANNE HASSELMANN  MPI     JUNE 1986.
C
C     ANNEGRET SPEIDEL    MPI  OCTOBER 1988. MODFIED FOR CYCLE_2.
C
C     K. HUBBERT          POL     JUNE 1989  DEPTH AND CURRENT
C                                            REFRACTION.
C
C     H. GUNTHER   ECMWF/GKSS    APRIL 1990  LAND POINTS ARE REMOVED
C                                            FROM BLOCKS AND THE CODE
C                                            HAS BEEN RESTRUCTURED.
C
C     R. PORTZ     MPI         JANUARY 1991  NESTED GRID OPTION.
C
C     H. GUNTHER   ECMWF/GKSS    APRIL 1991  CYCLE_4 MODIFICATIONS.
C                                            MULTI-PART REMOVED.
C                                            NEW SOURCE FUNCTIONS.
C                                            LOG. DEPTH TABLE.
C     J. BIDLOT    ECMWF         SEPTEMBER 1996 REDUCED GRID.
C
C*    PURPOSE.
C     --------
C
C       TO ARRANGE A GRID FOR THE WAM WAVE MODEL AND COMPUTE
C       ALL FIXED MODEL PARAMETERS WHICH ARE STORED IN DIFFERENT
C       COMMON BLOCKS.
C
C     METHOD.
C     -------
C
C       A REPRESENTATIVE TOPOGRAPHIC DATA SET ON LAT-LONG
C       COORDINATES CONTAINING THE MODEL SQUARE BOX REGION IS
C       READ IN.THE MODEL REGION IS EXTRACTED AND INTERPOLATED
C       ONTO GIVEN LAT-LONG GRID INCREMENTS (SEE SUB TOPOAR).
C       THE PROGRAM CHECKS FOR A PERIODIC LATITUDE GRID. IF THE
C       GRID IS NOT PERIODIC A CLOSED BASIN IS ASSUMED.
C       THE PROGRAM DOES NOT DISTINGUISH BETWEEN DEEP AND SHALLOW
C       WATER.
C
C       -BLOCK STRUCTURE :
C        GRID POINTS ARE COLLECTED INTO A 1-DIMENSIONAL ARRAY,
C        BLOCKS OF MAXIMALLY NIBLO ELEMENTS,  GRID POINTS
C        (ONLY SEAPOINTS) ARE COUNTED ALONG LINES OF LATITUDES
C        FROM WEST TO EAST WORKING FROM SOUTH TO NORTH.
C        BLOCKS OVERLAP OVER TWO LATITUDE LINES,TO COMPUTE NORTH
C        -SOUTH ADVECTION TERMS.
C
C       -NESTED GRIDS: THE GRID GENERATED CAN BE A
C         - COARSE GRID WHICH MEANS OUTPUT OF SPECTRA
C                       FOR A FOLLOW UP FINE GRID RUN.
C         - FINE   GRID WHICH MEANS INPUT OF SPECTRA
C                       FROM  AN EARLIER COARSE GRID RUN.
C         - COARSE AND FINE GRID
C
C       - REFRACTION: CONTROLLED BY THE REFRACTION OPTION
C         A CURRENT FIELD IS READ, INTERPOLATED TO THE MODEL
C         GRID AND STORED IN THE GRID OUTPUT FILE.
C
C       - PARAMETERS FOR ARRAY DIMENSIONS: THE PRORAM CHECKS
C         ALL DIMENSIONS INTERNALLY. ONLY THE BLOCK LENGTH
C         (NIBLO) IS USED FOR THE SET UP OF THE GRID, ALL
C         THE OTHER PARAMETERS HAVE TO BE LARGE ENOUGH TO
C         GET A SUCCESFULL RUN OF THE JOB. AT THE END OF
C         THE OUTPUT PROTOCOLL A LIST IS PRINTED FOR THE
C         OPTIMAL SETTINGS OF THE DIMENSION.
C
C**   INTERFACE.
C     ----------
C
C       *PROGRAM* *PREPROC*
C
C       *IU01*   - LOGICAL UNIT FOR INPUT OF TOPOGRAPHIC DATA.
C                  (SEE SUB TOPOAR).
C       *IU02*   - LOGICAL UNIT FOR INPUT OF CURRENTS.
C                  (SEE SUB READCUR).
C       *IU03*   - LOGICAL UNIT FOR INPUT OF COARSE GRID
C                  BOUNDARY ORGANISATION (COMMON CBOUND).
C                  IF THIS IS A FINE GRID PREPROC.
C                  FORMATED IF IFORM = 2 OTHERWISE UNFORMATED.
C                  (SEE SUB MBOUNF).
C       *IU05*   - LOGICAL UNIT FOR USER INPUT. (SEE SUB UIPREP).
C       *IU06*   - LOGICAL UNIT FOR PRINTER OUTPUT UNIT
C       *IU07*   - LOGICAL UNIT FOR OUTPUT OF GRID ORGANISATION
C                  AND COMPUTED CONSTANTS. (UNFORMATED)
C                  (SEE SUB OUTCOM).
C       *IU08*   - LOGICAL UNIT FOR OUTPUT OF COMMON UBUF.
C                  (UNFORMATED) (SEE SUB OUTUBUF).
C       *IU09*   - LOGICAL UNIT FOR UNFORMATED OUTPUT OF COARSE
C                  GRID BOUNDARY ORGANISATION (COMMON CBOUND),
C                  IF THIS IS A COARSE GRID PREPROC.
C                  (SEE SUB MBOUNC).
C       *IU10*   - LOGICAl UNIT FOR UNFORMATED OUTPUT OF FINE
C                  GRID BOUNDARY ORGANISATION (COMMON CBOUND).
C                  IF THIS IS A FINE GRID PREPROC.
C                  (SEE SUB MBOUNF).
C       *IU17*   - SAME AS IU07 BUT FORMATED.
C       *IU18*   - SAME AS IU08 BUT FORMATED.
C       *IU19*   - SAME AS IU09 BUT FORMATED.
C       *IU20*   - SAME AS IU10 BUT FORMATED.
C
C       ALL UNITS ARE DEFINE IN SECTION 1. OF THIS PROGRAM.
C
C       COMMON BLOCKS COUPLE, CURRENT, FREDIR, INDNL, GRIDPAR, MAP,
C       COUT, TABLE, AND SHALLOW ARE WRITTEN TO UNIT IU07 AND/OR IU17.
C       ALL FREQUENCY AND DIRECTION DEPENDENT ARRAYS ARE WRITTEN FROM
C       1 TO THE USED NUMBER OF FREQUENCIES AND THE USED NUMBER OF
C       DIRECTIONS.
C       OTHER ARRAYS ARE WRITTEN ACCORDING TO THEIR DIMENSIONS.
C
C     EXTERNALS.
C     ----------
C
C       *ABORT*     - TERMINATES PROCESSING.
C       *ADJUST*    - CORRECTS LONGITUDE INPUT.
C       *AKI*       - COMPUTES WAVE NUMBER.
C       *CHECK*     - CHECKS CONSISTENCY OF BLOCK OVERLAPS.
C       *FINDB*     - FIND BLOCK AND GRID POINT NUMBERS.
C       *JAFU*      - ANGULAR INDEX OF NON LINEAR INTERACTION
C       *LOCINT*    - INTERPOLATE TO BLOCKS.
C       *MBLOCK*    - PREPARES ONE BLOCK
C       *MBOUNC*    - MAKE COARSE GRID BOUNDARY.
C       *MBOUNF*    - MAKE FINE   GRID BOUNDARY.
C       *MBOXB*     - MAKE BOX FOR FINE GRID IN COARSE GRID.
C       *MCOUT*     - PREPARES OUTPUT COMMON COUT
C       *MFREDIR*   - COMPUTES FREQUENCY/DIRECTION COMMON FREDIR
C       *MGRID*     - ARRANGES GRID FOR MODEL.
C       *MINTF*     - MAKE INTERPOLATION TABLES FOR BOUNDARY INPUT.
C       *MTABS*     - COMPUTES TABLES USED FOR SHALLOW WATER
C       *MUBUF*     - COMPUTES COMMON UBUF.
C       *NLWEIGT*   - COMPUTES NON LINEAR WEIGHTS IN COMMON INDNL
C       *OUTCOM*    - OUTPUT OF COMPUTED COMMONS
C       *OUTPP*     - OUTPUT OF GRID.
C       *OUTUBUF*   - OUTPUT OF BLOCK COMMON UBUF
C       *PACKI*     - PACKS AN INTEGER ARRAY.
C       *PACKR*     - PACKS A REAL ARRAY.
C       *READCUR*   - READ A CURRENT FIELD.
C       *STRESS*    - STRESS TABLE.
C       *TAUHF*     - HIGH FREQUENCY STRESS TABLE.
C       *TOPOAR*    - PREPARE TOPOGRAPHY FOR GRID.
C       *UIPREP*    - READS USER INPUT
C       *WAMCUR*    - PREPARE CURRENT FIELD.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "txtpara.h"
C
#include "parameter_preproc.h"
C
#include "txtpara_preproc.h"
C
#include "parameter_topo.h"
C
#include "txtpara_topo.h"
C
#include "parameter_current.h"
C
#include "txtpara_current.h"
C
#include "comcpbo.h"
C
#include "txtcpbo.h"
C
#include "comcinp.h"
C
#include "txtcinp.h"
C
#include "commesspass.h"
C
#include "comcoup.h"
C
#include "txtcoup.h"
C
#include "comcout.h"
C
#include "txtcout.h"
C
#include "comcurg.h"
C
#include "txtcurg.h"
C
#include "comcurr.h"
C
#include "txtcurr.h"
C
#include "comfpbo.h"
C
#include "txtfpbo.h"
C
#include "comfred.h"
C
#include "txtfred.h"
C
#include "comgrid.h"
C
#include "txtgrid.h"
C
#include "comindn.h"
C
#include "txtindn.h"
C
#include "commap.h"
C
#include "txtmap.h"
C
#include "comshal.h"
C
#include "txtshal.h"
C
#include "comtabl.h"
C
#include "txttabl.h"
C
#include "comtest.h"
C
#include "txttest.h"
C
#include "comubuf.h"
C
#include "txtubuf.h"
C
#include "parcons.h"
C
#include "txtcons.h"
C
      DIMENSION IA1(NGX, NGY)
C
C          *IA1*   INTEGER  TOPOGRAPHIC DATA IN MODEL GRID.
C
C ----------------------------------------------------------------------
C
      LMESSPASS=.FALSE.
C*    1. INITIALISATION OF INPUT/OUTPUT UNITS.
C        -------------------------------------
C
      IU01  = 1
      IU02  = 2
      IU03  = 3
      IU05  = 5
      IU06  = 6
      IU07  = 7
      IU08  = 8
      IU09  = 9
      IU10  = 10
      IU17  = 17
      IU18  = 18
      IU19  = 19
      IU20  = 20
C
C ----------------------------------------------------------------------
C
C*    2. USER INPUT AND LINEPRINTER PROTOCOL.
C        ------------------------------------
C
2000  CONTINUE
      CALL UIPREP (IU05, IFORM, IREFRA, ML, KL)
C
C ----------------------------------------------------------------------
C
C*    2.1 DETERMINE COORDINATES FOR GAUSSIAN GRID.
C     -------------------------------------------
C
2100  CONTINUE
C
C
C ----------------------------------------------------------------------
C
C*    3. INITIALISE TOTAL NUMBER OF BLOCKS,
C*       AND GRID INCREMENTS IN RADIANS AND METRES.
C        ------------------------------------------
C
 3000 CONTINUE
      IGL=0
      DELPHI = XDELLA*CIRC/360.
      DO 3001 K=1,NY
         XLAT       = (AMOSOP + REAL(K-1)*XDELLA)*RAD
         SINPH(K)   = SIN(XLAT)
         COSPH(K)   = COS(XLAT)
         IF (IRGG.EQ.1) THEN
            NLONRGG(K) = NINT(NX*COSPH(K))
            IF(MOD(NLONRGG(K),2).eq.1) NLONRGG(K) = NLONRGG(K)+1
         ELSE
            NLONRGG(K) = NX
         END IF      
         PLONS=(AMOEAP-AMOWEP) + XDELLO
         ZDELLO(K)  = PLONS/FLOAT(NLONRGG(K))
         DELLAM(K)  = ZDELLO(K)*CIRC/360.
 3001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    4. COMPUTE GRID INDEPENDENT COMMON BLOCKS.
C        ---------------------------------------
C
 4000 CONTINUE
C
C*    4.1 COMMON FREDIR (FREQUENCY/DIRECTION CONST).
C         ------------------------------------------
C
 4100 CONTINUE
      CALL MFREDIR (ML, KL)
C
C*    4.2 COMMON INDNL (WEIGHT OF NON-LINEAR INTERACTION).
C         ------------------------------------------------
C
 4200 CONTINUE
      CALL NLWEIGT (ML, KL)
C
C*    4.3 COMMON SHALLOW (SHALLOW WATER TABLES).
C         --------------------------------------
C
      CALL MTABS (ML, KL)
C
C*    4.4 COMMON COUPLE.
C         --------------
C
      BETAMAX = 1.20
      ZALP    = 0.0110
      ALPHA   = 0.0090
      XKAPPA  = 0.40
C
C*    UP TO JPLEVC LEVELS CAN BE SPECIFIED.
C
      XNLEV(1)= 10.0
      XNLEV(2)= 31.0
C
C*    4.4 COMMON TABLE (STRESS TABLES).
C         -----------------------------
C
      CALL STRESS(IU06)
      CALL TAUHF (FR(ML))
C
C ----------------------------------------------------------------------
C
C*    5. GENERATE OUTPUT GRID INFORMATION.
C        ---------------------------------
C
 5000 CONTINUE
C
C*    5.1 READ IN TOPOGRAPHY AND ARRANGE ON REQUESTED MODEL AREA OF
C*        REQUESTED RESOLUTION.
C         ---------------------------------------------------------
C
 5100 CONTINUE
      CALL TOPOAR (IU01, IA1)
C
C*    5.2 COMPUTATION OF BLOCKS.
C         ----------------------
C
 5200 CONTINUE
      CALL MGRID (IA1)
      IF (ITEST.GT.0) WRITE (IU06,*) ' SUB MGRID DONE'
C
C*    5.3 COMPUTE OUTPUT POINT INDICES (COMMON COUT).
C         -------------------------------------------
C
 5300 CONTINUE
      CALL MCOUT
      IF (ITEST.GT.0) WRITE (IU06,*) ' SUB MCOUT DONE'
C
C ----------------------------------------------------------------------
C
C*    6. COMPUTE NEST INFORMATION.
C        -------------------------
C
 6000 CONTINUE
C
C*    6.1 COMPUTE FINE GRID NEST INFORMATION (COMMON FBOUND).
C         ---------------------------------------------------
C
 6100 CONTINUE
      IF (IBOUNF.EQ.1) THEN
         CALL MBOUNF (IU03, IU10, IU20, IFORM, IINPC)
         IF (ITEST.GT.0) WRITE (IU06,*) ' SUB MBOUNF DONE'
      ELSE
         IINPC  = 0
         NBOUNF = 0
      ENDIF
C
C*    6.2 COMPUTE COARSE GRID NEST INFORMATION (COMMON CBOUND).
C         -----------------------------------------------------
C
 6200 CONTINUE
      IF (IBOUNC.EQ.1) THEN
         CALL MBOUNC (IU09, IU19, IFORM)
         IF (ITEST.GT.0) WRITE (IU06,*) ' SUB MBOUNC DONE'
      ELSE
         NBOUNC = 0
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    7. GENERATE CURRENT FIELD.
C        -----------------------
C
 7000 CONTINUE
      IF (IREFRA.EQ.2) THEN
         CALL WAMCUR (IU02, NCC, NRC)
         IF (ITEST.GT.0) WRITE (IU06,*) ' SUB WAMCUR DONE'
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    8. GENERATE AND WRITE COMMON UBUF.
C        -------------------------------
C
 8000 CONTINUE
      DO 8001 IG=1,IGL
         CALL MUBUF (IA1, IG, IU08, IU18, IFORM)
         IF (ITEST.GT.0) THEN
            IF (IG.LE.ITESTB) WRITE (IU06,*)
     1        ' SUB MUBUF DONE FOR BLOCK ',IG
         ENDIF
 8001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    9. OUTPUT OF COMMON BLOCKS.
C        ------------------------
C
 9000 CONTINUE
      CALL OUTCOM (IU07, IU17, IFORM, IREFRA, ML, KL)
      IF (ITEST.GT.0) WRITE (IU06,*) ' SUB OUTCOM DONE'
C
C ----------------------------------------------------------------------
C
C*    10. CONSISTENCY CHECK OF COMPUTED BLOCK PARAMETERS AND
C*        OUTPUT OF NECESSARY DIMENSIONS.
C         --------------------------------------------------
C
 9100 CONTINUE
      CALL CHECK (IREFRA, ML, KL, IINPC)

      STOP
      END
