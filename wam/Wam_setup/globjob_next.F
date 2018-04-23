#  THIS JOB EXECUTES THE GLOBAL VERSION OF THE WAMODEL,
#  DRIVEN BY WINDS REQUESTED FROM MARS.
#  RUNING PREPROC IF PREPRO=.T.
#  CREATING STARTING FILES IF .COLDST.=.T.
#  PRODUCING WAVE HEIGHT AND DIRECTION PLOTS IF PLOT="WAVE"
#  ****WARNING**** FDB FILES ARE CREATED depending on the option
#                  when making the binary (default is fdb=off)
#  SAVES FDB FILES INTO MARS IF ARCHIVE=.T.

# QSUB -r an9209
# QSUB -q rdanexp
# QSUB -lM 20Mw
# QSUB -lT 2000
# QSUB -lt 1200
# QSUB -eo
# QSUB

set +v
set -a
PATH=$PATH:/tmp/emos_sms/bin:/ec/emos_sms/bin:.
cd $TMPDIR

#======================================================================
#   THE USER HAS TO DEFINE THE FOLLOWING VARIABLES !!!
#======================================================================

EXPVER='"xxx"'
EXPER=glo3ice
USEID=$USER

YYBEGA=1992
MMBEGA=09
DDBEGA=13
HHFORC=000
HHANAL=024

PREPRO=.T.
SARCOL=.F.
WAMRUN=.T.
COLDST=.F.
ARCHIVE=.F.

ASSIMILATION=YES   # ( YES or NO )
DOUB=.T.
PER=.T.

PLOT=              # WAVE
PAPER=a3           # SIZE OF THE PAPER  a3  OR a4  OR t for transparencies.
COL=-c             # put nothing after eq-sign for B&W or -c for colour plots.
# to continue a run take begin date = old begin date + hhanal

YMPPATH=$TEMP/$EXPER       # CRAY DISK CATALOG (WORKING CATALOG)
EXPPATH=/$USEID/$EXPER     # IBM TAPE NODE (DATA AND BINARY BACKUP)
SARPAT=/$USER/glo3ice      # PATH FOR SAR DATA (will be extended
                           #   automaticly by /cos/

#======================================================================
#   CREATE THESE IF FIRST RUNE !!!
#======================================================================

mkdir $TEMP/$EXPER
ecfile -p /$USEID/$EXPER add
ecfile -p /$USEID/cos/$EXPER/UWATEST add

#======================================================================
#   END OF USER RELATED CHANGES !!!
#======================================================================

EMOS_VERSION=11   # THIS IS A FIXED VARIABLE
 if [ "$COL" = "-c" ] ; then
   DEVICE=1
 else
   DEVICE=0
 fi

#======================================================================
#   Create WAMINFO file.     resp
#======================================================================

ymd=`substring $YYBEGA$MMBEGA$DDBEGA 3 8`

ymdm1=`newdate $YYBEGA$MMBEGA${DDBEGA}12 +$HHANAL`
ymdm1=`substring $ymdm1 3 8`

ymd1=`newdate 19${ymdm1}12 -24`
ymd1=`substring $ymd1 3 8`

ymdm2=`newdate 19${ymdm1}12 +$HHFORC`
ymdm2=`substring $ymdm2 3 8`

cat > WAMINFO << **EOF**
RUN MODEL FROM ${ymd}1200 TO ${ymdm2}1200
PARAM=10U/10V,
REPRES=GG,
ANALYSIS PERIOD = $HHANAL
FORECAST PERIOD = $HHFORC
WIND TIME STEP IN HOURS =  6
PLOT FROM ${ymdm1}1200 TO ${ymdm2}1200
ASSIMILATION $ASSIMILATION
NEW FORECAST YES
STATUS STORM1= UNFINISHED
STATUS STORM2= UNFINISHED
STATUS STORM3= UNFINISHED
**EOF** 
#
#======================================================================
#
# CRAY DISK CATALOGS (WORKING CATALOGS):
  YMPCY4=/tmp/dax/hrs                                                                                                                                                                                                                                           
#                                                                                                                                                                                                                                                               
# ECFILE NODE FOR CYCLE4 PERMANENT DATA:                                                                                                                                                                                                                        
  CFSCY4=/dax/c90/hrs                                                                                                                                                                                                                                           
#                                                                                                                                                                                                                                                               
  _WAVELIB=$YMPCY4/wavelib.a                                                                                                                                                                                                                                    
  _E_WAVELIB=$CFSCY4/wavelib.a                                                                                                                                                                                                                                  
#                                                                                                                                                                                                                                                               
  _ALTLIB=/tmp/emos_sms/lib/next/altlib.a                                                                                                                                                                                                                       
  _E_ALTLIB=/wavearch/0011/ers1/alt/global/altlib.a                                                                                                                                                                                                             
                                                                                                                                                                                                                                                                
#======================================================================                                                                                                                                                                                         
#  Load executable preproc.                                                                                                                                                                                                                                     
#======================================================================                                                                                                                                                                                         
#                                                                                                                                                                                                                                                               
emos_ecfile -p $_E_WAVELIB  -en get $_wavelib 
#                                                                                                                                                                                                                 
#======================================================================
#  General grid information and model constants are fetched(PREPRO=.F.)
#                                                or created(PREPRO=.T.)
#======================================================================

if [ "$PREPRO" = ".F." ] ; then
  #
  #   INPUT FILE
  #
  #
  #         GRID
  #
  test -r $YMPPATH/gridglou
  if [ $? != 0 ] ; then
  ecfile -p $EXPPATH/gridglou -eG get $YMPPATH/gridglou
    if [ $? != 0 ] ; then
      ls -lsa
      echo "\n\n\t\tPREPROC file gridglou was not found\n\n"
      echo "\n\n\t\tPROGRAM IS TERMINATED\n\n"
      exit 1
    fi
  fi
  cp $YMPPATH/gridglou GRID

  #
  #         UBUF
  #
  test -r $YMPPATH/ubufglou
  if [ $? != 0 ] ; then
    ecfile -p $EXPPATH/ubufglou -eG get $YMPPATH/ubufglou
    if [ $? != 0 ] ; then
      ls -lsa
      echo "\n\n\t\tPREPROC file ubufglou was not found\n\n"
      echo "\n\n\t\tPROGRAM IS TERMINATED\n\n"
      exit 2
    fi
  fi
  cp $YMPPATH/ubufglou UBUF

else
    segldr -ePREPROC -l $_WAVELIB $ECLIB $EMOSLIB -o preproc                                                                                                                                                                                                    
    if [ $? != 0 ] ; then
      echo "\n\n\t\tPREPROC NOT LOADED  \n\n"
      exit 4
    fi
  fi
  #
  #   INPUT FILES.
  #
  #      TOPOGRAFIC DATA.
  #
  test -r $YMPPATH/topoglo2
  if [ $? != 0 ] ; then
    ecfile -p /waj/cos/wamoper/topoglo2 -enG get $YMPPATH/topoglo2
    if [ $? != 0 ] ; then
      ls -lsa
      echo "\n\n\t\tTOTPOGLO2 absolute preproca is not found\n\n"
      echo "\n\n\t\tPROGRAM IS TERMINATED\n\n"
      exit 5
    fi
  fi
  cp $YMPPATH/topoglo2 fort.1

  #
  #   RUN PREPROC:
  #
  #
  cat > procin <<EOD
C
C       USER INPUT OF PREPROC.           CREATED 03/11/92
C       ======================
C
C       SET-UP OF GLOBAL 3 DEGREE GRID.
C
C        USER INFORMATION IS BEING READ BY SUB. UIPREP
C        WITH THE PRESUMPTIONS THAT:
C         1. EVERY LINE STARTING WITH 'C' IS A COMMENT LINE
C         2. VALUES ARE PUT IN BELOW POSITIONS INDICATED WITH '-'
C            (RIGHT-JUSTIFIED)
C
C ----------------------------------------------------------------------
C
C           **** HEADER OF MODEL RUN   ****
C
C     HEADER  CHARACTER*70  TEXT TO BE PRINTED
C
CHEADER
C---------+---------+---------+---------+---------+---------+---------+
 GLOBAL 3 DEGREE GRID                             PREPROC VERSION 4.0
C
C           **** FREQUENCY AND DIRECTION GRID ****
C
C     ML        INTEGER*5    NUMBER OF FREQUENCIES.
C     KL        INTEGER*5    NUMBER OF DIRECTIONS.
C     FR(1)     REAL*10      LOWEST FREQUENCY IN HERTZ.
C
C   ML    KL      FR(1)
C----- ----- ----------
    25    12  .04177248
C
C           **** OUTPUT GRID DEFINIONS ****
C
C     XDELLA    REAL       GRID INCREMENT FOR LATITUDE (DEGREE)
C     XDELLO    REAL       GRID INCREMENT FOR LONGITUDE (DEGREE)
C     AMOSOP    REAL       MOST SOUTHERN LATITUDE IN GRID MATRIX
C     AMONOP    REAL       MOST NORTHERN LATITUDE IN GRID MATRIX
C     AMOWEP    REAL       MOST WESTERN LONGITUDE IN GRID MATRIX
C     AMOEAP    REAL       MOST EASTERN LONGITUDE IN GRID MATRIX
C
C    XDELLA     XDELLO     AMOSOP     AMONOP     AMOWEP     AMOEAP
C---------- ---------- ---------- ---------- ---------- ----------
      3.000      3.000    -78.000     78.000      0.000    357.000
C
C           **** AREAS TO BE CHANGED IN GRID ****
C
C           SPECIFY A LIST WITH UP TO 80 AREAS.
C
C     XOUTS     REAL       SOUTHERN-MOST LATITUDE (DEGREE).
C     XOUTN     REAL       NORTHERN-MOST LATITUDE (DEGREE).
C     XOUTW     REAL       WESTERN-MOST LONGITUDE (DEGREE).
C     XOUTE     REAL       EASTERN-MOST LONGITUDE (DEGREE).
C     NOUTD     INTEGER    NEW DEPTH IN METRES (-999 IS LAND)
C
C     XOUTS      XOUTN      XOUTW      XOUTE      NOUTD
C---------- ---------- ---------- ---------- ----------
C
C   THE LIST HAS TO BE CLOSED BY  'END'
C
 END
C
C           **** OUTPUT LOCATIONS OF SPECTRA ****
C
C           SPECIFY A LIST WITH UP TO 100 SIDES.
C
C     OUTLAT    REAL       LATITUDE (DEGREE).
C     OUTLONG   REAL       LONGITUDE (DEGREE).
C
C    OUTLAT    OUTLONG
C---------- ----------
     57.000   -148.000    46001
     43.000   -130.000    46002
     52.000   -156.000    46003
     51.000   -136.000    46004
     46.000   -131.000    46005
     41.000   -138.000    46006
     23.000   -163.000    51001
     18.000   -158.000    51002
     20.000   -161.000    51003
     18.000   -153.000    51004
     40.000    -70.000    44008
     41.000    -67.000    44011
     35.000    -73.000    41001
     32.000    -75.000    41002
     30.000    -77.000    41006
     60.000      0.000    ERS-1 NORWEGEN
     60.000      3.000    ERS-1 NORWEGEN
     63.000      0.000    ERS-1 NORWEGEN
     63.000      3.000    ERS-1 NORWEGEN
     63.000      6.000    ERS-1 NORWEGEN
     66.000      0.000    ERS-1 NORWEGEN
     66.000      3.000    ERS-1 NORWEGEN
     66.000      6.000    ERS-1 NORWEGEN
     66.000      9.000    ERS-1 NORWEGEN
     45.000    -54.000    ERS-1 CANADA
     45.000    -51.000    ERS-1 CANADA
     45.000    -48.000    ERS-1 CANADA
     48.000    -51.000    ERS-1 CANADA
     48.000    -48.000    ERS-1 CANADA
C
C   THE LIST HAS TO BE CLOSED BY  'END'
C
 END
C
C           **** MODEL OPTIONS ****
C
C     IFORM     INTEGER    = 1   UNFORMATED WRITE ONLY.
C                          = 2   FORMATED WRITE ONLY.
C                          OTHERWISE BOTH OUTPUTS.
C     IREFRA    INTEGER    = 2   A CURRENT FIELD IS PROCESSED.
C                          OTHERWISE NOTHING IS DONE.
C     ITEST     INTEGER    TEST OUTPUT LEVEL..
C                          > 0 FOR OUPUT UP TO SUBROUTINE LEVEL
C                          = 0  NO TEST OUTPUT
C     ITESTB    INTEGER    TEST OUTPUT LEVEL BLOCKS  IF ITEST > 0.
C                          > 0 FOR BLOCK LEVEL TEST OUTPUT
C
C IFORM  IREFRA   ITEST  ITESTB
C------  ------  ------  ------
      1       0       0       4
C
C           **** NESTED GRIDS ****
C
C     IBOUNC  INTEGER   FLAG FOR THE COURSE GRID
C                       = 1 INFORMATION FOR A NESTED (FINE)
C                           GRID WILL BE GENERATED.
C     IBOUNF  INTEGER   FLAG FOR THE FINE GRID
C                       = 1  THIS IS A FINE GRID RUN, BOUNDARY
C                            INFORMATION IS EXPECTED FROM A
C                            PREVIOUS COARSE GRID RUN.
C     AMOSOC  REAL      SOUTHERN-MOST LONGITUDE OF NEST (DEGREE)
C     AMONOC  REAL      NORTHERN-MOST LONGITUDE OF NEST (DEGREE)
C     AMOWEC  REAL      WESTERN-MOST LATITUDE OF NEST (DEGREE)
C     AMOEAC  REAL      EASTERN-MOST LATITUDE OF NEST (DEGREE)
C                       (THESE NEST BOUNDARIES HAVE TO BE GIVEN
C                        IF IBOUNC = 1)
C
CIBOUNC IBOUNF     AMOSOC     AMONOC     AMOWEC     AMOEAC
C------ ------ ---------- ---------- ---------- ----------
      0      0      0.000      0.000      0.000      0.000
C
C ----------------------------------------------------------------------
C        END OF USER INPUT OF PREPROC
C ----------------------------------------------------------------------
EOD
  cat > fort.55 <<\EOFF
 SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
 SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
 SSSSSSSSSSSSSSSSSSLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
 SSSSSSLLLLLSSSSSSSSSLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLSSSSLLLLLLLLLLLSSSSSSS
 SSSSLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLSSLLLLLLLSSSSSSSSSSS
 SSSLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLSSLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLSSSSLLLLSSSSSSSSSSSSS
 SSLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLSSSSSSSSLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLSSSSSSSLSSSSSSSSSSSSSS
 SSSSLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLSSSSSLLSSSSSSSSSSSSSSLSSSSSSSLLLLLLLLLLLLLLLLLLLLLLLLSSSSSSSSSSSSSSSSSSLL
 LSSSLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLSSSSSSLLSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLLLLLLLLLLLLLLLLLSSSSSSSSSSSSSSSSLLL
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLSSSSLSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLLLLLLLLLLLLLLLLLSSSSSSSSSSSSSSSSSL
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLLLLLLLLLLLLLLLLLSSSSSSSSSSSSSSSSL
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLLLLLLLLLLLLLLSSSSSSSSSSSSSSSSSSSS
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLSSSLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLLLLLLLLLLLSSSSSSSSSSSSSSSSSSSSLLL
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLSSSSLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLLLLLLLLLLSSSSSSSSSSSSSSSSSSSSSLLL
 SLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLSSLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLLLLLLLLSSSSSSSSSSSSSSSSSSSSSSSLL
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLSLLLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLLLLLLSSSSSSSSSSSSSSSSSSSSSSSSLL
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLLLSSSSSSSSSSSSSSSSSSSSSSSSLLL
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLSSSSSLSSSSSSSSSSSSSSSSSSSSSSLLLL
 LLLLLLLLLLLLLLLLLLLSSSSLLLLLLLLLLLLLLLLLLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLL
 LLLLLLLLLLLLLLLLLLLLSSSLLLLLLLSLLLLLLLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLSSSLSSLSSSSSSSSSSSSSSSSSSSSLLLLL
 LLLLLLLLLLLLLLLLLLLSSSSSSLLLSSSSLLLLLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLSSSSSLSSSSSSSSSSSSSSSSSSLLLLL
 LLLLLLLLLLLLLLLLLLSSSSSSSLLSSSSSSLLLLSSSLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLSSSSSSSSSSSSSSSSSSSSSSLLLLL
 LLLLLLLLLLLLLLLLSSSSSSSSSLLSSSSSSLLLLSSSSLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLSSSSSSSSSSSSSSSSSSSSSSLLLLL
 LLLLLLLLLLLLLLLLLSSSSSSSSSLSSSSSSLSLSSSSSLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLSSSSSSSSSSSSSSSSLLLL
 LLLLLLLLLLLLLLLLLSSSSSSSSSSLSSSSSLSSSSSLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLSSSSSSSSSSSSSSSSLLL
 SSSSLLLLLLLLLLLLSSSSSSSSSSSSSSSSSLLSSSLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLSSSSSSSSSSSSSSSSS
 SSSLLLLLLLLLLLLSSSSSSSSSSSSSSSSSSLLSSLLLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLLLSSSSSSSSSSSSSSSS
 SSSSLLLLLLLLLLSSSSSSSSSSSSSSSSSSSSLLSLLSLSSSLSLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLLLLLLSSSSSSSSSSSSSS
 SSSSLLLLLLLLLSSSSSSSSSSSSSSSSSSSSSSLSSSSSSSSSSLLLSLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLLLLLLLLLSSSSSSSSSSS
 SSSSSLLLLLLLLSSSSSSSSSSSSSSSSSSSSSSSSSLLLSSSSSSLSLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLLLLLLLLSSSSSSSSSSS
 SSSSSLLLLLLLLLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLLLLLLLSSSSSSSSSSSS
 SSSSSLLLLLLLLLSSLSSSSSSSSSSSSSSSSSSSSSSSSSLLLLSLLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLLLLLSSSSSSSSSSSSS
 SSSSLLLLLLLLLSSLLSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLLLLSSSSSSSSSSSSS
 SSSSSLLLLLLLSSSLLSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLLLLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLLLSSSSSSSSSSSSS
 SSSSSLLLLLLLSSSLSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLLLLLLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLSSSSSSSSSSSSSSS
 SSSSSLLLLLLSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLLLLLLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLSSSSSSSSSSSSSSSS
 SSSSSSLLLLLSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLLLLLLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLSSSSSSSSSSSSSSSS
 SSSSSSLLLLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLSSSLLLLLLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLLLSSSSSSSSSSSSSSSSS
 SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLSSSSSSSSLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLSSSSSSSSSSSSSSSSSSS
 SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLLSSSSSSSSSSSSSSSSSSS
 SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLSSSSSSSLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLSSSSSSSSSSSSSSSSSSSSS
 SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLLSSSSSSSSSSSSSSSSSSSSS
 SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLLSSSSSSSSSSSSSSSSSSSSSS
 SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLSSSSSSSSSSSSSSSSSSSSSS
 SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLLSSSSSSSSSSSSSSSSSSSSSS
 SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
 SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
 SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
 SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLSSSSSSSSSSSSSSSSSSSS
 SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSLSSSSSSSSSSSSSSSSSSSSS
 SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
 SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
 SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
EOFF

  preproc < procin > procout ||  { ls -l ; asa procout ; exit 6 ; }
  asa procout

  #
  #    SAVE BACK-UP FILES IN ECFILE
  #

  mv fort.7 GRID || exit 7
  mv fort.8 UBUF || exit 8
  cp GRID $YMPPATH/gridglou   # because of script inconsistence

  ecfile -p $EXPPATH/gridglou -hO -uD save GRID   || exit 9
  ecfile -p $EXPPATH/ubufglou -hO -uD save UBUF   || exit 10

  echo "\n\n\t\tPREPROC ALL DONE\n\n"

  #=====================================================================
  #  PREPROC FILES GENERATED AND SAVED.
  #  END OF JOB PREPROC.
  #=====================================================================
  rm fort.55
fi 
#
#======================================================================
#
# TEMPORARY EXIT
#
#======================================================================
#
exit
#
#======================================================================
#   EXECUTE PROGRAM PREUWA IF SARCOL=.T.
#======================================================================

if [ "$SARCOL" = ".T." ] ; then

  # ====================================================================
  #  Copy WAMODEL PARAMETER FILE
  # ====================================================================

  cp $SOURCE/Input_Files/PARWAM PARWAM

  test -r $YMPPATH/uwapre.x
  if [ $? != 0 ] ; then
    cft77 -esz -Rabc -b uwapre.o $SOURCE/Sar/uwapre_6.f
    if [ $? != 0 ] ; then
      ls -lsa
      echo "\n\n\t\tUWAPRE.F NOT COMPILED\n\n"
      echo "\n\n\t\tPROGRAM IS TERMINATED\n\n"
      exit 11
    fi
    segldr uwapre.o -l $WAVELIB  $EMOSLIB $ECLIB -o uwapre.x
    if [ $? != 0 ] ; then
      ls -lsa
      echo "\n\n\t\tPREUWA.X NOT LOADED  \n\n"
      echo "\n\n\t\tPROGRAM IS TERMINATED\n\n"
      exit 11
    fi
    cp uwapre.x $YMPPATH
  else
    cp $YMPPATH/uwapre.x uwapre.x
  fi
  #
  #   RUN PREUWA:
  #
  DATE=${ymd}1800

  while [ $DATE -le ${ymdm1}1200  ] ; do
      if [ -s $YMPPATH/COL$DATE ] ; then
        cp $YMPPATH/COL$DATE COL$DATE
      else
        ecfile -p $EXPPATH/col${DATE} get $YMPPATH/COL${DATE}
        if [ $? != 0 ] ; then
           ecfile -p /oparch/ers1/uwax`da2op` -en get UWA$DATE
           uwapre.x
           if [ $? != 0 ] ; then
             ls -lsa
             echo "\n\n\t\tPREUWA FAILED\n\n"
             echo "\n\n\t\tPROGRAM IS TERMINATED\n\n"
             debug -s uwapre.x -B
             exit 12
           fi
           ecfile -p $EXPPATH/col${DATE} -uD -hO save  $YMPPATH/COL${DATE}
           cp COL$DATE $YMPPATH/COL$DATE
        else
           cp $YMPPATH/COL$DATE COL$DATE
        fi
      fi
      NDAT=`newdate 19$DATE +6`
      DATE=`substring $NDAT 3 10`00
  done

  #====================================================================
  #  SAR COLLOCATION FILES GENERATED AND SAVED.
  #  END OF JOB PREUWA.
  #====================================================================
fi


#======================================================================
#   EXECUTE PROGRAM URAPRE IF ASSIMILATION = YES
#======================================================================

if [ "$ASSIMILATION" = "YES" ] ; then
  test -r $BIN/urapre.x
  if [ $? != 0 ] ; then
    cft77 -esz $SOURCE/Work/urapre.f
    segldr urapre.o -l $MAGLIB $EMOSLIB $ECLIB -o $BIN/urapre.x
    if [ $? != 0 ] ; then
      asa urapre.l
      echo "\n\n\t\t  Load of urapre failed  \n\n"
      rm urapre.[ol]
      exit 2
    fi
    rm urapre.[ol]
  fi
  cp $BIN/urapre.x urapre.x

  # ====================================================================
  #  Copy WAMODEL PARAMETER FILE IF NOT ALREADY DONE.
  # ====================================================================

  [ -s PARWAM ] || cp $SOURCE/Input_Files/PARWAM PARWAM

  # ====================================================================
  #
  #  COPY USER INPUT FOR PREURA TO 'preinp'
  #
  #
  cat > preinp <<\***END***
C key  sub-routine
C----     --------
CSUVI     uraqwf
CTEST
***END***

  # ====================================================================
  #
  #  COPY PARAMETER NAMELIST FOR PREURA TO 'PARQC'
  #
  cat > PARQC   <<*EOF*
C ----------------------------------------------------------------------
C
C     Quality control parameters.
C     ---------------------------
C
        \$parqc
        pproz  =   0.5,     :  factor for variance test
        ppspi  =   3.0,     :  factor for spike test
        rmsf   =   2.0,     :  factor for rmsh
        rmsm   =   0.5,     :  maximum value for rms of hs
        hsdifm =   2.0,     :  maximum difference between two wave
                            :     heights before and behind a gap
        dtmax  =   3.0,     :  maximum time difference in seconds.
        hmin   =   0.441,   :  minimum waveheight
        hmax   =  17.479,   :  maximum waveheight
        nmax   =  30,       :  maximum number of records taken for
                            :     averaging
        nmin   =  20,       :  minimum number of records taken for
                            :     averaging
        rlamin = -90.0,     :  southern boundary
        rlamax =  90.0      :  northern boundary
        \$END
*EOF*

  #
  #   RUN URAPRE if there is no flagged file on $YMPPATH or $EXPPATH
  #
  DATE=${ymd}1800
  while [ $DATE -le ${ymdm1}1200  ]
    do
      if [ -s $YMPPATH/RFL$DATE ] ; then
        cp $YMPPATH/RFL$DATE RFL$DATE
      else
        ecfile -p $EXPPATH/rfl${DATE} -enG get RFL${DATE}
        if [ $? != 0 ] ; then
           urapre.x  > LOG < preinp  \
            || { debug -B -s urapre.x | cat ; asa LOG ; exit ; }
           ls -l
           asa LOG
           if [ -s RFL${DATE} ] ; then
              ecfile -p $EXPPATH/rfl${DATE} -hO -uD save RFL${DATE}
              cp RFL${DATE}  $YMPPATH/RFL$DATE
           fi
        fi
      fi
      NDAT=`newdate 19$DATE +6`
      DATE=`substring $NDAT 3 10`00
    done

  #====================================================================
  #
  #  RADAR ALTIMETER DATA PREPARED FOR USE BY THE WAModel ASSIMILATION
  #
  #====================================================================
fi


#======================================================================
#  START RUNNING WAM IF (WAMRUN=.T.)
#======================================================================

if [ "$WAMRUN" = ".T." ] ; then

  #====================================================================
  #   Prepare MARS files.
  #====================================================================

  ls -lt
cp $YMPPATH/m${ymd}a${HHANAL}f${HHFORC} fort.63
if [ $? != 0 ] ; then
  ecfile -enG -p $EXPPATH/m${ymd}a${HHANAL}f${HHFORC} get fort.63
  #
  #  If mars file not there prepare mars request
  #
  if [ $? != 0 ] ; then
    if [ "$ymd" = "$ymd1" ] ; then
      if [ "$HHFORC" = "000" ] ;  then
        cat > marsin << **EOF**
ret,
    date=$ymd,
    time=00/06/12/18,
    param=10u/10v,
    repres=gg,
    type=ia,
    levtype=sfc,
    levelist=off,
    format=packed,
    target="fort.63"
ret,
    date=$ymdm1,
    time=00/06/12
end
**EOF**

      else
        cat > marsin << **EOF**
ret,
    date=$ymd,
    time=00/06/12/18,
    param=10u/10v,
    repres=gg,
    type=ia,
    levtype=sfc,
    levelist=off,
    format=packed,
    target="fort.63"
ret,
    date=$ymdm1,
    time=00/06/12
ret,
    type=fc,
    date=$ymdm1,
    time=12,
    step=6/to/${HHFORC}/by/06
end
**EOF**

      fi
    else
      if [ "$HHFORC" = "000" ] ;  then
        cat > marsin << **EOF**
ret,
    date=$ymd/to/$ymd1,
    time=00/06/12/18,
    param=10u/10v,
    repres=gg,
    type=ia,
    levtype=sfc,
    levelist=off,
    format=packed,
    target="fort.63"
ret,
    date=$ymdm1,
    time=00/06/12
end
**EOF**

      else
        cat > marsin << **EOF**
ret,
    date=$ymd/to/$ymd1,
    time=00/06/12/18,
    param=10u/10v,
    repres=gg,
    type=ia,
    levtype=sfc,
    levelist=off,
    format=packed,
    target="fort.63"
ret,
    date=$ymdm1,
    time=00/06/12
ret,
    type=fc,
    date=$ymdm1,
    time=12,
    step=6/to/${HHFORC}/by/06
end
**EOF**

      fi
    fi
    mars < marsin || exit 13
    ecfile -p $EXPPATH/m${ymd}a${HHANAL}f${HHFORC} -hO -r20 save fort.63
    cp fort.63 $YMPPATH/m${ymd}a${HHANAL}f${HHFORC}
  fi
fi
  #=====================================================================
  #  MARS FILES PREPARED
  #=====================================================================

  #=====================================================================
  #   WAVE MODEL execution:
  #=====================================================================

  #=====================================================================
  #   CREATE STARTING FILES IF ASKED FOR (COLDST=.T.)
  #   OR FETCH THEM.
  #=====================================================================

  if [ "$COLDST" = ".T." ] ; then
    test -r $YMPPATH/preseta
    if [ $? != 0 ] ; then
      ecfile -p $EXPPATH/preset.o  -enG get preset.o
      if [ $? != 0 ] ; then
        ls -lsa
        echo "\n\n\t\tPRESET:  BINARY NOT FOUND \n\n"
        echo "\n\n\t\tJOB TERMINATED\n\n"
        exit 14
      fi
      segldr preset.o -g -V      \
           -l /ec/stg/wamcrlib2  \
              $ECLIB             \
              $EMOSLIB           \
           -o  $YMPPATH/preseta
      if [ $? != 0 ] ; then
        echo "\n\n\t\t  PRESET NOT LOADED \n\n"
        exit 15
      fi
    fi
    cp $YMPPATH/preseta preseta

    cat > presetin << EOD
C
C       USER INPUT OF PRESET.            CREATED 19/4/91
C       =====================
C
C        USER INFORMATION IS BEING READ WITH THE PRESUMPTIONS THAT:
C         1. EVERY LINE STARTING WITH 'C' IS A COMMENT LINE
C         2. VALUES ARE PUT IN BELOW POSITIONS INDICATED WITH '-'
C            (RIGHT-JUSTIFIED, BUT CHARACTER LEFT-JUSTIFIED)
C
C ----------------------------------------------------------------------
C
C           **** HEADER OF MODEL RUN   ****
C
C     HEADER  CHARACTER*70  TEXT TO BE PRINTED
C
CHEADER
C---------+---------+---------+---------+---------+---------+---------+
 GLOBAL        WAM MODEL INITIALISATION            PRESET VERSION 4.0
C
C           **** OPTION FOR INITIAL SPECTRA   ****
C
C     IOPTI   = 0 WIND INDEPENDENT INITIAL VALUES.
C             = 1 WIND DEPENDENT INITIAL VALUES AND
C                 ENERGY EQUAL ZERO IF WINDSPEED IS ZERO
C             = 2 WIND DEPENDENT INITIAL VALUES AND
C                 ENERGY COMPUTED FROM GIVEN PARAMETERS IF
C                 WINDSPEED IS ZERO.
C     ITEST   = 0 NO TEST OUTPUT
C             > 0 FOR SUBROUTINE LEVEL
C     ITESTB  > 0 FOR BLOCK LEVEL TEST OUTPUT
C
C  IOPTI    ITEST   ITESTB
C-------  -------   ------
       2        0        4
C
C           **** PARAMETER OF INITIAL SPECTRUM   ****
C
C     ALPHA   = PHILLIPS' PARAMETER  (NOT USED IF IOPTI = 1)
C     FM      = PEAK FREQUENCY (HZ) AND/OR MAXIMUM FREQUENCY
C     GAMMA   = OVERSHOOT FACTOR
C     SIGMAA  = LEFT PEAK WIDTH
C     SIGMAB  = RIGHT PEAK WIDTH
C     THETAQ  = WAVE DIRECTION (DEG) (NOT USED IF IOPTI = 1)
C     FETCH   = FETCH IN METRES (IF ZERO THEN 0.5 OF THE
C               LATITUDE INCREMENT IS USED.)
C
C     ALPHA          FM       GAMMA      SIGMAA      SIGMAB      THETAQ
C----------  ----------  ----------  ----------  ----------  ----------
    0.01800     0.20000     3.00000     0.07000     0.09000     0.00000
C
C     FETCH
C----------
     50000.
C
C           **** FILE NAMES ****
C
C     UID      CHARACTER*3   USERID
C     RID      CHARACTER*3   RUN IDENTIFIER
C     PATH     CHARACTER*60  PATH NAME OF FILES
C
CUID  RID  PATH
C---  ---  ---------+---------+---------+---------+---------+---------+
 $USEID  GLO  $EXPER
C
C           **** WIND INPUT PARAMETERS ****
C
C     IDATEA   INTERGER*10   START DATE OF RUN  (YYMMDDHHMM)
C     IDELWI   INTEGER*7     TIMESTEP ON INPUT WIND FILE (SECONDS)
C     U        CHARACTER*1   TIME UNIT  S SECONDS  OR  H HOURS
C
C    IDATEA   IDELWI U
C----------  ------- -
 7809060600    21600 S
C!!!!!!!!!!  !!!!!!! !
C
C ----------------------------------------------------------------------
C                END OF USER INPUT PROG. PRESET.
C ----------------------------------------------------------------------
EOD

    preseta < presetin > presetout
    if [ $? != 0 ] ; then
#     flowview -L -u -l
      echo "\n\n\t\tPRESET FAILED\n\n"
      exit 16
    else
#     flowview -L -u -l
      echo "\n\n\t\tPRESET ALL DONE\n\n"
    fi

    ls -l

    ecfile -p /$USEID/cos/$EXPER/lawianal -p $EXPPATH/lawi${ymd}a -hO copy
    ecfile -p /$USEID/cos/$EXPER/blspanal -p $EXPPATH/blsp${ymd}a -hO copy
    ecfile -p /$USEID/cos/$EXPER/slatanal -p $EXPPATH/slat${ymd}a -hO copy
      ecfile -p /$USEID/$EXPER/lawi${ymd}a get LAWI
      ecfile -p /$USEID/$EXPER/blsp${ymd}a get BLSP1
      ecfile -p /$USEID/$EXPER/slat${ymd}a get SLAT1
      cp LAWI   $YMPPATH/lawi${ymd}a
      cp BLSP1  $YMPPATH/blsp${ymd}a
      cp SLAT1  $YMPPATH/slat${ymd}a

    #===================================================================
    #  START FIELDS GENERATED AND SAVED.
    #  END OF JOB PRESET.
    #===================================================================

  else

    cp $YMPPATH/lawi${ymd}a LAWI  || \
      ecfile -p /$USEID/$EXPER/lawi${ymd}a get LAWI
    cp $YMPPATH/blsp${ymd}a BLSP1 || \
      ecfile -p /$USEID/$EXPER/blsp${ymd}a get BLSP1
    cp $YMPPATH/slat${ymd}a SLAT1 || \
      ecfile -p /$USEID/$EXPER/slat${ymd}a get SLAT1
  fi

  #=====================================================================
  #   RUN THE MODEL WITH the following USER INPUT.
  #=====================================================================

  test -r $YMPPATH/wama.next
  if [ $? != 0 ] ; then
    cft77 -esz -Rabc $SOURCE/wavemdl.f
    if [ $? != 0 ] ; then
      ls -lsa
      echo "\n\n\t\t  wavemdl.f not compiled        \n\n"
      echo "\n\n\t\t  JOB TERMINATED           \n\n"
      exit 17
    fi
    segldr wavemdl.o                 \
          -l /ec/stg/wamcrlib2       \
             $ECLIB                  \
             $EMOSLIB                \
             /tmp/emos_sms/lib/fdb.a \
             /usr/lib/libsort.a      \
           -g -V -o wama.next
    if [ $? != 0 ] ; then
      echo "\n\n\t\t WAVEMODEL NOT LOADED   !! \n\n"
      exit 18
    fi
    cp wama.next $YMPPATH/wama.next
  else
    cp $YMPPATH/wama.next wama.next
  fi

  cat > WAMINFO << **EOF**
RUN MODEL FROM ${ymd}1200 TO ${ymdm2}1200
PARAM=10U/10V,
REPRES=GG,
ANALYSIS PERIOD = $HHANAL
FORECAST PERIOD = $HHFORC
WIND TIME STEP IN HOURS =  6
PLOT FROM ${ymdm1}1200 TO ${ymdm2}1200
ASSIMILATION $ASSIMILATION
NEW FORECAST YES
STATUS STORM1= UNFINISHED
STATUS STORM2= UNFINISHED
STATUS STORM3= UNFINISHED
**EOF**

  cat > wamin << EOD
C
C       USER INPUT OF CHIEF.             CREATED 19/4/91
C       ====================
C
C        USER INFORMATION IS BEING READ WITH THE PRESUMPTIONS THAT:
C         1. EVERY LINE STARTING WITH 'C' IS A COMMENT LINE
C         2. VALUES ARE PUT IN BELOW POSITIONS INDICATED WITH '-'
C            (RIGHT-JUSTIFIED, BUT CHARACTER LEFT-JUSTIFIED)
C
C ----------------------------------------------------------------------
C
C         **** HEADER OF MODEL RUN   ****
C
C     HEADER  CHARACTER*70  TEXT TO BE PRINTED
C
CHEADER
C1--------+---------+---------+---------+---------+---------+---------+
 GLOBAL WAVE MODEL RUN                             MODEL VERSION 4.0
C
C           **** PERIOD OF MODEL RUN   ****
C
C     IDATEA   INTERGER*10   START DATE OF  RUN  (YYMMDDHHMM)
C     IDATEE   INTERGER*10   END   DATE OF  RUN  (YYMMDDHHMM)
C     IDATEF   INTERGER*10   END   ANALYSIS RUN  (YYMMDDHHMM)
C
C    IDATEA      IDATEE      IDATEF
C2---------  ----------  ----------
 7809060600  7809080600  7809080600
C!!!!!!!!!!  !!!!!!!!!!  !!!!!!!!!!
C
C           **** MODEL TIME STEPS ****
C
C     IDELPRO  INTERGER*7    PROPAGATION TIMESTEP
C     IDELT    INTERGER*7    SOURCE TIME STEP
C     IDELWO   INTERGER*7    OUTPUT WIND TIMESTEP
C     IDELWI   INTERGER*7    TIMESTEP ON INPUT WIND FILE
C     U        CHARACTER*1   TIME UNIT  S SECONDS  OR  H HOURS
C
C           --------   --------------------------------------------
C           DEMAND :   ALL RATIO'S ARE N OR 1/N AND N IS AN INTEGER
C           --------   --------------------------------------------
C
CIDELPRO U    IDELT U   IDELWO U   IDELWI U
C3------ -  ------- -  ------- -  ------- -
       1 H     1200 S        6 H        6 H
C                                 !!!!!!!!!
C
C           **** OUTPUT TIME IN FIXED INTERVALS ****
C
C     ------------------------------------------------------------
C     DEMAND :   ALL OUTPUT TIMESTEPS MUST BE MULTIPLE OF IDELPRO.
C     ------------------------------------------------------------
C
C     IDELINT  INTERGER*7    INTEGRATED PARAMETER OF TOTAL SEA
C     IDELINS  INTERGER*7    INTEGRATED PARAMETER OF SWELL
C     IDELAGI  INTERGER*7    INTEGRATED PARAMETER OF AG SWELL
C     IDELSPT  INTERGER*7    SPECTRA OF TOTAL SEA
C     IDELSPS  INTERGER*7    SPECTRA OF SWELL
C     IDELAGS  INTERGER*7    SPECTRA OF AG SWELL
C     IDELRES  INTERGER*7    OUTPUT FILES AND RESTART FILES ARE SAVED
C     U        CHARACTER*1   TIME UNIT  S SECONDS  OR  H HOURS
C
CIDELINT U  IDELINS U  IDELAGI U  IDELSPT U  IDELSPS U  IDELAGS U  IDELRES U
C4------ -  ------- -  ------- -  ------- -  ------- -  ------- -  ------- -
       6 H        6 H        6 H        6 H        6 H        6 H       24 H
C
C           **** OUTPUT TIME AT SPECIFIED TIMES ****
C
C     ----------------------------------------------------------------
C     DEMAND :   ALL OUTPUT TIMES MUST BE AT THE END OF A PROPAGATION.
C     ----------------------------------------------------------------
C   SPECIFY A LIST WITH UP TO 100 DATES (YYMMDDHHMM) AND CLOSE THE
C   LIST BY 'END       '. IF THE LIST INCLUDES ANYTHING ELSE THAN 'END'
C   ALL OUTPUT TIME STEPS SPECIFIED ABOVE ARE IGNORED.
C
CYYMMDDHHMM  YYMMDDHHMM  YYMMDDHHMM  YYMMDDHHMM  YYMMDDHHMM  YYMMDDHHMM
C5---------  ----------  ----------  ----------  ----------  ----------
 END
C           **** SELECTION OF OUTPUT DATA ****
C
C     P   CHARACTER*1   PRINTER OUTPUT..
C          Y  DATA ARE PRINTED, OTHERWISE NO OUTPUT
C     F   CHARACTER*1   FILE OUTPUT..
C          Y  DATA ARE WRITTEN TO FILE, OTHERWISE NO OUTPUT
C     G   CHARACTER*1   GRIB PACKING..
C          Y  DATA ARE PACKED AND WRITTEN TO FILE, OTHERWISE NO PACKING
C
CP F G                                 P F G
C6 - -                                 - - -
 N NYY SIGNIFICANT WAVE HEIGHT         N NYY MEAN WAVE DIRECTION
 N NYY WAVE MEAN FREQUENCY             N NYN FRICTION VELOCITY
 N NYY WIND DIRECTION                  N NYY WAVE PEAK FREQUENCY
 N NYY DRAG COEFFICIENT                N NYN NORMALISED WAVE STRESS
 N N N SWELL WAVE HEIGHT               N N Y SWELL DIRECTION
 N N N WIND WAVE DIRECTION             N N Y SWELL MEAN FREQUENCY
 N NYY SPECTRA OF TOTAL SEA            N N N SWELL SPECTRA
C
C
C           **** SELECTION OF OUTPUT OF WAVE TRAIN SEPERATION ****
C
C     P   CHARACTER*1   PRINTER OUTPUT..
C          Y  DATA ARE PRINTED, OTHERWISE NO OUTPUT
C     F   CHARACTER*1   FILE OUTPUT..
C          Y  DATA ARE WRITTEN TO FILE, OTHERWISE NO OUTPUT
C     G   CHARACTER*1   GRIB PACKING..
C          Y  DATA ARE PACKED AND WRITTEN TO FILE, OTHERWISE NO PACKING
C
CP F G                                 P F G
C7 - -                                 - - -
 N NYY SIGNIF. HEIGHT OF WIND WAVES    N NYY MEAN DIRECTION OF WIND WAVES
 N NYY MEAN PERIOD OF WIND WAVES       N NYY SIGNIF. HEIGHT OF PRIMARY SWELL
 N NYY MEAN DIRECTION OF PRIM. SWELL   N NYY MEAN PERIOD OF ORIMARY SWELL
C
C
C           **** TOGGLE FDB OPTION        ****
C
C     F   CHARACTER*1   ACCESS TO FIELDS DATA BASE..
C          Y  ACCESS IS ENABLED, OTHERWISE ACCESS IS DISABLED
C   F
C8  -
    N
C
C
C           **** TOGGLE PACKING OF 2D SPECTRA GLOBALY ***
C
C     F   CHARACTER*1   PACKING OF 2D SPECTRA..
C          Y  PACKING IS ENABLED, OTHERWISE NO PACKING
C   F
C9  -
    Y
C
C           **** MODEL OPTIONS  ****
C
C     ICASE    INTEGER*7   PROPAGATION OPTION..
C        1  SPHERICAL,  OTHERWISE  CARTESIAN PROPAGATION
C     ISHALLO  INTEGER*7   SHALLOW WATER FLAG..
C        1  DEEP,  OTHERWISE SHALLOW WATER MODEL
C     IREFRA   INTEGER*7  REFRACTION OPTION..
C        0 IF REFRACTION IS NOT USED
C        1 IF DEPTH REFRACTION IS USED
C        2 IF DEPTH AND CURRENT REFRACTION IS USED
C     ITEST    INTEGER*7  TEST OUTPUT LEVEL..
C        > 0 FOR OUPUT UPTO SUBROUTINE LEVEL
C        0  NO TEST OUTPUT
C     ITESTB   INTEGER*7  TEST OUTPUT LEVEL BLOCKS..
C        > 0 FOR BLOCK LEVEL TEST OUTPUT IF ITEST > 0
C     IREST    INTEGER*7  RESTART OPTION..
C        1 RESTART FILES ARE SAVED,  OTHERWISE FILES ARE NOT SAVED.
C     IASSI    INTEGER*7  DATA ASSIMILATION OPTION..
C        1 ASSIMILATION IS DONE,  OTHERWISE IS ARE NOT DONESAVED.
C     IFORCA   INTEGER*7  FORECAST OPTION..
C        1 MODEL STARTS FORM FORECAST FIELDS,  OTHERWISE FROM ANALYSIS.
C
C  ICASE  ISHALLO   IREFRA    ITEST   ITESTB    IREST    IASSI   IFORCA
C10-----  -------  -------  -------  -------  -------  -------  -------
       1        1        0        1        1        1        0        0
C                                                       !!!!!!  !!!!!!!
C
C           **** BOUNDARY POINTS *****
C
C           IBOUNC = COURSE GRID OPTION
C           IBOUNF = FINE GRID OPTION
C
C           IBOUNC MUST BE: 1 FOR COURSE GRID OUTPUT
C                           OTHERWISE NO BOUNDARY POINTS
C           IBOUNF MUST BE: 1 FOR FINE GRID INPUT
C                           OTHERWISE NO BOUNDARY POINTS
C
C IBOUNC   IBOUNF
C11-----  -------
       0        0
C
C           **** FILE NAMES ****
C
C           UID     = USERID              CHARACTER*3
C           RID     = RUN IDENTIFIER      CHARACTER*3
C           PATH    = PATH NAME OF FILES  CHARACTER*60
C
CUSD  RID  PATH
C12-  ---  ---------+---------+---------+---------+---------+---------+
 $USEID  GLO  $EXPER
C
C ----------------------------------------------------------------------
C                END OF USER INPUT PROG CHIEF.
C ----------------------------------------------------------------------
C
EOD

  #=====================================================================
  # PREPARE ICE BOUNDARIES.
  #=====================================================================

cp $YMPPATH/kelv$ymd  ICE
if [ $? != 0 ] ; then
  ecfile -p $EXPPATH/kelv$ymd -enG get ICE
  #
  #  If mars file not there prepare mars request
  #
  if [ $? != 0 ] ; then

    #===================================================================
    # GET SURFACE TEMPERATURES FROM MARS
    #===================================================================

    cat >myrequest <<EOF
    RETRIEVE,LEVTYPE=SFC,TYPE=AN,
      CLASS=OD,
      EXPVER=1,
      LEVEL=OFF,
      PARAM=ST,
      REPRES=GG,
      STEP=0,
      TIME=12,
      DATE=$ymd,
      GRID=3.0/3.0,
      AR= 78./ 0./ -78./357.,
      TARGET="ICE"
      END
EOF

    #===================================================================
    # DO THE MARS REQUEST.
    #===================================================================

    mars myrequest
    if [ $? != 0 ] ; then
      echo "\n\n\t\t  NO SURFACE TEMPS FROM MARS,  NOTHING IS DONE \n\n"
      exit 19
    fi
    ecfile -p $EXPPATH/kelv${ymd} -hO -uD save ICE
    cp ICE $YMPPATH/kelv${ymd}
  fi
fi

  wama.next < wamin > wamout  2> wamerr
  if [ $? != 0 ] ; then
    echo "\n\n\t\t  WAVE MODEL FAILED  \n\n"
    echo "\n\n\t\t  CORE FILE DUMPED   \n\n"
    ls -lsa
    asa wamout
    asa wamerr
    debug -B -s wama.next
    exit 22
  fi
  ls -lsa
  asa wamout
  cat wamerr
# flowview -L -u -l
# flowview -LA
# perfview -LA
# prof -x $YMPPATH/wama.next > my.rawfile
# profview -LA my.rawfile > profview_out
# ecfile -p $EXPPATH/prof.data_next    -hO -uD save prof.data
# ecfile -p $EXPPATH/perf.data_next    -hO -uD save perf.data
# ecfile -p $EXPPATH/flow.data_next    -hO -uD save flow.data
# ecfile -p $EXPPATH/my.rawfile_next    -hO -uD save my.rawfile
# ecfile -p $EXPPATH/profview_out_next    -hO -uD save profview_out
# cp prof.data $YMPPATH/prof.data
# cp perf.data $YMPPATH/perf.data
# cp flow.data $YMPPATH/flow.data

  #=====================================================================
  #   END OF WAM MODEL
  #   SAVE RESTART FILES
  #=====================================================================

  if [ "$HHANAL" != "000" ] ;  then
    ecfile -p /$USEID/cos/$EXPER/lawianal -p $EXPPATH/lawi${ymdm1}a -hO copy
    ecfile -p /$USEID/cos/$EXPER/blspanal -p $EXPPATH/blsp${ymdm1}a -hO copy
    ecfile -p /$USEID/cos/$EXPER/slatanal -p $EXPPATH/slat${ymdm1}a -hO copy
  fi
  if [ "$HHFORC" != "000" ] ;  then
    ecfile -p /$USEID/cos/$EXPER/lawiforc \
            -p $EXPPATH/lawi${ymdm1}f${HHFORC} -hO copy
    ecfile -p /$USEID/cos/$EXPER/blspforc \
            -p $EXPPATH/blsp${ymdm1}f${HHFORC} -hO copy
    ecfile -p /$USEID/cos/$EXPER/slatforc \
            -p $EXPPATH/slat${ymdm1}f${HHFORC} -hO copy
  fi

  #=====================================================================
  #  RESTART FILES ARE SAVED.
  #  END OF JOB WAM.
  #=====================================================================
fi
rm fort.*
ls -lt

#=====================================================================
#  PUSH THE RESULTS FROM FDB INTO MARS IF ASKED FOR.
#=====================================================================

if [ "$ARCHIVE" = ".T." ] ; then
  # ONLY ANALYSIS FIELDS ARE ARCHIVED.

  YYARCH=$YYBEGA
  MMARCH=$MMBEGA
  DDARCH=$DDBEGA
  HHARCH=12

  IDTPRO="${YYARCH}${MMARCH}${DDARCH}${HHARCH}"

  IDTPRO=`newdate "$IDTPRO" +6`
  while [ $IDTPRO -le "19${ymdm1}12" ] ; do
    DATEOFA=`substring $IDTPRO 3 8`
    TIMEOFA=`substring $IDTPRO 9 10`

    cat > marsin << EOFF
      retrieve,
        date=${DATEOFA},
        time=${TIMEOFA},
        step=0,
        expver=11,
        pass="runda",
        param=229/230/231/232/233/234/235/236/237/238/239,
        repres=ll,
        stream=wave,
        levtype=sfc,
        resol=030/030,
        target="fort.63",
        type=analysis
      end
EOFF

    mars -b always marsin || exit 23
    uscpblock -u -f BB fort.63 fort.64

    cat > marsin2 << EOFF
      archive,
        date=${DATEOFA},
        time=${TIMEOFA},
        step=0,
        expver=$EXPVER,
        pass="${USER}xxx",
        param=229/230/231/232/233/234/235/236/237/238/239,
        class=rd,
        repres=ll,
        levtype=sfc,
        resol=030/030,
        source="fort.64",
        disp=old,
        type=analysis
      end
EOFF

    mars < marsin2 || exit 24
    rm fort.6?
    IDTPRO=`newdate $IDTPRO +6`
  done
fi

#=======================================================================
#  PLOT THE RESULTS FROM FDB IF ASKED FOR.
#=======================================================================
if [ "$PLOT" = "WAVE" ] ; then

  #=====================================================================
  #  SECURE FOR ABSOLUT PLTURA.X
  #=====================================================================
  test -r $YMPPATH/pltura.x
  if [ $? != 0 ] ; then
    ecfile -p $EXPPATH/pltura.o  -enG get pltura.o
    if [ $? != 0 ] ; then
      ls -lsa
      echo "\n\n\t\tPLTURA:  BINARY NOT FOUND \n\n"
      echo "\n\n\t\tJOB TERMINATED\n\n"
      exit 25
    fi
    segldr pltura.o -g -V      \
         -l $ECLIB             \
            $EMOSLIB           \
            $MAGLIB            \
         -o pltura.x
    if [ $? != 0 ] ; then
      echo "\n\n\t\t  PLTURA NOT LOADED \n\n"
      exit 26
    fi
    cp pltura.x $YMPPATH/pltura.x
  else
    cp $YMPPATH/pltura.x pltura.x
  fi

  #=====================================================================
  #  SET ENVIRONMENT VARIABLES.
  #=====================================================================

  YEAR=`substring $YYBEGA 3 4`
  MONTH=$MMBEGA
  DAY=$DDBEGA
  HOUR=18
  MIN=00
  ANZAHL=4

  REG=/$EXPER
  PAT=$TEMP$REG
  UN=H
  LENGTH=000006
  for TYPE in FG AN ; do
  START="$YEAR$MONTH$DAY$HOUR$MIN"
  PRODUCT=0
  while [ $PRODUCT -lt "$ANZAHL" ] ; do
    PRODUCT=`expr $PRODUCT + 1`
    DATE=$START
    KONTIME=`newdate 19$DATE +00`
    KONTIME=`substring $KONTIME 3 10`00
    MAT1=`substring $KONTIME 1 6`
    HOU1=`substring $KONTIME 7 8`

# ===================================================================
#    DO MARS REQUEST FOR WAVE HEIGHT.
# ===================================================================

    cat > mcw << eoi
      RET,
       DATE=${MAT1},
       TIME=${HOU1}00,
       PARAM=229,
       DOMAIN=G,
       REPRES=LL,
       CLASS=OD,
       STREAM=WV,
       EXPVER=11,
       TYPE=$TYPE,
       LEVTYPE=SFC,
       LEVELIST=OFF,
       FORMAT=PACKED,
       TARGET="fort.71"
     END
eoi

    /usr/local/bin/mars mcw
    if [ $? != 0 ] ; then
      echo "\n\n\t\tMARS FAILED NOTHING IS DONE\n\n"
      exit 27
    fi

# ===================================================================
#    DO MARS REQUEST FOR WAVE DIRECTIONS
# ===================================================================

    cat > mcw << eoi
      RET,
       DATE=${MAT1},
       TIME=${HOU1}00,
       PARAM=230,
       DOMAIN=G,
       REPRES=LL,
       CLASS=OD,
       STREAM=WV,
       EXPVER=11,
       TYPE=$TYPE,
       LEVTYPE=SFC,
       LEVELIST=OFF,
       FORMAT=PACKED,
       TARGET="fort.72"
     END
eoi

    /usr/local/bin/mars mcw
    if [ $? != 0 ] ; then
      echo "\n\n\t\tMARS FAILED NOTHING IS DONE\n\n"
      exit 28
    fi

    pltura.x << ***FINI*** > OUT
C
C        INPUT OF PROGRAM PLTURA.
C        ========================
C
C        USER INFORMATION IS BEING READ WITH THE PRESUMPTIONS THAT:
C         1. EVERY LINE STARTING WITH 'C ' OR 'C-' IS A COMMENT LINE
C         2. LINE 1 IS A COMMENT LINE
C         3. VALUES ARE PUT IN BELOW POSITIONS INDICATED WITH '-'
C            (RIGHT-JUSTIFIED, BUT CHARACTER LEFT JUSTIFIED)
C ===
C  1.       **** KIND OF PLOTTER TO BE USED ****
C ===       **** INSERT KEYWORD 'DEVICE' ****
C
C           CODES FOR KDEVICE
C           1       = COLOUR PLOT
C           0       = BLACK & WHITE PLOT
C
C           CODES FOR SIZE
C           A4      = A4 PAPER WILL BE USED
C           A3      = A3 PAPER WILL BE USED
C
C key  kdevice  size
C----- -        --
DEVICE $DEVICE        $PAPER
C
C  ===      **** TYPE OF PARAMETER TO BE USED ****
C   2.      **** INSERT KEYWORD 'TYPE' ****
C  ===
C        IT IS TO SPECIFY WHICH PARAMETER SHOULD BE PLOTTED
C        TAKE ONE OF TWO POSSIBILITIES:
C        FIRST:   WAVE HIGHTS  (INSERT WAVE)
C        SECOND:  WIND SPEEDS  (INSERT WIND)
C        FOR THE PARAMETER THERE MUST BE A LINE BEGINNING WITH THE
C        KEYWORD 'TYPE' FOLLOWED BY EITHER 'WAVE' OR 'WIND' AT THE POSITIONS
C        GIVEN AT THE FOLLOWING.
C
C key     type of plot
C---      ----
TYPE      $PLOT
C
C ===       **** THE PERIOD OF INTEREST ****
C  3.       **** INSERT KEYWORD 'TIME' ****
C ===
C           THE PERIOD OF INTERST IS FIXED BY THE TIME OF
C           THE FIRST DATA FIELD USED AND THE LENGTH OF THE
C           PERIOD IN HOURS OR DAYS.  ADDITIONALLY THE TIME OF THE
C           BACKGROUND CONTOUR FIELD WHICH SHOULD APPEAR ON THE PLOT
C           HAS TO BE INSERTED AD CODE kontime.
C
C           CODES FOR start and kontime:
C             ANY 10 DIGIT NUMBER REPRESENTING A VALID DATE/TIME GROUP
C             FORMAT: YYMMDDhhmm
C
C           CODES FOR length:
C             NUMBER OF UP TO 6 DIGITS
C
C           CODES FOR unit:
C             D: THE LENGTH IS COUNTED IN DAYS
C             H: THE LENGTH IS COUNTED IN HOURS
C
C key       start  length  unit      kontime
C---   ----------  ------     -   ----------
TIME   $START  $LENGTH     $UN   $KONTIME
C
C
C ===       **** POSITION OF TIME CORRESPONDING TO A SWATH ****
C  4.       **** INSERT KEYWORD 'PWHERE' ****
C ===
C           THE DATE AND TIME OF A SWATH WILL BE PLOTTED
C           TOGETHER WITH THE SWATH. THE LATITUDE POSITION
C           OF THAT TIME ON THE MAP CAN BY CHOSEN HEREWITH.
C
C           CODES FOR LATITUDE
C             ANY f7.2 REAL NUMBER VALID FOR A LATITUDE   OR
C             -100.00 : THE TIME IS PLOTTED BELOW THE MOST SOUTHERN
C                       COLLOCATION OF THE SWATH.
C
C key    latitude
C-----    ----.--
PWHERE    -100.00
C
C ===       **** PLOT REQUESTS ****
C  5.       **** INSERT KEYWORD "PLOT**
C ===
C           CODES FOR map:
C             T: TO PLOT A MAP
C             F: NO MAP
C           CODES FOR blo:
C             T: TO PLOT BLOCK DIAGRAMS
C             F: NO BLOCK DIAGRAMS
C           CODES FOR sha:
C             T: SHADING OF CONTOUR FIELD
C             F: NO SHADING
C           CODES FOR cus:
C             T: TO PLOT A CUSTER DIAGRAM
C             F: NO CUSTER DIAGRAM
C
C key  map blo sha cus
C---     -   -   -   -
PLOT     T   F   T   T
C
C
C ===       *      ONLY IF MAP OF PLOT IS "T" !! *
C  6.       **** LAYOUT OF THE MAP ****
C ===       **** ADAPT THE MAGICS PARAMETERS AS APROPRIATE ****
C
C                            PSETC
C                            -------------------
SUBPAGE_MAP_PROJECTION       CYLINDRICAL
SUBPAGE_MAP_AREA_DEFINITION  CENTRE
C
C                              PSETR
C                              -----.---
SUBPAGE_MAP_CENTRE_LONGITUDE      15.0
SUBPAGE_MAP_CENTRE_LATITUDE       38.0
SUBPAGE_MAP_VERTICAL_LONGITUDE    15.0
SUBPAGE_MAP_SCALE                 18.0E6
SUBPAGE_LOWER_LEFT_LONGITUDE     120.0
SUBPAGE_LOWER_LEFT_LATITUDE      -90.0
SUBPAGE_UPPER_RIGHT_LONGITUDE    120.0
SUBPAGE_UPPER_RIGHT_LATITUDE      90.0
MAP_GRID_LONGITUDE_INCREMENT      20.0
MAP_GRID_LATITUDE_INCREMENT       10.0
MAP_GRID_LONGITUDE_REFERENCE     120.0
MAP_GRID_LATITUDE_REFERENCE      -90.0
WIND_THINNING_FACTOR               1.0          2.0
WIND_ARROW_UNIT_VELOCITY          10.0         10.0
NUMBER_HIGHT                       0.35        0.35
C                              PSETI
C                              -----
MAP_LABEL_LONGITUDE_FREQUENCY     20
MAP_LABEL_LATITUDE_FREQUENCY      10
NUMBER_THINNING_FACTOR             1            2
NUMBER_OF_CONTOUR_LEVELS           9
C
C                  REAL REAL REAL REAL REAL REAL REAL REAL REAL REAL REAL
C                  ---. ---. ---. ---. ---. ---. ---. ---. ---. ---. ---.
CONTOUR_LEVEL_LIST   0.   1.0  2.   3.0  4.   5.   6.   8.  25.    .    .
CONTOUR_DENSITIES  100. 300. 600. 300. 600. 300. 600. 600.    .    .
C
C                  INT  INT  INT  INT  INT  INT  INT  INT  INT  INT  INT
C                  ---. ---. ---. ---. ---. ---. ---. ---. ---. ---. ---.
CONTOUR_COLOURS      4    4    4    5    5    3    3    6
NUMBER_COLOURS      56   56   56    5    5    3    3    6
C
C
C ===       **** PRINTER OUTPUT REQUEST ****
C  7.
C ===
C
C   THERE WILL BE STANDARD OUTPUT OF EACH SUB-ROUTINE DOWN
C   TO THE LEVEL "lsuvi"     **** INSERT KEYWORD 'LSUVI' ****
C
C   THERE WILL BE TEST OUTPUT OF EACH SUB-ROUTINE DOWN
C   TO THE LEVEL "ltest"     **** INSERT KEYWORD 'LTEST' ****
C
C   THERE WILL BE STANDARD OUTPUT OF  EACH
C   SUB-ROUTINE MENTIONED IN THE LIST "csuvi"
C                            **** INSERT KEYWORD 'CSUVI' ****
C   THERE WILL BE TEST OUTPUT OF  EACH
C   SUB-ROUTINE MENTIONED IN THE LIST "ctest"
C                            **** INSERT KEYWORD 'CTEST' ****
C
C key  level
C----     --
LSUVI      1
LTEST      0
C
C key  sub-routine
C----     --------
CSUVI
CTEST
C
C ===       **** LIST OF USER INPUT FINISHED ****
C  8.       **** INSERT KEYWORD 'END' ****
C ===
C key
C--
END
***FINI***

    if [ $? = 0 ] ; then
      mv ps $YMPPATH/$PLOT$PRODUCT.ps
      if [ $? = 0 ] ; then
        echo 'hi'
        gksplot $COL  -f ps -d ps -n -k -p $PAPER \
              -j $PLOT -u $USER $YMPPATH/${PLOT}${PRODUCT}.ps
      fi
    else
      debug -B -s pltura.x | cat
      PRODUCT=200
    fi
    cat OUT
    asa LOG
    rm fort.7? fort.1
    START=`newdate 19$START +6`
    START=`substring $START 3 10`00
  done
  ls -l
  done
fi
env
exit 29
