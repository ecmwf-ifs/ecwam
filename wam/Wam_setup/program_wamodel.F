#ifdef JOB
#.
#. THIS JOB EXECUTES THE SUBROUTINE VERSION OF THE WAMODEL,
#. DRIVEN BY WINDS REQUESTED FROM MARS.
#.
#QSUB -r wamodel
#QSUB -q vpp700.normal
#QSUB -lP 1
# # #QSUB -a "tomorrow 4:00am"
#QSUB -lt 2500
#QSUB -lT 3000
#QSUB -lM 500mb
# # #QSUB -lM 1200mb
#QSUB -eo
#QSUB -mb -mi
#QSUB
#.
######################################################################
#.
PATH=$PATH:.
trap '{ echo "Signal SIGHUP received and trapped" ; }' 1
trap '{ echo "Signal SIGINT received and trapped" ; }' 2
set +v
#.
#.   Specify whether grib spectra are provided at input and output
#.
lgribin=T
lgribout=T

#.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#.   Specify whether input grib CD file is NOT provided
#.   It should only be the case for cold start !!!!!!!!!!!!!!!
#.
lnocdin=T
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#.
#.=====================================================================
#.  ALL PATHES, TIMES AND THE PRECISION CAN BE DEFINED IN path.h
#.  COMMONLY FOR ALL JOBS.
#.=====================================================================
#.
#include "path.h"
#.
# move to working directiory
############################
cd $WDIR
mkdir $WORKDIR
cd $WORKDIR
#.
#.=====================================================================
#. General grid information and model constants.
#.=====================================================================
#.
if [[ ! -r $ADIR/gridglou_${grid} ]] ; then
  ls -lsa
  echo "\n\n\t\tPREPROC file gridglou was not found\n\n"
  echo "\n\n\t\tPROGRAM IS TERMINATED\n\n"
  exit 1
else
  ln -s $ADIR/gridglou_${grid} fort.7
fi
#if region=='s'
#.
#. THE SWAMP CASE DOES NOT NEED EXTERNAL WINDS.
#.
#else
#.
#.
#if iceedge=='y'
#.check if SST file does exist
#.
if [[ ! -r $ADIR/kelv${begofrn}_${grid} ]] ; then
  print - "\n\n\t\t NO SURFACE TEMPS , NOTHING IS DONE \n\n"
  exit 2
else
  ln -s $ADIR/kelv${begofrn}_${grid} sfctempin
fi
#endif
#.
# ALTIMETER DATA IF NEEDED
############################
if [[ $ASSIMILATION = "YES" ]] ; then
  AS=1
  #
  ANALYSIS_DATE=$(substring $begofrn 1 10)
  ANALYSIS_DATE=$(newdate ${ANALYSIS_DATE} 6)
  LAST_ANALYSIS_DATE=`substring $endofrn 1 10`
  #
  until [[ $ANALYSIS_DATE -gt $LAST_ANALYSIS_DATE ]] ; do
    if [[ !  -r ${ADIR}/RFL${ANALYSIS_DATE}00 ]] ; then
      print - "\n\n\t\t SATELLITE DATA NOT ONLINE  \n\n"
      print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
      exit 3
    else
      ln -s ${ADIR}/RFL${ANALYSIS_DATE}00 .
    fi
    ANALYSIS_DATE=`newdate ${ANALYSIS_DATE} 6`
  done
else
  AS=0
fi
#.
# SURFACE WINDS
################
wind_filename=sfcwind_${grid}_${begofrn}_${begoffo}_${endofrn}
if [[ ! -r $WDIR/$wind_filename ]] ; then
  print - "\n\n\t\tWIND DATA NOT ONLINE\n\n"
  print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
  exit 4
else
  ln -s $WDIR/$wind_filename sfcwindin
  print - "\n\n\t\tSTART WAVE MODEL \n\n"
fi
#endif
#.
#.
#.=====================================================================
#.  WAVE MODEL execution.
#.=====================================================================
#.
#. Common UBUF.
#.
if [[ ! -r $ADIR/ubufglou_${grid} ]] ; then
  print - "\n\n\t\tPREPROC file ubufglou was not found\n\n"
  print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
  exit 5
else
  ln -s $ADIR/ubufglou_${grid} fort.8 
fi
#.
#. Coarse grid output points, if this is a coarse grid run.
#.
#. assign -a $ADIR/bouccatu_${grid} fort.9
#.
#. Fine grid input points, if this is a fine grid run.
#.
#. assign -a $ADIR/boufcatu_${grid} fort.10
#.
#.======================================================================
#. GET RESTART FILES
#.======================================================================
#. !!!! at this time we assume that the restart files have been saved
#. on disk in directory $WDIR and $ADIR
#.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#.
# Get last wind.
################
#.
if [[ $lgribin != "T" ]] ; then
  if [[ -r $WDIR/LAW${begofrn}00 ]] ; then
    ln -s $WDIR/LAW${begofrn}00 LAW${begofrn}00 
  else
    ls -lsa
    print - "LAW"${begofrn}"00 file was not found"
    print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
    exit 5
  fi
else
  if [[  -r $ADIR/CD${begofrn}_${grid} ]] ; then
    ln -s $ADIR/CD${begofrn}_${grid} cdwavein
  elif [[ $lnocdin = "T" ]] ; then
    print - "##################################################"
    print - " "
    print - "the program should start wihtout an input CD file"
    print - "as none is provided !!!!!"
    print - " "
    print - "##################################################"
  else
    ls -lsa
    print - "CD"${begofrn}_${grid}" file was not found"
    print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
    exit 5
  fi
  if [[  -r $ADIR/U10${begofrn}_${grid} ]] ; then
    ln -s $ADIR/U10${begofrn}_${grid} uwavein
  else
    print - "U10"${begofrn}" file was not found"
    print - "\n\n\t\tPROGRAM will continue as this can happen\n\n"
  fi
fi
#.
# Get initial spectrum.
#######################
#.
if [[ $lgribin != "T" ]] ; then
  if [[ -r $WDIR/BLS${begofrn}00 ]] ; then
    ln -s $WDIR/BLS${begofrn}00 BLS${begofrn}00 
  else
    ls -ls
    print - "BLS"${begofrn}"00 file was not found"
    print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
    exit 7
  fi
else
  if [[ -r $WDIR/SGS${begofrn}00 ]] ; then
    ln -s $WDIR/SGS${begofrn}00 specwavein 
  else
    ls -ls
    print - "SGS"${begofrn}"00 file was not found"
    print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
    exit 8
  fi
fi
#.
#.
#.
#.
#.======================================================================
#. load the WAM MODEL.
#.======================================================================
#.
if [[ ! -r $LIBS/$WAMLIB.a ]] ; then
  print - "\n\n\t\tWAMLIB file was not found\n\n"
  print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
  exit 12
fi
#.
#. LOAD LIBRARY
#.
ar vx $aux ifssig.o
ar vx $dum ishell.o
ar xv $LIBS/$WAMLIB.a chief.o
#.
frtcnt=$(cat frtcnt)
if [ $frtcnt = 0 -a -f $BIN/chief -a "$REMAKE" != YES  ] ; then
  print - WAM Model ready
else
  frt $rp -Wl,-m,-P,-uMAIN__ -X9 \
    -o $BIN/chief  chief.o ishell.o ifssig.o \
     -lwam  \
     $ECLIB $EMOSLIB $FDBLIB  \
     -Wl,-dy -lsocket $MPELIB  $MPLIB || {
    #.
    #. Check return code.
    #.
     print - "\n\n\t\tLOAD FAILED\n\n"
     print - "\n\n\t\tPROGRAM WILL TERMINATE\n\n"
     exit 13
    }
  #.
fi
#.
ls -ltr
###if [[ $NPES -ne 4 && $NPES -ne 8 && $NPES -ne 16 && $NPES -ne 24 ]] ; then
###  print - "\n\n\t\tWAVE MODEL IS COMPILED. NOW SET THE -lP TO 4, 8, 16 or 24 AND SUBMIT AGAIN\n\n"
###  cd $WDIR
###  \rm -rf $WORKDIR
###  exit
###else
  print - "\n\n\t\tWAVE MODEL WILL START\n\n"
###fi
#.
cd ${WDIR}/${WORKDIR}
#.
#. User input file
#.
cat > fort.3 << EOD
#if region=='s'
#include "input_wamodel_rswamp2.F"
#elif region=='m'
#if resolution == 50
#include "input_wamodel_rmedite.F"
#elif resolution == 25
#include "input_wamodel_rmedite_25.F"
#endif
#else
#if resolution == 150
#include "input_wamodel_rglobal150.F"
#elif resolution == 50
#include "input_wamodel_rglobal50.F"
#else
#include "input_wamodel_rglobal.F"
#endif
#endif
EOD
#.
#. RUN MODEL
#.
#. Increase size of the mail box.
#.
# specify mail box size:
VPP_MBX_SIZE=16000000 ; export VPP_MBX_SIZE
#.
#. Need to define those 2 variables for profiling : 
#. and link with libmp_hybrid.a
#####VPP_STATS=8 ; export VPP_STATS 
#####VPP_DETAIL=-1 ; export VPP_DETAIL
#.
if [[ $lgribin != "T" ]] ; then
  #modify runtime buffers
  #######################
  fu12bf=4096 ; export fu12bf
  fu15bf=344 ; export fu15bf
  fu63bf=524 ; export fu63bf
  fu90bf=53 ; export fu90bf
fi
#.
$BIN/chief 1> stdout 2>stderr || {
  #.
  #.  Abort on wave model failure.
  #.
  #.
  #.  Error: Copy program output to stdout and terminate.
  #.
   echo $?
   print - "\n\n\t\tWAVE MODEL FAILED\n\n"
   print - "\n\n\t\tOUTPUTS CAN BE FOUND ON\n\n"
   print - $WDIR"/"$WORKDIR
   ls -l
   ls -l  logfile.*
   exit 14
  }
ls -lsa
#.
#. Dump logfiles.
#################
#.
ls -l  logfile.*
#.
# stdout
########
cat stdout
# stderr
########
cat stderr
# logfile.1
############
cat logfile.1
#.
cd $WDIR
\rm -rf $WORKDIR
#.
#.
#.======================================================================
#.END OF JOB.
#.======================================================================
#.
#endif
