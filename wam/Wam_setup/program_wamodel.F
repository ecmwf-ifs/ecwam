#ifdef JOB
#.
#. THIS SCRIPT EXECUTES THE WAMODEL,
#.
#. All necessary input initial data and forcing wind fields
#. must have been retrieved before hand with the preset task
#. All grid informations and tables must also have been created with
#. task preproc.
#. 
#.
#QSUB -r wamodel
#QSUB -q vpp700.normal
#QSUB -lPv 8
#QSUB -lt 1000
#QSUB -lT 1000
#QSUB -lM 700mb
#QSUB -eo
#QSUB -mb -mi
#QSUB -lV 10Mb
#QSUB -ls 100Mb
#.
######################################################################
#.
PATH=$PATH:.
trap '{ echo "Signal SIGHUP received and trapped" ; }' 1
trap '{ echo "Signal SIGINT received and trapped" ; }' 2
set +v
#.
#.   Specify whether grib spectra are provided at input (see preset)
#.
lgribin=T
#.
#.   Specify whether spectra are output in grib format (or binary) 
#.
lgribout=T
#.
#.   Specify whether output is directed to the FDB or to files.
#.   when lfdb=F then grib output of integrated parameters will end up in
#.   fort.30 and the grib spectra in SGSYYYYMMDDHH00 in directory specified
#.   by $STORAGE_PATH. In the case of binary file outputs, the integrated
#.   will end up in fort.20 and the restart files (stress and spectra) in
#.   LAWYYYYMMDDHH00 and BLSYYYYMMDDHH00 in directory specified by
#.   $STORAGE_PATH.
#.   If lfdb=T you also need to specify FDB_ROOT (see below)
#.
lfdb=F
#.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#.   Specify whether input grib CD file is NOT provided
#.   It should only be the case for cold start (i.e. when preset has created
#.   the spectra) 
#.
lnocdin=F
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#.
#include "wkstation_vpp_share.h"
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
#.
if [[ ! -r $ADIR/ubufglou_${grid} ]] ; then
  print - "\n\n\t\tPREPROC file ubufglou was not found\n\n"
  print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
  exit 5
else
  ln -s $ADIR/ubufglou_${grid} fort.8 
fi
#.
#.======================================================================
#. GET RESTART FILES
#.======================================================================
#. !!!! at this time we assume that the restart files have been saved
#. on disk in directory $WDIR and $ADIR
#.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#.
# Get last stress.
##################
#.
if [[ $lgribin != "T" ]] ; then
  #. if binary files are used as input
  if [[ -r $WDIR/LAW${begofrn}00 ]] ; then
    ln -s $WDIR/LAW${begofrn}00 LAW${begofrn}00 
    ln -s $WDIR/LAW${begofrn}00 fort.15
  else
    ls -lsa
    print - "LAW"${begofrn}"00 file was not found"
    print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
    exit 5
  fi
else
  #. if grib files are used as input
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
  elif [[ $lnocdin = "T" ]] ; then
    print - "##################################################"
    print - " "
    print - "the program should start wihtout an input CD file"
    print - "as none is provided !!!!!"
    print - " "
    print - "##################################################"
  else
    ls -lsa
    print - "U10"${begofrn}_${grid}" file was not found"
    print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
    exit 5
  fi
fi
#.
# Get initial spectrum.
#######################
#.
if [[ $lgribin != "T" ]] ; then
  #. if binary files are used as input
  if [[ -r $WDIR/BLS${begofrn}00 ]] ; then
    ln -s $WDIR/BLS${begofrn}00 BLS${begofrn}00 
    ln -s $WDIR/BLS${begofrn}00 fort.12
  else
    ls -ls
    print - "BLS"${begofrn}"00 file was not found"
    print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
    exit 7
  fi
else
  #. if grib files are used as input
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
#.=====================================================================
# SURFACE WINDS
#.=====================================================================
#.
#if region=='s'
#.
#. THE SWAMP CASE DOES NOT NEED EXTERNAL WINDS.
#.
#else
#.
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

#if iceedge=='y'
#.=====================================================================
# SST for sea ice mask 
#.=====================================================================
#.check if SST file does exist
#.
if [[ ! -r $ADIR/SSTDATA${begofrn}_${grid} ]] ; then
  print - "\n\n\t\t NO SURFACE TEMPS , NOTHING IS DONE \n\n"
  exit 2
else
  ln -s $ADIR/SSTDATA${begofrn}_${grid} sfctempin
fi
#endif
#.
#.
#.======================================================================
#. load the WAM MODEL.
#.======================================================================
#.
#. copy library from workstation

cp $WKLIBS/$WAMLIB.a $LIBS/$WAMLIB.a

if [[ ! -r $LIBS/$WAMLIB.a ]] ; then
  print - "\n\n\t\tWAMLIB file was not found\n\n"
  print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
  exit 12
fi
#.
#. LOAD LIBRARY
#.
ar xv $LIBS/$WAMLIB.a chief.o
#.
frtcnt=$(cat frtcnt)
if [ $frtcnt = 0 -a -f $BIN/chief -a "$REMAKE" != YES  ] ; then
  print - WAM Model ready
else
#.
  frt $rp -Wl,-m,-P,-uMAIN__ -X9 \
    -o $BIN/chief  chief.o \
     -lwam  \
     ${mpllib} ${emoslib} ${eclib} ${fdblib} ${naglib} ${lapacklib} \
     -Wl,-dy -lsocket ${blaslib} ${mpilib} ${cvplib} ${mp2lib} || {
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
#.
#. !!! the model is run in WORKDIR !!!
cd ${WDIR}/${WORKDIR}
#.
#. User input file (SEE SUBROUTINE USERIN FOR A DESCRIPTION OF THE NAMELISTS)
#.
if [[ $lgribout != "T" ]] ; then
  lfdbioout=F
else
  lfdbioout=$lfdb
fi
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
# specify mail box size (for message passing):
VPP_MBX_SIZE=64000000 ; export VPP_MBX_SIZE
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
