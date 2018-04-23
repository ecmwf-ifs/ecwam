#ifdef JOB
# @ shell = /usr/bin/ksh
# @ core_limit = 4096
# @ notification = error
# @ restart = no
# @ cpu_limit = 4000,2000
# @ job_cpu_limit = 7600,5600
# @ class = any 
# @ job_name = wam_wamodel
# @ output = /scratch/rd/wab/wam_wamodel.$(jobid).out
# @ error  = /scratch/rd/wab/wam_wamodel.$(jobid).err
# @ account_no = ecrmwa
# @ job_type = parallel
# @ total_tasks = 2
# @ queue

#.
#. THIS SCRIPT EXECUTES THE WAMODEL,
#.
#. All necessary input initial data and forcing wind fields
#. must have been retrieved before hand with the preset task
#. All grid informations and tables must also have been created with
#. task preproc.
#. 
#. WAM can in principle run in parallel with any number of processors
#. The type of model decomposition is controlled by LL1D in the input
#. namelist.
#.
######################################################################
#.
PATH=$PATH:.
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
#.   parameters will end up in fort.20 and the restart files
#.   (stress and spectra) in
#.   LAWYYYYMMDDHH00 and BLSYYYYMMDDHH00 in directory specified by
#.   $STORAGE_PATH. (=$WDIR)
#.   If lfdb=T you also need to specify FDB_ROOT
#.   OUTSIDE ECMWF you should have lfdb=F
#.
lfdb=F
#.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#.   Specify whether input grib CD file is NOT provided
#.   It should only be the case for cold start (i.e. when preset has created
#.   the spectra). Normally grib CD is provided and lnocdin=F
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
mkdir tmp$$
cd tmp$$
#.
#.=====================================================================
#. General grid information and model constants.
#.=====================================================================
#.
if [[ ! -r $WDIR/gridglou_${grid} ]] ; then
  ls -lsa
  echo "\n\n\t\tPREPROC file gridglou was not found\n\n"
  echo "\n\n\t\tPROGRAM IS TERMINATED\n\n"
  exit 1
else
  ln -s $WDIR/gridglou_${grid} fort.7
fi
#.
if [[ ! -r $WDIR/ubufglou_${grid} ]] ; then
  print - "\n\n\t\tPREPROC file ubufglou was not found\n\n"
  print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
  exit 5
else
  ln -s $WDIR/ubufglou_${grid} fort.8 
fi
#.
#.======================================================================
#. GET RESTART FILES
#.======================================================================
#. !!!! at this time we assume that the restart files have been saved
#. on disk in directory $WDIR
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
  if [[  -r $WDIR/CD${begofrn}_${grid} ]] ; then
    ln -s $WDIR/CD${begofrn}_${grid} cdwavein
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
  if [[  -r $WDIR/U10${begofrn}_${grid} ]] ; then
    ln -s $WDIR/U10${begofrn}_${grid} uwavein
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
# SST or sea ice fraction for sea ice mask 
#.=====================================================================
#.check if SST file does exist
#.
if [[ ! -r $WDIR/SSTDATA${begofrn}_${grid} ]] ; then
  print - "\n\n\t\t NO SURFACE TEMPS , WAM CAN RUN WITHOUT IT\n\n"
else
  ln -s $WDIR/SSTDATA${begofrn}_${grid} sfctempin
fi
#endif
#.
#.
#. User input file (SEE SUBROUTINE USERIN FOR A DESCRIPTION OF THE NAMELISTS)
#############################################################################
#.
if [[ $lgribout != "T" ]] ; then
  lfdbioout=F
else
  lfdbioout=$lfdb
fi
#.
# WAM model namelists (see subroutine userin):
##############################################
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
############
#.
#. Increase size of the mail box.
#.
# specify mail box size (for message passing):
MBX_SIZE=64000000 ; export MBX_SIZE
export MPL_MBX_SIZE=$MBX_SIZE
#

${BINS}/chief 1> stdout 2>stderr || {
  #.
  #.  Abort on wave model failure.
  #.
  #.
  #.  Error: Copy program output to stdout and terminate.
  #.
   echo $?
   print - "\n\n\t\tWAVE MODEL FAILED\n\n"
   ls -l  logfile.*
   exit 14
  }
ls -lsa
#.
#. Dump logfiles.
#################
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
\rm -rf tmp$$

#.
#.
#.======================================================================
#.END OF JOB.
#.======================================================================
#.
#endif
