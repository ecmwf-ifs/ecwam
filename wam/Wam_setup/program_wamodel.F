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
#QSUB -lM 1200mb
# # #QSUB -lM 500mb
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
#.   Specify whether grib split spectra are provided
#.   (it's the operstional setting
#.
lgribin=T
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
#.  COMONLY FOR ALL JOBS.
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
if [ $COMPILE =  1 ] ; then
  #.
  cd $LIBS
  \rm -rf tmpdir || true
  mkdir tmpdir
  cd tmpdir
  #.
  #.
  #.======================================================================
  #. Copy the WAM MODEL to the compiler input.
  #.======================================================================
  #.
  #.
  cat > wamodel.src <<\EOFPRO
#endif
#include "wave_module.F"
#include "chief.F"
#include "buildstress.F"
#include "closend.F"
#include "airsea.F"
#include "bouinpt.F"
#include "dotdc.F"
#include "femean.F"
#include "fillbl.F"
#include "file_transfer.F"
#include "getspec.F"
#include "getstress.F"
#include "getwnd.F"
#include "gradi.F"
#include "grstname.F"
#include "headbc.F"
#include "igs2grs.F"
#include "implsch.F"
#include "incint.F"
#include "initmdl.F"
#include "intpol.F"
#include "locint.F"
#include "makeblo.F"
#include "makegrid.F"
#include "meansqs.F"
#include "newwind.F"
#
#include "mpbcastgrdfld.F"
#include "mpbcastgrid.F"
#include "mpbcastscfld.F"
#include "mpcrtbl.F"
#include "mpdecomp.F"
#include "mpdistribfl.F"
#include "mpdistribintfld.F"
#include "mpdistribscfld.F"
#include "mpexchng.F"
#include "mpgatherersfile.F"
#include "mpgatherfl.F"
#include "mpgathergrdfld.F"
#include "mpgatheroifld.F"
#include "mpgatherscfld.F"
#include "mpgatherspp.F"
#include "mpuserin.F"
#
#include "mt_allot.F"
#
#include "notim.F"
C #ifdef 
C #include "openfil_unicos.F"
C #endif
#include "outbeta.F"
#include "outbc.F"
#include "outbs.F"
#include "outers.F"
#include "outgrid.F"
#include "outint.F"
#include "outwnorm.F"
#include "outspp.F"
#include "outxt.F"
#include "peakfr.F"
#include "phys.F"
#include "prealloc_file.F"
#include "prewind.F"
#include "propags.F"
#include "propdot.F"
#include "readbou.F"
#include "readpre.F"
#include "readwgrib.F"
#include "readfl.F"
#include "readstress.F"
#include "writefl.F"
#include "writestress.F"

#if region=='s'
#ifndef READWND
#include "readwnd_rswamp2.F"
#define READWND
#endif
#endif

#ifdef SEASAT
#ifndef READWND
#include "readwnd_seasat.F"
#define READWND
#define READWND
#endif
#endif

#ifndef READWND
#include "readwnd_mars.F"
#define READWND
#endif

#include "savstress.F"
#include "savspec.F"
#include "sbottom.F"
#include "sdissip.F"
#include "semean.F"
#include "sepwisw.F"
#include "sinput.F"
#include "snonlin.F"
#include "spec2fdb.F"
#include "splitbl.F"
#include "spr.F"
#include "sthq.F"
#include "stresso.F"
#include "timin.F"
#include "userin.F"
#include "wamodel.F"
#include "wamwnd.F"
#include "wavemdl.F"
#include "wsmfen.F"
#include "wvalloc.F"
#include "wvdealloc.F"
#include "analyse.F"
#include "f4spec.F"
#include "fdur.F"
#include "fustar.F"
#include "fwsea.F"
#include "grdata.F"
#include "grfield.F"
#include "gribpac.F"
#include "oifield.F"
#include "readice.F"
#include "micep.F"
#include "readt.F"
#include "findb.F"
#include "confile.F"

#ifdef GEOASAT
#include "readsat_geosat.F"
#define READSAT
#endif

#ifdef SEASAT
#ifndef READSAT
#include "readsat_seasat.F"
#define READSAT
#endif
#endif

#ifdef ERS
#ifndef READSAT
#include "readsat.F"
#define READSAT
#endif
#endif

#include "setmat.F"
#include "setice.F"
#include "update.F"
#include "upwspec.F"
#include "wamassi.F"
#include "wamnorm.F"
#ifdef JOB
EOFPRO
  #.
  #.======================================================================
  #. Compile and load the WAM MODEL.
  #.======================================================================
  #.
  #. COMPILE SOURCE CODE
  #.
  /home/rd/sth/bin/f90split < wamodel.src
  #.
  [ ! -d tmp ] && mkdir tmp
  #.
  [ -f ${LIBS}/libwam.f.tar ] && {
    dd if=${LIBS}/libwam.f.tar of=libwam.f.tar bs=400k
    cd tmp
    tar xf ../libwam.f.tar
    cd ..
  }
  #.
  [ -f ${LIBS}/libwam.a ] && {
    dd if=${LIBS}/libwam.a of=libwam.a bs=400k
    ar xv libwam.a wave_module
  }
  #.
  integer frtcnt
  frtcnt=0
  export frtcnt
  print $frtcnt > frtcnt
  #.
  (
    if [ "$REMAKE" != YES ] ; then
      if [ -f tmp/wave_module.fm ] ; then
        if cmp wave_module.f tmp/wave_module.fm ] ; then
          if ar tv libwam.a | grep wave_module ; then
            exit 9
          fi
        fi
      fi
    fi
    frt $rp -Of -M. -Wv,-Of -Pos -Z wave_module.l -X9 \
        -c wave_module.f
    (( frtcnt += 1))
    print $frtcnt > frtcnt
  )
  mv wave_module.f wave_module.fm
  #.
  i=0
  atst=0
  #.
  for file in $(ls -1 *.f) ; do
    export file
    (
    if [ "$REMAKE" != YES ] ; then
      if [ -f tmp/$file ] ; then
        if cmp $file tmp/$file ; then
          if ar tv libwam.a | grep ${file%.f}.o ; then
            exit 10
          fi
        fi
      fi
    fi
    if [   $file = spec2fdb.f  -a $rp = "-Ad" ] ; then
      rp="-CcdRR8"
    elif [   $file = getspec.f  -a $rp = "-Ad" ] ; then
      rp="-CcdRR8"
    else
      rp=$rp
    fi
    frt $rp -Of -M. -Wv,-Of -Pos -Z ${file%.f}.l -X9 \
        -c -o ${file%.f}.o ${file} || print - ${file%.f}.l >> frt_error
    (( frtcnt += 1))
    print $frtcnt > frtcnt
    )&
    (( i += 1 ))
    [[ $i -gt $atst ]] && {
      atst=0
      i=0
      wait
    }
  done
  wait
  print - "compile done"
  #.
  #. BUILD LIBRARY
  #.
  if [[ $(ls -1 *.o | wc -l ) -gt 0 ]] ; then
  ar rv libwam.a *.o wave_module
  fi
  #.
  tar cf libwam.l.tar *.l
  tar cf libwam.f.tar *.f *.fm
  #.
  #. SAVE LIBRARY
  #.
  ############ecp -o libwam.* ec:$DHSPATH/.
  #.
  ls -la
  #.
  cp libwam.a $LIBS/libwam.a 
  cp libwam.f.tar $LIBS/libwam.f.tar
  cp libwam.l.tar $LIBS/libwam.l.tar
  #.
  chmod go+rx $LIBS/libwam.*
  #.
  if [ -r frt_error ] ; then
    cat $(cat frt_error)
    exit 11
  fi
  #.
  if [[ ! -r $LIBCR/lib$WAMCRLIB.a ]] ; then
    print - "\n\n\t\tWAMCRLIB file was not found\n\n"
    print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
    exit 12
  fi
  #.
  #. LOAD LIBRARY
  #.
  ar vx $aux ifssig.o
  ar vx $dum ishell.o
  [ ! -a chief.o ] &&  ar xv libwam.a chief.o
  #.
  frtcnt=$(cat frtcnt)
  if [ $frtcnt = 0 -a -f $BIN/chief -a "$REMAKE" != YES  ] ; then
    print - WAM Model ready
  else
    frt $rp -Wl,-m,-P,-uMAIN__ -X9 \
      -o $BIN/chief  chief.o ishell.o ifssig.o \
      libwam.a -lwamcr libwam.a \
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
  cd $LIBS 
  \rm -rf tmpdir
fi
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
