#ifdef JOB
#.
#. THIS JOB EXECUTES THE SUBROUTINE VERSION OF THE WAMODEL,
#. DRIVEN BY WINDS REQUESTED FROM MARS.
#.
#. LIANA ZAMBRESKY   GKSS/ECMWF    JUNE  1989
#. HEINZ GUNTHER     GKSS/ECMWF    JUNE  1990
#. PETER JANSSEN     KNMI          OCT   1990   (MIGRATION TO UNICOS)
#. HEINZ GUNTHER     GKSS/ECMWF    APRIL 1991   (CYCLE_3.5)
#. B HANSEN          ECMWF         NOVEMBER 1995 MIGRATION TO CLEARCASE
#. J BIDLOT          ECMWF         MAY 1996 MESSAGE PASSING 
#. B HANSEN          ECMWF         APRIL 1997 RESTART FACILITY
#.
#QSUB -r wamodel_CY18R4_test
#QSUB -q vpp700.normal
#QSUB -lP 1
# # #QSUB -a "tomorrow 4:00am"
#QSUB -lt 1800
#QSUB -lT 2000
#QSUB -lM 700mb
#QSUB -eo
#QSUB -mb -me -mi
#QSUB
#.
######################################################################
#.
PATH=$PATH:.
trap '{ echo "Signal SIGHUP received and trapped" ; }' 1
trap '{ echo "Signal SIGINT received and trapped" ; }' 2
cd $TMPDIR
set +v
#.
ASSIMILATION=YES
#.
#.=====================================================================
#.  ALL PATHES, TIMES AND THE PRECISION CAN BE DEFINED IN path.h
#.  COMONLY FOR ALL JOBS.
#.=====================================================================
#.
#include "path.h"
#.
CLASS=RD
NENS=000
TNE=000
typeset -l XID=${USER}a
#.
#.
#.====================================================================
# DEFINE FBD SERVER ON VPP-700
#.=====================================================================
#
#.export FDB_CONFIG_FILE=/vpp700/mrfs/fdb/FDbConfig_vpp700
export FDB_ROOT=/vpp700/fdb18           # /vpp700/fdb22
export FDB_CONFIG_MODE=async            # eps
export FDB_SERVER_HOST=vpp700-x18       # vpp700-x22   vpp700-x08  == wavefdb
export FDB_SIGNALS=no
#.
#.=====================================================================
#. General grid information and model constants.
#.=====================================================================
#.
if [[ ! -r $VPPPATH/gridglou ]] ; then
  ecp ec:$DHSPATH/gridglou $VPPPATH/gridglou || {
    ls -lsa
    echo "\n\n\t\tPREPROC file gridglou was not found\n\n"
    echo "\n\n\t\tPROGRAM IS TERMINATED\n\n"
    exit
  }
fi
ln -s $VPPPATH/gridglou fort.7
#if region=='s'
#.
#. THE SWAMP CASE DOES NOT NEED EXTERNAL WINDS.
#.
#else
#.
#######################################
#.check if need SST for ice boundaries 
#.
#if iceedge=='y'
#if resolution == 150
  GRID="1.5/1.5"
  AREA="81./ 0./ -81./358.5"
#elif region == 'g' &&  resolution == 50
  GRID="0.5/0.5"
  AREA="81./ 0./ -81./359.5"
#elif region == 'm' && resolution == 25
  GRID="0.25/0.25"
  AREA="66./ -6./ 30./42."
#endif
#.
#.check if SST file does exist
#.
if [[ !  -r $VPPPATH/kelv$begofan ]] ; then
  # in case file is not there do a mars request
  # GET SURFACE TEMPERATURES FROM MARS
  #===================================
  ICEDATE=`substring $begofan 1 6`
  cat > myrequest <<-EOF
	RETRIEVE,LEVTYPE=SFC,TYPE=AN,
	CLASS=OD,
	EXPVER=1,
	LEVEL=OFF,
	PARAM=ST,
	REPRES=GG,
	STEP=0,
	TIME=12,
	DATE=$ICEDATE,
	GRID=$GRID,
	AR=$AREA,
	TARGET="fort.90"
	END
	EOF
  #.
  #.DO THE MARS REQUEST.
  #=====================
  mars myrequest || {
    print - "\n\n\t\t NO SURFACE TEMPS FROM MARS, NOTHING IS DONE \n\n"
    exit 19
  }
  cp fort.90 $VPPPATH/kelv$begofan
else
  ln -s $VPPPATH/kelv$begofan fort.90
fi
#endif
#.
if [[ $ASSIMILATION = "YES" ]] ; then
  #.
  #.====================================================================
  #.
  #. Get quality controlled ERS-2 altimeter data if assimilation is
  #. required
  #.
  #. Note the location of the data in CFS is liable to change !!!!
  #.
  #. as well as the way day are kept in the archive (see Bjorn)
  #.
  #.====================================================================
  ##
  AS=1
  CFSALT=/emos/wave/0012/ers2_b/alt/global
  CFSALT=/emos/wave/0001/ers2_s/alt/global
  #
  # get 64 bit decoder from prepifs since the operational one is 32 bit
   if [ $PRECISION = DOUBLE ] ; then
      ALTLIB=/home/rd/rdx/lib/16r4/libalt.a
      ar xv $ALTLIB enderac.o
      frt -Ad -o enderac enderac.o $ALTLIB $EMOSLIB $ECLIB
   fi
  #
  CURRENT_DATE=$begofan
  CURRENT_DATE=$(substring $CURRENT_DATE 1 8)
  ANALYSIS_DATE=$(newdate 19${CURRENT_DATE} 6)
  ANALYSIS_DATE=$(substring $ANALYSIS_DATE 3 10)
  LAST_ANALYSIS_DATE=$endofan
  LAST_ANALYSIS_DATE=`substring $LAST_ANALYSIS_DATE 1 8`
  until [[ $ANALYSIS_DATE -gt $LAST_ANALYSIS_DATE ]] ; do
    if [[ !  -r ${VPPPATH}/RFL${ANALYSIS_DATE}00 ]] ; then
      # in case file is not there do a ecfile request
      ecfile -p $CFSALT/BFL${ANALYSIS_DATE}00 \
             -en get BFL${ANALYSIS_DATE}00
  # use appropriate decoder
      if [ $PRECISION = SINGLE ] ; then
         /home/ma/emos/wave/bin/enderac \
         -i RFL${ANALYSIS_DATE}00 -b BFL${ANALYSIS_DATE}00
      else
         ./enderac -i RFL${ANALYSIS_DATE}00 -b BFL${ANALYSIS_DATE}00
      fi

      /bin/rm BFL${ANALYSIS_DATE}00
      cp RFL${ANALYSIS_DATE}00 ${VPPPATH}/RFL${ANALYSIS_DATE}00
    else
      ln -s ${VPPPATH}/RFL${ANALYSIS_DATE}00 .
    fi
    CURRENT_DATE=$ANALYSIS_DATE
    CURRENT_DATE=`substring $CURRENT_DATE 1 8`
    ANALYSIS_DATE=`newdate 19${CURRENT_DATE} 6`
    ANALYSIS_DATE=`substring $ANALYSIS_DATE 3 10`
  done
else
  AS=0
fi
#.
#.=====================================================================
#.  Prepare MARS requests. Load and execute program GETINPT.
#.
# use MRFS for compilation
# cd $MRFSDIR
#.
cat >getinptb.f <<\EOFPRO
#include "getinptb.F"
#include "file_transfer.F"
EOFPRO
#.
#.compile getinptb
#.
frt $rp -X9 -c -Ps -o getinptb.o getinptb.f || exit
#.
frt $rp -Wl,-dy -lsocket -Of -Wv,-Of -Ps -X9 -o getinptb  getinptb.o \
    -l$WAMCRLIB $EMOSLIB $ECLIB $FDBLIB $MPELIB $MPLIB
#.
cat > fort.4 << EOF
RUN MODEL FROM $begofan TO $endoffo
PARAM=10U/10V,
REPRES=GG,
ANALYSIS PERIOD = $antime
FORECAST PERIOD = $fctime
WIND TIME STEP IN SECONDS = 0021600
ASSIMILATION FROM $begofan TO $endofan
ASSIMILATION YES
NEW FORECAST YES
STATUS STORM1= UNFINISHED
STATUS STORM2= UNFINISHED
STATUS STORM3= UNFINISHED
DATE FOR OUTPUT OF BOTH RESTART FILES = 0000000000
LAST DATE FOR SPECTRA FILE OUTPUT = 0000000000
EOF
#.
getinptb || {
  #.
  #. Abort on error.
  #.
   ls -lsa
   print - "\n\n\t\tGETINPT FAILED\n\n"
   print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
   exit 6
}
/bin/rm fort.4
#.
echo "\n\n\t\tSTART POSSIBLE MARS RETRIEVAL\n\n"
#.
rm tempout getinptb.* || true
#.
#.=====================================================================
#.  MARS retrieval.
#.=====================================================================
#.
#. Create filename for mars retrieval output
wind_filename=`awk '{print$1}' windfilename`
#.
echo $wind_filename
#.
#.check if mars retieval needed
#.
if [[ ! -r $VPPPATH/$wind_filename ]] ; then
  echo "mars wind file NOT online, start mars retrieval"
  mars fort.22 || {
    #.
    #. Abort on error.
    #.
      ls -lsa
      print - "\n\n\t\tMARS FAILED\n\n"
      print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
      exit 7
  }
  mv fort.80 fort.63
  cp fort.63 $VPPPATH/$wind_filename
  print - "\n\n\t\tSTART WAVE MODEL \n\n"
else
  print - "mars wind file already online"
  ln -s $VPPPATH/$wind_filename fort.63
  print - "\n\n\t\tSTART WAVE MODEL \n\n"
fi
#.
rm fort.22
#.
# cp $VPPPATH/WAMINFO fort.4
#endif
#.
#.
#.=====================================================================
#.  WAVE MODEL execution.
#.=====================================================================
#.
#. Common UBUF.
#.
if [[ ! -r $VPPPATH/ubufglou ]] ; then
  ecp ec:$DHSPATH/ubufglou $VPPPATH/ubufglou || {
    ls -lsa
    print - "\n\n\t\tPREPROC file ubufglou was not found\n\n"
    print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
    exit
  }
fi
ln -s $VPPPATH/ubufglou fort.8 
#.
#. Coarse grid output points, if this is a coarse grid run.
#.
#. assign -a $VPPPATH/bouccatu fort.9
#.
#. Fine grid input points, if this is a fine grid run.
#.
#. assign -a $VPPPAT1/boufcatu fort.10
#.
#.======================================================================
#. GET RESTART FILES
#.======================================================================
#. !!!! at this time we assume that the restart files have been saved
#. on disk in directory $VPPPATH
#.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#.
#.
#. Get last wind.
#.
if [[ !  -r $VPPPATH/LAW${begofan}00 ]] ; then
  ecp ec:$DHSPATH/LAW${begofan}00 $VPPPATH/LAW${begofan}00
  if [[ !  -r $VPPPATH/LAW${begofan}00 ]] ; then
    ls -lsa
    print - "LAW"${begofan}"00 file was not found"
    print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
    exit
  fi
fi
ln -s $VPPPATH/LAW${begofan}00 fort.15
ln -s $VPPPATH/LAW${begofan}00 LAW${begofan}00 
#.
#. Get initial spectrum.
#.
if [[ ! -r $VPPPATH/BLS${begofan}00 ]] ; then
  ecp ec:$DHSPATH/BLS${begofan}00 $VPPPATH/BLS${begofan}00
  if [[ ! -r $VPPPATH/BLS${begofan}00 ]] ; then
    ls -ls
    print - "BLS"${begofan}"00 file was not found"
    print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
    exit
  fi
fi
ln -s $VPPPATH/BLS${begofan}00 fort.12
ln -s $VPPPATH/BLS${begofan}00 BLS${begofan}00 
#.
#.======================================================================
#. Copy the WAM MODEL to the compiler input.
#.======================================================================
#.
#. Use MRFS for compilation.
#.
# cd $MRFSDIR
#.
cat > wamodel.src <<\EOFPRO
#endif
#include "wave_module.F"
#include "chief.F"
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
[ -f ${VPPPATH}/libwam.f.tar ] && {
  dd if=${VPPPATH}/libwam.f.tar of=libwam.f.tar bs=400k
  cd tmp
  tar xf ../libwam.f.tar
  cd ..
}
#.
[ -f ${VPPPATH}/libwam.a ] && {
  dd if=${VPPPATH}/libwam.a of=libwam.a bs=400k
  ar xv libwam.a wave_module
}
#.
export frtcntr=0
#.
(
  if [ "$REMAKE" != YES ] ; then
    if [ -f tmp/wave_module.fm ] ; then
      if cmp wave_module.f tmp/wave_module.fm ] ; then
        if ar tv libwam.a | grep wave_module ; then
          exit
        fi
      fi
    fi
  fi
  frt $rp -Of -M. -Wv,-Of -Pos -Z wave_module.l -X9 \
      -c wave_module.f
  mv wave_module.f wave_module.fm
  (( frtcnt += 1))
)
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
          exit
        fi
      fi
    fi
  fi
  frt $rp -Of -M. -Wv,-Of -Pos -Z ${file%.f}.l -X9 \
      -c -o ${file%.f}.o ${file} || print - ${file%.f}.l >> frt_error
  (( frtcnt += 1))
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
ecp -o libwam.* ec:$DHSPATH/.
#.
ls -la
#.
dd if=libwam.a     of=$VPPPATH/libwam.a     bs=400k
dd if=libwam.f.tar of=$VPPPATH/libwam.f.tar bs=400k
dd if=libwam.l.tar of=$VPPPATH/libwam.l.tar bs=400k
#.
chmod go+rx $VPPPATH/libwam.*
#.
if [ -r frt_error ] ; then
  cat $(cat frt_error)
  exit
fi
#.
if [[ ! -r $VPPLPATH/lib$WAMCRLIB.a ]] ; then
  ecp ec:$DHSLPATH/lib$WAMCRLIB.a $VPPLPATH/lib$WAMCRLIB.a || {
    ls -lsa
    print - "\n\n\t\tWAMCRLIB file was not found\n\n"
    print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
    exit 3
  }
fi
#.
#. LOAD LIBRARY
#.
ar vx $aux ifssig.o
ar vx $dum ishell.o
[ ! -a chief.o ] &&  ar xv libwam.a chief.o
#.
if [ $frtcntr = 0 -a -f $LTEMP/chief.coupled -a "$REMAKE" != YES  ] ; then
  ln -s  $LTEMP/chief.coupled chief
else
  frt $rp -Wl,-m,-P,-uMAIN__ -X9 \
    -o chief  chief.o ishell.o ifssig.o \
    libwam.a -lwamcr libwam.a \
     $ECLIB $EMOSLIB $FDBLIB  \
     -Wl,-dy -lsocket $MPELIB  $MPLIB || {
    #.
    #. Check return code.
    #.
     print - "\n\n\t\tLOAD FAILED\n\n"
     print - "\n\n\t\tPROGRAM WILL TERMINATE\n\n"
     exit 4
    }
  #.
  cp chief $LTEMP/chief.coupled
fi
#.
print - "\n\n\t\tWAVE MODEL WILL START\n\n"
#. Move executable back to disk.
# mv chief $TMPDIR 
# cd $TMPDIR 
#.
##################################################
# !!! the user input file has been moved to unit 3
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##################################################
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
#modify runtime buffers
#######################
fu12bf=4096 ; export fu12bf
fu15bf=344 ; export fu15bf
fu63bf=524 ; export fu63bf
fu90bf=53 ; export fu90bf
#.
chief || {
  #.
  #.  Abort on wave model failure.
  #.
  #.
  #.  Error: Copy program output to stdout and terminate.
  #.
   echo $?
   print - "\n\n\t\tWAVE MODEL FAILED\n\n"
   ls -l
   more logfile.*
   test -r split_table
   if [[ $? = 0 ]] ; then
     cp split_table $VPPPATH/split_table 
   fi
   exit
}
ls -lsa
cp fort.2* $VPPPATH/.
#.
#. Save split_table.
####################
#.
if [[ -r split_table ]] ; then
  cp split_table $VPPPATH/split_table 
fi
#.
#. Dump logfiles.
#################
#.
more logfile.*
#.
#.
#.======================================================================
#.END OF JOB.
#.======================================================================
#.
#endif
