#ifdef JOB
#.
#. PRESET TO GENERATE START FILES
#. FOR USE BY WAMODEL.
#.
#. HEINZ GUNTHER        ECMWF         JUNE 1990
#. B HANSEN             ECMWF         NOVEMBER 1995 MIRAGTION TO CLEARCASE
#.
#QSUB -r preset_CY18R4_test
#QSUB -q vpp700.normal
#QSUB -lP 0
#QSUB -lt 3000
#QSUB -lV 0mb
#QSUB -lM 700mb
#QSUB -eo
#QSUB -mb -me -mi
#QSUB
#.
######################################################################
# MAKE USE OF MEMORY RESIDENT FILE SYSTEM for TMPDIR
# TMPDIR=$MRFSDIR
######################################################################
#.
#!!! set variable COLDSTART !!!. By default it is set to yes, if set to no
# it means that the inital condition will be reconstructed from archived
# spectra and drag coefficient plus the usual surface wind. In that case,
# you must also specified the experiment version and class  of the archived
# data as well as the grid option regular/irregular grid (irgg=0/1). By default
# ALL but 0.5 grid should use irgg=0. Note that you can select irgg=0 for
# the 0.5 irregular grid, but this will result in unnecessary interpolation of
# the grib fields and it is best to retrieve the field on the same grid as
# created with.
#.
#.
COLDSTART=yes
expver_archive=1
class_archive=OD
#if resolution == 50
irgg=1
#else
irgg=0
#endif

#.
PATH=$PATH:.
cd $TMPDIR
set +v
#.
#include "path.h"
#.
if [ $PRECISION = SINGLE ] ; then
  export rp=-AD
  libselect -r 32 -d 64
  MPLIB=/opt/tools/lib/libmp2x.a
  MPELIB=/usr/local/lib/libmpe32.a
  #. !!!! MPELIB is selected for 32 bit real !!!!!!
  #################################################
elif [ $PRECISION = DOUBLE ] ; then
  libselect -r 64
  export rp=-Ad
  MPLIB=/opt/tools/lib/libmp2x.a
  MPELIB=/usr/local/lib/libmpe2.a
  #. !!!! MPELIB is selected for 64 bit real !!!!!!
  #################################################
else
  ls -lsa
  print - "\n\n\t\tPRECISION $PRECISION not valid\n\n"
  print - "\n\n\t\tUSE SINGLE or DOUBLE \n\n"
  print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
  exit 1
fi
#.
#.=====================================================================
#. Restore wamcrlib from DHS if lost.
#.=====================================================================
#.
if [[ ! -a  $VPPLPATH/lib$WAMCRLIB.a ]] ; then
  ecp ec:$DHSLPATH/lib$WAMCRLIB.a $VPPLPATH/lib$WAMCRLIB.a || {
    ls -lsa
    print - "\n\n\t\tWAMCRLIB file was not found\n\n"
    print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
    exit 3
  }
fi
#.
#.=====================================================================
#. General grid information and model constants.
#.=====================================================================
#.
if [[ ! -r $VPPPATH/gridglou ]] ; then
  ecp ec:$DHSPATH/gridglou $VPPPATH/gridglou || {
      ls -lsa
      print - "\n\n\t\tPREPROC file gridglou was not found\n\n"
      print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
      exit 2
  }
fi
cp $VPPPATH/gridglou fort.7 
#.
#. Common UBUF.
#.
if [[ "$COLDSTART" = "no" ]] ; then
  if [[ ! -r $VPPPATH/ubufglou ]] ; then
    ecp ec:$DHSPATH/ubufglou $VPPPATH/ubufglou || {
      ls -lsa
      print - "\n\n\t\tPREPROC file ubufglou was not found\n\n"
      print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
      exit
    }
  fi
  ln -s $VPPPATH/ubufglou fort.8
fi

# set restart file creation option
#
#.=====================================================================
#
if [[ "$COLDSTART" = "no" ]] ; then
### we will need to refine this as -1 should be use for the mediter.
### and for all OD global field prior to 961205
   date=`substring $begofan 1 6`
   if [[ $date -lt 961205 ]] ; then
     OPTI="-1"
   else
     OPTI="-2"
   fi
else
   OPTI=" 1"
fi

#if region=='s'
#. THE SWAMP CASE DOES NOT NEED EXTERNAL WINDS.
#else
#.
#.=====================================================================
#.  Prepare MARS requests. Load and execute program GETINPT.
#.
cat >getinptb.f <<\EOFPRO
#include "getinptb.F"
#include "file_transfer.F"
EOFPRO
# compile getinptb
#.
frt $rp -X9 -c -Z getinptb.lis -o getinptb.o getinptb.f
#.
frt $rp -Wl,-dy -lsocket -Ps -X9 -o getinptb getinptb.o -l$WAMCRLIB \
$EMOSLIB $ECLIB $FDBLIB $MPELIB $MPLIB
#.
cat > fort.4 <<- EOF
RUN MODEL FROM $begofan TO $endoffo 
PARAM=10U/10V,
REPRES=GG,
ANALYSIS PERIOD = $antime
FORECAST PERIOD = $fctime
WIND TIME STEP IN SECONDS = 21600
ANALYSIS FROM $begofan TO $endofan 
ASSIMILATION NO
NEW FORECAST YES
STATUS STORM1= UNFINISHED
STATUS STORM2= UNFINISHED
STATUS STORM3= UNFINISHED
DATE FOR OUTPUT OF BOTH RESTART FILES = $outofrf
LAST DATE FOR SPECTRA FILE OUTPUT = $outof2d
EOF
#.
getinptb || {
  #.
  #. Abort on error.
  #.
   ls -lsa
   print - "\n\n\t\tGETINPT FAILED\n\n"
   print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
   cat getinptb.lis
   exit 6
  }
#.
print - "\n\n\t\tSTART POSSIBLE MARS RETRIEVAL\n\n"
#.
rm tempout
rm getinptb.*
#.
#.=====================================================================
#.  MARS retrieval.
#.=====================================================================
#.
#.
#.create filename for mars retrieval output
#.
wind_filename=`awk '{print$1}' windfilename`
#.
print - $wind_filename
#.
#.check if mars retieval needed
#.
if [[ ! -r $VPPPATH/$wind_filename ]] ; then
   print -  "mars wind file NOT online, start mars retrieval"
   mars fort.22 || {
    #.
    #. Check MARS error code.
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
 cp  $VPPPATH/$wind_filename fort.63
 print - "\n\n\t\tSTART PRESET MODEL \n\n"
fi
#.
rm fort.22
#.
#endif
#
# Retrieve the drag coefficient if the initial data are contructed
# from known archive spectra , drag coefficient and surface wind speed
# as analysed by the wave model (if available, new parameter245)
# and link an existing BLS file to fort.12. This will change once we
# retrieve the spectra straight out of the archive.
#################################################################### 

if [[ "$COLDSTART" = "no" ]] ; then

#if resolution == 150
  GRID="1.5/1.5"
  AREA="81./ 0./ -81./358.5"
#elif region == 'g' &&  resolution == 50
  GRID="0.5/0.5"
  AREA="81./ 0./ -81./359.5"
#elif region == 'g' &&  resolution == 300 
  GRID="3.0/3.0"
  AREA="72./ 0./ -69./357.0"
#elif region == 'm' && resolution == 25
  GRID="0.25/0.25"
  AREA="66./ -6./ 30./42."
#endif

#
# Get the WAVE spectra 
# ====================
#
SPECDATE=`substring $begofan 1 6`
SPECTIME=`substring $begofan 7 8`
if [[ $SPECTIME -ne 12 && $SPECDATE -lt 980120 ]] ; then
   echo 'THERE IS NO SPECTRA IN THE ARCHIVE OUTSIDE 12Z'
   echo 'PLEASE CHANGE INITIAL DATE ACCORDINGLY'
   exit 8
fi
if [[ "$irgg" = "0" && $SPECDATE -ge 961205 ]] ; then
   echo 'THERE IS NO SPECTRA IN THE ARCHIVE ON A REGULAR GRID  '
   echo 'AFTER 961205'
   echo 'INTERPOLATION WILL BE USED TO CREATE ONE.'
fi
if [[ "$irgg" = "1" && $SPECDATE -lt 961205 ]] ; then
   echo 'THERE IS NO SPECTRA IN THE ARCHIVE ON AN IRREGULAR GRID  '
   echo 'BEFORE 961205'
   echo 'INTERPOLATION WILL BE USED TO CREATE ONE.'
fi

sfilename=GRS${begofan}00

if [[ !  -r $VPPPATH/${sfilename} ]] ; then
  # in case file is not there do a mars request
  # GET WAVE SPECTRA 
  #=================
  SPECDATE=`substring $begofan 1 6`
  SPECTIME=`substring $begofan 7 8`

  if [[ $OPTI -eq -2 ]] ; then
  cat > myrequest <<-EOF
	RETRIEVE,LEVTYPE=SFC,TYPE=AN,
	STREAM=WAVE,
	CLASS=$class_archive,
	EXPVER=$expver_archive,
	LEVEL=OFF,
	PARAM=250,
	REPRES=LL,
	STEP=0,
	TIME=$SPECTIME,
	DATE=$SPECDATE,
	DOMAIN=N,
	TARGET="${sfilename}"

	RETRIEVE,
	DOMAIN=S

	EOF
  else
  cat > myrequest <<-EOF
	RETRIEVE,LEVTYPE=SFC,TYPE=AN,
	STREAM=WAVE,
	CLASS=$class_archive,
	EXPVER=$expver_archive,
	LEVEL=OFF,
	PARAM=250,
	REPRES=LL,
	STEP=0,
	TIME=$SPECTIME,
	DATE=$SPECDATE,
	DOMAIN=G,
	TARGET="${sfilename}"
	EOF
  fi
  #.
  #.DO THE MARS REQUEST.
  #=====================
  mars myrequest || {
    print - "\n\n\t\t NO WAVE SPECTRA 250 FROM MARS, NOTHING IS DONE \n\n"
    exit 9
  }
  cp ${sfilename} $VPPPATH/${sfilename}
else
  ln -s $VPPPATH/${sfilename} ${sfilename}
fi
#
# Get drag coefficient with waves 
# ===============================
#
if [[ !  -r $VPPPATH/CD$begofan ]] ; then
  # in case file is not there do a mars request
  # GET DRAG COEFFICIENT WITH WAVE FROM MARS
  #i========================================
  CDDATE=`substring $begofan 1 6`
  CDTIME=`substring $begofan 7 8`
if [[ "$irgg" = "1" ]] ; then
  cat > myrequest <<-EOF
	RETRIEVE,LEVTYPE=SFC,TYPE=AN,
	STREAM=WAVE,
	CLASS=$class_archive,
	EXPVER=$expver_archive,
	LEVEL=OFF,
	PARAM=233,
	REPRES=LL,
	STEP=0,
	TIME=$CDTIME,
	DATE=$CDDATE,
	AR=$AREA,
	TARGET="fort.91"
	EOF
else
  cat > myrequest <<-EOF
	RETRIEVE,LEVTYPE=SFC,TYPE=AN,
	STREAM=WAVE,
	CLASS=$class_archive,
	EXPVER=$expver_archive,
	LEVEL=OFF,
	PARAM=233,
	REPRES=LL,
	STEP=0,
	TIME=$CDTIME,
	DATE=$CDDATE,
	GRID=$GRID,
	AR=$AREA,
	TARGET="fort.91"
	EOF
fi
  #.
  #.DO THE MARS REQUEST.
  #=====================
  mars myrequest || {
    print - "\n\n\t\t NO DRAG COEFFICIENT 233 FROM MARS, NOTHING IS DONE \n\n"
    exit 8
  }
  cp fort.91 $VPPPATH/CD$begofan
else
  ln -s $VPPPATH/CD$begofan fort.91
fi

#
# Get surface wind speed as modified by the wave model analysis
# =============================================================
#
if [[ !  -r $VPPPATH/U10$begofan ]] ; then
  # in case file is not there do a mars request
  # GET SURFACE WIND WITH WAVE FROM MARS
  #i=====================================
  CDDATE=`substring $begofan 1 6`
  CDTIME=`substring $begofan 7 8`
if [[ "$irgg" = "1" ]] ; then
  cat > myrequest <<-EOF
	RETRIEVE,LEVTYPE=SFC,TYPE=AN,
	STREAM=WAVE,
	CLASS=$class_archive,
	EXPVER=$expver_archive,
	LEVEL=OFF,
	PARAM=245,
	REPRES=LL,
	STEP=0,
	TIME=$CDTIME,
	DATE=$CDDATE,
	AR=$AREA,
	TARGET="fort.92"
	EOF
else
  cat > myrequest <<-EOF
	RETRIEVE,LEVTYPE=SFC,TYPE=AN,
	STREAM=WAVE,
	CLASS=$class_archive,
	EXPVER=$expver_archive,
	LEVEL=OFF,
	PARAM=245,
	REPRES=LL,
	STEP=0,
	TIME=$CDTIME,
	DATE=$CDDATE,
	GRID=$GRID,
	AR=$AREA,
	TARGET="fort.92"
	EOF
fi
  #.
  #.DO THE MARS REQUEST.
  #=====================
  mars myrequest || {
    print - "+++++++++++++++++++++++++++++++++++++++++++++++++++"
    print - "NO WIND SPEED 245 FROM MARS, program will continue "
    print - "and will use the atmospheric model winds instead   "
    print - "+++++++++++++++++++++++++++++++++++++++++++++++++++"
    \rm fort.92
  }
  if [[ -r $VPPPATH/U10$begofan ]] ; then
    cp fort.92 $VPPPATH/U10$begofan
  fi
else
  ln -s $VPPPATH/U10$begofan fort.92
fi


# end coldstart=no
fi

#.
#.======================================================================
#.  Copy the source code to the compiler input.
#.======================================================================
#.
cat >preset.f <<\EOFPRO
#endif
#include "wave_module.F"
#include "preset.F"
#include "airsea.F"
#include "file_transfer.F"
#include "grstname.F"
#include "getspec.F"
#include "getwnd.F"
#include "../Alt/i_get_unit.F"
#include "locint.F"
#include "makeblos.F"
#include "mpdecomp.F"
#include "mpdistribfl.F"
#include "mpgatherfl.F"
#include "mpgatherscfld.F"
#include "mstart.F"
#include "notim.F"
C  #ifdef 
C  #include "openfil_unicos.F"
C  #endif
#include "peak.F"
#include "prealloc_file.F"
#include "prewind.F"
#include "readfl.F"
#include "readgrib.F"
#include "readgrs.F"
#include "readpre.F"
#if region=='s'
#include "readwnd_rswamp2.F"
#else
#include "readwnd_mars.F"
#endif
#include "savspec.F"
#include "savstress.F"
#include "spectra.F"
#include "splitbl.F"
#include "spr.F"
#include "timin.F"
#include "../Alt/uraopc.F"
#include "wamwnd.F"
#include "writefl.F"
#include "writestress.F"
#ifdef JOB
EOFPRO
#.
#.======================================================================
#. COMPILE AND LOAD.
#.======================================================================
#.
#. COMPILE SOURCE CODE
#.
frt $rp -Of -Wv,-Of -Pos -Z preset.lis -X9 -c -o preset.o preset.f || {
  #.
  #. Abort on error.
  #.
    print - "\n\n\t\tCOMPILE FAILED\n\n"
    print - "\n\n\t\tPROGRAM WILL TERMINATE\n\n"
    cat preset.lis
    print - "\n\n\t\tCOMPILE FAILED\n\n"
    print - "\n\n\t\tPROGRAM WILL TERMINATE\n\n"
    exit 10
  }
#.
ar vx $dum ishell.o
frt $rp -Wl,-dy -lsocket -Of -Wv,-Of -Ps -X9  -o preset preset.o ishell.o \
    -l$WAMCRLIB $EMOSLIB $ECLIB $FDBLIB $MPELIB $MPLIB || {
  #.
  #.
  #. Abort on error.
  #.
    print - "\n\n\t\tLOAD FAILED\n\n"
    print - "\n\n\t\tPROGRAM WILL TERMINATE\n\n"
    exit 11
  }
#.
#.======================================================================
#. RUN THE PROGRAM AND FORMAT THE OUTPUT.
#.======================================================================
#.
preset <<EOD_for_preset || { print - "\n\n\t\tPRESET FAILED\n\n" ; ls -ltra ; exit ; }
#if region=='s'
#include "input_preset_rswamp2.F"
#elif region=='m'
#include "input_preset_rmedite.F"
#else
#include "input_preset_rglobal.F"
#endif
EOD_for_preset
#.
/bin/rm fort.4
#.
ls -la
nobls=0
nolaw=0
[[ ! -f $VPPPATH/BLS${begofan}00 ]] && nobls=1
[[ ! -f $VPPPATH/LAW${begofan}00 ]] && nolaw=1

if [[ -f fort.15 ]] ; then
  ecp fort.15 ec:$DHSPATH/.
  [[ $nolaw = 1 ]] && cp fort.15 $VPPPATH/LAW${begofan}00
elif [[ -f LAW${begofan}00 ]] ; then
  ecp -o LAW${begofan}00 ec:$DHSPATH/.
  [[ $nolaw = 1 ]] && cp LAW${begofan}00 $VPPPATH/.
endif
#.
if [[ -f fort.12 ]] ; then
  ecp fort.12 ec:$DHSPATH/.
  [[ $nobls = 1 ]] && cp fort.12 $VPPPATH/BLS${begofan}00
elif [[ -f BLS${begofan}00 ]] ; then
  ecp -o BLS${begofan}00 ec:$DHSPATH/.
  [[ $nobls = 1 ]] && cp BLS${begofan}00 $VPPPATH/.
endif
#.
print - "\n\n\t\tPRESET ALL DONE\n\n"
#.
#.======================================================================
#. START FIELDS GENERATED AND SAVED.
#. END OF JOB PRESET.
#.======================================================================
#.
#endif
