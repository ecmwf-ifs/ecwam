#ifdef JOB
#.
#. PRESET TO GENERATE and/or RETRIEVE  START FILES
#. FOR USE BY WAMODEL.
#.
#. HEINZ GUNTHER        ECMWF         JUNE 1990
#. B HANSEN             ECMWF         NOVEMBER 1995 MIRAGTION TO CLEARCASE
#. J. BIDLOT            ECMWF         1998 : use of grib spectra
#.                                           all data retrieval for wam run
#.                                           must be done here.
#. 
#QSUB -q vpp700.rdx04_io
#QSUB -r preset
#QSUB -lP 0
#QSUB -lt 3000
#QSUB -lV 0mb
#QSUB -lM 700mb
#QSUB -eo
#QSUB -mb -mi
#QSUB
#.
#!!! set variable COLDSTART !!!. 
# COLDSTART will be set to yes most of the time as preset should be used
# to generate all initial conditions for the wave model.
# When COLDSTART is set to no,
# it means that the inital condition will be reconstructed from archived
# spectra and drag coefficient plus the usual surface winds. In that case,
# you must also specified the experiment version and class of the archived
# data as well as the grid option regular/irregular grid (irgg=0/1). By default
# ALL but 0.5 and 0.25 grids should use irgg=0. Note that you can select irgg=1
# for the irregular grid, but this will result in unnecessary interpolation of
# the grib fields and it is best to retrieve the field on the same grid as
# created with. YOU must also specify whether the wave spectra are of the
# new type (specparam=251) or the old type (specparam=250). Finally set 
# luwave=yes if wind speed analysed by the wave model (parameter 245) is
# available (for OD, there was none until the coupled global model became
# operational (29 June 1998) and until the introduction of the new limited
# area model (27 October 1998).
#.
#
COLDSTART=yes
#
expver_archive=0001
class_archive=OD
specparam=251
luwave=yes
#.
#if region=='m'
domain=M
#if resolution == 25
irgg=1
#else
irgg=0
#endif
#else
domain=G
#if resolution == 50
irgg=1
#else
irgg=0
#endif
#endif
#.
PATH=$PATH:.
#.
#include "path.h"
#.
# move to working directiory
############################
cd $WDIR
mkdir $WORKDIR 
cd $WORKDIR 

#. determine if preset has to be run or only data retrieval is necessary
#. 
if [[ $COLDSTART = yes || ($COLDSTART != yes && $specparam = 250) ]] ; then
 print - " PRESET WILL BE RUN "
 lpreset=yes
  #.
  #. Restore WAMLIB from DHS if lost.
  #.=====================================================================
  #.
  if [[ ! -a  $LIBS/$WAMLIB.a ]] ; then
    ecp ec:$DHSLPATH/$WAMLIB.a $LIBS/$WAMLIB.a || {
      ls -lsa
      print - "\n\n\t\tWAMLIB file was not found\n\n"
      print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
      exit 3
    }
  fi
  #.
  #.=====================================================================
  #. General grid information and model constants.
  #.=====================================================================
  #.
  if [[ ! -r $ADIR/gridglou_${grid} ]] ; then
    ecp ec:$DHSPATH/gridglou_${grid} $ADIR/gridglou_${grid} || {
        ls -lsa
        print - "\n\n\t\tPREPROC file gridglou was not found\n\n"
        print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
        exit 2
    }
  fi
  ln -s  $ADIR/gridglou_${grid} fort.7 
  #.
  #. Common UBUF.
  #.
  if [[ ! -r $ADIR/ubufglou_${grid} ]] ; then
    ecp ec:$DHSPATH/ubufglou_${grid} $ADIR/ubufglou_${grid} || {
      ls -lsa
      print - "\n\n\t\tPREPROC file ubufglou was not found\n\n"
      print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
      exit
    }
  fi
  ln -s $ADIR/ubufglou_${grid} fort.8

fi
#.
# set restart file creation option
#.================================
#
if [[ "$COLDSTART" = "no" ]] ; then
### we will need to refine this as 1 should be use for the limited area model
### and for all OD global field prior to 19961205
   date=`substring $begofrn 1 8`
   if [[ $date -lt 19961205 ]] ; then
     OPTI="1"
   else
     OPTI="2"
   fi
   lgribin=T
else
   OPTI=" 1"
   lgribin=F
fi

#if region=='s'
#. THE SWAMP CASE DOES NOT NEED EXTERNAL WINDS.
#else
#.
#.=====================================================================
#.  Prepare MARS requests for winds. Load and execute program GETINPT.
#.
if [ ! -f $LIBS/getinptb ] ; then
  # extract
  ar xv $LIBS/$WAMLIB.a getinptb.o 

  # load
  frt $rp -Wl,-dy -lsocket -Ps -X9 -o $LIBS/getinptb getinptb.o \
  -lwam $EMOSLIB $ECLIB $FDBLIB $MPELIB $MPLIB || {
  #.
  #.
  #. Abort on error.
  #.
    print - "\n\n\t\tgetinptb LOAD FAILED\n\n"
    print - "\n\n\t\tPROGRAM WILL TERMINATE\n\n"
    exit 2
  }
  #.
  rm getinptb.*
  cd $WDIR/$WORKDIR 
fi
#.
cat > fort.4 <<- EOF
RUN MODEL FROM $begofrn TO $endofrn 
PARAM=10U/10V, CLASS=OD,
REPRES=GG,
ANALYSIS PERIOD = $antime
FORECAST PERIOD = $fctime
WIND TIME STEP IN SECONDS = 21600
ANALYSIS FROM $begofrn TO $begoffo 
ASSIMILATION NO
NEW FORECAST YES
STATUS STORM1= UNFINISHED
STATUS STORM2= UNFINISHED
STATUS STORM3= UNFINISHED
DATE FOR OUTPUT OF BOTH RESTART FILES = $outofrf
LAST DATE FOR SPECTRA FILE OUTPUT = $outof2d
EOF
#.
# gridglou is used by getinptb
  if [[ ! -r $ADIR/gridglou_${grid} ]] ; then
    ecp ec:$DHSPATH/gridglou_${grid} $ADIR/gridglou_${grid} || {
        ls -lsa
        print - "\n\n\t\tPREPROC file gridglou was not found\n\n"
        print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
        exit 2
    }
  fi
  rm -f fort.7
  ln -s  $ADIR/gridglou_${grid} fort.7 

$LIBS/getinptb || {
  #.
  #. Abort on error.
  #.
   ls -lsa
   print - "\n\n\t\tGETINPT FAILED\n\n"
   print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
   exit 6
  }
#.
print - "\n\n\t\tSTART POSSIBLE MARS RETRIEVAL\n\n"
#.
#.
#.=====================================================================
#.  MARS retrieval.
#.=====================================================================
#.
#.
#.create filename for mars retrieval output
#.
wind_filename=sfcwind_${grid}_${begofrn}_${begoffo}_${endofrn}
#.
print - $wind_filename
#.
#.check if mars retieval needed
#.
if [[ ! -r $WDIR/$wind_filename ]] ; then
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
   mv fort.80 sfcwindin 
   cp sfcwindin $WDIR/$wind_filename
   print - "\n\n\t\tSTART WAVE MODEL \n\n"
else
 print - "mars wind file already online"
 ln -s $WDIR/$wind_filename sfcwindin 
 print - "\n\n\t\tSTART PRESET MODEL \n\n"
fi
#.
rm fort.22
#.
#endif
#.
#.
#.check if need SST for ice boundaries (!!!! we are still using the OD SST!!!)
#######################################
#if iceedge=='y'
#.
#.check if SST file does exist
#.
if [[ !  -r $ADIR/kelv${begofrn}_${grid} ]] ; then
  # in case file is not there do a mars request
  # GET SURFACE TEMPERATURES FROM MARS
  #===================================
  ICEDATE=`substring $begofrn 1 8`
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
	TARGET="sfctempin"
	EOF
  #.
  #.DO THE MARS REQUEST.
  #=====================
  mars myrequest || {
    print - "\n\n\t\t NO SURFACE TEMPS FROM MARS, NOTHING IS DONE \n\n"
    exit 19
  }
  mv sfctempin $ADIR/kelv${begofrn}_${grid}
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
  #. Note the location of the data in ECFS has recently changed (June 1998) 
  #.
  #. as well as the way they are kept in the archive (see Bjorn)
  #.
  #.====================================================================
  ##
  AS=1
  DHSWAVE=/emos/0001/wave
  #
  # get 64 bit decoder from ops ( there is not one readily available
  # for 32 bit)
  #
   if [ $PRECISION = DOUBLE ] ; then
      cp /vpp700/emos_data/0001/bin/global50/enderac .
   else
      echo 'the following decoder does not work '
      echo '/home/ma/emos/wave/bin/enderac'
      exit 1
   fi
  #
  ANALYSIS_DATE=$(substring $begofrn 1 10)
  ANALYSIS_DATE=$(newdate ${ANALYSIS_DATE} 6)
  LAST_ANALYSIS_DATE=`substring $endofrn 1 10`
  #
  until [[ $ANALYSIS_DATE -gt $LAST_ANALYSIS_DATE ]] ; do
    if [[ !  -r ${ADIR}/RFL${ANALYSIS_DATE}00 ]] ; then
      # in case file is not there do a ecfs request
      yyyymm=$(substring ${ANALYSIS_DATE} 1 6)
      ddhh=$(substring ${ANALYSIS_DATE} 7 10)
      ecp ec:$DHSWAVE/$yyyymm/$ddhh/BFL BFL${ANALYSIS_DATE}00 
  #
  # use appropriate decoder
  #
      if [ $PRECISION = SINGLE ] ; then
           echo 'the following decoder does not work '
           echo '/home/ma/emos/wave/bin/enderac'
           exit 1
##         /home/ma/emos/wave/bin/enderac \
##         -i RFL${ANALYSIS_DATE}00 -b BFL${ANALYSIS_DATE}00
      else
         ./enderac -i RFL${ANALYSIS_DATE}00 -b BFL${ANALYSIS_DATE}00
      fi

      rm -f BFL${ANALYSIS_DATE}00
      mv RFL${ANALYSIS_DATE}00 ${ADIR}/RFL${ANALYSIS_DATE}00
    fi
    ANALYSIS_DATE=`newdate ${ANALYSIS_DATE} 6`
  done
else
  AS=0
fi

if [[ "$COLDSTART" = "no" ]] ; then
#
# Retrieve the drag coefficient if the initial data are contructed
# from known archive spectra , drag coefficient and surface wind speed
# as analysed by the wave model (if available, new parameter245)
#################################################################### 
#
# Get the WAVE spectra 
# ====================
#
SPECDATE=`substring $begofrn 1 8`
SPECTIME=`substring $begofrn 9 10`
if [[ $SPECTIME -ne 12 && $SPECDATE -lt 19980120 ]] ; then
   echo 'THERE IS NO SPECTRA IN THE ARCHIVE OUTSIDE 12Z'
   echo 'PLEASE CHANGE INITIAL DATE ACCORDINGLY'
   exit 8
fi
if [[ "$irgg" = "0" && $SPECDATE -ge 19961205 ]] ; then
   echo 'THERE IS NO SPECTRA IN THE ARCHIVE ON A REGULAR GRID  '
   echo 'AFTER 19961205'
   echo 'INTERPOLATION WILL BE USED TO CREATE ONE.'
fi
if [[ "$irgg" = "1" && $SPECDATE -lt 19961205 ]] ; then
   echo 'THERE IS NO SPECTRA IN THE ARCHIVE ON AN IRREGULAR GRID  '
   echo 'BEFORE 19961205'
   echo 'INTERPOLATION WILL BE USED TO CREATE ONE.'
fi

if [[ $specparam = 250 ]] ; then
  sfilename=GRS${begofrn}00
else
  sfilename=SGS${begofrn}00
fi

if [[ !  -r $WDIR/${sfilename} ]] ; then
  # in case file is not there do a mars request
  # GET WAVE SPECTRA 
  #=================
  SPECDATE=`substring $begofrn 1 8`
  SPECTIME=`substring $begofrn 9 10`

  if [[ $specparam = 250 ]] ; then

    if [[ $OPTI -eq 2 ]] ; then
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
		DOMAIN=${domain},
		STEP=0,
		TIME=$SPECTIME,
		DATE=$SPECDATE,
		DOMAIN=${domain},
		TARGET="${sfilename}"
		EOF
    fi
  else
    if [[ "$irgg" = "1" ]] ; then
      cat > myrequest <<-EOF
		RETRIEVE,LEVTYPE=SFC,TYPE=AN,
		STREAM=WAVE,
		CLASS=$class_archive,
		EXPVER=$expver_archive,
		LEVEL=OFF,
		PARAM=251,
		REPRES=LL,
		DOMAIN=${domain},
		STEP=0,
		TIME=$SPECTIME,
		DATE=$SPECDATE,
		DIRECTION=1/to/${nang},
		FREQUENCY=1/to/${nfre},
		TARGET="${sfilename}"
		EOF
    else
      cat > myrequest <<-EOF
		RETRIEVE,LEVTYPE=SFC,TYPE=AN,
		STREAM=WAVE,
		CLASS=$class_archive,
		EXPVER=$expver_archive,
		LEVEL=OFF,
		PARAM=251,
		REPRES=LL,
		DOMAIN=${domain},
		STEP=0,
		TIME=$SPECTIME,
		DATE=$SPECDATE,
		DIRECTION=1/to/${nang},
		FREQUENCY=1/to/${nfre},
		GRID=$GRID,
		AR=$AREA,
		TARGET="${sfilename}"
		EOF
     fi
  fi
  #.
  #.DO THE MARS REQUEST.
  #=====================
  mars myrequest || {
    print - "\n\n\t\t NO WAVE SPECTRA FROM MARS, NOTHING IS DONE \n\n"
    exit 9
  }
  mv ${sfilename} $WDIR/${sfilename}
  ln -s $WDIR/${sfilename} ${sfilename}
else
  ln -s $WDIR/${sfilename} ${sfilename}
fi
#
# Get drag coefficient with waves for wam model run 
# ===============================
#
if [[ !  -r $ADIR/CD${begofrn}_${grid} ]] ; then
  # in case file is not there do a mars request
  # GET DRAG COEFFICIENT WITH WAVE FROM MARS
  #i========================================
  CDDATE=`substring $begofrn 1 8`
  CDTIME=`substring $begofrn 9 10`
if [[ "$irgg" = "1" ]] ; then
  cat > myrequest <<-EOF
	RETRIEVE,LEVTYPE=SFC,TYPE=AN,
	STREAM=WAVE,
	CLASS=$class_archive,
	EXPVER=$expver_archive,
	LEVEL=OFF,
	PARAM=233,
	REPRES=LL,
	DOMAIN=${domain},
	STEP=0,
	TIME=$CDTIME,
	DATE=$CDDATE,
	TARGET="cdwavein"
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
	DOMAIN=${domain},
	STEP=0,
	TIME=$CDTIME,
	DATE=$CDDATE,
	GRID=$GRID,
	AR=$AREA,
	TARGET="cdwavein"
	EOF
fi
  #.
  #.DO THE MARS REQUEST.
  #=====================
  mars myrequest || {
    print - "\n\n\t\t NO DRAG COEFFICIENT 233 FROM MARS, NOTHING IS DONE \n\n"
    exit 8
  }
  mv cdwavein $ADIR/CD${begofrn}_${grid}
fi

rm -f cdwavein
ln -s $ADIR/CD${begofrn}_${grid} cdwavein

#
# Get surface wind speed as modified by the wave model analysis for wam run
# =============================================================
#
if [[ !  -r $ADIR/U10${begofrn}_${grid} ]] ; then
  if [[ $luwave  = yes ]] ; then
    # in case file is not there do a mars request
    # GET SURFACE WIND WITH WAVE FROM MARS
    #i=====================================
    CDDATE=`substring $begofrn 1 8`
    CDTIME=`substring $begofrn 9 10`
    if [[ "$irgg" = "1" ]] ; then
      cat > myrequest <<-EOF
		RETRIEVE,LEVTYPE=SFC,TYPE=AN,
		STREAM=WAVE,
		CLASS=$class_archive,
		EXPVER=$expver_archive,
		LEVEL=OFF,
		PARAM=245,
		REPRES=LL,
		DOMAIN=${domain},
		STEP=0,
		TIME=$CDTIME,
		DATE=$CDDATE,
		TARGET="uwavein"
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
		DOMAIN=${domain},
		STEP=0,
		TIME=$CDTIME,
		DATE=$CDDATE,
		GRID=$GRID,
		AR=$AREA,
		TARGET="uwavein"
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
    }
  fi
  cp uwavein $ADIR/U10${begofrn}_${grid}
fi
#.
# end coldstart=no
fi

if [[ $lpreset = yes ]] ; then
#.
#.
#.===========
#. LOAD.
#.===========
#.
  # extract
  ar xv $LIBS/$WAMLIB.a preset.o 
  ar vx $dum ishell.o

frt $rp -Wl,-dy -lsocket -Of -Wv,-Of -Ps -X9  -o preset preset.o ishell.o \
    -lwam $EMOSLIB $ECLIB $FDBLIB $MPELIB $MPLIB || {
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
#.
# end on lpreset
fi
#.
print - "\n\n\t\tPRESET ALL DONE\n\n"
#.
cd $WDIR
\rm -rf $WORKDIR 
#.
#.======================================================================
#. START FIELDS GENERATED AND SAVED.
#. END OF JOB PRESET.
#.======================================================================
#.
# cd $VPPPATH
#endif
