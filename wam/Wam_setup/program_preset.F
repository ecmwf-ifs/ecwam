#ifdef JOB
#.
#. PRESET IS USED TO GENERATE AND/OR RETRIEVE START FILES
#. FOR USE BY WAMODEL.
#.
#QSUB -q vpp700.normal
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
#################################
#.
COLDSTART=no
#.
# when COLDSTART is set to yes, preset will be used to generate all
# initial conditions for a run (usually used when first starting a set
# of runs with WAM). The initial spectra are obtained from the Jonswap
# parametrisation. The input namelist NALINE (see below) determines how
# the parameters for the Jonswap parametrisation are determined. When
# IOPTI is 1 or 2 then the parameters are specified using a fetch law and
# the local wind provided as input as well. When IOPTI is 0 then the
# parameters are as specified in the namelist (see ALFA, FM, GAMMA, SA, SB,
# THETA). The fetch used is also specified as input (FETCH). The difference
# between IOPTI=1 and IOPTI=2 is no longer inforced since there is a 
# requirement that the minimum wind speed is 1m/s.
#.
# When COLDSTART is set to no and LGRIBOUT to T (see below),
# it means that the inital condition will be reconstructed from archived
# spectra and drag coefficient plus the usual surface winds. In that case,
# you must also specified the experiment version and class of the archived
# data and its domain (G or M).
# It is assumed that the data were created with this latest version
# of the software and therefore all grib data fields are available.
# Note that this job does not bring online the binary files which are used
# instead of grib data as input to the wave model run (i.e when WAM is
# run with LGRIBIN=T see below)
#.
expver_archive=0001
class_archive=OD
domain=G
#.
#.
#. Type of input files !!!!
###########################
#.
#. The wave spectra and field needed to specify the initial wind and sfc stress
#  fields can either be in grib format or in binary forms. The binary format
#  can either be the pure binary form as created by a write statement in 
#  fortran or it can be created using the pbio software in use at ECMWF.
#  For this reason, one has to specify what format will be used to create
#  the initial data with preset and one has also to make sure that the same
#  format is specified in the namelist for the input fields for a wave model
#  run. Usually the same format is also used for the output of the wave model
#  run as some of the output is generally the input the the next run.

#  So in summary the following setting should specified in the input namelist
#  for PRESET            and   for WAMODEL 
#.
#  LGRIBOUT=T                  LGRIBIN=T  and LGRIBOUT=T
#.
#  LGRIBOUT=F and LPBIOOUT=T   LGRIBIN=F, LPBIOIN=T, GRIBOUT=F and LPBIOOUT=T
#.
#  LGRIBOUT=F and LPBIOOUT=F   LGRIBIN=F, LPBIOIN=F, GRIBOUT=F and LPBIOOUT=F
#.
#. When the the grib format is used LPBIOOUT and LPBIOIN have no meaning.
#. The grib format is commonly used at ECMWF the other option are feature from
#. the past.
#.
#.
#.
PATH=$PATH:.
#.
#include "wkstation_vpp_share.h"
#include "path.h"
#.
# move to working directiory
############################
cd $WDIR
mkdir $WORKDIR 
cd $WORKDIR 

#. Get initial data files for preset if it is run 
#. ###############################################
#. 
if [[ $COLDSTART = yes ]] ; then
  print - " PRESET WILL BE RUN "

  #. check is library is present
  if [[ ! -a  $LIBS/$WAMLIB.a ]] ; then
      ls -lsa
      print - "\n\n\t\tWAMLIB file was not found\n\n"
      print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
      exit 3
  fi
  #.
  #. General grid information and model constants.
  #.
  if [[ ! -r $ADIR/gridglou_${grid} ]] ; then
       ls -lsa
       print - "\n\n\t\tPREPROC file gridglou was not found\n\n"
       print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
       exit 2
  fi
  ln -s  $ADIR/gridglou_${grid} fort.7 
  #.
  #. Advection tables.
  #.
  if [[ ! -r $ADIR/ubufglou_${grid} ]] ; then
      ls -lsa
      print - "\n\n\t\tPREPROC file ubufglou was not found\n\n"
      print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
      exit
  fi
  ln -s $ADIR/ubufglou_${grid} fort.8
fi
#.

#if region=='s'
#. THE SWAMP CASE DOES NOT NEED EXTERNAL WINDS.
#else
#.
#.=====================================================================
#.  Prepare MARS requests for winds. Load and execute program GETINPT.
#.  GETINPT is intended to be quite general and can be replaced
#.  by a mars request if the latter is always the same
#.
#. copy library from workstation

cp $WKLIBS/$WAMLIB.a $LIBS/$WAMLIB.a


if [ ! -f $LIBS/getinptb ] ; then
  # extract
  ar xv $LIBS/$WAMLIB.a getinptb.o 

  # load

  frt $rp -Wl,-dy -lsocket -Ps -X9 -o $LIBS/getinptb getinptb.o \
  -lwam ${mpllib} ${emoslib} ${eclib} ${fdblib} \
  ${mpilib} ${mp2lib} || {
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
#. input for GETINPT which is uded to determined the mars request
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
        ls -lsa
        print - "\n\n\t\tPREPROC file gridglou was not found\n\n"
        print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
        exit 2
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
#. Retrieve sea ice fraction or sea surface temperature field
#. for sea ice mask
ciparam=CI
if [[ ! -r $ADIR/SSTDATA${begofrn}_${grid} ]] ; then
    CIDATE=`substring $begoffo 1 8`
    CITIME=`substring $begoffo 9 10`
    cat > myrequest <<-EOF
	RETRIEVE,LEVTYPE=SFC,TYPE=AN,
	STREAM=DA,
	CLASS=$class_archive,
	EXPVER=$expver_archive,
	LEVEL=OFF,
	PARAM=$ciparam,
	REPRES=GG,
	STEP=0,
	TIME=$CITIME,
	DATE=$CIDATE,
	TARGET="sfctempin"
	EOF
    #.
    #.DO THE MARS REQUEST.
    #=====================
    mars myrequest || {
    print - "\n\n\t\t NO SEA ICE FRACTION FROM MARS, NOTHING IS DONE \n\n"
    exit 8
    }
    mv sfctempin $ADIR/SSTDATA${begofrn}_${grid}
else
  ln -s $ADIR/SSTDATA${begofrn}_${grid} sfctempin
fi

if [[ "$COLDSTART" = "yes" ]] ; then
  #.
  #. Run preset
  #.
  # extract preset and load preset
  ar xv $LIBS/$WAMLIB.a preset.o 


  frt $rp -Wl,-dy -lsocket -Of -Wv,-Of -Ps -X9  -o preset preset.o \
    -lwam ${mpllib} ${emoslib} ${eclib} ${fdblib} \
    ${mpilib} ${cvplib} ${mp2lib} || {
  #.
  #.
  #. Abort on error.
  #.
    print - "\n\n\t\tLOAD FAILED\n\n"
    print - "\n\n\t\tPROGRAM WILL TERMINATE\n\n"
    exit 11
  }
  #.
  #.====================================================================
  #. RUN THE PROGRAM AND FORMAT THE OUTPUT.
  #.====================================================================
  #.
  #. THE DESCRIPTION OF THE INPUT NAMELIST IS DESCRIBED IN PROGRAM PRESET
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

else
  #.
  #.
  #. GET DATA FROM A PREVIOUS RUN
  #. ############################
  #.
  #. it it only set up to run with grib data which are available with
  #  a mars request (archived data or data still on ther FDB).
  #. If one wants to use input/output of binary
  #. restart files (the LAW and BLS files) then the user
  #. will have to properly administrate those files from the previous
  #. run.
  #.
  # Get the WAVE spectra 
  # ====================
  #.
  SPECDATE=`substring $begofrn 1 8`
  SPECTIME=`substring $begofrn 9 10`
  sfilename=SGS${begofrn}00

  if [[ !  -r $WDIR/${sfilename} ]] ; then
    # in case file is not there do a mars request
    # GET WAVE SPECTRA 
    #=================
    SPECDATE=`substring $begofrn 1 8`
    SPECTIME=`substring $begofrn 9 10`
 
    if [[ $grid -eq 050 ]] ; then
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
    elif [[ $grid -eq 100 ]] ; then
      echo "this option will not work if the archived fields"
      echo "not on that resolution "
      exit
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
  #.
  # Get drag coefficient with waves for wam model run 
  # ===============================
  #.
  if [[ !  -r $ADIR/CD${begofrn}_${grid} ]] ; then
    # in case file is not there do a mars request
    # GET DRAG COEFFICIENT WITH WAVE FROM MARS
    #=========================================
    CDDATE=`substring $begofrn 1 8`
    CDTIME=`substring $begofrn 9 10`
    if [[ $grid -eq 050 ]] ; then
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
    elif [[ $grid -eq 100 ]] ; then
      echo "this option will not work if the archived fields"
      echo "not on that resolution "
      exit
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

  #.
  # Get wind speed with waves for wam model run 
  # ===============================
  #.
  if [[ !  -r $ADIR/U10${begofrn}_${grid} ]] ; then
    # in case file is not there do a mars request
    # GET U10 WITH WAVE FROM MARS
    #=============================
    CDDATE=`substring $begofrn 1 8`
    CDTIME=`substring $begofrn 9 10`
    if [[ $grid -eq 050 ]] ; then
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
    elif [[ $grid -eq 100 ]] ; then
      echo "this option will not work if the archived fields"
      echo "not on that resolution "
      exit
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
    print - "\n\n\t\t U10 NO 245 FROM MARS, NOTHING IS DONE \n\n"
    exit 8
    }
    mv uwavein $ADIR/U10${begofrn}_${grid}
  fi

  rm -f uwavein
  ln -s $ADIR/U10${begofrn}_${grid} uwavein

# end test on COLDSTART
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
