#ifdef JOB
# @ shell = /usr/bin/ksh
# @ core_limit = 4096
# @ notification = error
# @ restart = no
# @ cpu_limit = 2200,200
# @ job_cpu_limit = 5800,3800
# @ class = any
# @ job_name = wam_preproc
# @ output = /scratch/rd/wab/wam_preproc.$(jobid).out
# @ error  = /scratch/rd/wab/wam_preproc.$(jobid).err
# @ account_no = ecrmwa
# @ job_type = parallel
# @ total_tasks = 1
# @ queue
#.
#. PREPROC GENERATES UNFORMATTED INPUT FILES (bathymetry, grid and tables)
#. FOR USE BY PRESET and WAMODEL
#. This example runs as a batch job on the ECMWF IBM servers.
#.
#.
PATH=$PATH:.
set +v
#.
#include "wkstation_vpp_share.h"

# necessary directories
[[ ! -d $ROOTWDIR  ]] && mkdir -p $ROOTWDIR
[[ ! -d $WDIR  ]] && mkdir -p $WDIR

#.
cd $WDIR
mkdir tmp$$
cd tmp$$

# check if excutable is available
#################################
if [[ ! -a  ${BINS}/preproc ]] ; then
      print - "\n\n\t\tpreproc was not found\n\n"
      print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
      exit 3
fi

#.
#.===================================================================
#.  INPUT FILES.
#.===================================================================
#.
#.  TOPOGRAPHIC DATA.
#.
#. Note : the topographic data were obtained from the ETOPO5 data set
#.        (Data announcement 88-MGG02, digital relief of the surface
#.         of the earth. National Geophysical data center, Boulder,
#.         Colorado 1988, http://www.ngdc.noaa.gov/mgg/global/etopos5.HTML)
#.         It was smoothed out to 0.5 or 0.25 degree and formatted to
#.         the appropriate format for preproc.

cat > fort.1 << EOF
#if region=='s'
#include "/ccvobs/wam/Wam_setup/Topographic_data/Swamp_test/topocat_999"
#elif region=='m' && resolution==50
#include "/ccvobs/wam/Wam_setup/Topographic_data/Medite/topomed"
#elif region=='m' && resolution==25
#include "/ccvobs/wam/Wam_setup/Topographic_data/Global/topoglo_25"
#else
#include "/ccvobs/wam/Wam_setup/Topographic_data/Global/topoglo_50"
#endif
EOF
#.
#.
#.===================================================================
#. RUN THE PROGRAM AND FORMAT THE OUTPUT.
#.===================================================================
#.
# input namelist for preproc (used to specify domain and resolutions
#                             as well as manual bathymetric corrections)
cat > fort.4 <<\EOD
#if region=='s'
#include "input_preproc_rswamp2.F"
#elif region=='m'
#if resolution == 25 
#include "input_preproc_rmedite25.F"
#else
#include "input_preproc_rmedite.F"
#endif
#else
#if resolution == 150
#include "input_preproc_rglobal150.F"
#elif resolution == 100
#include "input_preproc_rglobal100.F"
#elif resolution == 50
#include "input_preproc_rglobal50.F"
#elif resolution == 900
#include "input_preproc_rglobal900.F"
#else
#include "input_preproc_rglobal.F"
#endif
#endif
EOD

# run preproc:
${BINS}/preproc || { print - "PREPROC FAILED" ; ls -ltra ; exit ; }
#.
#.===================================================================
#. STORE TEMPORARILY THE UNFORMATED OUTPUT FILES.
#.===================================================================
#.
#. PARWAM FILE
#.
mv PARWAM $WDIR/PARWAM_${grid}
#.
#.  GRID FILE.
#.
mv fort.7 $WDIR/gridglou_${grid}
#.
#.  COMMON UBUF.
#.
mv fort.8 $WDIR/ubufglou_${grid}
#.
#.===================================================================
#. GRID FILES HAVE BEEN CREATED AND SAVED.
#. END OF JOB PREPROC.
#.===================================================================
#.
ls -l

cd ..
\rm -rf tmp$$

#endif
