#ifdef JOB
#.
#. PREPROC TO GENERATE UNFORMATTED INPUT FILES (bathymetry, grid and tables)
#. FOR USE BY PRESET and WAMODEL
#.
#QSUB -q vpp700.normal
#QSUB -eo
#QSUB -r preproc 
#QSUB -lP 0
#QSUB -lT 300
#QSUB -lV 60mb
#QSUB -lM 190mb
#QSUB -mb -mi
#QSUB
#.
PATH=$PATH:.
set +v
#.
#include "wkstation_vpp_share.h"
#include "path.h"
#.
#.
cd $TMPDIR
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
#include "/cc/rd/wam/Wam_setup/Topographic_data/Swamp_test/topocat_999"
#elif region=='m' && resolution==50
#include "/cc/rd/wam/Wam_setup/Topographic_data/Medite/topomed"
#elif region=='m' && resolution==25
#include "/cc/rd/wam/Wam_setup/Topographic_data/Global/topoglo_25"
#else
#include "/cc/rd/wam/Wam_setup/Topographic_data/Global/topoglo_50"
#endif
EOF
#.
#.========
#. LOAD
#.========
#.
#.
#. copy library from workstation

cp $WKLIBS/$WAMLIB.a $LIBS/$WAMLIB.a

if [[ ! -a $LIBS/$WAMLIB.a ]] ; then
    ls -lsa
    print - "\n\n\t\tWAMLIB file was not found\n\n"
    print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
    exit 2
fi
#.
#. extract preproc.o

ar xv $LIBS/$WAMLIB.a preproc.o

if [[ ! -a  preproc.o ]] ; then
      print - "\n\n\t\tpreproc.o file was not found\n\n"
      print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
      exit 3
fi
#.
frt $rp -Wl,-dy -lsocket -Of -Wv,-Of -Ps -X9  preproc.o -lwam \
${mpllib} ${emoslib} ${eclib} ${fdblib} ${mpilib} ${mp2lib} || {
   #.
   #. Check return code.
   #.
    print - "\n\n\t\tCOMPILE FAILED\n\n"
    print - "\n\n\t\tPROGRAM WILL TERMINATE\n\n"
    cat preproc.lis
    print - "\n\n\t\tCOMPILE FAILED\n\n"
    print - "\n\n\t\tPROGRAM WILL TERMINATE\n\n"
    exit 3
   }
#.
#.===================================================================
#. RUN THE PROGRAM AND FORMAT THE OUTPUT.
#.===================================================================
#.
a.out <<\EOD ||  { print - "\n\n\t\tPREPROC FAILED\n\n" ; ls -ltra ; exit ; }
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
#.
#.===================================================================
#. STORE TEMPORARILY THE UNFORMATED OUTPUT FILES.
#.===================================================================
#.
#. PARWAM FILE
#.
mv PARWAM $ADIR/PARWAM_${grid}
#.
#.  GRID FILE.
#.
mv fort.7 $ADIR/gridglou_${grid}
#.
#.  COMMON UBUF.
#.
mv fort.8 $ADIR/ubufglou_${grid}
#.
#.===================================================================
#. GRID FILES HAVE BEEN CREATED AND SAVED.
#. END OF JOB PREPROC.
#.===================================================================
#.
ls -l
#.
#endif
