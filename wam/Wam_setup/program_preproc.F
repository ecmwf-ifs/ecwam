#ifdef JOB
#.
#. PREPROC TO GENERATE UNFORMATTED INPUT FILE
#. FOR USE BY WAMODEL
#.
#QSUB -q vpp700.normal
#QSUB -eo
#QSUB -r preproc 
#QSUB -lP 0
#QSUB -lT 300
#QSUB -lV 10mb
#QSUB -lM 60mb
#QSUB -mb -mi
#QSUB
#.
######################################################################
# MAKE USE OF MEMORY RESIDENT FILE SYSTEM for TMPDIR
TMPDIR=$MRFSDIR
######################################################################
#.
PATH=$PATH:.
set +v
#.
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
cat > fort.1 << EOF
#if region=='s'
#include "/cc/rd/wam/Wam_setup/Topographic_data/Swamp_test/topocat_180"
#elif region=='m' && resolution==50
#include "/cc/rd/wam/Wam_setup/Topographic_data/Medite/topomed"
#elif region=='m' && resolution==25
#include "/cc/rd/wam/Wam_setup/Topographic_data/Global/topoglo_25"
#else
#include "/cc/rd/wam/Wam_setup/Topographic_data/Global/topoglo_50"
#endif
EOF
#.  ecp ec:$DHSPATH/topoglo2 fort.1
#.
#.  CURRENT DATA if refraction option in input is 2.
#.
#.ecp ec:$
#.
#.  UNFORMATED NEST FILE from coarse grid, if this is a fine grid run.
#.
#.ecp ec:$DHSPATH/BOUGLOU_${grid} fort.3
#.
#.======================================================================
#.  Copy the source code to the compiler input.
#.======================================================================
#.
cat >preproc.f <<\EOFPRO
#endif
#include "preproc.F"
#include "check.F"
#include "file_transfer.F"
#include "findb.F"
#include "jafu.F"
#include "locint.F"
#include "mblock.F"
#include "mbounc.F"
#include "mbounf.F"
#include "mboxb.F"
#include "mcout.F"
#include "mfredir.F"
#include "mgrid.F"
#include "mintf.F"
#include "mtabs.F"
#include "mubuf.F"
#include "nlweigt.F"
#include "outcom.F"
#include "outubuf.F"
#include "packi.F"
#include "packr.F"
#include "readcur.F"
#include "stress.F"
#include "tauhf.F"
#include "topoar.F"
#include "uiprep.F"
#include "wamcur.F"
#ifdef JOB
EOFPRO
#.
#.===================================================================
#. COMPILE AND LOAD
#.===================================================================
#.
frt $rp -Of -Wv,-Of -Posx -Z preproc.lis -X9 -c -o preproc.o preproc.f
#.
if [[ ! -a $LIBCR/lib$WAMCRLIB.a ]] ; then
  ecp  ec:$DHSLPATH/lib$WAMCRLIB.a $LIBCR/lib$WAMCRLIB.a || {
      ls -lsa
      print - "\n\n\t\tWAMCRLIB file was not found\n\n"
      print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
      exit 2
  }
fi
#.
frt $rp -Wl,-dy -lsocket -Of -Wv,-Of -Ps -X9  preproc.o -l$WAMCRLIB \
${EMOSLIB} ${ECLIB} ${FDBLIB} ${MPELIB} ${MPLIB} || {
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
#include "input_preproc_rmedit.F"
#endif
#else
#if resolution == 150
#include "input_preproc_rglobal150.F"
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
cp PARWAM $ADIR/PARWAM_${grid}
#.
#.  GRID FILE.
#.
cp fort.7 $ADIR/gridglou_${grid}
#.
#.  COMMON UBUF.
#.
cp fort.8 $ADIR/ubufglou_${grid}
#.
#.  COARSE GRID INFORMATION, if this is a coarse grid run.
#.
#.cp fort.9 $ADIR/boucglou${grid}
#.
#.  FINE GRID INFORMATION, if this is a fine grid run.
#.
#.cp fort.10 $ADIR/boufglou$_grid}
#.
#.
#.===================================================================
#.SAVE BACK-UP FILES IN ECFS 
#.===================================================================
#.
#.  Namelist-file parameter.
#.
####ecp -o $ADIR/PARWAM_${grid} ec:$DHSPATH/PARWAM_${grid}
#.
#.
#.  GRID FILE.
#.
#####ecp -o $ADIR/gridglou_${grid} ec:$DHSPATH/gridglou_${grid}
#.
#.  COMMON UBUF.
#.
#####ecp -o $ADIR/ubufglou_${grid} ec:$DHSPATH/ubufglou_${grid}
#.
#.  COARSE GRID INFORMATION, if this is a coarse grid run.
#.
#.ecp -o $ADIR/boucglou_${grid} ec:$DHSPATH/BOUCGLOU_${grid}
#.
#.  FINE GRID INFORMATION, if this is a fine grid run.
#.
#.ecp -o $ADIR/boufglou_${grid} ec:$ECPAT1/BOUFGLOU_${grid}
#.
#.===================================================================
#. GRID FILES HAVE BEEN CREATED AND SAVED.
#. END OF JOB PREPROC.
#.===================================================================
#.
ls -l
#.
#endif
