#
#  THIS JOB COMPILES, LOADS AND SAVES THE WAMODEL ABSOLUTE ON THE C90
#  AND ECFILE
#
#  ANNE GUILLAUME      ECMWF    NOVEMBER 1992
#
# QSUB -r WAM_TYE
# QSUB -lT 200
# QSUB -q high
# QSUB -eo
# QSUB

set +v
set -a
PATH=$PATH:/tmp/emos_sms/bin:/ec/emos_sms/bin:.
cd $TMPDIR

#======================================================================
#   Define variables
#======================================================================

EXPER=glo3d  
USEID=dax

#======================================================================

EXPPATH=/$USEID/$EXPER
YMPPATH=$TEMP/$EXPER

  ecfile -p /$USEID/$EXPER add
  ecfile -p /$USEID/$EXPER -v A/UV1/r/-/s modify
  mkdir $TEMP/$EXPER
  chmod go+rx $TEMP/$EXPER

#=======================================================================
#  Copy the WAM MODEL to the compiler input.
#=======================================================================

cat >wavemdl.f <<\EOFPRO
#include "readice.h"
#include "readt.h"
#include "inmarsb.h"
#include "micep.h"
#include "findb.h"
#include "msort.h"
#include "setice.h"
#include "findpic.h"
#include "regroup.h"
#include "trhou.h"
#include "vintpar.h"
#include "vtott.h"
#include "vagdirt.h"
#include "wtrain1.h"
#include "wtrain2.h"
#include "wtreorg.h"
#include "wtreosp.h"
#include "wtrain.h"
#include "wtraout.h"
#include "chief.h"
#include "airsea.h"
#include "bouinpt.h"
#include "dotdc.h"
#ifdef 
#include "ersfile_oper.h"
*elseif (ers=ersoff)
#include "ersfile.h"
#endif
#include "femean.h"
#include "fillbl.h"
#include "getwnd.h"
#include "gradi.h"
#ifdef 
#include "gsfile_unicos.h"
*else
#include "gsfile_oper_agswell.h"
#endif
#include "headbc.h"
#include "implsch.h"
#include "initmdl.h"
#include "intpol.h"
#include "locint.h"
#include "makegrid.h"
#include "notim.h"
#ifdef 
#include "openfil_unicos.h"
#endif
#include "outbc.h"
#include "outbs.h"
#include "outers.h"
#include "outgrid.h"
#include "outint.h"
#include "gribpac.h"
#include "fld2fdb.h"
#include "outspp.h"
#include "peakfr.h"
#include "prewind.h"
#include "propags.h"
#include "propdot.h"
#include "readbou.h"
#include "readpre.h"
#ifdef 
#include "readwnd_rswamp2.h"
*else
#ifdef 
#include "readwnd_seasat.h"
*else
#include "readwnd_mars.h"
#endif
#endif
#include "sbottom.h"
#include "sdissip.h"
#include "semean.h"
#include "sepwisw.h"
#include "sinput.h"
#include "snonlin.h"
#include "splitbl.h"
#include "sthq.h"
#include "stresso.h"
#include "timin.h"
CBH PACK: THIS ROUTINE TO GET USERINPUT IS NEW 2.12.92
#include "userin_pack.h"
#include "wamodel_agswell.h"
#include "wamwnd.h"
#ifdef 
#include "wavemdl_massi.h"
#include "analyse.h"
#include "f4spec.h"
#include "fdur.h"
#include "fustar.h"
#include "fwsea.h"
#include "grdata.h"
#include "makeblo.h"
#include "oifield.h"
#ifdef 
#include "readsat_geosat.h"
*elseif (v=seasat)
#include "readsat_seasat.h"
*elseif (v=ers1)
#include "readsat_ers1.h"
#endif
#include "setmat.h"
#include "spr.h"
#include "update.h"
#include "upspec.h"
#include "wamassi.h"
#include "wsmfen.h"
*else
#include "wavemdl_mwave.h"
#endif
EOFPRO

#=======================================================================
#  Compile and load the WAM MODEL.
#=======================================================================

NCPUS=14
export NCPUS
  cft77 -efszx -i 64 -Rabc -V wavemdl.f
# cft77   -efzs -V wavemdl.f
  if [ $? != 0 ] ; then
    echo "\n\n\t\t WAVEMODEL NOT COMPILED !! \n\n"
  else
    segldr wavemdl.o                 \
          -l /ec/stg/wamcrlib2       \
             $ECLIB                  \
             $EMOSLIB                \
             /usr/lib/libsort.a      \
             /tmp/emos_sms/lib/fdb.a \
             /lib/libprof.a          \
             /lib/libperf.a          \
           -g -V -o $YMPPATH/wama
    if [ $? != 0 ] ; then
      echo "\n\n\t\t WAVEMODEL NOT LOADED   !! \n\n"
    fi
    ecfile -p $EXPPATH/wama.o -hO -uD save wavemdl.o
  fi
  ls -l
  asa wavemdl.l

#=======================================================================
# END OF JOB.
#=======================================================================

