#
#  THIS JOB EXECUTES THE SUBROUTINE VERSION OF THE WAMODEL,
#  DRIVEN BY WINDS REQUESTED FROM MARS.
#
#  LIANA ZAMBRESKY     GKSS/ECMWF    JUNE  1989
#  HEINZ GUNTHER       GKSS/ECMWF    JUNE  1990
#  PETER JANSSEN       KNMI          OCT   1990   (MIGRATION TO UNICOS)
#  HEINZ GUNTHER       GKSS/ECMWF    APRIL 1991   (CYCLE_3.5)
#
# QSUB -r WAMODEL
# QSUB -q high
# QSUB -lT 400
# QSUB -eo
# QSUB
#
PATH=$PATH:.
cd $TMPDIR
set +v
#
#======================================================================
#   Define pathes for ECFILE and YMP
#======================================================================
#
#ifdef 
ECPATH=/nab/cos/swamp2
YMPPATH=$TEMP/swamp2
*elseif (r=medite)
ECPATH=/nab/cos/medite
YMPPATH=$TEMP/medite
*else
ECPATH=/nab/cos/wamoper
YMPPATH=$TEMP/wamoper
#endif
#
#======================================================================
#   Fetch WAMINFO file.
#======================================================================
#
test -r $YMPPATH/waminfo
if [ $? != 0 ]
then
   ecfile -p $ECPATH/waminfo -eG get $YMPPATH/waminfo
   if [ $? != 0 ]
   then
      ls -lsa
      echo "\n\n\t\tWAMINFO file was not found\n\n"
      echo "\n\n\t\tPROGRAM IS TERMINATED\n\n"
      exit 2
   fi
fi
cp $YMPPATH/waminfo WAMINFO
#
#======================================================================
#  General grid information and model constants.
#======================================================================
#
test -r $YMPPATH/gridglou
if [ $? != 0 ]
then
   ecfile -p $ECPATH/gridglou -eG get $YMPPATH/gridglou
   if [ $? != 0 ]
   then
      ls -lsa
      echo "\n\n\t\tPREPROC file gridglou was not found\n\n"
      echo "\n\n\t\tPROGRAM IS TERMINATED\n\n"
      exit 2
   fi
fi
assign -a $YMPPATH/gridglou GRID
#ifdef 
#
*else
#
#======================================================================
#   Prepare MARS requests. Execute program GETINPT.
#======================================================================
#
test -r $YMPPATH/getinptb
if [ $? != 0 ]
then
   ecfile -p $ECPATH/getinptb -eG get $YMPPATH/getinptb
   if [ $? != 0 ]
   then
      ls -lsa
      echo "\n\n\t\tGETINPTB WAS NOT FOUND\n\n"
      echo "\n\n\t\tPROGRAM IS TERMINATED\n\n"
      exit 5
   fi
fi
$YMPPATH/getinptb > tempout
#
#  Check error code.
#
if [ $? != 0 ]
then
   asa < tempout
   ls -lsa
   echo "\n\n\t\tGETINPT FAILED\n\n"
   echo "\n\n\t\tPROGRAM IS TERMINATED\n\n"
   exit 6
else
   asa < tempout
   echo "\n\n\t\tSTART MARS\n\n"
fi
rm tempout
#
#======================================================================
#   MARS retrieval.
#======================================================================
#
mars fort.22
#
#  Check MARS error code.
#
if [ $? != 0 ]
then
   ls -lsa
   echo "\n\n\t\tMARS FAILED\n\n"
   echo "\n\n\t\tPROGRAM IS TERMINATED\n\n"
   exit 7
else
   echo "\n\n\t\tSTART WAVE MODEL \n\n"
fi
#
rm fort.22
#endif
#
#======================================================================
#   WAVE MODEL execution.
#======================================================================
#
#  Common UBUF.
#
test -r $YMPPATH/ubufglou
if [ $? != 0 ]
then
   ecfile -p $ECPATH/ubufglou -eG get $YMPPATH/ubufglou
   if [ $? != 0 ]
   then
      ls -lsa
      echo "\n\n\t\tPREPROC file ubufglou was not found\n\n"
      echo "\n\n\t\tPROGRAM IS TERMINATED\n\n"
      exit 8
   fi
fi
assign -a $YMPPATH/ubufglou UBUF
#
#  Coarse grid output points, if this is a coarse grid run.
#
#assign -a $YMPPATH/bouccatu fort.9
#
#  Fine grid input points, if this is a fine grid run.
#
#assign -a $YMPPAT1/boufcatu fort.10
#
#=======================================================================
#  Copy the WAM MODEL to the compiler input.
#=======================================================================
#
cat >WAVEMDL.f <<\EOFPRO
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
#include "gsfile_oper.h"
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
#include "userin.h"
#include "wamodel.h"
#include "wamwnd.h"
#include "wposnam.h"
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
#
#=======================================================================
#  Compile and load the WAM MODEL.
#=======================================================================
#
NCPUS=8
export NCPUS
cf77 -l /ec/nab/wamoper/wamcrlib2,$ECLIB,$EMOSLIB WAVEMDL.f
if [ $? != 0 ]
then
   echo "\n\n\t\tCompile or load failed"
   echo "\t\tPROGRAM WILL TERMINATE\n\n"
   exit 3
else
   echo "\n\n\t\tSTART WAVE MODEL EXECUTION\n\n"
fi
rm WAVEMDL.f
pshell a.out <<\EOD > tempout
#ifdef 
#include "input_wamodel_rswamp2.h"
*elseif (r=medite)
#include "input_wamodel_rmedite.h"
*else
#include "input_wamodel_rglobal.h"
#endif
EOD
#
#   Check WAVE model error code.
#
if [ $? != 0 ]
then
#
#   Error: Copy program output to stdout and terminate.
#
   asa < tempout
   ls -lsa
   echo "\n\n\t\tWAVE MODEL FAILED\n\n"
   exit 10
fi
asa < tempout
#
#=======================================================================
# END OF JOB.
#=======================================================================
#
