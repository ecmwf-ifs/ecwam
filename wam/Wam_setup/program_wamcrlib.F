#ifdef JOB
#.
#. This version of WAM was extracted from the ECMWF software data base
#. from branch wam_${VERSION}
#.
#.
#. THIS IS USED TO CREATE THE WAVE MODEL LIBRARY (libwam.a)
#. using the VPP crosscompiler on workstation
#. Note that you will need to have available (i.e already compiled)
#. the modules from the MPL library (see ifsaux). In this job,
#. we assume that we can use those that were already created for
#. the default cycle (see ifs_cycle)
#.
#QSUB -q normal
#QSUB -eo
#QSUB -r WAMLIB
#QSUB -lt 1000
#QSUB -lV 20mb
#QSUB -lM 64mb
#QSUB -s  /bin/ksh
#QSUB -me
#QSUB
#.
PATH=$PATH:.
set -v
set -e
#.
#include "wkstation_vpp_share.h"
#.
mkdir $WKLIBS || true
cd $WKLIBS
mkdir tmp$$
cd tmp$$
#.
cat >WAMLIB.src <<\EOFPRO
#endif
#include "/cc/rd/wam/module/yowaltas.F"
#include "/cc/rd/wam/module/yowcard.F"
#include "/cc/rd/wam/module/yowcbou.F"
#include "/cc/rd/wam/module/yowcinp.F"
#include "/cc/rd/wam/module/yowcoer.F"
#include "/cc/rd/wam/module/yowcoup.F"
#include "/cc/rd/wam/module/yowcout.F"
#include "/cc/rd/wam/module/yowcpbo.F"
#include "/cc/rd/wam/module/yowcurg.F"
#include "/cc/rd/wam/module/yowcurr.F"
#include "/cc/rd/wam/module/yowdes.F"
#include "/cc/rd/wam/module/yowfbou.F"
#include "/cc/rd/wam/module/yowfpbo.F"
#include "/cc/rd/wam/module/yowfred.F"
#include "/cc/rd/wam/module/yowgap.F"
#include "/cc/rd/wam/module/yowgribhd.F"
#include "/cc/rd/wam/module/yowgrid.F"
#include "/cc/rd/wam/module/yowice.F"
#include "/cc/rd/wam/module/yowindn.F"
#include "/cc/rd/wam/module/yowintp.F"
#include "/cc/rd/wam/module/yowints.F"
#include "/cc/rd/wam/module/yowjons.F"
#include "/cc/rd/wam/module/yowmap.F"
#include "/cc/rd/wam/module/yowmean.F"
#include "/cc/rd/wam/module/yowmespas.F"
#include "/cc/rd/wam/module/yowmltsk.F"
#include "/cc/rd/wam/module/yowmpp.F"
#include "/cc/rd/wam/module/yowover.F"
#include "/cc/rd/wam/module/yowparam.F"
#include "/cc/rd/wam/module/yowpcons.F"
#include "/cc/rd/wam/module/yowprproc.F"
#include "/cc/rd/wam/module/yowrefd.F"
#include "/cc/rd/wam/module/yowsaras.F"
#include "/cc/rd/wam/module/yowshal.F"
#include "/cc/rd/wam/module/yowspec.F"
#include "/cc/rd/wam/module/yowstat.F"
#include "/cc/rd/wam/module/yowswel.F"
#include "/cc/rd/wam/module/yowtabl.F"
#include "/cc/rd/wam/module/yowtemp.F"
#include "/cc/rd/wam/module/yowtest.F"
#include "/cc/rd/wam/module/yowtext.F"
#include "/cc/rd/wam/module/yowubuf.F"
#include "/cc/rd/wam/module/yowunit.F"
#include "/cc/rd/wam/module/yowwami.F"
#include "/cc/rd/wam/module/yowwind.F"
#include "/cc/rd/wam/module/yowwndg.F"
#include "abort1.F"
#include "adjust.F"
#include "airsea.F"
#include "aki.F"
#include "analyse.F"
#include "altas.F"
#include "avoidpeaks.F"
#include "bouinpt.F"
#include "buildstress.F"
#include "check.F"
#include "chief.F"
#include "class_wgrib.F"
#include "closend.F"
#include "combinepeaks.F"
#include "command.F"
#include "confile.F"
#include "cormars.F"
#include "crossas.F"
#include "crossw.F"
#include "cpfdb.F"
#include "crewfn.F"
#include "dev.F"
#include "difdate.F"
#include "docombine.F"
#include "dotdc.F"
#include "dummy_alt.F"
#include "dummy_cray.F"
#include "ersfile.F"
#include "expand_string.F"
#include "f4spec.F"
#include "fdur.F"
#include "femean.F"
#include "file_transfer.F"
#include "fillbl.F"
#include "fillgaps.F"
#include "findb.F"
#include "fld2fdb.F"
#include "fldinter.F"
#include "fustar.F"
#include "fwsea.F"
#include "gapinterpol.F"
#include "getinptb.F"
#include "getrest.F"
#include "getspec.F"
#include "getstress.F"
#include "getwnd.F"
#include "getwspec.F"
#include "gradi.F"
#include "grb2wgrd.F"
#include "grdata.F"
#include "grfield.F"
#include "gribpac.F"
#include "grstname.F"
#include "gsfile_new.F"
#include "halfspread.F"
#include "headbc.F"
#include "implsch.F"
#include "incdate.F"
#include "incint.F"
#include "initialint.F"
#include "initmdl.F"
#include "iniwcst.F"
#include "inmars.F"
#include "inmarss.F"
#include "inmarssi.F"
#include "intgrs.F"
#include "intmars.F"
#include "intpol.F"
#include "intspec.F"
#include "intwaminput.F"
#include "ispoff.F"
#include "jafu.F"
#include "jonswap.F"
#include "kgribsize.F"
#include "lasttwobins.F"
#include "llalloc_grdfld.F"
#include "loc2glo.F"
#include "locint.F"
#include "lowenergy.F"
#include "makeblo.F"
#include "makeblos.F"
#include "makeframes.F"
#include "makegrid.F"
#include "mblock.F"
#include "mbounc.F"
#include "mbounf.F"
#include "mboxb.F"
#include "mcout.F"
#include "means.F"
#include "meansqs.F"
#include "mergesarcor.F"
#include "mfredir.F"
#include "mgrid.F"
#include "micep.F"
#include "mintf.F"
#include "mnintw.F"
#include "mpabort.F"
#include "mpbcastgrdfld.F"
#include "mpbcastgrid.F"
#include "mpbcastintfld.F"
#include "mpbcastsarin.F"
#include "mpbcastscfld.F"
#include "mpclose_unit.F"
#include "mpcrtbl.F"
#include "mpdecomp.F"
#include "mpdistribfl.F"
#include "mpdistribintfld.F"
#include "mpdistribscfld.F"
#include "mpexchng.F"
#include "mpexchngsarin.F"
#include "mpgatherersfile.F"
#include "mpgatherfl.F"
#include "mpgathergrdfld.F"
#include "mpgatheroifld.F"
#include "mpgatherscfld.F"
#include "mpgatherspp.F"
#include "mpuserin.F"
#include "mstart.F"
#include "mt_allot.F"
#include "mtabs.F"
#include "mubuf.F"
#include "myindex.F"
#include "mwp1.F"
#include "mwp2.F"
#include "newwind.F"
#include "nlweigt.F"
#include "notim.F"
#include "oifield.F"
#include "onestepapart.F"
#include "optint.F"
#include "outbc.F"
#include "outbeta.F"
#include "outbs.F"
#include "outcom.F"
#include "outers.F"
#include "outgrid.F"
#include "outint.F"
#include "outnam.F"
#include "outpp.F"
#include "outspp.F"
#include "outubuf.F"
#include "outwaminput.F"
#include "outwnorm.F"
#include "outxt.F"
#include "packi.F"
#include "packr.F"
#include "partitioning.F"
#include "peak.F"
#include "peakfr.F"
#include "peakfri.F"
#include "phys.F"
#include "preproc.F"
#include "preset.F"
#include "prewind.F"
#include "propags.F"
#include "propdot.F"
#include "prspp.F"
#include "prspps.F"
#include "readbou.F"
#include "readcur.F"
#include "readfl.F"
#include "readwgrib.F"
#include "readgrs.F"
#include "readice.F"
#include "readpre.F"
#include "readsat.F"
#include "readsarspec.F"
#include "readsta.F"
#include "readstress.F"
#include "readt.F"
#include "readwind.F"
#include "rearrngsar.F"
#include "resize_gap_array.F"
#include "resize_getwspec.F"
#include "resize_sarinput.F"
#include "rfl4wam.F"
#include "rotspec.F"
#include "saras.F"
#include "savrest.F"
#include "savspec.F"
#include "savstress.F"
#include "sbottom.F"
#include "sdissip.F"
#include "scosfl.F"
#include "secondhh.F"
#include "semean.F"
#include "sepwisw.F"
#include "setice.F"
#include "setwgribhd.F"
#include "setmarstype.F"
#include "setmat.F"
#include "sinput.F"
#include "skewness.F"
#include "smoothsarspec.F"
#include "snonlin.F"
#include "spec2fdb.F"
#include "spectra.F"
#include "splitbl.F"
#include "splitgrs.F"
#include "spr.F"
#include "stats.F"
#include "statsc.F"
#include "statse.F"
#include "statsi.F"
#include "statsp.F"
#include "sthq.F"
#include "stress.F"
#include "stresso.F"
#include "strspec.F"
#include "sumenergy.F"
#include "swellsep.F"
#include "tauhf.F"
#include "threshhold.F"
#include "timin.F"
#include "topoar.F"
#include "transmeans.F"
#include "transpart.F"
#include "tustreas.F"
#include "u10cor.F"
#include "uiprep.F"
#include "update.F"
#include "updatewd.F"
#include "upwspec.F"
#include "userin.F"
#include "vmin.F"
#include "vplus.F"
#include "wam_user_clock.F"
#include "wamassi.F"
#include "wamcur.F"
#include "wamnorm.F"
#include "wamodel.F"
#include "wamwnd.F"
#include "wavemdl.F"
#include "windcor.F"
#include "windsea.F"
#include "writefl.F"
#include "wposnam.F"
#include "writestress.F"
#include "writsta.F"
#include "wsmfen.F"
#include "wvalloc.F"
#include "wvdealloc.F"
#include "wdirspread.F"
#include "wstream_strg.F"
#include "../Alt/i_get_unit.F"
#include "../Alt/uraopc.F"
#include "../Sar/accuwa.F"
#include "../Sar/aftfft.F"
#include "../Sar/calibra.F"
#include "../Sar/chasar.F"
#include "../Sar/chpar.F"
#include "../Sar/compoutpar.F"
#include "../Sar/costp.F"
#include "../Sar/fftsub.F"
#include "../Sar/fimagep.F"
#include "../Sar/findbb.F"
#include "../Sar/fttok.F"
#include "../Sar/fullmap.F"
#include "../Sar/getsarpar.F"
#include "../Sar/iners1.F"
#include "../Sar/isoplot.F"
#include "../Sar/mpgathersarinv.F"
#include "../Sar/newspec.F"
#include "../Sar/prepfft.F"
#include "../Sar/rarmtf.F"
#include "../Sar/readwsp.F"
#include "../Sar/sarinv.F"
#include "../Sar/sarinvert.F"
#include "../Sar/sarpc.F"
#include "../Sar/sarpol.F"
#include "../Sar/select.F"
#include "../Sar/subs_check.F"
#include "../Sar/usage_sar.F"
#include "../Sar/uwaror.F"
#include "../Sar/vbmtf.F"
#include "../Sar/waven.F"
#include "../Sar/wrires.F"
#include "../Sar/wtail.F"

#ifdef JOB
EOFPRO
#.
######################################################################
#.create utility which will be used to split the source code into its
#.components
######################################################################
#.
cat > f90split.f90 <<EOF
#include "/cc/rd/wam/Wam_setup/f90split.f90"
EOF
#.
f90 -o f90split f90split.f90
#.
#.
./f90split < WAMLIB.src
rm -f f90split.f90 f90split
#.
#.
# COMPILE SOURCE CODE
#####################
#.
[ ! -d tmp ] && mkdir tmp
#.
#. compile WAM modules
#######################
use sgi_vpp
integer frtcnt
frtcnt=0
export frtcnt
print $frtcnt > frtcnt
export file
for file in $(ls -1 yow*.f) ; do
  (
  frtpx -CcdRR8 -Karraystack1 -AoQm -Oe,-e,-u,-p -Pstao -Wv,-m3 -X9 -c -v9 -M. -AD -Z ${file%.f}.l -c $file
  (( frtcnt += 1))
  print $frtcnt > frtcnt
  )
  mv $file ${file%.f}.fm
done
#.
#. Get MPL modules
###################
#. !!! we will need the module definition for the MPL routines
#. !!! we assumed that they were compiled by default. If you are using this
#. !!! outside ECMWF you need to get the mpl software (i.e. the ECMWF
#. !!! interface to the message passing software (MPI,...)).
#.
export file ifs_cycle
for file in $(ls -1 /cc/rd/module/frt/${ifs_cycle}/ifsaux/) ; do
   cp /cc/rd/module/frt/${ifs_cycle}/ifsaux/${file} .
done

#.
#.compile regular fortran file
##############################
i=0
export file
for file in $(ls -1 *.f) ; do
  frtpx -CcdRR8 -Karraystack1 -AoQm -Oe,-e,-u,-p -Pstao -Wv,-m3 -X9 -c -v9  -M. -AD -Z ${file%.f}.l -c -o ${file%.f}.o ${file}
done
#. View compile listing
cat *.l

print - "compile done"
#.
#. BUILD LIBRARY
#################
#.
if [[ $(ls -1 *.o | wc -l ) -gt 0 ]] ; then
  for file in $(ls -1 yow*.fm) ; do
    ar rv ${WKLIBS}/${WAMLIB}.a ${file%.fm}.mod
  done
  ar rv ${WKLIBS}/${WAMLIB}.a *.o 
fi
#.
if [[ -f ${WKLIBS}/${WAMLIB}.l.tar ]] ; then
  tar uf ${WKLIBS}/${WAMLIB}.l.tar *.l
else
  tar cf ${WKLIBS}/${WAMLIB}.l.tar *.l
fi

tar cf ${WKLIBS}/${WAMLIB}.f.tar *.f
tar uf ${WKLIBS}/${WAMLIB}.f.tar *.fm

chmod go+rx ${WKLIBS}/${WAMLIB}.*
#.
ls -la
#.
#.
#.
cd ..
\rm -rf tmp$$

#.
#.END OF JOB
#.
#endif
