#ifdef JOB
#!/bin/ksh
#.
#. This version of WAM was extracted from the ECMWF software data base
#. from branch view_${ECUSER}_${VERSION} (see below).
#.
#.===========================================================================
#. THIS IS USED TO CREATE THE WAVE MODEL LIBRARY (libwam.a)
#. AND THE EXCECUTABLES PROVIDED YOU HAVE THE NECESSARY LIBRARIES.

#. see below for the details but in the current configuration you will need
#. the ECLIB and the EMOSLIB from ECMWF
#. and the message passing library MPI for your computer.

#.===========================================================================
#.
#. This job run interactively on the ECMWF IBM servers and could also run
#. (in batch mode) on the IBM supercomputer. See the config.* and options.*
#. files below for the compiler options (those are the flags used at ECMWF
#. Note that the WAM code is still in fortran fixed format, however the
#. message passing code comes from the ifsaux library which is in free format.  
#. Note: this file contains all the different routines that are needed
#. and the job will cat each one of them as a separate file. For better
#. efficiency you might want to put all routines online some place and only
#. use this job for compilation and build.
#.
#. General information on the wave model at ECMWF can obtained at
#. http://ifsdoc.ecmwf.int/html/rd/
#. Then pick the relevent cycle , then WAM, then Documentation
#.
set +ve
#.
#include "wkstation_vpp_share.h"
#.

# necessary directories
[[ ! -d $TMPROOT  ]] && mkdir $TMPROOT
[[ ! -d $LIBS  ]] && mkdir $LIBS
[[ ! -d $BINS  ]] && mkdir $BINS

cd $LIBS

[[ ! -d build  ]] && mkdir build

\rm -rf module
mkdir module

mkdir tmp$$
cd tmp$$

#################
# WAM model code
#################

# Main programs
################
cat >chief.F <<\EOFPRO
#include "chief.F"
EOFPRO
cat >preproc.F <<\EOFPRO
#include "preproc.F"
EOFPRO
cat >preset.F <<\EOFPRO
#include "preset.F"
EOFPRO

# Modules
##########
cat >yowaltas.F <<\EOFPRO
#include "yowaltas.F"
EOFPRO
cat >yowcard.F <<\EOFPRO
#include "yowcard.F"
EOFPRO
cat >yowcbou.F <<\EOFPRO
#include "yowcbou.F"
EOFPRO
cat >yowcinp.F <<\EOFPRO
#include "yowcinp.F"
EOFPRO
cat >yowcoer.F <<\EOFPRO
#include "yowcoer.F"
EOFPRO
cat >yowcoup.F <<\EOFPRO
#include "yowcoup.F"
EOFPRO
cat >yowcout.F <<\EOFPRO
#include "yowcout.F"
EOFPRO
cat >yowcpbo.F <<\EOFPRO
#include "yowcpbo.F"
EOFPRO
cat >yowcurg.F <<\EOFPRO
#include "yowcurg.F"
EOFPRO
cat >yowcurr.F <<\EOFPRO
#include "yowcurr.F"
EOFPRO
cat >yowdes.F <<\EOFPRO
#include "yowdes.F"
EOFPRO
cat >yowfbou.F <<\EOFPRO
#include "yowfbou.F"
EOFPRO
cat >yowfpbo.F <<\EOFPRO
#include "yowfpbo.F"
EOFPRO
cat >yowfred.F <<\EOFPRO
#include "yowfred.F"
EOFPRO
cat >yowgap.F <<\EOFPRO
#include "yowgap.F"
EOFPRO
cat >yowgribhd.F <<\EOFPRO
#include "yowgribhd.F"
EOFPRO
cat >yowgrid.F <<\EOFPRO
#include "yowgrid.F"
EOFPRO
cat >yowice.F <<\EOFPRO
#include "yowice.F"
EOFPRO
cat >yowindn.F <<\EOFPRO
#include "yowindn.F"
EOFPRO
cat >yowintp.F <<\EOFPRO
#include "yowintp.F"
EOFPRO
cat >yowints.F <<\EOFPRO
#include "yowints.F"
EOFPRO
cat >yowjons.F <<\EOFPRO
#include "yowjons.F"
EOFPRO
cat >yowmap.F <<\EOFPRO
#include "yowmap.F"
EOFPRO
cat >yowmean.F <<\EOFPRO
#include "yowmean.F"
EOFPRO
cat >yowmespas.F <<\EOFPRO
#include "yowmespas.F"
EOFPRO
cat >yowmltsk.F <<\EOFPRO
#include "yowmltsk.F"
EOFPRO
cat >yowmpp.F <<\EOFPRO
#include "yowmpp.F"
EOFPRO
cat >yowover.F <<\EOFPRO
#include "yowover.F"
EOFPRO
cat >yowparam.F <<\EOFPRO
#include "yowparam.F"
EOFPRO
cat >yowpcons.F <<\EOFPRO
#include "yowpcons.F"
EOFPRO
cat >yowprproc.F <<\EOFPRO
#include "yowprproc.F"
EOFPRO
cat >yowrefd.F <<\EOFPRO
#include "yowrefd.F"
EOFPRO
cat >yowsaras.F <<\EOFPRO
#include "yowsaras.F"
EOFPRO
cat >yowshal.F <<\EOFPRO
#include "yowshal.F"
EOFPRO
cat >yowspec.F <<\EOFPRO
#include "yowspec.F"
EOFPRO
cat >yowstat.F <<\EOFPRO
#include "yowstat.F"
EOFPRO
cat >yowswel.F <<\EOFPRO
#include "yowswel.F"
EOFPRO
cat >yowtabl.F <<\EOFPRO
#include "yowtabl.F"
EOFPRO
cat >yowtemp.F <<\EOFPRO
#include "yowtemp.F"
EOFPRO
cat >yowtest.F <<\EOFPRO
#include "yowtest.F"
EOFPRO
cat >yowtext.F <<\EOFPRO
#include "yowtext.F"
EOFPRO
cat >yowubuf.F <<\EOFPRO
#include "yowubuf.F"
EOFPRO
cat >yowunit.F <<\EOFPRO
#include "yowunit.F"
EOFPRO
cat >yowwami.F <<\EOFPRO
#include "yowwami.F"
EOFPRO
cat >yowwind.F <<\EOFPRO
#include "yowwind.F"
EOFPRO
cat >yowwndg.F <<\EOFPRO
#include "yowwndg.F"
EOFPRO

# Subroutines
#############

cat >abort1.F <<\EOFPRO
#include "abort1.F"
EOFPRO
cat >adjust.F <<\EOFPRO
#include "adjust.F"
EOFPRO
cat >airsea.F <<\EOFPRO
#include "airsea.F"
EOFPRO
cat >aki.F <<\EOFPRO
#include "aki.F"
EOFPRO
cat >analyse.F <<\EOFPRO
#include "analyse.F"
EOFPRO
cat >altas.F <<\EOFPRO
#include "altas.F"
EOFPRO
cat >avoidpeaks.F <<\EOFPRO
#include "avoidpeaks.F"
EOFPRO
cat >bouinpt.F <<\EOFPRO
#include "bouinpt.F"
EOFPRO
cat >buildstress.F <<\EOFPRO
#include "buildstress.F"
EOFPRO
cat >check.F <<\EOFPRO
#include "check.F"
EOFPRO
cat >class_wgrib.F <<\EOFPRO
#include "class_wgrib.F"
EOFPRO
cat >closend.F <<\EOFPRO
#include "closend.F"
EOFPRO
cat >combinepeaks.F <<\EOFPRO
#include "combinepeaks.F"
EOFPRO
cat >confile.F <<\EOFPRO
#include "confile.F"
EOFPRO
cat >crossas.F <<\EOFPRO
#include "crossas.F"
EOFPRO
cat >crossw.F <<\EOFPRO
#include "crossw.F"
EOFPRO
cat >crewfn.F <<\EOFPRO
#include "crewfn.F"
EOFPRO
cat >dev.F <<\EOFPRO
#include "dev.F"
EOFPRO
cat >difdate.F <<\EOFPRO
#include "difdate.F"
EOFPRO
cat >docombine.F <<\EOFPRO
#include "docombine.F"
EOFPRO
cat >dotdc.F <<\EOFPRO
#include "dotdc.F"
EOFPRO
cat >dummy_alt.F <<\EOFPRO
#include "dummy_alt.F"
EOFPRO
cat >dummy_cray.F <<\EOFPRO
#include "dummy_cray.F"
EOFPRO
cat >dummy_eclib.F <<\EOFPRO
#include "dummy_eclib.F"
EOFPRO
cat >dummy_fdb.F <<\EOFPRO
#include "dummy_fdb.F"
EOFPRO
cat >ersfile.F <<\EOFPRO
#include "ersfile.F"
EOFPRO
cat >expand_string.F <<\EOFPRO
#include "expand_string.F"
EOFPRO
cat >fdur.F <<\EOFPRO
#include "fdur.F"
EOFPRO
cat >femean.F <<\EOFPRO
#include "femean.F"
EOFPRO
cat >file_transfer.F <<\EOFPRO
#include "file_transfer.F"
EOFPRO
cat >fillbl.F <<\EOFPRO
#include "fillbl.F"
EOFPRO
cat >fillgaps.F <<\EOFPRO
#include "fillgaps.F"
EOFPRO
cat >findb.F <<\EOFPRO
#include "findb.F"
EOFPRO
cat >fld2fdb.F <<\EOFPRO
#include "fld2fdb.F"
EOFPRO
cat >fldinter.F <<\EOFPRO
#include "fldinter.F"
EOFPRO
cat >fustar.F <<\EOFPRO
#include "fustar.F"
EOFPRO
cat >fwsea.F <<\EOFPRO
#include "fwsea.F"
EOFPRO
cat >gapinterpol.F <<\EOFPRO
#include "gapinterpol.F"
EOFPRO
cat >getinptb.F <<\EOFPRO
#include "getinptb.F"
EOFPRO
cat >getrest.F <<\EOFPRO
#include "getrest.F"
EOFPRO
cat >getspec.F <<\EOFPRO
#include "getspec.F"
EOFPRO
cat >getstress.F <<\EOFPRO
#include "getstress.F"
EOFPRO
cat >getwnd.F <<\EOFPRO
#include "getwnd.F"
EOFPRO
cat >getwspec.F <<\EOFPRO
#include "getwspec.F"
EOFPRO
cat >gradi.F <<\EOFPRO
#include "gradi.F"
EOFPRO
cat >grb2wgrd.F <<\EOFPRO
#include "grb2wgrd.F"
EOFPRO
cat >grdata.F <<\EOFPRO
#include "grdata.F"
EOFPRO
cat >grfield.F <<\EOFPRO
#include "grfield.F"
EOFPRO
cat >gribpac.F <<\EOFPRO
#include "gribpac.F"
EOFPRO
cat >grstname.F <<\EOFPRO
#include "grstname.F"
EOFPRO
cat >gsfile_new.F <<\EOFPRO
#include "gsfile_new.F"
EOFPRO
cat >halfspread.F <<\EOFPRO
#include "halfspread.F"
EOFPRO
cat >headbc.F <<\EOFPRO
#include "headbc.F"
EOFPRO
cat >implsch.F <<\EOFPRO
#include "implsch.F"
EOFPRO
cat >incdate.F <<\EOFPRO
#include "incdate.F"
EOFPRO
cat >incint.F <<\EOFPRO
#include "incint.F"
EOFPRO
cat >initialint.F <<\EOFPRO
#include "initialint.F"
EOFPRO
cat >initmdl.F <<\EOFPRO
#include "initmdl.F"
EOFPRO
cat >iniwcst.F <<\EOFPRO
#include "iniwcst.F"
EOFPRO
cat >inmars.F <<\EOFPRO
#include "inmars.F"
EOFPRO
cat >inmarssi.F <<\EOFPRO
#include "inmarssi.F"
EOFPRO
cat >inmarss.F <<\EOFPRO
#include "inmarss.F"
EOFPRO
cat >intgrs.F <<\EOFPRO
#include "intgrs.F"
EOFPRO
cat >intmars.F <<\EOFPRO
#include "intmars.F"
EOFPRO
cat >intpol.F <<\EOFPRO
#include "intpol.F"
EOFPRO
cat >intspec.F <<\EOFPRO
#include "intspec.F"
EOFPRO
cat >intwaminput.F <<\EOFPRO
#include "intwaminput.F"
EOFPRO
cat >ispoff.F <<\EOFPRO
#include "ispoff.F"
EOFPRO
cat >jafu.F <<\EOFPRO
#include "jafu.F"
EOFPRO
cat >jonswap.F <<\EOFPRO
#include "jonswap.F"
EOFPRO
cat >kgribsize.F <<\EOFPRO
#include "kgribsize.F"
EOFPRO
cat >kurtosis.F <<\EOFPRO
#include "kurtosis.F"
EOFPRO
cat >lasttwobins.F <<\EOFPRO
#include "lasttwobins.F"
EOFPRO
cat >llalloc_grdfld.F <<\EOFPRO
#include "llalloc_grdfld.F"
EOFPRO
cat >loc2glo.F <<\EOFPRO
#include "loc2glo.F"
EOFPRO
cat >locint.F <<\EOFPRO
#include "locint.F"
EOFPRO
cat >lowenergy.F <<\EOFPRO
#include "lowenergy.F"
EOFPRO
cat >makeblo.F <<\EOFPRO
#include "makeblo.F"
EOFPRO
cat >makeblos.F <<\EOFPRO
#include "makeblos.F"
EOFPRO
cat >makeframes.F <<\EOFPRO
#include "makeframes.F"
EOFPRO
cat >makegrid.F <<\EOFPRO
#include "makegrid.F"
EOFPRO
cat >mblock.F <<\EOFPRO
#include "mblock.F"
EOFPRO
cat >mbounc.F <<\EOFPRO
#include "mbounc.F"
EOFPRO
cat >mbounf.F <<\EOFPRO
#include "mbounf.F"
EOFPRO
cat >mboxb.F <<\EOFPRO
#include "mboxb.F"
EOFPRO
cat >mcout.F <<\EOFPRO
#include "mcout.F"
EOFPRO
cat >means.F <<\EOFPRO
#include "means.F"
EOFPRO
cat >meansqs.F <<\EOFPRO
#include "meansqs.F"
EOFPRO
cat >mergesarcor.F <<\EOFPRO
#include "mergesarcor.F"
EOFPRO
cat >mfredir.F <<\EOFPRO
#include "mfredir.F"
EOFPRO
cat >mgrid.F <<\EOFPRO
#include "mgrid.F"
EOFPRO
cat >micep.F <<\EOFPRO
#include "micep.F"
EOFPRO
cat >mintf.F <<\EOFPRO
#include "mintf.F"
EOFPRO
cat >mnintw.F <<\EOFPRO
#include "mnintw.F"
EOFPRO
cat >mpabort.F <<\EOFPRO
#include "mpabort.F"
EOFPRO
cat >mpbcastgrdfld.F <<\EOFPRO
#include "mpbcastgrdfld.F"
EOFPRO
cat >mpbcastgrid.F <<\EOFPRO
#include "mpbcastgrid.F"
EOFPRO
cat >mpbcastintfld.F <<\EOFPRO
#include "mpbcastintfld.F"
EOFPRO
cat >mpbcastsarin.F <<\EOFPRO
#include "mpbcastsarin.F"
EOFPRO
cat >mpbcastscfld.F <<\EOFPRO
#include "mpbcastscfld.F"
EOFPRO
cat >mpclose_unit.F <<\EOFPRO
#include "mpclose_unit.F"
EOFPRO
cat >mpcrtbl.F <<\EOFPRO
#include "mpcrtbl.F"
EOFPRO
cat >mpdecomp.F <<\EOFPRO
#include "mpdecomp.F"
EOFPRO
cat >mpdistribfl.F <<\EOFPRO
#include "mpdistribfl.F"
EOFPRO
cat >mpdistribintfld.F <<\EOFPRO
#include "mpdistribintfld.F"
EOFPRO
cat >mpdistribscfld.F <<\EOFPRO
#include "mpdistribscfld.F"
EOFPRO
cat >mpexchng.F <<\EOFPRO
#include "mpexchng.F"
EOFPRO
cat >mpexchngsarin.F <<\EOFPRO
#include "mpexchngsarin.F"
EOFPRO
cat >mpgatherersfile.F <<\EOFPRO
#include "mpgatherersfile.F"
EOFPRO
cat >mpgatherfl.F <<\EOFPRO
#include "mpgatherfl.F"
EOFPRO
cat >mpgathergrdfld.F <<\EOFPRO
#include "mpgathergrdfld.F"
EOFPRO
cat >mpgatheroifld.F <<\EOFPRO
#include "mpgatheroifld.F"
EOFPRO
cat >mpgatherscfld.F <<\EOFPRO
#include "mpgatherscfld.F"
EOFPRO
cat >mpgatherspp.F <<\EOFPRO
#include "mpgatherspp.F"
EOFPRO
cat >mpuserin.F <<\EOFPRO
#include "mpuserin.F"
EOFPRO
cat >mppewithindist.F <<\EOFPRO
#include "mppewithindist.F"
EOFPRO
cat >mstart.F <<\EOFPRO
#include "mstart.F"
EOFPRO
cat >mtabs.F <<\EOFPRO
#include "mtabs.F"
EOFPRO
cat >mubuf.F <<\EOFPRO
#include "mubuf.F"
EOFPRO
cat >mwp1.F <<\EOFPRO
#include "mwp1.F"
EOFPRO
cat >mwp2.F <<\EOFPRO
#include "mwp2.F"
EOFPRO
cat >newwind.F <<\EOFPRO
#include "newwind.F"
EOFPRO
cat >nlweigt.F <<\EOFPRO
#include "nlweigt.F"
EOFPRO
cat >notim.F <<\EOFPRO
#include "notim.F"
EOFPRO
cat >oifield.F <<\EOFPRO
#include "oifield.F"
EOFPRO
cat >onestepapart.F <<\EOFPRO
#include "onestepapart.F"
EOFPRO
cat >optint.F <<\EOFPRO
#include "optint.F"
EOFPRO
cat >outbc.F <<\EOFPRO
#include "outbc.F"
EOFPRO
cat >outbeta.F <<\EOFPRO
#include "outbeta.F"
EOFPRO
cat >outbs.F <<\EOFPRO
#include "outbs.F"
EOFPRO
cat >outcom.F <<\EOFPRO
#include "outcom.F"
EOFPRO
cat >outers.F <<\EOFPRO
#include "outers.F"
EOFPRO
cat >outgrid.F <<\EOFPRO
#include "outgrid.F"
EOFPRO
cat >outint.F <<\EOFPRO
#include "outint.F"
EOFPRO
cat >outnam.F <<\EOFPRO
#include "outnam.F"
EOFPRO
cat >outpp.F <<\EOFPRO
#include "outpp.F"
EOFPRO
cat >outspp.F <<\EOFPRO
#include "outspp.F"
EOFPRO
cat >outubuf.F <<\EOFPRO
#include "outubuf.F"
EOFPRO
cat >outwaminput.F <<\EOFPRO
#include "outwaminput.F"
EOFPRO
cat >outwnorm.F <<\EOFPRO
#include "outwnorm.F"
EOFPRO
cat >outxt.F <<\EOFPRO
#include "outxt.F"
EOFPRO
cat >packi.F <<\EOFPRO
#include "packi.F"
EOFPRO
cat >packr.F <<\EOFPRO
#include "packr.F"
EOFPRO
cat >partitioning.F <<\EOFPRO
#include "partitioning.F"
EOFPRO
cat >peak.F <<\EOFPRO
#include "peak.F"
EOFPRO
cat >peakfr.F <<\EOFPRO
#include "peakfr.F"
EOFPRO
cat >peakfri.F <<\EOFPRO
#include "peakfri.F"
EOFPRO
cat >phys.F <<\EOFPRO
#include "phys.F"
EOFPRO
cat >prewind.F <<\EOFPRO
#include "prewind.F"
EOFPRO
cat >propags.F <<\EOFPRO
#include "propags.F"
EOFPRO
cat >propdot.F <<\EOFPRO
#include "propdot.F"
EOFPRO
cat > prspp.F<<\EOFPRO
#include "prspp.F"
EOFPRO
cat >prspps.F <<\EOFPRO
#include "prspps.F"
EOFPRO
cat >readbou.F <<\EOFPRO
#include "readbou.F"
EOFPRO
cat >readcur.F <<\EOFPRO
#include "readcur.F"
EOFPRO
cat >readfl.F <<\EOFPRO
#include "readfl.F"
EOFPRO
cat >readwgrib.F <<\EOFPRO
#include "readwgrib.F"
EOFPRO
cat >readgrs.F <<\EOFPRO
#include "readgrs.F"
EOFPRO
cat >readice.F <<\EOFPRO
#include "readice.F"
EOFPRO
cat >readpre.F <<\EOFPRO
#include "readpre.F"
EOFPRO
cat >readsat.F <<\EOFPRO
#include "readsat.F"
EOFPRO
cat >readsarspec.F <<\EOFPRO
#include "readsarspec.F"
EOFPRO
cat >readsta.F <<\EOFPRO
#include "readsta.F"
EOFPRO
cat >readstress.F <<\EOFPRO
#include "readstress.F"
EOFPRO
cat >readt.F <<\EOFPRO
#include "readt.F"
EOFPRO
cat >readwind.F <<\EOFPRO
#include "readwind.F"
EOFPRO
cat >rearrngsar.F <<\EOFPRO
#include "rearrngsar.F"
EOFPRO
cat >resize_gap_array.F <<\EOFPRO
#include "resize_gap_array.F"
EOFPRO
cat >resize_getwspec.F <<\EOFPRO
#include "resize_getwspec.F"
EOFPRO
cat >resize_sarinput.F <<\EOFPRO
#include "resize_sarinput.F"
EOFPRO
cat >rotspec.F <<\EOFPRO
#include "rotspec.F"
EOFPRO
cat >saras.F <<\EOFPRO
#include "saras.F"
EOFPRO
cat >savrest.F <<\EOFPRO
#include "savrest.F"
EOFPRO
cat >savspec.F <<\EOFPRO
#include "savspec.F"
EOFPRO
cat >savstress.F <<\EOFPRO
#include "savstress.F"
EOFPRO
cat >sbottom.F <<\EOFPRO
#include "sbottom.F"
EOFPRO
cat >sdissip.F <<\EOFPRO
#include "sdissip.F"
EOFPRO
cat >scosfl.F <<\EOFPRO
#include "scosfl.F"
EOFPRO
cat >secondhh.F <<\EOFPRO
#include "secondhh.F"
EOFPRO
cat >semean.F <<\EOFPRO
#include "semean.F"
EOFPRO
cat >sepwisw.F <<\EOFPRO
#include "sepwisw.F"
EOFPRO
cat >setice.F <<\EOFPRO
#include "setice.F"
EOFPRO
cat >setwgribhd.F <<\EOFPRO
#include "setwgribhd.F"
EOFPRO
cat >setmarstype.F <<\EOFPRO
#include "setmarstype.F"
EOFPRO
cat >sinput.F <<\EOFPRO
#include "sinput.F"
EOFPRO
cat >skewness.F <<\EOFPRO
#include "skewness.F"
EOFPRO
cat >smoothsarspec.F <<\EOFPRO
#include "smoothsarspec.F"
EOFPRO
cat >snonlin.F <<\EOFPRO
#include "snonlin.F"
EOFPRO
cat >spec2fdb.F <<\EOFPRO
#include "spec2fdb.F"
EOFPRO
cat >spectra.F <<\EOFPRO
#include "spectra.F"
EOFPRO
cat >splitbl.F <<\EOFPRO
#include "splitbl.F"
EOFPRO
cat >splitgrs.F <<\EOFPRO
#include "splitgrs.F"
EOFPRO
cat >spr.F <<\EOFPRO
#include "spr.F"
EOFPRO
cat >sthq.F <<\EOFPRO
#include "sthq.F"
EOFPRO
cat >stress.F <<\EOFPRO
#include "stress.F"
EOFPRO
cat >stresso.F <<\EOFPRO
#include "stresso.F"
EOFPRO
cat >strspec.F <<\EOFPRO
#include "strspec.F"
EOFPRO
cat >sumenergy.F <<\EOFPRO
#include "sumenergy.F"
EOFPRO
cat >swellsep.F <<\EOFPRO
#include "swellsep.F"
EOFPRO
cat >tauhf.F <<\EOFPRO
#include "tauhf.F"
EOFPRO
cat >threshhold.F <<\EOFPRO
#include "threshhold.F"
EOFPRO
cat >timin.F <<\EOFPRO
#include "timin.F"
EOFPRO
cat >topoar.F <<\EOFPRO
#include "topoar.F"
EOFPRO
cat >transmeans.F <<\EOFPRO
#include "transmeans.F"
EOFPRO
cat >transpart.F <<\EOFPRO
#include "transpart.F"
EOFPRO
cat >tustreas.F <<\EOFPRO
#include "tustreas.F"
EOFPRO
cat >u10cor.F <<\EOFPRO
#include "u10cor.F"
EOFPRO
cat >uiprep.F <<\EOFPRO
#include "uiprep.F"
EOFPRO
cat >update.F <<\EOFPRO
#include "update.F"
EOFPRO
cat >updatewd.F <<\EOFPRO
#include "updatewd.F"
EOFPRO
cat >upwspec.F <<\EOFPRO
#include "upwspec.F"
EOFPRO
cat >userin.F <<\EOFPRO
#include "userin.F"
EOFPRO
cat >vmin.F <<\EOFPRO
#include "vmin.F"
EOFPRO
cat >vplus.F <<\EOFPRO
#include "vplus.F"
EOFPRO
cat >wam_user_clock.F <<\EOFPRO
#include "wam_user_clock.F"
EOFPRO
cat >wamassi.F <<\EOFPRO
#include "wamassi.F"
EOFPRO
cat >wamcur.F <<\EOFPRO
#include "wamcur.F"
EOFPRO
cat >wamnorm.F <<\EOFPRO
#include "wamnorm.F"
EOFPRO
cat >wamodel.F <<\EOFPRO
#include "wamodel.F"
EOFPRO
cat >wamwnd.F <<\EOFPRO
#include "wamwnd.F"
EOFPRO
cat >wavemdl.F <<\EOFPRO
#include "wavemdl.F"
EOFPRO
cat >windsea.F <<\EOFPRO
#include "windsea.F"
EOFPRO
cat >writefl.F <<\EOFPRO
#include "writefl.F"
EOFPRO
cat >wposnam.F <<\EOFPRO
#include "wposnam.F"
EOFPRO
cat >writestress.F <<\EOFPRO
#include "writestress.F"
EOFPRO
cat >writsta.F <<\EOFPRO
#include "writsta.F"
EOFPRO
cat >wsmfen.F <<\EOFPRO
#include "wsmfen.F"
EOFPRO
cat >wvalloc.F <<\EOFPRO
#include "wvalloc.F"
EOFPRO
cat >wvdealloc.F <<\EOFPRO
#include "wvdealloc.F"
EOFPRO
cat >wdirspread.F <<\EOFPRO
#include "wdirspread.F"
EOFPRO
cat >wstream_strg.F <<\EOFPRO
#include "wstream_strg.F"
EOFPRO
cat >i_get_unit.F <<\EOFPRO
#include "../Alt/i_get_unit.F"
EOFPRO
cat >uraopc.F <<\EOFPRO
#include "../Alt/uraopc.F"
EOFPRO
cat >sortini.F <<\EOFPRO
#include "../Alt/sortini.F"
EOFPRO
cat >sorti.F <<\EOFPRO
#include "../Alt/sorti.F"
EOFPRO
cat >accuwa.F <<\EOFPRO
#include "../Sar/accuwa.F"
EOFPRO
cat >aftfft.F <<\EOFPRO
#include "../Sar/aftfft.F"
EOFPRO
cat >calibra.F <<\EOFPRO
#include "../Sar/calibra.F"
EOFPRO
cat >chasar.F <<\EOFPRO
#include "../Sar/chasar.F"
EOFPRO
cat >chpar.F <<\EOFPRO
#include "../Sar/chpar.F"
EOFPRO
cat >compoutpar.F <<\EOFPRO
#include "../Sar/compoutpar.F"
EOFPRO
cat >costp.F <<\EOFPRO
#include "../Sar/costp.F"
EOFPRO
cat >fftsub.F <<\EOFPRO
#include "../Sar/fftsub.F"
EOFPRO
cat >fimagep.F <<\EOFPRO
#include "../Sar/fimagep.F"
EOFPRO
cat >findbb.F <<\EOFPRO
#include "../Sar/findbb.F"
EOFPRO
cat >fttok.F <<\EOFPRO
#include "../Sar/fttok.F"
EOFPRO
cat >fullmap.F <<\EOFPRO
#include "../Sar/fullmap.F"
EOFPRO
cat >getsarpar.F <<\EOFPRO
#include "../Sar/getsarpar.F"
EOFPRO
cat >iners1.F <<\EOFPRO
#include "../Sar/iners1.F"
EOFPRO
cat >isoplot.F <<\EOFPRO
#include "../Sar/isoplot.F"
EOFPRO
cat >mpgathersarinv.F <<\EOFPRO
#include "../Sar/mpgathersarinv.F"
EOFPRO
cat >newspec.F <<\EOFPRO
#include "../Sar/newspec.F"
EOFPRO
cat >prepfft.F <<\EOFPRO
#include "../Sar/prepfft.F"
EOFPRO
cat >rarmtf.F <<\EOFPRO
#include "../Sar/rarmtf.F"
EOFPRO
cat >readwsp.F <<\EOFPRO
#include "../Sar/readwsp.F"
EOFPRO
cat >sarinvert.F <<\EOFPRO
#include "../Sar/sarinvert.F"
EOFPRO
cat >sarpc.F <<\EOFPRO
#include "../Sar/sarpc.F"
EOFPRO
cat >sarpol.F <<\EOFPRO
#include "../Sar/sarpol.F"
EOFPRO
cat >select.F <<\EOFPRO
#include "../Sar/select.F"
EOFPRO
cat >subs_check.F <<\EOFPRO
#include "../Sar/subs_check.F"
EOFPRO
cat >usage_sar.F <<\EOFPRO
#include "../Sar/usage_sar.F"
EOFPRO
cat >uwaror.F <<\EOFPRO
#include "../Sar/uwaror.F"
EOFPRO
cat >vbmtf.F <<\EOFPRO
#include "../Sar/vbmtf.F"
EOFPRO
cat >waven.F <<\EOFPRO
#include "../Sar/waven.F"
EOFPRO
cat >wrires.F <<\EOFPRO
#include "../Sar/wrires.F"
EOFPRO
cat >wtail.F <<\EOFPRO
#include "../Sar/wtail.F"
EOFPRO

#.#####################################################################
# get message passing routines used at ECMWF
# They are defined using fortran modules and are therefore needed
# for the compilation of the WAM code.
# They can be found in the ifsaux library
# They also use common utility routines from the ifs library which
# are are replaced by a similar routine in libwam.
# This is an interface to MPI.

# The mpl modules (and associated modules)
##############################################

cat >parkind1.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/parkind1.F90"
EOFPRO

cat >parkind2.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/parkind2.F90"
EOFPRO

cat >mpl_mpif.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_mpif.F90"
EOFPRO

cat >mpl_data_module.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_data_module.F90"
EOFPRO

cat >mpl_abort_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_abort_mod.F90"
EOFPRO

cat >mpl_message_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_message_mod.F90"
EOFPRO

cat >mpl_ioinit_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_ioinit_mod.F90"
EOFPRO

cat >mpl_allgather_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_allgather_mod.F90"
EOFPRO

cat >mpl_allgatherv_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_allgatherv_mod.F90"
EOFPRO

cat >mpl_allreduce_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_allreduce_mod.F90"
EOFPRO

cat >mpl_alltoallv_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_alltoallv_mod.F90"
EOFPRO

cat >mpl_barrier_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_barrier_mod.F90"
EOFPRO

cat >mpl_broadcast_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_broadcast_mod.F90"
EOFPRO

cat >mpl_buffer_method_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_buffer_method_mod.F90"
EOFPRO

cat >mpl_close_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_close_mod.F90"
EOFPRO

cat >mpl_comm_create_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_comm_create_mod.F90"
EOFPRO

cat >mpl_end_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_end_mod.F90"
EOFPRO

cat >mpl_gatherv_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_gatherv_mod.F90"
EOFPRO

cat >mpl_groups.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_groups.F90"
EOFPRO

cat >mpl_init_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_init_mod.F90"
EOFPRO

cat >mpl_module.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_module.F90"
EOFPRO

cat >mpl_myrank_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_myrank_mod.F90"
EOFPRO

cat >mpl_nproc_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_nproc_mod.F90"
EOFPRO

cat >mpl_open_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_open_mod.F90"
EOFPRO

cat >mpl_probe_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_probe_mod.F90"
EOFPRO

cat >mpl_read_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_read_mod.F90"
EOFPRO

cat >mpl_recv_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_recv_mod.F90"
EOFPRO

cat >mpl_scatterv_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_scatterv_mod.F90"
EOFPRO

cat >mpl_send_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_send_mod.F90"
EOFPRO

cat >mpl_tour_table_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_tour_table_mod.F90"
EOFPRO

cat >mpl_wait_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_wait_mod.F90"
EOFPRO

cat >mpl_write_mod.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_write_mod.F90"
EOFPRO

cat >mpl_module.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/mpl_module.F90"
EOFPRO

cat >yomgstats.F90 <<\EOFPRO
#include "/ccvobs/ifsaux/module/yomgstats.F90"
EOFPRO



# configuration and options files for IBM power3 (IBM servers)
# and IBM hpcc and sgi
# as used by ALL ECMWF project.
# So they are quite extensive (WAM might actually need a lot less !
#####################################################################

cat >../build/config.free.rs6000 <<\EOFPRO
#include "/ccvobs/ifsaux/build/config.free.rs6000"
EOFPRO

cat >../build/config.fixed.rs6000 <<\EOFPRO
#include "/ccvobs/ifsaux/build/config.fixed.rs6000"
EOFPRO

cat >../build/options.rs6000 <<\EOFPRO
#include "/ccvobs/ifsaux/build/options.rs6000"
EOFPRO

cat >../build/config.free.ibm_power4 <<\EOFPRO
#include "/ccvobs/ifsaux/build/config.free.rs6000"
EOFPRO

cat >../build/config.fixed.ibm_power4 <<\EOFPRO
#include "/ccvobs/ifsaux/build/config.fixed.rs6000"
EOFPRO

cat >../build/options.ibm_power4 <<\EOFPRO
#include "/ccvobs/ifsaux/build/options.rs6000"
EOFPRO

cat >../build/config.free.sgi <<\EOFPRO
#include "/ccvobs/ifsaux/build/config.free.rs6000"
EOFPRO

cat >../build/config.fixed.sgi <<\EOFPRO
#include "/ccvobs/ifsaux/build/config.fixed.rs6000"
EOFPRO

cat >../build/options.sgi <<\EOFPRO
#include "/ccvobs/ifsaux/build/options.rs6000"
EOFPRO


######################
#. compile MPL modules
#######################

# The makefile for the mpl modules (my own version)
####################################################
cat >makefile_mpl <<\EOFPRO
include ../build/config.free.$(ARCH)
include ../build/options.$(ARCH)
SOURCE	= parkind1.F90 parkind2.F90 mpl_abort_mod.F90 mpl_allgather_mod.F90  mpl_allgatherv_mod.F90 mpl_alltoallv_mod.F90 mpl_barrier_mod.F90 mpl_broadcast_mod.F90 mpl_buffer_method_mod.F90 mpl_comm_create_mod.F90 mpl_data_module.F90 mpl_end_mod.F90 mpl_groups.F90 mpl_init_mod.F90 mpl_message_mod.F90 mpl_module.F90 mpl_myrank_mod.F90 mpl_probe_mod.F90 mpl_recv_mod.F90 mpl_scatterv_mod.F90 mpl_send_mod.F90 mpl_tour_table_mod.F90 mpl_mpif.F90 mpl_allreduce_mod.F90 mpl_ioinit_mod.F90 mpl_open_mod.F90 mpl_close_mod.F90 mpl_read_mod.F90 mpl_write_mod.F90 mpl_gatherv_mod.F90 mpl_wait_mod.F90  mpl_nproc_mod.F90 yomgstats.F90
OBJECTS1	=	$(SOURCE:.F=.o)
OBJECTS2	=	$(SOURCE:.F90=.o)
OBJECTS = $(OBJECTS1) $(OBJECTS2)
.SUFFIXES:
.SUFFIXES:	.o .F .F90
.F.o :
	$(FC) -c $(FFLAGS) $(CP_OPT) $*.F
.F90.o :
	$(FC) -c $(FFLAGS) $(CP_OPT) $*.F90
all:	$(OBJECTS)
	touch module.a
clean:
	 \rm -f $(OBJECTS)
# DO NOT DELETE THIS LINE - used by make depend
mpl_abort_mod.o: mpl_data_module.mod mpl_mpif.mod parkind1.mod
mpl_allgather_mod.o: mpl_barrier_mod.mod mpl_data_module.mod mpl_message_mod.mod mpl_mpif.mod mpl_recv_mod.mod mpl_send_mod.mod parkind1.mod
mpl_allgatherv_mod.o: mpl_barrier_mod.mod mpl_data_module.mod mpl_message_mod.mod mpl_mpif.mod mpl_recv_mod.mod mpl_send_mod.mod parkind1.mod
mpl_allreduce_mod.o: mpl_broadcast_mod.mod mpl_data_module.mod mpl_message_mod.mod mpl_mpif.mod mpl_recv_mod.mod mpl_send_mod.mod parkind1.mod
mpl_alltoallv_mod.o: mpl_data_module.mod mpl_message_mod.mod mpl_mpif.mod parkind1.mod
mpl_barrier_mod.o: mpl_data_module.mod mpl_message_mod.mod parkind1.mod
mpl_broadcast_mod.o: mpl_barrier_mod.mod mpl_data_module.mod mpl_message_mod.mod mpl_mpif.mod mpl_recv_mod.mod mpl_send_mod.mod parkind1.mod mpl_myrank_mod.mod
mpl_buffer_method_mod.o: mpl_data_module.mod mpl_message_mod.mod parkind1.mod
mpl_close_mod.o: mpl_data_module.mod mpl_ioinit_mod.mod mpl_message_mod.mod parkind1.mod
mpl_comm_create_mod.o: mpl_data_module.mod mpl_message_mod.mod mpl_mpif.mod parkind1.mod
mpl_end_mod.o: mpl_data_module.mod mpl_message_mod.mod parkind1.mod
mpl_gatherv_mod.o: mpl_data_module.mod mpl_message_mod.mod mpl_mpif.mod parkind1.mod
mpl_groups.o: mpl_data_module.mod mpl_message_mod.mod mpl_mpif.mod parkind1.mod
mpl_init_mod.o: mpl_buffer_method_mod.mod mpl_data_module.mod mpl_message_mod.mod mpl_mpif.mod mpl_tour_table_mod.mod parkind1.mod
mpl_ioinit_mod.o: mpl_data_module.mod mpl_message_mod.mod mpl_mpif.mod parkind1.mod
mpl_message_mod.o: mpl_abort_mod.mod mpl_data_module.mod mpl_mpif.mod parkind1.mod
mpl_module.o: mpl_abort_mod.mod mpl_allgather_mod.mod mpl_allgatherv_mod.mod mpl_allreduce_mod.mod mpl_barrier_mod.mod mpl_broadcast_mod.mod mpl_buffer_method_mod.mod mpl_close_mod.mod mpl_comm_create_mod.mod mpl_end_mod.mod mpl_gatherv_mod.mod mpl_groups.mod mpl_init_mod.mod mpl_ioinit_mod.mod mpl_message_mod.mod mpl_myrank_mod.mod mpl_nproc_mod.mod mpl_open_mod.mod mpl_probe_mod.mod mpl_read_mod.mod mpl_recv_mod.mod mpl_scatterv_mod.mod mpl_send_mod.mod mpl_tour_table_mod.mod mpl_wait_mod.mod mpl_write_mod.mod
mpl_myrank_mod.o: mpl_data_module.mod mpl_message_mod.mod parkind1.mod
mpl_nproc_mod.o: mpl_data_module.mod mpl_message_mod.mod parkind1.mod
mpl_open_mod.o: mpl_data_module.mod mpl_ioinit_mod.mod mpl_message_mod.mod mpl_mpif.mod parkind1.mod
mpl_probe_mod.o: mpl_data_module.mod mpl_message_mod.mod mpl_mpif.mod parkind1.mod
mpl_read_mod.o: mpl_data_module.mod mpl_ioinit_mod.mod mpl_message_mod.mod mpl_mpif.mod parkind1.mod
mpl_recv_mod.o: mpl_data_module.mod mpl_message_mod.mod mpl_mpif.mod parkind1.mod
mpl_scatterv_mod.o: mpl_data_module.mod mpl_message_mod.mod mpl_mpif.mod mpl_myrank_mod.mod parkind1.mod
mpl_send_mod.o: mpl_data_module.mod mpl_message_mod.mod mpl_mpif.mod parkind1.mod
mpl_tour_table_mod.o: mpl_data_module.mod mpl_message_mod.mod parkind1.mod
mpl_wait_mod.o: mpl_data_module.mod mpl_message_mod.mod mpl_mpif.mod parkind1.mod
mpl_write_mod.o: mpl_data_module.mod mpl_ioinit_mod.mod mpl_message_mod.mod mpl_mpif.mod parkind1.mod

mpl_abort_mod.mod: ./mpl_abort_mod.o
mpl_allgather_mod.mod: ./mpl_allgather_mod.o
mpl_allgatherv_mod.mod: ./mpl_allgatherv_mod.o
mpl_allreduce_mod.mod: ./mpl_allreduce_mod.o
mpl_alltoallv_mod.mod: ./mpl_alltoallv_mod.o
mpl_barrier_mod.mod: ./mpl_barrier_mod.o
mpl_broadcast_mod.mod: ./mpl_broadcast_mod.o
mpl_buffer_method_mod.mod: ./mpl_buffer_method_mod.o
mpl_close_mod.mod: ./mpl_close_mod.o
mpl_comm_create_mod.mod: ./mpl_comm_create_mod.o
mpl_data_module.mod: ./mpl_data_module.o
mpl_end_mod.mod: ./mpl_end_mod.o
mpl_gatherv_mod.mod: ./mpl_gatherv_mod.o
mpl_groups.mod: ./mpl_groups.o
mpl_init_mod.mod: ./mpl_init_mod.o
mpl_ioinit_mod.mod: ./mpl_ioinit_mod.o
mpl_message_mod.mod: ./mpl_message_mod.o
mpl_module.mod: ./mpl_module.o
mpl_mpif.mod: ./mpl_mpif.o
mpl_myrank_mod.mod: ./mpl_myrank_mod.o
mpl_nproc_mod.mod: ./mpl_nproc_mod.o
mpl_open_mod.mod: ./mpl_open_mod.o
mpl_probe_mod.mod: ./mpl_probe_mod.o
mpl_read_mod.mod: ./mpl_read_mod.o
mpl_recv_mod.mod: ./mpl_recv_mod.o
mpl_scatterv_mod.mod: ./mpl_scatterv_mod.o
mpl_send_mod.mod: ./mpl_send_mod.o
mpl_tour_table_mod.mod: ./mpl_tour_table_mod.o
mpl_wait_mod.mod: ./mpl_wait_mod.o
mpl_write_mod.mod: ./mpl_write_mod.o
parkind1.mod: ./parkind1.o
parkind2.mod: ./parkind2.o
yomgstats.mod: ./yomgstats.o
EOFPRO


gmake -f makefile_mpl

#.
# The makefile for the wam code 
################################
cat >makefile_wam <<\EOFPRO
include ../build/config.fixed.$(ARCH)
include ../build/options.$(ARCH)
OBJECTS	=	$(SOURCE:.F=.o)
.SUFFIXES:
.SUFFIXES:	.o .F .f
.F.o :
	$(FC) -c $(FFLAGS) $(CP_OPT) $*.F
all:	$(OBJECTS)
	touch module.a
clean:
	 \rm -f $(OBJECTS)
EOFPRO


#######################
# compile wam modules
#######################
for file in $(ls -1 yow*.F)
do
  gmake -f makefile_wam ${file%.F}.o
done

#######################
# compile wam routines
#######################
for file in $(ls -1 *.F)
do
  gmake -f makefile_wam ${file%.F}.o
done

#################
#. build library
#################
#.
\rm -f ${LIBS}/libwam.a
if [[ $(ls -1 *.o | wc -l ) -gt 0 ]] ; then
  ar rv ${LIBS}/libwam.a *.o 
fi

chmod go+rx ${LIBS}/libwam.a

######################
# create excecutables
######################

### necessary libraries ###
###########################
# ECLIB: needed for the time handeling routines secincr and secdiff (you can replace them
#        by yours). It is also needed for matrix inversion routine syminv used in the data
#        assimilation. A dummy definition of it is given here since the example given
#        with this code is for a wave model run without data assimilation. If you want
#        to run data assimilation get the proper definiton from eclib or use your own.

\rm -f ${LIBS}/libec.a
ln -s /usr/local/lib/eclib/libec_FORTRAN_WITH_UNDERSCORE.R64.D64.I32.a ${LIBS}/libec.a

# EMOSLIB: needed for grib coding/decoding and manipulation of binary data. WAM can in
#          principle be run without grib access using pure fortran binary files. If it is
#          the case just create extra dummy routines for the definition of the emoslib
#          routines and refer to the run example to switch the appropriate flags to have
#          pure binary input/output.

\rm -f ${LIBS}/libemos.a 
ln -s /home/ma/emos/lib/rs6000/000200/libemos.R64.D64.I32.a ${LIBS}/libemos.a

# NOTE: if SAR data assimilation was required you will need the following external library
# NAG library for c06fuf anf the LAPACK library for dgelss. Dummy definitions of these two
# routines are provided here.
# Dummy definitions are also provided for the field data base (FDB) library used at ECMWF
# to output data from model runs instead the example provided here will write the output
# to normal files.  

# IMPORTANT NOTE: you will also need the message passing library (mpi). In the ECMWF setup
#                 the mpi library is automatically loaded by using the mpxlf90_r command.
#                 If you are running outside ECMWF you will need to get your own version.
#                 Note also that you might have to modify the content of mpl_mpif.F90 
#                 which contains information on how the compile mpi. The version provided
#                 here is for the IBM servers.


cd ${BINS}

for file in chief preproc preset getinptb
do
  ar xv ${LIBS}/libwam.a ${file}.o
  mpxlf90_r -o ${file} -bmaxdata:0x70000000 -qsmp=noauto ${file}.o \
  -L${LIBS} -lwam -lec -lemos
  \rm ${file}.o
done

cd ..
\rm -rf tmp$$

#endif
