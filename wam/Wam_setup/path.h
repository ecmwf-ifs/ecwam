#.
#.               USER DEPENDENT SETTINGS.
#.               TO BE MODIFIED AS NESSECARY.
#.
#.=====================================================================
#.  Define paths for ECFS and VPP
#.=====================================================================
#.
VERSION=CY19R2_modules
#.
DHSROOT=/$USER/vpp700/wam_$VERSION
VPPROOT=/vpp700/wavedata44/${USER}/wam_$VERSION
TMPROOT=/vpp700/wavedir44/${USER}/wam_$VERSION
#.
DHSLPATH=$DHSROOT/${VERSION}
LIBS=$TMPROOT/lib
WORKDIR=tmp$$
#.
#ifndef region
#define region 'z'
#endif
#.
#if region=='s'
DHSPATH=$DHSROOT/swamp2/coarse
ADIR=$TMPROOT/swamp2/coarse
BIN=$TMPROOT/swamp2/coarse/bin
ROOTWDIR=$VPPROOT/swamp2/coarse
#elif  region=='m'
DHSPATH=$DHSROOT/medite
ADIR=$TMPROOT/medite
BIN=$TMPROOT/medite/bin
ROOTWDIR=$VPPROOT/medite
#elif region=='g'
DHSPATH=$DHSROOT/global
ADIR=$TMPROOT/global
BIN=$TMPROOT/global/bin
ROOTWDIR=$VPPROOT/global
#else
DHSPATH=$DHSROOT
ADIR=$TMPROOT
BIN=$LIBS
ROOTWDIR=$VPPROOT
#endif
#.
#if resolution==150
GRID="1.5/1.5"
AREA="81./ 0./ -81./358.5"
grid="150"
nang=12
nfre=25
#elif region=='g' &&  resolution==50
GRID="0.5/0.5"
AREA="81./ 0./ -81./359.5"
nang=12
nfre=25
grid="050"
#elif region=='g' &&  resolution==300
GRID="3.0/3.0"
AREA="72./ 0./ -63./357.0"
grid="300"
nang=12
nfre=25
#elif region=='m' && resolution==25
GRID="0.25/0.25"
AREA="81./ -98./ 9./42."
grid="025"
nang=24
nfre=25
#elif region=='m' && resolution==50
GRID="0.5/0.5"
AREA="46./ -6./ 30./36."
grid="050"
nang=24
nfre=25
#elif region=='g' && resolution==900
GRID="9.0/9.0"
AREA="72.0/0.0/-63.0/351.0"
grid="900"
nang=12
nfre=25
#else
   banner QUATSCH
#endif
#.
#.
[[ ! -d $ADIR  ]] && mkdir -p $ADIR
[[ ! -d $ROOTWDIR  ]] && mkdir -p $ROOTWDIR
[[ ! -d $LIBS  ]] && mkdir -p $LIBS
[[ ! -d $BIN  ]] && mkdir -p $BIN

WDIR=${ROOTWDIR}/${grid}
[[ ! -d $WDIR  ]] && mkdir -p $WDIR
#.
#.
#.=====================================================================
#.  Define user input 
#.=====================================================================
#.
ASSIMILATION=NO                        # altimeter data assimilation 
#
typeset -Z12 begofrn=199901011200      # BEGin date OF RUn 
typeset -Z12 endofrn=199901071200      # END   date OF RUn 
typeset -Z12 begoffo=199901021200      # BEGin date OF FOrcast.
#                                        This date must equal to endofrn
#                                        when analysis is only required
typeset -Z12 outofrf=199901021200      # Date to output restart file(s).
                                       # Set to 0000000000 if determined
                                       # on userinput.
typeset -Z12 outof2d=199901021200      # Date up to which 2D-spectra are
                                       # saved. Set to 0000000000 if not
                                       # required.
typeset -Z7 antime=86400               # Length of analysis in seconds.
typeset -Z7 fctime=432000              # Length of forcast in seconds.
#.
CLASS=RD                               # user class
NENS=000                               # only used for ensemble run.
TNE=000                                # only used for ensemble run.

#.! expver is the experiment id. It is used when gribbing the output data. 
typeset -l expver=${USER}a
#.
#.define file storage directory, as needed for output destination
#.
export STORAGE_PATH=$WDIR
#.
#. define fbd server on vpp-700
#. ============================
#
#.export FDB_CONFIG_FILE=/vpp700/mrfs/fdb/FDbConfig_vpp700
export FDB_ROOT=/vpp700/fdb48
export FDB_CONFIG_MODE=async
export FDB_SERVER_HOST=vpp700-x18
export FDB_SIGNALS=no
export FDB_DEBUG=no
#.
cfdb2dsp=$FDB_ROOT                      # fdb for 2d spectra.
cfdbsf=$FDB_ROOT                        # fdb for scalar fields.
#.
#
# ====================================================================
# DEFINE PRECISION OF REALS (SINGLE or DOUBLE).
# FORCE REMAKE OF LIBWAM (default NO).
# =====================================================================
#.
export PRECISION=SINGLE
export REMAKE=NO
#.
#.======================
#.  Define libraries.
#.======================
#.
aux=/home/rd/rdx/lib/18r3/libifsaux.a
dum=/home/rd/rdx/lib/18r3/libdummy.a
WAMLIB=libwam
LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib:$LIBS
wam=$LIBS/$WAMLIB.a
FDBLIB=-lifsio
#.
#.     END OF    USER DEPENDENT SETTINGS:
#.
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
  export rp=-CcdRR8
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
