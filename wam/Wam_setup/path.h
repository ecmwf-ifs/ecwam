#.
#.               USER DEPENDENT SETTINGS.
#.               TO BE MODIFIED AS NESSECARY.
#.
#.=====================================================================
#.  Define pathes for ECFS and VPP
#.=====================================================================
#.
VERSION=CY19R1_U10
#.
DHSROOT=/$USER/vpp700/wam_$VERSION
VPPROOT=/vpp700/wavedata44/${USER}/wam_$VERSION
TMPROOT=/vpp700/wavedir44/${USER}/wam_$VERSION
#.
#.
DHSLPATH=$DHSROOT/${VERSION}
LIBCR=$TMPROOT/lib
WORKDIR=tmp$$
#.
#ifndef region
#define region 'z'
#endif
#.
#if region=='s'
DHSPATH=$DHSROOT/swamp2/coarse
ADIR=$TMPROOT/swamp2/coarse
LIBS=$TMPROOT/swamp2/coarse/lib
BIN=$TMPROOT/swamp2/coarse/bin
WDIR=$VPPROOT/swamp2/coarse
nang=12
nfre=25
#elif  region=='m'
DHSPATH=$DHSROOT/medite
ADIR=$TMPROOT/medite
LIBS=$TMPROOT/medite/lib
BIN=$TMPROOT/medite/bin
WDIR=$VPPROOT/medite
nang=24
nfre=25
#elif region=='g'
DHSPATH=$DHSROOT/global
ADIR=$TMPROOT/global
LIBS=$TMPROOT/global/lib
BIN=$TMPROOT/global/bin
WDIR=$VPPROOT/global
nang=12
nfre=25
#else
DHSPATH=$DHSROOT
ADIR=$TMPROOT
LIBS=$LIBCR
BIN=$LIBCR
WDIR=$VPPROOT
#endif
#.
#.
WAMCRLIB=wamcr
LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib:$LIBS:$LIBCR
#.
[[ ! -d $ADIR  ]] && mkdir -p $ADIR
[[ ! -d $WDIR  ]] && mkdir -p $WDIR
[[ ! -d $LIBS  ]] && mkdir -p $LIBS
[[ ! -d $BIN  ]] && mkdir -p $BIN
[[ ! -d $LIBCR  ]] && mkdir -p $LIBCR
#.
# ====================================================================
# DEFINE FILE STORAGE DIRECTORY, AS NEEDED FOR CALL TO GSFILE
# =====================================================================
#.
export STORAGE_PATH=$WDIR
#.
# ====================================================================
# DEFINE PRECISION OF REALS (SINGLE or DOUBLE).
# FORCE REMAKE OF LIBRARIES AND X IF MAKE_LIBS=YES (DEFAULT NO).
# THE WHOLE COMPILATION CAN BE SWITCHED OFF COMPELTELY BY SETTING
# COMPILE=0 
# =====================================================================
#.
export PRECISION=DOUBLE
export REMAKE=NO
export COMPILE=1
#.
#.=====================================================================
#.  Define all times.
#.=====================================================================
#.
ASSIMILATION=NO
typeset -Z12 begofrn=199810251200      # BEGin date OF RUn 
typeset -Z12 endofrn=199810251800      # END   date OF RUn 
typeset -Z12 begoffo=199810251800      # BEGin date OF FOrcast.
#                                        This date must equal to endofrn
#                                        when analysis is only required
typeset -Z12 outofrf=199810251800      # Date to output restart files.
                                       # Set to 0000000000 if determined
                                       # on uerinput.
typeset -Z12 outof2d=199810251800      # Date up to which 2D-spectra are
                                       # saved. Set to 0000000000 if not
                                       # required.
typeset -Z7 antime=21600               # Length of analysis in seconds.
typeset -Z7 fctime=0                   # Length of forcast in seconds.
#.
#.=====================================================================
#.  Define auxiliary and dummy libraries.
#.=====================================================================
#.
#.
aux=/home/rd/rdx/lib/18r3/libifsaux.a
dum=/home/rd/rdx/lib/18r3/libdummy.a
wamcr=$LIBS/lib$WAMCRLIB.a
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
#.
CLASS=RD
NENS=000
TNE=000
typeset -l XID=${USER}a
#.
#.
#.====================================================================
# DEFINE FBD SERVER ON VPP-700
#.=====================================================================
#
#.export FDB_CONFIG_FILE=/vpp700/mrfs/fdb/FDbConfig_vpp700
export FDB_ROOT=/vpp700/fdb48           # /vpp700/fdb22
export FDB_CONFIG_MODE=async            # eps
export FDB_SERVER_HOST=vpp700-x18       # vpp700-x22   vpp700-x08  == wavefdb
export FDB_SIGNALS=no
export FDB_DEBUG=no
#.
cfdb2dsp=$FDB_ROOT                      # fdb for scalar fields.
cfdbsf=$FDB_ROOT                        # fdb for non grib 2d spectra.
#.
FDBLIB=-lifsio
#####

#.
#if resolution==150
  GRID="1.5/1.5"
  AREA="81./ 0./ -81./358.5"
  grid="150"
#elif region=='g' &&  resolution==50
  GRID="0.5/0.5"
  AREA="81./ 0./ -81./359.5"
  grid="050"
#elif region=='g' &&  resolution==300
  GRID="3.0/3.0"
  AREA="72./ 0./ -63./357.0"
  grid="300"
#elif region=='m' && resolution==25
  GRID="0.25/0.25"
  AREA="81./ -98./ 9./42."
  grid="025"
#elif region=='g' && resolution==900
  GRID="9.0/9.0"
  AREA="72.0/0.0/-63.0/351.0"
  grid="900"
#else
   banner QUATSCH
#endif
