#.
#.               USER DEPENDENT SETTINGS.
#.               TO BE MODIFIED AS NESSECARY.
#.
#.=====================================================================
#.  Define pathes for ECFS and VPP
#.=====================================================================
#.
VERSION=4r7_CY18R4_test
#.
DHSROOT=/$USER/vpp700/wam_CY18R4_test
VPPROOT=$LTEMP/wam_CY18R4_test
#.
DHSLPATH=$DHSROOT/${VERSION}
VPPLPATH=$VPPROOT/${VERSION}
#.
#ifndef region
#define region 'z'
#endif
#.
#if region=='s'
DHSPATH=$DHSROOT/${VERSION}_swamp2/coarse
DHSPAT1=$DHSROOT/${VERSION}_swamp2/fine
VPPPATH=$VPPROOT/${VERSION}_swamp2/coarse
VPPPAT1=$VPPROOT/${VERSION}_swamp2/fine
[[ ! -d $VPPPAT1 ]] && mkdir $VPPPAT1
#elif  region=='m'
DHSPATH=$DHSROOT/${VERSION}_medite
VPPPATH=$VPPROOT/${VERSION}_medite
#elif region=='g'
DHSPATH=$DHSROOT/${VERSION}_global
VPPPATH=$VPPROOT/${VERSION}_global
#else
DHSPATH=$DHSROOT/${VERSION}
VPPPATH=$VPPROOT/${VERSION}
#endif
#.
#.
WAMCRLIB=wamcr
LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib:$VPPLPATH
#.
[[ ! -d $VPPPATH  ]] && mkdir -p $VPPPATH
[[ ! -d $VPPLPATH ]] && mkdir -p $VPPLPATH
#.
# ====================================================================
# DEFINE FILE STORAGE DIRECTORY, AS NEEDED FOR CALL TO GSFILE
# =====================================================================
#.
export STORAGE_PATH=$VPPPATH
#.
# ====================================================================
# DEFINE PRECISION OF REALS (SINGLE or DOUBLE).
# FORCE REMAKE OF LIBRARIES AND X IF MAKE_LIBS=YES (DEFAULT NO).
# =====================================================================
#.
export PRECISION=SINGLE
export REMAKE=YES
#.
#.=====================================================================
#.  Define all times.
#.=====================================================================
#.
typeset -Z10 begofan=9801241200      # BEGin date OF ANalysis.
typeset -Z10 endofan=9801241800      # END   date OF ANalysis.
typeset -Z10 endoffo=9801241800      # END   date OF FOrcast.
typeset -Z10 outofrf=9801241800      # Date to output restart files.
                                     # Set to 0000000000 if determined
                                     # on uerinput.
typeset -Z10 outof2d=0000000000      # Date up to which 2D-spectra are
                                     # saved. Set to 0000000000 if not
                                     # required.
typeset -Z7 antime=21600             # Length of analysis in seconds.
typeset -Z7 fctime=0             # Length of forcast in seconds.
#.
#.=====================================================================
#.  Define auxiliary and dummy libraries.
#.=====================================================================
#.
#.
aux=/home/rd/rdx/lib/18r3/libifsaux.a
dum=/home/rd/rdx/lib/18r3/libdummy.a
wamcr=$VPPLPATH/lib$WAMCRLIB.a
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
  export rp=-Ad
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

