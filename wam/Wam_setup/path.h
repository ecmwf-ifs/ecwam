#.
#.               USER DEPENDENT SETTINGS.
#.               TO BE MODIFIED AS NECESSARY.
#.
#.=====================================================================
#.  Define paths 
#.=====================================================================
#.
#. Top directories on VPP
VPPROOT=${TEMP}/wam_${USER}_$VERSION
TMPROOT=${TEMP}/wam_${USER}_$VERSION
#.
#. Location of the compiled WAM library on VPP
LIBS=$TMPROOT/lib
#. Location of the compiled WAM library on workstation
WKLIBS=/ws/scratch/rd/${USER}/lib
#. 
WORKDIR=tmp$$
#.
#. Different directories on VPP
#ifndef region
#define region 'z'
#endif
#.
#if region=='s'
ADIR=$TMPROOT/swamp2/coarse
BIN=$TMPROOT/swamp2/coarse/bin
ROOTWDIR=$VPPROOT/swamp2/coarse
#elif  region=='m'
ADIR=$TMPROOT/medite
BIN=$TMPROOT/medite/bin
ROOTWDIR=$VPPROOT/medite
#elif region=='g'
ADIR=$TMPROOT/global
BIN=$TMPROOT/global/bin
ROOTWDIR=$VPPROOT/global
#else
ADIR=$TMPROOT
BIN=$LIBS
ROOTWDIR=$VPPROOT
#endif
#.
#.=====================================================================
#.  Define grid characteristics
#.=====================================================================
#if resolution==150
GRID="1.5/1.5"
AREA="81./ 0./ -81./358.5"
grid="150"
nang=12
nfre=25
#elif region=='g' &&  resolution==50
GRID="0.5/0.5"
AREA="81./ 0./ -81./359.5"
nang=24
nfre=30
grid="050"
#elif region=='g' &&  resolution==100
GRID="1.0/1.0"
AREA="81./ 0./ -78./359.0"
nang=12
nfre=25
grid="100"
#elif region=='g' &&  resolution==300
GRID="3.0/3.0"
AREA="81./ 0./ -78./357.0"
grid="300"
nang=12
nfre=25
#elif region=='m' && resolution==25
GRID="0.25/0.25"
AREA="48./ -6./ 30./42."
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
#.
# This example runs 24 hrs with forecast winds
typeset -Z10 start_date=2002010112     # starting date (YYYYMMDDHH)
typeset -i anlength=0                  # Length of analysis in hours
typeset -i fclength=24                 # Length of forcast in hours
#.
#.
typeset -i anlen=anlength*3600
typeset -Z7 antime=$anlen
typeset -i fclen=fclength*3600
typeset -Z7 fctime=$fclen
beginfcdt=$(newdate $start_date $anlength)
enddt=$(newdate $beginfcdt $fclength)
typeset -Z12 begofrn=${start_date}00   # BEGin date OF RUn
typeset -Z12 endofrn=${enddt}00        # END   date OF RUn
typeset -Z12 begoffo=${beginfcdt}00    # BEGin date OF FOrcast.
#                                        This date must equal to endofrn
#                                        when analysis is only required
typeset -Z12 outofrf=${endofrn}        # Date to output restart binary files
                                       # or the gribbed spectra if grib output
                                       # was requested. This date is on top of
                                       # all other dates that might be set by
                                       # cycling from the initial date with
                                       # step given by IDELRES.
                                       # Set it to 0000000000 if determined
                                       # by namelist NAOS (not used in 
                                       # this example)
#.
CLASS=RD                               # user class (used when gribbing data)
#.! expver is the experiment id. It is used when gribbing the output data. 
typeset -l expver=wave
#.
#if region=='s'
ASSIMILATION=NO                        # data assimilation 
iassi=0
laltas=F
laltcor=T
lsarinv=F
lsaras=F
#else
ASSIMILATION=NO                        # data assimilation 
#. Specify what type of data will be used for the assimilation
iassi=0
laltas=T                               # altimeter
laltcor=T                              # altimeter data correction
lsarinv=F                              # SAR inversion
lsaras=F                               # SAR
#endif
#.
#.define file storage directory, as needed for output destination
##################################################################
#.
export STORAGE_PATH=$WDIR
#.
#.
#. define fbd server on vpp-700
#. ============================
#
export FDB_ROOT=/vpp700/fdb4c
export FDB_CONFIG_MODE=async
export FDB_SERVER_HOST=vpp700-x4c
export FDB_SIGNALS=no
export FDB_DEBUG=no
#.
# ====================================================================
# DEFINE PRECISION OF REALS (SINGLE or DOUBLE).
# FORCE REMAKE OF LIBWAM (default NO).
# =====================================================================
#.
export PRECISION=DOUBLE
export REMAKE=YES
#.
#.======================
#.  Define libraries.
#.======================
#.
LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib:$LIBS
wam=$LIBS/$WAMLIB.a
#.
#.     END OF    USER DEPENDENT SETTINGS:
#.
#.
if [ $PRECISION = SINGLE ] ; then
#. not currently used
  exit 1
elif [ $PRECISION = DOUBLE ] ; then
  libselect -r 64
  export rp=-Ad
#.
#.the ECMWF message passing interface is part of auxiliary ifs library 
#.We use the defaults version corresponding to the extracted cycle
#.(see ifs_cycle)
#.If you are using WAM outside ECMWF you will need to get the
#.appropriate version.
  mpllib=/vpp700/rdx_dir04/xroot_rd/lib/${ifs_cycle}/libifsaux.a
#. Associated to the mpl routines are the MPI libraries
  mpilib=/usr/lang/mpi2/lib/libmpi.a
  mp2lib=/opt/tools/lib/libmp2.2.1x.a
#.The library EMOSLIB is used for coding and decoding data in GRIB and BUFR
 emoslib=$EMOSLIB
#.The library ECLIB is used for the routine syminv and date and time
#.operations
 eclib=$ECLIB
#.The library FDBLIB is used for writing grib data to the ECMWF field data
#.base 
 fdblib=$FDBLIB
#.NAGLIB and LAPACKLIB and BLASLIB are only used by code specific 
#.to SAR assimilation and to insure everything compile on VPP
#.The SAR assimilation software should not be used since it was not
#.yet operational.
 naglib=$NAGLIB
 lapacklib=/usr/local/lib/liblapack.a
 blaslib=/usr/local/lib/libblas.a
#. CVPLIB is only needed because the vpp will not load without a definition
#. for routines
 cvplib=/usr/lang/lib/libcvp.a
#.
else
  ls -lsa
  print - "\n\n\t\tPRECISION $PRECISION not valid\n\n"
  print - "\n\n\t\tUSE SINGLE or DOUBLE \n\n"
  print - "\n\n\t\tPROGRAM IS TERMINATED\n\n"
  exit 1
fi
