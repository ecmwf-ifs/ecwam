#.
#.=====================================================================
#.  Define version
#.=====================================================================
#.
ECUSER=wab
VERSION=CY26R3_extract_WAM
#.
#. Top directory
TMPROOT=/scratch/rd/${USER}/wam_$VERSION

#. Location of the compiled WAM library
LIBS=${TMPROOT}/lib

#. Location of the excecutables
BINS=${TMPROOT}/bin
#.
#.Define the spatial and spectral resolution
#.############################################
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
nfre=30
#elif region=='m' && resolution==50
GRID="0.5/0.5"
AREA="46./ -6./ 30./36."
grid="050"
nang=24
nfre=30
#elif region=='g' && resolution==900
GRID="9.0/9.0"
AREA="72.0/0.0/-63.0/351.0"
grid="900"
nang=12
nfre=25
#endif

#. Running directories
#if region=='s'
ROOTWDIR=$TMPROOT/swamp2
WDIR=${ROOTWDIR}/coarse
#elif  region=='m'
ROOTWDIR=$TMPROOT/medite
WDIR=${ROOTWDIR}/${grid}
#elif region=='g'
ROOTWDIR=$TMPROOT/global
WDIR=${ROOTWDIR}/${grid}
#else
ROOTWDIR=$TMPROOT
WDIR=${ROOTWDIR}/${grid}
#endif

