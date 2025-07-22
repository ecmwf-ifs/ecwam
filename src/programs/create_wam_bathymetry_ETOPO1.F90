! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

PROGRAM CREATE_BATHY_ETOPO1

!     CREATES THE WAM BATHYMETRY AND THE REDUCTION FACTORS
!     DUE TO SUB-GRID BATHYMETRIC FEATURES USING
!     THE ETOPO1 DATA SET,

!     JEAN BIDLOT JUNE 2017

!     DESCRIPTION OF ETOPO1:      
!     THE DATA FILE CAN BE FOUND ON TAPE:
!     ec:/rdx/prepdata/etopo1/ETOPO1_Ice_g_int.xyz.gz
!     It was obtained from the NCEI web server
!     https://www.ngdc.noaa.gov/mgg/global/
!     GUNZIP THE FILE AND MAKE SURE THAT ETOPO1_Ice_g_int.xyz IS ONLINE
!
!     DEPTH AND ELEVATIONS ARE IN METERS


!!!!! BECAUSE ETOPO1 HAS A RESOLUTION At BEST OF 1/60 degrees, IT DOES NOT
!!!!! MAKE ANY SENSE TO COMPUTE THE SMALL ISLANDS OBSTRUCTIONS IF THE RESOLUTION
!!!!! IS OF THAT ORDER OR SMALLER !!!!
!!!!! SO IF DXDELLA < 2*1/60, THEN THE OBSTRUCTIONS WILL ALL BE SET TO NON BLOCKING !!!


!     INPUT FILES:
!     -----------
!     ETOPO1_Ice_g_int.xyz 

!     USER INPUT WHICH SPECIFIES WHAT CONFIGURATION IS REQUIRED
!     input_to_wam_bathymetry...
!     if the file grid_description is also present, then the domain
!     and grid definition will be obtained from that file and not from
!     input_to_wam_bathymetry.... The reason is that grid_description
!     is used to supply information about the number of latitudes and
!     the number of points per latitudes in an attempt to mimic what
!     is done for gaussian grids.

!     REFERENCE LEVELS:
!     reference_levels : USED TO DEFINED AREAS FOR WHICH THE REFERENCE
!                        LEVEL IS NOT MEAN SEA LEVEL. IT ALSO CONTAINS
!                        THE NEW DEPTH TO BE IMPOSED ON THE DEFINED
!                        AREAS PROVIDED THE NEW DEPTH.
!                        THE FILE  MAY CONTAINS UP TO NREF ENTRIES

!     CORRECTION TO WAM POINTS (optional).

!     THEY WERE OBTAINED FROM THE PREVIOUS GRID SETUP
!     correction_to_wam_grid_xxxxx
!     where xxxxx=int(xdella*100) if grid_description is not present
!     else  xxxxx=the gaussian number of latitude - 1 (NGY)

!     OUTPUT FILES:

!     wam_topo_xxxxx
!     where xxxxx=int(xdella*100) if grid_description is not present
!     else  xxxxx=the gaussian spectral truncation ISPECTRUNC

!     AND ALL THE FILES NECESSARY FOR PLOTTING WITH METVIEW
!     IF REQUESTED IN THE USER INPUT !
!     meandepth.dat
!     and
!     obstructions_*.dat
 
!**************************************************************************

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : KCOUSTEP
      USE YOWCOUT  , ONLY : LFDB, LRSTST0
      USE YOWFRED  , ONLY : IFRE1, FRATIO, DELTH, FR, TH
      USE YOWGRIBHD, ONLY : LGRHDIFS ,LNEWLVTP, CEXPVERCLIM, NDATE_TIME_WINDOW_END, CDATECLIM, KPARAM_SUBGRIG
      USE YOWGRIB_HANDLES , ONLY : NGRIB_HANDLE_WAM_I,NGRIB_HANDLE_WAM_S
      USE YOWMAP   , ONLY : IPER, IRGG, IQGAUSS, NGX, NGY, NLONRGG, CLDOMAIN,    &
     &                      AMOWEP , AMOSOP , AMOEAP , AMONOP , XDELLA,  XDELLO, &
     &                      DAMOWEP, DAMOSOP, DAMOEAP, DAMONOP, DXDELLA, DXDELLO
      USE YOWPARAM , ONLY : NANG, NFRE_RED
      USE YOWPCONS , ONLY : PI, ZPI, G, ZMISS
      USE YOWSTAT  , ONLY : MARSTYPE ,YCLASS   ,YEXPVER  ,           &
     &                      NENSFNB  ,NTOTENS  ,NSYSNB   ,NMETNB   , &
     &                      IREFDATE ,ISTREAM  ,NLOCGRB
      USE YOWUBUF  , ONLY : NPROPAGS, NANG_OBS, KTOIS, KTOOBSTRUCT
      USE YOWUNIT  , ONLY : IU08

      USE YOWGRIB  , ONLY : IGRIB_OPEN_FILE, IGRIB_CLOSE_FILE

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "aki.intfb.h"
#include "iniwcst.intfb.h"
#include "iwam_get_unit.intfb.h"
#include "ktoobs.intfb.h"
#include "mfr.intfb.h"
#include "preset_wgrib_template.intfb.h"
#include "wgribenout.intfb.h"

!!    Parameters that can be adapted to tune the mean bathymetry
!!    **********************************************************
!     IF RATIOLAND_THRESHOLD OR MORE LAND OR THE CENTER OF THE GRID BOX IS LAND, THEN AVERAGE OVER LAND POINTS
      REAL(KIND=JWRU) ::  RATIOLAND_THRESHOLD
!     IF THERE IS A RATIO OF SHALLOWER POINTS LARGER THAN RATIOSHALLOW_THRESHOLD
!     THEN THE AVERAGE IS TAKEN OVER THOSE POINTS ALONE.
      REAL(KIND=JWRU) ::  RATIOSHALLOW_THRESHOLD


!!    Parameters that can be adapted to tune the obstruction scheme
!!    *************************************************************
!!    XKDMAX controls the overall impact of subnerged subgrid points in blocking waves (fully or locally)
      REAL(KIND=JWRU), PARAMETER :: XKDMAX=1.5_JWRU

!!    ALPR_DEEP controls the impact of submerged subgrid points in blocking waves as if they were subgrid land points
      REAL(KIND=JWRU), PARAMETER :: ALPR_DEEP=0.025_JWRU

!!    IREINF is used to reinforce land obstructions for small grid spacing (see below as it depends on DXDELLA)
!!    It works by artificially increasing the number of sub grid points detected as land
      INTEGER(KIND=JWIM) :: IREINF
!!    IREINF will also be used to reinforce fully blocking submerged obstructions for small grid spacings,
!!    only if the relative count of subgrid submerged points in a grid box is less than PSHALLOWTRHS
      REAL(KIND=JWRU) :: PSHALLOWTRHS = 0.8_JWRU
!!    and
!!    If the relative count of subgrid land points in a grid box is less than PLANDTRHS then IREINF reinforcement can be applied
      REAL(KIND=JWRU) :: PLANDTRHS = 0.3_JWRU


!!    For a subgrid submerged feature to be blocking, the grid box mean depth need to be at least  XKEXTHRS_DEEP * blocking depth
      REAL(KIND=JWRU), PARAMETER :: XKEXTHRS_DEEP=100.0_JWRU

!!    ISWTHRS is used to compute a depth dependent linear reduction factor for ALPR_DEEP
!!    i.e. ALPR_DEEP is linearly reduced for depth less than ISWTHRS to limit the impact of subgrid points in shallow waters.
      INTEGER(KIND=JWIM), PARAMETER :: ISWTHRS=200

!!    PENHCOR is used to adjust the corner obstructructions for IPROPAGS=2
      REAL(KIND=JWRU), PARAMETER :: PENHCOR = 1.0_JWRU


      INTEGER(KIND=JWIM), PARAMETER :: ILON=21601
      INTEGER(KIND=JWIM), PARAMETER :: ILAT=10801
      INTEGER(KIND=JWIM), PARAMETER :: NREF=500
      INTEGER(KIND=JWIM), PARAMETER :: NDPT=1000

      INTEGER(KIND=JWIM), PARAMETER :: NOOBSTRT=1000

      INTEGER(KIND=JWIM) :: IU01, IU06, IUGRD, IUNIT
      INTEGER(KIND=JWIM) :: I, J, IJ, K, KSN, KNS, M, IANG, IP
      INTEGER(KIND=JWIM) :: ISPECTRUNC
      INTEGER(KIND=JWIM) :: NLANDCENTREPM, NLANDCENTREMAX, NLANDCENTRE, NIOBSLAT
      INTEGER(KIND=JWIM) :: NSEA, NLAND, NSEASH
      INTEGER(KIND=JWIM) :: ILONL, ILONR, ILATB, ILATT
      INTEGER(KIND=JWIM) :: NTOT, ICOUNT, IC, IR
      INTEGER(KIND=JWIM) :: NREFERENCE
      INTEGER(KIND=JWIM) :: IX, IXLP, NJM, NJP, NIM, NIP, IH
      INTEGER(KIND=JWIM) :: II, JJ, IK, NPTS, IDPT
      INTEGER(KIND=JWIM) :: ITEMPEW
      INTEGER(KIND=JWIM) :: IS, KT, KB, IOBSRT
      INTEGER(KIND=JWIM) :: NOBSTRCT, NIOBSLON, NBLOCKLAND, NTOTPTS
      INTEGER(KIND=JWIM) :: INVRES
      INTEGER(KIND=JWIM) :: ITABLE, IPARAM, IZLEV, IFCST, ITEST, ITMIN, ITMAX

      INTEGER(KIND=JWIM) :: IDUM(15)
      INTEGER(KIND=JWIM), DIMENSION(NREF) :: LEVEL, NDEPTH
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:) :: ITHRSHOLD
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:) :: IBLOCKDPT
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:) :: IDEPTH
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:) :: ILSM
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:,:) :: IOBSLAT, IOBSLON
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:,:) :: IOBSCOR
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:,:) :: IOBSRLAT, IOBSRLON

      REAL(KIND=JWRU) :: DRAD, RESOL
      REAL(KIND=JWRU), ALLOCATABLE, DIMENSION(:) :: ZDELLO

      REAL(KIND=JWRU), PARAMETER :: RMIN_DEPTH = -0.3_JWRU
      REAL(KIND=JWRU), PARAMETER :: RMIN_DEPTH_SMOOTH = RMIN_DEPTH-0.01_JWRU
      REAL(KIND=JWRU) :: X60
      REAL(KIND=JWRU) :: ALONL, ALONR, ALATB, ALATT, XLON
      REAL(KIND=JWRU) :: REXCLTHRSHOLD
      REAL(KIND=JWRU) :: XLO, XLA, XI, YJ
      REAL(KIND=JWRU) :: SEA, XLAND, SEASH
      REAL(KIND=JWRU) :: XX 
      REAL(KIND=JWRU) :: STEPT, STEPB, XLATT, XLATB, XLONL, XLONR
      REAL(KIND=JWRU) :: STEPLAT, STEPLON
      REAL(KIND=JWRU) :: RR, XKEXTHRS, ALPR
      REAL(KIND=JWRU), DIMENSION(ILON) :: ALON
      REAL(KIND=JWRU), DIMENSION(ILAT) :: ALAT
      REAL(KIND=JWRU), DIMENSION(NDPT) :: XK 
      REAL(KIND=JWRU), DIMENSION(NREF) :: XINF, XSUP, YINF, YSUP
      REAL(KIND=JWRU), ALLOCATABLE, DIMENSION(:) :: COSPH
      REAL(KIND=JWRU), ALLOCATABLE, DIMENSION(:) :: XLAT
      REAL(KIND=JWRU), ALLOCATABLE, DIMENSION(:,:) :: WAMDEPTH
      REAL(KIND=JWRU), ALLOCATABLE, DIMENSION(:,:) :: PERCENTLAND, PERCENTSHALLOW


      REAL(KIND=JWRB) :: PRPLRADI, ZCONV, FR1, OMEGA, DEPTH
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:) :: FIELD  !!! because it will be used for grib encoding,
                                                             !!! it needs to be defined from North to South
      CHARACTER(LEN=  1) :: C1
      CHARACTER(LEN=  2) :: CFR
      CHARACTER(LEN=  5) :: CWAMRESOL
      CHARACTER(LEN=  5) :: CX
      CHARACTER(LEN= 11) :: FORMAT
      CHARACTER(LEN= 14) :: CDATE
      CHARACTER(LEN= 32) :: FILENM
      CHARACTER(LEN= 72) :: LOCATION(NREF)
      CHARACTER(LEN=144) :: CLINE, FILENAME

      LOGICAL :: LORIGINAL, LLPRINT
      LOGICAL :: LLAND, LREALLAND, L1ST, LNSW
      LOGICAL :: LLGRID
      LOGICAL :: LLGRIBIN !! funtionality not yet fully coded
      LOGICAL :: LLOBSTROUT, LLGRIBOUT
      LOGICAL :: LLOBSTRON
      LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: LLEXCLTHRSHOLD
      LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: LLSM

!----------------------------------------------------------------------

      PRPLRADI=1.0_JWRB
      CALL INIWCST(PRPLRADI)

      DRAD=4.0_JWRU*ATAN(1.0_JWRU)/180.0_JWRU 

!     ETOPO1 RESOLUTION
      INVRES=60
      X60=60.0_JWRU
      RESOL=1.0_JWRU/REAL(INVRES,JWRU)

      IU01=1
      IU06=6

      NANG=NANG_OBS !!! PSEUDO value used to encode obstruction coefficients in grib

      LFDB = .FALSE.

      LLOBSTRON = .TRUE.

!     READ INPUT SELECTION
!     --------------------
      READ(5,*) CLINE
      READ(5,*) CLINE
      READ(5,*) LLGRIBIN
      READ(5,*) CLINE
      READ(5,*) LLOBSTROUT
      READ(5,*) CLINE
      READ(5,*) LLGRIBOUT
      READ(5,*) CLINE
      READ(5,*) DXDELLA
      READ(5,*) CLINE
      READ(5,*) DAMOSOP, DAMONOP, DAMOWEP, DAMOEAP
      READ(5,*) CLINE
      READ(5,*) IPER
      READ(5,*) CLINE
      READ(5,*) IRGG
      READ(5,*) CLINE
      READ(5,*) FR1 
      READ(5,*) CLINE
      READ(5,*) NFRE_RED, IFRE1
      READ(5,*) CLINE
      READ(5,*) LLPRINT
      READ(5,*) CLINE
      READ(5,*) ALONL, ALONR, ALATB, ALATT
      READ(5,*) CLINE
      READ(5,*) LORIGINAL 

!     CHECK IF FILE grid_description IS PRESENT
!     IF IT IS THERE IT WILL SUPERSEDE THE OTHER INPUT

      FILENAME='grid_description'
!!!!!!!!!! grid_description is also read in uiprep  !!!!!!!!!!
      INQUIRE(FILE=FILENAME,EXIST=LLGRID)
      IF (LLGRID) THEN
        IUGRD=IWAM_GET_UNIT(IU06,FILENAME,'S','F',0,'READWRITE')
        OPEN(IUGRD,FILE=FILENAME,STATUS='OLD', FORM='FORMATTED')
        READ (IUGRD,*) ISPECTRUNC
        READ (IUGRD,*) DAMONOP
        READ (IUGRD,*) DAMOSOP
        READ (IUGRD,*) DAMOWEP
        READ (IUGRD,*) DAMOEAP
        READ (IUGRD,*) IPER
        READ (IUGRD,*) IRGG
        READ (IUGRD,*) NGY

        WRITE(IU06,*) "DAMONOP = ",DAMONOP
        WRITE(IU06,*) "DAMOSOP = ",DAMOSOP
        WRITE(IU06,*) "DAMOWEP = ",DAMOWEP
        WRITE(IU06,*) "DAMOEAP = ",DAMOEAP

        IF (ISPECTRUNC > 0) THEN
          IQGAUSS=1
        ELSE
          IQGAUSS=0
        ENDIF

        DXDELLA = (DAMONOP-DAMOSOP)/REAL(NGY-1,JWRU)
        ALLOCATE(NLONRGG(NGY))

        NGX = 0
        DO K=1,NGY
          KSN=NGY-K+1
          READ(IUGRD,*) NLONRGG(KSN)
          NGX = MAX(NGX,NLONRGG(KSN))
        ENDDO

        IF (IPER == 1) THEN
          DXDELLO  = 360._JWRU/REAL(NGX,JWRU)
          DAMOEAP = DAMOWEP + 360._JWRU - DXDELLO
          CLDOMAIN = 'g'
        ELSE
          DXDELLO = (DAMOEAP-DAMOWEP)/(NGX-1)
          CLDOMAIN = 'm'
        ENDIF

        CLOSE(IUGRD)

      ELSE
        DXDELLO=DXDELLA     
        NGX=NINT((DAMOEAP-DAMOWEP)/DXDELLO)+1
        NGY=NINT((DAMONOP-DAMOSOP)/DXDELLA)+1

        WRITE(IU06,*) "NGX = ",NGX
        WRITE(IU06,*) "NGY = ",NGY

        ALLOCATE(NLONRGG(NGY))

        IF (IPER == 1) THEN
          CLDOMAIN = 'g'
        ELSE
          CLDOMAIN = 'm'
        ENDIF
      ENDIF

      AMONOP = REAL(DAMONOP,JWRB)
      AMOSOP = REAL(DAMOSOP,JWRB)
      AMOWEP = REAL(DAMOWEP,JWRB)
      AMOEAP = REAL(DAMOEAP,JWRB)
      XDELLA = REAL(DXDELLA,JWRB)
      XDELLO = REAL(DXDELLO,JWRB)

    
      IF ( DXDELLA < 0.125_JWRU) THEN
        RATIOLAND_THRESHOLD = 0.5_JWRU
        RATIOSHALLOW_THRESHOLD = 1.0_JWRU
      ELSE
        RATIOLAND_THRESHOLD = 0.6_JWRU
        RATIOSHALLOW_THRESHOLD = 0.3_JWRU
      ENDIF
     
  
      NLANDCENTREPM=(NINT(0.2_JWRU*DXDELLA*INVRES)-1)/2
      NLANDCENTREPM=MAX(NLANDCENTREPM,1)
      NLANDCENTREMAX=(2*NLANDCENTREPM+1)**2


      ALLOCATE(ZDELLO(NGY))
      ALLOCATE(COSPH(NGY))
      ALLOCATE(XLAT(0:NGY+1))

      DO K=0,NGY+1
        XLAT(K) = (DAMOSOP + REAL(K-1,JWRU)*DXDELLA)
      ENDDO

      DO K=1,NGY
!       !!! from south to north !!!!
        COSPH(K) = COS(XLAT(K)*DRAD)
        IF (.NOT.LLGRID) THEN
          IF (IRGG == 1) THEN
!            The silly division by cos(x60*DRAD) is an attempt at making sure
!            that exactly 0.5 is used for cosine of 60 degrees.
             NLONRGG(K)=MAX(NINT(NGX*(COS(XLAT(K)*DRAD)/(2._JWRU*COS(X60*DRAD)))),2)
            IF (MOD(NLONRGG(K),2) == 1) NLONRGG(K) = NLONRGG(K)+1
          ELSE
            NLONRGG(K) = NGX
          ENDIF
          WRITE(IU06,*) 'POINTS PER LATITUDES: ',K, XLAT(K), NLONRGG(K)
        ENDIF

        IF (IPER == 1) THEN
          ZDELLO(K)  = 360._JWRU/REAL(NLONRGG(K),JWRU)
        ELSE
          ZDELLO(K)  = (DAMOEAP-DAMOWEP)/REAL(NLONRGG(K)-1,JWRU)
        ENDIF
      ENDDO

      ALLOCATE(FR(NFRE_RED))

      CALL MFR(NFRE_RED, IFRE1, FR1, FRATIO, FR)


      MARSTYPE = 'an'
      YCLASS   = 'od'
      YEXPVER = CEXPVERCLIM
      NENSFNB = 0  
      NTOTENS = 0
      ISTREAM = 1045 !! if changed to an ifs stream also change LNEWLVTP
      NLOCGRB = 1
      NSYSNB  = -1
      NMETNB  = -1
      IREFDATE = 0
      LGRHDIFS =.FALSE.
      LNEWLVTP =.FALSE.
      NDATE_TIME_WINDOW_END = 0
      KCOUSTEP = .FALSE.
      LRSTST0 = .FALSE.


      IF ( LLOBSTROUT .AND. LLGRIBOUT ) THEN

        WRITE(IU06,*) ''
        WRITE(IU06,*) 'OUTPUT IN GRIB '
        WRITE(IU06,*) ''

!       PREPARE OUTPUT

!       Use direction dimension to save obstruction coefficients
        DELTH = ZPI/REAL(NANG,JWRB)
        ALLOCATE(TH(NANG))
        DO IANG  = 1, NANG
          TH(IANG) = REAL(IANG-1,JWRB)*DELTH
        ENDDO 

!       FOR INTEGRATED PARAMETERS
        CALL PRESET_WGRIB_TEMPLATE("I",NGRIB_HANDLE_WAM_I,NGRIBV=2,LLCREATE=.true.,NBITSPERVALUE=24)
!       FOR SPECTRA
        CALL PRESET_WGRIB_TEMPLATE("S",NGRIB_HANDLE_WAM_S,NGRIBV=2,LLCREATE=.true.,NBITSPERVALUE=12)
        
        DO IP = 0, NPROPAGS 
          WRITE(C1,'(I1)') IP
          FILENAME='wam_grib_subgrid_'//C1
          CALL IGRIB_OPEN_FILE(IU08(IP),FILENAME,'w')
        ENDDO
      ENDIF

!     DATASET:
!     --------
 
      IF (ALONL < -180._JWRU .OR. ALONR > 180._JWRU) THEN
        WRITE(*,*) ' LONGITUDE SPECIFICATION ERROR +- 180'
        WRITE(*,*) ' ALONL, ALONR : ',ALONL,ALONR
        CALL ABORT1
      ENDIF
      IF (ALATT > 90.0_JWRU .OR. ALATB < -90.0_JWRU) THEN
        WRITE(*,*) ' LATITUDE SPECIFICATION ERROR +- 90'
        WRITE(*,*) ' ALATT, ALATB : ',ALATT,ALATB
        CALL ABORT1
      ENDIF
     
      ILONL = NINT((ALONL + 180._JWRU)*INVRES) + 1
      ILONR = NINT((ALONR + 180._JWRU)*INVRES) + 1
      ILATB = NINT((90.0_JWRU- ALATB)*INVRES) + 1
      ILATB = MAX(1,MIN(ILATB,ILAT))
      ILATT = NINT((90.0_JWRU- ALATT)*INVRES) + 1
      ILATT = MAX(1,MIN(ILATT,ILAT))
      IF (ILONR == ILON+1) ILONR=ILON
      IF (ILATB == ILAT+1) ILATB=ILAT

      ALLOCATE(IDEPTH(ILON,ILAT))
      ALLOCATE(WAMDEPTH(NGX,NGY))
      ALLOCATE(PERCENTLAND(NGX,NGY))
      ALLOCATE(PERCENTSHALLOW(NGX,NGY))

      OPEN(20,FILE='ETOPO1_Ice_g_int.xyz')

!     READ INPUT DATA
      NTOT=ILON*ILAT
      DO J = 1,ILAT
        DO I = 1,ILON
          READ(20,*,END=101) ALON(I),ALAT(J),IDEPTH(I,J)
        ENDDO
      ENDDO
101   CONTINUE

!     OUTPUT ORIGINAL DATA SET ON SUB AREA

      IF (LORIGINAL) THEN
        OPEN(11,file='depth.dat')
        WRITE(11,'(a4)') '#GEO'
        WRITE(11,'(a11)') '#FORMAT LLV'
        WRITE(11,'(a5)') '#DATA'
        write(*,*) 'output of original data set for indices:'
        write(*,*) ILATT,ILATB
        write(*,*) ILONL,ILONR
        DO J=ILATT,ILATB
          DO I=ILONL,ILONR
            IF (IDEPTH(I,J) >=  -300 .AND. IDEPTH(I,J) <= 2000 ) THEN
            WRITE(11,'(2(1X,F8.3),1X,I4)')ALON(I),ALAT(J),IDEPTH(I,J)
            ENDIF
          ENDDO
        ENDDO
      ENDIF


!     COMPUTE THE MEAN BATHYMETRY
!     ---------------------------
!     READ REFERENCE LEVEL DEFINITION FILE ( DEFINING REGIONS WHERE
!     MEAN SEA LEVEL SHOULD NOT BE THE REFERENCE LEVEL)

      OPEN(15,FILE='reference_levels',STATUS='OLD')

      DO IR=1,NREF
         READ(15,*,END=1000,ERR=1000)                                              &
     &        XINF(IR),YINF(IR),XSUP(IR),YSUP(IR),LEVEL(IR),NDEPTH(IR),LOCATION(IR)
      ENDDO


1000  NREFERENCE=IR-1
      WRITE(IU06,*) 'READ ',NREFERENCE,' NEW REFERENCE LEVELS'

      DO IR=1,NREFERENCE
        WRITE(IU06,*)                                                              &
     &        XINF(IR),YINF(IR),XSUP(IR),YSUP(IR),LEVEL(IR),NDEPTH(IR),LOCATION(IR)
        DO J=1,ILAT
          YJ=ALAT(J)
          IF (YJ >= YINF(IR) .AND. YJ <= YSUP(IR)) THEN
            DO I=1,ILON
              XI=ALON(I)
              IF (XI >= XINF(IR) .AND. XI <= XSUP(IR)) THEN
                 IF (IDEPTH(I,J) <= LEVEL(IR)) THEN
                   IF (NDEPTH(IR) /= 0) THEN
                     IDEPTH(I,J)=NDEPTH(IR)
                   ELSE
                     IDEPTH(I,J)=IDEPTH(I,J)-LEVEL(IR)
                   ENDIF
                 ELSE
                   IDEPTH(I,J)=IDEPTH(I,J)-LEVEL(IR)
                 ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      CALL FLUSH(IU06)

!     COMPUTE MEAN DEPTH 

      DO K=1,NGY
         DO IX=1,NGX
           WAMDEPTH(IX,K)=0.0_JWRU
         ENDDO
      ENDDO


      NJM=INT(0.5_JWRU*DXDELLA*INVRES)
      NJP=NINT(0.5_JWRU*DXDELLA*INVRES)
      IF ( DXDELLA < 0.125_JWRU) THEN
!       Allow a bit of extra smoothing
        NJM=NJM+1
        NJP=NJP+1
      ENDIF

      DO K=1,NGY
!        WE ASSUME THAT WAMGRID IS ALWAYS WITHIN ETOPO1
!        DETERMINE CLOSEST ETOPO1 J INDEX TO WAM POINT
         DO J=ILAT-1,1,-1
           IF (ALAT(J+1) < XLAT(K) .AND. XLAT(K) <= ALAT(J) ) EXIT
         ENDDO
         J=MIN(MAX(J,1),ILAT)

         DO IX=1,NLONRGG(K)

!          ETOPO1 STARTS AT -180
           XLON=DAMOWEP + REAL(IX-1,JWRU)*ZDELLO(K)
           IF (XLON > 180._JWRU) THEN
             XLON=XLON-360._JWRU
           ENDIF

!          DETERMINE CLOSEST ETOPO1 I INDEX TO WAM POINT
           DO I=1,ILON-1
             IF (ALON(I) <= XLON .AND. XLON < ALON(I+1) ) EXIT
           ENDDO

           NIM=INT(0.5_JWRU*ZDELLO(K)*REAL(INVRES,JWRU))
           NIP=NINT(0.5_JWRU*ZDELLO(K)*REAL(INVRES,JWRU))
           IF ( DXDELLA < 0.125_JWRU) THEN
!            Allow a bit of extra smoothing
             NIM=NIM+1
             NIP=NIP+1
           ENDIF

           NSEA=0
           SEA=0._JWRU
           NLAND=0
           XLAND=0._JWRU
           NSEASH=0
           SEASH=0._JWRU

!          AVERAGE OVER LAND AND SEA SEPARATELY
!          AROUND POINT I,J
           DO JJ=J-NJM,J+NJP
             IF (JJ >= 1 .AND. JJ <= ILAT) THEN
               DO II=I-NIM,I+NIP
                 IK=II
                 IF (II < 1) IK=ILON+II
                 IF (II > ILON) IK=II-ILON
                 IF (IDEPTH(IK,JJ) <= 0) THEN
                   NSEA=NSEA+1
!                  IN WAM 999m IS THE MAXIMUM DEPTH
                   SEA=SEA+MAX(-999,IDEPTH(IK,JJ))

!                  FIND SHALLOWER AREAS
                   IF (IDEPTH(IK,JJ) > -500) THEN
                      NSEASH=NSEASH+1
                      SEASH=SEASH+IDEPTH(IK,JJ)
                   ENDIF

                 ELSE
                   NLAND=NLAND+1
                   XLAND=XLAND+IDEPTH(IK,JJ)
                 ENDIF
               ENDDO
             ENDIF
           ENDDO

!          SEARCH FOR LAND AT THE CENTER OF THE GRID BOX
           NLANDCENTRE=0
           DO JJ=J-NLANDCENTREPM,J+NLANDCENTREPM
             IF (JJ >= 1 .AND. JJ <= ILAT) THEN
               DO II=I-NLANDCENTREPM,I+NLANDCENTREPM
                 IK=II
                 IF (II < 1) IK=ILON+II
                 IF (II > ILON) IK=II-ILON
                 IF (IDEPTH(IK,JJ) > 0) THEN
                   NLANDCENTRE=NLANDCENTRE+1
                 ENDIF
               ENDDO
             ENDIF
           ENDDO

!          IF RATIOLAND_THRESHOLD OR MORE LAND OR THE CENTER OF THE GRID BOX IS LAND, THEN AVERAGE OVER LAND POINTS
!          ELSE AVERAGE OVER SEA POINTS
           PERCENTLAND(IX,K)=REAL(NLAND,JWRU)/REAL((NLAND+NSEA),JWRU)
           IF (PERCENTLAND(IX,K) >  RATIOLAND_THRESHOLD .OR. NLANDCENTRE >= NLANDCENTREMAX ) THEN
             WAMDEPTH(IX,K)=XLAND/NLAND
           ELSE
!            IF THERE IS A PERCENTAGE OF SHALLOWER POINTS THEN
!            THE AVERAGE IS TAKEN OVER THOSE POINTS ALONE.
             PERCENTSHALLOW(IX,K)=REAL(NSEASH,JWRU)/REAL(NSEA,JWRU)
             IF (PERCENTSHALLOW(IX,K) >= RATIOSHALLOW_THRESHOLD) THEN
               WAMDEPTH(IX,K)=SEASH/NSEASH
               IF (PERCENTLAND(IX,K) < 0.10_JWRU) THEN
!                IF MOSTLY SEA THEN IT SHOULD BE SEA AND NOT 0 
                 WAMDEPTH(IX,K)=MIN(WAMDEPTH(IX,K),-1.0_JWRU)
               ENDIF
             ELSE 
               WAMDEPTH(IX,K)=SEA/NSEA
             ENDIF
           ENDIF

         ENDDO
      ENDDO

!     RESET TO LAND SEA POINTS THAT ARE NOT DEEP ENOUGH
      NPTS=0
      DO K=1,NGY
        DO IX=1,NLONRGG(K)
          IF ( WAMDEPTH(IX,K) > RMIN_DEPTH .AND. WAMDEPTH(IX,K) < 0.0_JWRU ) THEN
            NPTS=NPTS+1
            WAMDEPTH(IX,K) = -WAMDEPTH(IX,K)
          ENDIF
          IF ( WAMDEPTH(IX,K) > RMIN_DEPTH_SMOOTH .AND. WAMDEPTH(IX,K) < RMIN_DEPTH ) THEN
             WAMDEPTH(IX,K) = RMIN_DEPTH
          ENDIF
        ENDDO
      ENDDO
      WRITE(IU06,*) 'NUMBER OF POINTS RESET TO LAND AS NOT DEEP ENOUGH = ',NPTS
      WRITE(IU06,*) ' '


      NPTS=0
      DO K=1,NGY
        DO IX=1,NLONRGG(K)
          IF (WAMDEPTH(IX,K) < 0.0_JWRU) NPTS=NPTS+1
        ENDDO
      ENDDO
      WRITE(IU06,*) 'TOTAL NUMBER OF SEA POINTS BEFORE CORRECTIONS = ',NPTS


!     INTRODUCE CORRECTION TO WAM POINTS.
!     THEY WERE OBTAINED FROM THE PREVIOUS GRID SETUP
!     !!! SEA POINT DEPTHS ARE ALREADY POSITIVE !!!
!     IT's a bit obsolete !!!

      WRITE(IU06,*) 'CORRECTION TO WAM GRID'

      IF (LLGRID) THEN
        WRITE(CWAMRESOL,'(I5.5)') NGY-1
      ELSE
        WRITE(CWAMRESOL,'(I5.5)') INT(1000*DXDELLA)
      ENDIF
      FILENAME='correction_to_wam_grid_'//CWAMRESOL

      OPEN(35,FILE=FILENAME)

      DO WHILE(.TRUE.)
         READ(35,*,END=111,ERR=111) XLO,XLA,IX,K,IDPT
         IDPT=-IDPT
         XLON=DAMOWEP + REAL(IX-1,JWRU)*ZDELLO(K)
         IF (ABS(XLON-XLO) > REAL(ZDELLO(K),JWRU) .OR. ABS(XLAT(K)-XLA) > REAL(DXDELLA,JWRU) ) THEN
           WRITE(*,*) 'PROBLEM !!!!'
           WRITE(*,*) 'THE CORRECTION TO WAM GRID IS NOT A WAM POINT'
           WRITE(*,*) XLO,XLA,IDPT 
           WRITE(*,*) XLON,XLAT(K),WAMDEPTH(IX,K)
           WRITE(*,*) ''
         ELSE
           IF (WAMDEPTH(IX,K) >= 0.0_JWRU .AND. IDPT < 0) THEN
             NPTS=NPTS+1
           ELSEIF (WAMDEPTH(IX,K) < 0.0_JWRU .AND. IDPT >= 0) THEN
             NPTS=NPTS-1
           ENDIF  
           WAMDEPTH(IX,K)=REAL(IDPT,JWRU)
         ENDIF
      ENDDO
111   CONTINUE


!     SET MIN AND MAX TO +-999
      DO K=1,NGY
         DO IX=1,NLONRGG(K)
           WAMDEPTH(IX,K)=MIN(999._JWRU,MAX(-999._JWRU,WAMDEPTH(IX,K)))
        ENDDO
      ENDDO

      WRITE(IU06,*) 'TOTAL NUMBER OF SEA POINTS = ',NPTS
 
!     OUTPUT THE MEAN DATA SET ON SUB AREA
!     OMIT LAND POINTS
!     ------------------------------------
      IF (LLPRINT) THEN
        OPEN(12,file='meandepth.dat')
        WRITE(12,'(A4)') '#GEO'
        WRITE(12,'(A11)') '#FORMAT LLV'
        WRITE(12,'(A5)') '#DATA'
        DO K=1,NGY
           DO IX=1,NLONRGG(K)
             XLON=DAMOWEP + REAL(IX-1,JWRU)*ZDELLO(K)
             IF (XLON > 180._JWRU) THEN
               XLON=XLON-360._JWRU
             ENDIF
            IF (ALATB <= XLAT(K) .AND. XLAT(K) <= ALATT .AND.           &
     &          ALONL <= XLON .AND. XLON <= ALONR ) THEN
               IF (WAMDEPTH(IX,K) < 0.0_JWRU .AND. WAMDEPTH(IX,K) > -999._JWRU ) THEN
                  WRITE(12,'(3(1X,F8.3))') XLON,XLAT(K),WAMDEPTH(IX,K)
               ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF


!     OUTPUT THE MEAN DATA AS THE INPUT FILE FOR WAM
!     ----------------------------------------------
      IF (LLGRID) THEN
        WRITE(CWAMRESOL,'(I5.5)') ISPECTRUNC
      ELSE
        WRITE(CWAMRESOL,'(I5.5)') INT(1000*DXDELLA)
      ENDIF
      FILENAME='wam_topo_'//CWAMRESOL

      OPEN(IU01,FILE=FILENAME,FORM='FORMATTED')
      IF (LLGRID) THEN
        WRITE(IU01,'(A)') 'WAM BATHYMETRY'
        WRITE(IU01,'(6F13.8)') DXDELLA, DXDELLO, DAMOSOP, DAMONOP, DAMOWEP, DAMOEAP
      ELSE
        WRITE(IU01,'(6F10.5)') DXDELLA, DXDELLO, DAMOSOP, DAMONOP, DAMOWEP, DAMOEAP
      ENDIF

      CX='     '
      FORMAT='          '
      IF (LLGRID) THEN
        DO K=1,NGY
           WRITE(IU01,'(I5.5)') NLONRGG(K) 
        ENDDO
        WRITE(CX,'(I5.5)') NLONRGG(1)
        FORMAT='('//CX//'F9.2)'
        DO K=1,NGY
          DO IS = 1,NLONRGG(K),NLONRGG(1)
            WRITE(IU01,FORMAT) (WAMDEPTH(IX,K),IX=IS,MIN(IS+NLONRGG(1)-1,NLONRGG(K))) 
          ENDDO
        ENDDO
      ELSE
        DO K=1,NGY
           WRITE(IU01,'(I4.4)') NLONRGG(K) 
        ENDDO
        DO K=1,NGY
          WRITE(CX,'(I4.4)') NLONRGG(K) 
          FORMAT='('//CX//'I4)'
          WRITE(IU01,FORMAT) (NINT(WAMDEPTH(IX,K)),IX=1,NLONRGG(K)) 
        ENDDO
      ENDIF




IF ( LLOBSTROUT ) THEN

!     CREATE OBSTRUCTIONS COEFFICIENTS:
!     --------------------------------

      ALLOCATE(IOBSLAT(NGX,NGY,2))
      ALLOCATE(IOBSLON(NGX,NGY,2))
      ALLOCATE(IOBSCOR(NGX,NGY,4))
      ALLOCATE(IOBSRLAT(NGX,NGY,2))
      ALLOCATE(IOBSRLON(NGX,NGY,2))
      ALLOCATE(FIELD(NGX,NGY))

      WRITE(1,'(I4)') NFRE_RED

      IF (DXDELLA <= 0.125_JWRU) THEN
        IREINF=4
      ELSEIF (DXDELLA <= 0.5_JWRU) THEN
        IREINF=2
      ELSE
        IREINF=1
      ENDIF

      ALLOCATE(ITHRSHOLD(NGX,NGY))
      ALLOCATE(IBLOCKDPT(NGX,NGY))
      ALLOCATE(LLEXCLTHRSHOLD(NGX,NGY))
      ALLOCATE(LLSM(NGX,NGY))
      ALLOCATE(ILSM(NGX,NGY))

      IOBSLAT(:,:,:)  = NOOBSTRT
      IOBSLON(:,:,:)  = NOOBSTRT 
      IOBSCOR(:,:,:)  = NOOBSTRT 
      IOBSRLAT(:,:,:) = NOOBSTRT 
      IOBSRLON(:,:,:) = NOOBSTRT
      FIELD(:,:) = ZMISS 
      LLSM(:,:) = .FALSE.
      ILSM(:,:) = 0 

      ZCONV = 1.0_JWRU/REAL(NOOBSTRT,JWRB)



!     LOOP OVER ALL FREQUENCIES
!     -------------------------

      DO M=1,NFRE_RED

!!!!    COMPUTE THE OBSTRUCTIONS ONLY WHEN IT IS MEANINGFUL
        IF (DXDELLA <= 2.0_JWRU*RESOL ) THEN
          IF (M == 1) THEN
            WRITE(*,*) ''
            WRITE(*,*) '*********************************************'
            WRITE(*,*) 'THE REQUESTED RESOLUTION IS SMALL ENOUGH WITH'
            WRITE(*,*) 'RESPECT TO THE INPUT BATHYMETRY DATA !'
            WRITE(*,*) 'NO OBSTRUCTIONS WILL BE COMPUTED !'
            WRITE(*,*) '*********************************************'
            WRITE(*,*) ''
          ENDIF
          LLOBSTRON = .FALSE.
        ELSE 

!       COMPUTE WAVE NUMBERS
        OMEGA=ZPI*FR(M)
        DO IDPT=1,NDPT
          DEPTH=REAL(IDPT,JWRB)
          XK(IDPT)=REAL(AKI(OMEGA,DEPTH),JWRU)
        ENDDO


!       COMPUTE THE THRESHOLD AT WHICH THE WAVES ARE PARTIALLY OBSTRUCTED BY THE BOTTOM (ITHRSHOLD),
!       EXCEPT IF WAM DEPTH OF THE SAME ORDER OF MAGITUDE (.NOT. LLEXCLTHRSHOLD) .
!       ALSO COMPUTE THE DEPTH THAT IS CONSIDERED TO BE FULLY BLOCKING AS IF IT WAS LAND (IBLOCKDPT).
!       ALSO SET THE LAND SEA MASK
        DO K=1,NGY
          DO IX=1,NLONRGG(K)
            IF (WAMDEPTH(IX,K) < 0.0_JWRU) THEN
              XX=XKDMAX/XK(MAX(MIN(-NINT(WAMDEPTH(IX,K)),NDPT),1))
              ITHRSHOLD(IX,K)=NINT(-XX)
              RR=MAX(REAL((ISWTHRS-ABS(NINT(WAMDEPTH(IX,K)))),JWRU)/ISWTHRS,0.0_JWRU)
              XKEXTHRS=XKEXTHRS_DEEP*(1.0_JWRU+RR)
              ALPR=MAX(ALPR_DEEP*(1.0_JWRU-RR),0.0_JWRU)
              REXCLTHRSHOLD=MAX(XKEXTHRS*ITHRSHOLD(IX,K),-998._JWRU)
              LLEXCLTHRSHOLD(IX,K)=(WAMDEPTH(IX,K) < REXCLTHRSHOLD)
              IBLOCKDPT(IX,K)=INT(-ALPR*XX)
              LLSM(IX,K) = .TRUE.
              ILSM(IX,K) = 1
            ENDIF
          ENDDO
        ENDDO


!       NORTH-SOUTH OBSTRUCTIONS
!       -----------------------
!       LOOP OVER ADEVECTION DIRECTIONS
!       IS=1 is for the south-north advection
!       IS=2 is for the north-south advection
        WRITE(IU06,*) 'CREATE NORTH-SOUTH OBSTRUCTIONS '
        DO IS=1,2
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) &
!$OMP& PRIVATE(K,KT,KB,STEPT,STEPB,XLATT,XLATB,ILATT,ILATB,IX) &
!$OMP& PRIVATE(XLONL,XLONR,ILONL,ILONR,NOBSTRCT,NBLOCKLAND) &
!$OMP& PRIVATE(I,NIOBSLON,LLAND,LREALLAND,J,LNSW,L1ST) &
!$OMP& PRIVATE(NTOTPTS)
!         LOOP OVER MODEL LATITUDES
          DO K=1,NGY
            IF (IS == 1) THEN
              KT=K
              KB=K-1
              STEPT=-RESOL
              STEPB=0._JWRU
            ELSE
              KT=K+1
              KB=K
              STEPT=0._JWRU
              STEPB=RESOL
            ENDIF
!           LATIDUNAL INDEX OF THE OF THE SUBGRID POINTS THAT ARE INSIDE THE MODEL GRID BOX:  
            XLATT=XLAT(KT)+STEPT
            XLATB=XLAT(KB)+STEPB
            ILATT = NINT((90.0_JWRU- XLATT)*INVRES) + 1
            ILATT = MAX(1,MIN(ILATT,ILAT))
            ILATB = NINT((90.0_JWRU- XLATB)*INVRES) + 1
            ILATB = MAX(1,MIN(ILATB,ILAT))
            IF (ILATB == ILAT+1) ILATB=ILAT

!           LOOP OVER ALL MODEL POINTS FOR A GIVEN LATITUDE
            DO IX=1,NLONRGG(K)
              IF (LLSM(IX,K)) THEN
!               SEA POINT GRID BOX LATITUNAL EXTEND :
                XLONL=DAMOWEP + (REAL(IX-1,JWRU)-0.5_JWRU)*ZDELLO(K)
                IF (XLONL > 180._JWRU) THEN
                  XLONL=XLONL-360._JWRU
                ENDIF
                XLONR=DAMOWEP + (REAL(IX-1,JWRU)+0.5_JWRU)*ZDELLO(K)
                IF (XLONR > 180._JWRU) THEN
                  XLONR=XLONR-360._JWRU
                ENDIF

!               LONGITUDINAL INDEX OF THE OF THE SUBGRID POINTS THAT ARE INSIDE THE MODEL GRID BOX:  
                IF ( DXDELLA < 0.125_JWRU) THEN
                  ILONL = INT((XLONL + 180._JWRU)*INVRES) + 1
                  ILONR = INT((XLONR + 180._JWRU)*INVRES) + 2
                ELSE
!                 It was decided to not correct the double counting for low resolution (should be removed in future)
                  ILONL = NINT((XLONL + 180._JWRU)*INVRES) + 1
                  ILONR = NINT((XLONR + 180._JWRU)*INVRES) + 1
                ENDIF

!               COMPUTE THE OBSTRUCTIONS:
!               TALLY THE NUMBER OF SUB GRID POINTS THAT ARE POTENTIALLY BLOCKING WAVE PROPAGATION (NOBSTRCT)
                NOBSTRCT=0

!               AWAY FROM THE DATELINE
                IF (ILONL <= ILONR) THEN
                  NBLOCKLAND=0                  
!                 LOOP OVER SUBGRID LONGITUDE LINE:
                  DO I=ILONL,ILONR
                    NIOBSLON=0
                    LLAND=.FALSE.
                    LREALLAND=.FALSE.
!                   SCAN EACH SUBGRID LATTUDE:
                    DO J=ILATT,ILATB
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K) ) THEN
!                       IF THE LONGITUDE LINE CONTAINS ACTUAL LAND (> 0) THEN THE FULL LONGITUDE WILL BLOCK ONLY IF
!                       THERE IS A SWITCH BACK TO DEPTH < IBLOCKDPT OR VICE VERSA (see below).
!                       THIS IS TO AVOID CREATING FULL OBSTRUCTION WHEN APPROACHING THE COASTLINE.
!                       ELSE IF THE LINE CONTAINS PSEUDO LAND AS DEFINED AS ANYTHING ABOVE IBLOCKDPT(IX,K) (IBLOCKDPT is negative)
!                       THEN THE FULL LINE WILL BLOCK IF THERE IS NOT TOO MUCH LAND (see below).
!                       AGAIN THIS IS TO AVOID CREATING FULL OBSTRUCTION WHEN APPROACHING THE COASTLINE.

                        IF (IDEPTH(I,J) > 0 ) LREALLAND=.TRUE. 
                        LLAND=.TRUE.
                        NIOBSLON=NIOBSLON+1

                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND.  LLEXCLTHRSHOLD(IX,K)) THEN
!                       IF SEA ABOVE THE THRESHOLD THEN ONLY THAT SUBGRID POINT BLOCKS
                        NIOBSLON=NIOBSLON+1
                      ENDIF
                    ENDDO

!                   REVISIT THE LINE IF ANY SUBGRID LAND WAS DETECTED
                    IF (LLAND) THEN
                      IF (LREALLAND) THEN
!                       LINE CONTAINS ACTUAL LAND SUBGRID POINT(S)
!                       SEARCH FOR A CHANGE SEA-LAND-SEA OR VICE VERSA.
                        LNSW=.TRUE.  
                        IF (IDEPTH(I,ILATT) >= IBLOCKDPT(IX,K)) THEN
                          L1ST=.TRUE.
                        ELSE
                          L1ST=.FALSE.
                        ENDIF
                        DO J=ILATT+1,ILATB
                          IF ( ((IDEPTH(I,J) >= IBLOCKDPT(IX,K)) .NEQV. L1ST) .AND. LNSW ) THEN
                            LNSW=.FALSE.
                          ENDIF
                          IF ( ((IDEPTH(I,J) >= IBLOCKDPT(IX,K)) .EQV. L1ST) .AND. .NOT. LNSW ) THEN
!                           LAND IS BLOCKING
                            NIOBSLON=IREINF*(ILATB-ILATT+1)
                            NBLOCKLAND=NBLOCKLAND+1
                            EXIT
                          ENDIF
                        ENDDO
                        IF (LNSW) NIOBSLON=ILATB-ILATT+1

                      ELSE 
!                       SPEUDO LAND
                        IF (PERCENTSHALLOW(IX,K) > PSHALLOWTRHS) THEN
!                         mostly shallow, do not enhance obstruction
                          NIOBSLON=ILATB-ILATT+1
                        ELSEIF (PERCENTLAND(IX,K) < PLANDTRHS) THEN
!                         does not contain too much land
                          NIOBSLON=IREINF*(ILATB-ILATT+1)
                          NBLOCKLAND=NBLOCKLAND+1
                        ELSE
                          NIOBSLON=0
                        ENDIF
                      ENDIF
                    ENDIF

                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO

!                 TOTAL NUMBER OF SUBGRID POINTS, INCLUDING THE ARTIFICIALLY ENHANCED BLOCKING LINE(S)
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR-ILONL+1) + (IREINF-1)*NBLOCKLAND*(ILATB-ILATT+1)

!                 WAVE COMPONENT WILL BE ATTENUATED BY THE RATIO OF ALL BLOCKING SUBGRID POINTS TO THE TOTAL NUMBER OF POINTS
                  IOBSLAT(IX,K,IS) = NINT((1._JWRU-REAL(NOBSTRCT,JWRU)/NTOTPTS)*NOOBSTRT)
                  IOBSLAT(IX,K,IS) = MAX(IOBSLAT(IX,K,IS), 0)


!               AT THE DATELINE, DEALING WITH THE PERIODICITY (simplified version)
                ELSE
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR+ILON-ILONL+1)

                  DO I=1,ILONR
                    NIOBSLON=0
                    DO J=ILATT,ILATB
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K)) THEN
                        NIOBSLON=ILATB-ILATT+1
                        EXIT
                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K) ) THEN
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO
                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO

                  DO I=ILONL,ILON
                    NIOBSLON=0
                    DO J=ILATT,ILATB
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K)) THEN
                        NIOBSLON=ILATB-ILATT+1
                        EXIT
                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K) ) THEN
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO
                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO
                  IOBSLAT(IX,K,IS) = NINT((1._JWRU-REAL(NOBSTRCT,JWRU)/NTOTPTS)*NOOBSTRT)
                  IOBSLAT(IX,K,IS) = MAX(IOBSLAT(IX,K,IS), 0)
                ENDIF

              ENDIF
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ENDDO


!       EAST-WEST OBSTRUCTIONS
!       -----------------------
!       IS=1 is for the west-east advection
!       IS=2 is for the east-west advection
        WRITE(IU06,*) 'CREATE EAST-WEST OBSTRUCTIONS '
        DO IS=1,2
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) &
!$OMP& PRIVATE(K,XLATT,XLATB,ILATT,ILATB,IX) &
!$OMP& PRIVATE(XLONL,XLONR,ILONL,ILONR,NOBSTRCT,NBLOCKLAND) &
!$OMP& PRIVATE(J,NIOBSLAT,LLAND,LREALLAND,I,LNSW,L1ST) &
!$OMP& PRIVATE(NTOTPTS)
!         LOOP OVER MODEL LATITUDES
          DO K=1,NGY
!           LATIDUNAL INDEX OF THE OF THE SUBGRID POINTS THAT ARE INSIDE THE MODEL GRID BOX:  
            XLATT=XLAT(K)+0.5_JWRU*DXDELLA
            XLATB=XLAT(K)-0.5_JWRU*DXDELLA
            IF ( DXDELLA < 0.125_JWRU) THEN
              ILATT = INT((90.0_JWRU- XLATT)*INVRES) + 1
              ILATB = INT((90.0_JWRU- XLATB)*INVRES) + 2
            ELSE
!             It was decided to not correct the double counting for low resolution (should be removed in future)
              ILATT = NINT((90.0_JWRU- XLATT)*INVRES) + 1
              ILATB = NINT((90.0_JWRU- XLATB)*INVRES) + 1
            ENDIF
            ILATT = MAX(1,MIN(ILATT,ILAT))
            ILATB = MAX(1,MIN(ILATB,ILAT))
            IF (ILATB == ILAT+1) ILATB=ILAT

!           LOOP OVER ALL MODEL POINTS FOR A GIVEN LATITUDE
            DO IX=1,NLONRGG(K)
              IF (LLSM(IX,K)) THEN
!               SEA POINT GRID BOX LONGITUDINAL EXTEND :
                IF (IS == 1) THEN
                  XLONL=DAMOWEP + (REAL(IX-2,JWRU))*ZDELLO(K)
                  XLONR=DAMOWEP + (REAL(IX-1,JWRU))*ZDELLO(K) -RESOL
                ELSE
                  XLONL=DAMOWEP + (REAL(IX-1,JWRU))*ZDELLO(K) +RESOL
                  XLONR=DAMOWEP + (REAL(IX,JWRU))*ZDELLO(K)
                ENDIF
                IF (XLONL > 180._JWRU) THEN
                  XLONL=XLONL-360._JWRU
                ENDIF
                IF (XLONR > 180._JWRU) THEN
                  XLONR=XLONR-360._JWRU
                ENDIF

!               LONGITUDINAL INDEX OF THE OF THE SUBGRID POINTS THAT ARE INSIDE THE MODEL GRID BOX:  
                ILONL = NINT((XLONL + 180._JWRU)*INVRES) + 1
                ILONR = NINT((XLONR + 180._JWRU)*INVRES) + 1

!               COMPUTE THE OBSTRUCTIONS:
!               TALLY THE NUMBER OF SUB GRID POINTS THAT ARE POTENTIALLY BLOCKING WAVE PROPAGATION (NOBSTRCT)
                NOBSTRCT=0

!               AWAY FROM THE DATELINE
                IF (ILONL <= ILONR) THEN
                  NBLOCKLAND=0
!                 LOOP OVER SUBGRID LATITUDE LINE:
                  DO J=ILATT,ILATB
                    NIOBSLAT=0
                    LLAND=.FALSE.
                    LREALLAND=.FALSE.

!                   SCAN EACH SUBGRID LONGITUDE:
                    DO I=ILONL,ILONR
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K) ) THEN
!                       IF THE LATITUDE LINE CONTAINS ACTUAL LAND (> 0) THEN THE FULL LONGITUDE WILL BLOCK ONLY IF
!                       THERE IS A SWITCH BACK TO DEPTH < IBLOCKDPT OR VICE VERSA (see below).
!                       THIS IS TO AVOID CREATING FULL OBSTRUCTION WHEN APPROACHING THE COASTLINE.
!                       ELSE IF THE LINE CONTAINS PSEUDO LAND AS DEFINED AS ANYTHING ABOVE IBLOCKDPT(IX,K) (IBLOCKDPT is negative)
!                       THEN THE FULL LINE WILL BLOCK IF THERE IS NOT TOO MUCH LAND (see below).
!                       AGAIN THIS IS TO AVOID CREATING FULL OBSTRUCTION WHEN APPROACHING THE COASTLINE.

                        LLAND=.TRUE.
                        IF (IDEPTH(I,J) > 0 ) LREALLAND=.TRUE. 
                        NIOBSLAT=NIOBSLAT+1 

                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K) ) THEN
!                       IF SEA ABOVE THE THRESHOLD THEN ONLY THAT SUBGRID POINT BLOCKS
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO

!                   REVISIT THE LINE IF ANY SUBGRID LAND WAS DETECTED
                    IF (LLAND) THEN
                      IF (LREALLAND) THEN
!                       LINE CONTAINS ACTUAL LAND SUBGRID POINT(S)
!                       SEARCH FOR A CHANGE SEA-LAND-SEA OR VICE VERSA.
                        LNSW=.TRUE.  
                        IF (IDEPTH(ILONL,J) >= IBLOCKDPT(IX,K) ) THEN
                          L1ST=.TRUE.
                        ELSE
                          L1ST=.FALSE.
                        ENDIF
                        DO I=ILONL+1,ILONR
                          IF ( ((IDEPTH(I,J) >= IBLOCKDPT(IX,K)) .NEQV. L1ST) .AND. LNSW ) THEN
                            LNSW=.FALSE.
                          ENDIF
                          IF ( ((IDEPTH(I,J) >= IBLOCKDPT(IX,K)) .EQV. L1ST) .AND. .NOT. LNSW ) THEN
!                           LAND IS BLOCKING
                            NIOBSLAT=IREINF*(ILONR-ILONL+1)
                            NBLOCKLAND=NBLOCKLAND+1
                            EXIT
                          ENDIF
                        ENDDO
                        IF (LNSW) NIOBSLAT=ILONR-ILONL+1

                      ELSE
!                       SPEUDO LAND
                        IF (PERCENTSHALLOW(IX,K) > PSHALLOWTRHS) THEN
                          NIOBSLAT=ILONR-ILONL+1
                        ELSEIF (PERCENTLAND(IX,K) < PLANDTRHS) THEN
                          NIOBSLAT=IREINF*(ILONR-ILONL+1)
                          NBLOCKLAND=NBLOCKLAND+1
                        ELSE
                          NIOBSLAT=0
                        ENDIF
                      ENDIF
                    ENDIF

                    NOBSTRCT=NOBSTRCT+NIOBSLAT
                  ENDDO

!                 TOTAL NUMBER OF SUBGRID POINTS, INCLUDING THE ARTIFICIALLY ENHANCED BLOCKING LINE(S)
                  NTOTPTS = (ILATB-ILATT+1)*(ILONR-ILONL+1) + (IREINF-1)*NBLOCKLAND*(ILONR-ILONL+1) 

!                 WAVE COMPONENT WILL BE ATTENUATED BY THE RATIO OF ALL BLOCKING SUBGRID POINTS TO THE TOTAL NUMBER OF POINTS
                  IOBSLON(IX,K,IS) = NINT((1._JWRU-REAL(NOBSTRCT,JWRU)/NTOTPTS)*NOOBSTRT)
                  IOBSLON(IX,K,IS) = MAX(IOBSLON(IX,K,IS), 0)


!               AT THE DATELINE, DEALING WITH THE PERIODICITY (simplified version)
                ELSE
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR+ILON-ILONL+1)
                  DO J=ILATT,ILATB
                    NIOBSLAT=0
                    DO I=1,ILONR
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K)) THEN
                        NIOBSLAT=ILONR+ILON-ILONL+1
                        GOTO 1111 
                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K) ) THEN
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO
                    DO I=ILONL,ILON
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K)) THEN
                        NIOBSLAT=ILONR+ILON-ILONL+1
                        EXIT
                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K) ) THEN
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO
1111                CONTINUE
                    NOBSTRCT=NOBSTRCT+NIOBSLAT
                  ENDDO

                  IOBSLON(IX,K,IS) = NINT((1._JWRU-REAL(NOBSTRCT,JWRU)/NTOTPTS)*NOOBSTRT)
                  IOBSLON(IX,K,IS) = MAX(IOBSLON(IX,K,IS), 0)

                ENDIF

              ENDIF
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ENDDO


!       NORTH-WEST-SOUTH-EAST OBSTRUCTIONS (for IPROPAGS = 1)
!       ----------------------------------
!       IS=1 is for the southeast-northwest advection
!       IS=2 is for the northwest-southeast advection

        WRITE(IU06,*) 'CREATE NORTH-WEST-SOUTH-EAST OBSTRUCTIONS '

!       first search north-south
        DO IS=1,2
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) &
!$OMP& PRIVATE(K,KT,KB,STEPT,STEPB,XLATT,XLATB,ILATT,ILATB,IX) &
!$OMP& PRIVATE(XLON,XLONL,XLONR,ILONL,ILONR,NOBSTRCT,NBLOCKLAND) &
!$OMP& PRIVATE(I,NIOBSLON,LLAND,LREALLAND,J,LNSW,L1ST) &
!$OMP& PRIVATE(NTOTPTS)
          DO K=1,NGY
            IF (IS == 1) THEN
              KT=K
              KB=K-1
              STEPT=-RESOL
              STEPB=0._JWRU
            ELSE
              KT=K+1
              KB=K
              STEPT=0._JWRU
              STEPB=RESOL
            ENDIF
            XLATT=XLAT(KT)+STEPT
            XLATB=XLAT(KB)+STEPB
            ILATT = NINT((90.0_JWRU- XLATT)*INVRES) + 1
            ILATT = MAX(1,MIN(ILATT,ILAT))
            ILATB = NINT((90.0_JWRU- XLATB)*INVRES) + 1
            ILATB = MAX(1,MIN(ILATB,ILAT))
            IF (ILATB == ILAT+1) ILATB=ILAT

            DO IX=1,NLONRGG(K)
              IF (LLSM(IX,K)) THEN
                XLON=DAMOWEP + REAL(IX-1,JWRU)*ZDELLO(K)
                XLONL=XLON -(IS-1)*DXDELLA
                IF (XLONL > 180._JWRU) THEN
                  XLONL=XLONL-360._JWRU
                ENDIF
                XLONR=XLON +(2-IS)*DXDELLA
                IF (XLONR > 180._JWRU) THEN
                  XLONR=XLONR-360._JWRU
                ENDIF

                IF ( DXDELLA < 0.125_JWRU) THEN
                  ILONL = INT((XLONL + 180._JWRU)*INVRES) + 1
                  ILONR = INT((XLONR + 180._JWRU)*INVRES) + 2
                ELSE
!                 It was decided to not correct the double counting for low resolution (should be removed in future)
                  ILONL = NINT((XLONL + 180._JWRU)*INVRES) + 1
                  ILONR = NINT((XLONR + 180._JWRU)*INVRES) + 1
                ENDIF

                NOBSTRCT=0

                IF (ILONL <= ILONR) THEN
                  NBLOCKLAND=0
                  DO I=ILONL,ILONR
                    NIOBSLON=0
                    LLAND=.FALSE.
                    LREALLAND=.FALSE.
                    DO J=ILATT,ILATB
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K) ) THEN
!                     IF LAND THEN THE FULL LONGITUDE IS BLOCKED
!                     IF THERE IS A SWITCH BACK TO SEA OR VICE VERSA
!                     (SEE BELOW)
!                     LAND IS DEFINED AS ANYTHING ABOVE IBLOCKDPT(IX,K)
!                     ------------------------------------------
                        IF (IDEPTH(I,J) > 0 ) LREALLAND=.TRUE. 
                        LLAND=.TRUE.
                        NIOBSLON=NIOBSLON+1 
                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K)) THEN
!                     IF SEA ABOVE THE THRESHOLD THEN ONLY THAT
!                     GRID POINTS BLOCKS
!                     ------------------------------------------
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO

                    IF (LLAND) THEN
                      IF (LREALLAND) THEN
                        LNSW=.TRUE.  
                        IF (IDEPTH(I,ILATT) >= IBLOCKDPT(IX,K)) THEN
                          L1ST=.TRUE.
                        ELSE
                          L1ST=.FALSE.
                        ENDIF
                        DO J=ILATT+1,ILATB
                          IF ( ((IDEPTH(I,J) >= IBLOCKDPT(IX,K)) .NEQV. L1ST) .AND. LNSW ) THEN
                            LNSW=.FALSE.
                          ENDIF
                          IF ( ((IDEPTH(I,J) >= IBLOCKDPT(IX,K)) .EQV. L1ST) .AND. .NOT. LNSW ) THEN
!                           LAND IS BLOCKING
                            NIOBSLON=IREINF*(ILATB-ILATT+1)
                            NBLOCKLAND=NBLOCKLAND+1
                            EXIT
                          ENDIF
                        ENDDO
                        IF (LNSW) NIOBSLON=ILATB-ILATT+1
                      ELSE
                        IF (PERCENTSHALLOW(IX,K) > PSHALLOWTRHS) THEN
                          NIOBSLON=ILATB-ILATT+1
                        ELSEIF (PERCENTLAND(IX,K) < PLANDTRHS) THEN
                          NIOBSLON=IREINF*(ILATB-ILATT+1)
                          NBLOCKLAND=NBLOCKLAND+1
                        ELSE
                          NIOBSLON=0
                        ENDIF
                      ENDIF
                    ENDIF

                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR-ILONL+1)+              &
     &                    (IREINF-1)*NBLOCKLAND*(ILATB-ILATT+1)

                  IOBSRLAT(IX,K,IS) = NINT((1._JWRU-REAL(NOBSTRCT,JWRU)/NTOTPTS)*NOOBSTRT)
                ELSE
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR+ILON-ILONL+1)
                  DO I=1,ILONR
                    NIOBSLON=0
                    DO J=ILATT,ILATB
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K)) THEN
                        NIOBSLON=ILATB-ILATT+1
                        EXIT
                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K)) THEN
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO
                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO
                  DO I=ILONL,ILON
                    NIOBSLON=0
                    DO J=ILATT,ILATB
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K)) THEN
                        NIOBSLON=ILATB-ILATT+1
                        EXIT
                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K)) THEN
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO
                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO
                  IOBSRLAT(IX,K,IS) = NINT((1._JWRU-REAL(NOBSTRCT,JWRU)/NTOTPTS)*NOOBSTRT) 
                ENDIF

              ENDIF
            ENDDO

          ENDDO
!$OMP END PARALLEL DO
        ENDDO

!       then search east-west and average 
        DO IS=1,2
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) &
!$OMP& PRIVATE(K,XLATT,XLATB,ILATT,ILATB,IX) &
!$OMP& PRIVATE(XLON,XLONL,XLONR,ILONL,ILONR,NOBSTRCT,NBLOCKLAND) &
!$OMP& PRIVATE(J,NIOBSLAT,LLAND,LREALLAND,I,LNSW,L1ST) &
!$OMP& PRIVATE(NTOTPTS,ITEMPEW,XX)
          DO K=1,NGY
            XLATT=XLAT(K)+(IS-1)*DXDELLA
            XLATB=XLAT(K)-(2-IS)*DXDELLA
            IF ( DXDELLA < 0.125_JWRU) THEN
              ILATT = INT((90.0_JWRU- XLATT)*INVRES) + 1
              ILATB = INT((90.0_JWRU- XLATB)*INVRES) + 2
            ELSE
!             It was decided to not correct the double counting for low resolution (should be removed in future)
              ILATT = NINT((90.0_JWRU- XLATT)*INVRES) + 1
              ILATB = NINT((90.0_JWRU- XLATB)*INVRES) + 1
            ENDIF
            ILATT = MAX(1,MIN(ILATT,ILAT))
            ILATB = MAX(1,MIN(ILATB,ILAT))
            IF (ILATB == ILAT+1) ILATB=ILAT

            DO IX=1,NLONRGG(K)
              IF (LLSM(IX,K)) THEN
                XLON=DAMOWEP + REAL(IX-1,JWRU)*ZDELLO(K)
                IF (IS == 1) THEN
                  XLONL=XLON + RESOL
                  XLONR=XLON + DXDELLA
                ELSE
                  XLONL=XLON -  DXDELLA
                  XLONR=XLON - RESOL
                ENDIF
                IF (XLONL > 180._JWRU) THEN
                  XLONL=XLONL-360._JWRU
                ENDIF
                IF (XLONR > 180._JWRU) THEN
                  XLONR=XLONR-360._JWRU
                ENDIF

                ILONL = NINT((XLONL + 180._JWRU)*INVRES) + 1
                ILONR = NINT((XLONR + 180._JWRU)*INVRES) + 1

                NOBSTRCT=0

                IF (ILONL <= ILONR) THEN
                  NBLOCKLAND=0
                  DO J=ILATT,ILATB
                    NIOBSLAT=0
                    LLAND=.FALSE.
                    LREALLAND=.FALSE.
                    DO I=ILONL,ILONR
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K) ) THEN
!                     IF LAND THEN THE FULL LONGITUDE IS BLOCKED
!                     IF THERE IS A SWITCH BACK TO SEA OR VICE VERSA
!                     (SEE BELOW)
!                     LAND IS DEFINED AS ANYTHING ABOVE IBLOCKDPT(IX,K)
!                     ------------------------------------------
                        LLAND=.TRUE.
                        IF (IDEPTH(I,J) > 0 ) LREALLAND=.TRUE. 
                        NIOBSLAT=NIOBSLAT+1 
                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K)) THEN
!                     IF SEA ABOVE THE THRESHOLD THEN ONLY THAT
!                     GRID POINTS BLOCKS
!                     ------------------------------------------
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO

                    IF (LLAND) THEN
                      IF (LREALLAND) THEN
                        LNSW=.TRUE.  
                        IF (IDEPTH(ILONL,J) >= IBLOCKDPT(IX,K) ) THEN
                          L1ST=.TRUE.
                        ELSE
                          L1ST=.FALSE.
                        ENDIF
                        DO I=ILONL+1,ILONR
                          IF ( ((IDEPTH(I,J) >= IBLOCKDPT(IX,K)) .NEQV. L1ST) .AND. LNSW ) THEN
                            LNSW=.FALSE.
                          ENDIF
                          IF ( ((IDEPTH(I,J) >= IBLOCKDPT(IX,K)) .EQV. L1ST) .AND. .NOT. LNSW ) THEN
!                           LAND IS BLOCKING
                            NIOBSLAT=IREINF*(ILONR-ILONL+1)
                            NBLOCKLAND=NBLOCKLAND+1
                            EXIT
                          ENDIF
                        ENDDO
                        IF (LNSW) NIOBSLAT=ILONR-ILONL+1
                      ELSE
                        IF (PERCENTSHALLOW(IX,K) > PSHALLOWTRHS) THEN
                          NIOBSLAT=ILONR-ILONL+1
                        ELSEIF (PERCENTLAND(IX,K) < PLANDTRHS) THEN
                          NIOBSLAT=IREINF*(ILONR-ILONL+1)
                          NBLOCKLAND=NBLOCKLAND+1
                        ELSE
                          NIOBSLAT=0
                        ENDIF
                      ENDIF
                    ENDIF

                    NOBSTRCT=NOBSTRCT+NIOBSLAT
                  ENDDO

                  NTOTPTS=(ILATB-ILATT+1)*(ILONR-ILONL+1)+              &
     &                    (IREINF-1)*NBLOCKLAND*(ILONR-ILONL+1)
                  ITEMPEW=NINT((1._JWRU-REAL(NOBSTRCT,JWRU)/NTOTPTS)*NOOBSTRT)
                ELSE
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR+ILON-ILONL+1)
                  DO J=ILATT,ILATB
                    NIOBSLAT=0
                    DO I=1,ILONR
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K)) THEN
                        NIOBSLAT=ILONR+ILON-ILONL+1
                        GOTO 2222 
                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K)) THEN
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO
                    DO I=ILONL,ILON
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K)) THEN
                        NIOBSLAT=ILONR+ILON-ILONL+1
                        EXIT
                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K)) THEN
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO
2222                CONTINUE
                    NOBSTRCT=NOBSTRCT+NIOBSLAT
                  ENDDO
                  ITEMPEW=NINT((1._JWRU-REAL(NOBSTRCT,JWRU)/NTOTPTS)*NOOBSTRT)
                ENDIF
                XX=REAL((IOBSRLAT(IX,K,IS)*ITEMPEW),JWRU)
                XX=SQRT(XX)
                IOBSRLAT(IX,K,IS)=MIN(NINT(XX),NOOBSTRT)

              ENDIF
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ENDDO


!       SOUTH-WEST-NORTH-EAST OBSTRUCTIONS (for IPROPAGS = 1)
!       ----------------------------------
!       IS=1 is for the southwest-northeast advection
!       IS=2 is for the northeast-southwest advection

        WRITE(IU06,*) 'CREATE SOUTH-WEST-NORTH-EAST OBSTRUCTIONS '

!       first search north-south
        DO IS=1,2
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) &
!$OMP& PRIVATE(K,KT,KB,STEPT,STEPB,XLATT,XLATB,ILATT,ILATB,IX) &
!$OMP& PRIVATE(XLON,XLONL,XLONR,ILONL,ILONR,NOBSTRCT,NBLOCKLAND) &
!$OMP& PRIVATE(I,NIOBSLON,LLAND,LREALLAND,J,LNSW,L1ST) &
!$OMP& PRIVATE(NTOTPTS)
          DO K=1,NGY
            IF (IS == 1) THEN
              KT=K
              KB=K-1
              STEPT=-RESOL
              STEPB=0._JWRU
            ELSE
              KT=K+1
              KB=K
              STEPT=0._JWRU
              STEPB=RESOL
            ENDIF
            XLATT=XLAT(KT)+STEPT
            XLATB=XLAT(KB)+STEPB
            ILATT = NINT((90.0_JWRU- XLATT)*INVRES) + 1
            ILATT = MAX(1,MIN(ILATT,ILAT))
            ILATB = NINT((90.0_JWRU- XLATB)*INVRES) + 1
            ILATB = MAX(1,MIN(ILATB,ILAT))
            IF (ILATB == ILAT+1) ILATB=ILAT

            DO IX=1,NLONRGG(K)
              IF (LLSM(IX,K)) THEN
                XLON=DAMOWEP + REAL(IX-1,JWRU)*ZDELLO(K)
                XLONL=XLON -(2-IS)*DXDELLA
                IF (XLONL > 180._JWRU) THEN
                  XLONL=XLONL-360._JWRU
                ENDIF
                XLONR=XLON +(IS-1)*DXDELLA
                IF (XLONR > 180._JWRU) THEN
                  XLONR=XLONR-360._JWRU
                ENDIF

                IF ( DXDELLA < 0.125_JWRU) THEN
                  ILONL = INT((XLONL + 180._JWRU)*INVRES) + 1
                  ILONR = INT((XLONR + 180._JWRU)*INVRES) + 2
                ELSE
!                 It was decided to not correct the double counting for low resolution (should be removed in future)
                  ILONL = NINT((XLONL + 180._JWRU)*INVRES) + 1
                  ILONR = NINT((XLONR + 180._JWRU)*INVRES) + 1
                ENDIF

                NOBSTRCT=0

                IF (ILONL <= ILONR) THEN
                  NBLOCKLAND=0
                  DO I=ILONL,ILONR
                    NIOBSLON=0
                    LLAND=.FALSE.
                    LREALLAND=.FALSE.
                    DO J=ILATT,ILATB
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K) ) THEN
!                     IF LAND THEN THE FULL LONGITUDE IS BLOCKED
!                     IF THERE IS A SWITCH BACK TO SEA OR VICE VERSA
!                     (SEE BELOW)
!                     LAND IS DEFINED AS ANYTHING ABOVE IBLOCKDPT(IX,K)
!                     ------------------------------------------
                        IF (IDEPTH(I,J) > 0 ) LREALLAND=.TRUE. 
                        LLAND=.TRUE.
                        NIOBSLON=NIOBSLON+1 
                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K)) THEN
!                     IF SEA ABOVE THE THRESHOLD THEN ONLY THAT
!                     GRID POINTS BLOCKS
!                     ------------------------------------------
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO

                    IF (LLAND) THEN
                      IF (LREALLAND) THEN
                        LNSW=.TRUE.  
                        IF (IDEPTH(I,ILATT) >= IBLOCKDPT(IX,K)) THEN
                          L1ST=.TRUE.
                        ELSE
                          L1ST=.FALSE.
                        ENDIF
                        DO J=ILATT+1,ILATB
                          IF ( ((IDEPTH(I,J) >= IBLOCKDPT(IX,K)) .NEQV. L1ST) .AND. LNSW ) THEN
                            LNSW=.FALSE.
                          ENDIF
                          IF ( ((IDEPTH(I,J) >= IBLOCKDPT(IX,K)) .EQV. L1ST) .AND. .NOT. LNSW ) THEN
!                           LAND IS BLOCKING
                            NIOBSLON=IREINF*(ILATB-ILATT+1)
                            NBLOCKLAND=NBLOCKLAND+1
                            EXIT
                          ENDIF
                        ENDDO
                        IF (LNSW) NIOBSLON=ILATB-ILATT+1
                      ELSE
                        IF (PERCENTSHALLOW(IX,K) > PSHALLOWTRHS) THEN
                          NIOBSLON=ILATB-ILATT+1
                        ELSEIF (PERCENTLAND(IX,K) < PLANDTRHS) THEN
                          NIOBSLON=IREINF*(ILATB-ILATT+1)
                          NBLOCKLAND=NBLOCKLAND+1
                        ELSE
                          NIOBSLON=0
                        ENDIF
                      ENDIF
                    ENDIF

                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR-ILONL+1)+              &
     &                    (IREINF-1)*NBLOCKLAND*(ILATB-ILATT+1)

                  IOBSRLON(IX,K,IS) = NINT((1._JWRU-REAL(NOBSTRCT,JWRU)/NTOTPTS)*NOOBSTRT)
                ELSE
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR+ILON-ILONL+1)
                  DO I=1,ILONR
                    NIOBSLON=0
                    DO J=ILATT,ILATB
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K)) THEN
                        NIOBSLON=ILATB-ILATT+1
                        EXIT
                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K)) THEN
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO
                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO
                  DO I=ILONL,ILON
                    NIOBSLON=0
                    DO J=ILATT,ILATB
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K)) THEN
                        NIOBSLON=ILATB-ILATT+1
                        EXIT
                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K)) THEN
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO
                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO
                  IOBSRLON(IX,K,IS) = NINT((1._JWRU-REAL(NOBSTRCT,JWRU)/NTOTPTS)*NOOBSTRT)
                ENDIF

              ENDIF
            ENDDO

          ENDDO
!$OMP END PARALLEL DO
        ENDDO

!       then search east-west and average 
        DO IS=1,2
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) &
!$OMP& PRIVATE(K,XLATT,XLATB,ILATT,ILATB,IX) &
!$OMP& PRIVATE(XLON,XLONL,XLONR,ILONL,ILONR,NOBSTRCT,NBLOCKLAND) &
!$OMP& PRIVATE(J,NIOBSLAT,LLAND,LREALLAND,I,LNSW,L1ST) &
!$OMP& PRIVATE(NTOTPTS,ITEMPEW,XX)
          DO K=1,NGY
            XLATT=XLAT(K)+(IS-1)*DXDELLA
            XLATB=XLAT(K)-(2-IS)*DXDELLA
            IF ( DXDELLA < 0.125_JWRU) THEN
              ILATT = INT((90.0_JWRU- XLATT)*INVRES) + 1
              ILATB = INT((90.0_JWRU- XLATB)*INVRES) + 2
            ELSE
!             It was decided to not correct the double counting for low resolution (should be removed in future)
              ILATT = NINT((90.0_JWRU- XLATT)*INVRES) + 1
              ILATB = NINT((90.0_JWRU- XLATB)*INVRES) + 1
            ENDIF
            ILATT = MAX(1,MIN(ILATT,ILAT))
            ILATB = MAX(1,MIN(ILATB,ILAT))
            IF (ILATB == ILAT+1) ILATB=ILAT

            DO IX=1,NLONRGG(K)
              IF (LLSM(IX,K)) THEN
                XLON=DAMOWEP + REAL(IX-1,JWRU)*ZDELLO(K)
                IF (IS == 1) THEN
                  XLONL=XLON + RESOL
                  XLONR=XLON + DXDELLA
                ELSE
                  XLONL=XLON - DXDELLA
                  XLONR=XLON - RESOL
                ENDIF
                IF (XLONL > 180._JWRU) THEN
                  XLONL=XLONL-360._JWRU
                ENDIF
                IF (XLONR > 180._JWRU) THEN
                  XLONR=XLONR-360._JWRU
                ENDIF

                ILONL = NINT((XLONL + 180._JWRU)*INVRES) + 1
                ILONR = NINT((XLONR + 180._JWRU)*INVRES) + 1

                NOBSTRCT=0

                IF (ILONL <= ILONR) THEN
                  NBLOCKLAND=0
                  DO J=ILATT,ILATB
                    NIOBSLAT=0
                    LLAND=.FALSE.
                    LREALLAND=.FALSE.
                    DO I=ILONL,ILONR
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K) ) THEN
!                     IF LAND THEN THE FULL LONGITUDE IS BLOCKED
!                     IF THERE IS A SWITCH BACK TO SEA OR VICE VERSA
!                     (SEE BELOW)
!                     LAND IS DEFINED AS ANYTHING ABOVE IBLOCKDPT(IX,K)
!                     ------------------------------------------
                        LLAND=.TRUE.
                        IF (IDEPTH(I,J) > 0 ) LREALLAND=.TRUE. 
                        NIOBSLAT=NIOBSLAT+1 
                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K)) THEN
!                     IF SEA ABOVE THE THRESHOLD THEN ONLY THAT
!                     GRID POINTS BLOCKS
!                     ------------------------------------------
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO

                    IF (LLAND) THEN
                      IF (LREALLAND) THEN
                        LNSW=.TRUE.  
                        IF (IDEPTH(ILONL,J) >= IBLOCKDPT(IX,K) ) THEN
                          L1ST=.TRUE.
                        ELSE
                          L1ST=.FALSE.
                        ENDIF
                        DO I=ILONL+1,ILONR
                          IF ( ((IDEPTH(I,J) >= IBLOCKDPT(IX,K)) .NEQV. L1ST) .AND. LNSW ) THEN
                            LNSW=.FALSE.
                          ENDIF
                          IF ( ((IDEPTH(I,J) >= IBLOCKDPT(IX,K)) .EQV. L1ST) .AND. .NOT. LNSW ) THEN
!                           LAND IS BLOCKING
                            NIOBSLAT=IREINF*(ILONR-ILONL+1)
                            NBLOCKLAND=NBLOCKLAND+1
                            EXIT
                          ENDIF
                        ENDDO
                        IF (LNSW) NIOBSLAT=ILONR-ILONL+1
                      ELSE
                        IF (PERCENTSHALLOW(IX,K) > PSHALLOWTRHS) THEN
                          NIOBSLAT=ILONR-ILONL+1
                        ELSEIF (PERCENTLAND(IX,K) < PLANDTRHS) THEN
                          NIOBSLAT=IREINF*(ILONR-ILONL+1)
                          NBLOCKLAND=NBLOCKLAND+1
                        ELSE
                          NIOBSLAT=0
                        ENDIF
                      ENDIF
                    ENDIF

                    NOBSTRCT=NOBSTRCT+NIOBSLAT
                  ENDDO

                  NTOTPTS=(ILATB-ILATT+1)*(ILONR-ILONL+1)+              &
     &                    (IREINF-1)*NBLOCKLAND*(ILONR-ILONL+1)
                  ITEMPEW=NINT((1._JWRU-REAL(NOBSTRCT,JWRU)/NTOTPTS)*NOOBSTRT)
                ELSE
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR+ILON-ILONL+1)
                  DO J=ILATT,ILATB
                    NIOBSLAT=0
                    DO I=1,ILONR
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K)) THEN
                        NIOBSLAT=ILONR+ILON-ILONL+1
                        GOTO 3333 
                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K)) THEN
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO
                    DO I=ILONL,ILON
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K)) THEN
                        NIOBSLAT=ILONR+ILON-ILONL+1
                        EXIT
                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K)) THEN
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO
3333                CONTINUE
                    NOBSTRCT=NOBSTRCT+NIOBSLAT
                  ENDDO
                  ITEMPEW=NINT((1._JWRU-REAL(NOBSTRCT,JWRU)/NTOTPTS)*NOOBSTRT)
                ENDIF
                XX=REAL((IOBSRLON(IX,K,IS)*ITEMPEW),JWRU)
                XX=SQRT(XX)
                IOBSRLON(IX,K,IS)=MIN(NINT(XX),NOOBSTRT)

              ENDIF
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ENDDO


!       GRID CORNER POINT OBSTRUCTIONS (for IPROPAGS = 2)
!       ------------------------------
!       IS=1 is for the northeast-southwest advection
!       IS=2 is for the southeast-northwest advection
!       IS=3 is for the southwest-northeast advection
!       IS=4 is for the northwest-southeast advection

        WRITE(IU06,*) 'CREATE GRID CORNER OBSTRUCTIONS '

!       first search north-south
        DO IS=1,4
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) &
!$OMP& PRIVATE(K,KT,KB,STEPT,STEPB,XLATT,XLATB,ILATT,ILATB,IX) &
!$OMP& PRIVATE(XLON,XLONL,XLONR,ILONL,ILONR,NOBSTRCT,NBLOCKLAND) &
!$OMP& PRIVATE(I,NIOBSLON,LLAND,LREALLAND,J,LNSW,L1ST) &
!$OMP& PRIVATE(NTOTPTS)
          DO K=1,NGY
            IF (IS == 1) THEN
              KT=K+1
              KB=K
              STEPT=0._JWRU
              STEPB=RESOL
            ELSE IF (IS == 2) THEN
              KT=K
              KB=K-1
              STEPT=-RESOL
              STEPB=0._JWRU
            ELSE IF (IS == 3) THEN
              KT=K
              KB=K-1
              STEPT=-RESOL
              STEPB=0._JWRU
            ELSE IF (IS == 4) THEN
              KT=K+1
              KB=K
              STEPT=0._JWRU
              STEPB=RESOL
            ENDIF

            XLATT=XLAT(KT)+STEPT
            XLATB=XLAT(KB)+STEPB
            ILATT = NINT((90.0_JWRU- XLATT)*INVRES) + 1
            ILATT = MAX(1,MIN(ILATT,ILAT))
            ILATB = NINT((90.0_JWRU- XLATB)*INVRES) + 1
            ILATB = MAX(1,MIN(ILATB,ILAT))
            IF (ILATB == ILAT+1) ILATB=ILAT

            DO IX=1,NLONRGG(K)
              IF (LLSM(IX,K)) THEN
                XLON=DAMOWEP + REAL(IX-1,JWRU)*ZDELLO(K)
                XLONL=XLON -((IS-1)/2)*ZDELLO(K)
                IF (XLONL > 180._JWRU) THEN
                  XLONL=XLONL-360._JWRU
                ENDIF
                XLONR=XLON +((4-IS)/2)*ZDELLO(K)
                IF (XLONR > 180._JWRU) THEN
                  XLONR=XLONR-360._JWRU
                ENDIF

                IF ( DXDELLA < 0.125_JWRU) THEN
                  ILONL = INT((XLONL + 180._JWRU)*INVRES) + 1
                  ILONR = INT((XLONR + 180._JWRU)*INVRES) + 2
                ELSE
!                 It was decided to not correct the double counting for low resolution (should be removed in future)
                  ILONL = NINT((XLONL + 180._JWRU)*INVRES) + 1
                  ILONR = NINT((XLONR + 180._JWRU)*INVRES) + 1
                ENDIF

                NOBSTRCT=0

                IF (ILONL <= ILONR) THEN
                  NBLOCKLAND=0
                  DO I=ILONL,ILONR
                    NIOBSLON=0
                    LLAND=.FALSE.
                    LREALLAND=.FALSE.
                    DO J=ILATT,ILATB
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K) ) THEN
!                     IF LAND THEN THE FULL LONGITUDE IS BLOCKED
!                     IF THERE IS A SWITCH BACK TO SEA OR VICE VERSA
!                     (SEE BELOW)
!                     LAND IS DEFINED AS ANYTHING ABOVE IBLOCKDPT(IX,K)
!                     ------------------------------------------
                        IF (IDEPTH(I,J) > 0 ) LREALLAND=.TRUE. 
                        LLAND=.TRUE.
                        NIOBSLON=NIOBSLON+1 
                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K)) THEN
!                     IF SEA ABOVE THE THRESHOLD THEN ONLY THAT
!                     GRID POINTS BLOCKS
!                     ------------------------------------------
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO

                    IF (LLAND) THEN
                      IF (LREALLAND) THEN
                        LNSW=.TRUE.  
                        IF (IDEPTH(I,ILATT) >= IBLOCKDPT(IX,K)) THEN
                          L1ST=.TRUE.
                        ELSE
                          L1ST=.FALSE.
                        ENDIF
                        DO J=ILATT+1,ILATB
                          IF ( ((IDEPTH(I,J) >= IBLOCKDPT(IX,K)) .NEQV. L1ST) .AND. LNSW ) THEN
                            LNSW=.FALSE.
                          ENDIF
                          IF ( ((IDEPTH(I,J) >= IBLOCKDPT(IX,K)) .EQV. L1ST) .AND. .NOT. LNSW ) THEN
!                           LAND IS BLOCKING
                            NIOBSLON=IREINF*(ILATB-ILATT+1)
                            NBLOCKLAND=NBLOCKLAND+1
                            EXIT
                          ENDIF
                        ENDDO
                        IF (LNSW) NIOBSLON=ILATB-ILATT+1
                      ELSE
                        IF (PERCENTSHALLOW(IX,K) > PSHALLOWTRHS) THEN
                          NIOBSLON=ILATB-ILATT+1
                        ELSEIF (PERCENTLAND(IX,K) < PLANDTRHS) THEN
                          NIOBSLON=IREINF*(ILATB-ILATT+1)
                          NBLOCKLAND=NBLOCKLAND+1
                        ELSE
                          NIOBSLON=0
                        ENDIF
                      ENDIF
                    ENDIF

                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR-ILONL+1)+              &
     &                    (IREINF-1)*NBLOCKLAND*(ILATB-ILATT+1)

                  IOBSCOR(IX,K,IS) = NINT((1._JWRU-REAL(NOBSTRCT,JWRU)/NTOTPTS)*NOOBSTRT)
                  IOBSCOR(IX,K,IS) = MAX(IOBSCOR(IX,K,IS), 0)
                ELSE
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR+ILON-ILONL+1)
                  DO I=1,ILONR
                    NIOBSLON=0
                    DO J=ILATT,ILATB
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K)) THEN
                        NIOBSLON=ILATB-ILATT+1
                        EXIT
                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K)) THEN
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO
                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO
                  DO I=ILONL,ILON
                    NIOBSLON=0
                    DO J=ILATT,ILATB
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K)) THEN
                        NIOBSLON=ILATB-ILATT+1
                        EXIT
                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K)) THEN
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO
                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO
                  IOBSCOR(IX,K,IS) = NINT((1._JWRU-REAL(NOBSTRCT,JWRU)/NTOTPTS)*NOOBSTRT)
                  IOBSCOR(IX,K,IS) = MAX(IOBSCOR(IX,K,IS), 0)
                ENDIF

              ENDIF
            ENDDO

          ENDDO
!$OMP END PARALLEL DO
        ENDDO

!       then search east-west and average 
        DO IS=1,4
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) &
!$OMP& PRIVATE(K,XLATT,XLATB,ILATT,ILATB,IX) &
!$OMP& PRIVATE(XLON,XLONL,XLONR,ILONL,ILONR,NOBSTRCT,NBLOCKLAND) &
!$OMP& PRIVATE(J,NIOBSLAT,LLAND,LREALLAND,I,LNSW,L1ST) &
!$OMP& PRIVATE(NTOTPTS,ITEMPEW,XX)
          DO K=1,NGY
            IF (IS == 1 .OR. IS == 4) THEN
              XLATT=XLAT(K)+DXDELLA
              XLATB=XLAT(K)
            ELSE
              XLATT=XLAT(K)
              XLATB=XLAT(K)-DXDELLA
            ENDIF
            IF ( DXDELLA < 0.125_JWRU) THEN
              ILATT = INT((90.0_JWRU- XLATT)*INVRES) + 1
              ILATB = INT((90.0_JWRU- XLATB)*INVRES) + 2
            ELSE
!             It was decided to not correct the double counting for low resolution (should be removed in future)
              ILATT = NINT((90.0_JWRU- XLATT)*INVRES) + 1
              ILATB = NINT((90.0_JWRU- XLATB)*INVRES) + 1
            ENDIF
            ILATT = MAX(1,MIN(ILATT,ILAT))
            ILATB = MAX(1,MIN(ILATB,ILAT))
            IF (ILATB == ILAT+1) ILATB=ILAT

            DO IX=1,NLONRGG(K)
              IF (LLSM(IX,K)) THEN
                XLON=DAMOWEP + REAL(IX-1,JWRU)*ZDELLO(K)
                IF (IS == 1 .OR. IS == 2) THEN
                  XLONL=XLON + RESOL
                  XLONR=XLON + ZDELLO(K) 
                ELSE
                  XLONL=XLON - ZDELLO(K) 
                  XLONR=XLON - RESOL
                ENDIF
                IF (XLONL > 180._JWRU) THEN
                  XLONL=XLONL-360._JWRU
                ENDIF
                IF (XLONR > 180._JWRU) THEN
                  XLONR=XLONR-360._JWRU
                ENDIF

                ILONL = NINT((XLONL + 180._JWRU)*INVRES) + 1
                ILONR = NINT((XLONR + 180._JWRU)*INVRES) + 1

                NOBSTRCT=0

                IF (ILONL <= ILONR) THEN
                  NBLOCKLAND=0
                  DO J=ILATT,ILATB
                    NIOBSLAT=0
                    LLAND=.FALSE.
                    LREALLAND=.FALSE.
                    DO I=ILONL,ILONR
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K) ) THEN
!                     IF LAND THEN THE FULL LONGITUDE IS BLOCKED
!                     IF THERE IS A SWITCH BACK TO SEA OR VICE VERSA
!                     (SEE BELOW)
!                     LAND IS DEFINED AS ANYTHING ABOVE IBLOCKDPT(IX,K)
!                     ------------------------------------------
                        LLAND=.TRUE.
                        IF (IDEPTH(I,J) > 0 ) LREALLAND=.TRUE. 
                        NIOBSLAT=NIOBSLAT+1 
                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K)) THEN
!                     IF SEA ABOVE THE THRESHOLD THEN ONLY THAT
!                     GRID POINTS BLOCKS
!                     ------------------------------------------
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO

                    IF (LLAND) THEN
                      IF (LREALLAND) THEN
                        LNSW=.TRUE.  
                        IF (IDEPTH(ILONL,J) >= IBLOCKDPT(IX,K) ) THEN
                          L1ST=.TRUE.
                        ELSE
                          L1ST=.FALSE.
                        ENDIF
                        DO I=ILONL+1,ILONR
                          IF ( ((IDEPTH(I,J) >= IBLOCKDPT(IX,K)) .NEQV. L1ST) .AND. LNSW ) THEN
                            LNSW=.FALSE.
                          ENDIF
                          IF ( ((IDEPTH(I,J) >= IBLOCKDPT(IX,K)) .EQV. L1ST) .AND. .NOT. LNSW ) THEN
!                           LAND IS BLOCKING
                            NIOBSLAT=IREINF*(ILONR-ILONL+1)
                            NBLOCKLAND=NBLOCKLAND+1
                            EXIT
                          ENDIF
                        ENDDO
                        IF (LNSW) NIOBSLAT=ILONR-ILONL+1
                      ELSE
                        IF (PERCENTSHALLOW(IX,K) > PSHALLOWTRHS) THEN
                          NIOBSLAT=ILONR-ILONL+1
                        ELSEIF (PERCENTLAND(IX,K) < PLANDTRHS) THEN
                          NIOBSLAT=IREINF*(ILONR-ILONL+1)
                          NBLOCKLAND=NBLOCKLAND+1
                        ELSE
                          NIOBSLAT=0
                        ENDIF
                      ENDIF
                    ENDIF

                    NOBSTRCT=NOBSTRCT+NIOBSLAT
                  ENDDO

                  NTOTPTS=(ILATB-ILATT+1)*(ILONR-ILONL+1)+              &
     &                    (IREINF-1)*NBLOCKLAND*(ILONR-ILONL+1)
                  ITEMPEW=NINT((1._JWRU-REAL(NOBSTRCT,JWRU)/NTOTPTS)*NOOBSTRT)
                ELSE
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR+ILON-ILONL+1)
                  DO J=ILATT,ILATB
                    NIOBSLAT=0
                    DO I=1,ILONR
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K)) THEN
                        NIOBSLAT=ILONR+ILON-ILONL+1
                        GOTO 4444 
                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K)) THEN
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO
                    DO I=ILONL,ILON
                      IF (IDEPTH(I,J) >= IBLOCKDPT(IX,K)) THEN
                        NIOBSLAT=ILONR+ILON-ILONL+1
                        EXIT
                      ELSEIF (IDEPTH(I,J) >= ITHRSHOLD(IX,K) .AND. LLEXCLTHRSHOLD(IX,K)) THEN
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO
4444                CONTINUE
                    NOBSTRCT=NOBSTRCT+NIOBSLAT
                  ENDDO
                  ITEMPEW = NINT((1._JWRU-REAL(NOBSTRCT,JWRU)/NTOTPTS)*NOOBSTRT)
                  ITEMPEW = MAX(ITEMPEW, 0)
                ENDIF
                XX=REAL((IOBSCOR(IX,K,IS)*ITEMPEW),JWRU)
                XX=PENHCOR*SQRT(XX)
                IOBSCOR(IX,K,IS)=MIN(NINT(XX),NOOBSTRT)

              ENDIF
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ENDDO

        ENDIF ! end if xdella too small, do not compute obstructions


!       OUTPUT OBSTRUCTIONS
!       FOR PLOTTING

        IF (LLPRINT) THEN

          WRITE(CFR,'(I2.2)') M

          FILENM='obstructions_S_N_'//CFR//'.dat'
          OPEN(22,file=FILENM)
          WRITE(22,'(a4)') '#GEO'
          WRITE(22,'(a11)') '#FORMAT LLV'
          WRITE(22,'(a5)') '#DATA'

          FILENM='obstructions_N_S_'//CFR//'.dat'
          OPEN(23,file=FILENM)
          WRITE(23,'(a4)') '#GEO'
          WRITE(23,'(a11)') '#FORMAT LLV'
          WRITE(23,'(a5)') '#DATA'

          FILENM='obstructions_W_E_'//CFR//'.dat'
          OPEN(24,file=FILENM)
          WRITE(24,'(a4)') '#GEO'
          WRITE(24,'(a11)') '#FORMAT LLV'
          WRITE(24,'(a5)') '#DATA'

          FILENM='obstructions_E_W_'//CFR//'.dat'
          OPEN(25,file=FILENM)
          WRITE(25,'(a4)') '#GEO'
          WRITE(25,'(a11)') '#FORMAT LLV'
          WRITE(25,'(a5)') '#DATA'

          FILENM='obstructions_SE_NW_'//CFR//'.dat'
          OPEN(26,file=FILENM)
          WRITE(26,'(a4)') '#GEO'
          WRITE(26,'(a11)') '#FORMAT LLV'
          WRITE(26,'(a5)') '#DATA'

          FILENM='obstructions_NW_SE_'//CFR//'.dat'
          OPEN(27,file=FILENM)
          WRITE(27,'(a4)') '#GEO'
          WRITE(27,'(a11)') '#FORMAT LLV'
          WRITE(27,'(a5)') '#DATA'

          FILENM='obstructions_SW_NE_'//CFR//'.dat'
          OPEN(28,file=FILENM)
          WRITE(28,'(a4)') '#GEO'
          WRITE(28,'(a11)') '#FORMAT LLV'
          WRITE(28,'(a5)') '#DATA'

          FILENM='obstructions_NE_SW_'//CFR//'.dat'
          OPEN(29,file=FILENM)
          WRITE(29,'(a4)') '#GEO'
          WRITE(29,'(a11)') '#FORMAT LLV'
          WRITE(29,'(a5)') '#DATA'

          FILENM='obstructions_CORNER_NE_SW_'//CFR//'.dat'
          OPEN(30,file=FILENM)
          WRITE(30,'(a4)') '#GEO'
          WRITE(30,'(a11)') '#FORMAT LLV'
          WRITE(30,'(a5)') '#DATA'

          FILENM='obstructions_CORNER_SE_NW_'//CFR//'.dat'
          OPEN(31,file=FILENM)
          WRITE(31,'(a4)') '#GEO'
          WRITE(31,'(a11)') '#FORMAT LLV'
          WRITE(31,'(a5)') '#DATA'

          FILENM='obstructions_CORNER_SW_NE'//CFR//'.dat'
          OPEN(32,file=FILENM)
          WRITE(32,'(a4)') '#GEO'
          WRITE(32,'(a11)') '#FORMAT LLV'
          WRITE(32,'(a5)') '#DATA'

          FILENM='obstructions_CORNER_NW_SE'//CFR//'.dat'
          OPEN(33,file=FILENM)
          WRITE(33,'(a4)') '#GEO'
          WRITE(33,'(a11)') '#FORMAT LLV'
          WRITE(33,'(a5)') '#DATA'

          IUNIT=21

          DO IS =1,2
            IUNIT=IUNIT+1
            IF (IS == 1) THEN
              STEPLAT=-0.25_JWRU*DXDELLA
            ELSE
              STEPLAT=0.25_JWRU*DXDELLA
            ENDIF 
            DO K=1,NGY
               DO IX=1,NLONRGG(K)
                 XLON=DAMOWEP + REAL(IX-1,JWRU)*ZDELLO(K)
                 IF (XLON > 180._JWRU) THEN
                   XLON=XLON-360._JWRU
                 ENDIF
                IF (ALATB <= XLAT(K) .AND. XLAT(K) <= ALATT .AND.       &
     &              ALONL <= XLON .AND. XLON <= ALONR ) THEN
                   IF (LLSM(IX,K) .AND. IOBSLAT(IX,K,IS) < NOOBSTRT ) THEN
                      WRITE(IUNIT,'(2(1X,F8.3),1X,I4)')                 &
     &                XLON,XLAT(K)+STEPLAT,IOBSLAT(IX,K,IS)
                   ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO

          DO IS =1,2
            IUNIT=IUNIT+1
            IF (IS == 1) THEN
              STEPLON=-0.25_JWRU*DXDELLO
            ELSE
              STEPLON=0.25_JWRU*DXDELLO
            ENDIF 
            DO K=1,NGY
               DO IX=1,NLONRGG(K)
                 XLON=DAMOWEP + REAL(IX-1,JWRU)*ZDELLO(K)
                 IF (XLON > 180._JWRU) THEN
                   XLON=XLON-360._JWRU
                 ENDIF
                IF (ALATB <= XLAT(K) .AND. XLAT(K) <= ALATT .AND.       &
     &              ALONL <= XLON .AND. XLON <= ALONR ) THEN
                   IF (LLSM(IX,K) .AND. IOBSLON(IX,K,IS) < NOOBSTRT ) THEN
                      WRITE(IUNIT,'(2(1X,F8.3),1X,I4)')                 &
     &                XLON+STEPLON,XLAT(K),IOBSLON(IX,K,IS)
                   ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO

          DO IS =1,2
            IUNIT=IUNIT+1
            IF (IS == 1) THEN
              STEPLAT=-0.25_JWRU*DXDELLA
            ELSE
              STEPLAT=0.25_JWRU*DXDELLA
            ENDIF 
            DO K=1,NGY
               DO IX=1,NLONRGG(K)
                 XLON=DAMOWEP + REAL(IX-1,JWRU)*ZDELLO(K)
                 IF (XLON > 180._JWRU) THEN
                   XLON=XLON-360._JWRU
                 ENDIF
                IF (ALATB <= XLAT(K) .AND. XLAT(K) <= ALATT .AND.       &
     &              ALONL <= XLON .AND. XLON <= ALONR ) THEN
                   IF (LLSM(IX,K) .AND. IOBSRLAT(IX,K,IS) < NOOBSTRT ) THEN
                      WRITE(IUNIT,'(2(1X,F8.3),1X,I4)')                 &
     &                XLON,XLAT(K)+STEPLAT,IOBSRLAT(IX,K,IS)
                   ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO

          DO IS =1,2
            IUNIT=IUNIT+1
            IF (IS == 1) THEN
              STEPLON=-0.25_JWRU*DXDELLO
            ELSE
              STEPLON=0.25_JWRU*DXDELLO
            ENDIF 
            DO K=1,NGY
               DO IX=1,NLONRGG(K)
                 XLON=DAMOWEP + REAL(IX-1,JWRU)*ZDELLO(K)
                 IF (XLON > 180._JWRU) THEN
                   XLON=XLON-360._JWRU
                 ENDIF
                IF (ALATB <= XLAT(K) .AND. XLAT(K) <= ALATT .AND.       &
     &             ALONL <= XLON .AND. XLON <= ALONR ) THEN
                   IF (LLSM(IX,K) .AND. IOBSRLON(IX,K,IS) < NOOBSTRT ) THEN
                      WRITE(IUNIT,'(2(1X,F8.3),1X,I4)')                 &
     &                XLON+STEPLON,XLAT(K),IOBSRLON(IX,K,IS)
                   ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO

          DO IS =1,4
            IUNIT=IUNIT+1
            IF (IS == 1) THEN
              STEPLON=-0.25_JWRU*DXDELLO
            ELSE
              STEPLON=0.25_JWRU*DXDELLO
            ENDIF 
            DO K=1,NGY
               DO IX=1,NLONRGG(K)
                 XLON=DAMOWEP + REAL(IX-1,JWRU)*ZDELLO(K)
                 IF (XLON > 180._JWRU) THEN
                   XLON=XLON-360._JWRU
                 ENDIF
                IF (ALATB <= XLAT(K) .AND. XLAT(K) <= ALATT .AND.       &
     &              ALONL <= XLON .AND. XLON <= ALONR ) THEN
                   IF (LLSM(IX,K) .AND. IOBSRLON(IX,K,IS) < NOOBSTRT ) THEN
                      WRITE(IUNIT,'(2(1X,F8.3),1X,I4)')                 &
     &                XLON+STEPLON,XLAT(K),IOBSCOR(IX,K,IS)
                   ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO

        ENDIF



!       OUTPUT OBSTRUCTIONS
!       -------------------

        IF ( LLGRIBOUT ) THEN
!         GRIB OUTPUT (convert them to values between 0 and 1)
          CALL KTOOBS(IU06)


          DO IP = 0, NPROPAGS 
            DO IANG = 1, NANG

              IS = KTOIS(IANG,IP)

              IF ( IS > 0 .AND. LLOBSTRON ) THEN
                IOBSRT = KTOOBSTRUCT(IANG,IP) 
                SELECT CASE(IOBSRT)

                CASE(1)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(K,KNS,IX)
                  DO K=1,NGY
                    KNS=NGY-K+1
                    DO IX=1,NLONRGG(K)
                      FIELD(IX,KNS) = REAL(ILSM(IX,K)*IOBSLAT(IX,K,IS),JWRB) * ZCONV  + (1-ILSM(IX,K))*ZMISS
                    ENDDO
                  ENDDO
!$OMP END PARALLEL DO

                CASE(2)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(K,KNS,IX)
                  DO K=1,NGY
                    KNS=NGY-K+1
                    DO IX=1,NLONRGG(K)
                      FIELD(IX,KNS) = REAL(ILSM(IX,K)*IOBSLON(IX,K,IS),JWRB) * ZCONV  + (1-ILSM(IX,K))*ZMISS
                    ENDDO
                  ENDDO
!$OMP END PARALLEL DO

                CASE(3)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(K,KNS,IX)
                  DO K=1,NGY
                    KNS=NGY-K+1
                    DO IX=1,NLONRGG(K)
                      FIELD(IX,KNS) = REAL(ILSM(IX,K)*IOBSRLAT(IX,K,IS),JWRB) * ZCONV  + (1-ILSM(IX,K))*ZMISS
                    ENDDO
                  ENDDO
!$OMP END PARALLEL DO

                CASE(4)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(K,KNS,IX)
                  DO K=1,NGY
                    KNS=NGY-K+1
                    DO IX=1,NLONRGG(K)
                      FIELD(IX,KNS) = REAL(ILSM(IX,K)*IOBSRLON(IX,K,IS),JWRB) * ZCONV  + (1-ILSM(IX,K))*ZMISS
                    ENDDO
                  ENDDO
!$OMP END PARALLEL DO

                CASE(5)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(K,KNS,IX)
                  DO K=1,NGY
                    KNS=NGY-K+1
                    DO IX=1,NLONRGG(K)
                      FIELD(IX,KNS) = REAL(ILSM(IX,K)*IOBSCOR(IX,K,IS),JWRB) * ZCONV  + (1-ILSM(IX,K))*ZMISS
                    ENDDO
                  ENDDO
!$OMP END PARALLEL DO

                END SELECT

              ELSEIF ( IS <= 0 ) THEN
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(K,KNS,IX)
                DO K=1,NGY
                  KNS=NGY-K+1
                  DO IX=1,NLONRGG(K)
                    FIELD(IX,KNS) = ZMISS
                  ENDDO
                ENDDO
!$OMP END PARALLEL DO

              ELSE
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(K,KNS,IX)
                DO K=1,NGY
                  KNS=NGY-K+1
                  DO IX=1,NLONRGG(K)
                    FIELD(IX,KNS) = 1.0_JWRB
                  ENDDO
                ENDDO
!$OMP END PARALLEL DO


              ENDIF 

              ITEST = 0
              ITABLE=140
              IPARAM=KPARAM_SUBGRIG !! use the parameter id of the model bathymetry to insure the same interpolation method between the 2
              IZLEV=0
              ITMIN=0
              ITMAX=0
              CDATE=CDATECLIM
              IFCST=0


              CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,                   &
     &                        ITABLE, IPARAM, IZLEV, ITMIN, ITMAX, IANG , M,  &
     &                        CDATE, IFCST, MARSTYPE, LFDB, IU08(IP))

            ENDDO
          ENDDO

        ELSE
!         BINARY OUTPUT
!         FOR GLOBAL FIELD (in the same file as mean bathymetry)

          WRITE(CX,'(I5.5)') NLONRGG(1)
          FORMAT='('//CX//'I4)'

          DO IS =1,2
            DO K=1,NGY
              DO IXLP = 1,NLONRGG(K),NLONRGG(1)
                WRITE(IU01,FORMAT) (IOBSLAT(IX,K,IS),IX=IXLP,MIN(IXLP+NLONRGG(1)-1,NLONRGG(K)))
              ENDDO
            ENDDO
          ENDDO
          DO IS =1,2
            DO K=1,NGY
              DO IXLP = 1,NLONRGG(K),NLONRGG(1)
                WRITE(IU01,FORMAT) (IOBSLON(IX,K,IS),IX=IXLP,MIN(IXLP+NLONRGG(1)-1,NLONRGG(K)))
              ENDDO
            ENDDO
          ENDDO
          DO IS =1,2
            DO K=1,NGY
              DO IXLP = 1,NLONRGG(K),NLONRGG(1)
                WRITE(IU01,FORMAT) (IOBSRLAT(IX,K,IS),IX=IXLP,MIN(IXLP+NLONRGG(1)-1,NLONRGG(K)))
              ENDDO
            ENDDO
          ENDDO
          DO IS =1,2
            DO K=1,NGY
              DO IXLP = 1,NLONRGG(K),NLONRGG(1)
                WRITE(IU01,FORMAT) (IOBSRLON(IX,K,IS),IX=IXLP,MIN(IXLP+NLONRGG(1)-1,NLONRGG(K)))
              ENDDO
            ENDDO
          ENDDO
          DO IS =1,4
            DO K=1,NGY
              DO IXLP = 1,NLONRGG(K),NLONRGG(1)
                WRITE(IU01,FORMAT) (IOBSCOR(IX,K,IS),IX=IXLP,MIN(IXLP+NLONRGG(1)-1,NLONRGG(K)))
              ENDDO
            ENDDO
          ENDDO

        ENDIF  ! LLGRIBOUT

      ENDDO ! END LOOP ON FREQUENCIES


      IF ( LLGRIBOUT ) THEN
        DO IP = 0, NPROPAGS 
          CALL IGRIB_CLOSE_FILE(IU08(IP))
        ENDDO
      ENDIF

      CALL FLUSH(IU06)

ENDIF  !!! LLOBSTROUT

END PROGRAM CREATE_BATHY_ETOPO1
