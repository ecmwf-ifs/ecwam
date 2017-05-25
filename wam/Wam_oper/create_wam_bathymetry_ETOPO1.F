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
!!!!! SO IF XDELLA < 2*1/60, THEN THE OBSTRUCTIONS WILL ALL BE SET TO NON BLOCKING !!!


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
!     correction_to_wam_grid_xxx
!     where xxxx=int(xdella*100) if grid_description is not present
!     else  xxxx=the gaussian number of latitude - 1 (NY)

!     OUTPUT FILES:

!     wam_topo_xxxx
!     where xxxx=int(xdella*100) if grid_description is not present
!     else  xxxx=the gaussian spectral truncation ISPECTRUNC

!     AND ALL THE FILES NECESSARY FOR PLOTTING WITH METVIEW
!     IF REQUESTED IN THE USER INPUT !
!     meandepth.dat
!     and
!     obstructions_*.dat
 
!**************************************************************************

      IMPLICIT NONE
 
      INTEGER, PARAMETER :: ILON=21601
      INTEGER, PARAMETER :: ILAT=10801
      INTEGER, PARAMETER :: NREF=500
      INTEGER, PARAMETER :: NDPT=1000

      INTEGER :: I_GET_UNIT
      INTEGER :: IU06, IU, IUNIT
      INTEGER :: I, J, IJ, K, KSN, M
      INTEGER :: NX, NY
      INTEGER :: IPER, IRGG, NFRE, IFRE1, ISPECTRUNC
      INTEGER :: NLANDCENTREPM, NLANDCENTREMAX, NLANDCENTRE, NIOBSLAT
      INTEGER :: NSEA, NLAND, NSEASH
      INTEGER :: ILONL, ILONR, ILATB, ILATT
      INTEGER :: NTOT, ICOUNT, IC, IR
      INTEGER :: NREFERENCE
      INTEGER :: IX, NJM, NJP, NIM, NIP, IH
      INTEGER :: II, JJ, IK, NPTS, IDPT
      INTEGER :: IREINF, ITEMPEW
      INTEGER :: IS, KT, KB
      INTEGER :: NOBSTRCT, NIOBSLON, NBLOCKLAND, NTOTPTS
      INTEGER :: INVRES

      INTEGER :: IDUM(15)
      INTEGER, DIMENSION(NREF) :: LEVEL, NDEPTH
      INTEGER, ALLOCATABLE :: NLONRGG(:)
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ITHRSHOLD, IEXCLTHRSHOLD
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IBLOCKDPT
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IDEPTH, WAMDEPTH
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: IOBSLAT, IOBSLON
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: IOBSCOR
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: IOBSRLAT, IOBSRLON

      REAL, PARAMETER :: OLDPI = 3.1415927
      REAL :: AKI
      REAL :: PI, RAD, OLDRAD, G, X60, FRATIO, FR1
      REAL :: XDELLA, XDELLO
      REAL :: AMOSOP, AMONOP, AMOWEP, AMOEAP
      REAL :: ALONL, ALONR, ALATB, ALATT, XLON
      REAL :: PLANDTRHS, PSHALLOWTRHS 
      REAL :: PLONS, XLO, XLA, XI, YJ 
      REAL :: SEA, XLAND, SEASH 
      REAL :: OMEGA, XKDEEP, XKDMAX, XX, DEPTH
      REAL :: STEPT, STEPB, XLATT, XLATB, XLONL, XLONR  
      REAL :: STEPLAT, STEPLON
      REAL :: RESOL
      REAL, DIMENSION(ILON) :: ALON
      REAL, DIMENSION(ILAT) :: ALAT
      REAL, DIMENSION(NDPT) :: XK 
      REAL, DIMENSION(NREF) :: XINF, XSUP, YINF, YSUP
      REAL, ALLOCATABLE, DIMENSION(:) :: ZDELLO,COSPH
      REAL, ALLOCATABLE, DIMENSION(:) :: XLAT
      REAL, ALLOCATABLE, DIMENSION(:) :: FR
      REAL, ALLOCATABLE, DIMENSION(:,:) :: PERCENTLAND, PERCENTSHALLOW

      CHARACTER(LEN=2) :: CFR
      CHARACTER(LEN=5) :: CWAMRESOL
      CHARACTER(LEN=4) :: CX
      CHARACTER(LEN=10) :: FORMAT
      CHARACTER(LEN=32) :: FILENM
      CHARACTER(LEN=72) :: LOCATION(NREF)
      CHARACTER(LEN=144) :: CLINE,FILENAME

      LOGICAL :: LORIGINAL, LLPRINT
      LOGICAL :: LLAND, LREALLAND, L1ST, LNSW
      LOGICAL :: LLGRID
      LOGICAL :: LLWAMDEPTH(NDPT)


      PI=4.*ATAN(1.)
      RAD=PI/180.
      OLDRAD=OLDPI/180.
      G=9.806
      X60=60.

!     ETOPO1 RESOLUTION
      INVRES=60
      RESOL=1./INVRES

      FRATIO=1.1

      IU06=6


!     READ INPUT SELECTION
!     --------------------
      READ(5,*) CLINE
      READ(5,*) CLINE
      READ(5,*) XDELLA
      READ(5,*) CLINE
      READ(5,*) AMOSOP, AMONOP, AMOWEP, AMOEAP
      READ(5,*) CLINE
      READ(5,*) IPER
      READ(5,*) CLINE
      READ(5,*) IRGG
      READ(5,*) CLINE
      READ(5,*) FR1 
      READ(5,*) CLINE
      READ(5,*) NFRE, IFRE1
      READ(5,*) CLINE
      READ(5,*) LLPRINT
      READ(5,*) CLINE
      READ(5,*) ALONL, ALONR, ALATB, ALATT
      READ(5,*) CLINE
      READ(5,*) LORIGINAL 

!     CHECK IF FILE grid_description IS PRESENT
!     IF IT IS THERE IT WILL SUPERSEDE THE OTHER INPUT

      FILENAME='grid_description'
      INQUIRE(FILE=FILENAME,EXIST=LLGRID)
      IF(LLGRID) THEN
        IU=I_GET_UNIT(IU06,FILENAME,'S','F',0)
        OPEN(IU,FILE=FILENAME,STATUS='OLD', FORM='FORMATTED')
        READ (IU,'(I4)') ISPECTRUNC
        READ (IU,'(F8.3)') AMONOP
        READ (IU,'(F8.3)') AMOSOP
        READ (IU,'(F8.3)') AMOWEP
        READ (IU,'(F8.3)') AMOEAP
        READ (IU,'(I4)') IPER
        READ (IU,'(I4)') IRGG
        READ (IU,'(I4)') NY
      ENDIF


      IF(LLGRID) THEN
        XDELLA = (AMONOP-AMOSOP)/(NY-1)
        ALLOCATE(NLONRGG(NY))

        NX = 0
        DO K=1,NY
          KSN=NY-K+1
          READ(IU,'(I4)') NLONRGG(KSN)
          NX = MAX(NX,NLONRGG(KSN))
        ENDDO

        XDELLO = (AMOEAP-AMOWEP)/(NX-1)

        CLOSE(IU)

      ELSE
        XDELLO=XDELLA     
        NX=NINT((AMOEAP-AMOWEP)/XDELLO)+1
        NY=NINT((AMONOP-AMOSOP)/XDELLA)+1

        ALLOCATE(NLONRGG(NY))
      ENDIF

     
  
      NLANDCENTREPM=(NINT(0.2*XDELLA*INVRES)-1)/2
      NLANDCENTREPM=MAX(NLANDCENTREPM,1)
      NLANDCENTREMAX=(2*NLANDCENTREPM+1)**2

      PLANDTRHS=0.3
      PSHALLOWTRHS=0.8

      ALLOCATE(ZDELLO(NY))
      ALLOCATE(COSPH(NY))
      ALLOCATE(XLAT(0:NY+1))

      DO K=0,NY+1
        XLAT(K) = (AMOSOP + REAL(K-1)*XDELLA)
      ENDDO

      DO K=1,NY
!       !!! from south to north !!!!
        COSPH(K)   = COS(XLAT(K)*RAD)
        IF(.NOT.LLGRID) THEN
          IF (IRGG.EQ.1) THEN
            IF(XDELLA.EQ.0.25.AND.AMOWEP.EQ.-98..AND.AMOSOP.EQ.9..AND.
     &         AMOEAP.EQ.42..AND.AMONOP.EQ.81.) THEN
!         the old value for pi has to be taken in order to reproduce
!         exactly the irregular grid of the operational LAW model 0.25
              NLONRGG(K)=NINT(NX*COS((AMOSOP+REAL(K-1)*XDELLA)*OLDRAD))
            ELSE
!            The silly division by cos(x60*RAD) is an attempt at making sure
!            that exactly 0.5 is used for cosine of 60 degrees.
             NLONRGG(K)=
     &         MAX(NINT(NX*(COS(XLAT(K)*RAD)/(2.*COS(X60*RAD)))),2)
            ENDIF
            IF(MOD(NLONRGG(K),2).EQ.1) NLONRGG(K) = NLONRGG(K)+1
          ELSE
            NLONRGG(K) = NX
          ENDIF
          WRITE(IU06,*) 'POINTS PER LATITUDES: ',K, XLAT(K), NLONRGG(K)
        ENDIF

        PLONS=(AMOEAP-AMOWEP) + IPER*XDELLO
        IF(IPER.EQ.1) THEN
          ZDELLO(K)  = PLONS/REAL(NLONRGG(K))
        ELSE
          ZDELLO(K)  = PLONS/REAL(NLONRGG(K)-1)
        ENDIF
      ENDDO

      ALLOCATE(FR(NFRE))

      FR(IFRE1) = FR1
      DO M=IFRE1-1,1,-1
        FR(M) = (FR(M+1)/FRATIO)
      ENDDO
      DO M=IFRE1+1,NFRE
        FR(M) = FRATIO*FR(M-1)
      ENDDO


!     DATASET:
!     --------
 
      IF(ALONL.LT.-180..OR.ALONR.GT.180.) THEN
        WRITE(*,*) ' LONGITUDE SPECIFICATION ERROR +- 180'
        WRITE(*,*) ' ALONL, ALONR : ',ALONL,ALONR
        STOP
      ENDIF
      IF(ALATT.GT.90..OR.ALATB.LT.-90.) THEN
        WRITE(*,*) ' LATITUDE SPECIFICATION ERROR +- 90'
        WRITE(*,*) ' ALATT, ALATB : ',ALATT,ALATB
        STOP
      ENDIF
     
      ILONL = NINT((ALONL + 180.)*INVRES) + 1
      ILONR = NINT((ALONR + 180.)*INVRES) + 1
      ILATB = NINT((90.- ALATB)*INVRES) + 1
      ILATB = MAX(1,MIN(ILATB,ILAT))
      ILATT = NINT((90.- ALATT)*INVRES) + 1
      ILATT = MAX(1,MIN(ILATT,ILAT))
      IF(ILONR.EQ.ILON+1)ILONR=ILON
      IF(ILATB.EQ.ILAT+1)ILATB=ILAT

      ALLOCATE(IDEPTH(ILON,ILAT))
      ALLOCATE(WAMDEPTH(NX,NY))
      ALLOCATE(PERCENTLAND(NX,NY))
      ALLOCATE(PERCENTSHALLOW(NX,NY))

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

      IF(LORIGINAL) THEN
        OPEN(11,file='depth.dat')
        WRITE(11,'(a4)') '#GEO'
        WRITE(11,'(a11)') '#FORMAT LLV'
        WRITE(11,'(a5)') '#DATA'
        write(*,*) 'output of original data set for indices:'
        write(*,*) ILATT,ILATB
        write(*,*) ILONL,ILONR
        DO J=ILATT,ILATB
          DO I=ILONL,ILONR
            IF(IDEPTH(I,J).GE. -300 .and. 
     &         IDEPTH(I,J).LE. 2000 ) THEN
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
         READ(15,*,END=1000,ERR=1000)
     &        XINF(IR),YINF(IR),XSUP(IR),YSUP(IR),LEVEL(IR),NDEPTH(IR),
     &        LOCATION(IR)
      ENDDO


1000  NREFERENCE=IR-1
      WRITE(6,*) 'READ ',NREFERENCE,' NEW REFERENCE LEVELS'

      DO IR=1,NREFERENCE
        WRITE(*,*)
     &        XINF(IR),YINF(IR),XSUP(IR),YSUP(IR),LEVEL(IR),NDEPTH(IR),
     &        LOCATION(IR)
        DO J=1,ILAT
          YJ=ALAT(J)
          IF(YJ.GE.YINF(IR).AND.YJ.LE.YSUP(IR)) THEN
            DO I=1,ILON
              XI=ALON(I)
              IF (XI.GE.XINF(IR).AND.XI.LE.XSUP(IR)) THEN
                 IF(IDEPTH(I,J).LE.LEVEL(IR)) THEN
                   IF(NDEPTH(IR).NE.0) THEN
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

!     SOUTH OF 64S ALL NON DEEP POINTS ARE SET TO LAND TO AVOID
!     THE PROBLEM WITH PERMANENT ICE SHEET
        DO J=1,ILAT
          YJ=ALAT(J)
          IF(YJ.LE.-64.0) THEN
            DO I=1,ILON
              IF(IDEPTH(I,J).GE.-250) THEN
                IDEPTH(I,J)=1
              ENDIF
            ENDDO
          ENDIF
        ENDDO

!     COMPUTE MEAN DEPTH 

      DO K=1,NY
         DO IX=1,NX
           WAMDEPTH(IX,K)=0
         ENDDO
      ENDDO


      NJM=INT(0.5*XDELLA*INVRES)
      NJP=NINT(0.5*XDELLA*INVRES)

      DO K=1,NY
!        WE ASSUME THAT WAMGRID IS ALWAYS WITHIN ETOPO1
!        DETERMINE CLOSEST ETOPO1 J INDEX TO WAM POINT
         DO J=ILAT-1,1,-1
           IF(ALAT(J+1).LT.XLAT(K).AND.
     &        XLAT(K).LE.ALAT(J) ) EXIT
         ENDDO
         J=MIN(MAX(J,1),ILAT)

         DO IX=1,NLONRGG(K)

!          ETOPO1 STARTS AT -180
           XLON=AMOWEP + REAL(IX-1)*ZDELLO(K)
           IF(XLON.GT.180.) THEN
             XLON=XLON-360.
           ENDIF

!          DETERMINE CLOSEST ETOPO1 I INDEX TO WAM POINT
           DO I=1,ILON-1
             IF(ALON(I).LE.XLON.AND.
     &          XLON.LT.ALON(I+1) ) EXIT
           ENDDO

           NIM=INT(0.5*ZDELLO(K)*INVRES)
           NIP=NINT(0.5*ZDELLO(K)*INVRES)

           NSEA=0
           SEA=0.
           NLAND=0
           XLAND=0.
           NSEASH=0
           SEASH=0.

!          AVERAGE OVER LAND AND SEA SEPARATELY
!          AROUND POINT I,J
           DO JJ=J-NJM,J+NJP
             IF(JJ.GE.1 .AND. JJ.LE.ILAT) THEN
               DO II=I-NIM,I+NIP
                 IK=II
                 IF(II.LT.1) IK=ILON+II
                 IF(II.GT.ILON) IK=II-ILON
                 IF(IDEPTH(IK,JJ).LE.0) THEN
                   NSEA=NSEA+1
!                  IN WAM 999M IS THE MAXIMUM DEPTH
                   SEA=SEA+MAX(-999,IDEPTH(IK,JJ))

!                  FIND SHALLOWER AREAS
                   IF(IDEPTH(IK,JJ).GT.-500) THEN
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
             IF(JJ.GE.1 .AND. JJ.LE.ILAT) THEN
               DO II=I-NLANDCENTREPM,I+NLANDCENTREPM
                 IK=II
                 IF(II.LT.1) IK=ILON+II
                 IF(II.GT.ILON) IK=II-ILON
                 IF(IDEPTH(IK,JJ).GT.0) THEN
                   NLANDCENTRE=NLANDCENTRE+1
                 ENDIF
               ENDDO
             ENDIF
           ENDDO

!          IF 40% OR MORE LAND, THEN AVERAGE OVER LAND POINTS
!          OR THE CENTER OF THE GRID BOX IS LAND.
!          ELSE AVERAGE OVER SEA POINTS
           PERCENTLAND(IX,K)=FLOAT(NLAND)/FLOAT(NLAND+NSEA)
           IF(PERCENTLAND(IX,K).GT.0.60 .OR.
     &        NLANDCENTRE.GE.NLANDCENTREMAX ) THEN
             WAMDEPTH(IX,K)=XLAND/NLAND
           ELSE
!            IF THERE IS A PERCENTAGE OF SHALLOWER POINTS THEN
!            THE AVERAGE IS TAKEN OVER THOSE POINTS ALONE.
             PERCENTSHALLOW(IX,K)=FLOAT(NSEASH)/FLOAT(NSEA)
             IF(PERCENTSHALLOW(IX,K).GE.0.3) THEN
               WAMDEPTH(IX,K)=SEASH/NSEASH
             ELSE 
               WAMDEPTH(IX,K)=SEA/NSEA
             ENDIF
           ENDIF

         ENDDO
      ENDDO

      NPTS=0
      DO K=1,NY
         DO IX=1,NLONRGG(K)
          IF(WAMDEPTH(IX,K).LT.0) NPTS=NPTS+1
        ENDDO
      ENDDO

      WRITE(*,*) 'TOTAL NUMBER OF SEA POINT = ',NPTS


!     INTRODUCE CORRECTION TO WAM POINTS.
!     THEY WERE OBTAINED FROM THE PREVIOUS GRID SETUP
!     !!! SEA POINT DEPTHS ARE ALREADY POSITIVE !!!
!     IT's a bit obsolete !!!

      WRITE(*,*) 'CORRECTION TO WAM GRID'

      IF(LLGRID) THEN
        WRITE(CWAMRESOL,'(I5.5)') NY-1
      ELSE
        WRITE(CWAMRESOL,'(I5.5)') INT(1000*XDELLA)
      ENDIF
      FILENAME='correction_to_wam_grid_'//CWAMRESOL

      OPEN(35,FILE=FILENAME)

      DO WHILE(.TRUE.)
         READ(35,*,END=111,ERR=111) XLO,XLA,IX,K,IDPT
         IDPT=-IDPT
         XLON=AMOWEP + REAL(IX-1)*ZDELLO(K)
         IF(ABS(XLON-XLO).GT.ZDELLO(K) .OR.
     &      ABS(XLAT(K)-XLA).GT.XDELLA ) THEN
           WRITE(*,*) 'PROBLEM !!!!'
           WRITE(*,*) 'THE CORRECTION TO WAM GRID IS NOT A WAM POINT'
           WRITE(*,*) XLO,XLA,IDPT 
           WRITE(*,*) XLON,XLAT(K),WAMDEPTH(IX,K)
           WRITE(*,*) ''
         ELSE
           IF(WAMDEPTH(IX,K).GE.0 .AND. IDPT.LT.0) THEN
             NPTS=NPTS+1
           ELSEIF(WAMDEPTH(IX,K).LT.0 .AND. IDPT.GE.0) THEN
             NPTS=NPTS-1
           ENDIF  
           WAMDEPTH(IX,K)=IDPT
         ENDIF
      ENDDO
111   CONTINUE


!     SET MIN AND MAX TO +-999
      DO K=1,NY
         DO IX=1,NLONRGG(K)
           WAMDEPTH(IX,K)=MAX(-999,WAMDEPTH(IX,K))
           WAMDEPTH(IX,K)=MIN(999,WAMDEPTH(IX,K))
        ENDDO
      ENDDO

      WRITE(*,*) 'TOTAL NUMBER OF SEA POINT = ',NPTS
 
!     OUTPUT THE MEAN DATA SET ON SUB AREA
!     OMIT LAND POINTS
!     ------------------------------------
      IF(LLPRINT) THEN
        OPEN(12,file='meandepth.dat')
        WRITE(12,'(a4)') '#GEO'
        WRITE(12,'(a11)') '#FORMAT LLV'
        WRITE(12,'(a5)') '#DATA'
        DO K=1,NY
           DO IX=1,NLONRGG(K)
             XLON=AMOWEP + REAL(IX-1)*ZDELLO(K)
             IF(XLON.GT.180.) then
               XLON=XLON-360.
             ENDIF
            IF(ALATB.LE.XLAT(K) .AND. XLAT(K).LE.ALATT .AND.
     &         ALONL.LE.XLON .AND. XLON.LE.ALONR ) THEN
               IF(WAMDEPTH(IX,K).LT.0 .AND.
     &            WAMDEPTH(IX,K).GT.-999 ) THEN
                  WRITE(12,'(2(1X,F8.3),1X,I4)') 
     &            XLON,XLAT(K),WAMDEPTH(IX,K)
               ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF


!     OUTPUT THE MEAN DATA AS THE INPUT FILE FOR WAM
!     ----------------------------------------------
      IF(LLGRID) THEN
        WRITE(CWAMRESOL,'(I5.5)') ISPECTRUNC
      ELSE
        WRITE(CWAMRESOL,'(I5.5)') INT(1000*XDELLA)
      ENDIF
      FILENAME='wam_topo_'//CWAMRESOL

      OPEN(1,FILE=FILENAME,FORM='FORMATTED')
      WRITE(1,'(6F10.5)')XDELLA,XDELLO,AMOSOP,AMONOP,AMOWEP,AMOEAP
      DO K=1,NY
         WRITE(1,'(I4.4)') NLONRGG(K) 
      ENDDO
      DO K=1,NY
         WRITE(CX,'(I4.4)') NLONRGG(K) 
         FORMAT='('//CX//'I4)'
         WRITE(1,FORMAT) (WAMDEPTH(IX,K),IX=1,NLONRGG(K)) 
      ENDDO

!     FIND WHICH MEAN DEPTH EXITS
      DO IDPT=1,NDPT
        LLWAMDEPTH(IDPT)=.FALSE.
      ENDDO
      DO K=1,NY
         DO IX=1,NLONRGG(K)
           IH=MIN(-WAMDEPTH(IX,K),NDPT)
           IF(IH.GT.0) LLWAMDEPTH(IH)=.TRUE.
         ENDDO
      ENDDO


!     CREATE OBSTRUCTIONS:
!     --------------------

      ALLOCATE(IOBSLAT(NX,NY,2))
      ALLOCATE(IOBSLON(NX,NY,2))
      ALLOCATE(IOBSCOR(NX,NY,4))
      ALLOCATE(IOBSRLAT(NX,NY,2))
      ALLOCATE(IOBSRLON(NX,NY,2))

      WRITE(1,'(I4)') NFRE

!     IREINF IS USED TO REINFORCE LAND OBSTRUCTIONS FOR
!     SMALL GRID SPACING.
      IF(XDELLA.LE.0.5) THEN
        IREINF=2
      ELSE
        IREINF=1
      ENDIF

      ALLOCATE(ITHRSHOLD(NX,NY))
      ALLOCATE(IBLOCKDPT(NX,NY))
      ALLOCATE(IEXCLTHRSHOLD(NX,NY))

      IOBSLAT(:,:,:)=1000
      IOBSLON(:,:,:)=1000
      IOBSCOR(:,:,:)=1000
      IOBSRLAT(:,:,:)=1000
      IOBSRLON(:,:,:)=1000

      DO M=1,NFRE

!!!!    COMPUTE THE OBSTRUCTIONS ONLY WHEN IT IS MEANINGFUL
        IF(XDELLA .LE. 2.*RESOL ) THEN
          IF(M.EQ.1) THEN
            WRITE(*,*) ''
            WRITE(*,*) '*********************************************'
            WRITE(*,*) 'THE REQUESTED RESOLUTION IS SMALL ENOUGH WITH'
            WRITE(*,*) 'RESPECT TO THE INPUT BATHYMETRY DATA !'
            WRITE(*,*) 'NO OBSTRUCTIONS WILL BE COMPUTED !'
            WRITE(*,*) '*********************************************'
            WRITE(*,*) ''
          ENDIF
        ELSE 

!       COMPUTE WAVE NUMBERS
        OMEGA=2*PI*FR(M)
        XKDEEP=OMEGA**2/G
        DO IDPT=1,NDPT
          IF(LLWAMDEPTH(IDPT)) THEN
            DEPTH=FLOAT(IDPT)
            XK(IDPT)=AKI(OMEGA,DEPTH)
          ELSE
            XK(IDPT)=XKDEEP
          ENDIF
        ENDDO


!       COMPUTE THE THRESHOLD AT WHICH THE WAVES
!       ARE OBSTRUCTED BY THE BOTTOM,
!       EXCEPT IF WAM DEPTH OF THE SAME ORDER OF MAGITUDE.
!       ALSO COMPUTE THE DEPTH THAT IS CONSIDERED TO BE FULLY BLOCKING
!       AS IF IT WAS LAND.
        XKDMAX=2.
        DO K=1,NY
          DO IX=1,NLONRGG(K)
            IF(WAMDEPTH(IX,K).LT.0) THEN
              XX=XKDMAX/XK(MIN(-WAMDEPTH(IX,K),NDPT))
              ITHRSHOLD(IX,K)=NINT(-XX)
              IEXCLTHRSHOLD(IX,K)=MAX(10*ITHRSHOLD(IX,K),-998)
              IBLOCKDPT(IX,K)=INT(-0.05*XX)
            ENDIF
          ENDDO
        ENDDO


!       NORTH-SOUTH OBSTRUCTIONS
!       -----------------------
!       IS=1 is for the south-north advection
!       IS=2 is for the north-south advection
        WRITE(*,*) 'CREATE NORTH-SOUTH OBSTRUCTIONS '
        DO IS=1,2
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1)
!$OMP+PRIVATE(K,KT,KB,STEPT,STEPB,XLATT,XLATB,ILATT,ILATB,IX)
!$OMP+PRIVATE(XLONL,XLONR,ILONL,ILONR,NOBSTRCT,NBLOCKLAND)
!$OMP+PRIVATE(I,NIOBSLON,LLAND,LREALLAND,J,LNSW,L1ST)
!$OMP+PRIVATE(NTOTPTS)
          DO K=1,NY
            IF(IS.EQ.1) THEN
              KT=K
              KB=K-1
              STEPT=-RESOL
              STEPB=0.
            ELSE
              KT=K+1
              KB=K
              STEPT=0.
              STEPB=RESOL
            ENDIF
            XLATT=XLAT(KT)+STEPT
            XLATB=XLAT(KB)+STEPB
            ILATT = NINT((90.- XLATT)*INVRES) + 1
            ILATT = MAX(1,MIN(ILATT,ILAT))
            ILATB = NINT((90.- XLATB)*INVRES) + 1
            ILATB = MAX(1,MIN(ILATB,ILAT))
            IF(ILATB.EQ.ILAT+1)ILATB=ILAT

            DO IX=1,NLONRGG(K)
              IF(WAMDEPTH(IX,K).LT.0) THEN
                XLONL=AMOWEP + (REAL(IX-1)-0.5)*ZDELLO(K)
                IF(XLONL.GT.180.) then
                  XLONL=XLONL-360.
                ENDIF
                XLONR=AMOWEP + (REAL(IX-1)+0.5)*ZDELLO(K)
                IF(XLONR.GT.180.) then
                  XLONR=XLONR-360.
                ENDIF


                ILONL = NINT((XLONL + 180.)*INVRES) + 1
                ILONR = NINT((XLONR + 180.)*INVRES) + 1

                NOBSTRCT=0

                IF(ILONL.LE.ILONR) THEN
                  NBLOCKLAND=0
                  DO I=ILONL,ILONR
                    NIOBSLON=0
                    LLAND=.FALSE.
                    LREALLAND=.FALSE.
                    DO J=ILATT,ILATB
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K) ) THEN
!                     IF LAND THEN THE FULL LONGITUDE IS BLOCKED
!                     IF THERE IS A SWITCH BACK TO SEA OR VICE VERSA
!                     (SEE BELOW)
!                     LAND IS DEFINED AS ANYTHING ABOVE IBLOCKDPT(IX,K)
!                     ------------------------------------------
                        IF(IDEPTH(I,J).GE.0 ) LREALLAND=.TRUE. 
                        LLAND=.TRUE.
                        NIOBSLON=NIOBSLON+1 
                      ELSEIF (IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                        WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
!                     IF SEA ABOVE THE THRESHOLD THEN ONLY THAT
!                     GRID POINTS BLOCKS
!                     ------------------------------------------
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO

                    IF(LLAND) THEN
                      IF(LREALLAND) THEN
                        LNSW=.TRUE.  
                        IF(IDEPTH(I,ILATT).GE.IBLOCKDPT(IX,K)) THEN
                          L1ST=.TRUE.
                        ELSE
                          L1ST=.FALSE.
                        ENDIF
                        DO J=ILATT+1,ILATB
                          IF( ((IDEPTH(I,J).GE.IBLOCKDPT(IX,K))
     &                         .NEQV.L1ST)
     &                         .AND.
     &                         LNSW ) THEN
                            LNSW=.FALSE.
                          ENDIF
                          IF(((IDEPTH(I,J).GE.IBLOCKDPT(IX,K)).EQV.L1ST)
     &                        .AND.
     &                       .NOT. LNSW ) THEN
!                           LAND IS BLOCKING
                            NIOBSLON=IREINF*(ILATB-ILATT+1)
                            NBLOCKLAND=NBLOCKLAND+1
                            EXIT
                          ENDIF
                        ENDDO
                        IF(LNSW) NIOBSLON=ILATB-ILATT+1
                      ELSE
                        IF(PERCENTSHALLOW(IX,K).GT.PSHALLOWTRHS) THEN
!                         mostly shallow, do not enhance obstruction
                          NIOBSLON=ILATB-ILATT+1
                        ELSEIF(PERCENTLAND(IX,K).LT.PLANDTRHS) THEN
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
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR-ILONL+1)+
     &                    (IREINF-1)*NBLOCKLAND*(ILATB-ILATT+1)

                  IOBSLAT(IX,K,IS)=
     &               NINT((1.-FLOAT(NOBSTRCT)/NTOTPTS)*1000)
                ELSE
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR+ILON-ILONL+1)
                  DO I=1,ILONR
                    NIOBSLON=0
                    DO J=ILATT,ILATB
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K)) THEN
                        NIOBSLON=ILATB-ILATT+1
                        EXIT
                      ELSEIF(IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                       WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO
                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO
                  DO I=ILONL,ILON
                    NIOBSLON=0
                    DO J=ILATT,ILATB
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K)) THEN
                        NIOBSLON=ILATB-ILATT+1
                        EXIT
                      ELSEIF(IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                       WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO
                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO
                  IOBSLAT(IX,K,IS)=
     &               NINT((1.-FLOAT(NOBSTRCT)/NTOTPTS)*1000)
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
        WRITE(*,*) 'CREATE EAST-WEST  OBSTRUCTIONS '
        DO IS=1,2
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1)
!$OMP+PRIVATE(K,XLATT,XLATB,ILATT,ILATB,IX)
!$OMP+PRIVATE(XLONL,XLONR,ILONL,ILONR,NOBSTRCT,NBLOCKLAND)
!$OMP+PRIVATE(J,NIOBSLAT,LLAND,LREALLAND,I,LNSW,L1ST)
!$OMP+PRIVATE(NTOTPTS)
          DO K=1,NY
            XLATT=XLAT(K)+0.5*XDELLA
            XLATB=XLAT(K)-0.5*XDELLA
            ILATT = NINT((90.- XLATT)*INVRES) + 1
            ILATT = MAX(1,MIN(ILATT,ILAT))
            ILATB = NINT((90.- XLATB)*INVRES) + 1
            ILATB = MAX(1,MIN(ILATB,ILAT))
            IF(ILATB.EQ.ILAT+1)ILATB=ILAT

            DO IX=1,NLONRGG(K)
              IF(WAMDEPTH(IX,K).LT.0) THEN
                IF(IS.EQ.1) THEN
                  XLONL=AMOWEP + (REAL(IX-2))*ZDELLO(K)
                  XLONR=AMOWEP + (REAL(IX-1))*ZDELLO(K) -RESOL
                ELSE
                  XLONL=AMOWEP + (REAL(IX-1))*ZDELLO(K) +RESOL
                  XLONR=AMOWEP + (REAL(IX))*ZDELLO(K)
                ENDIF
                IF(XLONL.GT.180.) THEN
                  XLONL=XLONL-360.
                ENDIF
                IF(XLONR.GT.180.) THEN
                  XLONR=XLONR-360.
                ENDIF

                ILONL = NINT((XLONL + 180.)*INVRES) + 1
                ILONR = NINT((XLONR + 180.)*INVRES) + 1

                NOBSTRCT=0

                IF(ILONL.LE.ILONR) THEN
                  NBLOCKLAND=0
                  DO J=ILATT,ILATB
                    NIOBSLAT=0
                    LLAND=.FALSE.
                    LREALLAND=.FALSE.
                    DO I=ILONL,ILONR
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K) ) THEN
!                     IF LAND THEN THE FULL LONGITUDE IS BLOCKED
!                     IF THERE IS A SWITCH BACK TO SEA OR VICE VERSA
!                     (SEE BELOW)
!                     LAND IS DEFINED AS ANYTHING ABOVE IBLOCKDPT(IX,K)
!                     ------------------------------------------
                        LLAND=.TRUE.
                        IF(IDEPTH(I,J).GE.0 ) LREALLAND=.TRUE. 
                        NIOBSLAT=NIOBSLAT+1 
                      ELSEIF (IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                        WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
!                     IF SEA ABOVE THE THRESHOLD THEN ONLY THAT
!                     GRID POINTS BLOCKS
!                     ------------------------------------------
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO

                    IF(LLAND) THEN
                      IF(LREALLAND) THEN
                        LNSW=.TRUE.  
                        IF(IDEPTH(ILONL,J).GE.IBLOCKDPT(IX,K) ) THEN
                          L1ST=.TRUE.
                        ELSE
                          L1ST=.FALSE.
                        ENDIF
                        DO I=ILONL+1,ILONR
                          IF( ((IDEPTH(I,J).GE.IBLOCKDPT(IX,K))
     &                        .NEQV.L1ST)
     &                         .AND.
     &                         LNSW ) THEN
                            LNSW=.FALSE.
                          ENDIF
                          IF(((IDEPTH(I,J).GE.IBLOCKDPT(IX,K)).EQV.L1ST)
     &                        .AND.
     &                       .NOT. LNSW ) THEN
!                           LAND IS BLOCKING
                            NIOBSLAT=IREINF*(ILONR-ILONL+1)
                            NBLOCKLAND=NBLOCKLAND+1
                            EXIT
                          ENDIF
                        ENDDO
                        IF(LNSW) NIOBSLAT=ILONR-ILONL+1
                      ELSE
                        IF(PERCENTSHALLOW(IX,K).GT.PSHALLOWTRHS) THEN
                          NIOBSLAT=ILONR-ILONL+1
                        ELSEIF(PERCENTLAND(IX,K).LT.PLANDTRHS) THEN
                          NIOBSLAT=IREINF*(ILONR-ILONL+1)
                          NBLOCKLAND=NBLOCKLAND+1
                        ELSE
                          NIOBSLAT=0
                        ENDIF
                      ENDIF
                    ENDIF

                    NOBSTRCT=NOBSTRCT+NIOBSLAT
                  ENDDO

                  NTOTPTS=(ILATB-ILATT+1)*(ILONR-ILONL+1)+
     &                    (IREINF-1)*NBLOCKLAND*(ILONR-ILONL+1)
                  IOBSLON(IX,K,IS)=
     &               NINT((1.-FLOAT(NOBSTRCT)/NTOTPTS)*1000)
                ELSE
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR+ILON-ILONL+1)
                  DO J=ILATT,ILATB
                    NIOBSLAT=0
                    DO I=1,ILONR
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K)) THEN
                        NIOBSLAT=ILONR+ILON-ILONL+1
                        GOTO 1111 
                      ELSEIF(IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                       WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO
                    DO I=ILONL,ILON
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K)) THEN
                        NIOBSLAT=ILONR+ILON-ILONL+1
                        EXIT
                      ELSEIF(IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                       WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO
1111                CONTINUE
                    NOBSTRCT=NOBSTRCT+NIOBSLAT
                  ENDDO
                  IOBSLON(IX,K,IS)=
     &               NINT((1.-FLOAT(NOBSTRCT)/NTOTPTS)*1000)
                ENDIF

              ENDIF
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ENDDO


!       NORTH-WEST-SOUTH-EAST OBSTRUCTIONS
!       ----------------------------------
!       IS=1 is for the southeast-northwest advection
!       IS=2 is for the northwest-southeast advection

        WRITE(*,*) 'CREATE NORTH-WEST-SOUTH-EAST OBSTRUCTIONS '

!       first search north-south
        DO IS=1,2
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1)
!$OMP+PRIVATE(K,KT,KB,STEPT,STEPB,XLATT,XLATB,ILATT,ILATB,IX)
!$OMP+PRIVATE(XLON,XLONL,XLONR,ILONL,ILONR,NOBSTRCT,NBLOCKLAND)
!$OMP+PRIVATE(I,NIOBSLON,LLAND,LREALLAND,J,LNSW,L1ST)
!$OMP+PRIVATE(NTOTPTS)
          DO K=1,NY
            IF(IS.EQ.1) THEN
              KT=K
              KB=K-1
              STEPT=-RESOL
              STEPB=0.
            ELSE
              KT=K+1
              KB=K
              STEPT=0.
              STEPB=RESOL
            ENDIF
            XLATT=XLAT(KT)+STEPT
            XLATB=XLAT(KB)+STEPB
            ILATT = NINT((90.- XLATT)*INVRES) + 1
            ILATT = MAX(1,MIN(ILATT,ILAT))
            ILATB = NINT((90.- XLATB)*INVRES) + 1
            ILATB = MAX(1,MIN(ILATB,ILAT))
            IF(ILATB.EQ.ILAT+1)ILATB=ILAT

            DO IX=1,NLONRGG(K)
              IF(WAMDEPTH(IX,K).LT.0) THEN
                XLON=AMOWEP + REAL(IX-1)*ZDELLO(K)
                XLONL=XLON -(IS-1)*XDELLA
                IF(XLONL.GT.180.) THEN
                  XLONL=XLONL-360.
                ENDIF
                XLONR=XLON +(2-IS)*XDELLA
                IF(XLONR.GT.180.) THEN
                  XLONR=XLONR-360.
                ENDIF

                ILONL = NINT((XLONL + 180.)*INVRES) + 1
                ILONR = NINT((XLONR + 180.)*INVRES) + 1

                NOBSTRCT=0

                IF(ILONL.LE.ILONR) THEN
                  NBLOCKLAND=0
                  DO I=ILONL,ILONR
                    NIOBSLON=0
                    LLAND=.FALSE.
                    LREALLAND=.FALSE.
                    DO J=ILATT,ILATB
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K) ) THEN
!                     IF LAND THEN THE FULL LONGITUDE IS BLOCKED
!                     IF THERE IS A SWITCH BACK TO SEA OR VICE VERSA
!                     (SEE BELOW)
!                     LAND IS DEFINED AS ANYTHING ABOVE IBLOCKDPT(IX,K)
!                     ------------------------------------------
                        IF(IDEPTH(I,J).GE.0 ) LREALLAND=.TRUE. 
                        LLAND=.TRUE.
                        NIOBSLON=NIOBSLON+1 
                      ELSEIF (IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                        WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
!                     IF SEA ABOVE THE THRESHOLD THEN ONLY THAT
!                     GRID POINTS BLOCKS
!                     ------------------------------------------
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO

                    IF(LLAND) THEN
                      IF(LREALLAND) THEN
                        LNSW=.TRUE.  
                        IF(IDEPTH(I,ILATT).GE.IBLOCKDPT(IX,K)) THEN
                          L1ST=.TRUE.
                        ELSE
                          L1ST=.FALSE.
                        ENDIF
                        DO J=ILATT+1,ILATB
                          IF( ((IDEPTH(I,J).GE.IBLOCKDPT(IX,K))
     &                        .NEQV.L1ST)
     &                         .AND.
     &                         LNSW ) THEN
                            LNSW=.FALSE.
                          ENDIF
                          IF(((IDEPTH(I,J).GE.IBLOCKDPT(IX,K)).EQV.L1ST)
     &                        .AND.
     &                       .NOT. LNSW ) THEN
!                           LAND IS BLOCKING
                            NIOBSLON=IREINF*(ILATB-ILATT+1)
                            NBLOCKLAND=NBLOCKLAND+1
                            EXIT
                          ENDIF
                        ENDDO
                        IF(LNSW) NIOBSLON=ILATB-ILATT+1
                      ELSE
                        IF(PERCENTSHALLOW(IX,K).GT.PSHALLOWTRHS) THEN
                          NIOBSLON=ILATB-ILATT+1
                        ELSEIF(PERCENTLAND(IX,K).LT.PLANDTRHS) THEN
                          NIOBSLON=IREINF*(ILATB-ILATT+1)
                          NBLOCKLAND=NBLOCKLAND+1
                        ELSE
                          NIOBSLON=0
                        ENDIF
                      ENDIF
                    ENDIF

                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR-ILONL+1)+
     &                    (IREINF-1)*NBLOCKLAND*(ILATB-ILATT+1)

                  IOBSRLAT(IX,K,IS)=
     &               NINT((1.-FLOAT(NOBSTRCT)/NTOTPTS)*1000)
                ELSE
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR+ILON-ILONL+1)
                  DO I=1,ILONR
                    NIOBSLON=0
                    DO J=ILATT,ILATB
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K)) THEN
                        NIOBSLON=ILATB-ILATT+1
                        EXIT
                      ELSEIF(IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                       WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO
                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO
                  DO I=ILONL,ILON
                    NIOBSLON=0
                    DO J=ILATT,ILATB
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K)) THEN
                        NIOBSLON=ILATB-ILATT+1
                        EXIT
                      ELSEIF(IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                       WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO
                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO
                  IOBSRLAT(IX,K,IS)=
     &               NINT((1.-FLOAT(NOBSTRCT)/NTOTPTS)*1000)
                ENDIF

              ENDIF
            ENDDO

          ENDDO
!$OMP END PARALLEL DO
        ENDDO

!       then search east-west and average 
        DO IS=1,2
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1)
!$OMP+PRIVATE(K,XLATT,XLATB,ILATT,ILATB,IX)
!$OMP+PRIVATE(XLON,XLONL,XLONR,ILONL,ILONR,NOBSTRCT,NBLOCKLAND)
!$OMP+PRIVATE(J,NIOBSLAT,LLAND,LREALLAND,I,LNSW,L1ST)
!$OMP+PRIVATE(NTOTPTS,ITEMPEW,XX)
          DO K=1,NY
            XLATT=XLAT(K)+(IS-1)*XDELLA
            XLATB=XLAT(K)-(2-IS)*XDELLA
            ILATT = NINT((90.- XLATT)*INVRES) + 1
            ILATT = MAX(1,MIN(ILATT,ILAT))
            ILATB = NINT((90.- XLATB)*INVRES) + 1
            ILATB = MAX(1,MIN(ILATB,ILAT))
            IF(ILATB.EQ.ILAT+1)ILATB=ILAT

            DO IX=1,NLONRGG(K)
              IF(WAMDEPTH(IX,K).LT.0) THEN
                XLON=AMOWEP + REAL(IX-1)*ZDELLO(K)
                IF(IS.EQ.1) THEN
                  XLONL=XLON + RESOL
                  XLONR=XLON + XDELLA
                ELSE
                  XLONL=XLON -  XDELLA
                  XLONR=XLON - RESOL
                ENDIF
                IF(XLONL.GT.180.) THEN
                  XLONL=XLONL-360.
                ENDIF
                IF(XLONR.GT.180.) THEN
                  XLONR=XLONR-360.
                ENDIF

                ILONL = NINT((XLONL + 180.)*INVRES) + 1
                ILONR = NINT((XLONR + 180.)*INVRES) + 1

                NOBSTRCT=0

                IF(ILONL.LE.ILONR) THEN
                  NBLOCKLAND=0
                  DO J=ILATT,ILATB
                    NIOBSLAT=0
                    LLAND=.FALSE.
                    LREALLAND=.FALSE.
                    DO I=ILONL,ILONR
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K) ) THEN
!                     IF LAND THEN THE FULL LONGITUDE IS BLOCKED
!                     IF THERE IS A SWITCH BACK TO SEA OR VICE VERSA
!                     (SEE BELOW)
!                     LAND IS DEFINED AS ANYTHING ABOVE IBLOCKDPT(IX,K)
!                     ------------------------------------------
                        LLAND=.TRUE.
                        IF(IDEPTH(I,J).GE.0 ) LREALLAND=.TRUE. 
                        NIOBSLAT=NIOBSLAT+1 
                      ELSEIF (IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                        WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
!                     IF SEA ABOVE THE THRESHOLD THEN ONLY THAT
!                     GRID POINTS BLOCKS
!                     ------------------------------------------
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO

                    IF(LLAND) THEN
                      IF(LREALLAND) THEN
                        LNSW=.TRUE.  
                        IF(IDEPTH(ILONL,J).GE.IBLOCKDPT(IX,K) ) THEN
                          L1ST=.TRUE.
                        ELSE
                          L1ST=.FALSE.
                        ENDIF
                        DO I=ILONL+1,ILONR
                          IF( ((IDEPTH(I,J).GE.IBLOCKDPT(IX,K))
     &                         .NEQV.L1ST)
     &                         .AND.
     &                         LNSW ) THEN
                            LNSW=.FALSE.
                          ENDIF
                          IF(((IDEPTH(I,J).GE.IBLOCKDPT(IX,K)).EQV.L1ST)
     &                        .AND.
     &                       .NOT. LNSW ) THEN
!                           LAND IS BLOCKING
                            NIOBSLAT=IREINF*(ILONR-ILONL+1)
                            NBLOCKLAND=NBLOCKLAND+1
                            EXIT
                          ENDIF
                        ENDDO
                        IF(LNSW) NIOBSLAT=ILONR-ILONL+1
                      ELSE
                        IF(PERCENTSHALLOW(IX,K).GT.PSHALLOWTRHS) THEN
                          NIOBSLAT=ILONR-ILONL+1
                        ELSEIF(PERCENTLAND(IX,K).LT.PLANDTRHS) THEN
                          NIOBSLAT=IREINF*(ILONR-ILONL+1)
                          NBLOCKLAND=NBLOCKLAND+1
                        ELSE
                          NIOBSLAT=0
                        ENDIF
                      ENDIF
                    ENDIF

                    NOBSTRCT=NOBSTRCT+NIOBSLAT
                  ENDDO

                  NTOTPTS=(ILATB-ILATT+1)*(ILONR-ILONL+1)+
     &                    (IREINF-1)*NBLOCKLAND*(ILONR-ILONL+1)
                  ITEMPEW=NINT((1.-FLOAT(NOBSTRCT)/NTOTPTS)*1000)
                ELSE
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR+ILON-ILONL+1)
                  DO J=ILATT,ILATB
                    NIOBSLAT=0
                    DO I=1,ILONR
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K)) THEN
                        NIOBSLAT=ILONR+ILON-ILONL+1
                        GOTO 2222 
                      ELSEIF(IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                       WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO
                    DO I=ILONL,ILON
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K)) THEN
                        NIOBSLAT=ILONR+ILON-ILONL+1
                        EXIT
                      ELSEIF(IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                       WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO
2222                CONTINUE
                    NOBSTRCT=NOBSTRCT+NIOBSLAT
                  ENDDO
                  ITEMPEW=NINT((1.-FLOAT(NOBSTRCT)/NTOTPTS)*1000)
                ENDIF
                XX=FLOAT(IOBSRLAT(IX,K,IS)*ITEMPEW)
                XX=SQRT(XX)
                IOBSRLAT(IX,K,IS)=MIN(NINT(XX),1000)

              ENDIF
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ENDDO


!       SOUTH-WEST-NORTH-EAST OBSTRUCTIONS
!       ----------------------------------
!       IS=1 is for the southwest-northeast advection
!       IS=2 is for the northeast-southwest advection

        WRITE(*,*) 'CREATE SOUTH-WEST-NORTH-EAST OBSTRUCTIONS '

!       first search north-south
        DO IS=1,2
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1)
!$OMP+PRIVATE(K,KT,KB,STEPT,STEPB,XLATT,XLATB,ILATT,ILATB,IX)
!$OMP+PRIVATE(XLON,XLONL,XLONR,ILONL,ILONR,NOBSTRCT,NBLOCKLAND)
!$OMP+PRIVATE(I,NIOBSLON,LLAND,LREALLAND,J,LNSW,L1ST)
!$OMP+PRIVATE(NTOTPTS)
          DO K=1,NY
            IF(IS.EQ.1) THEN
              KT=K
              KB=K-1
              STEPT=-RESOL
              STEPB=0.
            ELSE
              KT=K+1
              KB=K
              STEPT=0.
              STEPB=RESOL
            ENDIF
            XLATT=XLAT(KT)+STEPT
            XLATB=XLAT(KB)+STEPB
            ILATT = NINT((90.- XLATT)*INVRES) + 1
            ILATT = MAX(1,MIN(ILATT,ILAT))
            ILATB = NINT((90.- XLATB)*INVRES) + 1
            ILATB = MAX(1,MIN(ILATB,ILAT))
            IF(ILATB.EQ.ILAT+1)ILATB=ILAT

            DO IX=1,NLONRGG(K)
              IF(WAMDEPTH(IX,K).LT.0) THEN
                XLON=AMOWEP + REAL(IX-1)*ZDELLO(K)
                XLONL=XLON -(2-IS)*XDELLA
                IF(XLONL.GT.180.) THEN
                  XLONL=XLONL-360.
                ENDIF
                XLONR=XLON +(IS-1)*XDELLA
                IF(XLONR.GT.180.) THEN
                  XLONR=XLONR-360.
                ENDIF

                ILONL = NINT((XLONL + 180.)*INVRES) + 1
                ILONR = NINT((XLONR + 180.)*INVRES) + 1

                NOBSTRCT=0

                IF(ILONL.LE.ILONR) THEN
                  NBLOCKLAND=0
                  DO I=ILONL,ILONR
                    NIOBSLON=0
                    LLAND=.FALSE.
                    LREALLAND=.FALSE.
                    DO J=ILATT,ILATB
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K) ) THEN
!                     IF LAND THEN THE FULL LONGITUDE IS BLOCKED
!                     IF THERE IS A SWITCH BACK TO SEA OR VICE VERSA
!                     (SEE BELOW)
!                     LAND IS DEFINED AS ANYTHING ABOVE IBLOCKDPT(IX,K)
!                     ------------------------------------------
                        IF(IDEPTH(I,J).GE.0 ) LREALLAND=.TRUE. 
                        LLAND=.TRUE.
                        NIOBSLON=NIOBSLON+1 
                      ELSEIF (IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                        WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
!                     IF SEA ABOVE THE THRESHOLD THEN ONLY THAT
!                     GRID POINTS BLOCKS
!                     ------------------------------------------
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO

                    IF(LLAND) THEN
                      IF(LREALLAND) THEN
                        LNSW=.TRUE.  
                        IF(IDEPTH(I,ILATT).GE.IBLOCKDPT(IX,K)) THEN
                          L1ST=.TRUE.
                        ELSE
                          L1ST=.FALSE.
                        ENDIF
                        DO J=ILATT+1,ILATB
                          IF( ((IDEPTH(I,J).GE.IBLOCKDPT(IX,K))
     &                        .NEQV.L1ST)
     &                         .AND.
     &                         LNSW ) THEN
                            LNSW=.FALSE.
                          ENDIF
                          IF(((IDEPTH(I,J).GE.IBLOCKDPT(IX,K)).EQV.L1ST)
     &                        .AND.
     &                       .NOT. LNSW ) THEN
!                           LAND IS BLOCKING
                            NIOBSLON=IREINF*(ILATB-ILATT+1)
                            NBLOCKLAND=NBLOCKLAND+1
                            EXIT
                          ENDIF
                        ENDDO
                        IF(LNSW) NIOBSLON=ILATB-ILATT+1
                      ELSE
                        IF(PERCENTSHALLOW(IX,K).GT.PSHALLOWTRHS) THEN
                          NIOBSLON=ILATB-ILATT+1
                        ELSEIF(PERCENTLAND(IX,K).LT.PLANDTRHS) THEN
                          NIOBSLON=IREINF*(ILATB-ILATT+1)
                          NBLOCKLAND=NBLOCKLAND+1
                        ELSE
                          NIOBSLON=0
                        ENDIF
                      ENDIF
                    ENDIF

                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR-ILONL+1)+
     &                    (IREINF-1)*NBLOCKLAND*(ILATB-ILATT+1)

                  IOBSRLON(IX,K,IS)=
     &               NINT((1.-FLOAT(NOBSTRCT)/NTOTPTS)*1000)
                ELSE
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR+ILON-ILONL+1)
                  DO I=1,ILONR
                    NIOBSLON=0
                    DO J=ILATT,ILATB
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K)) THEN
                        NIOBSLON=ILATB-ILATT+1
                        EXIT
                      ELSEIF(IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                       WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO
                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO
                  DO I=ILONL,ILON
                    NIOBSLON=0
                    DO J=ILATT,ILATB
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K)) THEN
                        NIOBSLON=ILATB-ILATT+1
                        EXIT
                      ELSEIF(IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                       WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO
                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO
                  IOBSRLON(IX,K,IS)=
     &               NINT((1.-FLOAT(NOBSTRCT)/NTOTPTS)*1000)
                ENDIF

              ENDIF
            ENDDO

          ENDDO
!$OMP END PARALLEL DO
        ENDDO

!       then search east-west and average 
        DO IS=1,2
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1)
!$OMP+PRIVATE(K,XLATT,XLATB,ILATT,ILATB,IX)
!$OMP+PRIVATE(XLON,XLONL,XLONR,ILONL,ILONR,NOBSTRCT,NBLOCKLAND)
!$OMP+PRIVATE(J,NIOBSLAT,LLAND,LREALLAND,I,LNSW,L1ST)
!$OMP+PRIVATE(NTOTPTS,ITEMPEW,XX)
          DO K=1,NY
            XLATT=XLAT(K)+(IS-1)*XDELLA
            XLATB=XLAT(K)-(2-IS)*XDELLA
            ILATT = NINT((90.- XLATT)*INVRES) + 1
            ILATT = MAX(1,MIN(ILATT,ILAT))
            ILATB = NINT((90.- XLATB)*INVRES) + 1
            ILATB = MAX(1,MIN(ILATB,ILAT))
            IF(ILATB.EQ.ILAT+1)ILATB=ILAT

            DO IX=1,NLONRGG(K)
              IF(WAMDEPTH(IX,K).LT.0) THEN
                XLON=AMOWEP + REAL(IX-1)*ZDELLO(K)
                IF(IS.EQ.1) THEN
                  XLONL=XLON + RESOL
                  XLONR=XLON + XDELLA
                ELSE
                  XLONL=XLON - XDELLA
                  XLONR=XLON - RESOL
                ENDIF
                IF(XLONL.GT.180.) THEN
                  XLONL=XLONL-360.
                ENDIF
                IF(XLONR.GT.180.) THEN
                  XLONR=XLONR-360.
                ENDIF

                ILONL = NINT((XLONL + 180.)*INVRES) + 1
                ILONR = NINT((XLONR + 180.)*INVRES) + 1

                NOBSTRCT=0

                IF(ILONL.LE.ILONR) THEN
                  NBLOCKLAND=0
                  DO J=ILATT,ILATB
                    NIOBSLAT=0
                    LLAND=.FALSE.
                    LREALLAND=.FALSE.
                    DO I=ILONL,ILONR
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K) ) THEN
!                     IF LAND THEN THE FULL LONGITUDE IS BLOCKED
!                     IF THERE IS A SWITCH BACK TO SEA OR VICE VERSA
!                     (SEE BELOW)
!                     LAND IS DEFINED AS ANYTHING ABOVE IBLOCKDPT(IX,K)
!                     ------------------------------------------
                        LLAND=.TRUE.
                        IF(IDEPTH(I,J).GE.0 ) LREALLAND=.TRUE. 
                        NIOBSLAT=NIOBSLAT+1 
                      ELSEIF (IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                        WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
!                     IF SEA ABOVE THE THRESHOLD THEN ONLY THAT
!                     GRID POINTS BLOCKS
!                     ------------------------------------------
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO

                    IF(LLAND) THEN
                      IF(LREALLAND) THEN
                        LNSW=.TRUE.  
                        IF(IDEPTH(ILONL,J).GE.IBLOCKDPT(IX,K) ) THEN
                          L1ST=.TRUE.
                        ELSE
                          L1ST=.FALSE.
                        ENDIF
                        DO I=ILONL+1,ILONR
                          IF( ((IDEPTH(I,J).GE.IBLOCKDPT(IX,K))
     &                         .NEQV.L1ST)
     &                         .AND.
     &                         LNSW ) THEN
                            LNSW=.FALSE.
                          ENDIF
                          IF(((IDEPTH(I,J).GE.IBLOCKDPT(IX,K)).EQV.L1ST)
     &                        .AND.
     &                       .NOT. LNSW ) THEN
!                           LAND IS BLOCKING
                            NIOBSLAT=IREINF*(ILONR-ILONL+1)
                            NBLOCKLAND=NBLOCKLAND+1
                            EXIT
                          ENDIF
                        ENDDO
                        IF(LNSW) NIOBSLAT=ILONR-ILONL+1
                      ELSE
                        IF(PERCENTSHALLOW(IX,K).GT.PSHALLOWTRHS) THEN
                          NIOBSLAT=ILONR-ILONL+1
                        ELSEIF(PERCENTLAND(IX,K).LT.PLANDTRHS) THEN
                          NIOBSLAT=IREINF*(ILONR-ILONL+1)
                          NBLOCKLAND=NBLOCKLAND+1
                        ELSE
                          NIOBSLAT=0
                        ENDIF
                      ENDIF
                    ENDIF

                    NOBSTRCT=NOBSTRCT+NIOBSLAT
                  ENDDO

                  NTOTPTS=(ILATB-ILATT+1)*(ILONR-ILONL+1)+
     &                    (IREINF-1)*NBLOCKLAND*(ILONR-ILONL+1)
                  ITEMPEW=NINT((1.-FLOAT(NOBSTRCT)/NTOTPTS)*1000)
                ELSE
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR+ILON-ILONL+1)
                  DO J=ILATT,ILATB
                    NIOBSLAT=0
                    DO I=1,ILONR
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K)) THEN
                        NIOBSLAT=ILONR+ILON-ILONL+1
                        GOTO 3333 
                      ELSEIF(IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                       WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO
                    DO I=ILONL,ILON
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K)) THEN
                        NIOBSLAT=ILONR+ILON-ILONL+1
                        EXIT
                      ELSEIF(IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                       WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO
3333                CONTINUE
                    NOBSTRCT=NOBSTRCT+NIOBSLAT
                  ENDDO
                  ITEMPEW=NINT((1.-FLOAT(NOBSTRCT)/NTOTPTS)*1000)
                ENDIF
                XX=FLOAT(IOBSRLON(IX,K,IS)*ITEMPEW)
                XX=SQRT(XX)
                IOBSRLON(IX,K,IS)=MIN(NINT(XX),1000)

              ENDIF
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ENDDO


!       GRID CORNER POINT OBSTRUCTIONS
!       ------------------------------
!       IS=1 is for the northeast-southwest advection
!       IS=2 is for the southeast-northwest advection
!       IS=3 is for the southwest-northeast advection
!       IS=4 is for the northwest-southeast advection

        WRITE(*,*) 'CREATE GRID CORNER OBSTRUCTIONS '

!       first search north-south
        DO IS=1,4
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1)
!$OMP+PRIVATE(K,KT,KB,STEPT,STEPB,XLATT,XLATB,ILATT,ILATB,IX)
!$OMP+PRIVATE(XLON,XLONL,XLONR,ILONL,ILONR,NOBSTRCT,NBLOCKLAND)
!$OMP+PRIVATE(I,NIOBSLON,LLAND,LREALLAND,J,LNSW,L1ST)
!$OMP+PRIVATE(NTOTPTS)
          DO K=1,NY
            IF(IS.EQ.1) THEN 
              KT=K+1
              KB=K
              STEPT=0.
              STEPB=RESOL
            ELSE IF(IS.EQ.2) THEN
              KT=K
              KB=K-1
              STEPT=-RESOL
              STEPB=0.
            ELSE IF(IS.EQ.3) THEN
              KT=K
              KB=K-1
              STEPT=-RESOL
              STEPB=0.
            ELSE IF (IS.EQ.4) THEN
              KT=K+1
              KB=K
              STEPT=0.
              STEPB=RESOL
            ENDIF

            XLATT=XLAT(KT)+STEPT
            XLATB=XLAT(KB)+STEPB
            ILATT = NINT((90.- XLATT)*INVRES) + 1
            ILATT = MAX(1,MIN(ILATT,ILAT))
            ILATB = NINT((90.- XLATB)*INVRES) + 1
            ILATB = MAX(1,MIN(ILATB,ILAT))
            IF(ILATB.EQ.ILAT+1)ILATB=ILAT

            DO IX=1,NLONRGG(K)
              IF(WAMDEPTH(IX,K).LT.0) THEN
                XLON=AMOWEP + REAL(IX-1)*ZDELLO(K)
                XLONL=XLON -((IS-1)/2)*ZDELLO(K)
                IF(XLONL.GT.180.) THEN
                  XLONL=XLONL-360.
                ENDIF
                XLONR=XLON +((4-IS)/2)*ZDELLO(K)
                IF(XLONR.GT.180.) THEN
                  XLONR=XLONR-360.
                ENDIF

                ILONL = NINT((XLONL + 180.)*INVRES) + 1
                ILONR = NINT((XLONR + 180.)*INVRES) + 1

                NOBSTRCT=0

                IF(ILONL.LE.ILONR) THEN
                  NBLOCKLAND=0
                  DO I=ILONL,ILONR
                    NIOBSLON=0
                    LLAND=.FALSE.
                    LREALLAND=.FALSE.
                    DO J=ILATT,ILATB
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K) ) THEN
!                     IF LAND THEN THE FULL LONGITUDE IS BLOCKED
!                     IF THERE IS A SWITCH BACK TO SEA OR VICE VERSA
!                     (SEE BELOW)
!                     LAND IS DEFINED AS ANYTHING ABOVE IBLOCKDPT(IX,K)
!                     ------------------------------------------
                        IF(IDEPTH(I,J).GE.0 ) LREALLAND=.TRUE. 
                        LLAND=.TRUE.
                        NIOBSLON=NIOBSLON+1 
                      ELSEIF (IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                        WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
!                     IF SEA ABOVE THE THRESHOLD THEN ONLY THAT
!                     GRID POINTS BLOCKS
!                     ------------------------------------------
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO

                    IF(LLAND) THEN
                      IF(LREALLAND) THEN
                        LNSW=.TRUE.  
                        IF(IDEPTH(I,ILATT).GE.IBLOCKDPT(IX,K)) THEN
                          L1ST=.TRUE.
                        ELSE
                          L1ST=.FALSE.
                        ENDIF
                        DO J=ILATT+1,ILATB
                          IF( ((IDEPTH(I,J).GE.IBLOCKDPT(IX,K))
     &                          .NEQV.L1ST)
     &                         .AND.
     &                         LNSW ) THEN
                            LNSW=.FALSE.
                          ENDIF
                          IF(((IDEPTH(I,J).GE.IBLOCKDPT(IX,K)).EQV.L1ST)
     &                        .AND.
     &                       .NOT. LNSW ) THEN
!                           LAND IS BLOCKING
                            NIOBSLON=IREINF*(ILATB-ILATT+1)
                            NBLOCKLAND=NBLOCKLAND+1
                            EXIT
                          ENDIF
                        ENDDO
                        IF(LNSW) NIOBSLON=ILATB-ILATT+1
                      ELSE
                        IF(PERCENTSHALLOW(IX,K).GT.PSHALLOWTRHS) THEN
                          NIOBSLON=ILATB-ILATT+1
                        ELSEIF(PERCENTLAND(IX,K).LT.PLANDTRHS) THEN
                          NIOBSLON=IREINF*(ILATB-ILATT+1)
                          NBLOCKLAND=NBLOCKLAND+1
                        ELSE
                          NIOBSLON=0
                        ENDIF
                      ENDIF
                    ENDIF

                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR-ILONL+1)+
     &                    (IREINF-1)*NBLOCKLAND*(ILATB-ILATT+1)

                  IOBSCOR(IX,K,IS)=
     &               NINT((1.-FLOAT(NOBSTRCT)/NTOTPTS)*1000)
                ELSE
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR+ILON-ILONL+1)
                  DO I=1,ILONR
                    NIOBSLON=0
                    DO J=ILATT,ILATB
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K)) THEN
                        NIOBSLON=ILATB-ILATT+1
                        EXIT
                      ELSEIF(IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                       WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO
                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO
                  DO I=ILONL,ILON
                    NIOBSLON=0
                    DO J=ILATT,ILATB
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K)) THEN
                        NIOBSLON=ILATB-ILATT+1
                        EXIT
                      ELSEIF(IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                       WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO
                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO
                  IOBSCOR(IX,K,IS)=
     &               NINT((1.-FLOAT(NOBSTRCT)/NTOTPTS)*1000)
                ENDIF

              ENDIF
            ENDDO

          ENDDO
!$OMP END PARALLEL DO
        ENDDO

!       then search east-west and average 
        DO IS=1,4
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1)
!$OMP+PRIVATE(K,XLATT,XLATB,ILATT,ILATB,IX)
!$OMP+PRIVATE(XLON,XLONL,XLONR,ILONL,ILONR,NOBSTRCT,NBLOCKLAND)
!$OMP+PRIVATE(J,NIOBSLAT,LLAND,LREALLAND,I,LNSW,L1ST)
!$OMP+PRIVATE(NTOTPTS,ITEMPEW,XX)
          DO K=1,NY
            IF(IS.EQ.1 .OR. IS.EQ.4) THEN
              XLATT=XLAT(K)+XDELLA
              XLATB=XLAT(K)
            ELSE
              XLATT=XLAT(K)
              XLATB=XLAT(K)-XDELLA
            ENDIF
            ILATT = NINT((90.- XLATT)*INVRES) + 1
            ILATT = MAX(1,MIN(ILATT,ILAT))
            ILATB = NINT((90.- XLATB)*INVRES) + 1
            ILATB = MAX(1,MIN(ILATB,ILAT))
            IF(ILATB.EQ.ILAT+1)ILATB=ILAT

            DO IX=1,NLONRGG(K)
              IF(WAMDEPTH(IX,K).LT.0) THEN
                XLON=AMOWEP + REAL(IX-1)*ZDELLO(K)
                IF(IS.EQ.1 .OR. IS.EQ.2) THEN
                  XLONL=XLON + RESOL
                  XLONR=XLON + ZDELLO(K) 
                ELSE
                  XLONL=XLON - ZDELLO(K) 
                  XLONR=XLON - RESOL
                ENDIF
                IF(XLONL.GT.180.) THEN
                  XLONL=XLONL-360.
                ENDIF
                IF(XLONR.GT.180.) THEN
                  XLONR=XLONR-360.
                ENDIF

                ILONL = NINT((XLONL + 180.)*INVRES) + 1
                ILONR = NINT((XLONR + 180.)*INVRES) + 1

                NOBSTRCT=0

                IF(ILONL.LE.ILONR) THEN
                  NBLOCKLAND=0
                  DO J=ILATT,ILATB
                    NIOBSLAT=0
                    LLAND=.FALSE.
                    LREALLAND=.FALSE.
                    DO I=ILONL,ILONR
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K) ) THEN
!                     IF LAND THEN THE FULL LONGITUDE IS BLOCKED
!                     IF THERE IS A SWITCH BACK TO SEA OR VICE VERSA
!                     (SEE BELOW)
!                     LAND IS DEFINED AS ANYTHING ABOVE IBLOCKDPT(IX,K)
!                     ------------------------------------------
                        LLAND=.TRUE.
                        IF(IDEPTH(I,J).GE.0 ) LREALLAND=.TRUE. 
                        NIOBSLAT=NIOBSLAT+1 
                      ELSEIF (IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                        WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
!                     IF SEA ABOVE THE THRESHOLD THEN ONLY THAT
!                     GRID POINTS BLOCKS
!                     ------------------------------------------
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO

                    IF(LLAND) THEN
                      IF(LREALLAND) THEN
                        LNSW=.TRUE.  
                        IF(IDEPTH(ILONL,J).GE.IBLOCKDPT(IX,K) ) THEN
                          L1ST=.TRUE.
                        ELSE
                          L1ST=.FALSE.
                        ENDIF
                        DO I=ILONL+1,ILONR
                          IF( ((IDEPTH(I,J).GE.IBLOCKDPT(IX,K))
     &                         .NEQV.L1ST)
     &                         .AND.
     &                         LNSW ) THEN
                            LNSW=.FALSE.
                          ENDIF
                          IF(((IDEPTH(I,J).GE.IBLOCKDPT(IX,K)).EQV.L1ST)
     &                        .AND.
     &                       .NOT. LNSW ) THEN
!                           LAND IS BLOCKING
                            NIOBSLAT=IREINF*(ILONR-ILONL+1)
                            NBLOCKLAND=NBLOCKLAND+1
                            EXIT
                          ENDIF
                        ENDDO
                        IF(LNSW) NIOBSLAT=ILONR-ILONL+1
                      ELSE
                        IF(PERCENTSHALLOW(IX,K).GT.PSHALLOWTRHS) THEN
                          NIOBSLAT=ILONR-ILONL+1
                        ELSEIF(PERCENTLAND(IX,K).LT.PLANDTRHS) THEN
                          NIOBSLAT=IREINF*(ILONR-ILONL+1)
                          NBLOCKLAND=NBLOCKLAND+1
                        ELSE
                          NIOBSLAT=0
                        ENDIF
                      ENDIF
                    ENDIF

                    NOBSTRCT=NOBSTRCT+NIOBSLAT
                  ENDDO

                  NTOTPTS=(ILATB-ILATT+1)*(ILONR-ILONL+1)+
     &                    (IREINF-1)*NBLOCKLAND*(ILONR-ILONL+1)
                  ITEMPEW=NINT((1.-FLOAT(NOBSTRCT)/NTOTPTS)*1000)
                ELSE
                  NTOTPTS=(ILATB-ILATT+1)*(ILONR+ILON-ILONL+1)
                  DO J=ILATT,ILATB
                    NIOBSLAT=0
                    DO I=1,ILONR
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K)) THEN
                        NIOBSLAT=ILONR+ILON-ILONL+1
                        GOTO 4444 
                      ELSEIF(IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                       WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO
                    DO I=ILONL,ILON
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(IX,K)) THEN
                        NIOBSLAT=ILONR+ILON-ILONL+1
                        EXIT
                      ELSEIF(IDEPTH(I,J).GE.ITHRSHOLD(IX,K) .AND.
     &                       WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(IX,K))THEN
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO
4444                CONTINUE
                    NOBSTRCT=NOBSTRCT+NIOBSLAT
                  ENDDO
                  ITEMPEW=NINT((1.-FLOAT(NOBSTRCT)/NTOTPTS)*1000)
                ENDIF
                XX=FLOAT(IOBSCOR(IX,K,IS)*ITEMPEW)
                XX=SQRT(XX)
                IOBSCOR(IX,K,IS)=MIN(NINT(XX),1000)

              ENDIF
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ENDDO

        ENDIF ! end if xdella too small, do not compute obstructions


!       OUTPUT OBSTRUCTIONS
!       FOR PLOTTING

        IF(LLPRINT) THEN

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
            IF(IS.EQ.1) THEN
              STEPLAT=-0.25*XDELLA
            ELSE
              STEPLAT=0.25*XDELLA
            ENDIF 
            DO K=1,NY
               DO IX=1,NLONRGG(K)
                 XLON=AMOWEP + REAL(IX-1)*ZDELLO(K)
                 IF(XLON.GT.180.) then
                   XLON=XLON-360.
                 ENDIF
                IF(ALATB.LE.XLAT(K) .AND. XLAT(K).LE.ALATT .AND.
     &             ALONL.LE.XLON .AND. XLON.LE.ALONR ) THEN
                   IF(WAMDEPTH(IX,K).LT.0 .AND.
     &                IOBSLAT(IX,K,IS).LT.1000 ) THEN
                      WRITE(IUNIT,'(2(1X,F8.3),1X,I4)') 
     &                XLON,XLAT(K)+STEPLAT,IOBSLAT(IX,K,IS)
                   ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO

          DO IS =1,2
            IUNIT=IUNIT+1
            IF(IS.EQ.1) THEN
              STEPLON=-0.25*XDELLO
            ELSE
              STEPLON=0.25*XDELLO
            ENDIF 
            DO K=1,NY
               DO IX=1,NLONRGG(K)
                 XLON=AMOWEP + REAL(IX-1)*ZDELLO(K)
                 IF(XLON.GT.180.) then
                   XLON=XLON-360.
                 ENDIF
                IF(ALATB.LE.XLAT(K) .AND. XLAT(K).LE.ALATT .AND.
     &             ALONL.LE.XLON .AND. XLON.LE.ALONR ) THEN
                   IF(WAMDEPTH(IX,K).LT.0 .AND.
     &                IOBSLON(IX,K,IS).LT.1000 ) THEN
                      WRITE(IUNIT,'(2(1X,F8.3),1X,I4)') 
     &                XLON+STEPLON,XLAT(K),IOBSLON(IX,K,IS)
                   ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO

          DO IS =1,2
            IUNIT=IUNIT+1
            IF(IS.EQ.1) THEN
              STEPLAT=-0.25*XDELLA
            ELSE
              STEPLAT=0.25*XDELLA
            ENDIF 
            DO K=1,NY
               DO IX=1,NLONRGG(K)
                 XLON=AMOWEP + REAL(IX-1)*ZDELLO(K)
                 IF(XLON.GT.180.) then
                   XLON=XLON-360.
                 ENDIF
                IF(ALATB.LE.XLAT(K) .AND. XLAT(K).LE.ALATT .AND.
     &             ALONL.LE.XLON .AND. XLON.LE.ALONR ) THEN
                   IF(WAMDEPTH(IX,K).LT.0 .AND.
     &                IOBSRLAT(IX,K,IS).LT.1000 ) THEN
                      WRITE(IUNIT,'(2(1X,F8.3),1X,I4)') 
     &                XLON,XLAT(K)+STEPLAT,IOBSRLAT(IX,K,IS)
                   ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO

          DO IS =1,2
            IUNIT=IUNIT+1
            IF(IS.EQ.1) THEN
              STEPLON=-0.25*XDELLO
            ELSE
              STEPLON=0.25*XDELLO
            ENDIF 
            DO K=1,NY
               DO IX=1,NLONRGG(K)
                 XLON=AMOWEP + REAL(IX-1)*ZDELLO(K)
                 IF(XLON.GT.180.) then
                   XLON=XLON-360.
                 ENDIF
                IF(ALATB.LE.XLAT(K) .AND. XLAT(K).LE.ALATT .AND.
     &             ALONL.LE.XLON .AND. XLON.LE.ALONR ) THEN
                   IF(WAMDEPTH(IX,K).LT.0 .AND.
     &                IOBSRLON(IX,K,IS).LT.1000 ) THEN
                      WRITE(IUNIT,'(2(1X,F8.3),1X,I4)') 
     &                XLON+STEPLON,XLAT(K),IOBSRLON(IX,K,IS)
                   ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO

          DO IS =1,4
            IUNIT=IUNIT+1
            IF(IS.EQ.1) THEN
              STEPLON=-0.25*XDELLO
            ELSE
              STEPLON=0.25*XDELLO
            ENDIF 
            DO K=1,NY
               DO IX=1,NLONRGG(K)
                 XLON=AMOWEP + REAL(IX-1)*ZDELLO(K)
                 IF(XLON.GT.180.) then
                   XLON=XLON-360.
                 ENDIF
                IF(ALATB.LE.XLAT(K) .AND. XLAT(K).LE.ALATT .AND.
     &             ALONL.LE.XLON .AND. XLON.LE.ALONR ) THEN
                   IF(WAMDEPTH(IX,K).LT.0 .AND.
     &                IOBSRLON(IX,K,IS).LT.1000 ) THEN
                      WRITE(IUNIT,'(2(1X,F8.3),1X,I4)') 
     &                XLON+STEPLON,XLAT(K),IOBSCOR(IX,K,IS)
                   ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO

        ENDIF


!       OUTPUT OBSTRUCTIONS
!       FOR GLOBAL FIELD (in the same file as mean bathymetry)
        DO IS =1,2
          DO K=1,NY
           WRITE(CX,'(I4.4)') NLONRGG(K) 
           FORMAT='('//CX//'I4)'
           WRITE(1,FORMAT) (IOBSLAT(IX,K,IS),IX=1,NLONRGG(K)) 
          ENDDO
        ENDDO
        DO IS =1,2
          DO K=1,NY
           WRITE(CX,'(I4.4)') NLONRGG(K) 
           FORMAT='('//CX//'I4)'
           WRITE(1,FORMAT) (IOBSLON(IX,K,IS),IX=1,NLONRGG(K)) 
          ENDDO
        ENDDO
        DO IS =1,2
          DO K=1,NY
           WRITE(CX,'(I4.4)') NLONRGG(K) 
           FORMAT='('//CX//'I4)'
           WRITE(1,FORMAT) (IOBSRLAT(IX,K,IS),IX=1,NLONRGG(K)) 
          ENDDO
        ENDDO
        DO IS =1,2
          DO K=1,NY
           WRITE(CX,'(I4.4)') NLONRGG(K) 
           FORMAT='('//CX//'I4)'
           WRITE(1,FORMAT) (IOBSRLON(IX,K,IS),IX=1,NLONRGG(K)) 
          ENDDO
        ENDDO
        DO IS =1,4
          DO K=1,NY
           WRITE(CX,'(I4.4)') NLONRGG(K) 
           FORMAT='('//CX//'I4)'
           WRITE(1,FORMAT) (IOBSCOR(IX,K,IS),IX=1,NLONRGG(K)) 
          ENDDO
        ENDDO

      ENDDO ! END LOOP ON FREQUENCIES

      END PROGRAM
