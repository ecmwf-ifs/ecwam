      PROGRAM CREATE_BATHY 

!     CREATES THE WAM BATHYMETRY AND THE REDUCTION FACTORS
!     DUE TO SUB-GRID BATHYMETRIC FEATURES USING
!     THE ETOPO2 DATA SET.

!     JEAN BIDLOT   FEBRUARY 2003

!     DESCRIPTION OF ETOPO2:      
!     THE ORIGINAL DATA FILE CAN BE FOUND ON TAPE:
!     ec:/wab/topo_wam/etopo2.tar
!     UNTAR THE FILE AND MAKE SURE THAT etopo2.dat IS ONLINE
!     ROWS OF DATA (LATITUDE LINES) ARE READ ONE BY
!     ONE STARTING AT 89 58' NORTH AND STEPPING SOUTHWARD
!     BY 2'. THERE ARE 5400 LATITUDE LINES.  THE 89 58'
!     SOUTH LINE IS NOT GIVEN, ELEVATION 2810 METERS CAN
!     BE ASSUMED.
!     THERE ARE 10800 LONGITUDE POINTS ALONG EACH LATITUDE
!     LINE. THE FIRST IS 180 wEST, SECOND 179 58' WEST, ETC.
!
!     DEPTH AND ELEVATIONS ARE IN METERS


!     INPUT FILES:
!     -----------
!     etopo2.dat

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

!     CORRECTION TO WAM POINTS.
!     THEY WERE OBTAINED FROM THE PREVIOUS GRID SETUP
!     correction_to_wam_grid_xxx
!     where xxxx=int(xdella*100) if grid_description is not present
!     else  xxxx=the gaussian number of latitude - 1 (NY)

!     OUTPUT FILES:

!     wamtopo_xxxx
!     where xxxx=int(xdella*100) if grid_description is not present
!     else  xxxx=the gaussian spectral truncation ISPECTRUNC

!     AND ALL THE FILES NECESSARY FOR PLOTTING WITH METVIEW
!     IF REQUESTED IN THE USER INPUT !
!     meandepth.dat
!     and
!     obstructions_*.dat
 
!**************************************************************************
 
      INTEGER, PARAMETER :: ILON=10800
      INTEGER, PARAMETER :: ILAT=5400
      INTEGER, PARAMETER :: NREF=300

      INTEGER, ALLOCATABLE :: NLONRGG(:)
      INTEGER, ALLOCATABLE :: ITHRSHOLD(:),IEXCLTHRSHOLD(:),IBLOCKDPT(:)
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IDEPTH,WAMDEPTH
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: IOBSLAT,IOBSLON
      INTEGER :: LEVEL(NREF),NDEPTH(NREF)

      REAL, ALLOCATABLE, DIMENSION(:) :: ZDELLO,COSPH
      REAL, ALLOCATABLE, DIMENSION(:) :: XLAT
      REAL, ALLOCATABLE, DIMENSION(:) :: FR
      REAL ::  ALON(ILON),ALAT(ILAT)
      REAL :: XINF(NREF),XSUP(NREF),YINF(NREF),YSUP(NREF)

      CHARACTER*2 :: CFR
      CHARACTER*4 :: CWAMRESOL
      CHARACTER*4 :: CX
      CHARACTER*10 :: FORMAT
      CHARACTER*23 :: FILENM
      CHARACTER*72 :: LOCATION(NREF)
      CHARACTER*144 :: CLINE,FILENAME

      LOGICAL :: LORIGINAL, LLPRINT
      LOGICAL :: LLAND, LREALLAND, L1ST, LNSW
      LOGICAL :: LLGRID

      PI=4.*ATAN(1.)
      RAD=PI/180.
      G=9.806

!     FIRST OR THIRD FREQUENCY IN WAM
      FR1= 4.177248E-02
      FRATIO=1.1
!!!   the obstructions will be produced for 30 frequencies
!!!   if 25 frequencies are used then the irrelevant
!!!   frequencies will be skipped (see topoar). Any other
!!!   number will not work.

      NFRE=30


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
        NX=INT((AMOEAP-AMOWEP)/XDELLO)+1
        NY=INT((AMONOP-AMOSOP)/XDELLA)+1

        ALLOCATE(NLONRGG(NY))
      ENDIF

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
            NLONRGG(K) = NINT(NX*COSPH(K))
            IF(MOD(NLONRGG(K),2).EQ.1) NLONRGG(K) = NLONRGG(K)+1
          ELSE
            NLONRGG(K) = NX
          ENDIF
        ENDIF

        PLONS=(AMOEAP-AMOWEP) + IPER*XDELLO
        IF(IPER.EQ.1) THEN
          ZDELLO(K)  = PLONS/REAL(NLONRGG(K))
        ELSE
          ZDELLO(K)  = PLONS/REAL(NLONRGG(K)-1)
        ENDIF
      ENDDO

      ALLOCATE(FR(NFRE))
      IF(NFRE.EQ.30) THEN
        FR(3)=FR1
        FR(2)=FR(3)/FRATIO
        FR(1)=FR(2)/FRATIO
        DO M=4,NFRE
          FR(M) = FRATIO*FR(M-1)
        ENDDO
      ELSEIF(NFRE.EQ.25) THEN
        FR(1)=FR1
        DO M=2,NFRE
          FR(M) = FRATIO*FR(M-1)
        ENDDO
      ELSE
        WRITE(*,*) 'UNKNOWN NUMBER OF FREQUENCIES !!!!'
        CALL ABORT
      ENDIF

!     COMPUTE THE THRESHOLD AT WHICH THE WAVES
!     ARE OBSTRUCTED BY THE BOTTOM.
!     EXCEPT IF WAM DEPTH OF THE SAME ORDER OF MAGITUDE IF IT'S SEA
!     ALSO COMPUTE THE DEPTH THAT IS CONSIDERED TO BE FULLY BLOCKING
!     AS IF IT WAS LAND.
      ALLOCATE(ITHRSHOLD(NFRE))
      ALLOCATE(IBLOCKDPT(NFRE))
      ALLOCATE(IEXCLTHRSHOLD(NFRE))
      DO M=1,NFRE
        XX=G/(2*PI**2*FR(M)**2)
        ITHRSHOLD(M)=NINT(-XX)
        IEXCLTHRSHOLD(M)=MAX(10*ITHRSHOLD(M),-998)
        IBLOCKDPT(M)=INT(-0.05*XX)
      ENDDO


!     DATASET:
!     --------
!     COMPUTE THE LATITUDE AND LONGITUDE POINTS
      DO I = 1,ILON
        ALON(I) = -180.0 + (FLOAT(I)-1.0)*2.0/60.0
      ENDDO
      DO J = 1,ILAT
        ALAT(J) = 90.0 - (FLOAT(J) -1.0)*2.0/60.0
      ENDDO
 
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
      
      ILONL = NINT((ALONL + 180.)*30) + 1
      ILONR = NINT((ALONR + 180.)*30) + 1
      ILATB = NINT((90.- ALATB)*30) + 1
      ILATT = NINT((90.- ALATT)*30) + 1
      IF(ILONR.EQ.ILON+1)ILONR=ILON
      IF(ILATB.EQ.ILAT+1)ILATB=ILAT

      ALLOCATE(IDEPTH(ILON,ILAT))
      ALLOCATE(WAMDEPTH(NX,NY))

      OPEN(20,FILE='etopo2.dat')

!     READ INPUT DATA
      DO J=1,ILAT
         READ(20,'(15I6)')(IDEPTH(I,J),I=1,ILON)
      ENDDO

!     OUTPUT ORIGINAL DATA SET ON SUB AREA

      IF(LORIGINAL) THEN
        OPEN(11,file='depth.dat')
        WRITE(11,'(a4)') '#GEO'
        WRITE(11,'(a11)') '#FORMAT LLV'
        WRITE(11,'(a5)') '#DATA'
      
        DO J=ILATT,ILATB
          DO I=ILONL,ILONR
            IF(IDEPTH(I,J).GE. -300 .and. 
     &         IDEPTH(I,J).LE. 0 ) THEN
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
                   IDEPTH(I,J)=NDEPTH(IR)
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

      NJM=INT(0.5*XDELLA*(60/2))
      NJP=NJM

      DO K=1,NY
!        WE ASSUME THAT WAMGRID IS ALWAYS WITHIN ETO2
!        DETERMINE CLOSEST ETO2 J INDEX TO WAM POINT
         DO J=ILAT-1,1,-1
           IF(ALAT(J+1).LE.XLAT(K).AND.
     &        XLAT(K).LT.ALAT(J) ) EXIT
         ENDDO
         IF(J.EQ.ILAT-1 .OR. J.LE.1) THEN
            WRITE(*,*) 'PROBLEM WITH J !!!'
            CALL ABORT
         ENDIF

         DO IX=1,NLONRGG(K)

!          ETO2 STARTS AT -180
           XLON=AMOWEP + REAL(IX-1)*ZDELLO(K)
           IF(XLON.GT.180.) THEN
             XLON=XLON-360.
           ENDIF

!          DETERMINE CLOSEST ETO2 I INDEX TO WAM POINT
           DO I=1,ILON-1
             IF(ALON(I).LE.XLON.AND.
     &          XLON.LT.ALON(I+1) ) EXIT
           ENDDO

           NIM=INT(0.5*ZDELLO(K)*(60/2))
           NIP=NIM

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
                 IF(II.LT.1) IK=ILON-II
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
           DO JJ=J-1,J+1
             IF(JJ.GE.1 .AND. JJ.LE.ILAT) THEN
               DO II=I-1,I+1
                 IK=II
                 IF(II.LT.1) IK=ILON-II
                 IF(II.GT.ILON) IK=II-ILON
                 IF(IDEPTH(IK,JJ).GT.0) THEN
                   NLANDCENTRE=NLANDCENTRE+1
                 ENDIF
               ENDDO
             ENDIF
           ENDDO

!          IF 50% OR MORE LAND, THEN AVERAGE OVER LAND POINTS
!          OR THE CENTER OF THE GRID BOX IS LAND.
!          ELSE AVERAGE OVER SEA POINTS
           IF(FLOAT(NLAND)/FLOAT(NLAND+NSEA).GT.0.50 .OR.
     &        NLANDCENTRE.GE.9 ) THEN
             WAMDEPTH(IX,K)=XLAND/NLAND
           ELSE
!            IF THERE IS A PERCENTAGE OF SHALLOWER POINTS THEN
!            THE AVERAGE IS TAKEN OVER THOSE POINTS ALONE.
             IF(FLOAT(NSEASH)/FLOAT(NSEA).GE.0.3) THEN
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

      WRITE(*,*) 'CORRECTION TO WAM GRID'

      IF(LLGRID) THEN
        WRITE(CWAMRESOL,'(I4.4)') NY-1
      ELSE
        WRITE(CWAMRESOL,'(I4.4)') INT(100*XDELLA)
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
        WRITE(CWAMRESOL,'(I4.4)') ISPECTRUNC
      ELSE
        WRITE(CWAMRESOL,'(I4.4)') INT(100*XDELLA)
      ENDIF
      FILENAME='wamtopo_'//CWAMRESOL

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



!     CREATE OBSTRUCTIONS:
!     --------------------

      ALLOCATE(IOBSLAT(NX,NY,2))
      ALLOCATE(IOBSLON(NX,NY,2))

      WRITE(1,'(I4)') NFRE

!     IREINF IS USED TO REINFORCE LAND OBSTRUCTIONS FOR
!     SMALL GRID SPACING.
      IF(XDELLA.LE.0.5) THEN
        IREINF=2
      ELSE
        IREINF=1
      ENDIF

      DO M=1,NFRE

        IOBSLAT=1000
        IOBSLON=1000

!       NORTH-SOUTH OBSTRUCTIONS
!       -----------------------
!       IS=1 is for the south-north advection
!       IS=2 is for the north-south advection
        DO IS=1,2
          DO K=1,NY
            IF(IS.EQ.1) THEN
              KT=K
              KB=K-1
              STEPT=-1./30.
              STEPB=0.
            ELSE
              KT=K+1
              KB=K
              STEPT=0.
              STEPB=1./30.
            ENDIF
            XLATT=XLAT(KT)+STEPT
            XLATB=XLAT(KB)+STEPB
            ILATT = NINT((90.- XLATT)*30) + 1
            ILATB = NINT((90.- XLATB)*30) + 1
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

                ILONL = NINT((XLONL + 180.)*30) + 1
                ILONR = NINT((XLONR + 180.)*30) + 1

                NOBSTRCT=0

                IF(ILONL.LE.ILONR) THEN
                  NBLOCKLAND=0
                  DO I=ILONL,ILONR
                    NIOBSLON=0
                    LLAND=.FALSE.
                    LREALLAND=.FALSE.
                    DO J=ILATT,ILATB
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(M) ) THEN
!                     IF LAND THEN THE FULL LONGITUDE IS BLOCKED
!                     IF THERE IS A SWITCH BACK TO SEA OR VICE VERSA
!                     (SEE BELOW)
!                     LAND IS DEFINED AS ANYTHING ABOVE IBLOCKDPT(M)
!                     ------------------------------------------
                        IF(IDEPTH(I,J).GE.0 ) LREALLAND=.TRUE. 
                        LLAND=.TRUE.
                        NIOBSLON=NIOBSLON+1 
                      ELSEIF (IDEPTH(I,J).GE.ITHRSHOLD(M) .AND.
     &                        WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(M) ) THEN
!                     IF SEA ABOVE THE THRESHOLD THEN ONLY THAT
!                     GRID POINTS BLOCKS
!                     ------------------------------------------
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO

                    IF(LLAND) THEN
                      IF(LREALLAND) THEN
                        LNSW=.TRUE.  
                        IF(IDEPTH(I,ILATT).GE.IBLOCKDPT(M)) THEN
                          L1ST=.TRUE.
                        ELSE
                          L1ST=.FALSE.
                        ENDIF
                        DO J=ILATT+1,ILATB
                          IF( ((IDEPTH(I,J).GE.IBLOCKDPT(M)).NEQV.L1ST)
     &                         .AND.
     &                         LNSW ) THEN
                            LNSW=.FALSE.
                          ENDIF
                          IF(((IDEPTH(I,J).GE.IBLOCKDPT(M)).EQV.L1ST)
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
                        NIOBSLON=IREINF*(ILATB-ILATT+1)
                        NBLOCKLAND=NBLOCKLAND+1
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
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(M)) THEN
                        NIOBSLON=ILATB-ILATT+1
                        EXIT
                      ELSEIF(IDEPTH(I,J).GE.ITHRSHOLD(M) .AND.
     &                       WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(M) ) THEN
                        NIOBSLON=NIOBSLON+1 
                      ENDIF
                    ENDDO
                    NOBSTRCT=NOBSTRCT+NIOBSLON
                  ENDDO
                  DO I=ILONL,ILON
                    NIOBSLON=0
                    DO J=ILATT,ILATB
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(M)) THEN
                        NIOBSLON=ILATB-ILATT+1
                        EXIT
                      ELSEIF(IDEPTH(I,J).GE.ITHRSHOLD(M) .AND.
     &                       WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(M) ) THEN
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
        ENDDO

!       EAST-WEST OBSTRUCTIONS
!       -----------------------
!       IS=1 is for the west-east advection
!       IS=2 is for the east-west advection
        DO IS=1,2
          DO K=1,NY
            XLATT=XLAT(K)+0.5*XDELLA
            XLATB=XLAT(K)-0.5*XDELLA
            ILATT = NINT((90.- XLATT)*30) + 1
            ILATB = NINT((90.- XLATB)*30) + 1
            IF(ILATB.EQ.ILAT+1)ILATB=ILAT

            DO IX=1,NLONRGG(K)
              IF(WAMDEPTH(IX,K).LT.0) THEN
                IF(IS.EQ.1) THEN
                  XLONL=AMOWEP + (REAL(IX-2))*ZDELLO(K)
                  XLONR=AMOWEP + (REAL(IX-1))*ZDELLO(K) -1./30.
                ELSE
                  XLONL=AMOWEP + (REAL(IX-1))*ZDELLO(K) +1./30.
                  XLONR=AMOWEP + (REAL(IX))*ZDELLO(K)
                ENDIF
                IF(XLONL.GT.180.) THEN
                  XLONL=XLONL-360.
                ENDIF
                IF(XLONR.GT.180.) THEN
                  XLONR=XLONR-360.
                ENDIF

                ILONL = NINT((XLONL + 180.)*30) + 1
                ILONR = NINT((XLONR + 180.)*30) + 1

                NOBSTRCT=0

                IF(ILONL.LE.ILONR) THEN
                  NBLOCKLAND=0
                  DO J=ILATT,ILATB
                    NIOBSLAT=0
                    LLAND=.FALSE.
                    LREALLAND=.FALSE.
                    DO I=ILONL,ILONR
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(M) ) THEN
!                     IF LAND THEN THE FULL LONGITUDE IS BLOCKED
!                     IF THERE IS A SWITCH BACK TO SEA OR VICE VERSA
!                     (SEE BELOW)
!                     LAND IS DEFINED AS ANYTHING ABOVE IBLOCKDPT(M)
!                     ------------------------------------------
                        LLAND=.TRUE.
                        IF(IDEPTH(I,J).GE.0 ) LREALLAND=.TRUE. 
                        NIOBSLAT=NIOBSLAT+1 
                      ELSEIF (IDEPTH(I,J).GE.ITHRSHOLD(M) .AND.
     &                        WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(M) ) THEN
!                     IF SEA ABOVE THE THRESHOLD THEN ONLY THAT
!                     GRID POINTS BLOCKS
!                     ------------------------------------------
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO

                    IF(LLAND) THEN
                      IF(LREALLAND) THEN
                        LNSW=.TRUE.  
                        IF(IDEPTH(ILONL,J).GE.IBLOCKDPT(M) ) THEN
                          L1ST=.TRUE.
                        ELSE
                          L1ST=.FALSE.
                        ENDIF
                        DO I=ILONL+1,ILONR
                          IF( ((IDEPTH(I,J).GE.IBLOCKDPT(M)).NEQV.L1ST)
     &                         .AND.
     &                         LNSW ) THEN
                            LNSW=.FALSE.
                          ENDIF
                          IF(((IDEPTH(I,J).GE.IBLOCKDPT(M)).EQV.L1ST)
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
                        NIOBSLAT=IREINF*(ILONR-ILONL+1)
                        NBLOCKLAND=NBLOCKLAND+1
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
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(M)) THEN
                        NIOBSLAT=ILONR+ILON-ILONL+1
                        GOTO 1111 
                      ELSEIF(IDEPTH(I,J).GE.ITHRSHOLD(M) .AND.
     &                       WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(M) ) THEN
                        NIOBSLAT=NIOBSLAT+1 
                      ENDIF
                    ENDDO
                    DO I=ILONL,ILON
                      IF(IDEPTH(I,J).GE.IBLOCKDPT(M)) THEN
                        NIOBSLAT=ILONR+ILON-ILONL+1
                        EXIT
                      ELSEIF(IDEPTH(I,J).GE.ITHRSHOLD(M) .AND.
     &                       WAMDEPTH(IX,K).LT.IEXCLTHRSHOLD(M) ) THEN
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
        ENDDO

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

      ENDDO ! END LOOP ON FREQUENCIES

      END PROGRAM
