      SUBROUTINE MPDECOMP(NPR,MAXLEN,LLIRANK)

!****  *MPDECOMP* - EVEN DECOMPOSITION OF THE GRID DOMAIN
!****               AMONG PROCESSES AFTER INPUT OF PREPROC GRID
!                   AND UBUF FILES.

!     J. BIDLOT    ECMWF   MARCH 1996  MESSAGE PASSING
!     J. BIDLOT    ECMWF   JANUARY 1998 introduce NPR 
!     J. BIDLOT    ECMWF   OCTOBER 1998 COMPLETE READING OF
!                                       IU07 AND IU08 
!     J. BIDLOT    ECMWF   FEBRUARY 1999 TAUT --> SQRT(TAUT)
!     J. BIDLOT    ECMWF   OCTOBER 2000 NOW READING SQRT(TAUT)

!     J. BIDLOT    ECMWF   FEBRUARY 2002 NEW 2-D DECOMPOSITION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     J. BIDLOT    ECMWF       11/2003
!                  IF YOUR ARE RUNNING AT ECMWF:
!                  BE AWARE THAT IF YOU CHANGE ANYTHING TO THE
!                  STRUCTURE OF THE UBUF FILE YOU WILL HAVE TO
!                  MAKE SURE THAT IT IS CREATED FOR YOUR RUN,
!                  OTHERWISE IT MIGHT PICK UP THE DEFAULT ONE
!                  THAT IS ALREADY ON DISK.
!                  YOU ALSO HAVE TO CHANGE OUTUBUF.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     PURPOSE.
!     --------
!     IT WILL READ IU07 AND IU08 AND IF MESSAGE PASSING IT WILL
!     DETERMINE AN EVEN DECOMPOSITION OF THE GRID ARRAYS FOR USE ON
!     A DISTRIBUTED MEMORY COMPUTER USING A MESSAGE PASSING PROTOCOL
!     FOR THE EXCHANGE OF INFORMATION ACROSS THE DIFFERENT PE's
!     THE LENGTH OF THE MESSAGES BETWEEN DIFFERENT PE's IS ALSO
!     COMPUTED.

!*    INTERFACE.
!     ----------
!     CALL *MPDECOMP*(NPR,MAXLEN,LLIRANK)

!      *NPR*       NUMBER OF SUBDOMAINS (USUALLY THE NUMBER OF PE'S )
!      *MAXLEN*    MAXIMUM NUMBER OF POINTS IN ANY SUB DOMAIN
!      *LLIRANK*   IF TRUE THEN SOME OF THE ARRAYS THAT ARE ONLY USED
!                  IF LMESSPASS IS TRUE ARE GENERATED EVEN IF IT IS NOT
!                  A MESSAGE PASSING RUN (USED WITH WRITE_MPDECOMP). 

!     METHOD.
!     -------

!    THE DOMAIN WILL BE SUBDIVIDED INTO NPR SUBAREAS IN SUCH A WAY
!    THAT EACH ONE WILL CONTAIN THE SAME NUMBER OF SEA POINTS (+-1).
!    THE MODEL IS FIRST DIVIDED INTO NYDECOMP LATITUDINAL BANDS
!    THIS IS REFFERED TO AS THE FIRST 1D DECOMPOSITION AND IS
!    IDENTICAL TO THE DECOMPOSITION WHICH WAS USED BEFORE. 
!    THE SECOND 2D DECOMPOSITION USED A SIMILAR PROCEDURE TO SPLIT
!    EACH BANDS INTO SUBAREAS OF EQUAL NUMBER OF POINTS.
!    THERE WILL BE NXDECOMP SUBAREAS IN THE FIRST NYCUT LATITUDINAL
!    BANDS (starting from the southern boundary)
!    AND (NXDECOMP-1) SUBAREAS IN THE REMAINING (NYDECOMP-NYCUT) BANDS
!    IN SUCH A WAY THAT 
!    NPR=NXDECOMP*NYCUT+(NYDECOMP-NYCUT)*(NXDECOMP-1) 
 
!    NXDECOMP and NYDECOMP ARE DETERMINED IN SUCH A WAY THAT THE SUB
!    AREAS ARE AS SQUARE AS POSSIBLE BY ASSUMING THAT THE GLOBAL
!    EXTEND OF THE DOMAIN IS TWICE AS LONG IN THE LONGITUDINAL DIRECTION
!    THAN IN THE LATITUDINAL DIRECTION (as is the case for the globe).

!    BASED ON THE CURRENT ADVECTION SCHEME (SEE PROPAGS),
!    THE DECOMPOSITION THEN YIELDS FOR EACH PE
!    NGBTOPE  WHICH GIVES THE TOTAL NUMBER OF NEIGHBOURING PE'S TO
!    WHICH INFORMATION FROM THE LOCAL PE IS POTENTIALLY NEEDED.
!    NTOPELST(INGB) INGB=1,NGBTOPE  THE LIST OF NEIGHBOURING PE'S ,
!    NTOPE(IPROC) IPROC=1,NPR  THE NUMBER OF POINTS FOR WHICH INFORMATION
!    HAS TO BE SENT TO PE IPROC,
!    IJTOPE(IH,IPROC) IH=1,NTOPE(IPROC), IPROC=1,NPR THE IJ INDEX OF THOSE
!    POINTS FOR WHICH INFORMATION HAS TO BE SENT TO PE IPROC,
!
!    NGBFROMPE WHICH GIVES THE TOTAL NUMBER OF NEIGHBOURING PE'S FROM
!    WHICH INFORMATION IS POTENTIALLY NEEDED,
!    NFROMPELST(KNGB) KNGB=1,NGBFROMPE  THE LIST OF NEIGHBOURING PE'S ,
!    NFROMPE(IPROC) IPROC=1,NPR  THE NUMBER OF POINTS THE LOCAL PE HAS
!    TO RECEIVE INFORMATION FROM PE IPROC,
!    NIJSTART(IPROC) IPROC=1,NPR THE IJ INDEX OF THE FIRST HALO POINT
!    OBTAINED FROM PE IPROC IN THE BUFFERS THAT PADS BOTH SIDE OF THE
!    1-D SEA POINT ARRAY OF THE LOCAL PE. THE OTHER POINTS ARE STORED IN
!    SUCCESSIVE ORDER FROM THAT IJSTART. THE ARRAYS KLON, KLAT, KALT2,
!    DEPTH, DELLAM1, COSPHM1 ARE LOCALLY REDEFINED
!    TO CORRECTLY POINT TO THOSE HALO POINTS IN THE PADDING BUFFERS.
!    BY CONVENTION, CONTRIBUTIONS FROM PE WITH PE NUMBER LESS THAN
!    THE LOCAL ONE (IRANK) ARE PUT IN THE BOTTOM BUFFER
!    (I.E. THEIR IJ'S ARE LESS THAN NSTART) AND THOSE CONTRIBUTIONS WITH
!    PE NUMBER GREATER THAN THE LOCAL ONE (IRANK) ARE PUT IN TOP BUFFER 
!    (I.E. THEIR IJ'S ARE GREATER THAN NEND).

!     EXTERNALS.
!     ----------
!          NONE

!     REFERENCES.
!     -----------
!          NONE

! -------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWFRED  , ONLY : FR       ,DFIM      ,DFIMOFR  ,FR5     ,
     &             FRM5     ,COFRM4   ,COEF4    ,DFFR     ,DFFR2
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL      ,IJLT     ,
     &            DELLAM   ,DELLAM1  ,COSPH    ,COSPHM1,
     &            DELPHI
      USE YOWICE   , ONLY : FLMIN    ,FLMINFR
      USE YOWJONS  , ONLY : AJONS    ,BJONS    ,DJONS    ,EJONS
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,KXLTMIN  ,KXLTMAX  ,
     &            IPER     ,IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,LLBOUND  ,
     &            KMNOP    ,KMSOP
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NINF     ,NSUP     ,KTAG 
      USE YOWPARAM , ONLY : NFRE     ,NBLO     ,NIBLO    ,NGY      ,
     &            LL1D
      USE YOWPCONS , ONLY : G
      USE YOWSHAL  , ONLY : DEPTH
      USE YOWSTAT  , ONLY : IREFRA 
      USE YOWTABL  , ONLY : JUMAX    ,UMAX     ,DELU
      USE YOWTEST  , ONLY : IU06
      USE YOWSPEC  , ONLY : NSTART   ,NEND     ,KLENTOP  ,KLENBOT  ,
     &            NFROMPE  ,NFROMPEMAX,NTOPE   ,NTOPEMAX ,NIJSTART ,
     &            IJTOPE   ,NGBTOPE  ,NTOPELST ,NGBFROMPE,NFROMPELST,
     &            IJ2NEWIJ
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,KLAT2    ,WADV     ,
     &            OBSLAT   ,OBSLON
      USE YOWUNIT  , ONLY : IU07     ,IU08
      USE MPL_MODULE


!----------------------------------------------------------------------


      INTEGER,ALLOCATABLE :: ICOMBUF(:), KDUM(:,:) 
      INTEGER,ALLOCATABLE, DIMENSION(:) :: INDEX,IJNDEX,NTOTSUB
      INTEGER,ALLOCATABLE, DIMENSION(:) :: NEWIJ2IJ
      INTEGER,ALLOCATABLE, DIMENSION(:) :: NSTART1D,NEND1D
      INTEGER,ALLOCATABLE, DIMENSION(:) :: KSTART1,KEND1, NLON, ILON
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IXLON
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITEMP,KTEMP,NLENHALO,IJHALO
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IJFROMPE,IPROCFROM
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: KOBSLON, KOBSLAT

      REAL,ALLOCATABLE :: RCOMBUF(:), RDUM(:,:)

      CHARACTER*80 LOGFILENAME
      LOGICAL :: LLIRANK

!     STANDARD OUTPUT UNIT
!     --------------------
      IF(LMESSPASS.AND.LWCOU) THEN
        IU06=20
      ELSEIF (LMESSPASS.AND..NOT.LWCOU) THEN
        IU06=66
        LOGFILENAME='logfile.%p'
        CALL EXPAND_STRING(IRANK,NPR,0,0,LOGFILENAME,1)
        OPEN(IU06,FILE=LOGFILENAME,STATUS='UNKNOWN')
        WRITE(*,*) ' STDOUT FILE ', LOGFILENAME, ' OPENED ON PE ', IRANK
        WRITE(IU06,*) ' STDOUT OF PE ', IRANK
        WRITE(IU06,*) ' ============ '
        CALL FLUSH (IU06)
      ELSE
        IU06=6
      ENDIF

!     GRID AND UBUF FILES UNITS
!     -------------------------

      IF(LWCOU) THEN
        IU07=29
        IU08=31
      ELSE
        IU07=7
        IU08=8
      ENDIF

!     0. READ GRID INPUT FROM PREPROC 
!        ----------------------------

      IREAD=1
      KTAG=1

      IF((LMESSPASS.AND.IRANK.EQ.IREAD).OR.(.NOT.LMESSPASS)) THEN
        CALL READPRE (IU07, IREFRA)
        CLOSE (UNIT=IU07, STATUS='KEEP')
      ENDIF

!     SEND INFORMATION FROM READPRE TO ALL PE's
      IF(LMESSPASS) CALL MPBCASTGRID(IU06,IREAD,KTAG,IREFRA)

      WRITE(IU06,*) ' PREPROC GRID INFORMATION READ IN '

!     DEFINE COEFFICIENT FOR MEAN PERIODS CALCULATION
      ALLOCATE(DFIMOFR(NFRE))
      ALLOCATE(DFFR(NFRE))
      ALLOCATE(DFFR2(NFRE))
      DO M=1,NFRE
        DFIMOFR(M) = DFIM(M)/FR(M)
        DFFR(M) = DFIM(M)*FR(M)
        DFFR2(M) = DFIM(M)*FR(M)**2
      ENDDO

!     DEFINE A FEW CONSTANTS FOR USE IN IMPLSCH

      ALLOCATE(FR5(NFRE))
      ALLOCATE(FRM5(NFRE))
      ALLOCATE(COFRM4(NFRE))

      DO M=1,NFRE
        FR5(M) = FR(M)**5.
        FRM5(M) = (1./FR(M))**5.
        COFRM4(M) = COEF4*G*FR(M)**(-4.)
      ENDDO

!     FIND KNMOP AND KMSOP
      KMNOP=1
      KMSOP=NGY
      DO IG=1,IGL
        DO IJ = IJS(IG),IJL(IG)
          KMNOP=MAX(KXLT(IJ,IG),KMNOP)
          KMSOP=MIN(KXLT(IJ,IG),KMSOP)
        ENDDO
      ENDDO


!     MINIMUM ENERGY IN SPECTRAL BINS 

      ALLOCATE(FLMINFR(JUMAX,NFRE))

      GAMMA=3.000000
      SA=7.000000E-02
      SB=9.000000E-02
      FETCH=MIN(0.5*DELPHI,25000.)
      GX = G * FETCH

!     THE ACTUAL VALUE OF ALPHJO IS NOT IMPORTANT. IT IS RESCALED
      ALPHJO = 1.

      DO J = 1,JUMAX
        U10=REAL(J)*DELU
!       find peak period  from fetch law
!       the same formulation as in subroutine peak is used.
        GXU = GX/(U10*U10)
        UG = G/U10
        FPK = AJONS * GXU ** DJONS
        FPK = MAX(0.13, FPK)
        FPK = FPK*UG
        FPK = MIN(FR(NFRE), FPK)
        CALL JONSWAP(ALPHJO, GAMMA, SA, SB, FPK, 1, 1, FLMINFR(J,:))
      ENDDO
!     RESCALE SPECTRUM SUCH THE MAX IS ALWAYS EQUAL to 4*FLMIN
!     (factor 4 because of the cos square spreading)
      DO J = 1,JUMAX
        FLMAX=0.
        DO M = 1,NFRE
          FLMAX=MAX(FLMAX,FLMINFR(J,M))
        ENDDO
        RSCALE=4*FLMIN/FLMAX
        DO M = 1,NFRE
          FLMINFR(J,M)=RSCALE*FLMINFR(J,M)
        ENDDO
      ENDDO


!*    1. INPUT NEIGHBOURING GRID POINT INDICES (UBUF) 
!        --------------------------------------------

      ALLOCATE(KLAT(NIBLO,2))  ! THE SIZE OF KLAT IS READJUSTED SEE BELOW
      ALLOCATE(KLON(NIBLO,2))  ! THE SIZE OF KLON IS READJUSTED SEE BELOW
      ALLOCATE(KLAT2(NIBLO,2)) ! THE SIZE OF KLAT2 IS READJUSTED SEE BELOW

      IF(LMESSPASS) THEN

!       READ IOU8 ON PE IREAD
        IF(IRANK.EQ.IREAD) THEN
          READ (IU08) KLAT, KLON
          READ (IU08) KLAT2
        ENDIF

!       SEND KLAT AND KLON TO OTHER PE'S
        IF(NPR.GT.1) THEN
          ITAG=KTAG
          MPLENGTH=6*NIBLO
          ALLOCATE(ICOMBUF(MPLENGTH))

          IF(IRANK.EQ.IREAD) THEN
            KCOUNT=0
            DO IC=1,2
              DO IJ=1,NIBLO
                KCOUNT=KCOUNT+1
                ICOMBUF(KCOUNT)=KLAT(IJ,IC)
                KCOUNT=KCOUNT+1
                ICOMBUF(KCOUNT)=KLAT2(IJ,IC)
              ENDDO
            ENDDO
            DO IC=1,2
              DO IJ=1,NIBLO
                KCOUNT=KCOUNT+1
                ICOMBUF(KCOUNT)=KLON(IJ,IC)
              ENDDO
            ENDDO
          ENDIF

          CALL MPL_BROADCAST(ICOMBUF(1:MPLENGTH),KROOT=IREAD,
     &     KTAG=ITAG,CDSTRING='MPDECOMP:')

          IF(IRANK.NE.IREAD) THEN
            KCOUNT=0
            DO IC=1,2
              DO IJ=1,NIBLO
                KCOUNT=KCOUNT+1
                KLAT(IJ,IC)=ICOMBUF(KCOUNT)
                KCOUNT=KCOUNT+1
                KLAT2(IJ,IC)=ICOMBUF(KCOUNT)
              ENDDO
            ENDDO
            DO IC=1,2
              DO IJ=1,NIBLO
                KCOUNT=KCOUNT+1
                KLON(IJ,IC)=ICOMBUF(KCOUNT)
              ENDDO
            ENDDO
          ENDIF

          DEALLOCATE(ICOMBUF)

        ENDIF

      ELSEIF(.NOT.LMESSPASS) THEN
        READ (IU08) KLAT, KLON
        READ (IU08) KLAT2
      ENDIF


!*    2. FIND THE NUMBER OF POINTS PER PE, THE START AND END INDEX
!        --------------------------------------------------------------

      ALLOCATE (NSTART(NPR),NEND(NPR),KLENBOT(NPR),KLENTOP(NPR))
      ALLOCATE (NFROMPE(NPR),NTOPE(NPR),NIJSTART(NPR))

!     DETERMINE THE SRUCTURE OF THE MODEL DECOMPOSITION

      IF(LL1D) THEN
!     1D DECOMPOSITION ONLY (old)
        NXDECOMP=1
        NYDECOMP=NPR
        NYCUT=NYDECOMP
	  ELSE
!     2D DECOMPOSITION (new)
        IF(NPR.EQ.1) THEN
          NXDECOMP=1
          NYDECOMP=1
          NYCUT=1
        ELSEIF(NPR.EQ.2) THEN
          NXDECOMP=2
          NYDECOMP=1
          NYCUT=1
        ELSE
          IPROC=0
          ICOUNT=0
!         find whether NPR can be expressed as 2*i**2 i=1,2,3,...
!         because in that case
!         NPR=NXDECOMP*NYCUT+(NYDECOMP-NYCUT)*(NXDECOMP-1)  (1) 
!         is satisfied with NXDECOMP=2*NYDECOMP and NYCUT=NYDECOMP
!         which is a perfect subdivision into identical squares 
!         of a rectangle of dimension 2L by L.
          DO WHILE (IPROC.LT.NPR)
            ICOUNT=ICOUNT+1
            IPROC=2*ICOUNT**2
          ENDDO
          IF(IPROC.EQ.NPR) THEN
            NYDECOMP=INT(SQRT(FLOAT(NPR)/2))
            NXDECOMP=2*NYDECOMP
            NYCUT=NYDECOMP
          ELSE
!         if the even decomposition into squares is not possible
!         start with the following approximation for NYDECOMP
!         found by setting NXDECOMP=2*NYDECOMP into (1) and
!         play around NXDECOMP=2*NYDECOMP,NYDECOMP,-1 and
!                     NYCUT=NYDECOMP,1,-1 until a solution to (1)
!         is reached. 
            IPROC=0
            NYDECOMP=INT(SQRT(FLOAT(NPR)/2))+1
            DO NXDECOMP=2*NYDECOMP,NYDECOMP,-1
              DO NYCUT=NYDECOMP,1,-1
                IPROC=NYDECOMP*(NXDECOMP-1)+NYCUT
                IF(IPROC.EQ.npr) EXIT
              ENDDO
              IF(IPROC.EQ.NPR) EXIT
            ENDDO
            IF(IPROC.NE.NPR) THEN
              WRITE(IU06,*) 'MPDECOMP :  decomposition problem !!!!'
              CALL ABORT1
            ENDIF
          ENDIF
        ENDIF
      ENDIF


!     FIRST 1-D DECOMPOSITION IN LATITUDINAL BANDS

      ALLOCATE(NSTART1D(NYDECOMP))
      ALLOCATE(NEND1D(NYDECOMP))

      IF(NYCUT.EQ.NYDECOMP) THEN
!       if the number of subareas per latitunal bands is the same
!       in all bands then the number of sea points in each band will
!       be determined to be as even as possible
        NMEAN=IJL(NBLO)/NYDECOMP
        NREST=IJL(NBLO)-NMEAN*NYDECOMP
     
        NSTART1D(1)=1
        IF(NREST.GT.0) THEN
          NPTS=NMEAN+1
          NREST=NREST-1 
        ELSE
          NPTS=NMEAN
        ENDIF
        NEND1D(1)=NSTART1D(1)+NPTS-1

        DO IP=2,NYDECOMP
          NSTART1D(IP)=NSTART1D(IP-1)+NPTS
          IF(NREST.GT.0) THEN
            NPTS=NMEAN+1
            NREST=NREST-1 
          ELSE
            NPTS=NMEAN
          ENDIF
          NEND1D(IP)=NSTART1D(IP)+NPTS-1
        ENDDO
      ELSE
!       if the number of subareas per latitunal bands is not the same
!       then the number of sea points per latitudinal bands will be
!       determined in such a way that number of points for the bands with
!       the least subareas (the top (nydecomp-nycut) bands) will
!       roughly scale like (nxdecomp-1)/nxdecomp the number of points
!       in the remaining bottom nycut bands

        NMEAN=NXDECOMP*IJL(NBLO)/((NXDECOMP-1)*NYDECOMP+NYCUT)
      
        NSTART1D(1)=1
        NPTS=NMEAN
        NEND1D(1)=NSTART1D(1)+NPTS-1

        DO IP=2,NYCUT
          NSTART1D(IP)=NSTART1D(IP-1)+NPTS
          NPTS=NMEAN
          NEND1D(IP)=NSTART1D(IP)+NPTS-1
        ENDDO

        NMEAN=(IJL(NBLO)-NEND1D(NYCUT))/(NYDECOMP-NYCUT)
        NREST=(IJL(NBLO)-NEND1D(NYCUT))-NMEAN*(NYDECOMP-NYCUT)

        DO IP=NYCUT+1,NYDECOMP
          NSTART1D(IP)=NSTART1D(IP-1)+NPTS
          IF(NREST.GT.0) THEN
            NPTS=NMEAN+1
            NREST=NREST-1 
          ELSE
            NPTS=NMEAN
          ENDIF
          NEND1D(IP)=NSTART1D(IP)+NPTS-1
        ENDDO
      ENDIF 


!     SECOND 1-D DECOMPOSITION IN EACH LATITUDINAL BAND

      IF(LL1D.OR.NPR.EQ.1) THEN
!       not needed
        DO IP=1,NYDECOMP
          NSTART(IP)=NSTART1D(IP)
          NEND(IP)=NEND1D(IP)
        ENDDO
      ELSE

        ALLOCATE(NEWIJ2IJ(NIBLO))
        ALLOCATE(IJ2NEWIJ(NIBLO))

        XDELLOINV=1.0/XDELLO

        STAGGER=0.5*(AMOEAP-AMOWEP+IPER*XDELLO)/NXDECOMP
        STAGGER=FLOAT(NINT(100*STAGGER))/100.
        ISTAGGER=NINT(STAGGER*XDELLOINV)

        IG=1
        IPROC=0
        NIJ=0
        DO IPR=1,NYDECOMP
          IPROC=IPROC+1
          NSTART(IPROC)=NIJ+1
          NTOT=NEND1D(IPR)-NSTART1D(IPR)+1

!         find number of points per subarea
          IF(IPR.LE.NYCUT) THEN
            NAREA=NXDECOMP
          ELSE
            NAREA=NXDECOMP-1
          ENDIF 
          ALLOCATE(NTOTSUB(NAREA))

          NMEAN=NTOT/NAREA
          NREST=NTOT-NMEAN*NAREA
          DO IAR=1,NAREA
            IF(NREST.GT.0) THEN
              NTOTSUB(IAR)=NMEAN+1
              NREST=NREST-1
            ELSE
              NTOTSUB(IAR)=NMEAN
            ENDIF
          ENDDO


!         sort sea points in latitudinal bands by increasing longitudes
!         and increasing latitude. Note that we use the fact that the
!         IJ's are already ordered for each latitude 
!         by determining the array IJNDEX which contain the IJ's with
!         increasing longitude and latitude.

          KLATBOT=KXLT(NSTART1D(IPR),IG)
          KLATTOP=KXLT(NEND1D(IPR),IG)
          ALLOCATE(KSTART1(KLATBOT:KLATTOP))
          ALLOCATE(KEND1(KLATBOT:KLATTOP))
          ALLOCATE(NLON(KLATBOT:KLATTOP))
          KXLAT=KLATBOT
          KSTART1(KXLAT) = NSTART1D(IPR)
          DO IJ=NSTART1D(IPR)+1,NEND1D(IPR)
            IF(KXLAT.LT.KXLT(IJ,IG)) THEN
              KXLAT = KXLT(IJ,IG)
              KSTART1(KXLAT) = IJ 
              KEND1(KXLAT-1) = IJ-1
            ENDIF
          ENDDO
          KEND1(KLATTOP)=NEND1D(IPR)

          NLONGMAX=0
          DO KXLAT=KLATBOT,KLATTOP
            NLONGMAX=MAX(KEND1(KXLAT)-KSTART1(KXLAT)+1,NLONGMAX)
          ENDDO

          ALLOCATE(IXLON(NLONGMAX,KLATBOT:KLATTOP))

          IXLONMAX=INT(AMOWEP*XDELLOINV)-1
          DO KXLAT=KLATBOT,KLATTOP
            NLON(KXLAT)=0
          ENDDO
          KXLAT=KLATBOT
          DO IJ=NSTART1D(IPR),NEND1D(IPR)
            IF(KXLAT.LT.KXLT(IJ,IG)) THEN
              KXLAT = KXLT(IJ,IG)
            ENDIF
            NLON(KXLAT)=NLON(KXLAT)+1
            IX = IXLG(IJ,IG)
            JSN= KXLT(IJ,IG)
            XLON=AMOWEP+(IX-1)*ZDELLO(JSN)
            XLON=FLOAT(NINT(100*XLON))/100.
            IXLON(NLON(KXLAT),KXLAT)=NINT(XLON*XDELLOINV)
            IXLONMAX=MAX(IXLONMAX,IXLON(NLON(KXLAT),KXLAT))
          ENDDO

          ALLOCATE(IJNDEX(NTOT))

          ALLOCATE(ILON(KLATBOT:KLATTOP))
          DO KXLAT=KLATBOT,KLATTOP
            ILON(KXLAT)=1
          ENDDO
          JC=0
          KMIN=KLATBOT
          DO WHILE(KMIN.GT.0)
            IXLONMIN=IXLONMAX+1
            KMIN=0
            DO KXLAT=KLATBOT,KLATTOP
              IF(ILON(KXLAT).LE.NLON(KXLAT)) THEN
                IF(IXLON(ILON(KXLAT),KXLAT).LT.IXLONMIN) THEN
                  KMIN=KXLAT
                  IXLONMIN=IXLON(ILON(KXLAT),KXLAT)
                ENDIF
              ENDIF
            ENDDO
            IF(KMIN.GT.0) THEN
              IJ=KSTART1(KMIN)+ILON(KMIN)-1 
              JC=JC+1
              IJNDEX(JC)=IJ
              ILON(KMIN)=ILON(KMIN)+1
            ENDIF
          ENDDO

!         find which points belong to a subarea

          JCS=1
          IF(MOD(IPR,2).EQ.0) THEN
!         staggering
            JCM=1
            DO KXLAT=KLATBOT,KLATTOP
              IIL=1
              DO WHILE (IXLON(MIN(IIL,NLON(KXLAT)),KXLAT).LT.ISTAGGER
     &                  .AND.   IIL.LE.NLON(KXLAT)
     &                  .AND.   NLON(KXLAT).GT.0 )
                IIL=IIL+1
                JCM=JCM+1 
              ENDDO
            ENDDO
          ELSE
            JCM=1
          ENDIF

          IAR=1
          IC=0
          DO JC=JCM,NTOT 
            NIJ=NIJ+1
            IC=IC+1 
            IF(IC.EQ.NTOTSUB(IAR)) THEN
              NEND(IPROC)=NIJ
            ELSEIF(IC.GT.NTOTSUB(IAR)) THEN
              IC=1
              IAR=IAR+1
              IPROC=IPROC+1
              NSTART(IPROC)=NIJ
            ENDIF
            IJ=IJNDEX(JC)
            NEWIJ2IJ(NIJ)=IJ
            IJ2NEWIJ(IJ)=NIJ
          ENDDO

          DO JC=JCS,JCM-1
            NIJ=NIJ+1
            IC=IC+1 
            IF(IC.EQ.NTOTSUB(IAR)) THEN
              NEND(IPROC)=NIJ
            ELSEIF(IC.GT.NTOTSUB(IAR)) THEN
              IC=1
              IAR=IAR+1
              IPROC=IPROC+1
              NSTART(IPROC)=NIJ
            ENDIF
            IJ=IJNDEX(JC)
            NEWIJ2IJ(NIJ)=IJ
            IJ2NEWIJ(IJ)=NIJ
          ENDDO


          DEALLOCATE(KSTART1)
          DEALLOCATE(KEND1)
          DEALLOCATE(NLON)
          DEALLOCATE(ILON)
          DEALLOCATE(IXLON)
          DEALLOCATE(IJNDEX)
          DEALLOCATE(NTOTSUB)

        ENDDO


!       RELABELLING OF THE ARRAYS KLAT KLAT2 KLON DEPTH
!       (also see below for WADV)

        DO IC=1,2
          DO IJ=NSTART(1),NEND(NPR)
            IF(KLON(IJ,IC).GT.0 .AND. KLON(IJ,IC).LE.NIBLO)
     &         KLON(IJ,IC) = IJ2NEWIJ(KLON(IJ,IC))
            IF(KLAT(IJ,IC).GT.0 .AND. KLAT(IJ,IC).LE.NIBLO)
     &         KLAT(IJ,IC) = IJ2NEWIJ(KLAT(IJ,IC))
            IF(KLAT2(IJ,IC).GT.0 .AND. KLAT2(IJ,IC).LE.NIBLO)
     &         KLAT2(IJ,IC) = IJ2NEWIJ(KLAT2(IJ,IC))
          ENDDO
        ENDDO

        ALLOCATE(KDUM(NIBLO,2))
        DO IC=1,2
          DO NIJ=NSTART(1),NEND(NPR)
            KDUM(NIJ,IC)=KLON(NEWIJ2IJ(NIJ),IC)
          ENDDO
          DO NIJ=NSTART(1),NEND(NPR)
            KLON(NIJ,IC)=KDUM(NIJ,IC)
          ENDDO
          DO NIJ=NSTART(1),NEND(NPR)
            KDUM(NIJ,IC)=KLAT(NEWIJ2IJ(NIJ),IC)
          ENDDO
          DO NIJ=NSTART(1),NEND(NPR)
            KLAT(NIJ,IC)=KDUM(NIJ,IC)
          ENDDO
          DO NIJ=NSTART(1),NEND(NPR)
            KDUM(NIJ,IC)=KLAT2(NEWIJ2IJ(NIJ),IC)
          ENDDO
          DO NIJ=NSTART(1),NEND(NPR)
            KLAT2(NIJ,IC)=KDUM(NIJ,IC)
          ENDDO
        ENDDO
        DEALLOCATE(KDUM)

        ALLOCATE(KDUM(NIBLO,IGL))
        DO IG=1,IGL
          DO NIJ=NSTART(1),NEND(NPR)
            KDUM(NIJ,IG)=IXLG(NEWIJ2IJ(NIJ),IG)
          ENDDO
          DO NIJ=NSTART(1),NEND(NPR)
            IXLG(NIJ,IG)=KDUM(NIJ,IG)
          ENDDO
          DO NIJ=NSTART(1),NEND(NPR)
            KDUM(NIJ,IG)=KXLT(NEWIJ2IJ(NIJ),IG)
          ENDDO
          DO NIJ=NSTART(1),NEND(NPR)
            KXLT(NIJ,IG)=KDUM(NIJ,IG)
          ENDDO
        ENDDO
        DEALLOCATE(KDUM)

        ALLOCATE(RDUM(NIBLO,IGL))
        DO IG=1,IGL
          DO NIJ=NSTART(1),NEND(NPR)
            RDUM(NIJ,IG)=DEPTH(NEWIJ2IJ(NIJ),IG)
          ENDDO
          DO NIJ=NSTART(1),NEND(NPR)
            DEPTH(NIJ,IG)=RDUM(NIJ,IG)
          ENDDO
        ENDDO
        DEALLOCATE(RDUM)

        IF (IREFRA.EQ.2) THEN
!!!       nothing has been done if currents are present
!!!       they will need to be reordered AND readdressed for the
!!!       halo points
          WRITE(IU06,*) 'MPDECOMP :  It is probable that  !!!!'
          WRITE(IU06,*) 'MPDECOMP :  U and V must be reordered !!!!'
          WRITE(IU06,*) 'MPDECOMP :  The program has aborted  !!!!'
          CALL ABORT1
        ENDIF

      ENDIF

      DEALLOCATE(NSTART1D)
      DEALLOCATE(NEND1D)


!     3. DETERMINE THE LENGTH OF THE MESSAGE THAT WILL BE EXCHANGED 
!        BETWEEN NEIGHBORING SUB GRID DOMAINS
!        -----------------------------------------------------------

      MAXLEN=0
      DO IP=1,NPR
        MAXLEN=MAX(MAXLEN,NEND(IP)-NSTART(IP)+1)
      ENDDO

!     FIND INDEX AND PE OF THE POINTS IN THE HALO
      MAXPERMLEN=2*MAXLEN+10
      MAXPERMLEN2=2*MAXPERMLEN

      ALLOCATE(IJFROMPE(MAXPERMLEN,NPR))
      ALLOCATE(IPROCFROM(MAXPERMLEN,NPR))
      ALLOCATE(NLENHALO(NPR))

      ALLOCATE(ITEMP(MAXPERMLEN2))

      DO IP=1,NPR
        DO IH=1,MAXPERMLEN
          IPROCFROM(IH,IP)=NPR+1
        ENDDO
      ENDDO

!     DETERMINE IJFROMPE
      DO IP=1,NPR
        IH=0
        DO IC=1,2
          DO IJ=NSTART(IP),NEND(IP)
            IF(  KLON(IJ,IC).GT.0 .AND. KLON(IJ,IC).LE.NIBLO .AND.
     &          (KLON(IJ,IC).LT.NSTART(IP).OR.
     &           KLON(IJ,IC).GT.NEND(IP)       ) ) THEN
               IH=IH+1
               IF(IH.GT.MAXPERMLEN2) THEN
                 WRITE(IU06,*) 'MPDECOMP :  decomposition problem !!!'
                 WRITE(IU06,*) 'MAXPERMLEN2 TOO SMALL !!!'
                 CALL ABORT1
               ENDIF
               ITEMP(IH)=KLON(IJ,IC)
            ENDIF
            IF(  KLAT(IJ,IC).GT.0 .AND. KLAT(IJ,IC).LE.NIBLO .AND.
     &          (KLAT(IJ,IC).LT.NSTART(IP).OR.
     &           KLAT(IJ,IC).GT.NEND(IP)       ) ) THEN 
               IH=IH+1
               IF(IH.GT.MAXPERMLEN2) THEN
                 WRITE(IU06,*) 'MPDECOMP :  decomposition problem !!!'
                 WRITE(IU06,*) 'MAXPERMLEN2 TOO SMALL !!!'
                 CALL ABORT1
               ENDIF
               ITEMP(IH)=KLAT(IJ,IC)
            ENDIF

            IF(  KLAT2(IJ,IC).GT.0 .AND. KLAT2(IJ,IC).LE.NIBLO .AND.
     &          (KLAT2(IJ,IC).LT.NSTART(IP).OR.
     &           KLAT2(IJ,IC).GT.NEND(IP)       ) ) THEN 
               IH=IH+1
               IF(IH.GT.MAXPERMLEN2) THEN
                 WRITE(IU06,*) 'MPDECOMP :  decomposition problem !!!'
                 WRITE(IU06,*) 'MAXPERMLEN2 TOO SMALL !!!'
                 CALL ABORT1
               ENDIF
               ITEMP(IH)=KLAT2(IJ,IC)
            ENDIF
          ENDDO
        ENDDO

        NH=IH
        IF(NH.GT.1) THEN
          ALLOCATE(INDEX(NH))
          CALL SORTINI(ITEMP,INDEX,NH)
          CALL SORTI(ITEMP,INDEX,NH)
          DEALLOCATE(INDEX)
        ENDIF

        JH=1

        IJFROMPE(1,IP)=ITEMP(1)
        DO IH=2,NH
          IF(ITEMP(IH).GT.ITEMP(IH-1))THEN
            JH=JH+1
            IF(JH.GT.MAXPERMLEN) THEN
               WRITE(IU06,*) 'MPDECOMP :  decomposition problem !!!'
               WRITE(IU06,*) 'MAXPERMLEN TOO SMALL !!!'
               WRITE(IU06,*) 'IH = ',IH,' JH = ',JH
               CALL ABORT1
            ENDIF
            IJFROMPE(JH,IP)=ITEMP(IH)
          ENDIF
        ENDDO
        NLENHALO(IP)=JH

      ENDDO
      DEALLOCATE(ITEMP)

!     DETERMINE IPROCFROM, KLENBOT and KLENTOP
      DO IP=1,NPR
        KLENBOT(IP)=0
        KLENTOP(IP)=0
        DO IH=1,NLENHALO(IP)
          DO IPROC=1,NPR
            IF( IJFROMPE(IH,IP).GE.NSTART(IPROC) .AND.
     &          IJFROMPE(IH,IP).LE. NEND(IPROC) ) THEN
              IPROCFROM(IH,IP)=IPROC
              IF(IPROC.LT.IP) THEN
                KLENBOT(IP)=KLENBOT(IP)+1
              ELSEIF(IPROC.GT.IP) THEN
                KLENTOP(IP)=KLENTOP(IP)+1
              ENDIF
              EXIT
            ENDIF
          ENDDO

        ENDDO

      ENDDO

      NINF=NSTART(IRANK)-KLENBOT(IRANK)
      NSUP=NEND(IRANK)+KLENTOP(IRANK)

!     FIND THE LOCAL NTOPE
      DO IP=1,NPR
        NTOPE(IP)=0
      ENDDO
      IF(LMESSPASS .OR. LLIRANK) THEN
        DO IP=1,NPR
          DO IH=1,NLENHALO(IP)
            IF(IPROCFROM(IH,IP).EQ.IRANK) THEN
              NTOPE(IP)=NTOPE(IP)+1
            ENDIF
          ENDDO
        ENDDO

        NTOPEMAX=0
        DO IP=1,NPR
          NTOPEMAX=MAX(NTOPEMAX,NTOPE(IP))
        ENDDO

      ENDIF


!     FIND THE LOCAL NFROMPE
      DO IP=1,NPR
        NFROMPE(IP)=0
      ENDDO
      IF(LMESSPASS .OR. LLIRANK) THEN
        DO IH=1,NLENHALO(IRANK)
          NFROMPE(IPROCFROM(IH,IRANK))=NFROMPE(IPROCFROM(IH,IRANK))+1
        ENDDO

        NFROMPEMAX=0
        DO IP=1,NPR
          NFROMPEMAX=MAX(NFROMPEMAX,NFROMPE(IP))
        ENDDO

      ENDIF


!     FIND NGBTOPE AND CREATE NTOPELST

      NGBTOPE=0
      DO IP=1,NPR
        IF(NTOPE(IP).GT.0) NGBTOPE=NGBTOPE+1
      ENDDO

      ALLOCATE(NTOPELST(NGBTOPE))
      INBNGH=0
      DO IP=1,NPR
        IF(NTOPE(IP).GT.0) THEN 
          INBNGH=INBNGH+1
          NTOPELST(INBNGH)=IP
        ENDIF
      ENDDO

!     FIND NGBFROMPE AND CREATE NFROMPELST 

      NGBFROMPE=0 
      DO IP=1,NPR
        IF(NFROMPE(IP).GT.0) NGBFROMPE=NGBFROMPE+1
      ENDDO

      ALLOCATE(NFROMPELST(NGBFROMPE))
      INBNGH=0
      DO IP=1,NPR
        IF(NFROMPE(IP).GT.0) THEN 
          INBNGH=INBNGH+1
          NFROMPELST(INBNGH)=IP
        ENDIF
      ENDDO


!     DETERMINE WHICH IJ's NEED TO BE SEND TO THE OTHER PE'S
      IF(LMESSPASS .OR. LLIRANK) THEN
        ALLOCATE(IJTOPE(NTOPEMAX,NPR))
        DO IP=1,NPR
          DO JH=1,NTOPEMAX
            IJTOPE(JH,IP)=NINF-1
          ENDDO
        ENDDO
        DO IP=1,NPR
          JH=0
          DO IH=1,NLENHALO(IP)
            IF(IPROCFROM(IH,IP).EQ.IRANK) THEN
              JH=JH+1
              IJTOPE(JH,IP)=IJFROMPE(IH,IP)
            ENDIF
          ENDDO
        ENDDO
      ENDIF


      ALLOCATE(IJHALO(NLENHALO(IRANK)))
      DO IH=1,NLENHALO(IRANK)
        IF(IPROCFROM(IH,IRANK).LT.IRANK) THEN
          IJHALO(IH)=NINF+IH-1
        ELSEIF(IPROCFROM(IH,IRANK).GT.IRANK) THEN
          IJHALO(IH)=NEND(IRANK)+IH-KLENBOT(IRANK)
        ENDIF
      ENDDO


!     CHANGE THE LOCAL ADDRESSING OF KLON,KLAT,KLAT2, DEPTH FOR POINTS
!     IN THE HALO
!     NOTE THAT THIS IMPLIES THAT THESE ARRAYS ARE LOCAL BECAUSE THEY
!     ARE DIFFERENT IN THE HALO REGIONS
      IF((LMESSPASS.AND.NPR.GT.1).OR. LLIRANK) THEN
        DO IC=1,2
          DO IJ=NSTART(IRANK),NEND(IRANK)
            DO IH=1,NLENHALO(IRANK)
              IF(KLON(IJ,IC).EQ.IJFROMPE(IH,IRANK)) THEN
                KLON(IJ,IC)=IJHALO(IH)
                EXIT
              ENDIF
            ENDDO
            DO IH=1,NLENHALO(IRANK)
              IF(KLAT(IJ,IC).EQ.IJFROMPE(IH,IRANK)) THEN
                KLAT(IJ,IC)=IJHALO(IH)
                EXIT
              ENDIF
            ENDDO
            DO IH=1,NLENHALO(IRANK)
              IF(KLAT2(IJ,IC).EQ.IJFROMPE(IH,IRANK)) THEN
                KLAT2(IJ,IC)=IJHALO(IH)
                EXIT
              ENDIF
            ENDDO
          ENDDO
        ENDDO

        ALLOCATE(RDUM(NLENHALO(IRANK),IGL))
        DO IG=1,IGL
          DO IH=1,NLENHALO(IRANK)
            IJ=IJFROMPE(IH,IRANK)
            RDUM(IH,IG)=DEPTH(IJ,IG)
          ENDDO
          DO IH=1,NLENHALO(IRANK)
            IJ=IJHALO(IH)
            DEPTH(IJ,IG)=RDUM(IH,IG)
          ENDDO
        ENDDO
        DEALLOCATE(RDUM)

      ENDIF

!     FIND NIJSTART

      DO IP=1,NPR
        NIJSTART(IP)=NINF-1
      ENDDO
      IF((LMESSPASS.AND.NPR.GT.1) .OR. LLIRANK) THEN
        IF(IPROCFROM(1,IRANK).LT.IRANK) THEN
          NIJSTART(IPROCFROM(1,IRANK))=NINF
        ELSEIF(IPROCFROM(1,IRANK).GT.IRANK) THEN
          NIJSTART(IPROCFROM(1,IRANK))=NEND(IRANK)+1
        ENDIF
        DO IH=2,NLENHALO(IRANK)
          IF(IPROCFROM(IH,IRANK).NE.IPROCFROM(IH-1,IRANK))THEN
            IF(IPROCFROM(IH,IRANK).LT.IRANK) THEN
              NIJSTART(IPROCFROM(IH,IRANK))=NINF+IH-1
            ELSEIF(IPROCFROM(IH,IRANK).GT.IRANK) THEN
              NIJSTART(IPROCFROM(IH,IRANK))=
     &                NEND(IRANK)+IH-KLENBOT(IRANK)
            ENDIF
          ENDIF
        ENDDO
      ENDIF

!     DETERMINE DELLAM1 AND COSPHM1 WHICH ARE USED IN PROPAGS.
!     Note that they will be defined locally for the halo points !!!

      ALLOCATE(DELLAM1(NINF-1:NSUP,IGL))
      ALLOCATE(COSPHM1(NINF-1:NSUP,IGL))
      DO IG=1,IGL
        DELLAM1(NINF-1,IG) = 0.
        IF((LMESSPASS.AND.NPR.GT.1)) THEN
          DO IJ=NSTART(IRANK),NEND(IRANK)
            JH = KXLT(IJ,IG)
            DELLAM1(IJ,IG) = 1./DELLAM(JH)
            COSPHM1(IJ,IG) = 1./COSPH(JH)
          ENDDO 
          DO IH=1,NLENHALO(IRANK)
            IJ=IJFROMPE(IH,IRANK)
            JH = KXLT(IJ,IG)
            IJ=IJHALO(IH)
            DELLAM1(IJ,IG) = 1./DELLAM(JH)
            COSPHM1(IJ,IG) = 1./COSPH(JH)
          ENDDO
        ELSE
          DO IJ=NINF,NSUP
            JH = KXLT(IJ,IG)
            DELLAM1(IJ,IG) = 1./DELLAM(JH)
            COSPHM1(IJ,IG) = 1./COSPH(JH)
          ENDDO 
        ENDIF
      ENDDO

      DEALLOCATE(IJFROMPE)
      DEALLOCATE(IPROCFROM)
      DEALLOCATE(NLENHALO)
      DEALLOCATE(IJHALO)


!     FIND BOUNDARY POINTS FOR EACH AREA
!     !!! it hads to run before 5. (see below) 

      ALLOCATE(LLBOUND(NIBLO))
      DO IJ=1,NIBLO
        LLBOUND(IJ)=.FALSE.
      ENDDO
      DO IP = 1, NPR
        DO IC=1,2
          DO IJ=NSTART(IP),NEND(IP)
            IF( KLON(IJ,IC).GT.0 .AND. KLON(IJ,IC).LE.NIBLO .AND.
     &         (KLON(IJ,IC).LT.NSTART(IP).OR.
     &          KLON(IJ,IC).GT.NEND(IP)       ) ) LLBOUND(IJ)=.TRUE.
           IF(  KLAT(IJ,IC).GT.0 .AND. KLAT(IJ,IC).LE.NIBLO .AND.
     &         (KLAT(IJ,IC).LT.NSTART(IP).OR.
     &          KLAT(IJ,IC).GT.NEND(IP)       ) ) LLBOUND(IJ)=.TRUE.
           IF(  KLAT2(IJ,IC).GT.0 .AND. KLAT2(IJ,IC).LE.NIBLO .AND.
     &         (KLAT2(IJ,IC).LT.NSTART(IP).OR.
     &          KLAT2(IJ,IC).GT.NEND(IP)       ) ) LLBOUND(IJ)=.TRUE.
          ENDDO
        ENDDO
      ENDDO

!     4. KEEP THE PART OF KLAT,KLAT2,KLON,DEPTH  WHICH IS NECESSARY
!        ----------------------------------------------------------

      IF(LMESSPASS.AND.NPR.GT.1) THEN

        ALLOCATE(KDUM(NSTART(IRANK):NEND(IRANK),2))

        DO IC=1,2
          DO IJ=NSTART(IRANK),NEND(IRANK)
            KDUM(IJ,IC)=KLAT(IJ,IC)
          ENDDO
        ENDDO

        DEALLOCATE(KLAT)
        ALLOCATE(KLAT(NSTART(IRANK):NEND(IRANK),2))

        DO IC=1,2
          DO IJ=NSTART(IRANK),NEND(IRANK)
            KLAT(IJ,IC)=KDUM(IJ,IC)
          ENDDO
        ENDDO

        DO IC=1,2
          DO IJ=NSTART(IRANK),NEND(IRANK)
            KDUM(IJ,IC)=KLON(IJ,IC)
          ENDDO
        ENDDO

        DEALLOCATE(KLON)
        ALLOCATE(KLON(NSTART(IRANK):NEND(IRANK),2))

        DO IC=1,2
          DO IJ=NSTART(IRANK),NEND(IRANK)
            KLON(IJ,IC)=KDUM(IJ,IC)
          ENDDO
        ENDDO

        DO IC=1,2
          DO IJ=NSTART(IRANK),NEND(IRANK)
            KDUM(IJ,IC)=KLAT2(IJ,IC)
          ENDDO
        ENDDO

        DEALLOCATE(KLAT2)
        ALLOCATE(KLAT2(NSTART(IRANK):NEND(IRANK),2))

        DO IC=1,2
          DO IJ=NSTART(IRANK),NEND(IRANK)
            KLAT2(IJ,IC)=KDUM(IJ,IC)
          ENDDO
        ENDDO

        DEALLOCATE(KDUM)

        ALLOCATE(RDUM(NINF:NSUP,IGL))
        DO IG=1,IGL
          DO IJ=NINF,NSUP
            RDUM(IJ,IG)=DEPTH(IJ,IG)
          ENDDO
        ENDDO
        DEALLOCATE(DEPTH)
        ALLOCATE(DEPTH(NINF:NSUP,IGL))
        DO IG=1,IGL
          DO IJ=NINF,NSUP
            DEPTH(IJ,IG)=RDUM(IJ,IG)
          ENDDO
        ENDDO
        DEALLOCATE(RDUM)

      ENDIF

!     5. MODIFY KLAT AND KLON SUCH THAT POINT INDICES FOR LAND IS
!        NINF-1.
!        ---------------------------------------------------------

      IF(.NOT.LMESSPASS) NSUP=NIBLO

      IF(LMESSPASS) THEN
        DO IC=1,2
          DO IJ = NSTART(IRANK),NEND(IRANK)
            IF(KLAT(IJ,IC).EQ.0) KLAT(IJ,IC) = NINF-1
            IF(KLAT2(IJ,IC).EQ.0) KLAT2(IJ,IC) = NINF-1
            IF(KLON(IJ,IC).EQ.0) KLON(IJ,IC) = NINF-1
          ENDDO
        ENDDO
      ENDIF


!*    6. INPUT THE WEIGHT FOR THE ADVECTION SCHEME 
!        --------------------------------------------

      ALLOCATE(WADV(NIBLO,2))   ! THE SIZE OF IS READJUSTED SEE BELOW
      ALLOCATE(KOBSLAT(NIBLO,NFRE,2))
      ALLOCATE(KOBSLON(NIBLO,NFRE,2))

      IF(LMESSPASS) THEN

!       READ IOU8 ON PE IREAD
        IF(IRANK.EQ.IREAD) THEN
          READ (IU08) WADV
          READ (IU08) KOBSLAT
          READ (IU08) KOBSLON
        ENDIF

!       SEND WADV TO OTHER PE'S
        IF(NPR.GT.1) THEN
          ITAG=KTAG+1
          MPLENGTH=2*NIBLO
          ALLOCATE(RCOMBUF(MPLENGTH))

          IF(IRANK.EQ.IREAD) THEN
            KCOUNT=0
            DO IC=1,2
              DO IJ=1,NIBLO
                KCOUNT=KCOUNT+1
                RCOMBUF(KCOUNT)=WADV(IJ,IC)
              ENDDO
            ENDDO
          ENDIF

          CALL MPL_BROADCAST(RCOMBUF(1:MPLENGTH),KROOT=IREAD,
     &     KTAG=ITAG,CDSTRING='MPDECOMP 1:')

          IF(IRANK.NE.IREAD) THEN
            KCOUNT=0
            DO IC=1,2
              DO IJ=1,NIBLO
                KCOUNT=KCOUNT+1
                WADV(IJ,IC)=RCOMBUF(KCOUNT)
              ENDDO
            ENDDO
          ENDIF

          DEALLOCATE(RCOMBUF)
        ENDIF

!       SEND KOBSLAT AND KOBSLON TO OTHER PE'S
        IF(NPR.GT.1) THEN
          ITAG=KTAG+2
          MPLENGTH=4*NFRE*NIBLO
          ALLOCATE(ICOMBUF(MPLENGTH))

          IF(IRANK.EQ.IREAD) THEN
            KCOUNT=0
            DO IC=1,2
              DO M=1,NFRE
                DO IJ=1,NIBLO
                  KCOUNT=KCOUNT+1
                  ICOMBUF(KCOUNT)=KOBSLAT(IJ,M,IC)
                ENDDO
                DO IJ=1,NIBLO
                  KCOUNT=KCOUNT+1
                  ICOMBUF(KCOUNT)=KOBSLON(IJ,M,IC)
                ENDDO
              ENDDO
            ENDDO
          ENDIF

          CALL MPL_BROADCAST(ICOMBUF(1:MPLENGTH),KROOT=IREAD,
     &     KTAG=ITAG,CDSTRING='MPDECOMP 2:')

          IF(IRANK.NE.IREAD) THEN
            KCOUNT=0
            DO IC=1,2
              DO M=1,NFRE
                DO IJ=1,NIBLO
                  KCOUNT=KCOUNT+1
                  KOBSLAT(IJ,M,IC)=ICOMBUF(KCOUNT)
                ENDDO
                DO IJ=1,NIBLO
                  KCOUNT=KCOUNT+1
                  KOBSLON(IJ,M,IC)=ICOMBUF(KCOUNT)
                ENDDO
              ENDDO
            ENDDO
          ENDIF

          DEALLOCATE(ICOMBUF)

        ENDIF

      ELSEIF(.NOT.LMESSPASS) THEN
        READ (IU08) WADV
        READ (IU08) KOBSLAT
        READ (IU08) KOBSLON
      ENDIF


!     RELABELLING OF THE ARRAYS

      IF(.NOT.LL1D.AND.NPR.GT.1) THEN
        ALLOCATE(RDUM(NIBLO,2))
        DO IC=1,2
          DO NIJ=NSTART(1),NEND(NPR)
            RDUM(NIJ,IC)=WADV(NEWIJ2IJ(NIJ),IC)
          ENDDO
          DO NIJ=NSTART(1),NEND(NPR)
           WADV(NIJ,IC)=RDUM(NIJ,IC)
          ENDDO
        ENDDO
        DEALLOCATE(RDUM)

        ALLOCATE(KDUM(NIBLO,2))
        DO IC=1,2
          DO M=1,NFRE
            DO NIJ=NSTART(1),NEND(NPR)
              KDUM(NIJ,IC)=KOBSLAT(NEWIJ2IJ(NIJ),M,IC)
            ENDDO
            DO NIJ=NSTART(1),NEND(NPR)
              KOBSLAT(NIJ,M,IC)=KDUM(NIJ,IC)
            ENDDO
          ENDDO
        ENDDO
        DO IC=1,2
          DO M=1,NFRE
            DO NIJ=NSTART(1),NEND(NPR)
              KDUM(NIJ,IC)=KOBSLON(NEWIJ2IJ(NIJ),M,IC)
            ENDDO
            DO NIJ=NSTART(1),NEND(NPR)
              KOBSLON(NIJ,M,IC)=KDUM(NIJ,IC)
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(KDUM)

        DEALLOCATE(NEWIJ2IJ)
      ENDIF

!     KEEP THE PART OF WADV, OBSLAT AND OBSLON WHICH ARE NECESSARY
!     ------------------------------------------------------------

      IF(LMESSPASS.AND.NPR.GT.1) THEN

        ALLOCATE(RDUM(NSTART(IRANK):NEND(IRANK),2))

        DO IC=1,2
          DO IJ=NSTART(IRANK),NEND(IRANK)
            RDUM(IJ,IC)=WADV(IJ,IC)
          ENDDO
        ENDDO
        DEALLOCATE(WADV)
        ALLOCATE(WADV(NSTART(IRANK):NEND(IRANK),2))
        DO IC=1,2
          DO IJ=NSTART(IRANK),NEND(IRANK)
            WADV(IJ,IC)=RDUM(IJ,IC)
          ENDDO
        ENDDO

        DEALLOCATE(RDUM)

      ENDIF

!     NOTE: OBSLON THE VALUE OF OBSLON WILL BE RESET IN THE FIRST
!     CALL TO PROPAGS TO CONTAIN THE OBSTRUCTION TIME THE GROUP
!     VELOCITY At THE INTERFACE !!!!!!!!
      ALLOCATE(OBSLON(NSTART(IRANK):NEND(IRANK),NFRE,2))

      DO IC=1,2
        DO M=1,NFRE
          DO IJ=NSTART(IRANK),NEND(IRANK)
            IF(KOBSLON(IJ,M,IC).EQ.0) THEN
              OBSLON(IJ,M,IC)=0.0
            ELSEIF(MOD(KOBSLON(IJ,M,IC),1000).EQ.0) THEN
              OBSLON(IJ,M,IC)=1.0
            ELSE
              OBSLON(IJ,M,IC)=FLOAT(KOBSLON(IJ,M,IC))*0.001
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(KOBSLON)

!     NOTE: OBSLON THE VALUE OF OBSLAT WILL BE RESET IN THE FIRST
!     CALL TO PROPAGS TO CONTAIN THE OBSTRUCTION TIME THE GROUP
!     VELOCITY At THE INTERFACE !!!!!!!!
      ALLOCATE(OBSLAT(NSTART(IRANK):NEND(IRANK),NFRE,2))
      DO IC=1,2
        DO M=1,NFRE
          DO IJ=NSTART(IRANK),NEND(IRANK)
            IF(KOBSLAT(IJ,M,IC).EQ.0) THEN
              OBSLAT(IJ,M,IC)=0.0
            ELSEIF(MOD(KOBSLAT(IJ,M,IC),1000).EQ.0) THEN
              OBSLAT(IJ,M,IC)=1.0
            ELSE
              OBSLAT(IJ,M,IC)=FLOAT(KOBSLAT(IJ,M,IC))*0.001
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(KOBSLAT)



!     7. DETERMINE KXLTMIN KXLTMAX

      ALLOCATE(KXLTMIN(NPR))
      ALLOCATE(KXLTMAX(NPR))

      DO IP=1,NPR
        KXLTMIN(IP)=NGY
        KXLTMAX(IP)=1
        DO IJ=NSTART(IP),NEND(IP)
           KXLTMIN(IP)=MIN(KXLTMIN(IP),KXLT(IJ,1))
           KXLTMAX(IP)=MAX(KXLTMAX(IP),KXLT(IJ,1))
        ENDDO
      ENDDO

!     8. TEST AND RESET
!        -------------- 
      IF(LMESSPASS) THEN
        IF(IGL.GT.1) THEN
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '*    FATAL ERROR                         *'
          WRITE (IU06,*) '*    ===========                         *'
          WRITE (IU06,*) '* A MULTI BLOCK VERSION RUN WAS REQUESTED*'
          WRITE (IU06,*) '* THIS OPTION IS NOT IMPLEMENTED FOR     *'
          WRITE (IU06,*) '* MESSAGE PASSING PROTOCOL               *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.      *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'
          CALL ABORT1
        ELSE
! SPECIFY THE LIMITS OF THE SEGMENT OVER WHICH THE PE HAS DIRECT ACCESS
          IJS(1)=NSTART(IRANK)
          IJL(1)=NEND(IRANK)
          IJLT(1)=NSUP
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE MPDECOMP
