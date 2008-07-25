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
!    KRLAT, KRLON, KRLAT2,KRLON2,
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
      USE YOWFRED  , ONLY : FR       ,DFIM      ,DFIMOFR  ,
     &             DFFR     ,DFFR2   ,COSTH     ,SINTH
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL      ,IJLT     ,
     &            DELLAM   ,DELLAM1  ,COSPH    ,COSPHM1  ,DELPHI   ,
     &            CDR      ,SDR      ,PRQRT    ,NLONRGG
      USE YOWICE   , ONLY : FLMIN    ,FLMINFR
      USE YOWJONS  , ONLY : AJONS    ,BJONS    ,DJONS    ,EJONS
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,KXLTMIN  ,KXLTMAX  ,
     &            IPER     ,IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,LLBOUND  ,
     &            KMNOP    ,KMSOP    ,IJFROMIK ,IFROMIJ  ,KFROMIJ
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NINF     ,NSUP     ,KTAG     ,
     &                      NPRECR   ,NPRECI   ,MBXSIZE
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO     ,NIBLO    ,
     &            NGX      ,NGY      ,LL1D     ,KWAMVER
      USE YOWPCONS , ONLY : G        ,PI       ,ZPI
      USE YOWSHAL  , ONLY : DEPTH
      USE YOWSTAT  , ONLY : IPROPAGS ,LSUBGRID
      USE YOWSPEC  , ONLY : NSTART   ,NEND     ,KLENTOP  ,KLENBOT  ,
     &            NFROMPE  ,NFROMPEMAX,NTOPE   ,NTOPEMAX ,NIJSTART ,
     &            IJTOPE   ,NGBTOPE  ,NTOPELST ,NGBFROMPE,NFROMPELST,
     &            IJ2NEWIJ
      USE YOWTABL  , ONLY : JUMAX    ,UMAX     ,DELU
      USE YOWTEST  , ONLY : IU06
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,KLAT2    ,KCOR      ,
     &            KRLAT    ,KRLON    ,KRLAT2   ,KRLON2   ,
     &            WLAT     ,WRLAT    ,WRLON    ,
     &            OBSLAT   ,OBSLON   ,OBSCOR   ,OBSRLAT  ,OBSRLON
      USE YOWUNIT  , ONLY : IU07     ,IU08
      USE MPL_MODULE
      USE MPL_DATA_MODULE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK


!----------------------------------------------------------------------


      INTEGER I4(2)
      INTEGER,ALLOCATABLE :: ICOMBUF(:)
      INTEGER,ALLOCATABLE :: KCOMBUF_S(:,:)
      INTEGER,ALLOCATABLE :: KCOMBUF_R(:)
      INTEGER,ALLOCATABLE :: ISENDREQ(:)
      INTEGER,ALLOCATABLE :: KDUM(:), KDUM2(:,:), KDUM3(:,:,:)
      INTEGER,ALLOCATABLE, DIMENSION(:) :: INDEX,IJNDEX,NTOTSUB
      INTEGER,ALLOCATABLE, DIMENSION(:) :: NEWIJ2IJ
      INTEGER,ALLOCATABLE, DIMENSION(:) :: NSTART1D,NEND1D
      INTEGER,ALLOCATABLE, DIMENSION(:) :: KSTART1,KEND1, NLON, ILON
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IXLON
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITEMP,KTEMP,NLENHALO,IJHALO
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IJFROMPE,IPROCFROM
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: KOBSLON, KOBSLAT
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: KOBSCOR
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: KOBSRLON, KOBSRLAT


      REAL X4(2)
      REAL,ALLOCATABLE :: RCOMBUF_S(:,:)
      REAL,ALLOCATABLE :: RCOMBUF_R(:)
      REAL,ALLOCATABLE :: RDUM(:), RDUM2(:,:)

      REAL ZHOOK_HANDLE

      CHARACTER*9 :: CMBXSIZE
      CHARACTER*80 LOGFILENAME
      LOGICAL :: LLIRANK
      LOGICAL :: LLEXIST

      IF (LHOOK) CALL DR_HOOK('MPDECOMP',0,ZHOOK_HANDLE)

!     STANDARD OUTPUT UNIT
!     --------------------
      IF(LMESSPASS.AND.LWCOU) THEN
        IU06=20
      ELSEIF (LMESSPASS.AND..NOT.LWCOU) THEN
        IU06=66
        LOGFILENAME='logfile.%p'
        CALL EXPAND_STRING(IRANK,NPR,0,0,LOGFILENAME,1)
        INQUIRE(UNIT=IU06, OPENED=LLEXIST)
        IF(.NOT. LLEXIST) THEN
          OPEN(IU06,FILE=LOGFILENAME,STATUS='UNKNOWN')
          WRITE(*,*) ' STDOUT FILE ', LOGFILENAME, ' OPENED ON PE ', IRANK
          WRITE(IU06,*) ' STDOUT OF PE ', IRANK
          WRITE(IU06,*) ' ============ '
        ENDIF
        CALL FLUSH (IU06)
        MPL_UNIT=IU06
      ELSE
        IU06=6
      ENDIF

      WRITE(IU06,*) ' WAM SOFTWARE VERSION: ', KWAMVER
      CALL FLUSH (IU06)


!     GRID AND UBUF FILES UNITS
!     -------------------------

      IF(LWCOU) THEN
        IU07=29
        IU08=31
      ELSE
        IU07=7
        IU08=8
      ENDIF

!     DETERMINE BYTE STORAGE REPRESENTATION OF REAL NUMBERS
!     -----------------------------------------------------

      X4=1.
      NPRECR = KIND(X4)
      I4=1
      NPRECI = KIND(I4)

!     GET MESSAGE PASSING BUFFER SIZE
      CALL GETENV('MBX_SIZE',CMBXSIZE)
      READ(CMBXSIZE,'(I9)') MBXSIZE

      ALLOCATE(ISENDREQ(MAX(1,NPR-1)))

!     0. READ GRID INPUT FROM PREPROC 
!        ----------------------------

      IREAD=1
      KTAG=1

      IF((LMESSPASS.AND.IRANK.EQ.IREAD).OR.(.NOT.LMESSPASS)) THEN
        CALL GSTATS(1771,0)
        CALL READPRE (IU07)
        CLOSE (UNIT=IU07, STATUS='KEEP')
        CALL GSTATS(1771,1)
      ENDIF

!     SEND INFORMATION FROM READPRE TO ALL PE's
      CALL GSTATS(694,0)
      IF(LMESSPASS) CALL MPBCASTGRID(IU06,IREAD,KTAG)
      CALL GSTATS(694,1)

      WRITE(IU06,*) ' PREPROC GRID INFORMATION READ IN '
      CALL FLUSH (IU06)

!     DEFINE COEFFICIENT FOR MEAN PERIODS CALCULATION
      IF(ALLOCATED(DFIMOFR)) DEALLOCATE(DFIMOFR)
      ALLOCATE(DFIMOFR(NFRE))
      IF(ALLOCATED(DFFR)) DEALLOCATE(DFFR)
      ALLOCATE(DFFR(NFRE))
      IF(ALLOCATED(DFFR2)) DEALLOCATE(DFFR2)
      ALLOCATE(DFFR2(NFRE))

      DO M=1,NFRE
        DFIMOFR(M) = DFIM(M)/FR(M)
        DFFR(M) = DFIM(M)*FR(M)
        DFFR2(M) = DFIM(M)*FR(M)**2
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

      FLMIN=0.0001

      IF(ALLOCATED(FLMINFR)) DEALLOCATE(FLMINFR)
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

      IF(ALLOCATED(KLAT)) DEALLOCATE(KLAT)
      ALLOCATE(KLAT(NIBLO,2))  ! THE SIZE OF KLAT IS READJUSTED SEE BELOW

      IF(ALLOCATED(KLON)) DEALLOCATE(KLON)
      ALLOCATE(KLON(NIBLO,2))  ! THE SIZE OF KLON IS READJUSTED SEE BELOW

      IF(ALLOCATED(KLAT2)) DEALLOCATE(KLAT2)
      ALLOCATE(KLAT2(NIBLO,2)) ! THE SIZE OF KLAT2 IS READJUSTED SEE BELOW

      IF(ALLOCATED(KCOR)) DEALLOCATE(KCOR)
      ALLOCATE(KCOR(NIBLO,4,-1:1)) ! THE SIZE OF KCOR IS READJUSTED SEE BELOW
                                   ! WILL ONLY BE KEPT IF IPROPAGS=2

      IF(ALLOCATED(KRLAT)) DEALLOCATE(KRLAT)
      ALLOCATE(KRLAT(NIBLO,2))  ! THE SIZE IS READJUSTED SEE BELOW
                                ! WILL ONLY BE KEPT IF IPROPAGS=1
      IF(ALLOCATED(KRLON)) DEALLOCATE(KRLON)
      ALLOCATE(KRLON(NIBLO,2))  ! THE SIZE IS READJUSTED SEE BELOW

      IF(ALLOCATED(KRLAT2)) DEALLOCATE(KRLAT2)
      ALLOCATE(KRLAT2(NIBLO,2)) ! THE SIZE IS READJUSTED SEE BELOW

      IF(ALLOCATED(KRLON2)) DEALLOCATE(KRLON2)
      ALLOCATE(KRLON2(NIBLO,2)) ! THE SIZE IS READJUSTED SEE BELOW


      IF(LMESSPASS) THEN

!       READ FIRST PART OF IOU8 ON PE IREAD
        IF(IRANK.EQ.IREAD) THEN
          READ (IU08) KLAT, KLON
          READ (IU08) KLAT2
          READ (IU08) KRLAT, KRLON
          READ (IU08) KRLAT2, KRLON2
          READ (IU08) KCOR
        ENDIF

!       SEND KLAT AND KLON TO OTHER PE'S
        IF(NPR.GT.1) THEN
          ITAG=KTAG
          IF(IPROPAGS.EQ.2) THEN
            MPLENGTH=(6+12)*NIBLO
          ELSE IF(IPROPAGS.EQ.1) THEN
            MPLENGTH=(6+8)*NIBLO
          ELSE
            MPLENGTH=6*NIBLO
          ENDIF
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

            IF(IPROPAGS.EQ.2) THEN
              DO ICL=-1,1
                DO IC=1,4
                  DO IJ=1,NIBLO
                    KCOUNT=KCOUNT+1
                    ICOMBUF(KCOUNT)=KCOR(IJ,IC,ICL)
                  ENDDO
                ENDDO
              ENDDO
            ELSE IF(IPROPAGS.EQ.1) THEN
              DO IC=1,2
                DO IJ=1,NIBLO
                  KCOUNT=KCOUNT+1
                  ICOMBUF(KCOUNT)=KRLAT(IJ,IC)
                  KCOUNT=KCOUNT+1
                  ICOMBUF(KCOUNT)=KRLAT2(IJ,IC)
                ENDDO
              ENDDO
              DO IC=1,2
                DO IJ=1,NIBLO
                  KCOUNT=KCOUNT+1
                  ICOMBUF(KCOUNT)=KRLON(IJ,IC)
                  KCOUNT=KCOUNT+1
                  ICOMBUF(KCOUNT)=KRLON2(IJ,IC)
                ENDDO
              ENDDO
            ENDIF
          ENDIF

          CALL GSTATS(694,0)
          CALL MPL_BROADCAST(ICOMBUF(1:MPLENGTH),KROOT=IREAD,
     &     KTAG=ITAG,CDSTRING='MPDECOMP:')
          CALL GSTATS(694,1)

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
            IF(IPROPAGS.EQ.2) THEN
              DO ICL=-1,1
                DO IC=1,4
                  DO IJ=1,NIBLO
                    KCOUNT=KCOUNT+1
                    KCOR(IJ,IC,ICL)=ICOMBUF(KCOUNT)
                  ENDDO
                ENDDO
              ENDDO
            ELSE IF(IPROPAGS.EQ.1) THEN
              DO IC=1,2
                DO IJ=1,NIBLO
                  KCOUNT=KCOUNT+1
                  KRLAT(IJ,IC)=ICOMBUF(KCOUNT)
                  KCOUNT=KCOUNT+1
                  KRLAT2(IJ,IC)=ICOMBUF(KCOUNT)
                ENDDO
              ENDDO
              DO IC=1,2
                DO IJ=1,NIBLO
                  KCOUNT=KCOUNT+1
                  KRLON(IJ,IC)=ICOMBUF(KCOUNT)
                  KCOUNT=KCOUNT+1
                  KRLON2(IJ,IC)=ICOMBUF(KCOUNT)
                ENDDO
              ENDDO
            ENDIF

          ENDIF

          DEALLOCATE(ICOMBUF)

        ENDIF

      ELSEIF(.NOT.LMESSPASS) THEN
        READ (IU08) KLAT, KLON
        READ (IU08) KLAT2
        READ (IU08) KRLAT, KRLON
        READ (IU08) KRLAT2, KRLON2
        READ (IU08) KCOR
      ENDIF

!     REMOVE ARRAYS THAT WILL NOT BE USED IN THE PROPAGATION SCHEME
      IF(IPROPAGS.NE.2) THEN
        DEALLOCATE(KCOR)
      ENDIF
      IF(IPROPAGS.NE.1) THEN
        DEALLOCATE(KRLAT)
        DEALLOCATE(KRLON)
        DEALLOCATE(KRLAT2)
        DEALLOCATE(KRLON2)
      ENDIF

!*    2. FIND THE NUMBER OF POINTS PER PE, THE START AND END INDEX
!        --------------------------------------------------------------

      IF (ALLOCATED(NSTART)) DEALLOCATE(NSTART)
      ALLOCATE (NSTART(NPR))
      IF (ALLOCATED(NEND)) DEALLOCATE(NEND)
      ALLOCATE (NEND(NPR))

      IF (ALLOCATED(KLENBOT)) DEALLOCATE(KLENBOT)
      ALLOCATE (KLENBOT(NPR))
      IF (ALLOCATED(KLENTOP)) DEALLOCATE(KLENTOP)
      ALLOCATE (KLENTOP(NPR))

      IF (ALLOCATED(NFROMPE)) DEALLOCATE(NFROMPE)
      ALLOCATE (NFROMPE(NPR))
      IF (ALLOCATED(NTOPE)) DEALLOCATE(NTOPE)
      ALLOCATE (NTOPE(NPR))
      IF (ALLOCATED(NIJSTART)) DEALLOCATE(NIJSTART)
      ALLOCATE (NIJSTART(NPR))

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

      IF (ALLOCATED(NSTART1D)) DEALLOCATE(NSTART1D)
      ALLOCATE(NSTART1D(NYDECOMP))
      IF (ALLOCATED(NEND1D)) DEALLOCATE(NEND1D)
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

        ALLOCATE(NEWIJ2IJ(0:NIBLO))

        IF(ALLOCATED(IJ2NEWIJ)) DEALLOCATE(IJ2NEWIJ)
        ALLOCATE(IJ2NEWIJ(0:NIBLO))
        NEWIJ2IJ(0)=0
        IJ2NEWIJ(0)=0

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


!       RELABELLING OF THE ARRAYS KLAT KLAT2 KLON KCOR KRLAT KRLAT2 
!       KRLON KRLON2 DEPTH IXLG KXLT
!       (also see below for WLAT)

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

        IF(IPROPAGS.EQ.1) THEN
          DO IC=1,2
            DO IJ=NSTART(1),NEND(NPR)
              IF(KRLON(IJ,IC).GT.0 .AND. KRLON(IJ,IC).LE.NIBLO)
     &           KRLON(IJ,IC) = IJ2NEWIJ(KRLON(IJ,IC))
              IF(KRLAT(IJ,IC).GT.0 .AND. KRLAT(IJ,IC).LE.NIBLO)
     &           KRLAT(IJ,IC) = IJ2NEWIJ(KRLAT(IJ,IC))
              IF(KRLAT2(IJ,IC).GT.0 .AND. KRLAT2(IJ,IC).LE.NIBLO)
     &           KRLAT2(IJ,IC) = IJ2NEWIJ(KRLAT2(IJ,IC))
              IF(KRLON2(IJ,IC).GT.0 .AND. KRLON2(IJ,IC).LE.NIBLO)
     &           KRLON2(IJ,IC) = IJ2NEWIJ(KRLON2(IJ,IC))
            ENDDO
          ENDDO
        ELSE IF(IPROPAGS.EQ.2) THEN
          DO ICL=-1,1
            DO IC=1,4
              DO IJ=NSTART(1),NEND(NPR)
                IF(KCOR(IJ,IC,ICL).GT.0 .AND. KCOR(IJ,IC,ICL).LE.NIBLO)
     &             KCOR(IJ,IC,ICL) = IJ2NEWIJ(KCOR(IJ,IC,ICL))
              ENDDO
            ENDDO
          ENDDO
        ENDIF

        ALLOCATE(KDUM(NIBLO))
        DO IC=1,2
          DO NIJ=NSTART(1),NEND(NPR)
            KDUM(NIJ)=KLON(NEWIJ2IJ(NIJ),IC)
          ENDDO
          DO NIJ=NSTART(1),NEND(NPR)
            KLON(NIJ,IC)=KDUM(NIJ)
          ENDDO
          DO NIJ=NSTART(1),NEND(NPR)
            KDUM(NIJ)=KLAT(NEWIJ2IJ(NIJ),IC)
          ENDDO
          DO NIJ=NSTART(1),NEND(NPR)
            KLAT(NIJ,IC)=KDUM(NIJ)
          ENDDO
          DO NIJ=NSTART(1),NEND(NPR)
            KDUM(NIJ)=KLAT2(NEWIJ2IJ(NIJ),IC)
          ENDDO
          DO NIJ=NSTART(1),NEND(NPR)
            KLAT2(NIJ,IC)=KDUM(NIJ)
          ENDDO
        ENDDO
        DEALLOCATE(KDUM)

        IF(IPROPAGS.EQ.1) THEN
          ALLOCATE(KDUM(NIBLO))
          DO IC=1,2
            DO NIJ=NSTART(1),NEND(NPR)
              KDUM(NIJ)=KRLON(NEWIJ2IJ(NIJ),IC)
            ENDDO
            DO NIJ=NSTART(1),NEND(NPR)
              KRLON(NIJ,IC)=KDUM(NIJ)
            ENDDO
            DO NIJ=NSTART(1),NEND(NPR)
              KDUM(NIJ)=KRLAT(NEWIJ2IJ(NIJ),IC)
            ENDDO
            DO NIJ=NSTART(1),NEND(NPR)
              KRLAT(NIJ,IC)=KDUM(NIJ)
            ENDDO
            DO NIJ=NSTART(1),NEND(NPR)
              KDUM(NIJ)=KRLON2(NEWIJ2IJ(NIJ),IC)
            ENDDO
            DO NIJ=NSTART(1),NEND(NPR)
              KRLON2(NIJ,IC)=KDUM(NIJ)
            ENDDO
            DO NIJ=NSTART(1),NEND(NPR)
              KDUM(NIJ)=KRLAT2(NEWIJ2IJ(NIJ),IC)
            ENDDO
            DO NIJ=NSTART(1),NEND(NPR)
              KRLAT2(NIJ,IC)=KDUM(NIJ)
            ENDDO
          ENDDO
          DEALLOCATE(KDUM)
        ELSE IF(IPROPAGS.EQ.2) THEN
          ALLOCATE(KDUM(NIBLO))
          DO ICL=-1,1
            DO IC=1,4
              DO NIJ=NSTART(1),NEND(NPR)
                KDUM(NIJ)=KCOR(NEWIJ2IJ(NIJ),IC,ICL)
              ENDDO
              DO NIJ=NSTART(1),NEND(NPR)
                KCOR(NIJ,IC,ICL)=KDUM(NIJ)
              ENDDO
            ENDDO
          ENDDO
          DEALLOCATE(KDUM)
        ENDIF

        ALLOCATE(KDUM(NIBLO))
        DO IG=1,IGL
          DO NIJ=NSTART(1),NEND(NPR)
            KDUM(NIJ)=IXLG(NEWIJ2IJ(NIJ),IG)
          ENDDO
          DO NIJ=NSTART(1),NEND(NPR)
            IXLG(NIJ,IG)=KDUM(NIJ)
          ENDDO
          DO NIJ=NSTART(1),NEND(NPR)
            KDUM(NIJ)=KXLT(NEWIJ2IJ(NIJ),IG)
          ENDDO
          DO NIJ=NSTART(1),NEND(NPR)
            KXLT(NIJ,IG)=KDUM(NIJ)
          ENDDO
        ENDDO
        DEALLOCATE(KDUM)

        ALLOCATE(RDUM(NIBLO))
        DO IG=1,IGL
          DO NIJ=NSTART(1),NEND(NPR)
            RDUM(NIJ)=DEPTH(NEWIJ2IJ(NIJ),IG)
          ENDDO
          DO NIJ=NSTART(1),NEND(NPR)
            DEPTH(NIJ,IG)=RDUM(NIJ)
          ENDDO
        ENDDO
        DEALLOCATE(RDUM)

      ENDIF ! END IF LL1D

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
      MAXPERMLEN=2*MAXLEN+12
      IF(IPROPAGS.EQ.0) THEN
        MAXPERMLEN2=6*MAXPERMLEN
      ELSEIF(IPROPAGS.EQ.1) THEN
        MAXPERMLEN2=8*MAXPERMLEN
      ELSE
        MAXPERMLEN2=12*MAXPERMLEN
      ENDIF

      ALLOCATE(IJFROMPE(MAXPERMLEN,NPR))
      ALLOCATE(IPROCFROM(MAXPERMLEN,NPR))
      ALLOCATE(NLENHALO(NPR))

      ALLOCATE(ITEMP(MAXPERMLEN2))
      DO IH=1,MAXPERMLEN2
        ITEMP(IH)=0
      ENDDO

      DO IP=1,NPR
        DO IH=1,MAXPERMLEN
          IJFROMPE(IH,IP)=0
          IPROCFROM(IH,IP)=NPR+1
        ENDDO
      ENDDO

!     DETERMINE IJFROMPE
      DO IP=1,NPR
        IH=0

!       CONTRIBUTION FROM NON-ROTATED GRID (NEEDED FOR ALL OPTIONS)
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

          ENDDO ! END DO ON IJ
        ENDDO ! END DO ON IC

        IF (IPROPAGS.EQ.2) THEN
!       CONTRIBUTION FROM CORNER GRID POINT
!       (NEEDED FOR CTU SCHEME: IPROPAGS=2)
          DO ICL=-1,1
            DO IC=1,4
              DO IJ=NSTART(IP),NEND(IP)

                IF(  KCOR(IJ,IC,ICL).GT.0 .AND.
     &               KCOR(IJ,IC,ICL).LE.NIBLO .AND.
     &              (KCOR(IJ,IC,ICL).LT.NSTART(IP).OR.
     &               KCOR(IJ,IC,ICL).GT.NEND(IP)       ) ) THEN
                   IH=IH+1
                   IF(IH.GT.MAXPERMLEN2) THEN
                     WRITE(IU06,*) 'MPDECOMP : DECOMPOSITION PROBLEM !'
                     WRITE(IU06,*) 'MAXPERMLEN2 TOO SMALL !!!'
                     CALL ABORT1
                   ENDIF
                   ITEMP(IH)=KCOR(IJ,IC,ICL)
                ENDIF

              ENDDO ! END DO ON IJ
            ENDDO ! END DO ON IC
          ENDDO ! END DO ON ICL
        ELSE IF (IPROPAGS.EQ.1) THEN
!       CONTRIBUTION FROM ROTATED GRID
!      (NEEDED FOR DUAL ROTATED SCHEME SCHEME: IPROPAGS=1)
          DO IC=1,2
            DO IJ=NSTART(IP),NEND(IP)

              IF(  KRLON(IJ,IC).GT.0 .AND. KRLON(IJ,IC).LE.NIBLO .AND.
     &            (KRLON(IJ,IC).LT.NSTART(IP).OR.
     &             KRLON(IJ,IC).GT.NEND(IP)       ) ) THEN
                 IH=IH+1
                 IF(IH.GT.MAXPERMLEN2) THEN
                   WRITE(IU06,*) 'MPDECOMP :  decomposition problem !!'
                   WRITE(IU06,*) 'MAXPERMLEN2 TOO SMALL !!!'
                   CALL ABORT1
                 ENDIF
                 ITEMP(IH)=KRLON(IJ,IC)
              ENDIF

              IF(  KRLAT(IJ,IC).GT.0 .AND. KRLAT(IJ,IC).LE.NIBLO .AND.
     &            (KRLAT(IJ,IC).LT.NSTART(IP).OR.
     &             KRLAT(IJ,IC).GT.NEND(IP)       ) ) THEN 
                 IH=IH+1
                 IF(IH.GT.MAXPERMLEN2) THEN
                   WRITE(IU06,*) 'MPDECOMP :  decomposition problem !!'
                   WRITE(IU06,*) 'MAXPERMLEN2 TOO SMALL !!!'
                   CALL ABORT1
                 ENDIF
                 ITEMP(IH)=KRLAT(IJ,IC)
              ENDIF

              IF(  KRLON2(IJ,IC).GT.0 .AND. KRLON2(IJ,IC).LE.NIBLO .AND.
     &            (KRLON2(IJ,IC).LT.NSTART(IP).OR.
     &             KRLON2(IJ,IC).GT.NEND(IP)       ) ) THEN 
                 IH=IH+1
                 IF(IH.GT.MAXPERMLEN2) THEN
                   WRITE(IU06,*) 'MPDECOMP :  decomposition problem !!'
                   WRITE(IU06,*) 'MAXPERMLEN2 TOO SMALL !!!'
                   CALL ABORT1
                 ENDIF
                 ITEMP(IH)=KRLON2(IJ,IC)
              ENDIF

              IF(  KRLAT2(IJ,IC).GT.0 .AND. KRLAT2(IJ,IC).LE.NIBLO .AND.
     &            (KRLAT2(IJ,IC).LT.NSTART(IP).OR.
     &             KRLAT2(IJ,IC).GT.NEND(IP)       ) ) THEN 
                 IH=IH+1
                 IF(IH.GT.MAXPERMLEN2) THEN
                   WRITE(IU06,*) 'MPDECOMP :  decomposition problem !!'
                   WRITE(IU06,*) 'MAXPERMLEN2 TOO SMALL !!!'
                   CALL ABORT1
                 ENDIF
                 ITEMP(IH)=KRLAT2(IJ,IC)
              ENDIF

            ENDDO ! END DO ON IJ
          ENDDO ! END DO ON IC
        ENDIF

        NH=IH

        IF(NH.GT.1) THEN
          ALLOCATE(INDEX(NH))
          CALL SORTINI(ITEMP,INDEX,NH)
          CALL SORTI(ITEMP,INDEX,NH)
          DEALLOCATE(INDEX)
        ENDIF

        JH=0
        IF(NH.GT.1) THEN
          JH=1
          IJFROMPE(JH,IP)=ITEMP(1)
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
        ENDIF
        NLENHALO(IP)=JH

      ENDDO ! END DO ON IP
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

      IF (ALLOCATED(NTOPELST)) DEALLOCATE(NTOPELST)
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

      IF (ALLOCATED(NFROMPELST)) DEALLOCATE(NFROMPELST)
      ALLOCATE(NFROMPELST(MAX(1,NGBFROMPE)))
      INBNGH=0
      DO IP=1,NPR
        IF(NFROMPE(IP).GT.0) THEN 
          INBNGH=INBNGH+1
          NFROMPELST(INBNGH)=IP
        ENDIF
      ENDDO


!     DETERMINE WHICH IJ's NEED TO BE SEND TO THE OTHER PE'S
      IF(LMESSPASS .OR. LLIRANK) THEN
        IF (ALLOCATED(IJTOPE)) DEALLOCATE(IJTOPE)
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

      ALLOCATE(IJHALO(MAX(1,NLENHALO(IRANK))))

      DO IH=1,NLENHALO(IRANK)
        IF(IPROCFROM(IH,IRANK).LT.IRANK) THEN
          IJHALO(IH)=NINF+IH-1
        ELSEIF(IPROCFROM(IH,IRANK).GT.IRANK) THEN
          IJHALO(IH)=NEND(IRANK)+IH-KLENBOT(IRANK)
        ENDIF
      ENDDO


!     CHANGE THE LOCAL ADDRESSING OF KLAT KLAT2 KLON KRLAT KRLAT2 
!     KRLON KRLON2 DEPTH FOR POINTS IN THE HALO
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

        IF(IPROPAGS.EQ.2) THEN
!       (NEEDED FOR CTU SCHEME: IPROPAGS=2)
          DO ICL=-1,1
            DO IC=1,4
              DO IJ=NSTART(IRANK),NEND(IRANK)
                DO IH=1,NLENHALO(IRANK)
                  IF(KCOR(IJ,IC,ICL).EQ.IJFROMPE(IH,IRANK)) THEN
                    KCOR(IJ,IC,ICL)=IJHALO(IH)
                    EXIT
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ELSE IF(IPROPAGS.EQ.1) THEN
!       (NEEDED FOR DUAL ROTATED SCHEME SCHEME: IPROPAGS=1)
          DO IC=1,2
            DO IJ=NSTART(IRANK),NEND(IRANK)
              DO IH=1,NLENHALO(IRANK)
                IF(KRLON(IJ,IC).EQ.IJFROMPE(IH,IRANK)) THEN
                  KRLON(IJ,IC)=IJHALO(IH)
                  EXIT
                ENDIF
              ENDDO
              DO IH=1,NLENHALO(IRANK)
                IF(KRLAT(IJ,IC).EQ.IJFROMPE(IH,IRANK)) THEN
                  KRLAT(IJ,IC)=IJHALO(IH)
                  EXIT
                ENDIF
              ENDDO
              DO IH=1,NLENHALO(IRANK)
                IF(KRLON2(IJ,IC).EQ.IJFROMPE(IH,IRANK)) THEN
                  KRLON2(IJ,IC)=IJHALO(IH)
                  EXIT
                ENDIF
              ENDDO
              DO IH=1,NLENHALO(IRANK)
                IF(KRLAT2(IJ,IC).EQ.IJFROMPE(IH,IRANK)) THEN
                  KRLAT2(IJ,IC)=IJHALO(IH)
                  EXIT
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF

        ALLOCATE(RDUM(MAX(1,NLENHALO(IRANK))))
        DO IG=1,IGL
          DO IH=1,NLENHALO(IRANK)
            IJ=IJFROMPE(IH,IRANK)
            RDUM(IH)=DEPTH(IJ,IG)
          ENDDO
          DO IH=1,NLENHALO(IRANK)
            IJ=IJHALO(IH)
            DEPTH(IJ,IG)=RDUM(IH)
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

!     DETERMINE DELLAM1 COSPHM1
!     CDR AND SDR ARE ALSO USED IN PROPAGS1.
!     Note that they all will be defined locally for the halo points !!!

      IF (ALLOCATED(DELLAM1)) DEALLOCATE(DELLAM1)
      ALLOCATE(DELLAM1(NINF-1:NSUP,IGL))
      IF (ALLOCATED(COSPHM1)) DEALLOCATE(COSPHM1)
      ALLOCATE(COSPHM1(NINF-1:NSUP,IGL))

      IF(IPROPAGS.EQ.1) THEN
        IF (ALLOCATED(CDR)) DEALLOCATE(CDR)
        ALLOCATE(CDR(NINF-1:NSUP,NANG,IGL))
        IF (ALLOCATED(SDR)) DEALLOCATE(SDR)
        ALLOCATE(SDR(NINF-1:NSUP,NANG,IGL))
!       PRQRT IS NOT DEFINED IN THE HALO !!!!
        IF (ALLOCATED(PRQRT)) DEALLOCATE(PRQRT)
        ALLOCATE(PRQRT(NSTART(IRANK):NEND(IRANK),IGL))
      ENDIF

      SQRT2O2=SIN(0.25*PI)

      DO IG=1,IGL
        DELLAM1(NINF-1,IG) = 0.
        IF((LMESSPASS.AND.NPR.GT.1)) THEN
          DO IJ=NSTART(IRANK),NEND(IRANK)
            JH = KXLT(IJ,IG)
            DELLAM1(IJ,IG)=1./DELLAM(JH)
            COSPHM1(IJ,IG)=1./COSPH(JH)
          ENDDO 
          DO IH=1,NLENHALO(IRANK)
            IJ=IJFROMPE(IH,IRANK)
            JH = KXLT(IJ,IG)
            IJ=IJHALO(IH)
            DELLAM1(IJ,IG)=1./DELLAM(JH)
            COSPHM1(IJ,IG)=1./COSPH(JH)
          ENDDO
        ELSE
          DO IJ=NINF,NSUP
            JH = KXLT(IJ,IG)
            DELLAM1(IJ,IG)=1./DELLAM(JH)
            COSPHM1(IJ,IG)=1./COSPH(JH)
          ENDDO 
        ENDIF
      ENDDO

      IF(IPROPAGS.EQ.1) THEN
        DO IG=1,IGL
          DO K=1,NANG
            CDR(NINF-1,K,IG) = 0.
            SDR(NINF-1,K,IG) = 0.
          ENDDO
          IF((LMESSPASS.AND.NPR.GT.1)) THEN
            DO IJ=NSTART(IRANK),NEND(IRANK)
              JH = KXLT(IJ,IG)
              DO K=1,NANG
                A=SQRT2O2*COSTH(K)*COSPH(JH)
                B=SQRT2O2*SINTH(K) 
                CDR(IJ,K,IG)=A-B
                SDR(IJ,K,IG)=A+B
              ENDDO
            ENDDO 
            IF(IRGG.EQ.1) THEN
              DO IJ=NSTART(IRANK),NEND(IRANK)
                JH = KXLT(IJ,IG)
                THETAMAX=ATAN2(1.,COSPH(JH))
                SINTHMAX=SIN(THETAMAX)
                DELTA=COSPH(JH)/(SINTHMAX*(1.+COSPH(JH)**2))
                PRQRT(IJ,IG)=MIN(0.5,DELTA)
              ENDDO
            ELSE
              DO IJ=NSTART(IRANK),NEND(IRANK)
                PRQRT(IJ,IG)=0.5
              ENDDO
            ENDIF
            DO IH=1,NLENHALO(IRANK)
              IJ=IJFROMPE(IH,IRANK)
              JH = KXLT(IJ,IG)
              IJ=IJHALO(IH)
              DO K=1,NANG
                A=SQRT2O2*COSTH(K)*COSPH(JH)
                B=SQRT2O2*SINTH(K) 
                CDR(IJ,K,IG)=A-B
                SDR(IJ,K,IG)=A+B
              ENDDO
            ENDDO
          ELSE
            DO IJ=NINF,NSUP
              JH = KXLT(IJ,IG)
              DO K=1,NANG
                A=SQRT2O2*COSTH(K)*COSPH(JH)
                B=SQRT2O2*SINTH(K) 
                CDR(IJ,K,IG)=A-B
                SDR(IJ,K,IG)=A+B
              ENDDO
            ENDDO 
            IF(IRGG.EQ.1) THEN
              DO IJ=NSTART(IRANK),NEND(IRANK)
                JH = KXLT(IJ,IG)
                THETAMAX=ATAN2(1.,COSPH(JH))
                SINTHMAX=SIN(THETAMAX)
                DELTA=COSPH(JH)/(SINTHMAX*(1.+COSPH(JH)**2))
                PRQRT(IJ,IG)=MIN(0.5,DELTA)
              ENDDO
            ELSE
              DO IJ=NSTART(IRANK),NEND(IRANK)
                PRQRT(IJ,IG)=0.5
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDIF

!     CREATE IJFROMIK, IFROMIJ, KFROMIJ (FOR CTUW)
!     !!!! IT IS ONLY DEFINED FOR GRID POINTS ON A GIVEN PE AND THEIR HALO !!!!

      IF(ALLOCATED(IJFROMIK)) DEALLOCATE(IJFROMIK)
      ALLOCATE(IJFROMIK(-1:NGX+2,0:NGY+1,IGL))
      IF(ALLOCATED(IFROMIJ)) DEALLOCATE(IFROMIJ)
      ALLOCATE(IFROMIJ(NINF-1:NSUP))
      IF(ALLOCATED(KFROMIJ)) DEALLOCATE(KFROMIJ)
      ALLOCATE(KFROMIJ(NINF-1:NSUP))

      IJFROMIK=NINF-1
      IFROMIJ(NINF-1)=0
      KFROMIJ(NINF-1)=0

      IG=1

      IF((LMESSPASS.AND.NPR.GT.1)) THEN
!       LOCAL POINTS
        DO IJ=NSTART(IRANK),NEND(IRANK)
          IJFROMIK(IXLG(IJ,IG),KXLT(IJ,IG),IG)=IJ
          IFROMIJ(IJ)=IXLG(IJ,IG)
          KFROMIJ(IJ)=KXLT(IJ,IG)
        ENDDO
!       HALO
        DO IH=1,NLENHALO(IRANK)
          IJ=IJFROMPE(IH,IRANK)
          IJFROMIK(IXLG(IJ,IG),KXLT(IJ,IG),IG)=IJHALO(IH)
          IFROMIJ(IJHALO(IH))=IXLG(IJ,IG)
          KFROMIJ(IJHALO(IH))=KXLT(IJ,IG)
        ENDDO
!       PERIODICITY
        IF(IPER.EQ.1) THEN
          DO IJ=NSTART(IRANK),NEND(IRANK)
            NL=NLONRGG(KXLT(IJ,IG))
            IF(IXLG(IJ,IG).EQ.1) THEN
              IJFROMIK(IXLG(IJ,IG)+NL,KXLT(IJ,IG),IG)=IJ
            ELSEIF(IXLG(IJ,IG).EQ.2) THEN
              IJFROMIK(IXLG(IJ,IG)+NL,KXLT(IJ,IG),IG)=IJ
            ELSE IF(IXLG(IJ,IG).EQ.NL) THEN
              IJFROMIK(IXLG(IJ,IG)-NL,KXLT(IJ,IG),IG)=IJ
            ELSE IF(IXLG(IJ,IG).EQ.NL-1) THEN
              IJFROMIK(IXLG(IJ,IG)-NL,KXLT(IJ,IG),IG)=IJ
            ENDIF
          ENDDO
          DO IH=1,NLENHALO(IRANK)
            IJ=IJFROMPE(IH,IRANK)
            NL=NLONRGG(KXLT(IJ,IG))
            IF(IXLG(IJ,IG).EQ.1) THEN
              IJFROMIK(IXLG(IJ,IG)+NL,KXLT(IJ,IG),IG)=IJHALO(IH)
            ELSEIF(IXLG(IJ,IG).EQ.2) THEN
              IJFROMIK(IXLG(IJ,IG)+NL,KXLT(IJ,IG),IG)=IJHALO(IH)
            ELSE IF(IXLG(IJ,IG).EQ.NL) THEN
              IJFROMIK(IXLG(IJ,IG)-NL,KXLT(IJ,IG),IG)=IJHALO(IH)
            ELSE IF(IXLG(IJ,IG).EQ.NL-1) THEN
              IJFROMIK(IXLG(IJ,IG)-NL,KXLT(IJ,IG),IG)=IJHALO(IH)
            ENDIF
          ENDDO
        ENDIF
      ELSE
        DO IJ=NINF,NSUP
          IJFROMIK(IXLG(IJ,IG),KXLT(IJ,IG),IG)=IJ
          IFROMIJ(IJ)=IXLG(IJ,IG)
          KFROMIJ(IJ)=KXLT(IJ,IG)
        ENDDO
!       PERIODICITY
        IF(IPER.EQ.1) THEN
          DO IJ=NINF,NSUP
            NL=NLONRGG(KXLT(IJ,IG))
            IF(IXLG(IJ,IG).EQ.1) THEN
              IJFROMIK(IXLG(IJ,IG)+NL,KXLT(IJ,IG),IG)=IJ
            ELSEIF(IXLG(IJ,IG).EQ.2) THEN
              IJFROMIK(IXLG(IJ,IG)+NL,KXLT(IJ,IG),IG)=IJ
            ELSE IF(IXLG(IJ,IG).EQ.NL) THEN
              IJFROMIK(IXLG(IJ,IG)-NL,KXLT(IJ,IG),IG)=IJ
            ELSE IF(IXLG(IJ,IG).EQ.NL-1) THEN
              IJFROMIK(IXLG(IJ,IG)-NL,KXLT(IJ,IG),IG)=IJ
            ENDIF
          ENDDO
        ENDIF
      ENDIF

      DEALLOCATE(IJFROMPE)
      DEALLOCATE(IPROCFROM)
      DEALLOCATE(NLENHALO)
      DEALLOCATE(IJHALO)


!     FIND BOUNDARY POINTS FOR EACH AREA
!     !!! it has to run before 5. (see below) 

      IF (ALLOCATED(LLBOUND)) DEALLOCATE(LLBOUND)
      ALLOCATE(LLBOUND(NIBLO))
      DO IJ=1,NIBLO
        LLBOUND(IJ)=.FALSE.
      ENDDO

      DO IP = 1, NPR

        DO IC=1,2
          DO IJ=NSTART(IP),NEND(IP)
           IF(  KLON(IJ,IC).GT.0 .AND. KLON(IJ,IC).LE.NIBLO .AND.
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

        IF(IPROPAGS.EQ.2) THEN
          DO ICL=-1,1
            DO IC=1,4
              DO IJ=NSTART(IP),NEND(IP)
               IF(  KCOR(IJ,IC,ICL).GT.0 .AND.
     &              KCOR(IJ,IC,ICL).LE.NIBLO .AND.
     &             (KCOR(IJ,IC,ICL).LT.NSTART(IP).OR.
     &              KCOR(IJ,IC,ICL).GT.NEND(IP) ) ) LLBOUND(IJ)=.TRUE.
              ENDDO
            ENDDO
          ENDDO
        ELSE IF(IPROPAGS.EQ.1) THEN
          DO IC=1,2
            DO IJ=NSTART(IP),NEND(IP)
             IF(  KRLON(IJ,IC).GT.0 .AND. KRLON(IJ,IC).LE.NIBLO .AND.
     &           (KRLON(IJ,IC).LT.NSTART(IP).OR.
     &            KRLON(IJ,IC).GT.NEND(IP)       ) ) LLBOUND(IJ)=.TRUE.
             IF(  KRLAT(IJ,IC).GT.0 .AND. KRLAT(IJ,IC).LE.NIBLO .AND.
     &           (KRLAT(IJ,IC).LT.NSTART(IP).OR.
     &            KRLAT(IJ,IC).GT.NEND(IP)       ) ) LLBOUND(IJ)=.TRUE.
             IF(  KRLON2(IJ,IC).GT.0 .AND. KRLON2(IJ,IC).LE.NIBLO .AND.
     &           (KRLON2(IJ,IC).LT.NSTART(IP).OR.
     &            KRLON2(IJ,IC).GT.NEND(IP)      ) ) LLBOUND(IJ)=.TRUE.
             IF(  KRLAT2(IJ,IC).GT.0 .AND. KRLAT2(IJ,IC).LE.NIBLO .AND.
     &           (KRLAT2(IJ,IC).LT.NSTART(IP).OR.
     &            KRLAT2(IJ,IC).GT.NEND(IP)      ) ) LLBOUND(IJ)=.TRUE.
            ENDDO
          ENDDO
        ENDIF

      ENDDO


!     4. KEEP THE PART OF KLAT,KLAT2,KLON,KCOR, DEPTH WHICH IS NECESSARY
!        ----------------------------------------------------------


      IF(LMESSPASS.AND.NPR.GT.1) THEN

        ALLOCATE(KDUM2(NSTART(IRANK):NEND(IRANK),2))

        DO IC=1,2
          DO IJ=NSTART(IRANK),NEND(IRANK)
            KDUM2(IJ,IC)=KLAT(IJ,IC)
          ENDDO
        ENDDO

        DEALLOCATE(KLAT)
        ALLOCATE(KLAT(NSTART(IRANK):NEND(IRANK),2))

        DO IC=1,2
          DO IJ=NSTART(IRANK),NEND(IRANK)
            KLAT(IJ,IC)=KDUM2(IJ,IC)
          ENDDO
        ENDDO

        DO IC=1,2
          DO IJ=NSTART(IRANK),NEND(IRANK)
            KDUM2(IJ,IC)=KLON(IJ,IC)
          ENDDO
        ENDDO

        DEALLOCATE(KLON)
        ALLOCATE(KLON(NSTART(IRANK):NEND(IRANK),2))

        DO IC=1,2
          DO IJ=NSTART(IRANK),NEND(IRANK)
            KLON(IJ,IC)=KDUM2(IJ,IC)
          ENDDO
        ENDDO

        DO IC=1,2
          DO IJ=NSTART(IRANK),NEND(IRANK)
            KDUM2(IJ,IC)=KLAT2(IJ,IC)
          ENDDO
        ENDDO

        DEALLOCATE(KLAT2)
        ALLOCATE(KLAT2(NSTART(IRANK):NEND(IRANK),2))

        DO IC=1,2
          DO IJ=NSTART(IRANK),NEND(IRANK)
            KLAT2(IJ,IC)=KDUM2(IJ,IC)
          ENDDO
        ENDDO

        DEALLOCATE(KDUM2)

        IF(IPROPAGS.EQ.1) THEN
          ALLOCATE(KDUM2(NSTART(IRANK):NEND(IRANK),2))
          DO IC=1,2
            DO IJ=NSTART(IRANK),NEND(IRANK)
              KDUM2(IJ,IC)=KRLAT(IJ,IC)
            ENDDO
          ENDDO

          DEALLOCATE(KRLAT)
          ALLOCATE(KRLAT(NSTART(IRANK):NEND(IRANK),2))

          DO IC=1,2
            DO IJ=NSTART(IRANK),NEND(IRANK)
              KRLAT(IJ,IC)=KDUM2(IJ,IC)
            ENDDO
          ENDDO

          DO IC=1,2
            DO IJ=NSTART(IRANK),NEND(IRANK)
              KDUM2(IJ,IC)=KRLON(IJ,IC)
            ENDDO
          ENDDO

          DEALLOCATE(KRLON)
          ALLOCATE(KRLON(NSTART(IRANK):NEND(IRANK),2))

          DO IC=1,2
            DO IJ=NSTART(IRANK),NEND(IRANK)
              KRLON(IJ,IC)=KDUM2(IJ,IC)
            ENDDO
          ENDDO

          DO IC=1,2
            DO IJ=NSTART(IRANK),NEND(IRANK)
              KDUM2(IJ,IC)=KRLAT2(IJ,IC)
           ENDDO
          ENDDO

          DEALLOCATE(KRLAT2)
          ALLOCATE(KRLAT2(NSTART(IRANK):NEND(IRANK),2))

          DO IC=1,2
            DO IJ=NSTART(IRANK),NEND(IRANK)
              KRLAT2(IJ,IC)=KDUM2(IJ,IC)
            ENDDO
          ENDDO

          DO IC=1,2
            DO IJ=NSTART(IRANK),NEND(IRANK)
              KDUM2(IJ,IC)=KRLON2(IJ,IC)
            ENDDO
          ENDDO

          DEALLOCATE(KRLON2)
          ALLOCATE(KRLON2(NSTART(IRANK):NEND(IRANK),2))

          DO IC=1,2
            DO IJ=NSTART(IRANK),NEND(IRANK)
              KRLON2(IJ,IC)=KDUM2(IJ,IC)
            ENDDO
          ENDDO
          DEALLOCATE(KDUM2)
        ELSE IF(IPROPAGS.EQ.2) THEN
          ALLOCATE(KDUM3(NSTART(IRANK):NEND(IRANK),4,-1:1))
          DO ICL=-1,1
            DO IC=1,4
              DO IJ=NSTART(IRANK),NEND(IRANK)
                KDUM3(IJ,IC,ICL)=KCOR(IJ,IC,ICL)
              ENDDO
            ENDDO
          ENDDO

          DEALLOCATE(KCOR)
          ALLOCATE(KCOR(NSTART(IRANK):NEND(IRANK),4,-1:1))

          DO ICL=-1,1
            DO IC=1,4
              DO IJ=NSTART(IRANK),NEND(IRANK)
                KCOR(IJ,IC,ICL)=KDUM3(IJ,IC,ICL)
              ENDDO
            ENDDO
          ENDDO
          DEALLOCATE(KDUM3)
        ENDIF

        ALLOCATE(RDUM2(NINF:NSUP,IGL))
        DO IG=1,IGL
          DO IJ=NINF,NSUP
            RDUM2(IJ,IG)=DEPTH(IJ,IG)
          ENDDO
        ENDDO
        DEALLOCATE(DEPTH)
        ALLOCATE(DEPTH(NINF:NSUP,IGL))
        DO IG=1,IGL
          DO IJ=NINF,NSUP
            DEPTH(IJ,IG)=RDUM2(IJ,IG)
          ENDDO
        ENDDO
        DEALLOCATE(RDUM2)

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
        IF(IPROPAGS.EQ.2) THEN
          DO ICL=-1,1
            DO IC=1,4
              DO IJ = NSTART(IRANK),NEND(IRANK)
                IF(KCOR(IJ,IC,ICL).EQ.0) KCOR(IJ,IC,ICL) = NINF-1
              ENDDO
            ENDDO
          ENDDO
        ELSE IF(IPROPAGS.EQ.1) THEN
          DO IC=1,2
            DO IJ = NSTART(IRANK),NEND(IRANK)
              IF(KRLAT(IJ,IC).EQ.0) KRLAT(IJ,IC) = NINF-1
              IF(KRLAT2(IJ,IC).EQ.0) KRLAT2(IJ,IC) = NINF-1
              IF(KRLON(IJ,IC).EQ.0) KRLON(IJ,IC) = NINF-1
              IF(KRLON2(IJ,IC).EQ.0) KRLON2(IJ,IC) = NINF-1
            ENDDO
          ENDDO
        ENDIF
      ENDIF



!*    6. INPUT THE WEIGHT FOR THE ADVECTION SCHEME 
!        --------------------------------------------

      IF (ALLOCATED(WLAT)) DEALLOCATE(WLAT)
      IF (ALLOCATED(WRLAT)) DEALLOCATE(WRLAT)
      IF (ALLOCATED(WRLON)) DEALLOCATE(WRLON)

      IF(LMESSPASS) THEN

!       READ SECOND PART OF IU08 ON PE IREAD

        IF(IRANK.EQ.IREAD) THEN
!         THEIR SIZE IS READJUSTED (SEE BELOW)
          ALLOCATE(WLAT(NIBLO,2))
          ALLOCATE(WRLAT(NIBLO,2))
          ALLOCATE(WRLON(NIBLO,2))
          CALL GSTATS(1771,0)
          READ (IU08) WLAT
          READ (IU08) WRLAT, WRLON
          CALL GSTATS(1771,1)

!         RELABELLING OF THE ARRAYS
          IF(.NOT.LL1D.AND.NPR.GT.1) THEN

            ALLOCATE(RDUM(NIBLO))
            DO IC=1,2
              DO NIJ=NSTART(1),NEND(NPR)
                RDUM(NIJ)=WLAT(NEWIJ2IJ(NIJ),IC)
              ENDDO
              DO NIJ=NSTART(1),NEND(NPR)
               WLAT(NIJ,IC)=RDUM(NIJ)
              ENDDO
            ENDDO

            DEALLOCATE(RDUM)

            IF(IPROPAGS.EQ.1) THEN
              ALLOCATE(RDUM(NIBLO))
              DO IC=1,2
                DO NIJ=NSTART(1),NEND(NPR)
                  RDUM(NIJ)=WRLAT(NEWIJ2IJ(NIJ),IC)
                ENDDO
                DO NIJ=NSTART(1),NEND(NPR)
                 WRLAT(NIJ,IC)=RDUM(NIJ)
                ENDDO
              ENDDO
              DO IC=1,2
                DO NIJ=NSTART(1),NEND(NPR)
                  RDUM(NIJ)=WRLON(NEWIJ2IJ(NIJ),IC)
                ENDDO
                DO NIJ=NSTART(1),NEND(NPR)
                 WRLON(NIJ,IC)=RDUM(NIJ)
                ENDDO
              ENDDO
              DEALLOCATE(RDUM)
            ENDIF

          ENDIF

        ENDIF

!       SEND WLAT WRLAT WRLON TO OTHER PE'S
        IF(NPR.GT.1) THEN
          ITAG=KTAG+1
          IF(IPROPAGS.EQ.1) THEN
            MPLENGTH=6*MAXLEN
          ELSE
            MPLENGTH=2*MAXLEN
          ENDIF

          IF(IRANK.EQ.IREAD) THEN
            ALLOCATE(RCOMBUF_S(MPLENGTH,NPR))
          ELSE
            ALLOCATE(RCOMBUF_R(MPLENGTH))
          ENDIF

          IF(IRANK.EQ.IREAD) THEN
!           SEND TO OTHER PE'S

!           FILL THE SEND BUFFER
            DO IP=1,NPR
              KCOUNT=0
              DO IC=1,2
                DO IJ=NSTART(IP),NEND(IP)
                  KCOUNT=KCOUNT+1
                  RCOMBUF_S(KCOUNT,IP)=WLAT(IJ,IC)
                ENDDO
              ENDDO

              IF(IPROPAGS.EQ.1) THEN
                DO IC=1,2
                  DO IJ=NSTART(IP),NEND(IP)
                    KCOUNT=KCOUNT+1
                    RCOMBUF_S(KCOUNT,IP)=WRLAT(IJ,IC)
                  ENDDO
                ENDDO
                DO IC=1,2
                  DO IJ=NSTART(IP),NEND(IP)
                    KCOUNT=KCOUNT+1
                    RCOMBUF_S(KCOUNT,IP)=WRLON(IJ,IC)
                  ENDDO
                ENDDO
              ENDIF
            ENDDO

            DEALLOCATE(WLAT)
            DEALLOCATE(WRLAT)
            DEALLOCATE(WRLON)

!           SEND NON BLOCKING THE BUFFERS
            CALL GSTATS(694,0)
            IPCOUNT=0
            DO IP=1,NPR
              IF(IP.NE.IREAD) THEN
                IPCOUNT=IPCOUNT+1
                CALL MPL_SEND(RCOMBUF_S(1:MPLENGTH,IP),
     &           KDEST=IP,KTAG=ITAG,
     &           KMP_TYPE=JP_NON_BLOCKING_STANDARD,
     &           KREQUEST=ISENDREQ(IPCOUNT),
     &           CDSTRING='MPDECOMP 1:')
              ENDIF
            ENDDO
            CALL GSTATS(694,1)

          ELSE
!           RECEIVE FROM PE IREAD 
            CALL GSTATS(694,0)
            CALL MPL_RECV(RCOMBUF_R(1:MPLENGTH),
     &       KSOURCE=IREAD,KTAG=ITAG,
     &       KOUNT=KRCOUNT,KRECVTAG=KRTAG,
     &       KMP_TYPE=JP_BLOCKING_STANDARD,
     &       CDSTRING='MPDECOMP 1:')
            IF(KRCOUNT.NE.MPLENGTH) CALL MPL_ABORT
     &       ('MPL_RECV ERROR in MPDECOMP 1:MISMATCHED MSG LENGTH')
            IF(KRTAG.NE.ITAG) CALL MPL_ABORT
     &      ('MPL_RECV ERROR in MPDECOMP 1 MISMATCHED TAGS' )
            CALL GSTATS(694,1)
          ENDIF

!         KEEP THE RELEVANT PART OF WLAT
          ALLOCATE(WLAT(NSTART(IRANK):NEND(IRANK),2))
          IF(IPROPAGS.EQ.1) THEN
            ALLOCATE(WRLAT(NSTART(IRANK):NEND(IRANK),2))
            ALLOCATE(WRLON(NSTART(IRANK):NEND(IRANK),2))
          ENDIF

          IF(IRANK.EQ.IREAD) THEN
            KCOUNT=0
            DO IC=1,2
              DO IJ=NSTART(IRANK),NEND(IRANK)
                KCOUNT=KCOUNT+1
                WLAT(IJ,IC)=RCOMBUF_S(KCOUNT,IRANK)
              ENDDO
            ENDDO

            IF(IPROPAGS.EQ.1) THEN
              DO IC=1,2
                DO IJ=NSTART(IRANK),NEND(IRANK)
                  KCOUNT=KCOUNT+1
                  WRLAT(IJ,IC)=RCOMBUF_S(KCOUNT,IRANK)
                ENDDO
              ENDDO
              DO IC=1,2
                DO IJ=NSTART(IRANK),NEND(IRANK)
                  KCOUNT=KCOUNT+1
                  WRLON(IJ,IC)=RCOMBUF_S(KCOUNT,IRANK)
                ENDDO
              ENDDO
            ENDIF

          ELSE
            KCOUNT=0
            DO IC=1,2
              DO IJ=NSTART(IRANK),NEND(IRANK)
                KCOUNT=KCOUNT+1
                WLAT(IJ,IC)=RCOMBUF_R(KCOUNT)
              ENDDO
            ENDDO
            IF(IPROPAGS.EQ.1) THEN
              DO IC=1,2
                DO IJ=NSTART(IRANK),NEND(IRANK)
                  KCOUNT=KCOUNT+1
                  WRLAT(IJ,IC)=RCOMBUF_R(KCOUNT)
                ENDDO
              ENDDO
              DO IC=1,2
                DO IJ=NSTART(IRANK),NEND(IRANK)
                  KCOUNT=KCOUNT+1
                  WRLON(IJ,IC)=RCOMBUF_R(KCOUNT)
                ENDDO
              ENDDO
            ENDIF
          ENDIF

!         WAIT ANY OUTSTANDING SENDS TO COMPLETE
!         --------------------------------------
          IF(IRANK.EQ.IREAD) THEN
            CALL GSTATS(694,0)
            CALL MPL_WAIT(RCOMBUF_S(1:MPLENGTH,1),KREQUEST=ISENDREQ,
     &              CDSTRING='WAIT IN MPDECOMP 1:')
            CALL GSTATS(694,1)
          ENDIF

          IF(ALLOCATED(RCOMBUF_S)) DEALLOCATE(RCOMBUF_S)
          IF(ALLOCATED(RCOMBUF_R)) DEALLOCATE(RCOMBUF_R)

        ENDIF


!       READ THE REST OF IU08
!       =====================

        IF(.NOT.LL1D.AND.NPR.GT.1) ALLOCATE(KDUM(NIBLO))

        IF(NPR.GT.1) THEN
          MPLENGTH=2*NFRE*MAXLEN
          IF(IRANK.EQ.IREAD) THEN
            ALLOCATE(KCOMBUF_S(MPLENGTH,NPR))
          ELSE
            ALLOCATE(KCOMBUF_R(MPLENGTH))
          ENDIF
        ENDIF

!       READING KOBSLAT
!       ---------------
        ITAG=KTAG+2
        IF(IRANK.EQ.IREAD) THEN
          CALL GSTATS(1771,0)
          ALLOCATE(KOBSLAT(NIBLO,NFRE,2))
          READ (IU08) KOBSLAT
          CALL GSTATS(1771,1)

!         RELABELLING OF THE ARRAY
          IF(.NOT.LL1D.AND.NPR.GT.1) THEN
            DO IC=1,2
              DO M=1,NFRE
                DO NIJ=NSTART(1),NEND(NPR)
                   KDUM(NIJ)=KOBSLAT(NEWIJ2IJ(NIJ),M,IC)
                ENDDO
                DO NIJ=NSTART(1),NEND(NPR)
                  KOBSLAT(NIJ,M,IC)=KDUM(NIJ)
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!         SEND TO OTHER PE'S
          IF(NPR.GT.1) THEN
!           FILL THE SEND BUFFER
            DO IP=1,NPR
              KCOUNT=0
              DO IC=1,2
                DO M=1,NFRE
                  DO IJ=NSTART(IP),NEND(IP)
                    KCOUNT=KCOUNT+1
                    KCOMBUF_S(KCOUNT,IP)=KOBSLAT(IJ,M,IC)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

            DEALLOCATE(KOBSLAT)

!           SEND NON BLOCKING THE BUFFERS
            CALL GSTATS(694,0)
            IPCOUNT=0
            DO IP=1,NPR
              IF(IP.NE.IREAD) THEN
                IPCOUNT=IPCOUNT+1
                CALL MPL_SEND(KCOMBUF_S(1:MPLENGTH,IP),
     &           KDEST=IP,KTAG=ITAG,
     &           KMP_TYPE=JP_NON_BLOCKING_STANDARD,
     &           KREQUEST=ISENDREQ(IPCOUNT),
     &           CDSTRING='MPDECOMP SENDING KOBSLAT:')
              ENDIF
            ENDDO
            CALL GSTATS(694,1)

          ENDIF

        ELSE
!         RECEIVE FROM PE IREAD 
          IF(NPR.GT.1) THEN
            CALL GSTATS(694,0)
            CALL MPL_RECV(KCOMBUF_R(1:MPLENGTH),
     &       KSOURCE=IREAD,KTAG=ITAG,
     &       KOUNT=KRCOUNT,KRECVTAG=KRTAG,
     &       KMP_TYPE=JP_BLOCKING_STANDARD,
     &       CDSTRING='MPDECOMP RECEIVING KOBSLAT:')
            IF(KRCOUNT.NE.MPLENGTH) CALL MPL_ABORT
     &       ('MPL_RECV ERROR in MPDECOMP 2:MISMATCHED MSG LENGTH')
            IF(KRTAG.NE.ITAG) CALL MPL_ABORT
     &      ('MPL_RECV ERROR in MPDECOMP 2 MISMATCHED TAGS' )
            CALL GSTATS(694,1)

          ENDIF
        ENDIF

        IF(NPR.GT.1) THEN
!       KEEP THE RELEVANT PART OF KOBSLAT
          ALLOCATE(KOBSLAT(NSTART(IRANK):NEND(IRANK),NFRE,2))
          IF(IRANK.EQ.IREAD) THEN
            KCOUNT=0
            DO IC=1,2
              DO M=1,NFRE
                DO IJ=NSTART(IRANK),NEND(IRANK)
                  KCOUNT=KCOUNT+1
                  KOBSLAT(IJ,M,IC)=KCOMBUF_S(KCOUNT,IRANK)
                ENDDO
              ENDDO
            ENDDO
          ELSE
            KCOUNT=0
            DO IC=1,2
              DO M=1,NFRE
                DO IJ=NSTART(IRANK),NEND(IRANK)
                  KCOUNT=KCOUNT+1
                  KOBSLAT(IJ,M,IC)=KCOMBUF_R(KCOUNT)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!       WAIT ANY OUTSTANDING SENDS TO COMPLETE
!       --------------------------------------
        IF(NPR.GT.1) THEN
          IF(IRANK.EQ.IREAD) THEN
            CALL GSTATS(694,0)
            CALL MPL_WAIT(KCOMBUF_S(1:MPLENGTH,1),KREQUEST=ISENDREQ,
     &              CDSTRING='WAIT IN MPDECOMP 2:')
            CALL GSTATS(694,1)
          ENDIF
        ENDIF


!       READING KOBSLON
!       ---------------
        ITAG=KTAG+3
        IF(IRANK.EQ.IREAD) THEN
          CALL GSTATS(1771,0)
          ALLOCATE(KOBSLON(NIBLO,NFRE,2))
          READ (IU08) KOBSLON
          CALL GSTATS(1771,1)

!         RELABELLING OF THE ARRAY
          IF(.NOT.LL1D.AND.NPR.GT.1) THEN
            DO IC=1,2
              DO M=1,NFRE
                DO NIJ=NSTART(1),NEND(NPR)
                   KDUM(NIJ)=KOBSLON(NEWIJ2IJ(NIJ),M,IC)
                ENDDO
                DO NIJ=NSTART(1),NEND(NPR)
                  KOBSLON(NIJ,M,IC)=KDUM(NIJ)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
!         SEND TO OTHER PE'S
          IF(NPR.GT.1) THEN
!           FILL THE SEND BUFFER
            DO IP=1,NPR
              KCOUNT=0
              DO IC=1,2
                DO M=1,NFRE
                  DO IJ=NSTART(IP),NEND(IP)
                    KCOUNT=KCOUNT+1
                    KCOMBUF_S(KCOUNT,IP)=KOBSLON(IJ,M,IC)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

            DEALLOCATE(KOBSLON)

!           SEND NON BLOCKING THE BUFFERS
            CALL GSTATS(694,0)
            IPCOUNT=0
            DO IP=1,NPR
              IF(IP.NE.IREAD) THEN
                IPCOUNT=IPCOUNT+1
                CALL MPL_SEND(KCOMBUF_S(1:MPLENGTH,IP),
     &           KDEST=IP,KTAG=ITAG,
     &           KMP_TYPE=JP_NON_BLOCKING_STANDARD,
     &           KREQUEST=ISENDREQ(IPCOUNT),
     &           CDSTRING='MPDECOMP SENDING KOBSLON:')
              ENDIF
            ENDDO
            CALL GSTATS(694,1)

          ENDIF

        ELSE
!         RECEIVE FROM PE IREAD 
          IF(NPR.GT.1) THEN
            CALL GSTATS(694,0)
            CALL MPL_RECV(KCOMBUF_R(1:MPLENGTH),
     &       KSOURCE=IREAD,KTAG=ITAG,
     &       KOUNT=KRCOUNT,KRECVTAG=KRTAG,
     &       KMP_TYPE=JP_BLOCKING_STANDARD,
     &       CDSTRING='MPDECOMP RECEIVING KOBSLON:')
            IF(KRCOUNT.NE.MPLENGTH) CALL MPL_ABORT
     &       ('MPL_RECV ERROR in MPDECOMP 3:MISMATCHED MSG LENGTH')
            IF(KRTAG.NE.ITAG) CALL MPL_ABORT
     &      ('MPL_RECV ERROR in MPDECOMP 3 MISMATCHED TAGS' )
            CALL GSTATS(694,1)

          ENDIF
        ENDIF

        IF(NPR.GT.1) THEN
!       KEEP THE RELEVANT PART OF KOBSLON
          ALLOCATE(KOBSLON(NSTART(IRANK):NEND(IRANK),NFRE,2))
          IF(IRANK.EQ.IREAD) THEN
            KCOUNT=0
            DO IC=1,2
              DO M=1,NFRE
                DO IJ=NSTART(IRANK),NEND(IRANK)
                  KCOUNT=KCOUNT+1
                  KOBSLON(IJ,M,IC)=KCOMBUF_S(KCOUNT,IRANK)
                ENDDO
              ENDDO
            ENDDO
          ELSE
            KCOUNT=0
            DO IC=1,2
              DO M=1,NFRE
                DO IJ=NSTART(IRANK),NEND(IRANK)
                  KCOUNT=KCOUNT+1
                  KOBSLON(IJ,M,IC)=KCOMBUF_R(KCOUNT)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!       WAIT ANY OUTSTANDING SENDS TO COMPLETE
!       --------------------------------------
        IF(NPR.GT.1) THEN
          IF(IRANK.EQ.IREAD) THEN
            CALL GSTATS(694,0)
            CALL MPL_WAIT(KCOMBUF_S(1:MPLENGTH,1),KREQUEST=ISENDREQ,
     &              CDSTRING='WAIT IN MPDECOMP 3:')
            CALL GSTATS(694,1)
          ENDIF
        ENDIF


!       READING KOBSRLAT
!       ----------------
        IF(IRANK.EQ.IREAD) THEN
          CALL GSTATS(1771,0)
          ALLOCATE(KOBSRLAT(NIBLO,NFRE,2))
          READ (IU08) KOBSRLAT
          CALL GSTATS(1771,1)
        ENDIF

!       DISTRIBUTE TO OTHER PE'S (IF NEEDED).
        IF(IPROPAGS.EQ.1) THEN
          ITAG=KTAG+2
          IF(IRANK.EQ.IREAD) THEN
!           RELABELLING OF THE ARRAY
            IF(.NOT.LL1D.AND.NPR.GT.1) THEN
              DO IC=1,2
                DO M=1,NFRE
                  DO NIJ=NSTART(1),NEND(NPR)
                     KDUM(NIJ)=KOBSRLAT(NEWIJ2IJ(NIJ),M,IC)
                  ENDDO
                  DO NIJ=NSTART(1),NEND(NPR)
                    KOBSRLAT(NIJ,M,IC)=KDUM(NIJ)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
!           SEND TO OTHER PE'S
            IF(NPR.GT.1) THEN
!             FILL THE SEND BUFFER
              DO IP=1,NPR
                KCOUNT=0
                DO IC=1,2
                  DO M=1,NFRE
                    DO IJ=NSTART(IP),NEND(IP)
                      KCOUNT=KCOUNT+1
                      KCOMBUF_S(KCOUNT,IP)=KOBSRLAT(IJ,M,IC)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO

              DEALLOCATE(KOBSRLAT)

!             SEND NON BLOCKING THE BUFFERS
              CALL GSTATS(694,0)
              IPCOUNT=0
              DO IP=1,NPR
                IF(IP.NE.IREAD) THEN
                  IPCOUNT=IPCOUNT+1
                  CALL MPL_SEND(KCOMBUF_S(1:MPLENGTH,IP),
     &             KDEST=IP,KTAG=ITAG,
     &             KMP_TYPE=JP_NON_BLOCKING_STANDARD,
     &             KREQUEST=ISENDREQ(IPCOUNT),
     &             CDSTRING='MPDECOMP SENDING KOBSRLAT:')
                ENDIF
              ENDDO
              CALL GSTATS(694,1)

            ENDIF

          ELSE
!           RECEIVE FROM PE IREAD 
            IF(NPR.GT.1) THEN
              CALL GSTATS(694,0)
              CALL MPL_RECV(KCOMBUF_R(1:MPLENGTH),
     &         KSOURCE=IREAD,KTAG=ITAG,
     &         KOUNT=KRCOUNT,KRECVTAG=KRTAG,
     &         KMP_TYPE=JP_BLOCKING_STANDARD,
     &         CDSTRING='MPDECOMP RECEIVING KOBSRLAT:')
              IF(KRCOUNT.NE.MPLENGTH) CALL MPL_ABORT
     &         ('MPL_RECV ERROR in MPDECOMP 4:MISMATCHED MSG LENGTH')
              IF(KRTAG.NE.ITAG) CALL MPL_ABORT
     &        ('MPL_RECV ERROR in MPDECOMP 4 MISMATCHED TAGS' )
              CALL GSTATS(694,1)

            ENDIF
          ENDIF

          IF(NPR.GT.1) THEN
!         KEEP THE RELEVANT PART OF KOBSRLAT
            ALLOCATE(KOBSRLAT(NSTART(IRANK):NEND(IRANK),NFRE,2))
            IF(IRANK.EQ.IREAD) THEN
              KCOUNT=0
              DO IC=1,2
                DO M=1,NFRE
                  DO IJ=NSTART(IRANK),NEND(IRANK)
                    KCOUNT=KCOUNT+1
                    KOBSRLAT(IJ,M,IC)=KCOMBUF_S(KCOUNT,IRANK)
                  ENDDO
                ENDDO
              ENDDO
           ELSE
              KCOUNT=0
              DO IC=1,2
                DO M=1,NFRE
                  DO IJ=NSTART(IRANK),NEND(IRANK)
                    KCOUNT=KCOUNT+1
                    KOBSRLAT(IJ,M,IC)=KCOMBUF_R(KCOUNT)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF

!         WAIT ANY OUTSTANDING SENDS TO COMPLETE
!         --------------------------------------
          IF(NPR.GT.1) THEN
            IF(IRANK.EQ.IREAD) THEN
              CALL GSTATS(694,0)
              CALL MPL_WAIT(KCOMBUF_S(1:MPLENGTH,1),KREQUEST=ISENDREQ,
     &                CDSTRING='WAIT IN MPDECOMP 4:')
              CALL GSTATS(694,1)
            ENDIF
          ENDIF

        ELSE
          IF(ALLOCATED(KOBSRLAT)) DEALLOCATE(KOBSRLAT)
        ENDIF


!       READING KOBSRLON
!       ----------------
        IF(IRANK.EQ.IREAD) THEN
          CALL GSTATS(1771,0)
          ALLOCATE(KOBSRLON(NIBLO,NFRE,2))
          READ (IU08) KOBSRLON
          CALL GSTATS(1771,1)
        ENDIF

!       DISTRIBUTE TO OTHER PE'S (IF NEEDED).
        IF(IPROPAGS.EQ.1) THEN
          ITAG=KTAG+5
          IF(IRANK.EQ.IREAD) THEN
!           RELABELLING OF THE ARRAY
            IF(.NOT.LL1D.AND.NPR.GT.1) THEN
              DO IC=1,2
                DO M=1,NFRE
                  DO NIJ=NSTART(1),NEND(NPR)
                     KDUM(NIJ)=KOBSRLON(NEWIJ2IJ(NIJ),M,IC)
                  ENDDO
                  DO NIJ=NSTART(1),NEND(NPR)
                    KOBSRLON(NIJ,M,IC)=KDUM(NIJ)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
!           SEND TO OTHER PE'S
            IF(NPR.GT.1) THEN
!             FILL THE SEND BUFFER
              DO IP=1,NPR
                KCOUNT=0
                DO IC=1,2
                  DO M=1,NFRE
                    DO IJ=NSTART(IP),NEND(IP)
                      KCOUNT=KCOUNT+1
                      KCOMBUF_S(KCOUNT,IP)=KOBSRLON(IJ,M,IC)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO

              DEALLOCATE(KOBSRLON)

!             SEND NON BLOCKING THE BUFFERS
              CALL GSTATS(694,0)
              IPCOUNT=0
              DO IP=1,NPR
                IF(IP.NE.IREAD) THEN
                  IPCOUNT=IPCOUNT+1
                  CALL MPL_SEND(KCOMBUF_S(1:MPLENGTH,IP),
     &            KDEST=IP,KTAG=ITAG,
     &             KMP_TYPE=JP_NON_BLOCKING_STANDARD,
     &             KREQUEST=ISENDREQ(IPCOUNT),
     &             CDSTRING='MPDECOMP SENDING KOBSRLON:')
                ENDIF
              ENDDO
              CALL GSTATS(694,1)

            ENDIF

          ELSE
!           RECEIVE FROM PE IREAD 
            IF(NPR.GT.1) THEN
              CALL GSTATS(694,0)
              CALL MPL_RECV(KCOMBUF_R(1:MPLENGTH),
     &         KSOURCE=IREAD,KTAG=ITAG,
     &         KOUNT=KRCOUNT,KRECVTAG=KRTAG,
     &         KMP_TYPE=JP_BLOCKING_STANDARD,
     &         CDSTRING='MPDECOMP RECEIVING KOBSRLON:')
              IF(KRCOUNT.NE.MPLENGTH) CALL MPL_ABORT
     &         ('MPL_RECV ERROR in MPDECOMP 5:MISMATCHED MSG LENGTH')
              IF(KRTAG.NE.ITAG) CALL MPL_ABORT
     &        ('MPL_RECV ERROR in MPDECOMP 5 MISMATCHED TAGS' )
              CALL GSTATS(694,1)

            ENDIF
          ENDIF

          IF(NPR.GT.1) THEN
!         KEEP THE RELEVANT PART OF KOBSRLON
            ALLOCATE(KOBSRLON(NSTART(IRANK):NEND(IRANK),NFRE,2))
            IF(IRANK.EQ.IREAD) THEN
              KCOUNT=0
              DO IC=1,2
                DO M=1,NFRE
                  DO IJ=NSTART(IRANK),NEND(IRANK)
                    KCOUNT=KCOUNT+1
                    KOBSRLON(IJ,M,IC)=KCOMBUF_S(KCOUNT,IRANK)
                  ENDDO
                ENDDO
              ENDDO
            ELSE
              KCOUNT=0
              DO IC=1,2
                DO M=1,NFRE
                  DO IJ=NSTART(IRANK),NEND(IRANK)
                    KCOUNT=KCOUNT+1
                    KOBSRLON(IJ,M,IC)=KCOMBUF_R(KCOUNT)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF

!         WAIT ANY OUTSTANDING SENDS TO COMPLETE
!         --------------------------------------
          IF(NPR.GT.1) THEN
            IF(IRANK.EQ.IREAD) THEN
              CALL GSTATS(694,0)
              CALL MPL_WAIT(KCOMBUF_S(1:MPLENGTH,1),KREQUEST=ISENDREQ,
     &                CDSTRING='WAIT IN MPDECOMP 5:')
              CALL GSTATS(694,1)
            ENDIF
          ENDIF

        ELSE
          IF(ALLOCATED(KOBSRLON)) DEALLOCATE(KOBSRLON)
        ENDIF


!       READING KOBSCOR
!       ---------------


        IF(IRANK.EQ.IREAD) THEN
          CALL GSTATS(1771,0)
          ALLOCATE(KOBSCOR(NIBLO,NFRE,4))
          READ (IU08) KOBSCOR
          CALL GSTATS(1771,1)
        ENDIF

!       DISTRIBUTE TO OTHER PE'S (IF NEEDED).
        ITAG=KTAG+6
        IF(IPROPAGS.EQ.2) THEN
          IF(NPR.GT.1) THEN
            MPLENGTH=4*NFRE*MAXLEN
            IF(IRANK.EQ.IREAD) THEN
              IF(ALLOCATED(KCOMBUF_S)) DEALLOCATE(KCOMBUF_S)
              ALLOCATE(KCOMBUF_S(MPLENGTH,NPR))
            ELSE
              IF(ALLOCATED(KCOMBUF_R)) DEALLOCATE(KCOMBUF_R)
              ALLOCATE(KCOMBUF_R(MPLENGTH))
            ENDIF
          ENDIF
          IF(IRANK.EQ.IREAD) THEN
!           RELABELLING OF THE ARRAY
            IF(.NOT.LL1D.AND.NPR.GT.1) THEN
              DO IC=1,4
                DO M=1,NFRE
                  DO NIJ=NSTART(1),NEND(NPR)
                     KDUM(NIJ)=KOBSCOR(NEWIJ2IJ(NIJ),M,IC)
                  ENDDO
                  DO NIJ=NSTART(1),NEND(NPR)
                    KOBSCOR(NIJ,M,IC)=KDUM(NIJ)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
!           SEND TO OTHER PE'S
            IF(NPR.GT.1) THEN
!             FILL THE SEND BUFFER
              DO IP=1,NPR
                KCOUNT=0
                DO IC=1,4
                  DO M=1,NFRE
                    DO IJ=NSTART(IP),NEND(IP)
                      KCOUNT=KCOUNT+1
                      KCOMBUF_S(KCOUNT,IP)=KOBSCOR(IJ,M,IC)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO

              DEALLOCATE(KOBSCOR)

!             SEND NON BLOCKING THE BUFFERS
              CALL GSTATS(694,0)
              IPCOUNT=0
              DO IP=1,NPR
                IF(IP.NE.IREAD) THEN
                  IPCOUNT=IPCOUNT+1
                  CALL MPL_SEND(KCOMBUF_S(1:MPLENGTH,IP),
     &             KDEST=IP,KTAG=ITAG,
     &             KMP_TYPE=JP_NON_BLOCKING_STANDARD,
     &             KREQUEST=ISENDREQ(IPCOUNT),
     &             CDSTRING='MPDECOMP SENDING KOBSCOR:')
                ENDIF
              ENDDO
              CALL GSTATS(694,1)

            ENDIF

          ELSE
!           RECEIVE FROM PE IREAD 
            IF(NPR.GT.1) THEN
              CALL GSTATS(694,0)
              CALL MPL_RECV(KCOMBUF_R(1:MPLENGTH),
     &         KSOURCE=IREAD,KTAG=ITAG,
     &         KOUNT=KRCOUNT,KRECVTAG=KRTAG,
     &         KMP_TYPE=JP_BLOCKING_STANDARD,
     &         CDSTRING='MPDECOMP RECEIVING KOBSCOR:')
              IF(KRCOUNT.NE.MPLENGTH) CALL MPL_ABORT
     &         ('MPL_RECV ERROR in MPDECOMP 6:MISMATCHED MSG LENGTH')
              IF(KRTAG.NE.ITAG) CALL MPL_ABORT
     &        ('MPL_RECV ERROR in MPDECOMP 6 MISMATCHED TAGS' )
              CALL GSTATS(694,1)

            ENDIF
          ENDIF

          IF(NPR.GT.1) THEN
!         KEEP THE RELEVANT PART OF KOBSCOR
            ALLOCATE(KOBSCOR(NSTART(IRANK):NEND(IRANK),NFRE,4))
            IF(IRANK.EQ.IREAD) THEN
              KCOUNT=0
              DO IC=1,4
                DO M=1,NFRE
                  DO IJ=NSTART(IRANK),NEND(IRANK)
                    KCOUNT=KCOUNT+1
                    KOBSCOR(IJ,M,IC)=KCOMBUF_S(KCOUNT,IRANK)
                  ENDDO
                ENDDO
              ENDDO
            ELSE
              KCOUNT=0
              DO IC=1,4
                DO M=1,NFRE
                  DO IJ=NSTART(IRANK),NEND(IRANK)
                    KCOUNT=KCOUNT+1
                    KOBSCOR(IJ,M,IC)=KCOMBUF_R(KCOUNT)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF

!         WAIT ANY OUTSTANDING SENDS TO COMPLETE
!         --------------------------------------
          IF(NPR.GT.1) THEN
            IF(IRANK.EQ.IREAD) THEN
              CALL GSTATS(694,0)
              CALL MPL_WAIT(KCOMBUF_S(1:MPLENGTH,1),KREQUEST=ISENDREQ,
     &                CDSTRING='WAIT IN MPDECOMP 5:')
              CALL GSTATS(694,1)
            ENDIF
          ENDIF

        ELSE
          IF(ALLOCATED(KOBSCOR)) DEALLOCATE(KOBSCOR)
        ENDIF


        IF(ALLOCATED(KCOMBUF_S)) DEALLOCATE(KCOMBUF_S)
        IF(ALLOCATED(KCOMBUF_R)) DEALLOCATE(KCOMBUF_R)

        IF(.NOT.LL1D.AND.NPR.GT.1) DEALLOCATE(KDUM)

      ELSEIF(.NOT.LMESSPASS) THEN
        ALLOCATE(WLAT(NIBLO,2))
        ALLOCATE(KOBSLAT(NIBLO,NFRE,2))
        ALLOCATE(KOBSLON(NIBLO,NFRE,2))
        ALLOCATE(WRLAT(NIBLO,2))
        ALLOCATE(WRLON(NIBLO,2))
        ALLOCATE(KOBSRLAT(NIBLO,NFRE,2))
        ALLOCATE(KOBSRLON(NIBLO,NFRE,2))
        ALLOCATE(KOBSCOR(NIBLO,NFRE,4))
        CALL GSTATS(1771,0)
        READ (IU08) WLAT
        READ (IU08) WRLAT, WRLON
        READ (IU08) KOBSLAT
        READ (IU08) KOBSLON
        READ (IU08) KOBSRLAT
        READ (IU08) KOBSRLON
        READ (IU08) KOBSCOR
        CALL GSTATS(1771,1)
!       REMOVE IF NOT USED IN PROPAGATION SCHEME
        IF(IPROPAGS.NE.1) THEN
          DEALLOCATE(WRLAT)
          DEALLOCATE(WRLON)
        ENDIF
        IF(IPROPAGS.NE.1) THEN
          DEALLOCATE(KOBSRLAT)
          DEALLOCATE(KOBSRLON)
        ENDIF
        IF(IPROPAGS.NE.2) THEN
          DEALLOCATE(KOBSCOR)
        ENDIF
      ENDIF


      IF(.NOT.LL1D.AND.NPR.GT.1) THEN
        DEALLOCATE(NEWIJ2IJ)
      ENDIF


!     OBSTRUCTION COEFFICIENTS

!     NOTE: THE VALUE OF OBSLON WILL BE RESET IN THE FIRST
!     CALL TO PROPAGS TO CONTAIN THE OBSTRUCTION TIME THE GROUP
!     VELOCITY At THE INTERFACE !!!!!!!!
      IF (ALLOCATED(OBSLON)) DEALLOCATE(OBSLON)
      ALLOCATE(OBSLON(NSTART(IRANK):NEND(IRANK),NFRE,2))

      DO IC=1,2
        DO M=1,NFRE
          DO IJ=NSTART(IRANK),NEND(IRANK)
            IF(.NOT. LSUBGRID) THEN
              OBSLON(IJ,M,IC)=1.0
            ELSEIF(KOBSLON(IJ,M,IC).EQ.0) THEN
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

!     NOTE: THE VALUE OF OBSLAT WILL BE RESET IN THE FIRST
!     CALL TO PROPAGS TO CONTAIN THE OBSTRUCTION TIME THE GROUP
!     VELOCITY At THE INTERFACE !!!!!!!!
      IF (ALLOCATED(OBSLAT)) DEALLOCATE(OBSLAT)
      ALLOCATE(OBSLAT(NSTART(IRANK):NEND(IRANK),NFRE,2))
      DO IC=1,2
        DO M=1,NFRE
          DO IJ=NSTART(IRANK),NEND(IRANK)
            IF(.NOT. LSUBGRID) THEN
              OBSLAT(IJ,M,IC)=1.0
            ELSE IF(KOBSLAT(IJ,M,IC).EQ.0) THEN
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

      IF(IPROPAGS.EQ.1) THEN
!       NOTE: THE VALUE OF OBSRLON WILL NOT BE RESET IN THE FIRST
        IF(ALLOCATED(OBSRLON)) DEALLOCATE(OBSRLON)
        ALLOCATE(OBSRLON(NSTART(IRANK):NEND(IRANK),NFRE,2))
        DO IC=1,2
          DO M=1,NFRE
            DO IJ=NSTART(IRANK),NEND(IRANK)
              IF(.NOT. LSUBGRID) THEN
                OBSRLON(IJ,M,IC)=1.0
              ELSEIF(KOBSRLON(IJ,M,IC).EQ.0) THEN
                OBSRLON(IJ,M,IC)=0.0
              ELSEIF(MOD(KOBSRLON(IJ,M,IC),1000).EQ.0) THEN
                OBSRLON(IJ,M,IC)=1.0
              ELSE
                OBSRLON(IJ,M,IC)=FLOAT(KOBSRLON(IJ,M,IC))*0.001
              ENDIF
            ENDDO
          ENDDO
        ENDDO

!       NOTE: THE VALUE OF OBSRLAT WILL NOT BE RESET IN THE FIRST
        IF(ALLOCATED(OBSRLAT)) DEALLOCATE(OBSRLAT)
        ALLOCATE(OBSRLAT(NSTART(IRANK):NEND(IRANK),NFRE,2))
        DO IC=1,2
          DO M=1,NFRE
            DO IJ=NSTART(IRANK),NEND(IRANK)
              IF(.NOT. LSUBGRID) THEN
                OBSRLAT(IJ,M,IC)=1.0
              ELSE IF(KOBSRLAT(IJ,M,IC).EQ.0) THEN
                OBSRLAT(IJ,M,IC)=0.0
               ELSEIF(MOD(KOBSRLAT(IJ,M,IC),1000).EQ.0) THEN
              OBSRLAT(IJ,M,IC)=1.0
              ELSE
                OBSRLAT(IJ,M,IC)=FLOAT(KOBSRLAT(IJ,M,IC))*0.001
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(KOBSRLON)
        DEALLOCATE(KOBSRLAT)
      ENDIF

      IF(IPROPAGS.EQ.2) THEN
        IF(ALLOCATED(OBSCOR)) DEALLOCATE(OBSCOR)
        ALLOCATE(OBSCOR(NSTART(IRANK):NEND(IRANK),NFRE,4))
        DO IC=1,4
          DO M=1,NFRE
            DO IJ=NSTART(IRANK),NEND(IRANK)
              IF(.NOT. LSUBGRID) THEN
                OBSCOR(IJ,M,IC)=1.0
              ELSEIF(KOBSCOR(IJ,M,IC).EQ.0) THEN
                OBSCOR(IJ,M,IC)=0.0
              ELSEIF(MOD(KOBSCOR(IJ,M,IC),1000).EQ.0) THEN
                OBSCOR(IJ,M,IC)=1.0
              ELSE
                OBSCOR(IJ,M,IC)=FLOAT(KOBSCOR(IJ,M,IC))*0.001
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(KOBSCOR)
      ENDIF

!     7. DETERMINE KXLTMIN, KXLTMAX

      IF (ALLOCATED(KXLTMIN)) DEALLOCATE(KXLTMIN)
      ALLOCATE(KXLTMIN(NPR))
      IF (ALLOCATED(KXLTMAX)) DEALLOCATE(KXLTMAX)
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

      KTAG=KTAG+1

      IF (LHOOK) CALL DR_HOOK('MPDECOMP',1,ZHOOK_HANDLE)

      CLOSE (UNIT=IU08, STATUS='KEEP')

      DEALLOCATE(ISENDREQ)

      WRITE(IU06,*) ' MODEL DECOMPOSITION FINISHED.'
      CALL FLUSH(IU06)

      RETURN
      END SUBROUTINE MPDECOMP
