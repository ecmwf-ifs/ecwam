! ----------------------------------------------------------------------

      PROGRAM WRITE_MPDECOMP

! ----------------------------------------------------------------------

!**** *WRITE_MPDECOMP* - WRITES AN ASCII TEXT FILE TO PLOT
!                        THE MODEL DECOMPOSITION FOR A GIVEN
!                        NUMBER OF PE and other stats etc..

!**   INTERFACE.
!     ----------

!          ---- FORMAL PARAMETERS ----

!     *WRITE_MPDECOMP*  -n nproc [-a] [-i irank] [-d]
!       *nproc*     INTEGER  NUMBER of SUB DOMAIN (i.e PE WHEN RUNNING) 
!       *-a*        all sea points will be written to output
!                   else only the boundary points.
!       *irank*     if present the output will only be limited to
!                   the area the pe irank would have if mpdecomp
!                   was run in the message passing environment
!                   As a debugging tool it will also output an extra
!                   file which contains the points that are to be
!                   exchanged from pe irank with the others.
!       *-d*        if present the domain decomposition will only be
!                   1-dimensional (i.e per latitudinal bands)

!          ---- INPUT/OUTPUT UNITS ---

!           *IU06*   - PRINTER OUTPUT.
!           *IU07*   - INPUT  UNIT OF PRECOMPUTED COMMON BLOCKS.
!                      (OUTPUT OF PREPROC). fort.7 must be there !
!           *IU08*   - INPUT  UNIT OF COMMON UBUF.
!                      (OUTPUT OF PREPROC). fort.8 must be there !
!           *IUOUT*  - ASCII TEXT FILE grid_decomposition.xxxx
!            where xxxx=nproc 

!          MESSAGE PASSING NOT CATERED FOR YET.
!          IF MESSAGE PASSING ENABLED THE UNITS WILL BE REASSIGNED.

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE
!     ---------


! ----------------------------------------------------------------------

      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL      ,IJLT
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,
     &            IPER     ,IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,LLBOUND
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPREVIOUS,NNEXT    ,
     &            MPMAXLENGTH,KTAG     ,NPRECR   ,
     &            NPRECI
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NBLO     ,NIBLO    ,
     &            LL1D
      USE YOWSPEC  , ONLY : NSTART   ,NEND     ,KLENTOP  ,KLENBOT  ,
     &            NFROMPE  ,NTOPE    ,NIJSTART ,IJTOPE
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,KLAT2

! -------------------------------------------------------------------

      INTEGER GETCLO, GETCLA, IOPTVAL
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJPROC, NPE
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NEIGHBOUR
      LOGICAL :: LLALLPTS
      CHARACTER*  3 CLL1
      CHARACTER*  8 CLOPTS
      CHARACTER* 12 CLFMT
      CHARACTER*128 CLARG
      CHARACTER*28 FILENAME
      CHARACTER     CLOPTLET
      DATA CLOPTS/'n;a;i;d;'/

      LOGICAL :: LLIRANK
      REAL Z4(2)
      INTEGER I4(2)

      CHARACTER (LEN=12),ALLOCATABLE :: CLMD(:)

! ----------------------------------------------------------------------

!*    1. INITIAL VALUES SET AND CRACK COMMAND LINE.
!        -----------------------------------------

      CALL INIWCST

      IU06 = 6
      NPROC = 0
      IRANK = 1
      LMESSPASS = .FALSE.
      LL1D = .FALSE.
      LLIRANK = .FALSE.
      LLALLPTS = .FALSE.

      CMDLINE: DO
        IOPTVAL=GETCLO(clopts,clarg)
        IF (IOPTVAL .LE. 0 )  THEN
          EXIT CMDLINE
        ENDIF

        CLOPTLET=CHAR(IOPTVAL)

!       GETS VARIABLE ARGUMENT FOR OPTION
        MORARG=GETCLA(clarg)
        IF (MORARG.NE.0) THEN
          IF ( CLOPTLET .EQ. 'n' ) THEN
            I1=LEN_TRIM(CLARG)
            WRITE (CLL1,'(I3)') I1
            CLFMT = '(I'//CLL1//')'
            READ (CLARG(1:I1),FMT=CLFMT) NPROC
          ELSEIF ( CLOPTLET .EQ. 'a' ) THEN
            LLALLPTS = .TRUE.
          ELSEIF ( CLOPTLET .EQ. 'i' ) THEN
            I1=LEN_TRIM(CLARG)
            WRITE (CLL1,'(I3)') I1
            CLFMT = '(I'//CLL1//')'
            READ (CLARG(1:I1),FMT=CLFMT) IRANK 
            LLIRANK = .TRUE.
          ELSEIF ( CLOPTLET .EQ. 'd' ) THEN
            LL1D = .TRUE.
          ENDIF
        ENDIF

      ENDDO CMDLINE

      IF(NPROC.LE.0) THEN
        WRITE(*,*) 'NPROC MUST BE POSITIVE !!!!'
        CALL ABORT
      ENDIF

      IF(IRANK.LE.0 .OR. IRANK.GT.NPROC) THEN
        WRITE(*,*) 'IRANK MUST BE >0 and <= NPROC !!!!'
        CALL ABORT
      ENDIF

      IF(LL1D) THEN
        WRITE(*,*) '1-d DECOMPOSITION HAS BEEN REQUESTED !'
      ENDIF

! ---------------------------------------------------------------------

!     2. GLOBAL SETTINGS.
!        ----------------


      Z4=1.
      NPRECR = KIND(Z4)
      I4=1
      NPRECI = KIND(I4)

      IREAD=1

      NPREVIOUS=IRANK-1
      IF(IRANK.EQ.NPROC) THEN
        NNEXT=0
      ELSE
        NNEXT=IRANK+1
      ENDIF

      NPR=NPROC

      CALL MPDECOMP(NPR,MAXLEN,LLIRANK)

      MPMAXLENGTH=MAXLEN

! ----------------------------------------------------------------------

      FILENAME='grid_decomposition_'
      WRITE (FILENAME(20:23),'(I4.4)') NPROC
      IF(LLIRANK) THEN
        FILENAME(24:24)='.'
        WRITE (FILENAME(25:28),'(I4.4)') IRANK
      ENDIF
      IUOUT=10
      OPEN(IUOUT,FILE=FILENAME,FORM='FORMATTED')

      WRITE(IUOUT,'(a4)') '#GEO'
      WRITE(IUOUT,'(a11)') '#FORMAT LLV'
      WRITE(IUOUT,'(a5)') '#DATA'

      IG=1
      IPROC=0

      ALLOCATE(IJPROC(NIBLO))
      
      IF(LLIRANK) THEN
        IPROC_S=IRANK
        IPROC_E=IRANK
      ELSE
        IPROC_S=1
        IPROC_E=NPROC
      ENDIF

      DO IPROC=IPROC_S,IPROC_E
        DO IJ=NSTART(IPROC),NEND(IPROC)
          IX = IXLG(IJ,IG)
          JSN= KXLT(IJ,IG)
          XLON=AMOWEP+(IX-1)*ZDELLO(JSN)
          XLAT=AMOSOP + (JSN-1)*XDELLA
          isymbol=mod(IPROC,9)
          if(isymbol.eq.0) isymbol=9

          IF(LLALLPTS) THEN
            WRITE(IUOUT,'(2(1X,F8.3),1X,I4)') XLON,XLAT,isymbol
          ELSE
            IF(LLBOUND(IJ))
     &        WRITE(IUOUT,'(2(1X,F8.3),1X,I4)') XLON,XLAT,isymbol
          ENDIF
        ENDDO
      ENDDO

      CLOSE(IUOUT)

      DO IPROC=1,NPROC
        DO IJ=NSTART(IPROC),NEND(IPROC)
          IJPROC(IJ)=IPROC
        ENDDO
      ENDDO


      IF(LLIRANK) THEN
        FILENAME='pts_in_message_'
        WRITE (FILENAME(16:19),'(I4.4)') NPROC
        FILENAME(20:20)='.'
        WRITE (FILENAME(21:24),'(I4.4)') IRANK
        IUOUT=10
        OPEN(IUOUT,FILE=FILENAME,FORM='FORMATTED')

        WRITE(IUOUT,'(a4)') '#GEO'
        WRITE(IUOUT,'(a11)') '#FORMAT LLV'
        WRITE(IUOUT,'(a5)') '#DATA'

        DO IPROC=1,NPROC
          DO IH=1,NTOPE(IPROC)

          IJ = IJTOPE(IH,IPROC)
          IX = IXLG(IJ,IG)
          JSN= KXLT(IJ,IG)
          XLON=AMOWEP+(IX-1)*ZDELLO(JSN)
          XLAT=AMOSOP + (JSN-1)*XDELLA
          ISYMBOL=MOD(IPROC,9)
          IF(ISYMBOL.EQ.0) ISYMBOL=9
          WRITE(IUOUT,'(2(1X,F8.3),1X,I4)') XLON,XLAT,isymbol

          ENDDO
        ENDDO
        CLOSE(IUOUT)
      ENDIF

!!!   find max number of neighbours
      ALLOCATE(NEIGHBOUR(NPROC,NPROC))
      ALLOCATE(NPE(NPROC))
      NEIGHBOUR=0
      NPE=0
      DO IPR=IPROC_S,IPROC_E
        DO IJ=NSTART(IPR),NEND(IPR)
          DO IC=1,2
            IF(KLON(IJ,IC).GT.0.AND.KLON(IJ,IC).LE.NIBLO)
     &        NEIGHBOUR(IJPROC(IJ),IJPROC(KLON(IJ,IC)))=1 
            IF(KLAT(IJ,IC).GT.0.AND.KLAT(IJ,IC).LE.NIBLO)
     &        NEIGHBOUR(IJPROC(IJ),IJPROC(KLAT(IJ,IC)))=1 
            IF(KLAT2(IJ,IC).GT.0.AND.KLAT2(IJ,IC).LE.NIBLO)
     &        NEIGHBOUR(IJPROC(IJ),IJPROC(KLAT2(IJ,IC)))=1 
          ENDDO
        ENDDO
        NPE(IPR)=NEND(IPR)-NSTART(IPR)+1
      ENDDO

      FILENAME='neighbours_'
      WRITE (FILENAME(12:15),'(I4.4)') NPROC 
      IF(LLIRANK) THEN
        FILENAME(16:16)='.'
        WRITE (FILENAME(17:20),'(I4.4)') IRANK
      ENDIF
      IUOUN=20
      OPEN(IUOUN,FILE=FILENAME,FORM='FORMATTED')

      XNGHB_MEAN=0.
      SNGHB=0.
      XNPE=0.
      DO IPR=IPROC_S,IPROC_E
        NGHB=0
        DO JPR=IPROC_S,IPROC_E
           NGHB=NGHB+NEIGHBOUR(IPR,JPR)
        ENDDO
        IF(NPE(IPR).gt.1) NGHB=NGHB-1
        XNGHB_MEAN=XNGHB_MEAN+NGHB
        SNGHB=SNGHB+NGHB*NGHB
         XNPE=XNPE+NPE(IPR)
        WRITE(IUOUN,'(i4.4,1x,i6,1x,i4)') IPR,NPE(IPR),NGHB
      ENDDO
      XNGHB_MEAN=XNGHB_MEAN/(IPROC_E-IPROC_S+1)
      SNGHB=SNGHB/(IPROC_E-IPROC_S+1)
      SNGHB=SNGHB-XNGHB_MEAN**2
      SNGHB=SQRT(MAX(SNGHB,0.))
      XNPE=XNPE/(IPROC_E-IPROC_S+1)

      WRITE(*,*) 'in write_mpdecomp ',NPROC,XNGHB_MEAN,SNGHB,XNPE

      CLOSE(IUOUN)
      DEALLOCATE(IJPROC)
      DEALLOCATE(NEIGHBOUR)
      DEALLOCATE(NPE)

      CLOSE(IUOUT)

      END
