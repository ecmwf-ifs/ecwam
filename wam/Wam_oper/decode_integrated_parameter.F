!
! ----------------------------------------------------------------------
!
      PROGRAM DECODE_INTEGRATED_PARAMETER
!
!     J BIDLOT ECMWF   OCTOBRE 1998
!
!     PURPOSE.
!     --------
!     decodes the grib integrated parameters as obtained from a mars request
!     Input parameter can be any integrated parameter as output from WAM.
!     or the atmospheric surface fields (winds) on a regular lat-lon grid.
!     The results will be saved as in ascii as a header and 2d matrix.
!     If more than one field is present, they will written out in the
!     same output file.
!     (the data will be written fron north to south)

!     USAGE : decode_integrated_parameter [-i input_filename]
!                                         [-o output_filename]
!                                         [-t itest] 
!    
!             default value for input_filename : input 
!                               output_filename : output 
!                               itest : 0 (no diagnostics)
!                               itest : 0 (no diagnostics)
!                               itest : 1 (some diagnostics)
!                               itest : 2 (some more diagnostics)
!                               itest : 3 (all diagnostics)

!            itest = 1 : prints built in diagnostics (only indication to calls)
!            itest = 2 : prints built in diagnostics
!            itest = 3 : prints built in diagnostics + gribex debugger 
!                        messages (!!! it can be big)

 
!     LIBRARY : FORTRAN 90, and EMOSLIB
!     -------

!     OUTPUT FILE FORMAT :
!     --------------------
!     2d matrix in binary format of dimension NGX x NGY. NGX number of
!     longitudes and NGY number of latitudes of the regular grid.
!     Missing data will be set to -999.


!     MARS REQUEST EXAMPLES :
!     -----------------------

!    RETRIEVE,
!      DATE=19960704,
!      TIME=12,
!      FORMAT=PACKED,
!      TARGET="input",
!      PARAM=229, !!!! could 220 --> 245
!      DOMAIN=G,    !!! could be M
!      REPRES=LL,
!      STREAM=WV,
!      CLASS=OD,
!      EXPVER=1,
!      TYPE=AN,
!      LEVTYPE=SFC,
!      GRID=1.0/1.0,
!      AREA=81/0/45/30,
!      LEVELIST=OFF


!     FOR DATA AFTER DECEMBER 4, 1996 18Z, and DOMAIN=G or
!     FOR DATA AFTER OCTOBER 27, 1998 0Z, and DOMAIN=M
!     if GRID and AREA are not specified then the data will be retrieved
!     on the original irregular lat lon grid. The output matrix will
!     still be NGX x NGY of the corresponding regular grid and only the
!     column 1 to column N will contain non missing data value per
!     row where N is the number of point on the corresponding latitude.  

! ----------------------------------------------------------------------

      PARAMETER(NFREMAX=36)
      PARAMETER (G = 9.806, PI = 3.1415927, CIRC = 40000000.,
     1           ZPI = 2.*PI, RAD = PI/180., DEG = 180./PI,
     2           R = CIRC/ZPI, EPSMIN=0.1E-32, ZMISS=-999.0)

      INTEGER I4
      INTEGER GETCLO, GETCLA, IOPTVAL
      INTEGER :: KSEC0(2),KSEC1(2048),KSEC2(2060),
     &           KSEC3(2),KSEC4(252)
      INTEGER, ALLOCATABLE :: INGRIB(:),OUTGRIB(:),KDOMRGG(:)
      REAL :: ONE,ZTHETA,ZFRE
      REAL :: PSEC2(96),PSEC3(2)
      REAL, ALLOCATABLE :: PSEC4OUT(:)
      REAL, ALLOCATABLE :: FIELD(:,:)
      CHARACTER     CLOPTLET
      CHARACTER*  3 CLL1
      CHARACTER*  6 CLOPTS
      CHARACTER* 12 CLFMT
      CHARACTER*40 ERRMSG(-4:6), MSG
      CHARACTER*128 CLARG, FNAMEIN, FNAMEOUT
      LOGICAL LLEXIST, LINTEGRATE, LWAVEHGT, LLASTFRE, llfrst
      LOGICAL :: LLSCANNS

      DATA CLOPTS/'i;o;t;'/
      DATA ERRMSG/
     4            ' DECODED WITH BIT MAP                  *',
     3            '                                       *',
     2            '                                       *',
     1            '                                       *',
     1            ' NO ERROR                              *',
     1            ' END OF FILE ENCOUNTED                 *',
     2            ' DECODING ERROR SEE GRIBEX DESCRIPTION *',
     3            ' SUSPICIOUS TIME UNIT IN BLOCK 1       *',
     4            ' DIMENSION 1   IS TOO SMALL            *',
     5            ' error  *',
     6            ' ENCODING ERROR SEE GRIBEX DESCRIPTION *'/

! ----------------------------------------------------------------------

!*    INITIAL VALUES SET AND CRACK COMMAND LINE.
!     -----------------------------------------

      I4=1
      NPRECI = KIND(I4) 
      FNAMEIN='input'
      FNAMEOUT='output'
      IFRE_FIRST=1
      IFRE_LAST=NFREMAX
      LINTEGRATE=.FALSE.
      LWAVEHGT=.FALSE.
      LLASTFRE=.FALSE.
      ITEST=0
      llfrst=.true.

      CMDLINE: DO
        IOPTVAL=GETCLO(CLOPTS,CLARG)
        IF (IOPTVAL .LE. 0 )  THEN
          EXIT CMDLINE
        ENDIF
        CLOPTLET=CHAR(IOPTVAL)
!       GETS VARIABLE ARGUMENT FOR OPTION
        MORARG=GETCLA(CLARG)
        IF (MORARG.NE.0) THEN
          IF ( CLOPTLET .EQ. 'i' ) THEN
            FNAMEIN=CLARG
          ELSE IF ( CLOPTLET .EQ. 'o' ) THEN
            FNAMEOUT=CLARG
          ELSE IF ( CLOPTLET .EQ. 't' ) THEN
            I1=LEN_TRIM(CLARG)
            WRITE (CLL1,'(I3)') I1
            CLFMT = '(I'//CLL1//')'
            READ (CLARG(1:I1),FMT=CLFMT) ITEST
          ENDIF
        ENDIF
      ENDDO CMDLINE

      KSEC1=0
      KSEC2=0
      KSEC4=0
      KSEC3(1)=0
      KSEC3(2)=0
      PSEC2=0.
      PSEC3(2)=ZMISS

!*    INPUT GRIB DATA FILE
!     --------------------
      LFILE=0
      LLEXIST=.FALSE.
      IF (FNAMEIN.NE. ' ') LFILE=LEN_TRIM(FNAMEIN)
      INQUIRE(FILE=FNAMEIN(1:LFILE),EXIST=LLEXIST)
      IF(LLEXIST) THEN
        CALL PBOPEN(IUGRS,FNAMEIN(1:LFILE),'r',KRET)
        IF(KRET.LT.0) THEN
          WRITE (*,*) '****************************************'
          WRITE (*,*) '*                                      *'
          WRITE (*,*) '*   ERROR FOLLOWING CALL TO PBOPEN     *'
          IF(KRET.EQ.-1)
     &      WRITE (*,*) 'COULD NOT OPEN FILE ',FNAMEIN
          IF(KRET.EQ.-2)
     &      WRITE (*,*) 'INVALID FILENAME ',FNAMEIN
          IF(KRET.EQ.-3) WRITE (*,*) 'INVALID OPEN MODE SPECIFIED'
          WRITE (*,*) '*                                      *'
          WRITE (*,*) '****************************************'
          CALL ABORT
        ENDIF
      ELSE
        WRITE(*,*)'****************************'
        WRITE(*,*)'*                          *'
        WRITE(*,*)'*GRIB DATA NOT FOUND IN *'
        WRITE(*,*)  FNAMEIN 
        WRITE(*,*)'*PROGRAM WILL ABORT        *'
        WRITE(*,*)'*                          *'
        WRITE(*,*)'****************************'
        CALL ABORT
      ENDIF
 
!     GET SIZE OF INGRIB


      CALL PBSIZE (IUGRS, IPLENG)
      ISIZE=(IPLENG+NPRECI-1)/NPRECI
1111  IPLENG=ISIZE*NPRECI
      IF(.NOT.ALLOCATED(INGRIB)) ALLOCATE(INGRIB(ISIZE))

!     GET DATA FILE

      CALL PBGRIB(IUGRS,INGRIB,IPLENG,ILENG,KRET)
      IF     (KRET.EQ.-1) THEN
        WRITE (*,*) ' REACHED EOF IN ',FNAMEIN 
        GOTO 2222
      ELSEIF (KRET .EQ. -2) THEN
        WRITE (*,*) ' ERROR IN FILE HANDLING IN ',FNAMEIN 
      ELSEIF (KRET .EQ. -3) THEN
        DEALLOCATE(INGRIB)

        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' SIZE  OF KGRIB IS NOT BIG ENOUGH.'
        WRITE(IU06,*) ' IT WAS ', ISIZE

        KKRET=0

        ISIZE=(ILENG+KIND(ISIZE)-1)/KIND(ISIZE)
        WRITE(IU06,*) ' IT SHOULD AT LEAST BE ', ISIZE
        WRITE(IU06,*) ' THE SIZE WAS RESET AUTOMATICALLY'
        WRITE(IU06,*) ' '
        CALL FLUSH(IU06)

!       RESET THE FILE POINTER TO READ FIELD AGAIN
        KOFFSET=-ILENG
        CALL PBSEEK(IUGRS,KOFFSET,1,KKRET)
        IF(KKRET.EQ.-1) THEN
          WRITE(IU06,*) '***********************************'
          WRITE(IU06,*) '*  PBSEEK : END OF FILE ENCOUNTED'
          WRITE(IU06,*) '***********************************'
          CALL ABORT
        ENDIF
        IF(KKRET.EQ.-2) THEN
          WRITE(IU06,*) '***********************************'
          WRITE(IU06,*) '*  PBSEEK : FILE HANDLING ERROR'
          WRITE(IU06,*) '***********************************'
          CALL ABORT
        ENDIF

        GOTO 1111

      ENDIF

!*    GET GRIB HEADERS 

      KRET  = 1
      ILEN1 = 1
      ALLOCATE(PSEC4OUT(ILEN1))
      IF (ITEST.gt.0)
     &     WRITE(*,*)' GETTING GRIB HEADER OF 1st INPUT FIELD' 
      IF (ITEST.GT.2) CALL GRSDBG (1)
      CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &             PSEC4OUT, ILEN1, INGRIB, ILENG , KWORD, 'J', KRET)
      IF (ITEST.GT.1) WRITE(*,*)' GRIBEX DONE status=' , KRET 
      IF(KRET.GT.0) THEN
        MSG = ERRMSG(2)
        CALL GRPRS0 (KSEC0)
        CALL GRPRS1 (KSEC0, KSEC1)
        WRITE (*,*) MSG
        CALL ABORT
      ENDIF
      IF(ALLOCATED(PSEC4OUT)) DEALLOCATE(PSEC4OUT)

!*    DETERMINE DATA FIELD CHARACTERISTICS 

      IPARAM = KSEC1(6)
      IRGG = KSEC2(17)

      IF(KSEC2(11).EQ.0) THEN
        LLSCANNS=.TRUE.
      ELSEIF(KSEC2(11).EQ.64 .OR. KSEC2(11).EQ.2) THEN
        LLSCANNS=.FALSE.
      ELSE
        WRITE(IU06,*) '***********************************'
        WRITE(IU06,*) '*  SCANNING MODE NOT RECOGNIZED !!!'
        WRITE(IU06,*) ' KSEC2(11) = ', KSEC2(11)
        WRITE(IU06,*) '***********************************'
        CALL ABORT
      ENDIF

      YFRST = FLOAT(KSEC2(4)/1000)+0.1*(MOD(KSEC2(4),1000)/100)+
     &         0.01*(MOD(KSEC2(4),100)/10)+0.001*MOD(KSEC2(4),10)
      YLAST = FLOAT(KSEC2(7)/1000)+0.1*(MOD(KSEC2(7),1000)/100)+
     &         0.01*(MOD(KSEC2(7),100)/10)+0.001*MOD(KSEC2(7),10)
      IF(LLSCANNS) THEN
        AMONOP = YFRST 
        AMOSOP = YLAST 
      ELSE
        AMONOP = YLAST 
        AMOSOP = YFRST 
      ENDIF
      AMOWEP = FLOAT(KSEC2(5)/1000)+0.1*(MOD(KSEC2(5),1000)/100)+
     &         0.01*(MOD(KSEC2(5),100)/10)+0.001*MOD(KSEC2(5),10)

      AMOEAP = FLOAT(KSEC2(8)/1000)+0.1*(MOD(KSEC2(8),1000)/100)+
     &         0.01*(MOD(KSEC2(8),100)/10)+0.001*MOD(KSEC2(8),10)

      XDELLA = FLOAT(KSEC2(10))/1000
      XDELLO = FLOAT(KSEC2(9))/1000

      IYYYY=(KSEC1(21)-1)*100+KSEC1(10)
      IMM=KSEC1(11)
      IDD=KSEC1(12)
      IHH=KSEC1(13)
      IMI=KSEC1(14)

!     FORECAST STEP (defined here in hours)
      IF(KSEC1(15) .EQ. 1 ) THEN
        FCST=KSEC1(16)
      ELSEIF(KSEC1(15) .EQ. 0 ) THEN
        FCST=KSEC1(16)/60.
      ELSEIF(KSEC1(15) .EQ. 2 ) THEN
        FCST=24.*KSEC1(16)
      ELSEIF(KSEC1(15) .EQ. 12 ) THEN
        FCST=12.*KSEC1(16)
      ELSEIF(KSEC1(15) .EQ. 11 ) THEN
        FCST=6.*KSEC1(16)
      ELSEIF(KSEC1(15) .EQ. 10 ) THEN
        FCST=3.*KSEC1(16)
      ELSE
        WRITE(*,*) 'UNKNOWN DEFINITION OF FORECAST STEP !!!'
        WRITE(*,*) 'KSEC1(15) = ',KSEC1(15)
        CALL ABORT
      ENDIF

      NGY = KSEC2(3)
      IF (ITEST.gt.0) WRITE(*,*) ' THE INPUT PARAMETER IS ',IPARAM 

      ALLOCATE(KDOMRGG(NGY))
      KDOMRGG=0
      IF(IRGG.EQ.0) THEN
        NGX = KSEC2(2)
        KDOMRGG = NGX
      ELSE
        ISTART=0
        DO WHILE(KSEC2(23+ISTART).EQ.0)
           ISTART=ISTART+1
        ENDDO
        NGX = 0
        DO J=1,NGY-ISTART
           KDOMRGG(J) = KSEC2(22+J+ISTART)
           NGX = MAX(NGX,KDOMRGG(J))
        ENDDO
      ENDIF

!     DECODE INPUT GRIB DATA
!     ----------------------


      IF(IRGG.EQ.0) THEN
        ILEN1 = NGY*NGX
      ELSE
        ILEN1=0
        DO J=1,NGY
          ILEN1 = ILEN1 + KDOMRGG(J) 
        ENDDO
      ENDIF

      ALLOCATE(PSEC4OUT(ILEN1))
      PSEC4OUT=ZMISS

      KSEC3(2)=0
      PSEC3(2)=ZMISS
      IF (ITEST.gt.0) 
     &   WRITE(*,*)' DECODING INPUT FIELD '
      CALL GRIBEX(KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &            PSEC4OUT, ILEN1, INGRIB, ILENG , KWORD, 'D', KRET)
      IF (ITEST.gt.1) THEN 
        WRITE(*,*)' GRIBEX DONE status=' , KRET 
        WRITE(*,*) ' KSEC0 : '
        WRITE(*,*) KSEC0
        WRITE(*,*) ' KSEC1 : '
        WRITE(*,*) KSEC1
        WRITE(*,*) ' KSEC2 : '
        WRITE(*,*) KSEC2
        WRITE(*,*) ' KSEC3 : '
        WRITE(*,*) KSEC3
        WRITE(*,*) ' KSEC4 : '
        WRITE(*,*) KSEC4
        WRITE(*,*) ' PSEC2 : '
        WRITE(*,*) PSEC2
      ENDIF
      IF(KRET.GT.0) THEN
        MSG = ERRMSG(2)
        CALL GRPRS0 (KSEC0)
        CALL GRPRS1 (KSEC0, KSEC1)
        WRITE (*,*) MSG
        CALL ABORT
      ENDIF


!     WRITES OUT DECODED DATA IN FORMATTED FORMAT
!     *******************************************

      ALLOCATE(FIELD(NGX,NGY))
!     BY DEFINITION FIELD IS FRON NORTH to SOUTH
      FIELD=ZMISS
      IF(LLSCANNS) THEN
        L=0
        DO K=1,NGY
          JSN=NGY-K+1
          DO I=1,KDOMRGG(JSN)
            L = L+1
            FIELD(I,K) = PSEC4OUT(L)
          ENDDO
        ENDDO
      ELSE
        L=0
        DO K=NGY,1,-1
          JSN=K
          DO I=1,KDOMRGG(JSN)
            L = L+1
            FIELD(I,K) = PSEC4OUT(L)
          ENDDO
        ENDDO
      ENDIF

      if(llfrst) then
        LFILE=0
        IF (FNAMEOUT.NE. ' ') LFILE=LEN_TRIM(FNAMEOUT)
        OPEN(10,FILE=FNAMEOUT(1:LFILE),STATUS='UNKNOWN',
     &       FORM='FORMATTED')
        llfrst=.false.
        WRITE(*,*) ' '
        WRITE(*,*) ' THE INPUT PARAMETER IS ',IPARAM 
        WRITE(*,*) ' FOR DATE ',IYYYY,IMM,IDD,IHH,IMI,
     &             ' STEP ',INT(FCST),' HOURS'
        WRITE(*,*) ' '
        WRITE(*,*) ' OUTPUT MATRIX DIMENSION : ',NGX,NGY
        WRITE(*,*) ' '
        WRITE(*,*) ' AMONOP = ',AMONOP
        WRITE(*,*) ' AMOWEP = ',AMOWEP
        WRITE(*,*) ' AMOSOP = ',AMOSOP
        WRITE(*,*) ' AMOEAP = ',AMOEAP
        IF(IRGG.EQ.0) THEN
          WRITE(*,*) ' REGULAR LAT LON GRID '
          WRITE(*,*) ' XDELLO = ',XDELLO
          WRITE(*,*) ' XDELLA = ',XDELLA
        ELSE
          WRITE(*,*) ' IRREGULAR LAT LON GRID '
          WRITE(*,*) ' XDELLA = ',XDELLA
        ENDIF
      else
        WRITE(*,*) ' '
        WRITE(*,*) ' THE INPUT PARAMETER IS ',IPARAM 
        WRITE(*,*) ' FOR DATE ',IYYYY,IMM,IDD,IHH,IMI,
     &             ' STEP ',INT(FCST),' HOURS'
      endif

      DO I=1,NGX
      DO K=1,NGY
        if(field(i,k).ne.zmiss) then
        WRITE(10,110) IPARAM,IYYYY,IMM,IDD,IHH,I,K,FIELD(I,K)
        endif
      ENDDO 
      ENDDO 

110   FORMAT(I3,1X,I4.4,3(I2.2),1x,2I4,1x,f12.8)


      DEALLOCATE(PSEC4OUT)
      DEALLOCATE(KDOMRGG)
      DEALLOCATE(FIELD)
      GOTO 1111


2222  CONTINUE
      CALL PBCLOSE(IUGRS,KRET)
      IF(KRET.LT.0) THEN
        WRITE (*,*) '************************************'
        WRITE (*,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
        WRITE (*,*) '* FILE ',FNAMEIN
        WRITE (*,*) '************************************'
        CALL ABORT
      ENDIF

      CLOSE(10)

      END

!#######################################################################
!z    FUNCTION getclo(yaoptions, yaargument)
!z    INTEGER getclo, getcla, rtb
!z    CHARACTER*   1 yolastarg
!z    CHARACTER* (*) yaoptions, yaargument
!z    CHARACTER* 120 arg

!z    INTEGER here, imorearg, ivarg
!z    DATA here, imorearg, ivarg, arg / 1, 0, 0, "  " /
!z    DATA yolastarg / " " /

!z    arg=' '
!z    CALL getarg(here,arg)
!-->  PRINT*,'-------------getclo---------------'
!-->  PRINT*,' ###',arg,'###'

!z    iol=rtb(arg)
!-->  PRINT*,' iol: ', iol
!z    IF (iol .EQ. 2 .AND. arg(1:1) .EQ. '-' .AND. ivarg .EQ. 0 ) THEN
!z      iol = rtb(yaoptions)
!-->    PRINT*,' iol 2',iol,' options: ', yaoptions
!z      DO jl=1,iol
!z        getclo = 0
!z        IF ( yaoptions(jl:jl) .EQ. arg(2:2) ) THEN
!z          getclo = ichar(arg(2:2))
!-->        PRINT *,' FOUND ',yaoptions(jl:jl), ' IN THE COMMAND LINE',
!--> .               yaoptions(jl+1:jl+1)
!z          IF (yaoptions(jl+1:jl+1) .EQ. ':' ) THEN
!-->          PRINT*, yaoptions(jl:jl),' requires arguments'
!z            yolastarg=yaoptions(jl:jl)
!z            ivarg=1
!z          ENDIF
!z          EXIT
!z        ENDIF
!z      ENDDO
!z    ELSEIF ( ivarg .EQ. 1 ) THEN
!z       WRITE(*,*) ' option -', yolastarg, ' requires arguments'
!z       getclo=-1
!z    ELSEIF (iol .EQ. 0) THEN
!z      getclo=0
!z    ELSE
!z       WRITE(*,*) 'illegal option: ',arg(1:iol)
!z       getclo=-1
!z    ENDIF
!-->  PRINT *,' HERE ins getclo', here
!z    here = here + 1
!z    RETURN

!z    ENTRY getcla(yaargument)
!-->  PRINT*,'-------------getcla--------------'
!-->  PRINT*, 'HERE ins getcla :', here,' options: ', yaoptions

!z    getcla = 1
!z    CALL getarg(here,arg)
!-->  PRINT*,' arg in getcla ', arg
!z    IF ( arg (1:1) .NE. '-' ) THEN
!z      here = here + 1
!z      yaargument=arg
!z    ELSE
!z      IF (ivarg.EQ.1) THEN
!z        WRITE(*,*)' refused to take ', arg (1:2) ,' as argument for',
!z   .    ' the option -',yolastarg
!z        getcla = -1
!z      ELSE
!z        getcla = 0
!z      ENDIF
!z    ENDIF
!z    ivarg=0
!-->  PRINT*,' getcla in getcla ', getcla

!z    RETURN
!z    END
