!
! ----------------------------------------------------------------------
!
      PROGRAM DECODE_INTEGRATED_PARAMETER
!
!     J BIDLOT ECMWF   OCTOBRE 1998
!     J BIDLOT ECMWF   MAY 2009: CONVERSION TO GRIB_API 
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

 
!     LIBRARY : FORTRAN 90, GRIB_API ($GRIB_API_INCLUDE $GRIB_API_LIB)
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

      USE GRIB_API

      IMPLICIT NONE

      REAL, PARAMETER :: ZMISS=-999.0

      INTEGER :: MORARG, IU06
      INTEGER :: GETCLO, GETCLA, IOPTVAL
      INTEGER :: I1, ITEST, LFILE
      INTEGER :: NGX, NGY, L, K, I, J, JSN, ISTART
      INTEGER :: IFILE, IRET, IGRIB 
      INTEGER :: ITABPAR, IPARAM, ITABLE, IRGG, ISCAN
      INTEGER :: IYYYYMMDD, IHHMM, IHH
      INTEGER :: NUMBEROFVALUES
      INTEGER :: PLPRESENT, NB_PL
      INTEGER, DIMENSION(:), ALLOCATABLE :: PL
      INTEGER, DIMENSION(:), ALLOCATABLE :: KDOMRGG

      REAL :: YFRST, YLAST, AMONOP, AMOSOP, AMOWEP, AMOEAP
      REAL :: XDELLO, XDELLA
      REAL :: FCST, STEP
      REAL(KIND=8), ALLOCATABLE :: VALUES(:)
      REAL, ALLOCATABLE :: FIELD(:,:)

      CHARACTER*1  :: CLOPTLET
      CHARACTER*3  :: CLL1
      CHARACTER*6  :: CLOPTS
      CHARACTER*8  :: CSTEPUNITS
      CHARACTER*8  :: CSTEPTYPE
      CHARACTER*12 :: CGRIDTYPE
      CHARACTER*12 :: CLFMT
      CHARACTER*128:: CLARG, FNAMEIN, FNAMEOUT

      LOGICAL :: LLEXIST, LLFRST
      LOGICAL :: LLSCANNS

      DATA CLOPTS/'i;o;t;'/

! ----------------------------------------------------------------------

!*    INITIAL VALUES SET AND CRACK COMMAND LINE.
!     -----------------------------------------

      FNAMEIN='input'
      FNAMEOUT='output'
      ITEST=0
      LLFRST=.TRUE.

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*    OPEN DATA FILE
!     --------------
      LFILE=0
      LLEXIST=.FALSE.
      IF (FNAMEIN.NE. ' ') LFILE=LEN_TRIM(FNAMEIN)
      INQUIRE(FILE=FNAMEIN(1:LFILE),EXIST=LLEXIST)
      IF(LLEXIST) THEN
        CALL GRIB_OPEN_FILE(IFILE,FNAMEIN(1:LFILE),'r')
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

!     LOOP ON ALL MESSAGES IN INPUT FILE
 
      IGRIB=-99
      CALL  GRIB_NEW_FROM_FILE(IFILE,IGRIB,IRET)

      LOOP: DO WHILE (IRET /= GRIB_END_OF_FILE)

!       GET THE SIZE OF THE VALUES ARRAY
        CALL GRIB_GET_SIZE(IGRIB,'values',NUMBEROFVALUES)

        ALLOCATE(VALUES(NUMBEROFVALUES))

!       SET THE MISSING DATA INDICATOR
        CALL GRIB_SET(IGRIB, 'missingValue',ZMISS)

!       GET DATA VALUES
        CALL GRIB_GET(IGRIB,'values',VALUES)

!       DETERMINE DATA FIELD CHARACTERISTICS 

        CALL GRIB_GET(IGRIB,'paramId',ITABPAR)
        IF (ITEST.gt.0) WRITE(*,*) ' THE INPUT PARAMETER IS ',ITABPAR
        ITABLE=ITABPAR/1000
        IPARAM=ITABPAR-ITABLE*1000

        CALL GRIB_GET(IGRIB,'gridType', CGRIDTYPE)
        IF(CGRIDTYPE(1:7) == 'regular') THEN
          IRGG=0
        ELSEIF(CGRIDTYPE(1:7) == 'reduced') THEN
          IRGG=1
        ELSE
          WRITE(IU06,*) '***********************************'
          WRITE(IU06,*) '*  GRID TYPE NOT RECOGNIZED !!!'
          WRITE(IU06,*) '   gridType = ', CGRIDTYPE 
          WRITE(IU06,*) '***********************************'
          CALL ABORT
        ENDIF

        CALL GRIB_GET(IGRIB,'iScansNegatively',ISCAN)
        IF(ISCAN.EQ.0) THEN
          LLSCANNS=.TRUE.
        ELSEIF(ISCAN.EQ.64 .OR. ISCAN.EQ.2) THEN
          LLSCANNS=.FALSE.
        ELSE
          WRITE(IU06,*) '***********************************'
          WRITE(IU06,*) '*  SCANNING MODE NOT RECOGNIZED !!!'
          WRITE(IU06,*) ' ISCAN = ', ISCAN 
          WRITE(IU06,*) '***********************************'
          CALL ABORT
        ENDIF

        CALL GRIB_GET(IGRIB,'latitudeOfFirstGridPointInDegrees',YFRST)
        CALL GRIB_GET(IGRIB,'latitudeOfLastGridPointInDegrees',YLAST)

        IF(LLSCANNS) THEN
          AMONOP = YFRST 
          AMOSOP = YLAST 
        ELSE
          AMONOP = YLAST 
          AMOSOP = YFRST 
        ENDIF

        CALL GRIB_GET(IGRIB,'longitudeOfFirstGridPointInDegrees',AMOWEP)
        CALL GRIB_GET(IGRIB,'longitudeOfLastGridPointInDegrees',AMOEAP)

        IF(CGRIDTYPE(9:10) == 'gg') THEN
          XDELLA=0.0
          XDELLO=0.0
        ELSEIF(CGRIDTYPE(1:7) == 'reduced') THEN
          CALL GRIB_GET(IGRIB,'jDirectionIncrementInDegrees',XDELLA)
          XDELLO=0.0
        ELSE
          CALL GRIB_GET(IGRIB,'jDirectionIncrementInDegrees',XDELLA)
          CALL GRIB_GET(IGRIB,'iDirectionIncrementInDegrees',XDELLO)
        ENDIF


        CALL GRIB_GET(IGRIB,'dataDate',IYYYYMMDD)
        CALL GRIB_GET(IGRIB,'time',IHHMM)
        IHH=IHHMM/100

        CALL GRIB_GET(IGRIB,'stepType',CSTEPTYPE)

!       FORECAST STEP (defined here in hours)
        IF(CSTEPTYPE(1:7) .EQ. 'instant') THEN
          CALL GRIB_SET(IGRIB,'stepUnits','h')
          CALL GRIB_GET(IGRIB,'step',STEP)
          FCST=STEP
        ELSE
          WRITE(*,*) 'UNKNOWN DEFINITION OF FORECAST STEP TYPE !!!'
          WRITE(*,*) 'stepType = ',CSTEPTYPE
          CALL ABORT
        ENDIF


        CALL GRIB_GET(IGRIB,'Ny',NGY)

        ALLOCATE(KDOMRGG(NGY))
        KDOMRGG=0
        IF(IRGG.EQ.0) THEN
          CALL GRIB_GET(IGRIB,'Nx',NGX)
          KDOMRGG = NGX
        ELSE
          CALL GRIB_GET(IGRIB,'PLPresent',PLPRESENT)
          IF (PLPRESENT == 1) THEN
            CALL GRIB_GET_SIZE(IGRIB,'pl',NB_PL)
            ALLOCATE(PL(NB_PL))
            CALL GRIB_GET(IGRIB,'pl',PL)
          ELSE
            WRITE(*,*) 'NUMBER OF POINTS PER LATITUDE MISSING !!!'
            CALL ABORT
          ENDIF

          ISTART=1
          DO WHILE(PL(ISTART).EQ.0 .AND. ISTART.LT.NB_PL)
             ISTART=ISTART+1
          ENDDO
          ISTART=ISTART-1
          NGX = 0
          DO J=1,NGY-ISTART
             KDOMRGG(J) = PL(J+ISTART)
             NGX = MAX(NGX,KDOMRGG(J))
          ENDDO
          DEALLOCATE(PL)
        ENDIF

!       WRITES OUT DECODED DATA IN FORMATTED FORMAT
!       *******************************************

        ALLOCATE(FIELD(NGX,NGY))
!       BY DEFINITION FIELD IS FRON NORTH to SOUTH
        FIELD=ZMISS
        IF(LLSCANNS) THEN
          L=0
          DO K=1,NGY
            JSN=NGY-K+1
            DO I=1,KDOMRGG(JSN)
              L = L+1
            FIELD(I,K) = VALUES(L)
            ENDDO
          ENDDO
        ELSE
          L=0
          DO K=NGY,1,-1
            JSN=K
            DO I=1,KDOMRGG(JSN)
              L = L+1
              FIELD(I,K) = VALUES(L)
            ENDDO
          ENDDO
        ENDIF

        IF(LLFRST) THEN
          LFILE=0
          IF (FNAMEOUT.NE. ' ') LFILE=LEN_TRIM(FNAMEOUT)
          OPEN(10,FILE=FNAMEOUT(1:LFILE),STATUS='UNKNOWN',
     &         FORM='FORMATTED')
          LLFRST=.FALSE.
          WRITE(*,*) ' '
          WRITE(*,*) ' THE INPUT PARAMETER IS ',IPARAM
          WRITE(*,*) ' FOR DATE ',IYYYYMMDD,IHHMM,
     &               ' STEP ',INT(FCST),' HOURS'
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
        ELSE
          WRITE(*,*) ' '
          WRITE(*,*) ' THE INPUT PARAMETER IS ',IPARAM 
          WRITE(*,*) ' FOR DATE ',IYYYYMMDD,IHHMM,
     &               ' STEP ',INT(FCST),' HOURS'
        ENDIF

        DO I=1,NGX
          DO K=1,NGY
            IF(FIELD(I,K).NE.ZMISS) THEN
            WRITE(10,110) IPARAM,IYYYYMMDD,IHH,I,K,FIELD(I,K)
            ENDIF
          ENDDO 
        ENDDO 

110     FORMAT(I3,1X,I8.8,i2.2,1x,2I4,1x,f12.8)

        DEALLOCATE(FIELD)

        CALL GRIB_RELEASE(IGRIB)

        IGRIB=-99
        CALL GRIB_NEW_FROM_FILE(IFILE,IGRIB, IRET)

      END DO LOOP

      CALL GRIB_CLOSE_FILE(IFILE)

      CLOSE(10)

      END

!#######################################################################
      FUNCTION getclo(yaoptions, yaargument)
      INTEGER getclo, getcla, my_rtb
      CHARACTER*   1 yolastarg
      CHARACTER* (*) yaoptions, yaargument
      CHARACTER* 120 arg

      INTEGER here, imorearg, ivarg
      DATA here, imorearg, ivarg, arg / 1, 0, 0, "  " /
      DATA yolastarg / " " /

      arg=' '
      CALL getarg(here,arg)
!-->  PRINT*,'-------------getclo---------------'
!-->  PRINT*,' ###',arg,'###'

      iol=my_rtb(arg)
!-->  PRINT*,' iol: ', iol
      IF (iol .EQ. 2 .AND. arg(1:1) .EQ. '-' .AND. ivarg .EQ. 0 ) THEN
        iol = my_rtb(yaoptions)
!-->    PRINT*,' iol 2',iol,' options: ', yaoptions
        DO jl=1,iol
          getclo = 0
          IF ( yaoptions(jl:jl) .EQ. arg(2:2) ) THEN
            getclo = ichar(arg(2:2))
!-->        PRINT *,' FOUND ',yaoptions(jl:jl), ' IN THE COMMAND LINE',
!--> .               yaoptions(jl+1:jl+1)
            IF (yaoptions(jl+1:jl+1) .EQ. ':' ) THEN
!-->          PRINT*, yaoptions(jl:jl),' requires arguments'
              yolastarg=yaoptions(jl:jl)
              ivarg=1
            ENDIF
            EXIT
          ENDIF
        ENDDO
      ELSEIF ( ivarg .EQ. 1 ) THEN
         WRITE(*,*) ' option -', yolastarg, ' requires arguments'
         getclo=-1
      ELSEIF (iol .EQ. 0) THEN
        getclo=0
      ELSE
         WRITE(*,*) 'illegal option: ',arg(1:iol)
         getclo=-1
      ENDIF
!-->  PRINT *,' HERE ins getclo', here
      here = here + 1
      RETURN

      ENTRY getcla(yaargument)
!-->  PRINT*,'-------------getcla--------------'
!-->  PRINT*, 'HERE ins getcla :', here,' options: ', yaoptions

      getcla = 1
      CALL getarg(here,arg)
!-->  PRINT*,' arg in getcla ', arg
      IF ( arg (1:1) .NE. '-' ) THEN
        here = here + 1
        yaargument=arg
      ELSE
        IF (ivarg.EQ.1) THEN
          WRITE(*,*)' refused to take ', arg (1:2) ,' as argument for',
     .    ' the option -',yolastarg
          getcla = -1
        ELSE
          getcla = 0
        ENDIF
      ENDIF
      ivarg=0
!-->  PRINT*,' getcla in getcla ', getcla

      RETURN
      END
C Copyright 1981-2007 ECMWF
C 
C Licensed under the GNU Lesser General Public License which
C incorporates the terms and conditions of version 3 of the GNU
C General Public License.
C See LICENSE and gpl-3.0.txt for details.
C


      INTEGER FUNCTION MY_RTB (CHAR)
C
C---->
C*
C*    NAME      : MY_RTB
C*
C*    FUNCTION  : COUNT THE NUMBER OF CHARACTERS IN A CHARACTER
C*                STRING, EXCLUDING TRAILING SPACES.
C*
C*    INPUT     : CHAR - CHARACTER STRING
C*
C*    OUTPUT    : MY_RTB  - NO OF CHARACTERS
C*
C*    GENERAL   : MY_RTB CALLS  -----
C*
C*    AUTHOR    : J.HENNESSY  15.4.85
C*
C     ---------------------------------------------------------------
C----<
C
      CHARACTER*(*) CHAR
C
      J = LEN (CHAR)
C
C     No trailing blanks.
C
      IF (CHAR(J:J).NE.' ')
     C     THEN
               MY_RTB = J
               RETURN
           ENDIF
C
C     String with trailing blanks.
C
      DO 100 I = J,1,-1
         IF (CHAR(I:I).NE.' ')
     C          THEN
                    MY_RTB = I
                    RETURN
                ENDIF
  100 CONTINUE
C
C     String consists of only blanks.
C
      MY_RTB = 0
C
      RETURN
      END


