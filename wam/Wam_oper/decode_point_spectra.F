!
! ----------------------------------------------------------------------
!
      PROGRAM DECODE_POINT_SPECTRA 
!
!     J BIDLOT ECMWF   OCTOBRE 1998
!     J BIDLOT ECMWF   MAY 2009: CONVERSION TO GRIB_API 
!
!     PURPOSE.
!     --------
!     decodes the grib spectrum for locations as obtained from a mars
!     request which has limited the area to a few points 
!     (it will be a little big otherwise) (see example below). 
!     Input parameter can be 250 (old) or 251 (NEW).
!     The results will be saved as an ascii text file
!     containing the location and date and time  and the decoded sprectrum
!     The 1d frequency spectrum will also be computed

!     USAGE : decode_point_spectra [-i input_filename] [-o output_filename]
!     ------                       [-t itest] 
!    
!             default value for input_filename : input_spectra 
!                               output_filename : output_spectra 
!                               the 1d spectrum will be saved in 
!                               output_spectra.1d
!                               itest : 0 (no diagnostics)
!                               itest : 0 (no diagnostics)
!                               itest : 1 (some diagnostics)
!                               itest : 2 (some more diagnostics)
!                               itest : 3 (all diagnostics)

!            itest = 1 : prints built in diagnostics (only indication to calls)
!            itest = 2 : prints built in diagnostics
!            itest = 3 : prints built in diagnostics + gribex debugger 
!                        messages (!!! it can be big)

!     INPUT FILE REQUIREMENT :
!     -----------------------
!     The input file can only contain one parameter (250 or 251)
 
!     If the input file contains parameter 251, it must have been obtained
!     for all directions and frequencies in the default order. Namelly, the
!     mars request should be done with DIRECTION=1/TO/nang, and
!     FREQUENCY=1/TO/nfre. At this time nfre=30 and nang=24 for global model
!     and for mediterranean data.
 
 
!     LIBRARY : FORTRAN 90, GRIB_API ($GRIB_API_INCLUDE $GRIB_API_LIB)
!     -------

!     OUTPUT FILE FORMAT :
!     --------------------

!     LAT LON DATETIME
!     0.000 direction1 direction2 ..... direction12
!     freq1 value      value      ......value
!     .
!     .
!     .
!     freq25

!     freq1 to freq25 are the frequencies (Hz)

!     and so on for the other locations

!     values are the spectral components (m**2 Hz-1)

!     The directions are defined such that 0. means that the waves are
!     propagating to the North, 90 to the East, etc...

!     MARS REQUEST EXAMPLES :
!     -----------------------

!     PLEASE NOTE THAT WITH THE INTRODUCTION OF THE COUPLED WAVE MODEL
!     ON June 29, 1998, THE "GLOBAL" WAVE MODEL SPECTRA ARE NOW SAVED IN
!     GRIB AS PARAMETER 251 AND ARE AVAILABLE FOR THE ANALYSIS 
!     AT 0, 6, 12, 18Z. THE DOMAIN is G for Global. The ATTRIBUTE 
!     FREQUENCY=1/TO/nfre AND DIRECTION=1/TO/nang MUST BE SPECIFIED (see below) 

!     From June 28, 1998 18Z until November 20,2000 12Z
!     nfre=25 and nang=12

!     Since November 20, 2000 18Z
!     nfre=30 and nang=24 

!     BEFORE THAT TIME, AND FROM February 1995, 26 12z,
!     PARAMETER 250 IS TO BE USED AND ANALYSIS ARE ONLY AVAILABLE FOR 12Z
!     (except after January 20 1998 12Z when WAVE MODEL
!     ANALYSIS SPECTRA BECAME AVAILABLE AT 0, 6, 12, 18Z). FOR SPECTRAL
!     DATA BETWEEN December 4, 1996, 18Z and June 28, 1998 18z, the DOMAIN
!     WAS SPLIT BETWEEN NORTH AND SOUTH HEMISPHERES (DOMAIN=N OR DOMAIN=S).
!     The EQUATOR IS PART OF THE NORTHERN HEMISPHERE. BEFORE that TIME,
!     THE DOMAIN WAS GLOBAL (DOMAIN=G). NOTE THAT THE DATA FILE FOR 250
!     CAN BE QUITE CUMBERSOM IF RETRIEVED FOR THE ALL GLOBE.

!     THE DATA FROM THE MEDITERRANEAN/BALTIC SEA MODEL ARE ARCHIVED AS
!     PARAMETER 250 FROM February 1995, 26 12z,
!     UNTIL THE INTRODUCTION OF THE NEW MODEL (on October 27 1998)
!     where IT BECAME 251 AS WELL.

!!!   NOTE: THE DISCRETISED DIRECTION STARTING POINT WAS 0 DEGREE BEFORE
!     May 14 1997 12Z. After THAT DATE, IT IS 15 DEGREES FOR THE GLOBAL
!     MODEL AND 7.5 FOR THE MEDITERRANEAN MODEL. SO ALWAYS GET THE
!     INFORMATION FOR DIRECTION (AND FREQUENCY) FROM THE GRIB HEADERS.

!     REANALYSIS DATA FROM ERA15 WERE ARCHIVED AS PARAMETER 250.

!     REANALYSIS DATA FROM ERA40 WERE ARCHIVED AS PARAMETER 251.


!     EXAMPLE:  spectra at (30W,59S) and (29W,59S) from the global model:

!     FOR DATA PRIOR TO DECEMBER 4, 1996 18Z :

!    RETRIEVE,
!      DATE=19960704,
!      TIME=12,
!      FORMAT=PACKED,
!      TARGET="input_spectra",
!      PARAM=250,
!      DOMAIN=G,
!      REPRES=LL,
!      STREAM=WV,
!      CLASS=OD,
!      EXPVER=1,
!      TYPE=AN,
!      LEVTYPE=SFC,
!      GRID=1.0/1.0,
!      AREA=-59/-30/-59/-29,
!      LEVELIST=OFF


!     FOR DATA AFTER DECEMBER 4, 1996 18Z, but PRIOR TO JUNE 29, 1998 0Z :

!    RETRIEVE,
!      DATE=19970704,
!      TIME=12,
!      FORMAT=PACKED,
!      TARGET="input_spectra",
!      PARAM=250,
!      DOMAIN=S,
!      REPRES=LL,
!      STREAM=WV,
!      CLASS=OD,
!      EXPVER=1,
!      TYPE=AN,
!      LEVTYPE=SFC,
!      GRID=1.0/1.0,
!      AREA=-59/-30/-59/-29,
!      LEVELIST=OFF


!     FOR DATA AFTER TO JUNE 28, 1998 18Z UNTIL NOVEMBER 20, 2000 12Z 

!    RETRIEVE,
!      DATE=19980704,
!      TIME=0/6/12/18,
!      FORMAT=PACKED,
!      TARGET="input_spectra",
!      PARAM=251,
!      DOMAIN=G,
!      REPRES=LL,
!      STREAM=WV,
!      CLASS=OD,
!      EXPVER=1,
!      TYPE=AN,
!      LEVTYPE=SFC,
!      GRID=1.0/1.0,
!      DIRECTION=1/to/12,
!      FREQUENCY=1/to/25,
!      AREA=-59/-30/-59/-29,
!      LEVELIST=OFF

!    SINCE NOVEMBER 20, 2000 18Z 

!    RETRIEVE,
!      DATE=20020521,
!      TIME=0/6/12/18,
!      FORMAT=PACKED,
!      TARGET="input_spectra",
!      PARAM=251,
!      DOMAIN=G,
!      REPRES=LL,
!      STREAM=WV,
!      CLASS=OD,
!      EXPVER=1,
!      TYPE=AN,
!      LEVTYPE=SFC,
!      GRID=1.0/1.0,
!      DIRECTION=1/to/24,
!      FREQUENCY=1/to/30,
!      AREA=-59/-30/-59/-29,
!      LEVELIST=OFF

!    !!! Since 26th January 2010 for the global model
!      DIRECTION=1/to/36,
!      FREQUENCY=1/to/36,

!    !!! if data come from the archive then one can also use
!      DIRECTION=all,
!      FREQUENCY=all,

! ----------------------------------------------------------------------

      USE GRIB_API

      IMPLICIT NONE

      REAL, PARAMETER :: ZMISS=-999.0
      REAL, PARAMETER :: PI=3.1415927
      REAL, PARAMETER :: ZPI=2.*PI

      INTEGER :: MORARG, IU06
      INTEGER :: I1, ITEST, LFILE
      INTEGER :: ITABPAR, IPARAM, ITABLE, IRGG, ISCAN
      INTEGER :: IYYYYMMDD, IHHMM, IHH
      INTEGER :: NGX, NGY, L, K, I, J, M, JSN, ISTART, IJ
      INTEGER :: NANG, NFRE, NFRANG, KK, MM
      INTEGER :: IDIRSCALING, IFRESCALING
      INTEGER :: ILON, ILAT
      INTEGER :: GETCLO, GETCLA, IOPTVAL
      INTEGER :: IFILE, IRET, IGRIB 
      INTEGER :: NUMBEROFVALUES, NBVAL
      INTEGER :: NGXLOC
      INTEGER :: PLPRESENT, NB_PL
      INTEGER, DIMENSION(:), ALLOCATABLE :: PL
      INTEGER, DIMENSION(:), ALLOCATABLE :: KDOMRGG

      REAL :: ONE,ZTHETA,ZFRE
      REAL :: DUM, ZDELLO
      REAL :: YFRST, YLAST, AMONOP, AMOSOP, AMOWEP, AMOEAP
      REAL :: XDELLO, XDELLA
      REAL :: DELLO
      REAL :: FCST, STEP
      REAL :: DELTH,PPREC,PPEPS,SUM
      REAL :: XLAT, XLON
      REAL(KIND=8), ALLOCATABLE :: VALUES(:)
      REAL, ALLOCATABLE :: SPEC(:) 
      REAL, ALLOCATABLE :: FR(:), TH(:), DFIM(:)
      REAL, ALLOCATABLE :: SCFR(:), SCTH(:)

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
      LOGICAL :: LPERIODIC

      DATA CLOPTS/'i;o;t;'/

! ----------------------------------------------------------------------

!*    INITIAL VALUES SET AND CRACK COMMAND LINE.
!     -----------------------------------------

      FNAMEOUT='output_spectra'
      ITEST=0
      IU06=6

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


!     GET FIRST DATA FILE

!     LOOP ON ALL MESSAGES IN INPUT FILE
 
      IGRIB=-99
      CALL  GRIB_NEW_FROM_FILE(IFILE,IGRIB,IRET)

      LOOP: DO WHILE (IRET /= GRIB_END_OF_FILE)

!*      DETERMINE DATA FIELD CHARACTERISTICS 

        CALL GRIB_GET(IGRIB,'paramId',ITABPAR)
        ITABLE=ITABPAR/1000
        IPARAM=ITABPAR-ITABLE*1000
        IF (ITEST.gt.0) WRITE(*,*) ' THE INPUT PARAMETER IS ',IPARAM 

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

        CALL GRIB_GET(IGRIB,'jDirectionIncrementInDegrees',XDELLA)
        CALL GRIB_GET(IGRIB,'iDirectionIncrementInDegrees',XDELLO)


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

        CALL ADJUST (AMOWEP, AMOEAP)
        LPERIODIC=.FALSE.
        DELLO=(AMOEAP-AMOWEP)/MAX(1,NGX-1)
        IF (AMOEAP-AMOWEP+1.5*DELLO.GE.360.) LPERIODIC=.TRUE.

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


        IF(IPARAM.EQ.250) THEN
!         obsolete !!!
          write(*,*) 'not ready yet for 250 !!!!'
          call abort
!!        probably the same
          CALL GRIB_GET(IGRIB,'numberOfDirections',NANG)
          CALL GRIB_GET(IGRIB,'numberOfFrequencies',NFRE)

        ELSE IF(IPARAM.EQ.251) THEN
          CALL GRIB_GET(IGRIB,'numberOfDirections',NANG)
          CALL GRIB_GET(IGRIB,'numberOfFrequencies',NFRE)
        ELSE
          WRITE(*,*) 'THE INPUT GRIB PARAMETER IS NOT 250 OR 251 BUT',
     &               IPARAM 
          WRITE(*,*) 'WHICH IS NOT A WAVE SPECTRUM PARAMETER !!!'
          WRITE(*,*) 'PROGRAM WILL ABORT'
          CALL ABORT
        ENDIF

        NFRANG=NANG*NFRE

        ALLOCATE(FR(NFRE))
        ALLOCATE(SCFR(NFRE))
        ALLOCATE(TH(NANG))
        ALLOCATE(SCTH(NANG))

        DELTH=ZPI/NANG

!       DECODE INPUT GRIB DATA
!       ----------------------

        IF(IPARAM.EQ.250) THEN

           write(*,*) ' not yet !!! '

        ELSE IF(IPARAM.EQ.251) THEN

!         GET THE SIZE OF THE VALUES ARRAY
          CALL GRIB_GET_SIZE(IGRIB,'values',NUMBEROFVALUES)
          ALLOCATE(VALUES(NUMBEROFVALUES))
          ALLOCATE(SPEC(NUMBEROFVALUES*NFRANG))

          DO M=1,NFRE
            DO K=1,NANG
            IF (ITEST.gt.0) 
     &          WRITE(*,*)' DECODING INPUT FIELD ',(M-1)*NANG+K 

!             SET THE MISSING DATA INDICATOR
              CALL GRIB_SET(IGRIB, 'missingValue',ZMISS)

!             CHECK ON NUMBER OF VALUES
              CALL GRIB_GET_SIZE(IGRIB,'values',NBVAL)
              IF(NBVAL.NE.NUMBEROFVALUES) THEN
                WRITE (*,*) '*************************************'
                WRITE (*,*) 'THE NUMBER OF VALUES SHOULD ALWAYS BE'
                WRITE (*,*) 'NUMBEROFVALUES = ',NUMBEROFVALUES
                WRITE (*,*) 'IT IS NOW ',NBVAL
                WRITE (*,*) 'PROBLEM !!!! ABORTING NOW !'
                WRITE (*,*) '*************************************'
                CALL ABORT
              ENDIF
!             GET DATA VALUES
              CALL GRIB_GET(IGRIB,'values',VALUES)

!             DETERMINE DATA FIELD CHARACTERISTICS 

              CALL GRIB_GET(IGRIB,'directionNumber',KK)
              CALL GRIB_GET(IGRIB,'directionScalingFactor',IDIRSCALING)
              CALL GRIB_GET(IGRIB,'scaledDirections',SCTH)
              TH(KK) = SCTH(KK)/IDIRSCALING

              CALL GRIB_GET(IGRIB,'frequencyNumber',MM)
              CALL GRIB_GET(IGRIB,'frequencyScalingFactor',IFRESCALING)
              CALL GRIB_GET(IGRIB,'scaledFrequencies',SCFR)
              FR(MM) =  SCFR(MM)/IFRESCALING

              IF(KK.NE.K.AND.MM.NE.M) THEN
                WRITE (*,*) '************************************'
                WRITE (*,*) '* DIRECTION AND FREQUENCY INDEX ARE'
                WRITE (*,*) '* OUT OF ORDER. '
                WRITE (*,*) '* IN FILE ',FNAMEIN
                WRITE (*,*) '* WE EXPECTED DIRECTION : ', K 
                WRITE (*,*) '* AND GOT ', KK
                WRITE (*,*) '* WE EXPECTED FREQUENCY : ', M 
                WRITE (*,*) '* AND GOT ', MM 
                WRITE (*,*) '* VERIFY THAT YOUR MARS REQUEST WAS'
                WRITE (*,*) '* DONE FOR ALL DIRECTIONS AND FREQUENCIES'
                WRITE (*,*) '* IN THE DEFAULT ORDER.'
                WRITE (*,*) '************************************'
                CALL ABORT
              ENDIF

              PPREC=0.0
              PPEPS=1.0e-10
              DO IJ=1,NUMBEROFVALUES
                IF (VALUES(IJ) .NE. ZMISS) THEN
!!                the use of variable DUM is to insure that we get same
!!                results as with gribex
                  DUM=VALUES(IJ)
                  SPEC(K+(M-1)*NANG+(IJ-1)*NFRANG) =
     &                      10.**(DUM-ABS(PPREC))- PPEPS
                ELSE
                  SPEC(K+(M-1)*NANG+(IJ-1)*NFRANG) = 0.0
                ENDIF
              ENDDO

!             GET NEXT FIELD

              IF(.NOT.(K.EQ.NANG.AND.M.EQ.NFRE)) THEN
                CALL GRIB_RELEASE(IGRIB)
                IGRIB=-99
                CALL GRIB_NEW_FROM_FILE(IFILE,IGRIB, IRET)
              ENDIF

            ENDDO
          ENDDO

        ENDIF
 
!       WRITES OUT THE SPECTRUM IN ASCCII FORMAT 
!       ****************************************

        LFILE=0
        IF (FNAMEOUT.NE. ' ') LFILE=LEN_TRIM(FNAMEOUT)
        OPEN(10,FILE=FNAMEOUT(1:LFILE))
        OPEN(11,FILE=FNAMEOUT(1:LFILE)//'.1d')
 
        dum=0.
        ilon=0
        ilat=1
        DO IJ=1,NUMBEROFVALUES
          ilon=ilon+1
          IF(IRGG.EQ.0) THEN
            if(ilon.gt.NGX) then
              ilon=1
              ilat=ilat+1
            endif
            zdello=xdello
          ELSE
            ngxloc = KDOMRGG(ILAT)
            if(ilon.gt.NGXLOC) then
              ilon=1
              ilat=ilat+1
              IF(LPERIODIC) THEN
                zdello=360./ngxloc
              ELSE
                zdello=(AMOEAP-AMOWEP)/(ngxloc-1)
              ENDIF
            endif
          ENDIF
          xlon=AMOWEP+(ilon-1)*zdello
          xlat=AMONOP-(ilat-1)*xdella

          WRITE(10,'(1X,F7.3,1X,F8.4,1x,i8.8,i4.4,1x,f4.0)')
     &          xlat,xlon,IYYYYMMDD,IHHMM,FCST 
          WRITE(11,'(1X,F7.3,1X,F8.4,1x,i8.8,i4.4,1x,f4.0)')
     &          xlat,xlon,IYYYYMMDD,IHHMM,FCST 

          WRITE(10,110) dum,(TH(K),K=1,NANG)
          DO M=1,NFRE
            WRITE(10,110) FR(M),
     &         (SPEC(K+(M-1)*NANG+(IJ-1)*NFRANG),K=1,NANG)
            SUM=0.
            DO K=1,NANG
               SUM=SUM+SPEC(K+(M-1)*NANG+(IJ-1)*NFRANG)
            ENDDO
            SUM=SUM*DELTH
            WRITE(11,*) FR(M),SUM
          ENDDO

        ENDDO

110   FORMAT(37(1X,F9.5))


        IF(ALLOCATED(KDOMRGG))DEALLOCATE(KDOMRGG)
        IF(ALLOCATED(TH))DEALLOCATE(TH)
        IF(ALLOCATED(SCTH))DEALLOCATE(SCTH)
        IF(ALLOCATED(SCFR))DEALLOCATE(SCFR)
        IF(ALLOCATED(VALUES)) DEALLOCATE(VALUES)
        IF(ALLOCATED(SPEC)) DEALLOCATE(SPEC)


        CALL GRIB_RELEASE(IGRIB)
        IGRIB=-99
        CALL GRIB_NEW_FROM_FILE(IFILE,IGRIB, IRET)

      END DO LOOP

      CALL GRIB_CLOSE_FILE(IFILE)

      CLOSE(10)
      CLOSE(11)

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

      SUBROUTINE ADJUST(WEST, EAST)

! ----------------------------------------------------------------------

!**** *ADJUST* - ROUTINE TO CORRECT BORDERS OF INTERVALS.

!     H.GUNTHER            ECMWF       04/04/1990

!*    PURPOSE.
!     -------

!       ADJUSTS INTERVAL BORDERS GIVEN IN DEGREE.

!**   INTERFACE.
!     ----------

!       *CALL* *ADJUST (WEST, EAST)*
!          *WEST*    - LEFT INTERVAL BORDER IN DEGREE.
!          *EAST*    - RIGHT INTERVAL BORDER IN DEGREE.

!     METHOD.
!     -------

!       THE INTERVAL BORDERS ARE CHANGED TO FULLFILL:
!         0. .LE. EAST  .AND. EAST .LT. 360. .AND. WEST .LE. EAST

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

!* 1. CORRECT BORDERS.
!     ----------------

      WEST = MOD(WEST+720.,360.)
      EAST = MOD(EAST+720.,360.)
      IF (WEST.GT.EAST) WEST = WEST-360.

      RETURN
      END SUBROUTINE ADJUST
