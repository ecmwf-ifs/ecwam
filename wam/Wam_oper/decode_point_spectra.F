!
! ----------------------------------------------------------------------
!
      PROGRAM DECODE_POINT_SPECTRA 
!
!     J BIDLOT ECMWF   OCTOBRE 1998
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
 
 
!     LIBRARY : FORTRAN90, and EMOSLIB
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
! ----------------------------------------------------------------------

      PARAMETER(NFREMAX=36)
      PARAMETER (NBIT = 50000)
      PARAMETER (G = 9.806, PI = 3.1415927, CIRC = 40000000.,
     1           ZPI = 2.*PI, RAD = PI/180., DEG = 180./PI,
     2           R = CIRC/ZPI, EPSMIN=0.1E-32, ZMISS=-999.0)

      INTEGER I4
      INTEGER GETCLO, GETCLA, IOPTVAL
      INTEGER :: KSEC0(2),KSEC1(255),KSEC2(2060),
     &           KSEC3(2),KSEC4(128)
      INTEGER :: ISEC0(2),ISEC1(255),ISEC2(2060),
     &           ISEC3(2),ISEC4(128)
      INTEGER, ALLOCATABLE :: INGRIB(:),OUTGRIB(:),KDOMRGG(:)
      REAL :: ONE,ZTHETA,ZFRE
      REAL :: PSEC2(96),PSEC3(2)
      REAL :: ZSEC2(96),ZSEC3(2)
      REAL, ALLOCATABLE :: PSEC4(:),PSEC4OUT(:)
      REAL, ALLOCATABLE :: FR(:), TH(:), DFIM(:)
      CHARACTER     CLOPTLET
      CHARACTER*  3 CLL1
      CHARACTER*  6 CLOPTS
      CHARACTER* 12 CLFMT
      CHARACTER*40 ERRMSG(-4:6), MSG
      CHARACTER*128 CLARG, FNAMEIN, FNAMEOUT
      LOGICAL LLEXIST, LINTEGRATE, LWAVEHGT, LLASTFRE

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
     5            ' NO 2D SPECTRA IN INPUT USE "INMARSB"  *',
     6            ' ENCODING ERROR SEE GRIBEX DESCRIPTION *'/

! ----------------------------------------------------------------------

!*    INITIAL VALUES SET AND CRACK COMMAND LINE.
!     -----------------------------------------

      I4=1
      NPRECI = KIND(I4) 
      FNAMEIN='input_spectra'
      FNAMEOUT='output_spectra'
      IFRE_FIRST=1
      IFRE_LAST=NFREMAX
      LINTEGRATE=.FALSE.
      LWAVEHGT=.FALSE.
      LLASTFRE=.FALSE.
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

      KSEC1=0
      KSEC2=0
      KSEC4=0
      KSEC3(1)=0
      KSEC3(2)=0
      PSEC2=0.
      PSEC3(2)=ZMISS

!*    INPUT FIRST GRIB DATA FILE
!     --------------------------
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
        WRITE(*,*)'*GRIB SPECTRA NOT FOUND IN *'
        WRITE(*,*)  FNAMEIN 
        WRITE(*,*)'*PROGRAM WILL ABORT        *'
        WRITE(*,*)'*                          *'
        WRITE(*,*)'****************************'
        CALL ABORT
      ENDIF
 
!     GET SIZE OF INGRIB
 
      CALL PBSIZE (IUGRS, IPLENG)
      ISIZE=(IPLENG+NPRECI-1)/NPRECI
      ISIZEMAX=ISIZE
      ALLOCATE(INGRIB(ISIZE))
      ISIZE=NBIT
1111  IPLENG=ISIZE*NPRECI
      IF(.NOT.ALLOCATED(INGRIB)) ALLOCATE(INGRIB(ISIZE))


!     GET FIRST DATA FILE

      CALL PBGRIB(IUGRS,INGRIB,IPLENG,ILENG,KRET)
      IF     (KRET.EQ.-1) THEN
        WRITE (*,*) ' REACHED EOF IN ',FNAMEIN 
        GOTO 3333 
      ELSEIF (KRET .EQ. -2) THEN
        WRITE (*,*) ' ERROR IN FILE HANDLING IN ',FNAMEIN 
      ELSEIF (KRET .EQ. -3) THEN
        DEALLOCATE(INGRIB)

        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' ***** WARNING IN '
        WRITE(IU06,*) ' SIZE  OF KGRIB IS NOT BIG ENOUGH.'
        WRITE(IU06,*) ' IT WAS ', ISIZE

        KKRET=0

        ISIZE=(ILENG+KIND(ISIZE)-1)/KIND(ISIZE)
        WRITE(IU06,*) ' IT SHOULD AT LEAST BE ', ISIZE
        WRITE(IU06,*) ' THE SIZE WAS RESET AUTOMATICALLY'
        WRITE(IU06,*) ' AND THE FIELD READ WITH THE NEW SIZE'
        WRITE(IU06,*) ' IF THIS PROBLEM OCCURS TOO OFTEN'
        WRITE(IU06,*) ' MODIFY THE VALUE OF NBIT IN SOURCE'
        WRITE(IU06,*) ' ***** WARNING ****** WARNING *****'
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
      ILENP = 1
      ALLOCATE(PSEC4(ILENP))
      IF (ITEST.gt.0)
     &     WRITE(*,*)' GETTING GRIB HEADER OF 1st INPUT FIELD' 
      IF (ITEST.GT.2) CALL GRSDBG (1)
      CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &             PSEC4, ILENP, INGRIB, ILENG , KWORD, 'J', KRET)
      IF (ITEST.GT.1) WRITE(*,*)' GRIBEX DONE status=' , KRET 
!     the error code 811: 'Cannot handle 2ndary bitmaps for J option'
!     is actually a warning message and should not cause a failure.
      IF(KRET.GT.0.and.kret.ne.811) THEN
        MSG = ERRMSG(2)
        CALL GRPRS0 (KSEC0)
        CALL GRPRS1 (KSEC0, KSEC1)
        WRITE (*,*) MSG
        CALL ABORT
      ENDIF
      IF(ALLOCATED(PSEC4)) DEALLOCATE(PSEC4)

!*    DETERMINE DATA FIELD CHARACTERISTICS 

      IPARAM = KSEC1(6)
      IRGG = KSEC2(17)
      AMONOP = FLOAT(KSEC2(4)/1000)+0.1*(MOD(KSEC2(4),1000)/100)+
     &         0.01*(MOD(KSEC2(4),100)/10)+0.001*MOD(KSEC2(4),10)
      AMOSOP = FLOAT(KSEC2(7)/1000)+0.1*(MOD(KSEC2(7),1000)/100)+
     &         0.01*(MOD(KSEC2(7),100)/10)+0.001*MOD(KSEC2(7),10)
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
      NGYN = NGY 
      IF (ITEST.gt.0) WRITE(*,*) ' THE INPUT PARAMETER IS ',IPARAM 
      IF(IPARAM.EQ.250) THEN
        IF(IRGG.EQ.0) THEN
          NGYMAX = NGY
        ELSE
          IF(AMONOP.GT.0..AND.AMOSOP.LE.0.) THEN
            NGYMAX = 2*KSEC2(3)-1
          ELSE IF(AMONOP.GT.0..AND.AMOSOP.GT.0.) THEN
            NGYMAX = NGY 
          ELSE
            WRITE (*,*) 'THE GRIB DATA SHOULD BE INPUT SUCH THAT' 
            WRITE (*,*) 'THE NORTHERN HEMISPHERE WAS FIRST RETRIEVED'
            WRITE (*,*) 'AND THEN THE SOUTHERN PART.'
            WRITE (*,*) 'THE PROGRAM WILL ABORT'
            CALL ABORT
          ENDIF
        ENDIF
        NANG = KSEC4(53) 
        NFRE = KSEC4(55) 
      ELSE IF(IPARAM.EQ.251) THEN
        NGYMAX = NGY
        NANG = KSEC1(46)
        NFRE = KSEC1(47)
      ELSE
        WRITE(*,*) 'THE INPUT GRIB PARAMETER IS NOT 250 OR 251 BUT',
     &               IPARAM 
        WRITE(*,*) 'WHICH IS NOT A WAVE SPECTRUM PARAMETER !!!'
        WRITE(*,*) 'PROGRAM WILL ABORT'
        CALL ABORT
      ENDIF

      ALLOCATE(KDOMRGG(NGYMAX))
      KDOMRGG=0
      NFRANG=NANG*NFRE
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

      ALLOCATE(FR(NFRE))
      ALLOCATE(TH(NANG))

      DELTH=ZPI/NANG

!     DECODE INPUT GRIB DATA
!     ----------------------

      IF(IPARAM.EQ.250) THEN

        IF(IRGG.EQ.0) THEN
          ILENP = NGY*NGX*NFRANG
          ILEN1 = NGY*NGX
        ELSE
          ILENP=0
          DO J=1,NGY-1
            ILENP = ILENP + KDOMRGG(J) 
          ENDDO
          ILENPS = ILENP
          IOFF = ILENP + NGX
          ILENP = 2*ILENP
          ILENP = ILENP + NGX
          ILEN1 = ILENP
          ILENP = ILENP*NFRANG
          ILENPS = ILENPS*NFRANG
          IOFF = IOFF*NFRANG
        ENDIF

        ALLOCATE(PSEC4(ILENP))
        PSEC4=ZMISS
        KSEC3(2)=0
        PSEC3(2)=ZMISS

        IF (ITEST.gt.0) WRITE(*,*)' DECODING INPUT FIELD' 
        CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &               PSEC4, ILENP, INGRIB, ILENG , KWORD, 'D', KRET)
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

        DO IC=1,NFRE
          FR(IC) = TRANSFER(KSEC4(59+NANG+IC),ONE) 
        ENDDO
        DO IC=1,NANG
          TH(IC) = TRANSFER(KSEC4(59+IC),ONE) 
        ENDDO

        IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)

        PPREC=0.0
        PPEPS=1.0e-10
        DO IJ=1,ILENP
          IF (PSEC4(IJ) .NE. ZMISS) THEN
            PSEC4(IJ) = 10.**(PSEC4(IJ)-ABS(PPREC))- PPEPS
          ELSE
            PSEC4(IJ) = 0.
          ENDIF
        ENDDO

      ELSE IF(IPARAM.EQ.251) THEN

        IF(IRGG.EQ.0) THEN
          ILENP = NGY*NGX*NFRANG
          ILEN1 = NGY*NGX
        ELSE
          ILENP=0
          DO J=1,NGY
            ILENP = ILENP + KDOMRGG(J) 
          ENDDO
          ILEN1 = ILENP
          ILENP = ILENP*NFRANG
        ENDIF

        ALLOCATE(PSEC4OUT(ILEN1))
        PSEC4OUT=ZMISS

        ALLOCATE(PSEC4(ILENP))

        DO M=1,NFRE
          DO K=1,NANG

            KSEC3(2)=0
            PSEC3(2)=ZMISS
            IF (ITEST.gt.0) 
     &          WRITE(*,*)' DECODING INPUT FIELD ',(M-1)*NANG+K 
            CALL GRIBEX(KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &               PSEC4OUT, ILEN1, INGRIB, ILENG , KWORD, 'D', KRET)
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

            IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)

            KK = KSEC1(44)
            MM = KSEC1(45)
            FR(MM) =  FLOAT(KSEC1(49+NANG+MM))/KSEC1(49)
            TH(KK) = FLOAT(KSEC1(49+KK))/KSEC1(48)
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
            DO IJ=1,ILEN1
              IF (PSEC4OUT(IJ) .NE. ZMISS) THEN
                PSEC4OUT(IJ) = 10.**(PSEC4OUT(IJ)-ABS(PPREC))- PPEPS
              ELSE
                PSEC4OUT(IJ) = 0.0
              ENDIF
            ENDDO

            DO IJ=1,ILEN1
              PSEC4(K+(M-1)*NANG+(IJ-1)*NFRANG) =  PSEC4OUT(IJ)
            ENDDO

!           GET NEXT FIELD

            IF(.NOT.(K.EQ.NANG.AND.M.EQ.NFRE)) THEN
2222          IPLENG=ISIZE*NPRECI
              IF(.NOT.ALLOCATED(INGRIB)) ALLOCATE(INGRIB(ISIZE))

              IF (ITEST.gt.1) WRITE (*,*) ' CALLING PBGRIB ' 
              CALL PBGRIB(IUGRS,INGRIB,IPLENG,ILENG,KRET)
              IF     (KRET.EQ.-1) THEN
                WRITE (*,*) ' REACHED EOF IN ',FNAMEIN 
                WRITE (*,*) ' CHECK IF YOU HAVE PROVIDED THE INPUT' 
                WRITE (*,*) ' SPECTRA FOR ALL DIRECTIONS ie: ',NANG
                WRITE (*,*) ' AND FOR ALL FREQUENCIES ie: ',NFRE 
              ELSEIF (KRET .EQ. -2) THEN
                WRITE (*,*) ' ERROR IN FILE HANDLING IN ',FNAMEIN 
              ELSEIF (KRET .EQ. -3) THEN
                DEALLOCATE(INGRIB)

                WRITE(IU06,*) ' '
                WRITE(IU06,*) ' ***** WARNING IN '
                WRITE(IU06,*) ' SIZE  OF KGRIB IS NOT BIG ENOUGH.'
                WRITE(IU06,*) ' IT WAS ', ISIZE

                KKRET=0

                ISIZE=(ILENG+KIND(ISIZE)-1)/KIND(ISIZE)
                WRITE(IU06,*) ' IT SHOULD AT LEAST BE ', ISIZE
                WRITE(IU06,*) ' THE SIZE WAS RESET AUTOMATICALLY'
                WRITE(IU06,*) ' AND THE FIELD READ WITH THE NEW SIZE'
                WRITE(IU06,*) ' IF THIS PROBLEM OCCURS TOO OFTEN'
                WRITE(IU06,*) ' MODIFY THE VALUE OF NBIT IN SOURCE'
                WRITE(IU06,*) ' ***** WARNING ****** WARNING *****'
                WRITE(IU06,*) ' '
                CALL FLUSH(IU06)

!               RESET THE FILE POINTER TO READ FIELD AGAIN
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

                GOTO 2222

              ENDIF
            ENDIF

          ENDDO
        ENDDO

      ENDIF

!     WRITES OUT THE SPECTRUM IN ASCCII FORMAT 
!     ****************************************

      LFILE=0
      IF (FNAMEOUT.NE. ' ') LFILE=LEN_TRIM(FNAMEOUT)
      OPEN(10,FILE=FNAMEOUT(1:LFILE))
      OPEN(11,FILE=FNAMEOUT(1:LFILE)//'.1d')
 
      dum=0.
      ilon=0
      ilat=1
      DO IJ=1,ILEN1
        ilon=ilon+1
        IF(IRGG.EQ.0) THEN
          if(ilon.gt.NGX) then
            ilon=1
            ilat=ilat+1
          endif
          zdello=xdello
        ELSE
          ngxloc = KSEC2(22+ilat+ISTART)
          if(ilon.gt.NGXLOC) then
            ilon=1
            ilat=ilat+1
            zdello=(AMOEAP-AMOWEP)/(ngxloc-1)
          endif
        ENDIF
        xlon=AMOWEP+(ilon-1)*zdello
        xlat=AMONOP-(ilat-1)*xdella

        WRITE(10,'(2(1x,F8.4),1x,i4.4,4(i2.2),1x,f4.0)')
     &        xlat,xlon,IYYYY,IMM,IDD,IHH,IMI,FCST 
        WRITE(11,'(2(1x,F8.4),1x,i4.4,4(i2.2)1x,f4.0)')
     &        xlat,xlon,IYYYY,IMM,IDD,IHH,IMI,FCST 

        IF(NANG.EQ.12) THEN
          WRITE(10,100) dum,(TH(K),K=1,NANG)
          DO M=1,NFRE
            WRITE(10,100) FR(M),
     &         (PSEC4(K+(M-1)*NANG+(IJ-1)*NFRANG),K=1,NANG)
            SUM=0.
            DO K=1,NANG
               SUM=SUM+PSEC4(K+(M-1)*NANG+(IJ-1)*NFRANG)
            ENDDO
            SUM=SUM*DELTH
            WRITE(11,*) FR(M),SUM
          ENDDO
        ELSE IF(NANG.EQ.24) THEN
          WRITE(10,110) FR(M),(TH(K),K=1,NANG)
          DO M=1,NFRE
            WRITE(10,110) FR(M),
     &         (PSEC4(K+(M-1)*NANG+(IJ-1)*NFRANG),K=1,NANG)
            SUM=0.
            DO K=1,NANG
               SUM=SUM+PSEC4(K+(M-1)*NANG+(IJ-1)*NFRANG)
            ENDDO
            SUM=SUM*DELTH
            WRITE(11,*) FR(M),SUM
          ENDDO
        ELSE
          WRITE(*,*) 'NO FORMAT WRITE POSSIBLE'
          WRITE(*,*) 'THE NUMBER OF DIRECTIONS IS NOT 12 OR 24'
        ENDIF
      ENDDO

100   FORMAT(13(1X,F8.4))
110   FORMAT(25(1X,F8.4))


      IF(ALLOCATED(KDOMRGG))DEALLOCATE(KDOMRGG)
      IF(ALLOCATED(TH))DEALLOCATE(TH)
      IF(ALLOCATED(FR))DEALLOCATE(FR)
      IF(ALLOCATED(PSEC4)) DEALLOCATE(PSEC4)
      IF(ALLOCATED(PSEC4OUT)) DEALLOCATE(PSEC4OUT)
      GOTO 1111


3333  CALL PBCLOSE(IUGRS,KRET)
      IF(KRET.LT.0) THEN
        WRITE (*,*) '************************************'
        WRITE (*,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
        WRITE (*,*) '* FILE ',FNAMEIN
        WRITE (*,*) '************************************'
        CALL ABORT
      ENDIF

      CLOSE(10)
      CLOSE(11)

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
