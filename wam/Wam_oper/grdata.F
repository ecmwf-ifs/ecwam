      SUBROUTINE GRDATA(NWH, NSTART, NEND, LAVERAGE)

!--------------------------------------------------------------------

!**** *GRDATA* - TRANSFERS MEASUREMENTS TO MODEL GLOBAL GRID,
!****            DISCARDING UNRELIABLE DATA.

!     P.LIONELLO     ECMWF       APRIL 1990

!     PURPOSE.
!     --------

!       PREPARE THE DATA FOR OPTIMAL INTERPOLATION.

!**   INTERFACE.
!     ----------

!       *CALL* *GRDATA(NWH NSTART, NEND)*
!        *NWH*  NUMBER OF RELEVANT GRIDDED OBSERVATIONS.
!        *NSTART* INTEGER INDEX OF THE FIRST POINT OF SUB GRID DOMAIN
!        *NEND*   INTEGER INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!        *LAVERAGE* LOGICAL TRUE IF SUPER-OBSERVATIONS ARE TO BE 
!                                   ASSIMILATED

!     METHOD.
!     -------

!        A BOX OF SIZE DELLA*DELLO IS CENTERED IN EACH GRID POINT.
!        THE AVERAGE VALUE OF THE SATELLITE MEASUREMENTS INSIDE THE
!        BOX IS TAKEN TO UPDATE THE WAM MODEL. THE VALUE IS DISCARDED
!        IF THE ROOT MEAN SQUARE ERROR IS TOO BIG.
!        MOREOVER THE VALUE IS DISCARDED IF THERE ARE LESS THEN FOUR
!        MEASUREMENTS INSIDE THE BOX.

!      EXTERNALS.
!      ----------

!         READSAT  - READ A MEASUREMENT.

!      MODIFICATONS.
!      -------------
!         FEB 92  Bjoern Hansen :
!                 Interface to readsat extended in order to transfer
!                 the flags added to the satellite information by the
!                 quality control preura.
!                 The as unreliable flagged satellite measurements
!                 are rejected.
!                 The input unit of the measurements is automatically
!                 looked up each time grdata is called. It is assumed
!                 that the needed input file is present in the working
!                 directory, contains data of the same period of time
!                 as the output time step and is named RFLyyyymmddhhmm,
!                 where yyyymmddhhmm is a DATE/TIME group.

!          Aug 97 Jean Bidlot : in coupled mode IDELWO is not the
!                 appropriate time span for the selection of the data
!                 to be used but rather IDELALT, A NEW VARIABLE WHICH
!                 IS SUPPLIED IN WAMINPUT

!          Oct 97 Jean Bidlot : use of the new ice mask instead of 
!                 SST field to determine whether an observation which
!                 falls on a model sea point.

!          May 01 Jean Bidlot : do not use wind speed observations

!                 Jean Bidlot : only use grid information to thin
!                               and smooth the data.

!          J. BIDLOT       *ECMWF*       SEP 2001
!            INTRODUCE THE CONCEPT OF A BLACKLIST TO REMOVE DATA FROM
!            BEING PROCESSED. IT IS TRIGGERED IF THE FILE
!            wave_data_blacklist IS PRESENT.
!            THE BLACKLIST WILL CONTAIN THE FOLLOWING INFORMATION
!            STATID : THE SATELLITE IDENTIFICATION NUMBER (AS IN BUFR)
!            SENSOR : THE INSTRUMENT IDENTIFICATION NUMBER (AS IN BUFR)
!            THE DATA WILL BE BLACKLISTED BETWEEN THE DATE CDATESTART
!            AND CDATEEND AND BETWEEN LATMI<= LATITUDE <= LATMA
!            AND BETWEEN LONGMI <= LONGITUDE <= LONGMA
!            FOR the OBSERVED PARAMETER PRESCRIBED BY PARAM (AS IN BUFR)
!            PROVIDED IT IS BETWEEN THE VALUES VALUE_MIN AND VALUE_MAX

!-----------------------------------------------------------------------

      USE YOWALTAS , ONLY : NUMALT, NIJALT   ,NALTDT   ,IJALT    ,
     &              ALTDATA,IBUFRSAT 
      USE YOWGRID  , ONLY : NLONRGG  ,IGL      ,IJS      ,IJL
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,IPER     ,AMOWEP   ,
     &            AMOSOP   ,XDELLA   ,ZDELLO
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : NPROC
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWSTAT  , ONLY : CDTPRO   ,IDELALT
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------

      LOGICAL LAVERAGE

      INTEGER, DIMENSION(NUMALT) :: ICOUNT, ICOUNTT, ICOUNTD,
     &                              ICOUNTB, ICOUNTQCF,
     &                              ICOUNTBL, IOLD, JOLD
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: NUMBWH, NR
      REAL, DIMENSION(NUMALT) :: SWHO
      REAL, ALLOCATABLE, DIMENSION(:,:,:) ::  WHME, WHSE
      CHARACTER*12 CBEGINDT, CENDDT
      CHARACTER*14 CDATE
      CHARACTER*14 :: CDATEO(NUMALT)
      LOGICAL EOFD

!*     VARIABLE     TYPE     PURPOSE.
!      --------     ----     --------

!     *CBEGINDT*    CHAR*12  BEGINNING DATE OF THE TIME WINDOW
!                            DURING WHICH DATA ARE USED FOR ANALYSIS.
!     *CENDDT*      CHAR*12  END DATE OF THE TIME WINDOW
!                            DURING WHICH DATA ARE SAVE FOR ANALYSIS.
!     *WHME*        REAL     GRID WAVE HEIGHT.
!     *WHSE*        REAL     WAVE HEIGHT ROOT MEAN SQUARE ERROR.
!     *NUMBWH*      INTEGER  NUMBER OF HS MEASUREMENTS INSIDE A BOX.

! ----------------------------------------------------------------------

      INTEGER,DIMENSION(NPROC) :: NSTART, NEND
      INTEGER,ALLOCATABLE,DIMENSION(:) :: STATID, SENSOR, PARAM

      REAL,ALLOCATABLE,DIMENSION(:) :: LATMI, LONGMI, LATMA, LONGMA
      REAL,ALLOCATABLE,DIMENSION(:) :: VALUE_MIN, VALUE_MAX, DLON

      CHARACTER*19 YOFILE
      CHARACTER*100 CDUM
      CHARACTER (LEN=12),ALLOCATABLE :: CDATESTART(:), CDATEEND(:)

      LOGICAL :: LLBLEXIST
      LOGICAL :: OPENED

      PARAMETER (MAXSATBC=10)
      INTEGER, DIMENSION(MAXSATBC) :: ISATBC, NBCT
      REAL,ALLOCATABLE,DIMENSION(:,:) :: BCINC
      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: BCTABLE
      CHARACTER (LEN=12),ALLOCATABLE,DIMENSION(:,:) :: CBCDATESTART
      CHARACTER (LEN=12),ALLOCATABLE,DIMENSION(:,:) :: CBCDATEEND
      LOGICAL :: LLBCTEXIST, LLDOBC, LNEW

      CHARACTER* 3 YOFID

! ----------------------------------------------------------------------


!*    1. FIX TIME WINDOW EQUAL TO ONE WIND TIME STEP.
!        --------------------------------------------

      CBEGINDT = CDTPRO
      CALL INCDATE (CBEGINDT,-IDELALT/2)
      CENDDT = CBEGINDT
      CALL INCDATE (CENDDT,IDELALT)

      WRITE(IU06,*) ' GRDATA : SELECTING ALTIMETER MEASUREMENTS '
      WRITE(IU06,*) '          PERIOD FROM ', CBEGINDT,' TO ', CENDDT

!     FIX THE MINIMUM NUMBER OF OBSERVATIONS PER GRID BOX REQUIRED
!     IN ORDER TO ACCEPT THE DATA AT THE GRIB POINTS

      IF (LAVERAGE) THEN
!       if superobservations fall within the same grid box, 
!       they will be averaged.
        NNMIN=1
      ELSEIF(XDELLA.GE.0.5) THEN
        NNMIN=4
      ELSEIF(XDELLA.GT.0.25) THEN
        NNMIN=3
      ELSE
        NNMIN=2
      ENDIF
      WRITE(IU06,*) '   GRIDDED DATA WILL BE KEPT WHEN THE NUMBER OF'
      WRITE(IU06,*) '   OBSERVATIONS PER GRID BOX IS LARGER THAN ',
     &  NNMIN-1

      IF (LAVERAGE) THEN
        YOFID='AFA'
      ELSE 
        YOFID='QCP'
      ENDIF
! ----------------------------------------------------------------------

!*    2. INITIALIZE ARRAYS, counters and units.
!        --------------------------------------

      ALLOCATE(WHME(NGX,NGY,NUMALT))
      ALLOCATE(WHSE(NGX,NGY,NUMALT))
      ALLOCATE(NUMBWH(NGX,NGY,NUMALT))
      ALLOCATE(NR(NGX,NGY,NUMALT))

      DO ISAT=1,NUMALT
        DO J=1,NGY
          DO I=1,NGX
            WHME(I,J,ISAT) = 0.
            WHSE(I,J,ISAT) = 0.
            NUMBWH(I,J,ISAT) = 0
            NR(I,J,ISAT)=0
          ENDDO
        ENDDO
        ICOUNT(ISAT)=0
        ICOUNTT(ISAT)=0
        ICOUNTD(ISAT)=0
        ICOUNTB(ISAT)=0
        ICOUNTQCF(ISAT)=0
        ICOUNTBL(ISAT)=0
        IOLD(ISAT) = 0
        JOLD(ISAT) = 0
        IBUFRSAT(ISAT)=0 
        CDATEO(ISAT) = CBEGINDT(1:12)//'00'
      ENDDO


      DO IUME=33,99,1
        INQUIRE ( UNIT=IUME, OPENED=opened)
        IF ( .NOT. opened) THEN
          WRITE(iu06,*)
     &     ' GRDATA : TAKES UNIT ',iume,' AS ALTIMETER DATA INPUT UNIT'
          GOTO 2000
        ENDIF
      ENDDO

!     CHECK IF BLACKLISTING FILE IS PRESENT
2000  CONTINUE
      NBL=0
      YOFILE='wave_data_blacklist'
      INQUIRE (FILE=YOFILE,EXIST=LLBLEXIST)
      IF(LLBLEXIST) THEN
        IOBL = I_GET_UNIT(IU06, YOFILE, 'r', 'f', 0)
        WRITE(IU06,*) '   GRDATA : READING BLACKLISTING FILE ', YOFILE
        WRITE(IU06,*) ' '

        DO IC=1,21
          READ(IOBL,'(A100)',END=200,ERR=200) CDUM
        ENDDO

        DO WHILE(.TRUE.)
          READ(IOBL,'(A100)',END=100,ERR=100) CDUM
          NBL=NBL+1
        ENDDO
100     WRITE(IU06,*) '            TOTAL NUMBER OF ENTRIES :',NBL
        REWIND(IOBL)

        IF(NBL.GT.0) THEN
          ALLOCATE(STATID(NBL))
          ALLOCATE(SENSOR(NBL))
          ALLOCATE(CDATESTART(NBL))
          ALLOCATE(CDATEEND(NBL))
          ALLOCATE(LATMI(NBL))
          ALLOCATE(LONGMI(NBL))
          ALLOCATE(LATMA(NBL))
          ALLOCATE(LONGMA(NBL))
          ALLOCATE(DLON(NBL))
          ALLOCATE(PARAM(NBL))
          ALLOCATE(VALUE_MIN(NBL))
          ALLOCATE(VALUE_MAX(NBL))
        ENDIF

        DO IC=1,21
          READ(IOBL,'(A100)',END=200,ERR=200) CDUM
        ENDDO

        DO IBL=1,NBL
          READ(IOBL,111,ERR=200) STATID(IBL), SENSOR(IBL),
     &         CDATESTART(IBL), CDATEEND(IBL), LATMI(IBL), LONGMI(IBL),
     &         LATMA(IBL), LONGMA(IBL), PARAM(IBL), VALUE_MIN(IBL),
     &         VALUE_MAX(IBL)
        ENDDO
        GOTO 300
200     CONTINUE
        WRITE(IU06,*) '************************************************'
        WRITE(IU06,*) '   GRDATA : ERROR READING ', YOFILE
        WRITE(IU06,*) '   THE PROGRAM WILL ABORT !!!!'
        WRITE(IU06,*) '************************************************'
        CALL ABORT

300     CLOSE(IOBL)
        DO IBL=1,NBL
          IF(ABS(LONGMI(IBL)-LONGMA(IBL)) .NE. 360.) THEN
            LONGMI(IBL) = MOD(LONGMI(IBL)+720.,360.)
            LONGMA(IBL) = MOD(LONGMA(IBL)+720.,360.)
          ENDIF
          IF(LONGMI(IBL).GE.LONGMA(IBL)) LONGMI(IBL)=LONGMI(IBL)-360.
          DLON(IBL) =  LONGMA(IBL) - LONGMI(IBL)

        ENDDO
      ELSE
        WRITE(IU06,*) '************************************************'
        WRITE(IU06,*) '   GRDATA : THE BLACK LISTING FILE ',YOFILE
        WRITE(IU06,*) '   WAS NOT FOUND !!!'
        WRITE(IU06,*) '   THE PROCESSING WILL CONTINUE WITHOUT IT !!!!'
        WRITE(IU06,*) '************************************************'
      ENDIF
111   FORMAT(2(I8,1x),2(A12,1x),2(f5.1,1x,f6.1,1x),I6,1x,2(f9.3,1x))

! ----------------------------------------------------------------------
!     CHECK IF THE FILE WITH "BIAS CORRECTION TABLES" IS PRESENT
2500  CONTINUE
      NBCT=0
      YOFILE='biascorrection.swh'
      INQUIRE (FILE=YOFILE,EXIST=LLBCTEXIST)
      IF(LLBCTEXIST) THEN
        IOBCT = I_GET_UNIT(IU06, YOFILE, 'r', 'f', 0)
        WRITE(IU06,*) '   GRDATA : READING FILE WITH BIAS CORRECTION '
     &                    //'TABLES  ', YOFILE

        WRITE(IU06,*) ' '

        DO IC=1,21
          READ(IOBCT,'(A)',END=2800,ERR=2800) CDUM
        ENDDO

        NSATBC=0
        MAXNBCINC=0
        DO WHILE(.TRUE.)
          READ(IOBCT,'(I4,1X,I4,A27,I4,1X,F5.2)',END=2600,ERR=2600)
     &             ISAT, IDUM, CDUM, NBCINC, BCINC0
          MAXNBCINC=MAX(MAXNBCINC,NBCINC)
          LNEW=.TRUE.
          DO IBCT=1, NSATBC
            IF (ISAT.EQ.ISATBC(IBCT)) THEN
              LNEW=.FALSE.
              NBCT(IBCT) = NBCT(IBCT) + 1
              EXIT
            ENDIF
          ENDDO
          IF (LNEW) THEN
            NSATBC=NSATBC+1
            IF (NSATBC.GT.MAXSATBC) THEN
              WRITE(IU06,*) '******************************************'
              WRITE(IU06,*) '   GRDATA : ERROR:  NUMBER OF SATELLITES '
     &             //'IN BIAS CORRECTION FILE EXCEEDS MAXSATBC =', 
     &                                                MAXSATBC 
              WRITE(IU06,*) '       INCREASE MAXSATBC IN GRDATA'
              WRITE(IU06,*) '   THE PROGRAM WILL ABORT !!!!'
              WRITE(IU06,*) '******************************************'
              CALL ABORT
            ENDIF
            ISATBC(NSATBC)=ISAT
            NBCT(NSATBC) = 1
          ENDIF
          READ(IOBCT,'(A)',END=2800,ERR=2800) CDUM
        ENDDO
2600    WRITE(IU06,*) '            NUMBER OF SATELLITE ENTRIES :',NSATBC
        MAXBCT=NBCT(1)
        DO IBCT=1, NSATBC
          WRITE(IU06,*) '              SATELLITE ',ISATBC(IBCT),
     &                                  ' HAS ', NBCT(IBCT), ' TABLES'
          MAXBCT=MAX(MAXBCT,NBCT(IBCT))
        ENDDO
        REWIND(IOBCT)

        IF(NSATBC.GT.0) THEN
          ALLOCATE(BCINC(NSATBC,MAXBCT))
          ALLOCATE(BCTABLE(NSATBC,MAXBCT,MAXNBCINC))
          ALLOCATE(CBCDATESTART(NSATBC,MAXBCT))
          ALLOCATE(CBCDATEEND(NSATBC,MAXBCT))
        ENDIF

        DO IC=1,21
          READ(IOBCT,'(A)',END=2800,ERR=2800) CDUM
        ENDDO

        DO IBCT=1, NSATBC
          NBCT(IBCT) = 0
        ENDDO
        DO WHILE(.TRUE.)
          READ(IOBCT,'(I4,1X,I4,A27,I4,1X,F5.2)',END=2900,ERR=2900)
     &          ISAT, IDUM, CDUM, NBCINC, BCINC0
          DO IBCT=1, NSATBC
            IF (ISAT.EQ.ISATBC(IBCT)) THEN
              NBCT(IBCT) = NBCT(IBCT) + 1
              BCINC(IBCT,NBCT(IBCT))=BCINC0
              CBCDATESTART(IBCT,NBCT(IBCT))=CDUM(2:13)
              CBCDATEEND  (IBCT,NBCT(IBCT))=CDUM(15:26)
              READ(IOBCT,'(9999F5.2)',END=2800,ERR=2800)
     &             (BCTABLE(IBCT,NBCT(IBCT),I),I=1,NBCINC)
              IF (ITEST.GE.1) THEN
                WRITE(IU06,'(A,I4,A)') 'BIAS CORRECTION TABLE FOR '//
     &             'SATELLITE ',ISAT, ' BETWEEN '//
     &             CBCDATESTART(IBCT,NBCT(IBCT))//' AND '//
     &             CBCDATEEND(IBCT,NBCT(IBCT))//' IS:'
                WRITE(IU06,'(5X,A)') 'WAVE HEIGHT     BIAS CORRECTION'
                DO I=1,NBCINC
                  WRITE(IU06,'(5X,F5.2,11X,F5.2)') (I-1.)*BCINC0, 
     &                                    BCTABLE(IBCT,NBCT(IBCT),I)
                ENDDO
                WRITE(IU06,'(A/)') 'END OF TABLE'
              ENDIF
              EXIT
            ENDIF
          ENDDO
        ENDDO

2800    CONTINUE
        WRITE(IU06,*) '************************************************'
        WRITE(IU06,*) '   GRDATA : ERROR READING ', YOFILE
        WRITE(IU06,*) '   THE PROGRAM WILL ABORT !!!!'
        WRITE(IU06,*) '************************************************'
        CALL ABORT

2900    CONTINUE
        CLOSE (IOBCT)
        WRITE(IU06,*) '            BIAS CORRECTION TABLES READ IN'

      ELSE
        WRITE(IU06,*) '************************************************'
        WRITE(IU06,*) '   GRDATA : FILE WITH BIAS CORRECTION TABLES ',
     &                                                           YOFILE
        WRITE(IU06,*) '   WAS NOT FOUND !!!'
        WRITE(IU06,*) '   THE PROCESSING WILL CONTINUE WITHOUT IT !!!!'
        WRITE(IU06,*) '************************************************'
      ENDIF
! ----------------------------------------------------------------------

!*    3. READ THE MEASUREMENT TAPE AND CUMULATE THE DATA ON GRID POINTS.
!        ---------------------------------------------------------------

 3000 CONTINUE
      CALL READSAT (IU06, IUME, CDTPRO, IDENTI, ISENSOR,
     &              CDATE, RLAT, RLON, IDES_HS, SWH,
     &              IDES_WS, WS, EOFD, YOFID)

      IF(ITEST.GE.2)
     & WRITE(IU06,'(A,A,2I6,2(I6,F6.2))') 
     &  'READSAT READ: ', CDATE,IDENTI,ISENSOR,
     &  IDES_WS, WS, IDES_HS, SWH
!     which satellite is it
      DO ISAT=1,NUMALT
        IF(IBUFRSAT(ISAT).EQ.0) THEN
          IF(ITEST.GE.1)
     &      WRITE(IU06,*) 'NEW SATELLITE DETECTED WITH ID ', IDENTI
          IBUFRSAT(ISAT)=IDENTI
          KSAT=ISAT
          EXIT
        ENDIF
        IF(IBUFRSAT(ISAT).EQ.IDENTI) THEN
          KSAT=ISAT
          EXIT
        ENDIF
      ENDDO
      IF(ISAT.GT.NUMALT) THEN
        WRITE(IU06,*) ' GRDATA : WARNING MORE SATELLITE IN INPUT'
        WRITE(IU06,*) '          THAN ALLOWED BY NUMALT !!!! '
        WRITE(IU06,*) '          THE SATELLITE IS ',IDENTI
        WRITE(IU06,*) '          NEW DATA WILL NOT BE ALLOWED IN ' 
        GOTO 3000
      ENDIF

      IF (EOFD) GOTO 4000
      ICOUNTD(KSAT) = ICOUNTD(KSAT)+1

!     CHECK DATE.

      IF (CDATE(1:12).LT.CBEGINDT)  then
        ICOUNTB(KSAT) = ICOUNTB(KSAT) + 1
        GOTO 3000
      ENDIF

      IF (CDATE(1:12).LE.CENDDT) THEN

!       BLACKLISTING

        IF(LLBLEXIST) THEN
          DO IBL=1,NBL
            IF(STATID(IBL).EQ.IDENTI .AND.
     &         SENSOR(IBL).EQ.ISENSOR .AND.
     &         CDATESTART(IBL).LE.CDATE(1:12) .AND. 
     &         CDATEEND(IBL).GE.CDATE(1:12) .AND. 
!           BLACKLISTING OF WAVE HEIGHT ONLY !!!
     &         PARAM(IBL).EQ.IDES_HS .AND. 
     &         LATMI(IBL).LE.RLAT .AND.
     &         LATMA(IBL).GE.RLAT .AND.
     &         VALUE_MIN(IBL).LE.SWH .AND.
     &         VALUE_MAX(IBL).GE.SWH
     &        ) THEN
              DX = MOD(RLON-LONGMI(IBL)+720.,360.)
              IF(DX.LE.DLON(IBL)) THEN
                ICOUNTBL(KSAT)=ICOUNTBL(KSAT)+1
                GOTO 3000
              ENDIF
            ENDIF
          ENDDO
        ENDIF

!       DATE IS INSIDE THE TIME WINDOW AND RECORD FLAGGED AS RELIABLE.
!       THE INDICES ARE COMPUTED

        ICOUNTT(KSAT) = ICOUNTT(KSAT)+1
        J       = NGY-NINT((RLAT-AMOSOP)/XDELLA)
        XI      = MOD(RLON-AMOWEP+720.,360.)

        JSN     = NGY-J+1
        IF((J.LE.0).OR.(J.GT.NGY)) JSN=1
        NRGG    = NLONRGG(JSN)
        I       = NINT(XI/ZDELLO(JSN)) + 1
        IF (IPER.EQ.1.AND.I.EQ.NRGG+1) I = 1

!       IF THE INDICES ARE INSIDE THE GRID THE MEASUREMENT IS CUMULATED

        IF((I.LE.NRGG).AND.(I.GE.1) .AND. (J.LE.NGY).AND.(J.GE.1)) THEN

          ICOUNT(KSAT)=ICOUNT(KSAT)+1
          IF (LAVERAGE) THEN
!           the superobservation is assigned to the closest grid point.
            IF (SWH.GE.0.) THEN
              WHME(I,J,KSAT) = WHME(I,J,KSAT) + SWH
              NUMBWH(I,J,KSAT) = NUMBWH(I,J,KSAT) + 1
            ENDIF
          ELSE
            IF (I.NE.IOLD(KSAT) .OR. J.NE.JOLD(KSAT)) THEN

!               NEW BOX

              CALL DIFDATE(CDATEO(KSAT),CDATE,ISHIFT)
              IF (ISHIFT.LT.3.AND.ICOUNT(KSAT).GT.1) THEN
                IF (SWH-SWHO(KSAT).GT.2.0) THEN
                  NR(I,J,KSAT)=NR(I,J,KSAT)+1
                  GO TO 3000
                ENDIF
              ENDIF
              IOLD(KSAT) = I
              JOLD(KSAT) = J
              SWHO(KSAT) = SWH
            ELSE

!             SUDDEN POSITVE STEPS ARE REFUSED.

              IF (SWH-SWHO(KSAT).GT.1.0) THEN
                NR(I,J,KSAT) = NR(I,J,KSAT)+1
                GO TO 3000
              ENDIF
            ENDIF

            IF (SWH.GE.0.) THEN
              WHME(I,J,KSAT) = WHME(I,J,KSAT) + SWH
              WHSE (I,J,KSAT) = WHSE(I,J,KSAT) + SWH**2
              NUMBWH(I,J,KSAT) = NUMBWH(I,J,KSAT) + 1
            ENDIF

!           A PREVIOUS SUDDEN NEGATIVE STEP IS COMPENSATED.

            IF (SWHO(KSAT)-SWH.GT.2.0) THEN
              NR(I,J,KSAT) = NR(I,J,KSAT)+1
              WHME(I,J,KSAT) = WHME(I,J,KSAT) - SWHO(KSAT)
              WHSE (I,J,KSAT) = WHSE(I,J,KSAT) - SWHO(KSAT)**2
              NUMBWH(I,J,KSAT) = NUMBWH(I,J,KSAT) - 1
            ENDIF

            SWHO(KSAT) = SWH
            CDATEO(KSAT) = CDATE
          ENDIF
        ENDIF
!*

        GOTO 3000
      ENDIF


! ----------------------------------------------------------------------

!*    4. DETERMINE AVERAGE VALUE AND ERROR OF THE MEASUREMENTS IN THE
!*       BOX AROUND GRIDPOINT I,J.
!        ------------------------------------------------------------

 4000 CONTINUE

      NWH = 0

      DO ISAT=1,NUMALT
        IF(IBUFRSAT(ISAT).NE.0) THEN
          NUWH=0
          DO J2=1,NGY
            J2SN = NGY-J2+1     
            DO J1=1,NLONRGG(J2SN)
              NUWH=NUWH+NUMBWH(J1,J2,ISAT)
            ENDDO
          ENDDO

          WRITE (IU06,'(" GRDATA: INPUT STATISTICS:")')
          WRITE (IU06,'("         THE SATELLITE IS",
     &                  " ......:", I8)') IBUFRSAT(ISAT)
          WRITE (IU06,'("         NUMBER OF RECORDS READ FROM DATA ",
     &                  "FILE IS ......:", I8)') ICOUNTD(ISAT)
          WRITE (IU06,'("         THEREOF NUMBER OF BLACKLISTED ",
     &                  "RECORDS IS ......:", I8)') ICOUNTBL(ISAT)
          WRITE (IU06,'("         NUMBER OF DATA IN SELECTED PERIOD ",
     &                  "IS ..........:", I8)') ICOUNTT(ISAT)
          WRITE (IU06,'("         NUMBER OF DATA BEFORE SELECTED ",
     &                  "PERIOD IS .....:", I8)') ICOUNTB(ISAT)
          WRITE (IU06,'("         NUMBER OF DATA IN BOXES AROUND ",
     &                  "SEA POINTS IS ..:", I8)') ICOUNT(ISAT)
          WRITE (IU06,'("         NUMBER OF GOOD DATA ABOVE LAND IS ",
     &                  ".............:", I8)') ICOUNTQCF(ISAT)
          WRITE (iu06,'("         NUMBER OF DATA ACCEPTED BY ",
     &                  "''GRDATA'' IS FOR HS :", I8)') NUWH

          DO J=1,NGY
            JSN = NGY-J+1
            DO I=1,NLONRGG(JSN)
              NN = NUMBWH(I,J,ISAT)
              IF (NN.GE.NNMIN) THEN
!               MEAN VALUES
                WHME(I,J,ISAT) = WHME(I,J,ISAT)/REAL(NN)
!               COMPUTE THE STANDARD DEVIATION.
                IF (.NOT. LAVERAGE) THEN
                  WHSE2 =  (WHSE(I,J,ISAT)-WHME(I,J,ISAT)**2*REAL(NN))
     &                     / (REAL(NN)-1.)
                  WHSE(I,J,ISAT) = SQRT( MAX(WHSE2,0.))
                ENDIF
              ENDIF
            ENDDO
          ENDDO

! ----------------------------------------------------------------------

!*    5. FLAG GRID POINTS, IF THERE ARE FEW OR STRONGLY SCATTERED DATA.
!        --------------------------------------------------------------

!         ONLY KEEP THE RELEVANT INFORMATION
!         ALL DATA THAT ARE NOT OVER SEA ARE DISCARDED HERE !

          IF (LAVERAGE) THEN
            DO IJ = NSTART(1), NEND(NPROC)
              I = IXLG(IJ,1)
              J = NGY- KXLT(IJ,1) +1
              NN = NUMBWH(I,J,ISAT)
              IF (NN.EQ.0) THEN
!               IF THERE ARE NO MEASUREMENTS
                WHME(I,J,ISAT)=-3.
              ELSEIF (WHME(I,J,ISAT).LT.0.50) THEN
                WHME(I,J,ISAT) = -WHME(I,J,ISAT)
              ELSE
                NWH = NWH+1
              ENDIF
            ENDDO
          ELSE

           DO IJ = NSTART(1), NEND(NPROC)
            I = IXLG(IJ,1)
            J = NGY- KXLT(IJ,1) +1
            NN = NUMBWH(I,J,ISAT)
            IF (NN.EQ.0) THEN
!             IF THERE ARE NO MEASUREMENTS
              WHME(I,J,ISAT)=-3.
            ELSEIF (NN.GT.0. AND. NN.LT.NNMIN) THEN
!             IF THERE ARE FEW MEASUREMENTS.
              WHME(I,J,ISAT)=-2.
            ELSEIF (NN.GE.NNMIN) THEN
!             IF THE VARIANCE IS TOO LARGE.
!             OR SWH IS TOO SMALL  OR THERE ARE TOO MANY SPIKES
              XX = REAL(NR(I,J,ISAT))/REAL(NN)
              RFWM = MAX(0.5,0.25*WHME(I,J,ISAT))
              IF (WHSE(I,J,ISAT).GT.RFWM .OR.
     &            WHME(I,J,ISAT).LT.0.50 .OR.
     &            XX.GT.0.1) THEN
                WHME(I,J,ISAT) = -WHME(I,J,ISAT)
              ELSE
                NWH = NWH+1
              ENDIF
            ENDIF
           ENDDO

          ENDIF
        ENDIF
      ENDDO  ! ON ISAT

!*
      DEALLOCATE(WHSE)
      DEALLOCATE(NUMBWH)
      DEALLOCATE(NR)

      IF(NWH.GT.0) THEN
        ALLOCATE(IJALT(NWH,NIJALT))
        ALLOCATE(ALTDATA(NWH,NALTDT))

        IF(ITEST.GE.2)
     &     WRITE(IU06,'(A)') 'S.ID. ORIG. CORR. INCR.'
        IOBS=0
        DO ISAT=1,NUMALT
          LLDOBC=.FALSE.
          IF(IBUFRSAT(ISAT).NE.0) THEN
            IF(LLBCTEXIST) THEN
              DO IBCS=1, NSATBC
                IF (IBUFRSAT(ISAT).EQ.ISATBC(IBCS)) THEN
                  DO IBCT=1,NBCT(IBCS)
                    IF (CDTPRO.GE.CBCDATESTART(IBCS,IBCT) .AND. 
     &                  CDTPRO.LE.CBCDATEEND  (IBCS,IBCT)      ) THEN
                      LLDOBC=.TRUE.
                      EXIT
                    ENDIF
                  ENDDO
                  EXIT
                ENDIF
              ENDDO
            ENDIF

            DO IJ = NSTART(1), NEND(NPROC)
              IX = IXLG(IJ,1)
              IY = NGY- KXLT(IJ,1) +1
              IF(WHME(IX,IY,ISAT) .GT. 0. ) THEN
                IOBS=IOBS+1
                ALTDATA(IOBS,3) = WHME(IX,IY,ISAT)  !  1 REPLACED BY 3
                IJALT(IOBS,1) = IJ
                IJALT(IOBS,2) = IBUFRSAT(ISAT)

                ! BIAS CORRECTION IS DONE IF NEEDED
                IF(LLDOBC) THEN
                  ZI=MIN(FLOAT(MAXNBCINC), WHME(IX,IY,ISAT)/
     &                                             BCINC(IBCS,IBCT))
                  I=INT(ZI)+1
                  IP1=MIN(MAXNBCINC, I+1)
                  BCIV=(I+1-ZI)*BCTABLE(IBCS,IBCT,I)+
     &                 (ZI-I)  *BCTABLE(IBCS,IBCT,IP1) 
                  ALTDATA(IOBS,1)=WHME(IX,IY,ISAT)-BCIV
                  IF(ITEST.GE.2)
     &              WRITE(IU06,'(I5,3(1X,F5.2))') ISAT,
     &                     ALTDATA(IOBS,3), ALTDATA(IOBS,1), BCIV
                ELSE
                  ALTDATA(IOBS,1)=WHME(IX,IY,ISAT)
                ENDIF

              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF

      DEALLOCATE(WHME)

      IF(NBL.GT.0) THEN
        DEALLOCATE(STATID)
        DEALLOCATE(SENSOR)
        DEALLOCATE(CDATESTART)
        DEALLOCATE(CDATEEND)
        DEALLOCATE(LATMI)
        DEALLOCATE(LONGMI)
        DEALLOCATE(LATMA)
        DEALLOCATE(LONGMA)
        DEALLOCATE(DLON)
        DEALLOCATE(PARAM)
        DEALLOCATE(VALUE_MIN)
        DEALLOCATE(VALUE_MAX)
      ENDIF

      WRITE (IU06,*) ' GRDATA: GRIDDED MEASUREMENTS ARE SELECTED:'
      WRITE (IU06,*) '         WAVE HEIGHT AT NWH = ', NWH,' GRIDPOINTS'

      CLOSE (iume,IOSTAT=ioszx)
      IF(ITEST.GT.0)
     & WRITE(iu06,*) ' UNIT ',iume,' CLOSED WITH IOSTAT= ',ioszx

      RETURN
      END SUBROUTINE GRDATA
