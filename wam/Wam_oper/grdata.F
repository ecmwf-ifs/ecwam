      SUBROUTINE GRDATA(NWH, NSTART, NEND)

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

      USE YOWALTAS , ONLY : NIJALT   ,NALTDT   ,IJALT    ,ALTDATA 
      USE YOWGRID  , ONLY : NLONRGG  ,IGL      ,IJS      ,IJL
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,IPER     ,AMOWEP   ,
     &            AMOSOP   ,XDELLA   ,ZDELLO
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : NPROC
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWSTAT  , ONLY : CDTPRO   ,IDELALT
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------

      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NUMBWH, NR
      REAL, ALLOCATABLE, DIMENSION(:,:) ::  WHME, WHSE
      CHARACTER*12 CBEGINDT, CENDDT
      CHARACTER*14 CDATE,CDATEO
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
      LOGICAL flagged, opened

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

      IF(XDELLA.GE.0.5) THEN
        NNMIN=4
      ELSE
        NNMIN=2
      ENDIF
      WRITE(IU06,*) '   GRIDDED DATA WILL BE KEPT WHEN THE NUMBER OF'
      WRITE(IU06,*) '   OBSERVATIONS PER GRID BOX IS LARGER THAN ',
     &  NNMIN-1

! ----------------------------------------------------------------------

!*    2. INITIALIZE ARRAYS, counters and units.
!        --------------------------------------

      ICOUNT=0
      ICOUNTT=0
      ICOUNTD=0
      icountf=0
      icountb=0
      icountqcf=0
      icountbl=0
      ALLOCATE(WHME(NGX,NGY))
      ALLOCATE(WHSE(NGX,NGY))
      ALLOCATE(NUMBWH(NGX,NGY))
      ALLOCATE(NR(NGX,NGY))
      DO J=1,NGY
        DO I=1,NGX
          WHME(I,J) = 0.
          WHSE(I,J) = 0.
          NUMBWH(I,J) = 0
          NR(I,J)=0
        ENDDO
      ENDDO

      IOLD = 0
      JOLD = 0
      CDATEO = CBEGINDT(1:12)//'00'

      DO iume=33,99,1
        INQUIRE ( UNIT=iume, OPENED=opened)
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
          READ(IOBL,'(A100)') CDUM
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
          READ(IOBL,'(A100)') CDUM
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

!*    3. READ THE MEASUREMENT TAPE AND CUMULATE THE DATA ON GRID POINTS.
!        ---------------------------------------------------------------

 3000 CONTINUE
      flagged = .TRUE.
      CALL READSAT (IU06, iume, CDTPRO, IDENTI, ISENSOR,
     &              CDATE, RLAT, RLON, IDES_HS, SWH,
     &              IDES_WS, WS, flagged, EOFD)
!!! note this will need to be adapted when there is more than one 
!!! satellite
      IF (EOFD) GOTO 4000
      ICOUNTD = ICOUNTD+1

!     CHECK DATE.

      IF (CDATE(1:12).LT.CBEGINDT)  then
        icountb = icountb + 1
        GOTO 3000
      endif

      IF (CDATE(1:12).LE.CENDDT) THEN

!     CHECK QUALITY.

        IF (flagged) THEN
          icountf=icountf+1
          GOTO 3000
        ENDIF

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
                icountbl=icountbl+1
                GOTO 3000
              ENDIF
            ENDIF
          ENDDO
        ENDIF

!       DATE IS INSIDE THE TIME WINDOW AND RECORD FLAGGED AS RELIABLE.
!       THE INDICES ARE COMPUTED

        ICOUNTT = ICOUNTT+1
        J       = NGY-NINT((RLAT-AMOSOP)/XDELLA)
        XI      = MOD(RLON-AMOWEP+720.,360.)

        JSN     = NGY-J+1
        IF((J.LE.0).OR.(J.GT.NGY)) JSN=1
        NRGG    = NLONRGG(JSN)
        I       = NINT(XI/ZDELLO(JSN)) + 1
        IF (IPER.EQ.1.AND.I.EQ.NRGG+1) I = 1

!     IF THE INDICES ARE INSIDE THE GRID THE MEASUREMENT IS CUMULATED

        IF((I.LE.NRGG).AND.(I.GE.1) .AND. (J.LE.NGY).AND.(J.GE.1)) THEN

          ICOUNT=ICOUNT+1
          IF (I.NE.IOLD .OR. J.NE.JOLD) THEN

!               NEW BOX

            CALL DIFDATE(CDATEO,CDATE,ISHIFT)
            IF (ISHIFT.LT.3.AND.ICOUNT.GT.1) THEN
              IF (SWH-SWHO.gt.2.0) THEN
                NR(I,J)=NR(I,J)+1
                GO TO 3000
              ENDIF
            ENDIF
            IOLD = I
            JOLD = J
            SWHO = SWH
          ELSE

!     SUDDEN POSITVE STEPS ARE REFUSED.

            IF (SWH-SWHO.GT.1.0) THEN
              NR(I,J) = NR(I,J)+1
              GO TO 3000
            ENDIF
          ENDIF

          IF (SWH.GE.0.) THEN
            WHME(I,J) = WHME(I,J) + SWH
            WHSE (I,J) = WHSE(I,J) + SWH**2
            NUMBWH(I,J) = NUMBWH(I,J) + 1
          ENDIF

!     A PREVIOUS SUDDEN NEGATIVE STEP IS COMPENSATED.

          IF (SWHO-SWH.GT.2.0) then
            NR(I,J) = NR(I,J)+1
            WHME(I,J) = WHME(I,J) - SWHO
            WHSE (I,J) = WHSE(I,J) - SWHO**2
            NUMBWH(I,J) = NUMBWH(I,J) - 1
          ENDIF

          SWHO = SWH
          CDATEO = CDATE
        ENDIF
!*

        GOTO 3000
      ENDIF


! ----------------------------------------------------------------------

!*    4. DETERMINE AVERAGE VALUE AND ERROR OF THE MEASUREMENTS IN THE
!*       BOX AROUND GRIDPOINT I,J.
!        ------------------------------------------------------------

 4000 CONTINUE
      nuwh=0
      DO j2=1,ngy
        j2sn = ngy-j2+1     
        DO j1=1,nlonrgg(j2sn)
          nuwh=nuwh+numbwh(j1,j2)
        ENDDO
      ENDDO

      WRITE (IU06,'(" GRDATA: INPUT STATISTICS:")')
      WRITE (IU06,'("         NUMBER OF RECORDS READ FROM DATA FILE IS",
     &              "......:", I8)') ICOUNTD
      WRITE (IU06,'("         THEREOF NUMBER OF RECORDS WITH FLAGS IS ",
     &              "......:", I8)') ICOUNTF
      WRITE (IU06,'("         THEREOF NUMBER OF BLACKLISTED RECORDS IS",
     &              "......:", I8)') ICOUNTBL
      WRITE (IU06,'("         NUMBER OF DATA IN SELECTED PERIOD IS ...",
     &              "......:", I8)') ICOUNTT
      WRITE (IU06,'("         NUMBER OF DATA BEFORE SELECTED PERIOD IS",
     &              "......:", I8)') ICOUNTB
      WRITE (IU06,'("         NUMBER OF DATA IN BOXES AROUND SEA POINT",
     &              "S IS..:", I8)') ICOUNT
      WRITE (IU06,'("         NUMBER OF GOOD DATA ABOVE LAND IS.......",
     &              "......:", I8)') ICOUNTQCF
      WRITE (iu06,'("         NUMBER OF DATA ACCEPTED BY ''GRDATA'' IS",
     &              " FOR wh :", I8)') nuwh

      DO J=1,NGY
        JSN = NGY-J+1
        DO I=1,NLONRGG(JSN)
          NN = NUMBWH(I,J)
          IF (NN.GE.NNMIN) THEN

!           MEAN VALUES

            WHME(I,J) = WHME(I,J)/REAL(NN)

!           COMPUTE THE STANDARD DEVIATION.

            WHSE2 =  (WHSE(I,J)-WHME(I,J)**2*REAL(NN))
     &       / (REAL(NN)-1.)
            WHSE(I,J) = SQRT( MAX(WHSE2,0.))
          ENDIF
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    5. FLAG GRID POINTS, IF THERE ARE FEW OR STRONGLY SCATTERED DATA.
!        --------------------------------------------------------------

!     ONLY KEEP THE RELEVANT INFORMATION
!     ALL DATA THAT ARE NOT OVER SEA ARE DISCARDED HERE !

      NWH = 0
      DO IJ = NSTART(1), NEND(NPROC)
        I = IXLG(IJ,1)
        J = NGY- KXLT(IJ,1) +1
        NN = NUMBWH(I,J)
        IF (NN.EQ.0) THEN

!         IF THERE ARE NO MEASUREMENTS

          WHME(I,J)=-3.
        ELSEIF (NN.GT.0. AND. NN.LT.NNMIN) THEN

!         IF THERE ARE FEW MEASUREMENTS.

          WHME(I,J)=-2.
        ELSEIF (NN.GE.NNMIN) THEN

!         IF THE VARIANCE IS TOO LARGE.
!         OR SWH IS TOO SMALL  OR THERE ARE TOO MANY SPIKES

          XX = REAL(NR(I,J))/REAL(NN)
          RFWM = MAX(0.5,0.25*WHME(I,J))
          IF (WHSE(I,J).GT.RFWM .OR. WHME(I,J).LT.0.50 .OR.
     &     XX.GT.0.1) THEN
            WHME(I,J) = -WHME(i,j)
          ELSE
            NWH = NWH+1
          ENDIF
        ENDIF
      ENDDO

!*
      DEALLOCATE(WHSE)
      DEALLOCATE(NUMBWH)
      DEALLOCATE(NR)

      IF(NWH.GT.0) THEN
        ALLOCATE(IJALT(NWH,NIJALT))
        ALLOCATE(ALTDATA(NWH,NALTDT))

        IOBS=0
        DO IJ = NSTART(1), NEND(NPROC)
          IX = IXLG(IJ,1)
          IY = NGY- KXLT(IJ,1) +1
          IF(WHME(IX,IY) .GT. 0. ) THEN
            IOBS=IOBS+1
            ALTDATA(IOBS,1) = WHME(IX,IY)
            IJALT(IOBS,1) = IJ
            IJALT(IOBS,2) = IDENTI 
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
