      SUBROUTINE GRDATA (WHME, USME, NSTART, NEND)

!--------------------------------------------------------------------

!**** *GRDATA* - TRANSFERS MEASUREMENTS TO MODEL GLOBAL GRID,
!****            DISCARDING UNRELIABLE DATA.

!     P.LIONELLO     ECMWF       APRIL 1990

!     PURPOSE.
!     --------

!       PREPARE THE DATA FOR OPTIMAL INTERPOLATION.

!**   INTERFACE.
!     ----------

!       *CALL* *GRDATA (WHME, USME)*

!        *WHME*   REAL    OUTPUT FIELD OF MEASURED WAVE HEIGHTS.
!        *USME*   REAL    OUTPUT FIELD OF MEASURED WIND SPEEDS.
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

!-----------------------------------------------------------------------

      USE YOWGRID  , ONLY : NLONRGG  ,IGL      ,IJS      ,IJL
      USE YOWICE   , ONLY : ICEMASK
      USE YOWINTP  , ONLY : CDG
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,IPER     ,AMOWEP   ,
     &            AMOSOP   ,XDELLA   ,ZDELLO
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : NPROC
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWSTAT  , ONLY : CDTPRO   ,IDELALT
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------

      LOGICAL EOFD
      CHARACTER*12 CBEGINDT, CENDDT
      CHARACTER*14 CDATE,CDATEO
      DIMENSION WHME(NGX,NGY), USME(NGX,NGY),
     &          WHSE(NGX,NGY), USSE(NGX,NGY),
     &          NUMBWH(NGX,NGY), NUMBUS(NGX,NGY),NR(NGX,NGY)

!*     VARIABLE     TYPE     PURPOSE.
!      --------     ----     --------

!     *CBEGINDT*    CHAR*12  BEGINNING DATE OF THE TIME WINDOW
!                            DURING WHICH DATA ARE USED FOR ANALYSIS.
!     *CENDDT*      CHAR*12  END DATE OF THE TIME WINDOW
!                            DURING WHICH DATA ARE SAVE FOR ANALYSIS.
!     *WHME*        REAL     GRID WAVE HEIGHT.
!     *USME*        REAL     GRID WIND SPEED.
!     *WHSE*        REAL     WAVE HEIGHT ROOT MEAN SQUARE ERROR.
!     *USSE*        REAL     WIND SPEED ROOT MEAN SQUARE ERROR.
!     *NUMBWH*      INTEGER  NUMBER OF HS MEASUREMENTS INSIDE A BOX.
!     *NUMBUS*      INTEGER  NUMBER OF US MEASUREMENTS INSIDE A BOX.

! ----------------------------------------------------------------------

      INTEGER,DIMENSION(NPROC) :: NSTART, NEND
! ----------------------------------------------------------------------
      INTEGER,ALLOCATABLE :: ICEMASKG(:,:)

      INTEGER, PARAMETER :: jpice=1000

      REAL,ALLOCATABLE :: RLAICE(:),RLOICE(:),SWHICE(:),WSICE(:)
      CHARACTER*12,ALLOCATABLE :: CDAICE(:)
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

      ALLOCATE(RLAICE(JPICE),RLOICE(JPICE),SWHICE(JPICE),WSICE(JPICE))
      ALLOCATE(CDAICE(JPICE))
      ICOUNT=0
      ICOUNTT=0
      ICOUNTD=0
      icountf=0
      icountb=0
      icountice=0
      icountqcf=0
      DO J=1,NGY
        DO I=1,NGX
          WHME(I,J) = 0.
          USME(I,J) = 0.
          WHSE(I,J) = 0.
          USSE(I,J) = 0.
          NUMBWH(I,J) = 0
          NUMBUS(I,J) = 0
          NR(I,J)=0
        ENDDO
      ENDDO

!     INITIALISE THE ICE MASK ARRAY

      ALLOCATE(ICEMASKG(NGX,NGY))

      ICEMASKG=1

      IF(LMESSPASS) THEN
        DO IJ = NSTART(1), NEND(NPROC)
          IF (ICEMASK(IJ,1).EQ.0) THEN
            IX = IXLG(IJ,1)
            IY = NGY- KXLT(IJ,1) +1
            ICEMASKG(IX,IY)=0
          ENDIF
        ENDDO
      ELSE
        DO IG = 1, IGL 
          DO IJ = IJS(IG),IJL(IG) 
            IF (ICEMASK(IJ,IG).EQ.0) THEN
              IX = IXLG(IJ,IG)
              IY = NGY- KXLT(IJ,IG) +1
              ICEMASKG(IX,IY)=0
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      IOLD = 0
      JOLD = 0
      CDATEO = CBEGINDT(1:12)//'00'

      DO iume=33,99,1
        INQUIRE ( UNIT=iume, OPENED=opened)
        IF ( .NOT. opened) THEN
          WRITE(iu06,*)
     &     ' GRDATA : TAKES UNIT ',iume,' AS ALTIMETER DATA INPUT UNIT'
          GOTO 3000
        ENDIF
      ENDDO

! ----------------------------------------------------------------------

!*    3. READ THE MEASUREMENT TAPE AND CUMULATE THE DATA ON GRID POINTS.
!        ---------------------------------------------------------------

 3000 CONTINUE
      flagged = .TRUE.
      CALL READSAT (IU06, iume, CDTPRO,
     &              CDATE, RLAT, RLON, SWH, WS, flagged, EOFD)
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

!     DATE IS INSIDE THE TIME WINDOW and record flagged as reliable.
!     THE INDECES ARE COMPUTED

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

!     IF THE GRID POINT IS SET TO BE ICE TAKE NEXT MEASUREMENT.

          IF (ICEMASKG(I,J).EQ.0) THEN
            icountice = icountice + 1
            if (icountice .le. jpice) then
              cdaice(icountice)=cdate(1:12)
              rlaice(icountice)=rlat
              rloice(icountice)=rlon
              swhice(icountice)=swh
              wsice(icountice)=ws
            endif
            GOTO 3000
          ENDIF

!     SOME MEASURENTS ARE VERY CLOSE TO LAND, PASSED THE QC AS
!     BEING ON SEA BUT REGARDED AS BEING ON LAND BY GRDATA.
!     THESE ARE ELIMINATED QUIETLY HERE.
!     NOTE: CDG CAN BE = 0.0 FOR SEA POINTS WHEN GRDATA IS USED 
!     BY *RFL4WAM*

          IF (CDG(I,J).LT.0.000) THEN
            icountQCF = icountQCF + 1
            GOTO 3000
          ENDIF


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
            WSO  = WS
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
            IF (WS.GE.0.) THEN
              USME(I,J) = USME(I,J) + WS
              USSE (I,J) = USSE(I,J) + WS**2
              NUMBUS(I,J) = NUMBUS(I,J) + 1
            ENDIF
          ENDIF

!     A PREVIOUS SUDDEN NEGATIVE STEP IS COMPENSATED.

          IF (SWHO-SWH.GT.2.0) then
            NR(I,J) = NR(I,J)+1
            WHME(I,J) = WHME(I,J) - SWHO
            WHSE (I,J) = WHSE(I,J) - SWHO**2
            NUMBWH(I,J) = NUMBWH(I,J) - 1
            IF (WSO.GE.0.) THEN
              USME(I,J) = USME(I,J) - WSO
              USSE (I,J) = USSE(I,J) - WSO**2
              NUMBUS(I,J) = NUMBUS(I,J) - 1
            ENDIF
          ENDIF

          SWHO = SWH
          WSO = WS
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
      nuus=0
      DO j2=1,ngy
        j2sn = ngy-j2+1     
        DO j1=1,nlonrgg(j2sn)
          nuwh=nuwh+numbwh(j1,j2)
          nuus=nuus+numbus(j1,j2)
        ENDDO
      ENDDO

      WRITE (IU06,'(" GRDATA: INPUT STATISTICS:")')
      WRITE (IU06,'("         NUMBER OF RECORDS READ FROM DATA FILE IS",
     &              "......:", I8)') ICOUNTD
      WRITE (IU06,'("         THEREOF NUMBER OF RECORDS WITH FLAGS IS ",
     &              "......:", I8)') ICOUNTF
      WRITE (IU06,'("         NUMBER OF DATA IN SELECTED PERIOD IS ...",
     &              "......:", I8)') ICOUNTT
      WRITE (IU06,'("         NUMBER OF DATA BEFORE SELECTED PERIOD IS",
     &              "......:", I8)') ICOUNTB
      WRITE (IU06,'("         NUMBER OF DATA IN BOXES AROUND SEA POINT",
     &              "S IS..:", I8)') ICOUNT
      WRITE (IU06,'("         NUMBER OF GOOD DATA ABOVE LAND IS.......",
     &              "......:", I8)') ICOUNTQCF
      WRITE (IU06,'("         NUMBER OF GOOD DATA ABOVE ICE IS........",
     &              "......:", I8)') ICOUNTICE
      WRITE (iu06,'("         NUMBER OF DATA ACCEPTED BY ''GRDATA'' IS",
     &              " FOR wh :", I8, /, 
     &"                                            AND FOR us :", I8)')
     &    nuwh, nuus

      DO iubl=30,99,1
        INQUIRE ( UNIT=iubl, OPENED=opened)
        IF ( .NOT. opened.and.ITEST.GT.0.and.ITEST.GT.0) THEN
          WRITE(IU06,*)
     &     ' GRDATA: TAKES UNIT ',iubl,' AS BLACKLIST OUTPUT UNIT'
          CALL FLUSH(IU06)
          GOTO 4160
        ENDIF
      ENDDO
 4160 CONTINUE
      OPEN (IUBL,FILE='BLL'//CDTPRO)
      WRITE(IUBL,'(A12)') CDTPRO
      DO j=1,MIN(icountice,jpice)
        WRITE(IUBL,'(a12,4f8.2)') cdaice(j),rlaice(j),rloice(j),
     &                            swhice(j),wsice(j)
      ENDDO

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

      DO J=1,NGY
        JSN = NGY-J+1
        DO I=1,NLONRGG(JSN)
          NN = NUMBUS(I,J)
          IF (NN.GE.NNMIN) THEN

!           MEAN VALUES

            USME(I,J) = USME(I,J)/REAL(NN)

!           COMPUTE THE STANDARD DEVIATION.

            USSE2 =  (USSE(I,J)-USME(I,J)**2*REAL(NN))
     &       / (REAL(NN)-1.)
            USSE(I,J) = SQRT( MAX(USSE2,0.))
          ENDIF
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    5. FLAG GRID POINTS, IF THERE ARE FEW OR STRONGLY SCATTERED DATA.
!        --------------------------------------------------------------

      NWH = 0
      DO J=1,NGY
        JSN = NGY-J+1
        DO I=1,NLONRGG(JSN)
          NN = NUMBWH(I,J)
          IF (NN.EQ.0) THEN

!           IF THERE ARE NO MEASUREMENTS

            WHME(I,J)=-3.
          ELSEIF (NN.GT.0. AND. NN.LT.NNMIN) THEN

!           IF THERE ARE FEW MEASUREMENTS.

            WHME(I,J)=-2.
          ELSEIF (NN.GE.NNMIN) THEN

!           IF THE VARIANCE IS TOO LARGE.
!           OR SWH IS TOO SMALL  OR THERE ARE TOO MANY SPIKES

            XX = REAL(NR(I,J))/REAL(NN)
            RFWM = MAX(0.5,0.25*WHME(I,J))
            IF (WHSE(I,J).GT.RFWM .OR. WHME(I,J).LT.0.50 .OR.
     &       XX.GT.0.1) THEN
              WHME(I,J) = -WHME(i,j)
            ELSE
              NWH = NWH+1
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      NUS = 0
      DO J=1,NGY
        JSN = NGY-J+1
        DO I=1,NLONRGG(JSN)
          NN = NUMBUS(I,J)
          IF (NN.EQ.0) THEN

!           IF THERE ARE NO MEASUREMENTS

            USME(I,J)=-3.
          ELSEIF (NN.GT.0. AND. NN.LT.NNMIN) THEN

!           IF THERE ARE FEW MEASUREMENTS.

            USME(I,J)=-2.
          ELSEIF (NN.GE.NNMIN)THEN

!           IF THE VARIANCE IS TOO LARGE.
!           OR THERE IS NO CORRESPONDING SWH

            RFU = 0.0
            IF (USME(I,J).GT.0.0001) RFU = USSE(I,J)/USME(I,J)
            IF (RFU.GE.0.5.OR.WHME(I,J).LE.0.) THEN
              USME(I,J) = -USME(I,J)
            ELSE
              NUS = NUS + 1
            ENDIF
          ENDIF
        ENDDO
      ENDDO

!*

      WRITE (IU06,*) ' GRDATA: GRIDDED MEASUREMENTS ARE SELECTED:'
      WRITE (IU06,*) '         WAVE HEIGHT AT NWH = ', NWH,' GRIDPOINTS'
      WRITE (IU06,*) '         WIND SPEED  AT NUS = ', NUS,' GRIDPOINTS'

      CLOSE (iume,IOSTAT=ioszx)
      IF(ITEST.GT.0)
     & WRITE(iu06,*) ' UNIT ',iume,' CLOSED WITH IOSTAT= ',ioszx
      CLOSE (iubl,IOSTAT=ioszx)
      IF(ITEST.GT.0)
     & WRITE(iu06,*) ' UNIT ',iubl,' CLOSED WITH IOSTAT= ',ioszx

      IF(ALLOCATED(ICEMASKG)) DEALLOCATE(ICEMASKG)
      DEALLOCATE(RLAICE,RLOICE,SWHICE,WSICE,CDAICE)

      RETURN
      END SUBROUTINE GRDATA
