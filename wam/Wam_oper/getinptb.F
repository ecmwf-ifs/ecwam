      PROGRAM GETINPTB

! --------------------------------------------------------------------

!     PETER JANSSEN     ECMWF        FEBRUARI 1987
!     JEAN BIDLOT       ECMWF        MAY 1995  MIGRATION TO FUJITSU

!     PURPOSE
!     -------

!          THIS PROGRAM WRITES THE INPUT DATA CARDS FOR MARS.
!          THE PARAMETER SPECIFICATIONS FOR MARS ARE READ FROM FILE
!          WAMINFO.
!          ALSO THE FILE NAME FOR THE RESULTS OF THE MARS RETRIEVAL
!          IS OUTPUT IN WINDFILENAME AND SHOULD BE USED TO SAVE MARS
!          RETRIEVAL ON LINE.

!     INTERFACE.
!     ---------

!     *FILE DEFINITION*

!      IU04  - INFO FROM STATUSFILE WAMINFO
!      IU06  - PRINTER OUTPUT UNIT
!      IU22  - CONTAINS RETRIEVAL COMMANDS FOR MARS

!     PARAMETERS.
!     ----------

!     VARIABLE    TYPE      PURPOSE
!     --------    ----      -------

!     CDAT1       CHAR*12   BEGINDATE OF RUN
!     CDAT2       CHAR*12   ENDDATE OF RUN
!     IDATEFC     INTEGER   START DATE FORECAST RUN
!     IANALPD     INTEGER   PERIOD OF ANALYSIS
!     IFOREPD     INTEGER   FORECAST PERIOD
!     ISHIFT      INTEGER   TIMESTEP OF WINDFIELDS FROM MARS


!     EXTERNALS.
!     ---------

!     *INCDATE*  -  FROM WAMCRLIB INCREASES DATE TIME GROUP
!     *READSTA*  -  READ WAMINFO FILE.

!----------------------------------------------------------------------

      USE YOWCARD  , ONLY : CARD     ,JPCL
      USE YOWMAP   , ONLY : AMOWEP   ,AMOSOP   ,AMOEAP   ,AMONOP   ,
     &            XDELLA   ,XDELLO



      CHARACTER*12 CDAT1, CDAT2, CPBDT, CPEDT , CDATF, CRSDT, CLSDT
      CHARACTER*12 CDUMMY
      CHARACTER CDATS*8,CDATEFC*8,CDE*8,CGRID*3,CHOUS*2,CHOUE*2
      CHARACTER CANPERIOD*3,CFCPERIOD*3,CSHIFT*2
      CHARACTER MARS(21)*80,CFILENAME*44
      INTEGER ISTAT(3)


! ---------------------------------------------------------------------

!*    1. PRINT WAMINFO FILE AND THE GRID INFORMATION
!        -------------------------------------------

      IU04 = 4 
      IU06 = 6
      IU07 = 7
      IU22 = 22
      CARD = ' '
      CALL READSTA (IU04, CDAT1, CDAT2, IANALPD, IFOREPD, ISHIFT,
     &                    CRSDT, CLSDT, CPBDT, CPEDT, IRST, INF, ISTAT)
      WRITE(IU06,'(1H1,'' START OF GETINPT'',//,
     &       '' THE WAMINFO FILE IS AS FOLLOWS:'',/)')
      DO I=1,JPCL
        WRITE(IU06,'(4X,A72)') CARD(I)
      ENDDO
      ISHIFT = ISHIFT / 3600


      CALL READPRE(IU07,IREFRA)
      WRITE(IU06,*) ' '  
      WRITE(IU06,*) '    THE GRID IS AS FOLLOWS : '  
      WRITE(IU06,*) '    NORTH : ',AMONOP 
      WRITE(IU06,*) '    WEST  : ',AMOWEP
      WRITE(IU06,*) '    SOUTH : ',AMOSOP 
      WRITE(IU06,*) '    EAST  : ',AMOEAP
      WRITE(IU06,*) '    X-INCREMENT ; ',XDELLO 
      WRITE(IU06,*) '    Y-INCREMENT ; ',XDELLA 
      WRITE(IU06,*) ' ' 

! ---------------------------------------------------------------------

!*    2. DETERMINE DATES AND HOURS.
!        --------------------------

      CDATS=CDAT1(1:8)
      READ(CDATS,'(I8)') IDATS
      CHOUS=CDAT1(9:10)
      READ(CHOUS,'(I2)') IHOUS
      IF      (IHOUS.GE.0  .AND. IHOUS.LT. 6) THEN
        IHOUS = 0
      ELSEIF (IHOUS.GE.6  .AND. IHOUS.LT.12) THEN
        IHOUS = 6
      ELSEIF (IHOUS.GE.12 .AND. IHOUS.LT.18) THEN
        IHOUS = 12
      ELSEIF (IHOUS.GE.18) THEN
        IHOUS = 18
      ENDIF

      CDATF=CDAT1
      CALL INCDATE(CDATF,IANALPD)
      CDATEFC=CDATF(1:8)
      READ(CDATEFC,'(I8)') IDATEFC
      CHOUE=CDATF(9:10)
      READ(CHOUE,'(I2)') IHOUE
      IF      (IHOUE.GT.0  .AND. IHOUE.LE. 6) THEN
        IHOUE = 6
      ELSEIF (IHOUE.GT.6  .AND. IHOUE.LE.12) THEN
        IHOUE = 12
      ELSEIF (IHOUE.GT.12 .AND. IHOUE.LE.18) THEN
        IHOUE = 18
      ELSEIF (IHOUE.GE.18) THEN
        IHOUE = 0
        CDUMMY=CDATEFC//'0000'
        CALL INCDATE (CDUMMY, 24*3600)
        CDATEFC = CDUMMY(1:8)
        READ(CDATEFC,'(I8)') IDATEFC
      ENDIF

! ---------------------------------------------------------------------

!*    3. INITIALISE UNIT TO RETRIEVE MARS DATA AND LINE COUNTER.
!        -------------------------------------------------------

      IUNIT = 80
      ILINE = 0

! ---------------------------------------------------------------------

!*    3. WRITE MARS COMMANDS.
!        --------------------

      MARS(2) = CARD(2)
      MARS(3) = CARD(3)
      WRITE(MARS(4),'("TYPE=AN,LEVTYPE=SFC,LEVELIST=OFF,")')

      WRITE(MARS(5),'("GRID=",F5.2,"/",F5.2,",FORMAT=PACKED,")')
     &      XDELLO, XDELLA
      WRITE(MARS(6),'("AREA=",f7.2,"/",f7.2,"/",f7.2,"/",f7.2,",")')
     &      AMONOP, AMOWEP, AMOSOP, AMOEAP

      IGRID=INT(100*XDELLA) 
      WRITE(CGRID,'(I3.3)') IGRID

      IF (IANALPD.GE.0) THEN

!*    3.1 RETRIEVAL OF ANALYSIS.
!         ----------------------


!    ALL FIELDS ARE ON ONE DAY.

        IF (IDATS.EQ.IDATEFC) THEN
          ILINE = ILINE + 1
          IF (ILINE.EQ.1) THEN
            IF (IHOUS.EQ.IHOUE) THEN
              WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     &         ",TIME=",I2.2,",")')
     &         IDATS, IHOUS
            ELSE
              WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     &         ",TIME=",I2.2,"/TO/",I2.2,",")')
     &         IDATS, IHOUS, IHOUE
            ENDIF
            WRITE(MARS(7), '("TARGET=""fort.",I2,"""")') IUNIT
            ILINE = 7
          ELSE
            IF (IHOUS.EQ.IHOUE) THEN
              WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     &         ",TIME=",I2.2)')
     &         IDATS, IHOUS
            ELSE
              WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     &         ",TIME=",I2.2,"/TO/",I2.2)')
     &         IDATS, IHOUS, IHOUE
            ENDIF
          ENDIF
        ELSE

!    MORE THAN ONE DAY

          IF (IHOUS.NE.0) THEN
            ILINE = ILINE + 1
            IF (ILINE.EQ.1) THEN
              IF (IHOUS.EQ.24-ISHIFT) THEN
                WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     &           ",TIME=",I2.2,",")')
     &           IDATS, IHOUS
              ELSE
                WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     &           ",TIME=",I2.2,"/TO/",I2.2,",")')
     &           IDATS, IHOUS, 24-ISHIFT
              ENDIF
              ILINE = 7
              WRITE(MARS(ILINE), '("TARGET=""fort.",I2,"""")')
     &         IUNIT
            ELSE
              IF (IHOUS.EQ.24-ISHIFT) THEN
                WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     &           ",TIME=",I2.2)')
     &           IDATS, IHOUS
              ELSE
                WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     &           ",TIME=",I2.2,"/TO/",I2.2)')
     &           IDATS, IHOUS, 24-ISHIFT
              ENDIF
            ENDIF
            CDUMMY=CDATS//'0000'
            CALL INCDATE (CDUMMY, 24*3600)
            CDATS = CDUMMY(1:8)
            READ(CDATS,'(I8)') IDATS
          ENDIF
          IF (IHOUE.NE.24-ISHIFT) THEN
            CDUMMY=CDATEFC//'0000'
            CALL INCDATE(CDUMMY,-24*3600)
            CDE = CDUMMY(1:8)
            READ(CDE,'(I8)') IDE
          ELSE
            IDE = IDATEFC
          ENDIF
          IF (IDATS.LE.IDE) THEN
            ILINE = ILINE + 1
            IF (ILINE.EQ.1) THEN
              IF (IDATS.LT.IDE) THEN
                WRITE(MARS(ILINE),'("RET,DATE=",I8.8,"/TO/",I8.8,
     &           ",TIME=00/TO/18,")')
     &           IDATS, IDE
              ELSE
                WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     &           ",TIME=00/TO/18,")')
     &           IDATS
              ENDIF
              ILINE = 7
              WRITE(MARS(ILINE), '("TARGET=""fort.",I2,"""")')
     &         IUNIT
            ELSE
              IF (IDATS.LT.IDE) THEN
                WRITE(MARS(ILINE),'("RET,DATE=",I8.8,"/TO/",I8.8,
     &           ",TIME=00/TO/18")')
     &           IDATS, IDE
              ELSE
                WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     &           ",TIME=00/TO/18")')
     &           IDATS
              ENDIF
            ENDIF
          ENDIF
          IF (IDATEFC.NE.IDE) THEN
            ILINE = ILINE + 1
            IF (ILINE.EQ.1) THEN
              IF (IHOUE.EQ.0) THEN
                WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     &           ",TIME=00,")') IDATEFC
              ELSE
                WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     &           ",TIME=00/TO/",I2.2,",")')
     &           IDATEFC, IHOUE
              ENDIF
              ILINE = 7
              WRITE(MARS(ILINE), '("TARGET=""fort.",I2,"""")')
     &         IUNIT
            ELSE
              IF (IHOUE.EQ.0) THEN
                WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     &           ",TIME=00")') IDATEFC
              ELSE
                WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     &           ",TIME=00/TO/",I2.2)')
     &           IDATEFC, IHOUE
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF

!*    3.2 RETRIEVAL OF FORECASTED FIELDS.
!         -------------------------------

      IF (IFOREPD.GT.0) THEN
        ILINE = ILINE + 1
        WRITE(MARS(ILINE),'("RET,DATE=",I8,",STEP=6/TO/",I3,
     &   "/BY/",I2,",")') IDATEFC,IFOREPD/3600,
     &   ISHIFT
        IF (ILINE.EQ.1) THEN
          WRITE(MARS(4),'("LEVTYPE=SFC,LEVELIST=OFF,")')
          ILINE = 7
          WRITE(MARS(ILINE),
     &     '("TIME=12,TYPE=FC,TARGET=""fort.",I2,"""")') IUNIT
        ELSE
          ILINE = ILINE + 1
          WRITE(MARS(ILINE),'("TIME=12,TYPE=FC")')
        ENDIF
      ENDIF

!*    4. OUTPUT MARS TO TAPE22 AND PRINTER.
!        ----------------------------------

      WRITE(IU06,*) ' '
      WRITE(IU06,*) ' THE MARS INPUT CARDS ARE AS FOLLOWS --'
      WRITE(IU06,*) '  '

      DO I=1,ILINE
        WRITE(IU22,'(A80)') MARS(I)
        WRITE(IU06,*) MARS(I)
      ENDDO
      CLOSE(UNIT=IU22,STATUS='KEEP')

!     5. OUTPUT FILENAME TO BE USED TO KEEP THE RESULT OF THE MARS
!        REQUEST
!     ------------------------------------------------------------

      WRITE(CANPERIOD,'(I3.3)') IANALPD/3600
      WRITE(CFCPERIOD,'(I3.3)') IFOREPD/3600
      WRITE(CSHIFT,'(I2.2)') ISHIFT

      CFILENAME='windfile_'//CGRID//'_'//CDAT1(1:10)//'_'//CDAT2(1:10)
     &//'_'//CANPERIOD//'_'//CFCPERIOD//'_'//CSHIFT 

      open(17,file='windfilename')
      write(17,*) CFILENAME 
      close(17)

      STOP
      END


