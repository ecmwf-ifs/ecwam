      PROGRAM GETINPTB
C
C --------------------------------------------------------------------
C
C     PETER JANSSEN     ECMWF        FEBRUARI 1987
C     JEAN BIDLOT       ECMWF        MAY 1995  MIGRATION TO FUJITSU
C
C     PURPOSE
C     -------
C
C          THIS PROGRAM WRITES THE INPUT DATA CARDS FOR MARS.
C          THE PARAMETER SPECIFICATIONS FOR MARS ARE READ FROM FILE
C          WAMINFO.
C          ALSO THE FILE NAME FOR THE RESULTS OF THE MARS RETRIEVAL
C          IS OUTPUT IN WINDFILENAME AND SHOULD BE USED TO SAVE MARS
C          RETRIEVAL ON LINE.
C
C     INTERFACE.
C     ---------
C
C     *FILE DEFINITION*
C
C      IU04  - INFO FROM STATUSFILE WAMINFO
C      IU06  - PRINTER OUTPUT UNIT
C      IU22  - CONTAINS RETRIEVAL COMMANDS FOR MARS
C
C     PARAMETERS.
C     ----------
C
C     VARIABLE    TYPE      PURPOSE
C     --------    ----      -------
C
C     CDAT1       CHAR*12   BEGINDATE OF RUN
C     CDAT2       CHAR*12   ENDDATE OF RUN
C     IDATEFC     INTEGER   START DATE FORECAST RUN
C     IANALPD     INTEGER   PERIOD OF ANALYSIS
C     IFOREPD     INTEGER   FORECAST PERIOD
C     ISHIFT      INTEGER   TIMESTEP OF WINDFIELDS FROM MARS
C
C
C     EXTERNALS.
C     ---------
C
C     *INCDATE*  -  FROM WAMCRLIB INCREASES DATE TIME GROUP
C     *READSTA*  -  READ WAMINFO FILE.
C
C----------------------------------------------------------------------

      CHARACTER*12 CDAT1, CDAT2, CPBDT, CPEDT , CDATF, CRSDT, CLSDT
      CHARACTER*12 CDUMMY
      CHARACTER CDATS*8,CDATEFC*8,CDE*8,CGRID*3,CHOUS*2,CHOUE*2
      CHARACTER CANPERIOD*3,CFCPERIOD*3,CSHIFT*2
      CHARACTER MARS(21)*80,CFILENAME*44
      INTEGER ISTAT(3)
C
#include "comcard.h"
C
C ---------------------------------------------------------------------
C
C*    1. PRINT WAMINFO FILE
C        ------------------
C
      IU04 = 4 
      IU06 = 6
      IU22 = 22
      CARD = ' '
      CALL READSTA (IU04, CDAT1, CDAT2, IANALPD, IFOREPD, ISHIFT,
     1                    CRSDT, CLSDT, CPBDT, CPEDT, IRST, INF, ISTAT)
      WRITE(IU06,'(1H1,'' START OF GETINPT'',//,
     1       '' THE WAMINFO FILE IS AS FOLLOWS:'',/)')
      DO 1000 I=1,14
         WRITE(IU06,'(4X,A72)') CARD(I)
 1000 CONTINUE
       ISHIFT = ISHIFT / 3600
C
C ---------------------------------------------------------------------
C
C*    2. DETERMINE DATES AND HOURS.
C        --------------------------
C
      CDATS=CDAT1(1:8)
      READ(CDATS,'(I8)') IDATS
      CHOUS=CDAT1(9:10)
      READ(CHOUS,'(I2)') IHOUS
      IF      (IHOUS.GE.0  .AND. IHOUS.LT. 6) THEN
         IHOUS = 0
      ELSE IF (IHOUS.GE.6  .AND. IHOUS.LT.12) THEN
         IHOUS = 6
      ELSE IF (IHOUS.GE.12 .AND. IHOUS.LT.18) THEN
         IHOUS = 12
      ELSE IF (IHOUS.GE.18) THEN
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
      ELSE IF (IHOUE.GT.6  .AND. IHOUE.LE.12) THEN
         IHOUE = 12
      ELSE IF (IHOUE.GT.12 .AND. IHOUE.LE.18) THEN
         IHOUE = 18
      ELSE IF (IHOUE.GE.18) THEN
         IHOUE = 0
         CDUMMY=CDATEFC//'0000'
         CALL INCDATE (CDUMMY, 24*3600)
         CDATEFC = CDUMMY(1:8)
         READ(CDATEFC,'(I8)') IDATEFC
      ENDIF
C
C ---------------------------------------------------------------------
C
C*    3. INITIALISE UNIT TO RETRIEVE MARS DATA AND LINE COUNTER.
C        -------------------------------------------------------
C
      IUNIT = 80
      ILINE = 0
C
C ---------------------------------------------------------------------
C
C*    3. WRITE MARS COMMANDS.
C        --------------------
C
      MARS(2) = CARD(2)
      MARS(3) = CARD(3)
      WRITE(MARS(4),'("TYPE=AN,LEVTYPE=SFC,LEVELIST=OFF,")')

#if   resolution==100 && region=='m'
       WRITE(MARS(5),'("GRID=1.0/1.0,FORMAT=PACKED,")')
       WRITE(MARS(6),'("AREA=46./-6./30./42.0,")')
       CGRID='100'
#elif resolution==50 && region=='m'
       WRITE(MARS(5),'("GRID=0.5/0.5,FORMAT=PACKED,")')
       WRITE(MARS(6),'("AREA=66./-6./30./42.,")')
       CGRID='050'
#elif resolution==25 && region=='m'
       WRITE(MARS(5),'("GRID=0.25/0.25,FORMAT=PACKED,")')
       WRITE(MARS(6),'("AREA=81./-98./9./42.,")')
       CGRID='025'
#elif resolution==900 && region=='g'
       WRITE(MARS(5),'("GRID=9.0/9.0,FORMAT=PACKED,")')
       WRITE(MARS(6),'("AREA=72./0./-72./351.0,")')
       CGRID='900'
#elif resolution==300 && region=='g'
       WRITE(MARS(5),'("GRID=3.0/3.0,FORMAT=PACKED,")')
       WRITE(MARS(6),'("AREA=72./0./-63./357.0,")')
       CGRID='300'
#elif resolution==150 && region=='g'
       WRITE(MARS(5),'("GRID=1.5/1.5,FORMAT=PACKED,")')
       WRITE(MARS(6),'("AREA=81./0./-81./358.5,")')
       CGRID='150'
#elif resolution==50 && region=='g'
       WRITE(MARS(5),'("GRID=0.5/0.5,FORMAT=PACKED,")')
       WRITE(MARS(6),'("AREA=81./0./-81./359.5,")')
       CGRID='050'
#elif resolution==25 && region=='g'
       WRITE(MARS(5),'("GRID=0.25/0.25,FORMAT=PACKED,")')
       WRITE(MARS(6),'("AREA=81./0./-81./359.75,")')
       CGRID='025'
#else
C  NO PARAMETER FILE SELECTED
C "ICEDGE"     = iceedge
C "REGION"     = region
C "BLOCKS"     = blocks
C "RESOLUTION" = resolution
C "DIRECTIONS" = directions
#endif

      IF (IANALPD.GE.0) THEN
C
C*    3.1 RETRIEVAL OF ANALYSIS.
C         ----------------------
C

C    ALL FIELDS ARE ON ONE DAY.

         IF (IDATS.EQ.IDATEFC) THEN
            ILINE = ILINE + 1
            IF (ILINE.EQ.1) THEN
               IF (IHOUS.EQ.IHOUE) THEN
                  WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     1                                ",TIME=",I2.2,",")')
     2                                  IDATS, IHOUS
               ELSE
                  WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     1                                ",TIME=",I2.2,"/TO/",I2.2,",")')
     2                                  IDATS, IHOUS, IHOUE
               ENDIF
               WRITE(MARS(7), '("TARGET=""fort.",I2,"""")') IUNIT
               ILINE = 7
            ELSE
               IF (IHOUS.EQ.IHOUE) THEN
                  WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     1                                ",TIME=",I2.2)')
     2                                  IDATS, IHOUS
               ELSE
                  WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     1                                ",TIME=",I2.2,"/TO/",I2.2)')
     2                                  IDATS, IHOUS, IHOUE
               ENDIF
            ENDIF
         ELSE

C    MORE THAN ONE DAY

            IF (IHOUS.NE.0) THEN
               ILINE = ILINE + 1
               IF (ILINE.EQ.1) THEN
                  IF (IHOUS.EQ.24-ISHIFT) THEN
                     WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     1                                   ",TIME=",I2.2,",")')
     2                                     IDATS, IHOUS
                  ELSE
                     WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     1                          ",TIME=",I2.2,"/TO/",I2.2,",")')
     2                                     IDATS, IHOUS, 24-ISHIFT
                  ENDIF
                  ILINE = 7
                  WRITE(MARS(ILINE), '("TARGET=""fort.",I2,"""")')
     1                               IUNIT
               ELSE
                  IF (IHOUS.EQ.24-ISHIFT) THEN
                     WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     1                                   ",TIME=",I2.2)')
     2                                     IDATS, IHOUS
                  ELSE
                     WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     1                          ",TIME=",I2.2,"/TO/",I2.2)')
     2                                     IDATS, IHOUS, 24-ISHIFT
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
     1                                ",TIME=00/TO/18,")')
     2                                  IDATS, IDE
                  ELSE
                     WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     1                                ",TIME=00/TO/18,")')
     2                                  IDATS
                  ENDIF
                  ILINE = 7
                  WRITE(MARS(ILINE), '("TARGET=""fort.",I2,"""")')
     1                               IUNIT
               ELSE
                  IF (IDATS.LT.IDE) THEN
                     WRITE(MARS(ILINE),'("RET,DATE=",I8.8,"/TO/",I8.8,
     1                                ",TIME=00/TO/18")')
     2                                  IDATS, IDE
                  ELSE
                     WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     1                                ",TIME=00/TO/18")')
     2                                  IDATS
                  ENDIF
               ENDIF
            ENDIF
            IF (IDATEFC.NE.IDE) THEN
               ILINE = ILINE + 1
               IF (ILINE.EQ.1) THEN
                  IF (IHOUE.EQ.0) THEN
                     WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     1                                   ",TIME=00,")') IDATEFC
                  ELSE
                     WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     1                                   ",TIME=00/TO/",I2.2,",")')
     2                                    IDATEFC, IHOUE
                  ENDIF
                  ILINE = 7
                  WRITE(MARS(ILINE), '("TARGET=""fort.",I2,"""")')
     1                   IUNIT
               ELSE
                  IF (IHOUE.EQ.0) THEN
                     WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     1                                   ",TIME=00")') IDATEFC
                  ELSE
                     WRITE(MARS(ILINE),'("RET,DATE=",I8.8,
     1                                   ",TIME=00/TO/",I2.2)')
     2                                    IDATEFC, IHOUE
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
C
C*    3.2 RETRIEVAL OF FORECASTED FIELDS.
C         -------------------------------
C
      IF (IFOREPD.GT.0) THEN
         ILINE = ILINE + 1
         WRITE(MARS(ILINE),'("RET,DATE=",I8,",STEP=6/TO/",I3,
     1                       "/BY/",I2,",")') IDATEFC,IFOREPD/3600,
     2                        ISHIFT
         IF (ILINE.EQ.1) THEN
             WRITE(MARS(4),'("LEVTYPE=SFC,LEVELIST=OFF,")')
             ILINE = 7
             WRITE(MARS(ILINE),
     1          '("TIME=12,TYPE=FC,TARGET=""fort.",I2,"""")') IUNIT
         ELSE
            ILINE = ILINE + 1
            WRITE(MARS(ILINE),'("TIME=12,TYPE=FC")')
         ENDIF
      ENDIF
C
C*    4. OUTPUT MARS TO TAPE22 AND PRINTER.
C        ----------------------------------
C
      ILINE = ILINE + 1
      MARS(ILINE)='END'

      WRITE(IU06,*) ' '
      WRITE(IU06,*) ' THE MARS INPUT CARDS ARE AS FOLLOWS --'
      WRITE(IU06,*) '  '

      DO 77 I=1,ILINE
         WRITE(IU22,'(A80)') MARS(I)
         WRITE(IU06,*) MARS(I)
  77  CONTINUE
      CLOSE(UNIT=IU22,STATUS='KEEP')
C
C     5. OUTPUT FILENAME TO BE USED TO KEEP THE RESULT OF THE MARS
C        REQUEST
C     ------------------------------------------------------------

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
