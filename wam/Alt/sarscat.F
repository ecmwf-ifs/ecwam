C ------------------------------------------------------------------
C
      PROGRAM SARSCAT
C
C ------------------------------------------------------------------
C
C**** *sarscat*
C
C
C      b. hansen         ecmwf      july    1992
C
C     PURPOSE.
C     --------
C     COLLECTION OF PERIODICAL STATISTICS ON ERS-1 DATA.
C     1. PRODUCT URA:
C        SCATTER PLOTS OF ERS-1 WAVE HEIGHTS OVER WAM-MODEL WAVE HEIGHTS
C        AND OF ERS-1 WINDSPEEDS OVER Txxx-MODEL WINDSPEEDS.
C        DATA INPUT IS THE COLLOCATION FILE RACyymmddhh, ALL FILES
C        HAVE TO BE ONLINE i.e. ECFILE CALL (IF REQUIRED HAVE
C        TO BE DONE BEFORE. TO RUN THIS PROGRAM USE THE FOLLOWING
C        COMMAND:
C        ls RAC* | urascat   
C        WHERE urascat IS THE NAME OF THE EXECUTABLE OF THIS SOURCE.
C        INTERMEDIATE DATA PROCESSING IS PROJECTED ON A GRID. THERE-
C        FORE THE GRID INFORMATION OF THE CORRESPONDING WAModel GRID
C        NEEDS TO BE PROVIDED i.e. gridglou MUST BE PRESENT AS FILE
C        GRID.
C
C**   INTERFACE.
C     ----------
C          NONE
C
C**   LIBRARIES.
C     ----------
C          $ECLIB
C          $MAGLIB
C          $EMOSLIB
C
C**   EXTERNALS.
C     ----------
C          POPEN      INITIATES MAGISC.
C          PSETC      SINGLE CHARACTER PARAMETER SETTING.
C          SCATIN     READS AND EVALUATES THE USERINPUT.
C          ABORT      TERMINATION ROUTINE IN ERROR CASE.
C          INDLOLE    CONVERTS CORNER COORDINATIONS TO GRID INDICES.
C          INDSEEK    CONVERTS LAT/LONG COORDINATES TO GRID INDICES.
C          STATSI     INITIALIZES SUMS.
C          STATSE     EVALUATES STATISTICS.
C          STATSP     PRINTS STATISTICS.
C          PLTCONTW   PERFORMS SCATTER PLOTT.
C          TITMAK     GENERATES HEADLINES FOR MAGICS.
C          TEXTBOX    GENERATES ADJUSTED TWO COLUMNS OF TEXT.
C          PNEW       EJECTS PAGE INSIDE MAGICS.
C          PCLOSE     TERMINATES MAGICS.
C
C**   INTERFACE.
C     ----------
C          STANDARD OUT ON REQUEST IS WRITTEN TO jpuso.
C          ALL DATA ARE READ FROM THE FILENAME GIVEN TO THIS
C          PROGRAM INTO STANDARD INPUT VIA A PIPE OR A SEPERATE
C          LIST OF FILENAMES.
C          
C
C     METHOD.
C     -------
C          NONE.
C
C     REFERENCE.
C     ----------
C          NONE.
C
C ----------------------------------------------------------------------
C
C
C*    0. DEFINITIONS.
C        ------------
      PARAMETER (
     . jpuso    =   6,   !  log unit for standard output
     . jpuin    = 'INP', !  log unit for user input     
     . jpumer   =  23,   !  log unit for file: RAC
     . jpuwast  =  46,   !  log unit for statistics on file: RAC wave
     . jpuwist  =  47,   !  log unit for statistics on file: RAC wind
     . jpuphs   =  48,   !  log unit for statistics to be plotted wave
     . jpupu10  =  49,   !  log unit for statistics to be plotted wind
     . jpust    =  50,   !  log unit to write namelists to
     . jpuari   =  51,   !  log unit to write scatterfields to
     . jpusum   =  52,   !  log unit to write gridpoint sums to
     . jpngx    = 120,   !  size of model grid along latitudes
     . jpngy    =  61,   !  size of model grid along longitudes
     . jpns     =  10,   !  dimension of array sum and sumwi
     . jps      =  10,   !  dimension of array ph and wiscfi
     . jpar     =  10,   !  dimension of array ka1 - ppo
     . jpnx     = 101,   !  x-size of scatterfield of waveheights
     . jpny     = 101,   !  y-size of scatterfield of waveheights
     . ppdx     =  0.25, !  x-increment  -----   "    ----------
     . ppdy     =  0.25, !  y-increment  -----   "    ----------
     . ppxa     =  0.0 , !  x-origen     -----   "    ----------
     . ppya     =  0.0 , !  y-orgen      -----   "    ----------
     . ppymaxwa =  8.0 , !  maximum y-value of scatterdiagram ----
     . ppamaxwa =  8.0 , !  maximum x-value of scatterdiagram ----
     . ppatintwa = 2.0 , !  x-tick-interval of scatterplot   -----
     . ppytintwa = 2.0 , !  y-tick-interval of scatterplot   -----
     . jpwindnx = 81   , !  x-size of scatterfield of windspeeds
     . jpwindny = 81   , !  y-size of scatterfield of windspeeds
     . ppwinddx =  0.5 , !  x-increment  -----   "    --------
     . ppwinddy =  0.5 , !  y-increment  -----   "    --------
     . ppwindxa =  0.0 , !  x-origen     -----   "    --------
     . ppwindya =  0.0 , !  y-orgen      -----   "    --------
     . ppaminwi =  0.0 , !  x-minimum value of scatterdiagram --
     . ppyminwi =  0.0 , !  y-minimum value of scatterdiagram --
     . ppamaxwi = 25.0 , !  x-maximum value of scatterdiagram --
     . ppymaxwi = 25.0 , !  y-maximum value of scatterdiagram --
     . ppatintwi = 5.0 , !  x-tick-interval of scatterplot  ----
     . ppytintwi = 5.0 , !  y-tick-interval of scatterplot  ----
     . jpkl     = 10   , !  dimension of headlines array yoline
     . jpkv     =  8     !  dimension of values array pval
     .          )

      INTEGER    ka1(jpar), ka2(jpar),
     .           ko1(jpar), ko2(jpar),
     .           la1(jpar), la2(jpar),
     .           lo1(jpar), lo2(jpar)

      REAL       pwa(jpar), pwo(jpar), ppa(jpar), ppo(jpar), pval(jpkv)

      REAL       swhmean(jps),swherro(jps),swhbias(jps)
      REAL       u10mean(jps),u10erro(jps),u10bias(jps)

      REAL       sum  (jpngx, jpngy, jpns)
      REAL       ph   (jpnx,  jpny,  jps )
      REAL       sumwi (jpngx,    jpngy,    jpns)
      REAL       wiscfi(jpwindnx, jpwindny, jps )

      DIMENSION  isym(2)
      CHARACTER*  1 yooption
      CHARACTER*  2 yoyear
      CHARACTER*  3 yousd
      CHARACTER*  3 yomonth(14)
      CHARACTER*  4 yofid
      CHARACTER*  7 yor(4)
      CHARACTER* 30 yoaxtite(2)
      CHARACTER* 30 yocom          ! CURRENT DATA FILE NAME
      CHARACTER* 36 yotitle1(2)
      CHARACTER* 36 yotitle2(2)
      CHARACTER* 60 youserli       ! PAGE_ID_LINE_USER_TEXT
      CHARACTER* 60 yoregion(10)
      CHARACTER*120 yohead, yosubhead, yoperiod, yosubreg,
     .              yoextra, yoreg, yortxt
      CHARACTER*130 yoline(jpkl)
      LOGICAL    loabort, lotwice(jps), lofailed, lohead
C
C ----------------------------------------------------------------------
C
      INTEGER rtb
C
      INTEGER f_getenv
      CHARACTER*256 YOISCOM
      CHARACTER*  3 YOENV  
C
C ----------------------------------------------------------------------
C
C*       PARAMETER NAMELIST AS IN COMMON PARAM OF THE WAM-MODEL.
C        -------------------------------------------------------
C
      LOGICAL loper          !  the land sea mask is periodic (T)
      NAMELIST /PARWAM/
     .  NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER, MOUTP, MOUTT,
     .  NMAXC, NMAXF, NBINP, NIBL1, NIBLD, NBLD, NIBLC, NBLC ,
     .  ITAUMAX, JUMAX, IUSTAR, IALPHA, NDEPTH, IREFRA, loper
C
C ----------------------------------------------------------------------
C
      NAMELIST /wave/yotitle1, yoperiod, yoreg,
     .             NTOT,
     .             XMAX, YMAX, XMEAN, YMEAN, STDEVX, STDEVY,
     .             DMAX, BIAS, RMSE, STDEVD, SID, RV, CORR,
     .             SLOPES, RMSESS, RMSEUS, SLOPE1, RMSES1, RMSEU1,
     .             SLOPE2, RINTR2, RMSES2, RMSEU2

      NAMELIST /wind/yotitle2, yoperiod, yoreg,
     .             nwindtot,
     .             xwindmax, ywindmax, XMEAN, YMEAN, STDEVX, STDEVY,
     .             dwindmax, BIAS, RMSE, STDEVD, SID, RV, CORR,
     .             SLOPES, RMSESS, RMSEUS, SLOPE1, RMSES1, RMSEU1,
     .             SLOPE2, RINTR2, RMSES2, RMSEU2
C            
C ----------------------------------------------------------------------
C
      DATA yomonth /'jan','feb','mar','apr',
     .              'may','jun','jul','aug',
     .              'sep','oct','nov','dec',
     .              '_Hs','u10' /
C
C ----------------------------------------------------------------------
C
C*    1. INITIALIZATION.
C        ---------------
C
 1000 CONTINUE
      isym(1)=18
      isym(2)=2
      yoline(1) =' \UL\STATISTICS '
      yoline(2) ='             '
      yoline(3) =' ENTRIES          '
      yoline(4) =' MEAN WAM         '
      yoline(5) =' MEAN ERS-1       '
      yoline(6) =' BIAS (ERS-1 - WAM) '
      yoline(7) =' STANDARD DEVIATION '
      yoline(8) =' SCATTER INDEX    '
      yoline(9) =' CORRELATION      '
      yoline(10)=' SYMMETRIC SLOPE  '
      youserli  = 'H>ansen,<B.'
      yosubhead = ' SAR - <W>aveheights        '
      yoextra   = ' '
      yohead    =  ' WAM-M>odel  -  <ERS-1  -  C>omparison'
      yotitle1(1) = 'WAM-Model - ERS-1  -  Comparison'
      yotitle1(2) = 'SAR - Waveheights'
      yotitle2(1) = ' '
      yotitle2(2) = ' '
      yosubreg  = ' first guess ocean wave spectra from <ECMWF>'//
     .            ' e-suite <WAM>odel (<ERS1>-ALT-DATA-ASSIMILATION) '
      yortxt      ='REGION: '//
     .'>lower left longitude   lower left latitude   upper right lo'//
     .'ngitude   upper right latitude'

      LOFAILED  = .FALSE.
      jpfail    =   0
      ient      =   0
      katime    = -10
C
      DO 1100 jngx=1,jpngx
        DO 1100 jngy=1,jpngy
          DO 1100 jx=1,jpns
            sum  (jngx,jngy,jx) = 0.0
            sumwi(jngx,jngy,jx) = 0.0
 1100 CONTINUE
C
      DO 1200 jx=6,8
        DO 1200 jngy=1,jpngy
          DO 1200 jngx=1,jpngx
            sum  (jngx,jngy,jx) = -999999.
            sumwi(jngx,jngy,jx) = -999999.
 1200 CONTINUE
C
      DO 1300 jx=1,jps
        DO 1300 jngy=1,jpny
          DO 1300 jngx=1,jpnx
            ph(jngx,jngy,jx) = 0.0
 1300 CONTINUE
C
      DO 1400 jx=1,jps
        DO 1400 jngy=1,jpwindny
          DO 1400 jngx=1,jpwindnx
            wiscfi(jngx,jngy,jx) = 0.0
 1400 CONTINUE
C
C* 1.5   READ WAMODEL-PARAMETERS FROM FILE PARWAM.
C  -----------------------------------------------
C
      rewind 'PARWAM'
      READ('PARWAM',PARWAM)
C
C  1.6 DATA FILES WITH HEADER.
C  ---------------------------
C
      lohead=.F.
      yoiscom='HEADER'
      ir=f_getenv(yoiscom,yoenv)
      IF ( YOENV.EQ.".T.")lohead=.T.
      WRITE(jpuso,*)'  HEADER ?  >',lohead,'<'
C
C ----------------------------------------------------------------------
C
C*    2. READ USER INPUT.
C        ----------------
C
 2000 CONTINUE
      loabort = .FALSE.
      CALL scatin (jpuso, jpuin, iareas, ppa, pwo, pwa, ppo, yoregion,
     .             yooption)
C
C ----------------------------------------------------------------------
C
C*    3. CONVERT CORNER COORDINATES TO INDICES.
C        --------------------------------------
C
 3000 CONTINUE
      DO 3100 jre = 1,iareas
        CALL indlole(jpuso,  
     .     nang, nfre, ngx, ngy, nblo, niblo, nover, moutp, moutt,
     .     nmaxc, nmaxf, nbinp, nibl1, nibld, nbld, niblc, nblc ,
     .     itaumax, jumax, iustar, ialpha, ndepth, irefra, loper,
     .     ppa(jre), pwo(jre), pwa(jre), ppo(jre),
     .     ka1(jre), ko1(jre), la1(jre), lo1(jre), lotwice(jre),
     .     ka2(jre), ko2(jre), la2(jre), lo2(jre), loabort,
     .     amonop, amoeap, amosop, amowep,
     .     delat, delong, ngx, ngy, iper)
        IF (loabort) CALL ABORT
        write(jpuso,'(a20)') ' '
 3100 CONTINUE
      PRINT *, ' PASSED INDLOLE ', JRE-1 ,' TIMES'
      PRINT *, ' AMOWEP ', AMOWEP
      PRINT *, ' AMOSOP ', AMOSOP
      PRINT *, ' AMONOP ', AMONOP
      PRINT *, ' AMOEAP ', AMOEAP
      PRINT *, ' DELAT  ', DELAT
      PRINT *, ' DELONG ', DELONG
      PRINT *, ' NGX    ', NGX
      PRINT *, ' NGY    ', NGY
      PRINT *, ' IPER   ', IPER
C
C ----------------------------------------------------------------------
C
C*    4. LOOP OVER PRODUCTS TO CALCULATE GRID-POINT-SUMS.
C        ------------------------------------------------
C
 4000 CONTINUE
      CALL popen
      CALL psetc ('PAGE_FRAME',                  'ON' )
      CALL psetc ('PAGE_ID_LINE_USER_TEXT', youserli  )
      CALL psetc ('SUBPAGE_FRAME',              'OFF' )
C
C*   4.1 OPEN THE NEXT DATA FILE AS FROM STDIN ON UNIT 'jpumer'.
C        -------------------------------------------------------
C
 4100 CONTINUE
      DO 4501 jv = 1,1000                 ! dummy loop 
       lofailed = .FALSE.                 !
       READ (*,'(a30)',END=5000) yocom    ! processing of all data 
C                                           is finished when all
C                                           lines are read from stdin.
       lfp=rtb(yocom)
       kios = 0
       OPEN(jpumer,FILE=yocom(1:lfp),STATUS='old',
     .    IOSTAT=kios, FORM='unformatted')
       IF ( kios .NE. 0) THEN
         WRITE(jpuso,'(a50)')' product '//yocom(1:lfp)//
     .                       ' not opened  '
         WRITE(jpuso,'(a40,i10)')' reason is ',kios
         lofailed = .TRUE.
       ELSE
         IF (lohead) THEN
            READ (jpumer,END=4500) ztime
            WRITE(jpuso,*)' ZTIME = ',ztime
            ktime=NINT(ztime)
         ENDIF
       ENDIF
C
C*   4.2 READ FILE TO ADD SUMMS.
C        -----------------------
C
 4200 CONTINUE
       IF (.NOT. lofailed ) THEN
        DO 4300 jm = 1,1000
C         READ(jpumer,END=4400)ktime, inul1, rlat, rlot, rnul1,
C    .                         rnul2, inul2, rwhe, rwhw,
C    .                         zwinders,zwindwam,xt

          READ(jpumer,END=4400)rlot,rlat,xdatw
          READ(jpumer,END=4400)rwhw
          READ(jpumer,END=4400)xlons,xlats,xdats
          READ(jpumer,END=4400)rwhe
          READ(jpumer,END=4400)xvssi
          READ(jpumer,END=4400)xvsbe
          READ(jpumer,END=4400)xvsob
          READ(jpumer,END=4400)track
          READ(jpumer,END=4400)jpksp
          READ(jpumer,END=4400)afl
          READ(jpumer,END=4400)xefre
          READ(jpumer,END=4400)xsdata
C
C    rwhw --> x-axis
C    rwhe --> y-axis
C
          IF(rwhw .GE.  0.0 .AND. rwhe .GE.  0.0) THEN
           LOFAILED=.FALSE.
           CALL indseek(rlat,rlot,ila,ilo,lofailed,
     .                  amonop, amoeap, amosop, amowep,
     .                  delat, delong, ngx, ngy, iper)
C          PRINT4201,' HS-M: ', rwhw,' HS-S: ', rwhe,
C    .               '  LAT: ',rlat,' INDEX: ',ila,
C    .               '  LON: ',rlot,' INDEX: ',ilo
 4201      FORMAT(2(a7,f5.2),2(a7,f8.3,a8,i5))

           IF (.NOT.LOFAILED) THEN
             IF (katime .EQ. -10)THEN
                katime=NINT(xdats)*100
             ENDIF
C
C*   4.2.1 WAVE SUMS.
C          ----------
C
 4210 CONTINUE
             sum(ilo,ila,1)=sum    (ilo,ila,1) + rwhw
             sum(ilo,ila,2)=sum    (ilo,ila,2) + rwhe
             sum(ilo,ila,3)=sum    (ilo,ila,3) + rwhw * rwhw
             sum(ilo,ila,4)=sum    (ilo,ila,4) + rwhe * rwhe
             sum(ilo,ila,5)=sum    (ilo,ila,5) + rwhw * rwhe
             sum(ilo,ila,6)=MAX(sum(ilo,ila,6), rwhw)
             sum(ilo,ila,7)=MAX(sum(ilo,ila,7), rwhe)
             sum(ilo,ila,8)=MAX(sum(ilo,ila,8), ABS(rwhe-rwhw))
             sum(ilo,ila,9)=sum    (ilo,ila,9)+(rwhe-rwhw)*(rwhe-rwhw)
             sum(ilo,ila,10)=sum   (ilo,ila,10)+1.
C
             ii = nint (( rwhw - ppxa )/ ppdx + 1.)
             ik = nint (( rwhe - ppya )/ ppdy + 1.)
             IF ( ii .GT. jpnx  .OR. ik .GT. jpny .OR.
     .            ii .EQ. 0 .OR. ik .EQ. 0 ) THEN
                WRITE(jpuso,*) RWHW,PPXA,PPYA,RWHE,PPDY,PPDX,II,IK
             ELSE
                DO 4211 jars=1,iareas
                   IF (lotwice(jars) ) THEN
                     IF ( ila .LE. la1(jars) .AND.
     .                    ila .GE. ka1(jars) .AND.
     .                   (ilo .GE. lo1(jars) .OR.
     .                    ilo .LE. ko2(jars))) THEN
                          ph(ii,ik,jars) = 1. + ph(ii,ik,jars)
                     ENDIF
                   ELSE
                     IF ( ila .LE. la1(jars) .AND.
     .                    ila .GE. ka1(jars) .AND.
     .                    ilo .GE. lo1(jars) .AND.
     .                    ilo .LE. ko1(jars)) THEN
                          ph(ii,ik,jars) = 1. + ph(ii,ik,jars)
                     ENDIF
                   ENDIF
 4211           CONTINUE
             ENDIF
             GOTO 4225 
C
C*   4.2.2 WIND SUMS.
C          ----------
C
 4220 CONTINUE
             sumwi(ilo,ila,1)=sumwi    (ilo,ila,1) + zwindwam
             sumwi(ilo,ila,2)=sumwi    (ilo,ila,2) + zwinders
             sumwi(ilo,ila,3)=sumwi    (ilo,ila,3) + zwindwam * zwindwam
             sumwi(ilo,ila,4)=sumwi    (ilo,ila,4) + zwinders * zwinders
             sumwi(ilo,ila,5)=sumwi    (ilo,ila,5) + zwindwam * zwinders
             sumwi(ilo,ila,6)=MAX(sumwi(ilo,ila,6), zwindwam)
             sumwi(ilo,ila,7)=MAX(sumwi(ilo,ila,7), zwinders)
             sumwi(ilo,ila,8)=
     .        MAX(sumwi(ilo,ila,8), ABS(zwinders-zwindwam))
             sumwi(ilo,ila,9)=sumwi(ilo,ila,9)+(zwinders-zwindwam)
     .        *(zwinders-zwindwam)
             sumwi(ilo,ila,10)=sumwi   (ilo,ila,10)+1.
C
             ii = nint (( zwindwam - ppwindxa )/ ppwinddx + 1.)
             ik = nint (( zwinders - ppwindya )/ ppwinddy + 1.)
C
             WRITE(1,8000)ktime, inul1, rlat, rlot, rnul1,
     .                    rnul2, inul2, rwhe, rwhw,
     .                    zwinders, zwindwam, xt, ii, ik
C
             IF ( ii .GT. jpwindnx  .OR. ik .GT. jpwindny .OR.
     .            ii .EQ. 0 .OR. ik .EQ. 0 ) THEN
                WRITE(jpuso,*) zwindwam,ppwindxa,ppwindya,
     .                         zwinders,ppwinddy,ppwinddx,II,IK
             ELSE
                DO 4221 jars=1,iareas
                   IF (lotwice(jars) ) THEN
                     IF ( ila .LE. la1(jars) .AND.
     .                    ila .GE. ka1(jars) .AND.
     .                   (ilo .GE. lo1(jars) .OR.
     .                    ilo .LE. ko2(jars))) THEN
                          wiscfi(ii,ik,jars) = 1. + wiscfi(ii,ik,jars)
                     ENDIF
                   ELSE
                     IF ( ila .LE. la1(jars) .AND.
     .                    ila .GE. ka1(jars) .AND.
     .                    ilo .GE. lo1(jars) .AND.
     .                    ilo .LE. ko1(jars)) THEN
                         wiscfi(ii,ik,jars) = 1. + wiscfi(ii,ik,jars)
                     ENDIF
                   ENDIF
 4221           CONTINUE
             ENDIF
 4225        CONTINUE
           ELSE
             WRITE(jpuso,*)'  found values at Latitude: ',rlat,
     .       ' and Longitude: ',rlot,' ',yocom(1:lfp), jm
           ENDIF
          ENDIF
C
C*   4.3 FIN TO READ FILE 'jpumer' IF LONGER THAN 'jm'.
C        ----------------------------------------------
 4300    CONTINUE
         WRITE(jpuso,'(a55)')' not all data read in '
         GOTO 4500
C
C*   4.4 FIN TO READ FILE 'jpumer' IF SHORTER THAN 'jm'.
C        -----------------------------------------------
C
 4400    CONTINUE
         WRITE(jpuso,'(a55)')'     all data read in '
       ENDIF
C
C*   4.5 END OF PROCESSING OF ONE PRODUCT.
C        ---------------------------------
C
 4500 CONTINUE
      CLOSE(jpumer,STATUS='keep')
 4501 CONTINUE
C-    WRITE(jpuso,'(50f3.0)')wiscfi
C
C ----------------------------------------------------------------------
C
C*    5. LOOP OVER REQUESTED AREAS.
C        --------------------------
C
 5000 CONTINUE
      ktime=NINT(xdats)*100
      CALL titmak(katime,ktime,imon,iyear,yoperiod)
      WRITE(yoyear,'(i2)')MOD(iyear,100)
      open(jpuwast,FILE='RAC'//yomonth(imon)//yoyear//'WAVE')
C     open(jpuwist,FILE='RAC'//yomonth(imon)//yoyear//'WIND')
      open(jpuphs ,FILE='RAC'//yomonth(  13)//'PLOT',position='append')
C     open(jpupu10,FILE='RAC'//yomonth(  14)//'PLOT',position='append')
      open(jpust,  FILE='RAC'//yomonth(imon)//yoyear//'LIST')
      open(jpuari, FILE='RAC'//yomonth(imon)//yoyear//'SCAT')
      open(jpusum, FILE='RAC'//yomonth(imon)//yoyear//'SUMS',
     .     form='unformatted')
      BUFFER OUT (jpusum,0) (sum(1,1,1),sum(jpngx,jpngy,jpns))
      BUFFER OUT (jpusum,0) (sumwi(1,1,1),sumwi(jpngx,jpngy,jpns))

      DO 6000 jre = 1, iareas
C       print*,' THIS IS LOOP NO <',jre,'>'
        yoreg=' '
        istrl=rtb(yoregion(jre))
        IF ( yoregion(jre) .EQ. ' ') THEN
          write(yor(1),'(f7.2)')pwo(jre)
          write(yor(2),'(f7.2)')pwa(jre)
          write(yor(3),'(f7.2)')ppo(jre)
          write(yor(4),'(f7.2)')ppa(jre)
          yoreg = '                    '  // yor(1) //
     .                    '            '  // yor(2) //
     .                '                '  // yor(3) //
     .                 '               '  // yor(4)
        ELSE
          yoreg(1:istrl)=yoregion(jre)
        ENDIF
C
C*   5.2 ADD SUMS FOR REGION.
C        --------------------
C
 5200   CONTINUE
C
C*   5.2.1 WAVE-SUMS.
C          ----------
C
 5210   CONTINUE
C
C*   5.2.1.1 INITIALIZE SUMS.
C            ----------------
C
 5211   CONTINUE
        CALL STATSI (NTOT, XMAX, YMAX, DMAX,
     .               SUMX, SUMY, SUMXX, SUMYY, SUMXY, SUMXY2)
        xtot = 0.0
C
C*   5.2.1.2 BUILD SUMS.
C            -----------
C
 5212   CONTINUE
        DO 5213 jy = ka1(jre),la1(jre)
        DO 5213 jx = lo1(jre),ko1(jre)
          xtot = sum(jx,jy,10) + xtot
          xmax = MAX(sum(jx,jy, 6), xmax)
          ymax = MAX(sum(jx,jy, 7), ymax)
          dmax = MAX(sum(jx,jy, 8), dmax)
          sumx = sum(jx,jy, 1) + sumx
          sumy = sum(jx,jy, 2) + sumy
          sumxx = sum(jx,jy, 3) + sumxx
          sumyy = sum(jx,jy, 4) + sumyy
          sumxy = sum(jx,jy, 5) + sumxy
          sumxy2 = sum(jx,jy, 9) + sumxy2
 5213   CONTINUE
C
        IF (lotwice(jre)) THEN
          DO 5214 jy = ka2(jre),la2(jre)
          DO 5214 jx = lo2(jre),ko2(jre)
            xtot = sum(jx,jy,10) + xtot
            xmax = MAX(sum(jx,jy, 6), xmax)
            ymax = MAX(sum(jx,jy, 7), ymax)
            dmax = MAX(sum(jx,jy, 8), dmax)
            sumx = sum(jx,jy, 1) + sumx
            sumy = sum(jx,jy, 2) + sumy
            sumxx = sum(jx,jy, 3) + sumxx
            sumyy = sum(jx,jy, 4) + sumyy
            sumxy = sum(jx,jy, 5) + sumxy
            sumxy2 = sum(jx,jy, 9) + sumxy2
 5214     CONTINUE
        ENDIF
C
C*   5.2.2 WIND-SUMS.
C          ----------
C
 5220   CONTINUE
C
C*   5.2.2.1 INITIALIZE SUMS.
C            ----------------
C
 5221   CONTINUE
        CALL STATSI (nwindtot, xwindmax, ywindmax, dwindmax,
     .             sumwiX, sumwiY, sumwiXX, sumwiYY, sumwiXY, sumwiXY2)
        xwindtot = 0.0
C
C*   5.2.2.2 BUILD SUMS.
C            -----------
C
 5222   CONTINUE
        DO 5223 jy = ka1(jre),la1(jre)
        DO 5223 jx = lo1(jre),ko1(jre)
          xwindtot = sumwi(jx,jy,10) + xwindtot
          xwindmax = MAX(sumwi(jx,jy, 6), xwindmax)
          ywindmax = MAX(sumwi(jx,jy, 7), ywindmax)
          dwindmax = MAX(sumwi(jx,jy, 8), dwindmax)
          sumwix = sumwi(jx,jy, 1) + sumwix
          sumwiy = sumwi(jx,jy, 2) + sumwiy
          sumwixx = sumwi(jx,jy, 3) + sumwixx
          sumwiyy = sumwi(jx,jy, 4) + sumwiyy
          sumwixy = sumwi(jx,jy, 5) + sumwixy
          sumwixy2 = sumwi(jx,jy, 9) + sumwixy2
 5223   CONTINUE
C
        IF (lotwice(jre)) THEN
          DO 5224 jy = ka2(jre),la2(jre)
          DO 5224 jx = lo2(jre),ko2(jre)
            xwindtot = sumwi(jx,jy,10) + xwindtot
            xwindmax = MAX(sumwi(jx,jy, 6), xwindmax)
            ywindmax = MAX(sumwi(jx,jy, 7), ywindmax)
            dwindmax = MAX(sumwi(jx,jy, 8), dwindmax)
            sumwix = sumwi(jx,jy, 1) + sumwix
            sumwiy = sumwi(jx,jy, 2) + sumwiy
            sumwixx = sumwi(jx,jy, 3) + sumwixx
            sumwiyy = sumwi(jx,jy, 4) + sumwiyy
            sumwiyy = sumwi(jx,jy, 4) + sumwiyy
            sumwixy = sumwi(jx,jy, 5) + sumwixy
            sumwixy2 = sumwi(jx,jy, 9) + sumwixy2
 5224     CONTINUE
        ENDIF
C
        ntot = NINT(xtot)
        nwindtot = NINT(xwindtot)
C
C*   5.3 PROCESS SUMS.
C        -------------
C
 5300   CONTINUE
C
C*  5.3.1 WAVE SUMS.
C        -----------
C
 5310   CONTINUE
C*  5.3.1.1 COMPUTE STATISTICS.
C           -------------------
C
 5311 CONTINUE
        CALL STATSE
     .        (NTOT, SUMX, SUMY, SUMXX, SUMYY, SUMXY, SUMXY2,
     .         XMAX, YMAX, XMEAN, YMEAN, STDEVX, STDEVY,
     .         DMAX, BIAS, RMSE, STDEVD, SID, RV, CORR,
     .         SLOPES, RMSESS, RMSEUS, SLOPE1, RMSES1, RMSEU1,
     .         SLOPE2, RINTR2, RMSES2, RMSEU2)
C
C*   5.3.1.2 PRINT WAVE STATISTICS.
C            ----------------------
C
 5312   CONTINUE
        ir=MOD(jre,2)
        IF (ir.EQ.1) THEN
          write(jpuwast,'(1h1)')
        ELSE
          write(jpuwast,'(a  )')' '
          write(jpuwast,'(a  )')' '
        ENDIF
        WRITE(jpust,wave)
        WRITE(jpuwast,'(1x,a33)') yotitle1(1)
        WRITE(jpuwast,'(1x,a32)') yotitle1(2)
        WRITE(jpuwast,'(2x,a79)') yoperiod(3:81)
        WRITE(jpuwast,'(2x,a79)') yoreg
        WRITE(jpuwast,'(2x,a79)') yosubreg
        CALL STATSP (jpuwast , NTOT,
     .               XMAX, YMAX, XMEAN, YMEAN, STDEVX, STDEVY,
     .               DMAX, BIAS, RMSE, STDEVD, SID, RV, CORR,
     .               SLOPES, RMSESS, RMSEUS, SLOPE1, RMSES1, RMSEU1,
     .               SLOPE2, RINTR2, RMSES2, RMSEU2)
C
C*   5.3.1.3 PERFORM WAVE SCATTER PLOT.
C            --------------------------
C
 5313   CONTINUE
        pval(1)= xtot
        pval(2)= xmean
        pval(3)= ymean
        pval(4)= bias
        pval(5)= stdevd
        pval(6)= sid
        pval(7)= corr
        pval(8)= slopes
        yoaxtite(1) = 'WAM  WAVEHEIGHTS (M)'
        yoaxtite(2) = 'ERS-1  WAVEHEIGHTS (M)'
        CALL pltcontw(jpuso, jpuari, ph(1,1,jre), jpnx,  jpny,
     .                yohead,   yosubhead, yoperiod, yoreg, yosubreg,
     .                ppxa,     ppxa,      ppamaxwa, ppamaxwa,
     .                ppxa,     ppya,      ppdx,     ppdy,
     .                ppxa,     ppamaxwa,  ppatintwa, yoaxtite,
     .                ppya,     ppymaxwa,  ppytintwa, yooption)
C
        yoline(6) =' BIAS (ERS-1 - WAM) '
        CALL textbox(jpkl,yoline,jpkv,pval)
        CALL plmesy (ph(1,1,jre),
     .               ppxa, ppya, ppdx, ppdy, jpnx, jpny, isym)
        swhmean(jre)=xmean
        swherro(jre)=stdevd
        swhbias(jre)=bias
        CALL pnew ('PAGE')
        GOTO 6000
C
C*  5.3.2 WIND SUMS.
C        -----------
C
 5320   CONTINUE
C
C*  5.3.2.1 COMPUTE STATISTICS.
C           -------------------
C
 5321   CONTINUE
        CALL STATSE
     .        (nwindtot,
     .         sumwiX, sumwiY, sumwiXX, sumwiYY, sumwiXY, sumwiXY2,
     .         xwindmax, ywindmax, XMEAN, YMEAN, STDEVX, STDEVY,
     .         dwindmax, BIAS, RMSE, STDEVD, SID, RV, CORR,
     .         SLOPES, RMSESS, RMSEUS, SLOPE1, RMSES1, RMSEU1,
     .         SLOPE2, RINTR2, RMSES2, RMSEU2)
C
C*   5.3.2.2 PRINT WIND STATISTICS.
C            ----------------------
C
 5322   CONTINUE

        ir=MOD(jre,2)
        IF (ir.EQ.1) THEN
          write(jpuwist,'(1h1)')
        ELSE
          write(jpuwist,'(a  )')' '
          write(jpuwist,'(a  )')' '
        ENDIF
        WRITE(jpust,wind)
        WRITE(jpuwist,'(1x,a33)') yotitle2(1)
        WRITE(jpuwist,'(1x,a32)') yotitle2(2)
        WRITE(jpuwist,'(2x,a79)') yoperiod(3:81)
        WRITE(jpuwist,'(2x,a79)') yoreg
        WRITE(jpuwist,'(2x,a79)') yosubreg
        CALL STATSP (jpuwist , nwindtot,
     .               xwindmax, ywindmax, XMEAN, YMEAN, STDEVX, STDEVY,
     .               dwindmax, BIAS, RMSE, STDEVD, SID, RV, CORR,
     .               SLOPES, RMSESS, RMSEUS, SLOPE1, RMSES1, RMSEU1,
     .               SLOPE2, RINTR2, RMSES2, RMSEU2)
C
C*   5.3.2.3 PERFORM WIND SCATTER PLOT.
C            --------------------------
 5323   CONTINUE
        pval(1)= xwindtot
        pval(2)= xmean
        pval(3)= ymean
        pval(4)= bias
        pval(5)= stdevd
        pval(6)= sid
        pval(7)= corr
        pval(8)= slopes
        yoaxtite(1) ='ECMWF  WINDSPEEDS (M/S)'
        yoaxtite(2) ='ERS-1  WINDSPEEDS (M/S)'

        CALL pltcontw(jpuso,jpuari,wiscfi(1,1,jre),jpwindnx,jpwindny,
     .                yohead,   yoextra,  yoperiod, yoreg, yosubreg,
     .                ppaminwi, ppaminwi, ppamaxwi, ppamaxwi,
     .                ppwindxa, ppwindya, ppwinddx, ppwinddy,
     .                ppaminwi, ppamaxwi, ppatintwi, yoaxtite,
     .                ppyminwi, ppymaxwi, ppytintwi, yooption)

        yoline(6) =' BIAS (ERS-1 - ECMWF)'
        CALL textbox(jpkl,yoline,jpkv,pval)
        CALL plmesy (wiscfi(1,1,jre), ppwindxa, ppwindya,
     .               ppwinddx, ppwinddy, jpwindnx, jpwindny, isym)
        u10mean(jre)=xmean
        u10erro(jre)=stdevd
        u10bias(jre)=bias
        CALL pnew ('PAGE')
C
C ----------------------------------------------------------------------
C
C*    6. END OF LOOP OVER AREAS.
C        -----------------------
 6000 CONTINUE
      WRITE(jpuphs ,8100)iyear, imon, iareas,
     .                   (swhmean(j),swherro(j),swhbias(j), j=1,iareas)
C     WRITE(jpupu10,8100)iyear, imon, iareas,
C    .                   (u10mean(j),u10erro(j),u10bias(j), j=1,iareas)
      CLOSE (jpuwast)
C     CLOSE (jpuwist)
      CLOSE (jpust)
      CLOSE (jpuari)
      CLOSE (jpusum)
      CALL pclose
      STOP ' NORMAL'
C
C ----------------------------------------------------------------------
C
C*    7. ABNORMAL TERMINATION.
C        ---------------------
 7000 CONTINUE
      STOP  ' ABNORMAL '
 8000 FORMAT(1x, 2i16, 4f10.2, i4, 4f8.2, f2.0, 2i6 )
 8100 FORMAT(1x, i4,i2.2,i4,60f5.1)
      END
