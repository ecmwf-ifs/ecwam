! ------------------------------------------------------------------

      PROGRAM BSDCOL 

! ------------------------------------------------------------------

!**** *BSDCOL* - EXTRACTS RESULTS FROM GRIDDED MODEL OUTPUT AND
!****            STORES THESE TOGETHER WITH BUOY DATA AS OBTAINED
!                FROM BUOY/SHIP/DRIFTER DATA. THESE DATA HAVE
!                PREPROCESSED AND ARCHIVED AS SUCH

!                MODEL DATA SHOULD BE ALL FOR THE SAME DATE AND TIME 
!                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                THE OBSERVATIONS WILL BE PICKED FOR THAT DATE.

 
!**   INTERFACE.
!     ----------

!          FILES       UNIT     CONTENT
!          -----       ----     -------

!          MAP          IUMAP   CONTAINS WAVE AND WIND FIELDS.
!          BUO          IUBUO   CONTAINS BUOY DATA.
!          CBSDM        IUCBSDM COLLOCATION FILE.
!          IU06         6       PRINTER OUTPUT.

!          LIBRARIES
!          ---------
!          LIBWAM AND EMOSLIB

!     METHOD.
!     -------

!          THIS PROGRAM TAKES THE MODEL OUTPUTS AS INPUT
!          EXTRACTS TIMESERIES AT SPECIFIED LOCATIONS.

! ----------------------------------------------------------------------

      USE MPL_MODULE,ONLY : MPL_INIT, MPL_END

      PARAMETER (ID1=3601, ID2=1801, IDP=14)
      PARAMETER ( JPLENG=5800000 )

      CHARACTER*5 :: BUOY
      CHARACTER*40 ::  MSG
      CHARACTER*12 :: CDATEM, CDATEB

      LOGICAL :: LLINTER, PER
      LOGICAL :: LLEXIST
      REAL, ALLOCATABLE, DIMENSION(:,:) :: WAVEHT, WAVEDIR, WAVEMPR,
     & HMAX, WINDSPD, WINDDIR, WAVEPPR, UWND, VWND,
     & TAIR, TSEA, PRS, WGUST, FIELD 
      REAL, ALLOCATABLE, DIMENSION(:) :: BLAT, BLON, MLAT, MLON, FLD
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SERIEB, SERIEM, R2DUM 
      CHARACTER*5, ALLOCATABLE, DIMENSION(:) :: NAME, C5DUM
      CHARACTER*12, ALLOCATABLE, DIMENSION(:) :: CDTBO, C12DUM
      CHARACTER*12 :: CDTPR(255)
      LOGICAL :: LFLDPRST(255)
      LOGICAL :: LMODLOC(255)
      LOGICAL :: LNOTFOUND
      LOGICAL :: LLNEW
      REAL :: BSDVAL(IDP)
!     BSDVAL(1) : LONGITUDE (degree)
!     BSDVAL(2) : LATITUDE (degree)
!     BSDVAL(3) : WIND SPEED (m/s)
!     BSDVAL(4) : WIND DIRECTION (degree)
!     BSDVAL(5) : WIND GUST (m/s)
!     BSDVAL(6) : SURFACE PRESSURE (hPa)
!     BSDVAL(7) : AIR TEMPERATURE AT 2M (Celsius)
!     BSDVAL(8) : SEA SURFACE TEMPERATURE (Celsius)
!     BSDVAL(9) : WAVE HEIGHT (m)
!     BSDVAL(10): WAVE MEAN PERIOD (TZ) (s)
!     BSDVAL(11): WAVE PEAK PERIOD (TP) (s)
!     BSDVAL(12): WAVE DIRECTION (degree)
!     BSDVAL(13): ANEMOMETER HEIGHT (m) 
!     BSDVAL(14): MAXIMIM WAVE HEIGHT (m) 

! ----------------------------------------------------------------------

      CALL MPL_INIT(KOUTPUT=1)

      PI=4*ATAN(1.)
      RAD = PI/180.
      DEG = 180./PI

!     1. INITALISATION.
!        --------------
      ALLOCATE(WAVEHT(ID1,ID2), WAVEDIR(ID1,ID2), WAVEMPR(ID1,ID2))
      ALLOCATE(HMAX(ID1,ID2))
      ALLOCATE(WINDSPD(ID1,ID2), WINDDIR(ID1,ID2), WAVEPPR(ID1,ID2))
      ALLOCATE(UWND(ID1,ID2), VWND(ID1,ID2), WGUST(ID1,ID2))
      ALLOCATE(TAIR(ID1,ID2), TSEA(ID1,ID2), PRS(ID1,ID2))
      ALLOCATE(FIELD(ID1,ID2))

      DUMMY = -999.00

      WINDSPD = DUMMY 
      WINDDIR = DUMMY 
      WGUST   = DUMMY
      PRS     = DUMMY 
      TAIR    = DUMMY 
      TSEA    = DUMMY 
      WAVEHT  = DUMMY 
      WAVEMPR = DUMMY 
      WAVEPPR = DUMMY 
      WAVEDIR = DUMMY 
      HMAX    = DUMMY

      UWND    = DUMMY 
      VWND    = DUMMY 

      FIELD   = DUMMY 

!     MAX NUMBER OF STATIONS
      IDB=10000

      ALLOCATE(NAME(IDB))
      ALLOCATE(CDTBO(IDB))
      ALLOCATE(SERIEB(IDB,IDP))

      LFLDPRST=.FALSE.
      LMODLOC=.FALSE.
      CDTPR='000000000000'
!*    1.1 INITALISE UNITS.
!         ----------------

      IU06  = 6

!     INPUT BUOY SHIP DRIFTER DATA
      INQUIRE(FILE='BUO',EXIST=LLEXIST)
      IF(.NOT.LLEXIST) THEN
        WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        WRITE(*,*) 'SOMETHING MUST HAVE GONE WRONG. '
        WRITE(*,*) 'FILE BUO IS NOT PRESENT !!!'
        WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        CALL ABORT
      ENDIF
      IUBUO = I_GET_UNIT(-1, 'BUO', 'S', 'F', 0)
!     test which format (old or new)
      READ (IUBUO,444,ERR=1111,END=1111) BUOY, CDATEB, BSDVAL
      LLNEW=.TRUE.
      GOTO 1112
1111  CONTINUE
      LLNEW=.FALSE.
1112  CONTINUE
      REWIND(IUBUO)


!     OUTPUT COLLOCATION
      IUCBSDM = I_GET_UNIT(-1, 'CBSDM', 'S', 'F', 0)

!*    1.2  GET START DATE AND LENGTH OF ANALYSIS PERIOD FROM ENV.
!          ------------------------------------------------------

      WRITE (IU06,*) ' ----- START OF PROG. BSDCOL -----'


!*    1.4 FIX PARAMETERS.
!         ---------------

      LLINTER = .TRUE.
      NB = 0
      IDATBO = 0
      IDATEL = 0
      NAME= ' '
      IDATBO = 0
      NB = 0
      CDTBO='000000000000'
      SERIEB=DUMMY
 
! ----------------------------------------------------------------------

!*    4. LOOP OVER OUTPUT FILE.
!        ----------------------

!*    4.2 READ IN WIND AND WAVE FIELDS.
!         -----------------------------
      CALL PBOPEN(IUMAP, 'MAP','r',KRET)
      IF(KRET.LT.0) THEN
        WRITE (IU06,*) '******************************************'
        WRITE (IU06,*) '*                                        *'
        WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO PBOPEN       *'
        WRITE (IU06,*) '*   IN BCDCOL                            *'
        IF(KRET.EQ.-1)
     1     WRITE (IU06,*) ' COULD NOT OPEN FILE MAP '
        IF(KRET.EQ.-2) WRITE (IU06,*) ' INVALID FILENAME MAP'
        IF(KRET.EQ.-3) WRITE (IU06,*) ' INVALID OPEN MODE SPECIFIED'
        WRITE (IU06,*) '*                                        *'
        WRITE (IU06,*) '******************************************'
        CALL ABORT
      ENDIF

!     LOOP OVER ALL INPUT FIELDS
      IPARAM=0
 4200 CONTINUE
      CALL INMARSB (6, IUMAP, ID1, ID2, JPLENG, CDATEM, IFORP, IPARAM,
     .              NGX, NGY, AMOWEP, AMOSOP, AMOEAP, AMONOP,
     .              FIELD, IERR, MSG)
      PER = .TRUE.
      CALL ADJUST (AMOWEP, AMOEAP)
      XDELLO = (AMOEAP-AMOWEP)/(NGX-1)
      IF (MOD(AMOEAP-AMOWEP+XDELLO+720., 360.).NE.0.) PER = .FALSE.

!*    4.2.1 END OF FILE ENCOUNTED OR ELSE.
!           ------------------------------

      IF (IERR.EQ.1) THEN
!       ALL DATA ARE ONLINE NOW.
        CALL PBCLOSE(IUMAP,KRET)
        IF(KRET.LT.0) THEN
          WRITE (IU06,*) '************************************'
          WRITE (IU06,*) '*                                  *'
          WRITE (IU06,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
          WRITE (IU06,*) '* IN BCDCOL                        *'
          WRITE (IU06,*) '* FILE MAP                         *'
          WRITE (IU06,*) '*                                  *'
          WRITE (IU06,*) '************************************'
        ENDIF
        GOTO 4300 
      ELSEIF (IERR .GT. 1) THEN
        GOTO 5005
      ELSEIF ( IPARAM .EQ. 218) THEN
        HMAX=FIELD
        LFLDPRST(IPARAM)=.TRUE.
        CDTPR(IPARAM)=CDATEM
      ELSEIF ( IPARAM .EQ. 221) THEN
        WAVEMPR=FIELD
        LFLDPRST(IPARAM)=.TRUE.
        CDTPR(IPARAM)=CDATEM
      ELSEIF ( IPARAM .EQ. 229) THEN
        WAVEHT=FIELD
        LFLDPRST(IPARAM)=.TRUE.
        CDTPR(IPARAM)=CDATEM
      ELSEIF ( IPARAM .EQ. 230) THEN
        WAVEDIR=FIELD
        LFLDPRST(IPARAM)=.TRUE.
        CDTPR(IPARAM)=CDATEM
      ELSEIF ( IPARAM .EQ. 231) THEN
        WAVEPPR=FIELD
        LFLDPRST(IPARAM)=.TRUE.
        CDTPR(IPARAM)=CDATEM
      ELSEIF ( IPARAM .EQ. 245) THEN
        WINDSPD=FIELD
        LFLDPRST(IPARAM)=.TRUE.
        CDTPR(IPARAM)=CDATEM
      ELSEIF ( IPARAM .EQ. 249) THEN
        WINDDIR=FIELD
        LFLDPRST(IPARAM)=.TRUE.
        CDTPR(IPARAM)=CDATEM
      ELSEIF ( IPARAM .EQ. 34) THEN
!       CONVERT TO CELSIUS
        DO J2=1,ID2
          DO J1=1,ID1
            IF(FIELD(J1,J2).NE.DUMMY) THEN
              TSEA(J1,J2)=FIELD(J1,J2)-273.15
            ELSE
              TSEA(J1,J2)=DUMMY
            ENDIF
          ENDDO
        ENDDO
        LFLDPRST(IPARAM)=.TRUE.
        CDTPR(IPARAM)=CDATEM
      ELSEIF ( IPARAM .EQ. 49) THEN
        WGUST=FIELD
        LFLDPRST(IPARAM)=.TRUE.
        CDTPR(IPARAM)=CDATEM
      ELSEIF ( IPARAM .EQ. 151) THEN
!       CONVERT TO hPa
        DO J2=1,ID2
          DO J1=1,ID1
            IF(FIELD(J1,J2).NE.DUMMY) THEN
              PRS(J1,J2)=0.01*FIELD(J1,J2)
            ELSE
              PRS(J1,J2)=DUMMY
            ENDIF
          ENDDO
        ENDDO
        LFLDPRST(IPARAM)=.TRUE.
        CDTPR(IPARAM)=CDATEM
      ELSEIF ( IPARAM .EQ. 165) THEN
        UWND=FIELD
        LFLDPRST(IPARAM)=.TRUE.
        CDTPR(IPARAM)=CDATEM
      ELSEIF ( IPARAM .EQ. 166) THEN
        VWND=FIELD
        LFLDPRST(IPARAM)=.TRUE.
        CDTPR(IPARAM)=CDATEM
      ELSEIF ( IPARAM .EQ. 167) THEN
        DO J2=1,ID2
          DO J1=1,ID1
            IF(FIELD(J1,J2).NE.DUMMY) THEN
              TAIR(J1,J2)=FIELD(J1,J2)-273.15
            ELSE
              TAIR(J1,J2)=DUMMY
            ENDIF
          ENDDO
        ENDDO
        LFLDPRST(IPARAM)=.TRUE.
        CDTPR(IPARAM)=CDATEM
      ENDIF
!     LOOP BACK
      GOTO 4200

4300  CONTINUE

!     COMPUTE WIND SPEED AND DIRECTION FROM U AND V
      IF(LFLDPRST(165) .AND. LFLDPRST(166) ) THEN
        LFLDPRST(245)=.TRUE.
        LFLDPRST(249)=.TRUE.
        CDTPR(245)=CDTPR(165)
        CDTPR(249)=CDTPR(166)
        DO J2=1,ID2
          DO J1=1,ID1
             WINDSPD(J1,J2)=SQRT(UWND(J1,J2)**2+VWND(J1,J2)**2)
             UU = UWND(J1,J2)
             VV = VWND(J1,J2)
             IF (WINDSPD(J1,J2).GT.0.) THEN
               WINDDIR(J1,J2) = ATAN2(UU,VV)
               IF (WINDDIR(J1,J2).LT.0.)
     &             WINDDIR(J1,J2) = WINDDIR(J1,J2)+ZPI
               WINDDIR(J1,J2)=MOD( (WINDDIR(J1,J2)*DEG)+180.,360.)
             ENDIF
          ENDDO
        ENDDO
      ENDIF

!     CHECK MODEL DATE
      CDATEM='000000000000'
      DO IPARAM=1,255
        IF(LFLDPRST(IPARAM)) THEN
          IF(CDATEM.NE.'000000000000') THEN
            IF(CDATEM.NE.CDTPR(IPARAM)) THEN
              WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
              WRITE(*,*) 'SOMETHING MUST HAVE GONE WRONG. '
              WRITE(*,*) 'THE MODEL DATES ARE NOT ALL THE SAME !!!'
              WRITE(*,*) CDATEM,' ',CDTPR(IPARAM)
              WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
              CALL ABORT
            ENDIF
          ENDIF
          CDATEM=CDTPR(IPARAM)
        ENDIF
      ENDDO
      IF(CDATEM.EQ.'000000000000') THEN
        WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        WRITE(*,*) 'SOMETHING MUST HAVE GONE WRONG. '
        WRITE(*,*) 'THE MODEL DATE IS NOT DEFINED !!!'
        WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        CALL ABORT
      ENDIF

!     FIND WHICH FIELD WILL BE USED TO KEEP THE MODEL OBSERVATION POSITION 
      LNOTFOUND=.TRUE.
      DO IPARAM=217,255
        IF(LFLDPRST(IPARAM)) THEN
          LMODLOC(IPARAM)=.TRUE.
          LNOTFOUND=.FALSE.
          EXIT
        ENDIF
      ENDDO
      IF(LNOTFOUND) THEN
        DO IPARAM=1,216
          IF(LFLDPRST(IPARAM)) THEN
            LMODLOC(IPARAM)=.TRUE.
            LNOTFOUND=.FALSE.
            EXIT
          ENDIF
        ENDDO
      ENDIF

      WRITE (IU06,*) '      MODEL DATA FOR ',CDATEM,' WERE READ IN ' 

!*    4.4 LOOP OVER BUOY DATES.
!         ---------------------

 4400 CONTINUE
      IF(LLNEW) THEN
        READ (IUBUO,444,ERR=5001,END=4500) BUOY, CDATEB, BSDVAL
      ELSE
        BSDVAL=DUMMY
        READ (IUBUO,444,ERR=5001,END=4500) BUOY, CDATEB, BSDVAL(1:13)
      ENDIF


!     PICK ONLY DATA THAT CORRESPOND TO MODEL TIME 
!     AND WITH VALID COORDINATES.
      IF (CDATEB.LT.CDATEM .OR. CDATEB.GT.CDATEM ) GOTO 4400
      IF (BSDVAL(1).EQ.-99. .AND. BSDVAL(2).EQ.-99 ) GOTO 4400  

444   FORMAT(A5,1X,A12,1X,F7.2,1X,F6.2,1X,F6.2,1X,F5.0,1X,F6.2,
     &       1X,F8.2,1X,5(F6.2,1X),F5.0,2(1X,F6.2))

 
!*    4.4.1.2 STORE BUOY REPORT.
!             ------------------

 4412 CONTINUE
      NB = NB+1
      if(NB.gt.IDB) then
        IDBOLD=IDB
        IDB=IDB+100
        WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        WRITE(*,*) 'THE SIZE OF THE ARRAYS HAVE BEEN ADJUSTED'
        WRITE(*,*) 'TO ALLOW FOR MORE BUOYS. IF THAT OFTEN OCCURS '
        WRITE(*,*) 'THEN CHANGE THE DEFAULT VALUE IN THE SOURCE !!'
        WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        IF(IDB.gt.200000) then
          WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
          WRITE(*,*) 'SOMETHING MUST HAVE GONE WRONG, THE NUMBER '
          WRITE(*,*) 'BUOYS HAS INCREASED TO 100000. THE PROGRAM '
          WRITE(*,*) 'WILL ABORT AS SAFETY PRECAUTION... '
          WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
          CALL ABORT
        ENDIF
        ALLOCATE(C5DUM(IDBOLD))
        C5DUM=NAME
        DEALLOCATE(NAME)
        ALLOCATE(NAME(IDB))
        NAME=' '
        DO IDUM=1,IDBOLD
          NAME(IDUM)=C5DUM(IDUM)
        ENDDO
        DEALLOCATE(C5DUM)

        ALLOCATE(C12DUM(IDBOLD))
        C12DUM=CDTBO
        DEALLOCATE(CDTBO)
        ALLOCATE(CDTBO(IDB))
        CDTBO='000000000000'
        DO IDUM=1,IDBOLD
          CDTBO(IDUM)=C12DUM(IDUM)
        ENDDO
        DEALLOCATE(C12DUM)

        ALLOCATE(R2DUM(IDBOLD,IDP))
        R2DUM=SERIEB
        DEALLOCATE(SERIEB)
        ALLOCATE(SERIEB(IDB,IDP))
        SERIEB=DUMMY
        DO IP=1,IDP
          DO IDUM=1,IDBOLD
            SERIEB(IDUM,IP)=R2DUM(IDUM,IP)
          ENDDO
        ENDDO
        DEALLOCATE(R2DUM)
      ENDIF

      NAME(NB) = BUOY
      CDTBO(NB) = CDATEB

      DO IP=1,2
        SERIEB(NB,IP) = BSDVAL(IP)
      ENDDO
      DO IP=3,IDP
         IF(BSDVAL(IP).GT.-99.) THEN
          SERIEB(NB,IP) = BSDVAL(IP)
         ELSE
          SERIEB(NB,IP) = DUMMY
         ENDIF
      ENDDO

!     BRANCH BACK FOR NEXT BUOY REPORT.

      GOTO 4400

!*    4.5 END OF DATA FOR THIS TIME WINDOW.
!         ---------------------------------

 4500 CONTINUE
      IF (NB.EQ.0) THEN
         WRITE (IU06,*) '      NO REPORT FOUND FOR CDATEM = ', CDATEM
         GOTO 4800
      ENDIF
      WRITE (IU06,*) NB,' REPORTS FOUND FOR CDATEM = ', CDATEM

 
!*    4.6 EXTRACT MODEL DATA.
!         ------------------

      ALLOCATE(BLAT(NB), BLON(NB), MLAT(NB), MLON(NB), FLD(NB))
      ALLOCATE(SERIEM(NB,IDP))
      SERIEM=DUMMY

      DO IS=1,NB
        BLON(IS)=SERIEB(IS,1)
        BLAT(IS)=SERIEB(IS,2)
      ENDDO
!     GET MODEL VALUES AT BLON, BLAT

!         WIND SPEED.

      IPARAM = 245
      IF(LFLDPRST(IPARAM)) THEN
        CALL EXPOINT (IU06, ID1, ID2, NGX, NGY, IPARAM, NB, NB, LLINTER,
     1                DUMMY, BLON, BLAT, MLON, MLAT, FLD, PER,
     2                AMOWEP, AMOSOP, AMOEAP, AMONOP, WINDSPD)
        DO IS=1,NB
          IF(LMODLOC(IPARAM)) THEN
            SERIEM(IS,1)=MLON(IS)
            SERIEM(IS,2)=MLAT(IS)
          ENDIF
          SERIEM(IS,3)=FLD(IS)
        ENDDO
      ENDIF

!         WIND DIRECTION.

      IPARAM = 249
      IF(LFLDPRST(IPARAM)) THEN
      CALL EXPOINT (IU06, ID1, ID2, NGX, NGY, IPARAM,NB, NB, LLINTER,
     1              DUMMY, BLON, BLAT, MLON, MLAT, FLD, PER,
     2              AMOWEP, AMOSOP, AMOEAP, AMONOP, WINDDIR)
        DO IS=1,NB
          IF(LMODLOC(IPARAM)) THEN
            SERIEM(IS,1)=MLON(IS)
            SERIEM(IS,2)=MLAT(IS)
          ENDIF
          SERIEM(IS,4)=FLD(IS)
        ENDDO
      ENDIF

!         WIND GUST 

      IPARAM = 49 
      IF(LFLDPRST(IPARAM)) THEN
      CALL EXPOINT (IU06, ID1, ID2, NGX, NGY, IPARAM,NB, NB, LLINTER,
     1              DUMMY, BLON, BLAT, MLON, MLAT, FLD, PER,
     2              AMOWEP, AMOSOP, AMOEAP, AMONOP, WGUST)
        DO IS=1,NB
          IF(LMODLOC(IPARAM)) THEN
            SERIEM(IS,1)=MLON(IS)
            SERIEM(IS,2)=MLAT(IS)
          ENDIF
          SERIEM(IS,5)=FLD(IS)
        ENDDO
      ENDIF

!        SURFACE PRESSURE 

      IPARAM = 151 
      IF(LFLDPRST(IPARAM)) THEN
      CALL EXPOINT (IU06, ID1, ID2, NGX, NGY, IPARAM,NB, NB, LLINTER,
     1              DUMMY, BLON, BLAT, MLON, MLAT, FLD, PER,
     2              AMOWEP, AMOSOP, AMOEAP, AMONOP, PRS)
        DO IS=1,NB
          IF(LMODLOC(IPARAM)) THEN
            SERIEM(IS,1)=MLON(IS)
            SERIEM(IS,2)=MLAT(IS)
          ENDIF
          SERIEM(IS,6)=FLD(IS)
        ENDDO
      ENDIF

!        TAIR 

      IPARAM = 167 
      IF(LFLDPRST(IPARAM)) THEN
      CALL EXPOINT (IU06, ID1, ID2, NGX, NGY, IPARAM,NB, NB, LLINTER,
     1              DUMMY, BLON, BLAT, MLON, MLAT, FLD, PER,
     2              AMOWEP, AMOSOP, AMOEAP, AMONOP, TAIR)
        DO IS=1,NB
          IF(LMODLOC(IPARAM)) THEN
            SERIEM(IS,1)=MLON(IS)
            SERIEM(IS,2)=MLAT(IS)
          ENDIF
          SERIEM(IS,7)=FLD(IS)
        ENDDO
      ENDIF

!        TSEA 

      IPARAM = 34 
      IF(LFLDPRST(IPARAM)) THEN
      CALL EXPOINT (IU06, ID1, ID2, NGX, NGY, IPARAM,NB, NB, LLINTER,
     1              DUMMY, BLON, BLAT, MLON, MLAT, FLD, PER,
     2              AMOWEP, AMOSOP, AMOEAP, AMONOP, TSEA)
        DO IS=1,NB
          IF(LMODLOC(IPARAM)) THEN
            SERIEM(IS,1)=MLON(IS)
            SERIEM(IS,2)=MLAT(IS)
          ENDIF
          SERIEM(IS,8)=FLD(IS)
        ENDDO
      ENDIF

!         WAVE HEIGHT.

      IPARAM = 229
      IF(LFLDPRST(IPARAM)) THEN
        CALL EXPOINT (IU06, ID1, ID2, NGX, NGY, IPARAM, NB, NB, LLINTER,
     1                DUMMY, BLON, BLAT, MLON, MLAT, FLD, PER,
     2                AMOWEP, AMOSOP, AMOEAP, AMONOP, WAVEHT)
        DO IS=1,NB
          IF(LMODLOC(IPARAM)) THEN
            SERIEM(IS,1)=MLON(IS)
            SERIEM(IS,2)=MLAT(IS)
          ENDIF
          SERIEM(IS,9)=FLD(IS)
        ENDDO
      ENDIF

!         WAVE MEAN PERIOD.

      IPARAM = 221
      IF(LFLDPRST(IPARAM)) THEN
        CALL EXPOINT (IU06, ID1, ID2, NGX, NGY, IPARAM, NB, NB, LLINTER,
     1                DUMMY, BLON, BLAT, MLON, MLAT, FLD, PER,
     2                AMOWEP, AMOSOP, AMOEAP, AMONOP, WAVEMPR)
        DO IS=1,NB
          IF(LMODLOC(IPARAM)) THEN
            SERIEM(IS,1)=MLON(IS)
            SERIEM(IS,2)=MLAT(IS)
          ENDIF
          SERIEM(IS,10)=FLD(IS)
        ENDDO
      ENDIF

!         WAVE PEAK PERIOD.

      IPARAM = 231
      IF(LFLDPRST(IPARAM)) THEN
        CALL EXPOINT (IU06, ID1, ID2, NGX, NGY, IPARAM, NB, NB, LLINTER,
     1                DUMMY, BLON, BLAT, MLON, MLAT, FLD, PER,
     2                AMOWEP, AMOSOP, AMOEAP, AMONOP, WAVEPPR)
        DO IS=1,NB
          IF(LMODLOC(IPARAM)) THEN
            SERIEM(IS,1)=MLON(IS)
            SERIEM(IS,2)=MLAT(IS)
          ENDIF
          SERIEM(IS,11)=FLD(IS)
        ENDDO
      ENDIF

!         WAVE DIRECTION.

      IPARAM = 230
      IF(LFLDPRST(IPARAM)) THEN
        CALL EXPOINT (IU06, ID1, ID2, NGX, NGY, IPARAM, NB, NB, LLINTER,
     1                DUMMY, BLON, BLAT, MLON, MLAT, FLD, PER,
     2                AMOWEP, AMOSOP, AMOEAP, AMONOP, WAVEDIR)
        DO IS=1,NB
          IF(LMODLOC(IPARAM)) THEN
            SERIEM(IS,1)=MLON(IS)
            SERIEM(IS,2)=MLAT(IS)
          ENDIF
          SERIEM(IS,12)=FLD(IS)
        ENDDO
      ENDIF

!     MODEL WIND SPEED REFERENCE HEIGHT 
      DO IS=1,NB
        SERIEM(IS,13)=10.
      ENDDO

!         MAXIMUM WAVE HEIGHT.

      IPARAM = 218
      IF(LFLDPRST(IPARAM)) THEN
        CALL EXPOINT (IU06, ID1, ID2, NGX, NGY, IPARAM, NB, NB, LLINTER,
     1                DUMMY, BLON, BLAT, MLON, MLAT, FLD, PER,
     2                AMOWEP, AMOSOP, AMOEAP, AMONOP, HMAX)
        DO IS=1,NB
          IF(LMODLOC(IPARAM)) THEN
            SERIEM(IS,1)=MLON(IS)
            SERIEM(IS,2)=MLAT(IS)
          ENDIF
          SERIEM(IS,14)=FLD(IS)
        ENDDO
      ENDIF

!     -180 < longitude < 180      
      DO IS=1,NB
        IF(SERIEM(IS,1).GT.180.) SERIEM(IS,1)=SERIEM(IS,1)-360.
      ENDDO

!*    4.7 WRITE TO COLLOCATION FILE.
!         --------------------------

      DO IB=1,NB
!       IF WAVE HEIGHT FIELD PRESENT, USE IT TO MASK LAND
        IF((LFLDPRST(229) .AND. SERIEM(IB,9).NE.DUMMY) .OR.
     &     .NOT. LFLDPRST(229) ) THEN
            WRITE (IUCBSDM,'(A5,1X,A12,28F8.2)')
     &         NAME(IB), CDTBO(IB),
     &         (SERIEB(IB,IP),SERIEM(IB,IP),IP=1,IDP)
        ENDIF
      ENDDO

!*        SAVE NEW COLLOCATION FILE TO ECFILE.
!         ------------------------------------

 4800 CONTINUE
      CLOSE (UNIT=IUCBSDM)
      

!*    5. TERMINATION.
!        ------------

 5000 CONTINUE
      WRITE (IU06,*) ' ----- NORMAL END OF PROG. BSDCOL -----'
      CALL FLUSH(IU06)
      CALL MPL_END()
      STOP 
 5001 CONTINUE
      WRITE(*,*) ' !!! PROBLEM ACCESSING FILE BUO !!!'
      CALL MPL_END()
      STOP
 5005 CONTINUE
      WRITE(*,*) ' !!! ERROR WHILE DECODING GRIB DATA !!!'
      WRITE(*,*) MSG
      CALL ABORT
      CALL MPL_END()

      END
