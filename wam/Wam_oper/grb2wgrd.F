      SUBROUTINE GRB2WGRD (IU06, ITEST, NKSEK0, NKSEK1, NKSEK2, NKSEK3,
     &                     NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
     &                     KGRIB, ISIZE,
     &                     NGX, NGY, KRGG, KLONRGG, XDELLA, ZDELLO,
     &                     AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &                     PPREC, PPEPS,
     &                     CDATE, IFORP, IPARAM, KZLEV, KKK, MMM, FIELD)
! ----------------------------------------------------------------------    

!!!!!
!     As it stands, it
!!!   should not be use to interpolate a directional field !!!!
!!!   see condition below which forces the program to abort !!!





!***  *GRB2WGRD* - UNPACKS GRIB DATA FIELD
!                  AND PUTS IT ON THE WAVE MODEL GRID USING BILINEAR
!                  INTERPOLATION IF NECESSARY.

!      J. BIDLOT  ECMWF  FEBRUARY 2000. 

!     PURPOSE.                                                          
!     --------                                                          

!     IT UNPACKS A GRIB DATA FIELD INCLUDING THE INFORMATION ON
!     ITS GRID AND USE IT TO PRODUCE A FIELD ON THE WAVE MODEL GRID
!     USING BILINEAR INTERPOLATION IF THE 2 GRIDS DO NOT MATCH.
!     THE CLOSEST GRID POINT IS USED INSTEAD IF FIELDS ARE
!     ARE ENCOUNTERED WITH MISSING DATA POINTS.
!     WAVE FIELDS ARE RECOGNISED BY THEIR VALUES OF KSEC1(40)
!     SEE WSTREAM_STRG
!     FOR NON WAVE FIELDS IF THE CLOSEST VALUE IS STILL MISSING THEN
!     THE AVERAGE OF THE NEIGHBORING NON MISSING POINTS IS TAKEN.
!     IT ONLY WORKS IF THE INPUT AND OUTPUT GRIDS ARE LATITUDE/LONGITUDE
!     (REGULAR OR IRREGULAR) OR GAUSSIAN GRIDS (FULL OR REDUCED) !!!!

!**   INTERFACE.                                                        
!     ----------                                                        

!      *CALL GRB2WGRD* (IU06, ITEST, NKSEK0, NKSEK1, NKSEK2,NKSEK3,
!    &                  NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
!    &                  KGRIB, ISIZE,
!    &                  NGX, NGY, KRGG, KLONRGG, XDELLA, ZDELLO,
!    &                  AMOWEP, AMOSOP, AMOEAP, AMONOP,
!    &                  PPREC, PPEPS,
!    &                  CDATE, IFORP, IPARAM, KZLEV, KKK, MMM, FIELD)

!        *IU06*   - OUTPUT UNIT.
!        *ITEST*  - TEST MESSAGE LEVEL.
!        *NKSEK0* - SIZE OF KSEC0
!        *NKSEK1* - SIZE OF KSEC1
!        *NKSEK2* - SIZE OF KSEC2
!        *NKSEK3* - SIZE OF KSEC3
!        *NKSEK4* - SIZE OF KSEC4
!        *NPSEK2* - SIZE OF PSEC2
!        *NPSEK3* - SIZE OF PSEC3
!        *KSEK3*  - DEFAULT VALUES FOR KSEC3
!        *PSEK3*  - DEFAULT VALUES FOR PSEC3
!        *KGRIB*  - GRIB CODED DATA ARRAY
!        *ISIZE*  - SIZE OF KGRIB
!        WAVE MODEL GRID SPECIFICATION:
!        *NGX  *  - NUMBER OF COLUMNS IN ARRAY FIELD USED.              
!        *NGY  *  - NUMBER OF ROWS    IN ARRAY FIELD USED.              
!        *KRGG*   - GRID DEFINITION PARAMETER (O=REGULAR, 1=IRREGULAR)
!        *KLONRGG - NUMBER OF GRID POINTS FOR EACH LATITUDE
!        *XDELLA* - GRID POINT SPACING BETWEEN LATITUDES.
!        *ZDELL0* - GRID POINT SPACING PER LATITUDES. 
!        *AMOWEP* - MOST WESTERN LONGITUDE IN GRID (  1, ? ).           
!        *AMOSOP* - MOST SOUTHERN LATITUDE IN GRID.( ? ,NGY).           
!        *AMOEAP* - MOST EASTERN LONGITUDE IN GRID (NGX, ? ).           
!        *AMONOP* - MOST NORTHERN LATITUDE IN GRID ( ? , 1 ).           
!        *PPREC*  - ONLY USED FOR WAVE SPECTRAL FIELD.
!                   SMALL NUMBER USED IN SPECTRAL PACKING OF 251.
!        *PPEPS*  - ONLY USED FOR WAVE SPECTRAL FIELD.
!                   REFERENCE VALUE FOR SPECTRAL PACKING OF 251.
!        OUTPUT:
!        *CDATE*  - DATE/TIME OF THE DATA READ.         
!        *IFORP*  - FORCAST PERIOD IN SECONDS.                          
!        *IPARAM* - DATA CODE
!        *KZLEV*  - REFERENCE LEVEL IN full METER
!                   SHOULD BE 0 EXCEPT FOR 233 AND 245 WHERE IT
!                   MIGHT BE DIFFERENT THAN 0 PROVIDED IT WAS
!                   CODED IN KSEC1(7) AND KSEC1(8) SEE *GRIBPAC*
!        *KKK*    - ONLY USED FOR WAVE SPECTRAL FIELD: DIRECTION NUMBER.
!        *MMM*    - ONLY USED FOR WAVE SPECTRAL FIELD: FREQUENCY NUMBER.
!        *FIELD*  - UNPACKED DATA ON WAVE MODEL GRID.

!     EXTERNALS.                                                        
!     ----------                                                        

!     EMOSLIB                                                           

!     *GRIBEX*         UNPACKS MARS DATA.                               

! ----------------------------------------------------------------------

      INTEGER KLONRGG(NGY)
      INTEGER,ALLOCATABLE :: RLONRGG(:)

      INTEGER KSEC0(NKSEK0), KSEC1(NKSEK1), KSEC2(NKSEK2),
     &        KSEC3(NKSEK3), KSEC4(NKSEK4)
      INTEGER KSEK3(NKSEK3)

      INTEGER KGRIB(ISIZE)

      REAL PSEC2(NPSEK2)
      REAL PSEK3(NPSEK3), PSEC3(NPSEK3)

      REAL FIELD(NGX,NGY)                                          
      REAL ZDELLO(NGY)
      REAL, ALLOCATABLE :: PSEC4(:), RDELLO(:), WORK(:,:)

      CHARACTER CSTREAM*4
      CHARACTER CDATE*12

      LOGICAL LLINTERPOL, LLNEAREST, LLNONWAVE

      DATA KPARFRSTT,KPARFRSTB /0,0/
                                                                        
! ----------------------------------------------------------------------

      CSTREAM='****'
      KSEC3 = KSEK3 
      PSEC3 = PSEK3 

!*    UNPACK MARS FIELDS.                                           
!     -------------------                                           

!     GET SIZE OF THE UNPACKED DATA
      ALLOCATE(PSEC4(1))
      KLENG = ISIZE 
      KLENP = 1
      KRET  = 1
!     CALL GRSDBG(1)
      CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &             PSEC4, KLENP, KGRIB, KLENG , KWORD, 'J', KRET)
      IF(KRET.GT.0) THEN                                                
        WRITE(IU06,*) '***************************************'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*  FATAL ERROR IN SUB. GRB2WGRD       *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*  GRIBEX FAILED WITH ERROR CODE ',KRET
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*     THE PROGRAM ABORTS              *'
        WRITE(IU06,*) '***************************************'
        CALL ABORT1
      ENDIF                                                             

      NR=KSEC2(3)
      NRFULL=NR
      IF(NR.GT.NKSEK2-22) THEN
        WRITE(IU06,*) '***************************************'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*  FATAL ERROR IN SUB. GRB2WGRD       *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*  THE DIMENSION FOR KSEC2 ',NR+22
        WRITE(IU06,*) '*  IS LARGER THAN THE PRESCRIBED ONE ',NKSEK2
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*     THE PROGRAM ABORTS              *'
        WRITE(IU06,*) '***************************************'
        CALL ABORT1
      ENDIF

      JRGG = KSEC2(17)

      IF(JRGG.EQ.1) THEN
        NC=0
        NTOT=0
        DO J=1,NR
          NC = MAX(NC,KSEC2(22+J))
          NTOT=NTOT+KSEC2(22+J)
        ENDDO

        IR=0
        DO J=1,NR
          IF(KSEC2(22+J).NE.0) IR=IR+1
        ENDDO
        NR=IR
        KLEN=NTOT
      ELSEIF(JRGG.EQ.0) THEN
        NC=KSEC2(2)
        KLEN = NC*NR
      ELSE
        WRITE(IU06,*)
     &  ' SUB GRB2WGRD : REPRESENTATION OF THE FIELD NOT KNOWN'
        CALL ABORT1
      ENDIF

      ISTREAM=KSEC1(40)
      CALL WSTREAM_STRG(ISTREAM,CSTREAM)
      IF(CSTREAM.NE.'****') THEN 
        LLNONWAVE=.FALSE.
      ELSE
        LLNONWAVE=.TRUE.
      ENDIF

      IREPR=KSEC2(1)

      IF(IREPR.NE.0.AND.IREPR.NE.4) THEN
        WRITE(IU06,*) '***************************************'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*  FATAL ERROR IN SUB. GRB2WGRD       *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*  UNKNOWN GRID REPRESENTATION = ',IREPR
        WRITE(IU06,*) '*  GRB2WGRD CAN ONLY DEAL WITH        *'
        WRITE(IU06,*) '*  LATITUDE/LONGITUDE GRID (IREPR=0)  *'
        WRITE(IU06,*) '*   OR GAUSSIAN (IREPR=4)              *' 
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*     THE PROGRAM ABORTS              *'
        WRITE(IU06,*) '***************************************'
        CALL ABORT1
      ENDIF

      DEALLOCATE(PSEC4)

!     GET THE DATA

      KSEC3 = KSEK3 
      PSEC3 = PSEK3 
      KLENG = ISIZE 
      KLENP = KLEN
      ALLOCATE(PSEC4(KLENP))
      KRET  = 1                                                         

      CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &             PSEC4, KLENP, KGRIB, KLENG, KWORD, 'D', KRET)
      IF(KRET.GT.0) THEN                                                
        WRITE(IU06,*) '***************************************'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*    FATAL ERROR IN SUB. GRB2WGRD     *'
        WRITE(IU06,*) '*    ============================     *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*  GRIBEX FAILED WITH ERROR CODE ',KRET
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*     THE PROGRAM ABORTS              *'
        WRITE(IU06,*) '***************************************'
        CALL ABORT1
      ENDIF                                                             

! ----------------------------------------------------------------------

!*    DETERMINE INFORMATION ABOUT THE DECODED DATA 
!     --------------------------------------------

!     START DATE. 

      IYYYY=(KSEC1(21)-1)*100+KSEC1(10)
      WRITE(CDATE,'(I4.4,4I2.2)') IYYYY,KSEC1(11),KSEC1(12),
     &                             KSEC1(13),KSEC1(14)

!*    DETERMINE TIME UNIT.

      IF (KSEC1(15).EQ.0) THEN                                          
        JCONS = 60                                                     
      ELSEIF (KSEC1(15).EQ.1) THEN                                     
        JCONS = 3600                                                   
      ELSEIF (KSEC1(15).EQ.2) THEN                                     
        JCONS = 86400                                                  
      ELSEIF (KSEC1(15).EQ.10) THEN                                     
        JCONS = 10800 
      ELSEIF (KSEC1(15).EQ.11) THEN                                     
        JCONS = 21600 
      ELSEIF (KSEC1(15).EQ.12) THEN                                     
        JCONS = 43200 
      ELSE                                                              
        WRITE(IU06,*) '***********************************'
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '*    PROBLEM IN SUB. GRB2WGRD     *'
        WRITE(IU06,*) '*    =========================    *'
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '*  UNABLE TO DECODE THE TIME UNIT *' 
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '***********************************'
        CALL FLUSH(IU06)
      ENDIF                                                             

!     DETERMINE FORECAST PERIOD.

      IFORP = JCONS*KSEC1(16)

!     GET DATE OF THE FORECAST INSTEAD OF STARTING DATE

      CALL INCDATE (CDATE,IFORP)

!*    DETERMINE CODE FOR DATA FIELD TYPE.                           

      IPARAM = KSEC1(6)                                                 

      IF(KSEC1(7).EQ.105) THEN
        KZLEV=KSEC1(8)
      ELSE
        KZLEV=0
      ENDIF

!*    DETERMINE GRID PARAMETERS.                                    
!     --------------------------                                    

      ALLOCATE(RLONRGG(NR))
      RLONRGG=0
      ALLOCATE(RDELLO(NR))

      RMONOP = REAL(KSEC2(4)/1000)+0.1*(MOD(KSEC2(4),1000)/100)+
     &         0.01*(MOD(KSEC2(4),100)/10)+0.001*MOD(KSEC2(4),10)

      RMOWEP = REAL(KSEC2(5)/1000)+0.1*(MOD(KSEC2(5),1000)/100)+
     &         0.01*(MOD(KSEC2(5),100)/10)+0.001*MOD(KSEC2(5),10)

      RMOSOP = REAL(KSEC2(7)/1000)+0.1*(MOD(KSEC2(7),1000)/100)+
     &         0.01*(MOD(KSEC2(7),100)/10)+0.001*MOD(KSEC2(7),10)

!!! THERE IS A DANGER THAT THE DEFINITON FOR RMOEAP MIGHT VARY DUE TO
!!! THE AMBIGOUS DEFINITION FOR IRREGULAR GRIDS. FOR NON WAVE FIELDS,
!!! A GAUSSIAN GRID IMPLIES THAT THE GRID IS GLOBAL, THEREFORE
!!! RMOEAP IS IMPLICITLY KNOWN.
      IF(IREPR.EQ.4 .AND. LLNONWAVE) THEN
        DELLO = 360./MAX(1,NC)
        RMOEAP = RMOWEP+360. - DELLO
        IPERIODIC = 1
      ELSE
        RMOEAP = REAL(KSEC2(8)/1000)+0.1*(MOD(KSEC2(8),1000)/100)+
     &           0.01*(MOD(KSEC2(8),100)/10)+0.001*MOD(KSEC2(8),10)
        IF(JRGG.EQ.1) CALL ADJUST (RMOWEP, RMOEAP)
        IPERIODIC = 0
        DELLO=(RMOEAP-RMOWEP)/MAX(1,NC-1)
        IF (RMOEAP-RMOWEP+1.5*DELLO.GE.360.) IPERIODIC = 1
      ENDIF


      IF(JRGG.EQ.1) THEN

        ISTART=0
        DO WHILE(KSEC2(23+ISTART).EQ.0)
          ISTART=ISTART+1 
        ENDDO

        ISTOP=0
        DO WHILE(KSEC2(23+NRFULL-1-ISTOP).EQ.0)
          ISTOP=ISTOP+1 
        ENDDO

        DO J=1,NR
          JSN=NR-J+1
          RLONRGG(JSN) = KSEC2(22+J+ISTART) 
        ENDDO

        IDUM = KSEC2(4)-ISTART*KSEC2(10)
        RMONOP = REAL(IDUM/1000)+0.1*(MOD(IDUM,1000)/100)+
     &           0.01*(MOD(IDUM,100)/10)+0.001*MOD(IDUM,10)

        IDUM = KSEC2(7)+ISTOP*KSEC2(10)
        RMOSOP = REAL(IDUM/1000)+0.1*(MOD(IDUM,1000)/100)+
     &         0.01*(MOD(IDUM,100)/10)+0.001*MOD(IDUM,10)

      ELSEIF(JRGG.EQ.0) THEN
        RLONRGG=NC
      ELSE
        WRITE(IU06,*)
     &  ' SUB GRB2WGRD : REPRESENTATION OF THE FIELD NOT KNOWN'
        CALL ABORT1
      ENDIF

      PMISS=PSEC3(2)

!     FIND WHETHER INTERPOLATION IS NEEDED

      DELLA=(RMONOP-RMOSOP)/MAX(1,NR-1)

      KAMOWEP=NINT(AMOWEP*100.)
      KAMOEAP=NINT(AMOEAP*100.)
      KAMONOP=NINT(AMONOP*100.)
      KAMOSOP=NINT(AMOSOP*100.)
      KRMOWEP=NINT(RMOWEP*100.)
      KRMOEAP=NINT(RMOEAP*100.)
      KRMONOP=NINT(RMONOP*100.)
      KRMOSOP=NINT(RMOSOP*100.)

      IF(IPERIODIC.EQ.1) THEN
        DELLO=360./MAX(1,NC)
        DO J=1,NR
          JSN=NR-J+1
          RDELLO(JSN) = 360./MAX(1,RLONRGG(JSN)) 
        ENDDO
      ELSE
        DELLO=(RMOEAP-RMOWEP)/MAX(1,NC-1)
        DO J=1,NR
          JSN=NR-J+1
          RDELLO(JSN) = (RMOEAP-RMOWEP)/MAX(1,RLONRGG(JSN)-1) 
        ENDDO
      ENDIF

      ITOP=0
      IF(KAMONOP.GT.KRMONOP) THEN
        IF(KPARFRSTT.NE.IPARAM) THEN
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*  WARNING SUB. GRB2WGRD                  *'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*  THE MODEL NORTHERN BOUNDARY IS OUTSIDE *'
          WRITE(IU06,*) '*  THE INPUT DOMAIN FOR PARAMETER         *'
          WRITE(IU06,*) '  ',IPARAM 
          WRITE(IU06,*) '*  MISSING DATA WILL BE ASSUMED FOR THE   *'
          WRITE(IU06,*) '*  AREA IN BETWEEN.'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*  AMONOP: ', AMONOP, 'RMONOP : ',RMONOP
          KPARFRSTT=IPARAM
        ENDIF
        ITOP=INT((AMONOP-RMONOP)/DELLA)+2
      ENDIF
      IBOT=0
      IF (KAMOSOP.LT.KRMOSOP) THEN
        IF(KPARFRSTB.NE.IPARAM) THEN
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*  WARNING SUB. GRB2WGRD                  *'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*  THE MODEL SOUTHERN BOUNDARY IS OUTSIDE *'
          WRITE(IU06,*) '*  THE INPUT DOMAIN FOR PARAMETER         *'
          WRITE(IU06,*) '  ',IPARAM 
          WRITE(IU06,*) '*  MISSING DATA WILL BE ASSUMED FOR THE   *'
          WRITE(IU06,*) '*  AREA IN BETWEEN.'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*  AMOSOP: ', AMOSOP, 'RMOSOP : ',RMOSOP
          WRITE(IU06,*) '*                                         *'
          KPARFRSTB=IPARAM
        ENDIF
        IBOT=INT((RMOSOP-AMOSOP)/DELLA)+2
      ENDIF

      IF(KAMONOP.LT.KRMOSOP.OR.
     &   KAMOSOP.GT.KRMONOP.OR.
     &   (JRGG.EQ.0.AND.IREPR.NE.4.AND.IPERIODIC.NE.1
     &    .AND.KAMOWEP.LT.KRMOWEP).OR.
     &   (JRGG.EQ.0.AND.IREPR.NE.4.AND.IPERIODIC.NE.1
     &    .AND.KAMOEAP.GT.NINT((RMOEAP+DELLO)*100.)) ) THEN

         WRITE(IU06,*) '***********************************'
         WRITE(IU06,*) '*                                 *'
         WRITE(IU06,*) '*  FATAL ERROR IN SUB. GRB2WGRD   *'
         WRITE(IU06,*) '*  ============================   *'
         WRITE(IU06,*) '*                                 *'
         WRITE(IU06,*) '*  THE MODEL DOMAIN IS OUTSIDE    *' 
         WRITE(IU06,*) '*  THE INPUT DOMAIN FOR PARAMETER *'
         WRITE(IU06,*) '  ',IPARAM 
         WRITE(IU06,*) '*                                 *'
         WRITE(IU06,*) '*  AMOSOP: ', AMOSOP, 'RMOSOP : ',RMOSOP
         WRITE(IU06,*) '*  AMONOP: ', AMONOP, 'RMONOP : ',RMONOP
         WRITE(IU06,*) '*  AMOWEP: ', AMOWEP, 'RMOWEP : ',RMOWEP
         WRITE(IU06,*) '*  AMOEAP: ', AMOEAP, 'RMOEAP : ',RMOEAP
         WRITE(IU06,*) '*  DELLO: ', DELLO
         WRITE(IU06,*) '*                                 *'
         WRITE(IU06,*) '*     THE PROGRAM ABORTS          *'
         WRITE(IU06,*) '***********************************'
         CALL ABORT1
      ENDIF

      LLINTERPOL=.TRUE.

      IF(KAMONOP.EQ.KRMONOP.AND.KAMOSOP.EQ.KRMOSOP.AND.
     &   KAMOWEP.EQ.KRMOWEP.AND.KAMOEAP.EQ.KRMOEAP ) THEN
         IF(JRGG.EQ.KRGG.AND.
     &      NC.EQ.NGX.AND.NR.EQ.NGY) THEN

            LLINTERPOL=.FALSE.

            IF(KRGG.EQ.1) THEN
              DO J=1,NGY
                IF(RLONRGG(J).NE.KLONRGG(J)) THEN
                  LLINTERPOL=.TRUE.
                  EXIT
                ENDIF
              ENDDO
            ENDIF

         ENDIF
      ENDIF

!     TRANSFORM WAVE SPECTRAL VALUE TO THEIR ACTUAL SCALE
      KKK=0
      MMM=0
      IF(IPARAM.EQ.251 .AND. CSTREAM.NE.'****') THEN
        KKK=KSEC1(44)
        MMM=KSEC1(45)
        DO L=1,KLENP
          IF(PSEC4(L).NE.PMISS) THEN
            PSEC4(L) = 10.**(PSEC4(L)-ABS(PPREC))- PPEPS
          ELSE
            PSEC4(L) = PMISS
          ENDIF
        ENDDO
      ENDIF


      IF(.NOT.LLINTERPOL) THEN

!       REARRANGE DATA FIELD.
!       --------------------

        IF(ITEST.GT.0) THEN
          WRITE(IU06,*) ' '
          WRITE(IU06,*) ' THE FIELD WITH GRIB PARAMETER ',IPARAM
          WRITE(IU06,*) ' IS TRANSFERRED ONTO THE WAVE MODEL GRID '
          WRITE(IU06,*) ' '
        ENDIF
        L = 0                                                          
        DO K=1,NGY                                                
          JSN=NGY-K+1
          DO I=1,KLONRGG(JSN)
            L = L+1                                                     
            FIELD(I,K) = PSEC4(L)                                       
          ENDDO
        ENDDO
        DEALLOCATE(PSEC4)

      ELSE

!       INTERPOLATE TO WAVE MODEL GRID
!       ------------------------------

        IF((IPARAM.EQ.230 .OR. IPARAM.EQ.235 .OR. IPARAM.EQ.238 .OR.
     &      IPARAM.EQ.242 .OR. IPARAM.EQ.249) .AND.
     &      CSTREAM.NE.'****') THEN
          WRITE(IU06,*) '****************************************'
          WRITE(IU06,*) '*    FATAL ERROR IN SUB. GRB2WGRD      *'
          WRITE(IU06,*) '*    ============================      *'
          WRITE(IU06,*) '*                                      *'
          WRITE(IU06,*) '* GRB2WGRD NOT YET READY TO INTERPOLATE*' 
          WRITE(IU06,*) '* DIRECTIONAL FIELDS !!                *'
          WRITE(IU06,*) '*                                      *'
          WRITE(IU06,*) '*     THE PROGRAM ABORTS               *'
          WRITE(IU06,*) '****************************************'
          CALL ABORT1
        ENDIF

        IF(.NOT. LLNONWAVE) THEN 
!         FOR WAVE DATA ALWAYS USE THE CLOSEST GRID POINT IF
!         CORNER DATA ARE MISSING.
          LLNEAREST=.TRUE.
        ELSE
!         FOR NON WAVE DATA USE THE CLOSEST GRID POINT IF
!         CORNER DATA ARE FOUND TO BE MISSING.
          NDECODED=ABS(KSEC4(1))
          L=1
          DO WHILE(PSEC4(L).NE.PMISS .AND. L.LE.NDECODED) 
            L=L+1
          ENDDO
          IF(L.GT.NDECODED) THEN
            LLNEAREST=.FALSE.
          ELSE
            LLNEAREST=.TRUE.
          ENDIF
        ENDIF


        IF(ITEST.GT.0) THEN
          WRITE(IU06,*) ' '
          WRITE(IU06,*) ' THE FIELD WITH GRIB PARAMETER ',IPARAM
          IF(IREPR.EQ.0.AND.JRGG.EQ.0) THEN
            WRITE(IU06,*) ' ON A REGULAR LATITUDE/LONGITUDE GRID '
          ELSEIF(IREPR.EQ.0.AND.JRGG.EQ.1) THEN
            WRITE(IU06,*) ' ON A REDUCED LATITUDE/LONGITUDE GRID '
          ELSEIF(IREPR.EQ.4.AND.JRGG.EQ.0) THEN
            WRITE(IU06,*) ' ON A REGULAR GAUSSIAN GRID '
          ELSE
            WRITE(IU06,*) ' ON A REDUCED GAUSSIAN GRID '
          ENDIF
          WRITE(IU06,*) ' IS INTERPOLATED ONTO THE WAVE MODEL GRID '
            IF(LLNEAREST) THEN
              WRITE(IU06,*) ' USING THE CLOSEST GRID POINT '
              WRITE(IU06,*) ' IF MISSING DATA ARE PRESENT '
              WRITE(IU06,*) ' OTHERWISE USING A BILINEAR INTERPOLATION '
            ELSE
              WRITE(IU06,*) ' USING A BILINEAR INTERPOLATION '
            ENDIF
          WRITE(IU06,*) ' '
          WRITE(IU06,*) ' '
        ENDIF


        KKMIN=1-ITOP
        KKMAX=NR+IBOT
        ALLOCATE(WORK(1-IPERIODIC:NC+IPERIODIC,KKMIN:KKMAX))

        DO KK=KKMIN,0
          DO II=1-IPERIODIC,NC+IPERIODIC
            WORK(II,KK) = PMISS
          ENDDO
        ENDDO

        DO KK=NR+1,KKMAX
          DO II=1-IPERIODIC,NC+IPERIODIC
            WORK(II,KK) = PMISS
          ENDDO
        ENDDO

        L = 0                                                          
        DO KK=1,NR                                                
          JSN=NR-KK+1
          DO II=1,RLONRGG(JSN)
            L = L+1                                                     
            WORK(II,KK) = PSEC4(L)                                       
          ENDDO
        ENDDO
        DEALLOCATE(PSEC4)

        IF(IPERIODIC.EQ.1) THEN
          DO KK=1,NR                                                
          JSN=NR-KK+1
          WORK(0,KK)= WORK(RLONRGG(JSN),KK)
          WORK(RLONRGG(JSN)+1,KK)= WORK(1,KK)
          ENDDO
        ENDIF

        DO K=1,NGY                                                

          JSN=NGY-K+1

          XK = RMONOP - (AMOSOP + (JSN-1)*XDELLA)
          XK = XK/DELLA+1.
          KK = MAX(KKMIN,INT(XK))
          KSN= NR-KK+1
          KK1=MIN(KK+1,KKMAX)
          KSN1=NR-KK1+1
          DK1=XK-REAL(KK)
          DK2=1.-DK1


          DO I=1,KLONRGG(JSN)

            XI = AMOWEP+(I-1)*ZDELLO(JSN) - RMOWEP
            XI = MOD(XI+720.,360.)

            KSNLIM=MIN(MAX(KSN,1),NR)
            XII = XI/RDELLO(KSNLIM)+1.
            II = MAX(1-IPERIODIC,INT(XII))
            II = MIN(II,RLONRGG(KSNLIM))
            II1 = MIN(II+1,RLONRGG(KSNLIM)+IPERIODIC)
            DII1=XII-REAL(II)
            DII2=1.-DII1

            KSN1LIM=MIN(MAX(KSN1,1),NR)
            XIIP = XI/RDELLO(KSN1LIM)+1.
            IIP = MAX(1-IPERIODIC,INT(XIIP))
            IIP = MIN(IIP,RLONRGG(KSN1LIM)) 
            IIP1 = MIN(IIP+1,RLONRGG(KSN1LIM)+IPERIODIC) 
            DIIP1=XIIP-REAL(IIP)
            DIIP2=1.-DIIP1

            IF(LLNEAREST) THEN
!             DETERMINE WHETHER ANY OF THE 4 CONERS HAS A MISSING DATA
!             IF SO USE THE CLOSEST GRID POINT VALUE.
              IF(WORK(II,KK).EQ.PMISS.OR.
     &           WORK(II1,KK).EQ.PMISS.OR.
     &           WORK(IIP,KK1).EQ.PMISS.OR.
     &           WORK(IIP1,KK1).EQ.PMISS) THEN
                IF(DK1.LE.0.5) THEN
                  KCL=KK
                  ICL=II1
                  IF(DII1.LE.0.5) ICL=II
                ELSE
                  KCL=KK1
                  ICL=IIP1
                  IF(DIIP1.LE.0.5) ICL=IIP
                ENDIF
                FIELD(I,K)=WORK(ICL,KCL)

!               FOR NON WAVE FIELD, NON MISSING VALUE OVER SEA IS NEEDED
                IF(LLNONWAVE) THEN
                  IF(FIELD(I,K).EQ.PMISS) THEN
                    NCOUNT=0
                    FIELD(I,K)=0.
                    IF(WORK(II,KK).NE.PMISS) THEN
                      FIELD(I,K)=FIELD(I,K)+WORK(II,KK)
                      NCOUNT=NCOUNT+1
                    ENDIF
                    IF(WORK(II1,KK).NE.PMISS) THEN
                      FIELD(I,K)=FIELD(I,K)+WORK(II1,KK)
                      NCOUNT=NCOUNT+1
                    ENDIF
                    IF(WORK(IIP,KK1).NE.PMISS) THEN
                      FIELD(I,K)=FIELD(I,K)+WORK(IIP,KK1)
                      NCOUNT=NCOUNT+1
                    ENDIF
                    IF(WORK(IIP1,KK1).NE.PMISS) THEN
                      FIELD(I,K)=FIELD(I,K)+WORK(IIP1,KK1)
                      NCOUNT=NCOUNT+1
                    ENDIF
                    IF(NCOUNT.GT.0) THEN
                      FIELD(I,K)=FIELD(I,K)/NCOUNT
                    ENDIF
                  ENDIF
                ENDIF

              ELSE
                FIELD(I,K)=DK2*(DII2*WORK(II,KK)+DII1*WORK(II1,KK)) +
     &                  DK1*(DIIP2*WORK(IIP,KK1)+DIIP1*WORK(IIP1,KK1))
              ENDIF
            ELSE
              FIELD(I,K)=DK2*(DII2*WORK(II,KK)+DII1*WORK(II1,KK)) +
     &                DK1*(DIIP2*WORK(IIP,KK1)+DIIP1*WORK(IIP1,KK1))
            ENDIF
            
          ENDDO
        ENDDO

        DEALLOCATE(WORK)
      ENDIF


      DEALLOCATE(RLONRGG)
      DEALLOCATE(RDELLO)

      RETURN
      END SUBROUTINE GRB2WGRD

