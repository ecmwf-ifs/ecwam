      SUBROUTINE GRB2WGRD (IU06, ITEST, NKSEK0, NKSEK1, NKSEK2, NKSEK3,
     &                     NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
     &                     KGRIB, ISIZE,
     &                     NGX, NGY, KRGG, KLONRGG, XDELLA, ZDELLO,
     &                     AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &                     CDATE, IFORP, IPARAM, KZLEV, FIELD )
! ----------------------------------------------------------------------    

!!!!!
!     As it stands, it
!!!   should not be use to interpolate a directional field !!!!





!***  *GRB2WGRD* - UNPACKS GRIB DATA FIELD
!                  AND PUTS IT ON THE WAVE MODEL GRID USING BILINEAR
!                  INTERPOLATION IF NECESSARY.

!      J. BIDLOT  ECMWF  FEBRUARY 2000. 

!     PURPOSE.                                                          
!     --------                                                          

!     IT UNPACKS A GRIB DATA FIELD INCLUDING THE INFORMATION ON
!     ITS GRID AND USE IT TO PRODUCE A FIELD ON THE WAVE MODEL GRID
!     USING BILINEAR INTERPOLATION IF THE 2 GRIDS DO NOT MATCH.
!     THE CLOSEST GRID POINT IS USED INSTEAD IF WAVE FIELDS ARE
!     ARE ENCOUNTERED WITH MISSING DATA POINTS.
!     WAVE FIELDS ARE RECOGNISED BY THEIR VALUES OF KSEC1(40)=1045 OR
!     1081 OR 1082 OR 1027.
!     IT ONLY WORKS IF THE INPUT AND OUTPUT GRIDS ARE LATITUDE/LONGITUDE
!     (REGULAR OR IRREGULAR) OR GAUSSIAN GRIDS (FULL OR REDUCED) !!!!

!**   INTERFACE.                                                        
!     ----------                                                        

!      *CALL GRB2WGRD* (IU06, ITEST, NKSEK0, NKSEK1, NKSEK2,NKSEK3,
!    &                  NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
!    &                  KGRIB, ISIZE,
!    &                  NGX, NGY, KRGG, KLONRGG, XDELLA, ZDELLO,
!    &                  AMOWEP, AMOSOP, AMOEAP, AMONOP,
!    &                  CDATE, IFORP, IPARAM, KZLEV, FIELD )

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
!        OUTPUT:
!        *CDATE*  - DATE/TIME OF THE DATA READ.         
!        *IFORP*  - FORCAST PERIOD IN SECONDS.                          
!        *IPARAM* - DATA CODE
!        *KZLEV*  - REFERENCE LEVEL IN full METER
!                   SHOULD BE 0 EXCEPT FOR 233 AND 245 WHERE IT
!                   MIGHT BE DIFFERENT THAN 0 PROVIDED IT WAS
!                   CODED IN KSEC1(7) AND KSEC1(8) SEE *GRIBPAC*
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

      CHARACTER CDATE*12

      LOGICAL LLINTERPOL, LLNEAREST
                                                                        
! ----------------------------------------------------------------------

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
     &             PSEC4, KLENP, KGRIB, KLENG , KWORD, 'I', KRET)
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

      ISTREAM=KSEC2(40)
      IF(ISTREAM.EQ.1045.OR.ISTREAM.EQ.1081.OR.
     &   ISTREAM.EQ.1082.OR.ISTREAM.EQ.1027) THEN
        LLNEAREST=.TRUE.
      ELSE
        LLNEAREST=.FALSE.
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

!         START DATE. 

      IYYYY=(KSEC1(21)-1)*100+KSEC1(10)
      WRITE(CDATE,'(I4.4,4I2.2)') IYYYY,KSEC1(11),KSEC1(12),
     &                             KSEC1(13),KSEC1(14)

!*         DETERMINE TIME UNIT.

      IF (KSEC1(15).EQ.0) THEN                                          
        JCONS = 60                                                     
      ELSEIF (KSEC1(15).EQ.1) THEN                                     
        JCONS = 3600                                                   
      ELSEIF (KSEC1(15).EQ.2) THEN                                     
        JCONS = 86400                                                  
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

!         DETERMINE FORECAST PERIOD.

      IFORP = JCONS*KSEC1(16)

!     GET DATE OF THE FORECAST INSTEAD OF STARTING DATE
      IF(KSEC1(39).EQ.9.OR.KSEC1(39).EQ.10.OR.KSEC1(39).EQ.11) THEN
        CALL INCDATE (CDATE,IFORP)
      ENDIF

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

      RMOEAP = REAL(KSEC2(8)/1000)+0.1*(MOD(KSEC2(8),1000)/100)+
     &         0.01*(MOD(KSEC2(8),100)/10)+0.001*MOD(KSEC2(8),10)


      IF(JRGG.EQ.1) THEN

        ISTART=0
        DO WHILE(KSEC2(23+ISTART).EQ.0)
          ISTART=ISTART+1 
        ENDDO

        ISTOP=0
        DO WHILE(KSEC2(23+NR-1-ISTOP).EQ.0)
          ISTOP=ISTOP+1 
        ENDDO

        IR=0
        DO J=1,NR
          IF(KSEC2(22+J).NE.0) IR=IR+1
        ENDDO
        NR=IR

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

      IF(JRGG.EQ.1) CALL ADJUST (RMOWEP, RMOEAP)

      IPERIODIC = 0
      DELLO=(RMOEAP-RMOWEP)/MAX(1,NC-1)
      IF (RMOEAP-RMOWEP+DELLO.GE.360.) IPERIODIC = 1

      DELLA=(RMONOP-RMOSOP)/MAX(1,NR-1)

      IF(IPERIODIC.EQ.1) THEN
        DELLO=(RMOEAP-RMOWEP)/MAX(1,NC)
        DO J=1,NR
          JSN=NR-J+1
          RDELLO(JSN) = (RMOEAP-RMOWEP+DELLO)/MAX(1,RLONRGG(JSN)) 
        ENDDO
      ELSE
        DELLO=(RMOEAP-RMOWEP)/MAX(1,NC-1)
        DO J=1,NR
          JSN=NR-J+1
          RDELLO(JSN) = (RMOEAP-RMOWEP)/MAX(1,RLONRGG(JSN)-1) 
        ENDDO
      ENDIF


      IF(AMONOP.GT.RMONOP.OR.AMONOP.LT.RMOSOP.OR.
     &   AMOSOP.LT.RMOSOP.OR.AMOSOP.GT.RMONOP.OR.
     &   (JRGG.EQ.0.AND.AMOWEP.LT.RMOWEP).OR.
     &   (JRGG.EQ.0.AND.AMOEAP.GT.RMOEAP+DELLO) ) THEN

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

      IF(AMONOP.EQ.RMONOP.AND.AMOSOP.EQ.RMOSOP.AND.
     &   AMOWEP.EQ.RMOWEP.AND.AMOEAP.EQ.RMOEAP ) THEN
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
          WRITE(IU06,*) ' '
        ENDIF

        IF(IPERIODIC.EQ.1) THEN
          ALLOCATE(WORK(0:NC+1,NR))
        ELSE
          ALLOCATE(WORK(NC,NR))
        ENDIF

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
          KK = MAX(1,INT(XK))
          KSN= NR-KK+1
          KK1=MIN(KK+1,NR)
          KSN1=NR-KK1+1
          DK1=XK-REAL(KK)
          DK2=1.-DK1


          DO I=1,KLONRGG(JSN)

            XI = AMOWEP+(I-1)*ZDELLO(JSN) - RMOWEP
            XI = MOD(XI+720.,360.)

            XII = XI/RDELLO(KSN)+1.
            II = MAX(1-IPERIODIC,INT(XII))
            II = MIN(II,RLONRGG(KSN))
            II1 = MIN(II+1,RLONRGG(KSN)+IPERIODIC)
            DII1=XII-REAL(II)
            DII2=1.-DII1

            XIIP = XI/RDELLO(KSN1)+1.
            IIP = MAX(1-IPERIODIC,INT(XIIP))
            IIP = MIN(IIP,RLONRGG(KSN1)) 
            IIP1 = MIN(IIP+1,RLONRGG(KSN1)+IPERIODIC) 
            DIIP1=XIIP-REAL(IIP)
            DIIP2=1.-DIIP1

            IF(LLNEAREST) THEN
!             DETERMINE WHETHER ANY OF THE 4 CONERS HAS A MISSING DATA
!             IF SO USE THE CLOSEST GRID POINT VALUE.
              IF(WORK(II,KK).EQ.PMISS.OR.
     &           WORK(II1,KK).EQ.PMISS.OR.
     &           WORK(IIP,KK1).EQ.PMISS.OR.
     &           WORK(IIP1,KK1).EQ.PMISS) THEN
                FIELD(I,K)=WORK(II,KK)
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

