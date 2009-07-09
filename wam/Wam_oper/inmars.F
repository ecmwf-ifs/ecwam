      SUBROUTINE INMARS (FILNM, ID1, ID2, PMISS, NKSEK0, NKSEK1, NKSEK2,
     &                   NKSEK3, NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
     &                   CDATE, IFORP, IPARAM, NGX, NGY,
     &                   AMOWEP, AMOSOP, AMOEAP, AMONOP, KRGG, KLONRGG,
     &                   KZLEV, FIELD, KRET, MSG, IU06,
     &                   IREAD, IRANK, NPROC)

! -----------------------------------------------------------------     

!***  *INMARS* - READS IN AND UNPACKS WAVE MODEL FIELDS FROM MARS       
!                AS WELL AS SURFACE TEMPERATURE.

!      H. GUNTHER   ECMWF/GKSS    JANUARY 1990                          
!                                 JUNE    1991  GRIB EDITION 1.         
!      J. BIDLOT    ECMWF         1996    MIGRATION TO FUJITSU
!      J. BIDLOT    ECMWF         1997   ADD READ OF IRREGULAR GRID DATA

!      J. BIDLOT    ECMWF         1998   REMOVE THE INTERCHANGE 
!                                        NORTH AND SOUTH FOR OLD DATA

!     PURPOSE.                                                          
!     --------                                                          

!        READS TARGET FILES FROM THE MARS ARCHIVE AND UNPACK THE DATA.  
!             A) THE DATE OF THE FIELD                                  
!             B) THE DEFINITION OF THE MARS GRID                        
!             C) THE FIELD OF DATA                                      
!             D) THE TYPE OF DATA                                       

!**   INTERFACE.                                                        
!     ----------                                                        

!      *CALL INMARS* (FILNM, ID1, ID2, PMISS, NKSEK0, NKSEK1, NKSEK2,
!    &                NKSEK3, NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
!    &                CDATE, IFORP, IPARAM, NGX, NGY,
!    &                AMOWEP, AMOSOP, AMOEAP, AMONOP, KRGG, KLONRGG,
!    &                KZLEV, FIELD, KRET, MSG, IU06, 
!    &                IREAD, IRANK, NPROC)

!        *FILNM*  - DATA INPUT FILENAME. 
!        *ID1  *  - NUMBER OF COLUMNS IN ARRAY FIELD (1 DIMENSION).     
!        *ID2  *  - NUMBER OF ROWS    IN ARRAY FIELD (2 DIMENSION).     
!        *PMISS*  - REAL VALUE USED TO INDICATE MISSING DATA IN FIELD.
!        *NKSEK0* - SIZE OF KSEC0
!        *NKSEK1* - SIZE OF KSEC1
!        *NKSEK2* - SIZE OF KSEC2
!        *NKSEK3* - SIZE OF KSEC3
!        *NKSEK4* - SIZE OF KSEC4
!        *NPSEK2* - SIZE OF PSEC2
!        *NPSEK3* - SIZE OF PSEC3
!        *KSEK3*  - DEFAULT VALUES FOR KSEC3
!        *PSEK3*  - DEFAULT VALUES FOR PSEC3
!        *CDATE*  - DATE/TIME OF THE DATA READ.         
!        *IFORP*  - FORCAST PERIOD IN SECONDS.                          
!        *IPARAM* - DATA CODE: 229 = WAVE HEIGHT    (METRES).           
!                              230 = WAVE DIRECTION (DEGREE).           
!                              231 = PEAK PERIOD    (SECONDS)           
!                              232 = MEAN PERIOD    (SECONDS).          
!                              139 = SURFACE TEMPERATURE(DEGREE CELCIUS)
!        *NGX  *  - NUMBER OF COLUMNS IN ARRAY FIELD USED.              
!        *NGY  *  - NUMBER OF ROWS    IN ARRAY FIELD USED.              
!        *AMOWEP* - MOST WESTERN LONGITUDE IN GRID (  1, ? ).           
!        *AMOSOP* - MOST SOUTHERN LATITUDE IN GRID.( ? ,NGY).           
!        *AMOEAP* - MOST EASTERN LONGITUDE IN GRID (NGX, ? ).           
!        *AMONOP* - MOST NORTHERN LATITUDE IN GRID ( ? , 1 ).           
!        *KRGG*   - GRID DEFINITION PARAMETER (O=REGULAR, 1=IRREGULAR)
!        *KLONRGG - NUMBER OF GRID POINTS FOR EACH LATITUDE
!        *KZLEV*  - REFERENCE LEVEL IN full METER
!                   SHOULD BE 0 EXCEPT FOR 233, 245 AND 249 WHERE IT
!                   MIGHT BE DIFFERENT THAN 0 PROVIDED IT WAS
!                   CODED IN KSEC1(7) AND KSEC1(8) SEE *GRIBPAC*
!        *FIELD*  - UNPACKED DATA.
!        *KRET*   - ERROR NUMBER (0 = NO ERROR).
!        *MSG*    - ERROR MESSAGE.
!        *IU06*   - OUTPUT UNIT.
!        *IREAD*  - ACCESS TO FILE ONLY FOR PE=IREAD.
!        *IRANK*  - PROCESSOR NUMBER.
!        *NPROC*  - TOTAL NUMBER OF PE'S. 


!       THIS FEATURE HAS BEEN REMOVED :
!       WAVE DATA ARCHIVED BEFORE 3. APRIL 1990 12.00 ARE CORRECTED,
!       THAT MEANS..                                                    
!            - NORTH AND SOUTH ARE INTERCHANGED                         
!            - NUMBER OF LONGITUDES IS SET TO 120                       
!              (THE DIMENSION ID1 HAS TO BE 121)                        

!     EXTERNALS.                                                        
!     ----------                                                        

!     EMOSLIB                                                           

!     *PBOPEN*         OPEN FILE TO BE READ BY PBGRIB
!     *PBGRIB*         READ GRIB FILE
!     *PBCLOSE*        CLOSE FILE
!     *GRIBEX*         UNPACKS MARS DATA.                               
!     *GRPRS1*         PRINTS PRODUCT DEFINITION BLOCK.                 
!     *GRPRS2*         PRINTS GRID  DEFINITION BLOCK.                   

! ----------------------------------------------------------------------

      USE MPL_MODULE

! ----------------------------------------------------------------------

      PARAMETER (NBIT = 430000)

      INTEGER KLONRGG(ID2)                                          
      INTEGER KSEC0(NKSEK0), KSEC1(NKSEK1), KSEC2(NKSEK2),
     &        KSEC3(NKSEK3), KSEC4(NKSEK4)
      INTEGER KSEK3(NKSEK3)
      INTEGER, ALLOCATABLE ::  KGRIB(:)
      INTEGER :: LNAME 
      INTEGER :: IBUF(2) 

      REAL FIELD(ID1,ID2)                                          
      REAL PSEC2(NPSEK2), PSEC3(NPSEK3), PSEC4(ID1*ID2)
      REAL PSEK3(NPSEK3)

      CHARACTER*40 ERRMSG(-5:4), MSG
      CHARACTER HOPER*1, CDATE*14, FILNM*24

      LOGICAL :: LLEXIST
                                                                        
      DATA HOPER / 'D' /                                                
      DATA ERRMSG/' -5                                    *',
     &            ' -4 BIT MAP FOUND                      *',
     &            ' -3                                    *',
     &            ' -2                                    *',
     &            ' -1                                    *',
     &            ' NO ERROR                              *',
     &            ' END OF FILE ENCOUNTED                 *',
     &            ' DECODING ERROR SEE GRIBEX DESCRIPTION *',
     &            ' SUSPICIOUS TIME UNIT IN BLOCK 1       *',
     &            ' DIMENSIONS ID1 OR ID2 ARE TOO SMALL   *'/

! ----------------------------------------------------------------------

      KSEC3 = KSEK3 
      PSEC3 = PSEK3 

!*    1.0 READ MARS DATA IN GRIB CODE FORMAT.                           
!         -----------------------------------                           
     
      LLEXIST=.FALSE.
      LNAME = LEN_TRIM(FILNM)
      INQUIRE(FILE=FILNM(1:LNAME),EXIST=LLEXIST)
      IF(.NOT. LLEXIST) THEN
        WRITE (IU06,*) '*************************************'
        WRITE (IU06,*) '*                                   *'
        WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO INQUIRE  *'
        WRITE (IU06,*) '*  IN INMARS:                       *'
        WRITE (IU06,*) '*  COULD NOT FIND FILE ',FILNM
        WRITE (IU06,*) '*                                   *'
        WRITE (IU06,*) '*************************************'
        MSG=' COULD NOT FIND FILE '//FILNM
        RETURN 
      ENDIF

      IF(IRANK.EQ.IREAD) THEN
        CALL PBOPEN(IUNIT1,FILNM,'r',KRET)                            
        IF(KRET.LT.0) THEN
          IF(KRET.EQ.-1) MSG=' COULD NOT OPEN '//FILNM
          IF(KRET.EQ.-2) MSG='INVALID FILENAME'//FILNM
          RETURN 
        ENDIF

        ISIZE=NBIT
1001    KPLENG=ISIZE*KIND(ISIZE)
        IF(.NOT.ALLOCATED(KGRIB)) ALLOCATE(KGRIB(ISIZE))
        CALL PBGRIB(IUNIT1,KGRIB,KPLENG,KLEN,KRET)
        IF(KRET.EQ.-1) MSG='PBGRIB : END OF FILE ENCOUNTED'
        IF(KRET.EQ.-2) MSG='PBGRIB : FILE HANDLING ERROR'
        IF(KRET.EQ.-3) THEN 
          DEALLOCATE(KGRIB)
          CALL KGRIBSIZE(IU06, IUNIT1, KLEN, ISIZE, 'INMARS')
          GOTO 1001
        ENDIF
        IF(KRET.LT.0) RETURN

        CALL PBCLOSE(IUNIT1,KRET)
        IF(KRET.EQ.-1) MSG='PBCLOSE :FILE HANDLING ERROR'
        IF(KRET.LT.0) RETURN
      ENDIF

      WRITE(IU06,*) ' SUB. INMARS - READ FROM ',FILNM

      CALL MPL_BARRIER(CDSTRING='INMARS: KGRIB')

!     SEND GRIB DATA TO THE OTHER PE'S
!     --------------------------------
      IF(NPROC.GT.1) THEN
        CALL GSTATS(619,0)
        IF(IRANK.EQ.IREAD) THEN
          IBUF(1)=ISIZE
          IBUF(2)=KLEN
        ENDIF
        CALL MPL_BROADCAST(IBUF(1:2),KROOT=IREAD,KTAG=1,
     &                     CDSTRING='INMARS IBUF:')
        IF(IRANK.NE.IREAD) THEN
          ISIZE=IBUF(1)
          KLEN=IBUF(2)
          ALLOCATE(KGRIB(ISIZE))
        ENDIF

        CALL MPL_BROADCAST(KGRIB(1:ISIZE),KROOT=IREAD,KTAG=2,
     &                     CDSTRING='INMARS: KGRIB')

        CALL GSTATS(619,1)
      ENDIF

! ----------------------------------------------------------------------

!*    2.0 UNPACK MARS FIELDS.                                           
!         -------------------                                           

      KLENG = KLEN 
      KLENP = ID1*ID2                                                   
      KRET  = 1                                                         
      CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,     
     &             PSEC4, KLENP, KGRIB, KLENG, KWORD, HOPER, KRET)      
      IF(ALLOCATED(KGRIB)) DEALLOCATE(KGRIB)

! ----------------------------------------------------------------------

!*    3.0 PRODUCT DEFINITION BLOCK.                                     
!         -------------------------                                     

!*    3.1 PRINT BLOCK.                                                  
!         ------------                                                  

!     CALL GRPRS1 (KSEC0, KSEC1)                                        

!*    3.2 DETERMINE DATE FROM KSEC1(10)-KSEC1(16).                      
!         ----------------------------------------                      

!         START DATE.                                                   

      IYYYY=(KSEC1(21)-1)*100+KSEC1(10)
      WRITE(CDATE(1:12),'(I4.4,4I2.2)') IYYYY,KSEC1(11),KSEC1(12),
     &                             KSEC1(13),KSEC1(14)
      CDATE(13:14)='00'

!*         DETERMINE TIME UNIT.                                         

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
        KRET = 3                                                       
        MSG = ERRMSG(KRET)                                             
      ENDIF                                                             

!         DETERMINE FORECAST PERIOD.                                    

      IFORP = JCONS*KSEC1(16)                                           

      IF(KSEC1(39).EQ.9.OR.KSEC1(39).EQ.6.OR.
     &   KSEC1(39).EQ.10.OR.KSEC1(39).EQ.11) THEN
        CALL INCDATE (CDATE,IFORP)
      ENDIF

!*    3.3 DETERMINE CODE FOR DATA FIELD TYPE.                           
!         -----------------------------------                           

      IPARAM = KSEC1(6)                                                 

      IF(KSEC1(7).EQ.105) THEN
        KZLEV=KSEC1(8)
      ELSE
        KZLEV=0
      ENDIF

! ----------------------------------------------------------------------

!*    4.0 GRID DEFINITION BLOCK.                                        
!         ----------------------                                        

!*    4.1 PRINT BLOCK.                                                  
!         ------------                                                  

!     CALL GRPRS2 (KSEC0, KSEC2, PSEC2)                                 

!*    4.2 DETERMINE GRID PARAMETERS.                                    
!         --------------------------                                    

      NGY = KSEC2(3)                                                    

      AMONOP = REAL(KSEC2(4)/1000)+0.1*(MOD(KSEC2(4),1000)/100)+
     &         0.01*(MOD(KSEC2(4),100)/10)+0.001*MOD(KSEC2(4),10)

      AMOWEP = REAL(KSEC2(5)/1000)+0.1*(MOD(KSEC2(5),1000)/100)+
     &         0.01*(MOD(KSEC2(5),100)/10)+0.001*MOD(KSEC2(5),10)

      AMOSOP = REAL(KSEC2(7)/1000)+0.1*(MOD(KSEC2(7),1000)/100)+
     &         0.01*(MOD(KSEC2(7),100)/10)+0.001*MOD(KSEC2(7),10)

      AMOEAP = REAL(KSEC2(8)/1000)+0.1*(MOD(KSEC2(8),1000)/100)+
     &         0.01*(MOD(KSEC2(8),100)/10)+0.001*MOD(KSEC2(8),10)


      KRGG = KSEC2(17)

      IF(KRGG.EQ.1) THEN

        ISTART=0
        DO WHILE(KSEC2(23+ISTART).EQ.0)
          ISTART=ISTART+1 
        ENDDO

        ISTOP=0
        DO WHILE(KSEC2(23+NGY-1-ISTOP).EQ.0)
          ISTOP=ISTOP+1 
        ENDDO

        IC=0
        DO J=1,NGY
          IF(KSEC2(22+J).NE.0) IC=IC+1
        ENDDO
        NGY=IC

        NGX = 0 
        DO J=1,NGY
          JSN=NGY-J+1
          KLONRGG(JSN) = KSEC2(22+J+ISTART) 
          NGX = MAX(NGX,KLONRGG(JSN))
        ENDDO

        IDUM = KSEC2(4)-ISTART*KSEC2(10)
        AMONOP = REAL(IDUM/1000)+0.1*(MOD(IDUM,1000)/100)+
     &           0.01*(MOD(IDUM,100)/10)+0.001*MOD(IDUM,10)

        IDUM = KSEC2(7)+ISTOP*KSEC2(10)
        AMOSOP = REAL(IDUM/1000)+0.1*(MOD(IDUM,1000)/100)+
     &         0.01*(MOD(IDUM,100)/10)+0.001*MOD(IDUM,10)

      ELSEIF(KRGG.EQ.0) THEN
        NGX = KSEC2(2)
        KLONRGG=NGX
      ELSE
        WRITE(IU06,*)
     &  ' SUB INMARS : REPRESENTATION OF THE FIELD NOT KNOWN'
        CALL ABORT1
      ENDIF

      IF ((NGX.GT.ID1).OR.(NGY.GT.ID2)) THEN                            
        KRET = 4                                                       
        MSG = ERRMSG(KRET)                                             
        RETURN                                                         
      ENDIF                                                            

! ----------------------------------------------------------------------

!*    5.0 REARANGE DATA FIELD.                                          
!         --------------------                                          

!*    5.2 NEW DATA.                                                     
!         ---------                                                     

      L = 0                                                          
      DO K=1,NGY                                                
        JSN=NGY-K+1
        DO I=1,KLONRGG(JSN)
          L = L+1                                                     
          FIELD(I,K) = PSEC4(L)                                       
        ENDDO
      ENDDO

      MSG = ERRMSG(KRET)                                                
      IF(KRET.LE.0) KRET=0

      RETURN
      END SUBROUTINE INMARS
