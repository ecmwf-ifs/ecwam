      SUBROUTINE INMARS (IUNIT, ID1, ID2, PMISS, CDATE, IFORP, IPARAM,
     1                   NGX, NGY, AMOWEP, AMOSOP, AMOEAP, AMONOP,
     2                   KRGG,KLONRGG,FIELD, KRET, MSG) 
C
C -----------------------------------------------------------------     
C                                                                       
C***  *INMARS* - READS IN AND UNPACKS WAVE MODEL FIELDS FROM MARS       
C                AS WELL AS SURFACE TEMPERATURE.
C                                                                       
C      H. GUNTHER   ECMWF/GKSS    JANUARY 1990                          
C                                 JUNE    1991  GRIB EDITION 1.         
C      J. BIDLOT    ECMWF         1996    MIGRATION TO FUJITSU
C      J. BIDLOT    ECMWF         1997    ADD READ OF IRREGULAR GRID DATA
C                                                                       
C     PURPOSE.                                                          
C     --------                                                          
C                                                                       
C        READS TARGET FILES FROM THE MARS ARCHIVE AND UNPACK THE DATA.  
C             A) THE DATE OF THE FIELD                                  
C             B) THE DEFINITION OF THE MARS GRID                        
C             C) THE FIELD OF DATA                                      
C             D) THE TYPE OF DATA                                       
C                                                                       
C**   INTERFACE.                                                        
C     ----------                                                        
C                                                                       
C      *CALL INMARS*                                                    
C                                                                       
C        *IUNIT*  - DATA INPUT UNIT.                                    
C        *ID1  *  - NUMBER OF COLUMNS IN ARRAY FIELD (1 DIMENSION).     
C        *ID2  *  - NUMBER OF ROWS    IN ARRAY FIELD (2 DIMENSION).     
C        *PMISS*  - REAL VALUE USED TO INDICATE MISSING DATA IN FIELD.
C        *CDATE*  - DATE/TIME OF THE DATA READ (ANALYSIS DATE).         
C        *IFORP*  - FORCAST PERIOD IN SECONDS.                          
C        *IPARAM* - DATA CODE: 229 = WAVE HEIGHT    (METRES).           
C                              230 = WAVE DIRECTION (DEGREE).           
C                              231 = PEAK PERIOD    (SECONDS)           
C                              232 = MEAN PERIOD    (SECONDS).          
C                              139 = SURFACE TEMPERATURE (DEGREE CELCIUS)
C        *NGX  *  - NUMBER OF COLUMNS IN ARRAY FIELD USED.              
C        *NGY  *  - NUMBER OF ROWS    IN ARRAY FIELD USED.              
C        *AMOWEP* - MOST WESTERN LONGITUDE IN GRID (  1, ? ).           
C        *AMOSOP* - MOST SOUTHERN LATITUDE IN GRID.( ? ,NGY).           
C        *AMOEAP* - MOST EASTERN LONGITUDE IN GRID (NGX, ? ).           
C        *KRGG*   - GRID DEFINITION PARAMETER (O=REGULAR, 1=IRREGULAR)
C        *KLONRGG - NUMBER OF GRID POINTS FOR EACH LATITUDE
C        *AMONOP* - MOST NORTHERN LATITUDE IN GRID ( ? , 1 ).           
C        *FIELD*  - UNPACKED DATA.                                      
C        *KRET*   - ERROR NUMBER.                                       
C        *MSG*    _ ERROR MESSAGE.                                      
C                                                                       
C       WAVE DATA ARCHIVED BEFORE 3. APRIL 1990 12.00 ARE CORRECTED,         
C       THAT MEANS..                                                    
C            - NORTH AND SOUTH ARE INTERCHANGED                         
C            - NUMBER OF LONGITUDES IS SET TO 120                       
C              (THE DIMENSION ID1 HAS TO BE 121)                        
C                                                                       
C     EXTERNALS.                                                        
C     ----------                                                        
C                                                                       
C     EMOSLIB                                                           
C                                                                       
C     *PBOPEN*         OPEN FILE TO BE READ BY PBGRIB
C     *PBGRIB*         READ GRIB FILE
C     *PBCLOSE*        CLOSE FILE
C     *GRIBEX*         UNPACKS MARS DATA.                               
C     *GRPRS1*         PRINTS PRODUCT DEFINITION BLOCK.                 
C     *GRPRS2*         PRINTS GRID  DEFINITION BLOCK.                   
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
      PARAMETER (NBIT = 470000)                                          
      DIMENSION FIELD(ID1,ID2)                                          
      DIMENSION KLONRGG(ID2)                                          
      INTEGER   KSEC1(60), KSEC2(22+361), KSEC3(2), KSEC4(42), KSEC0(2)
      INTEGER, ALLOCATABLE ::  KGRIB(:)
      DIMENSION PSEC2(10+ID2), PSEC3(2), PSEC4(ID1*ID2)                 
      CHARACTER*40 ERRMSG(-5:4), MSG                                     
      CHARACTER HOPER*1, CDATE*10, FILNM*24                                                 
                                                                        
      DATA HOPER / 'D' /                                                
      DATA KSEC3 /0,     0 /
      DATA PSEC3 /0., -999./
      DATA ERRMSG/' -5                                    *',
     4            ' -4 BIT MAP FOUND                      *',
     3            ' -3                                    *',
     2            ' -2                                    *',
     1            ' -1                                    *',
     &            ' NO ERROR                              *',
     1            ' END OF FILE ENCOUNTED                 *',
     2            ' DECODING ERROR SEE GRIBEX DESCRIPTION *',
     3            ' SUSPICIOUS TIME UNIT IN BLOCK 1       *',
     4            ' DIMENSIONS ID1 OR ID2 ARE TO SMALL    *'/
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    1.0 READ MARS DATA IN GRIB CODE FORMAT.                           
C         -----------------------------------                           
C                                                                       
C
      IF (IUNIT.LT.10) THEN                                                  
         WRITE(FILNM,'(A,I1)') 'fort.',IUNIT                  
      ELSE
         WRITE(FILNM,'(A,I2)') 'fort.',IUNIT                  
      ENDIF
C             
      CALL PBOPEN(IUNIT1,FILNM,'r',KRET)                            
C
      IF(.NOT.ALLOCATED(KGRIB)) ALLOCATE(KGRIB(NBIT))
      CALL PBGRIB(IUNIT1,KGRIB,NBIT,KLEN,KRET)
      IF(KRET.EQ.-1) GOTO 6000                                           
      IF(KRET.EQ.-2) GOTO 6000                                           
      IF(KRET.EQ.-3) GOTO 6000                                           
C
      CALL PBCLOSE(IUNIT1,KRET)      
C           
 1100 CONTINUE                                                          
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    2.0 UNPACK MARS FIELDS.                                           
C         -------------------                                           
C                                                                       
      KLENG = KLEN 
      KLENP = ID1*ID2                                                   
      KRET  = 1                                                         
      PSEC3(2) = PMISS
      CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,     
     1             PSEC4, KLENP, KGRIB, KLENG, KWORD, HOPER, KRET)      
      IF(ALLOCATED(KGRIB)) DEALLOCATE(KGRIB)
      IF(KRET.GT.0) THEN                                                
         MSG = ERRMSG(2)                                                
         RETURN                                                         
      ENDIF                                                             
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    3.0 PRODUCT DEFINITION BLOCK.                                     
C         -------------------------                                     
C                                                                       
C*    3.1 PRINT BLOCK.                                                  
C         ------------                                                  
C                                                                       
C     CALL GRPRS1 (KSEC0, KSEC1)                                        
C                                                                       
C*    3.2 DETERMINE DATE FROM KSEC1(10)-KSEC1(16).                      
C         ----------------------------------------                      
C                                                                       
C         START DATE.                                                   
C                                                                       
      WRITE(CDATE,'(5I2.2)') KSEC1(10),KSEC1(11),KSEC1(12),KSEC1(13),
     1                       KSEC1(14)                               
C                                                                       
C*         DETERMINE TIME UNIT.                                         
C                                                                       
      IF (KSEC1(15).EQ.0) THEN                                          
         JCONS = 60                                                     
      ELSE IF (KSEC1(15).EQ.1) THEN                                     
         JCONS = 3600                                                   
      ELSE IF (KSEC1(15).EQ.2) THEN                                     
         JCONS = 86400                                                  
      ELSE                                                              
         KRET = 3                                                       
         MSG = ERRMSG(KRET)                                             
      ENDIF                                                             
C                                                                       
C         DETERMINE FORECAST PERIOD.                                    
C                                                                       
      IFORP = JCONS*KSEC1(16)                                           
C                                                                       
C*    3.3 DETERMINE CODE FOR DATA FIELD TYPE.                           
C         -----------------------------------                           
C                                                                       
      IPARAM = KSEC1(6)                                                 
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    4.0 GRID DEFINITION BLOCK.                                        
C         ----------------------                                        
C                                                                       
C*    4.1 PRINT BLOCK.                                                  
C         ------------                                                  
C                                                                       
C     CALL GRPRS2 (KSEC0, KSEC2, PSEC2)                                 
C                                                                       
C*    4.2 DETERMINE GRID PARAMETERS.                                    
C         --------------------------                                    
C                                                                       
      NGY = KSEC2(3)                                                    

      AMONOP = FLOAT(KSEC2(4)/1000)+0.1*(MOD(KSEC2(4),1000)/100)+
     &         0.01*(MOD(KSEC2(4),100)/10)+0.001*MOD(KSEC2(4),10)

      AMOWEP = FLOAT(KSEC2(5)/1000)+0.1*(MOD(KSEC2(5),1000)/100)+
     &         0.01*(MOD(KSEC2(5),100)/10)+0.001*MOD(KSEC2(5),10)

      AMOSOP = FLOAT(KSEC2(7)/1000)+0.1*(MOD(KSEC2(7),1000)/100)+
     &         0.01*(MOD(KSEC2(7),100)/10)+0.001*MOD(KSEC2(7),10)

      AMOEAP = FLOAT(KSEC2(8)/1000)+0.1*(MOD(KSEC2(8),1000)/100)+
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
        AMONOP = FLOAT(IDUM/1000)+0.1*(MOD(IDUM,1000)/100)+
     &           0.01*(MOD(IDUM,100)/10)+0.001*MOD(IDUM,10)

        IDUM = KSEC2(7)+ISTOP*KSEC2(10)
        AMOSOP = FLOAT(IDUM/1000)+0.1*(MOD(IDUM,1000)/100)+
     &         0.01*(MOD(IDUM,100)/10)+0.001*MOD(IDUM,10)

      ELSE IF(KRGG.EQ.0) THEN
        NGX = KSEC2(2)                                                    
        KLONRGG=NGX
      ELSE
        WRITE(*,*) ' SUB INMARS : REPRESENTATION OF THE FIELD NOT KNOWN'
        CALL abort
      ENDIF

      IF ((NGX.GT.ID1).OR.(NGY.GT.ID2)) THEN                            
         KRET = 4                                                       
         MSG = ERRMSG(KRET)                                             
         RETURN                                                         
      END IF                                                            
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    5.0 REARANGE DATA FIELD.                                          
C         --------------------                                          
C                                                                       
      IF (CDATE.LE.'9004031200'.AND.IPARAM.NE.139) THEN
C                                                                       
C*    5.1 INTERCHANGE NORTH AND SOUTH FOR OLD DATA.                     
C         -----------------------------------------                     
C                                                                       
         L = 0                                                          
         DO K=NGY,1,-1                                             
           DO I=1,NGX                                                
             L = L+1                                                     
             FIELD(I,K) = PSEC4(L)                                       
           ENDDO
         ENDDO
         IF (CDATE.LE.'9001071200'.AND.CDATE.GE.'8910041800') THEN          
            CALL CORMARS (FIELD, ID1, ID2, NGX, NGY, IPARAM)            
         ENDIF                                                          
         NGX = 120                                                      
         AMONOP = .001*FLOAT(KSEC2(7))                                  
         AMOSOP = .001*FLOAT(KSEC2(4))                                  
      ELSE                                                              
C                                                                       
C*    5.2 NEW DATA.                                                     
C         ---------                                                     
C                                                                       
         L = 0                                                          
         DO K=1,NGY                                                
           JSN=NGY-K+1
           DO I=1,KLONRGG(JSN)
            L = L+1                                                     
            FIELD(I,K) = PSEC4(L)                                       
           ENDDO
         ENDDO
      ENDIF                                                             
                                                                        
      MSG = ERRMSG(KRET)                                                
                                                                        
      RETURN                                                            
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    6. END OF FILE BUFFER IN.                                         
C        ----------------------                                         
C                                                                       
 6000 CONTINUE                                                          
      KRET = 1                                                          
      MSG = ERRMSG(KRET)                                                
                                                                        
      RETURN                                                            
                                                                        
      END                                                               
