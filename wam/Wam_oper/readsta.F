      SUBROUTINE READSTA (IUIN, CBDT, CEDT, ANALPD, FOREPD, IS,           
     1                    CRSDT, CLSDT, CABDT, CAEDT, IASS, NF, ISTAT)
                                                                        
C -------------------------------------------------------------------   
C                                                                       
C****  *READSTA*  READS THE WAMINFO FILE.                               
C                                                                       
C     PURPOSE.                                                          
C     --------                                                          
C                                                                       
C         *READSTA* - READS THE WAMINFO FILE AND DETERMINES             
C                     BEGIN DATE TIME, END DATE TIME, ANALYSIS PERIOD,  
C                     FORECAST PERIOD, SHIFT,                           
C                     PLOTTING DATES, IASS, NF, ISTATUS                 
C                     FROM THE PREVIOUS RUN AND WRITES THIS             
C                     INFO ON ARRAY CARD                                
C                     AND THE LAST DATE FOR OUTPUT OF RESTART FILES
C                     IF SPECIFIED OTHERWISE IT WILL BE SET TO END DATE
C                                                                       
C     PROGRAMMER.                                                       
C     -----------                                                       
C                                                                       
C          P. JANSSEN                                                   
C          J. BIDLOT NOVEMBER 96  : add CRSDT and CLSDT definition
C          B. HANSEN APRIL    97  : move COMMON CARD INTO INCLUDE FILE.
C                                                                       
C**   INTERFACE.                                                        
C     ----------                                                        
C                                                                       
C          IUIN     INTEGER    DISC UNIT FOR WAMINFO.                   
C          CBDT     CHAR*10    BEGIN DATE-TIME-GROUP.                   
C          CEDT     CHAR*10    END DATE-TIME-GROUP.                     
C          ANALPD   INTEGER    ANALYSIS PERIOD IN SECONDS.
C          FOREPD   INTEGER    FORECAST PERIOD IN SECONDS.
C          IS       INTEGER    WIND INPUT TIME STEP IN SECONDS.           
C          CRSDT    CHAR*10    DATE FOR RESTART FILES OUTPUT
C          CLSDT    CHAR*10    LAST DATE FOR SPECTRUM FILE OUTPUT
C          CABDT    CHAR*10    BEGIN ANALYSIS DATE.                         
C          CAEDT    CHAR*10    END ANALYSIS DATE.                           
C          IASS     INTEGER    ASSIMILATION CONDITION.                  
C          NF       INTEGER    NEW FORECAST CONDITION.                  
C          ISTAT    INTEGER    STATUS FLAGS FOR JOBS.                   
C                                                                       
C     METHOD.                                                           
C     -------                                                           
C                                                                       
C          READ THE WAMINFO FILE                                        
C                                                                       
C     EXTERNALS.                                                        
C     ----------                                                        
C                                                                       
C          ABORT1    TERMINATES PROCESSING                               
C                                                                       
C     REFERENCES.                                                       
C     -----------                                                       
C                                                                       
C          NONE                                                         
C                                                                       
C---------------------------------------------------------------------- 
                                                                        
      INTEGER ANALPD, FOREPD                      
      INTEGER ISTAT(3)                                                  
      CHARACTER*10 CBDT, CEDT, CABDT, CAEDT, CRSDT, CLSDT
C
#include "comcard.h"
C                                                                       
C*    1.0  READ WAMINFO FILE FROM IUIN.                                 
C          ----------------------------                                 
C                                                                       
      CARD=' '
      CRSDT = ' '
      REWIND IUIN                                                       
      DO 1000 I=1,JPCL                                                    
         READ(IUIN,'(A72)', END=3000, ERR=3000) CARD(I)                 
 1000 CONTINUE                                                          
      REWIND IUIN                                                       
C                                                                       
C*    2.0  CONVERT CHARACTER TO INTEGER.                                
C          -----------------------------                                
C                                                                       
      READ (CARD(1), '(15X,A10,4X,A10)') CBDT, CEDT                       
      READ (CARD(4), '(18X,I7)') ANALPD                                 
      READ (CARD(5), '(18X,I7)') FOREPD                                 
      READ (CARD(6), '(28X,I7)') IS                                     
      READ (CARD(7), '(14X,A10,4X,A10)') CABDT, CAEDT                     
      IASS=0                                                            
      IF (CARD(8)(14:14).EQ.'Y') IASS=1                                 
      NF=0                                                              
      IF (CARD(9)(14:14).EQ.'Y') NF=1                                   
      ISTAT(1)=0                                                        
      IF (CARD(10)(16:16).EQ.'F') ISTAT(1)=1                            
      ISTAT(2)=0                                                        
      IF (CARD(11)(16:16).EQ.'F') ISTAT(2)=1                            
      ISTAT(3)=0                                                        
      IF (CARD(12)(16:16).EQ.'F') ISTAT(3)=1                            
                  
      IF(CARD(13).NE.' ') THEN
         READ (CARD(13), '(40X,A10)') CRSDT
      ELSE
         CRSDT = '0000000000'
      ENDIF

      IF(CARD(14).NE.' ') THEN
         READ (CARD(14), '(36X,A10)') CLSDT
      ELSE
         CLSDT = CEDT
      ENDIF

      RETURN                                                            
C                                                                       
C*    3.0  ERROR EXIT                                                   
C          -----------                                                  
C                                                                       
 3000 CONTINUE                                                          
      IU06 = 6                                                          
      WRITE(IU06,*) ' *************************************************'
      WRITE(IU06,*) ' *                                               *'
      WRITE(IU06,*) ' *        FATAL ERROR IN SUB. READSTA            *'
      WRITE(IU06,*) ' *        ===========================            *'
      WRITE(IU06,*) ' * READ ERROR ON UNIT  IUIN = ', IUIN              
      WRITE(IU06,*) ' * WAMINFO FILE.                                 *'
      WRITE(IU06,*) ' * RECORD NUMBER IS I = ', I                       
      WRITE(IU06,*) ' * WAMINFO READ BEFORE ERROR IS:                 *'
      WRITE(IU06,*) ' *                                               *'
      DO 3100 L=1,I                                                     
          WRITE(IU06,'(5X,A72)') CARD(L)                                
 3100 CONTINUE                                                          
      WRITE(IU06,*) ' *                                               *'
      WRITE(IU06,*) ' *  PROGRAMM ABORTS  PROGRAMM ABORTS             *'
      WRITE(IU06,*) ' *                                               *'
      WRITE(IU06,*) ' *************************************************'
      CALL ABORT1                                                        
C                                                                       
      RETURN                                                            
      END                                                               
