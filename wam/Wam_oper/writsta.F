      SUBROUTINE WRITSTA (IDIN, CBDT, CEDT, ANALPD, FOREPD, IS,
     1                    CRSDT, CLSDT, CABDT, CAEDT, IASS, NF, ISTAT)
                                                                        
C --------------------------------------------------------------------- 
C                                                                       
C***      *WRITSTA* - WRITES A NEW WAMINFO FILE                         
C                                                                       
C          P. JANSSEN                                                   
C          B. HANSEN APRIL 97  : add CRSDT and CLSDT definition
C                                move COMMON CARD INTO INCLUDE FILE.
C                                                                       
C     PURPOSE.                                                          
C     --------                                                          
C                                                                       
C          WRITES A NEW WAMINFO FILE TO UNIT IDIN                       
C                                                                       
C     INTERFACE.                                                        
C     ----------                                                        
C                                                                       
C         *IDIN*     INTEGER     UNIT NUMBER.                           
C         *CBDT*     CHAR*12     BEGIN DATE TIME GROUP.                 
C         *CEDT*     CHAR*12     END DATE TIME GROUP.                   
C         *ANALPD*   INTEGER     ANALYSIS PERIOD IN SECONDS.
C         *FOREPD*   INTEGER     FORECAST PERIOD IN SECONDS.
C         *IS*       INTEGER     WIND INPUT TIME STEP IN SECONDS.
C         *CRSDT*    CHAR*12     DATE FOR RESTART FILES OUTPUT
C         *CLSDT*    CHAR*12     LAST DATE FOR SPECTRUM FILE OUTPUT
C         *CABDT*    CHAR*12     BEGIN ANALYSIS DATE TIME GROUP.            
C         *CAEDT*    CHAR*12     END ANALYSIS DATE TIME GROUP.              
C         *IASS*     INTEGER     ASSIMILATION CONDITION.                
C         *NF*       INTEGER     NEW FORECAST CONDITION.                
C         *ISTAT*    INTEGER     JOB STATUS INDICATORS.                 
C                                                                       
C     METHOD.                                                           
C     -------                                                           
C                                                                       
C          NEW WAMINFO FILE IS WRITTEN TO UNIT IDNIN                    
C                                                                       
C     EXTERNALS.                                                        
C     ----------                                                        
C                                                                       
C          NONE                                                         
C                                                                       
C     REFERENCES.                                                       
C     -----------                                                       
C                                                                       
C          NONE                                                         
C                                                                       
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
      INTEGER ANALPD, FOREPD                      
      INTEGER ISTAT(3)                                                  
      CHARACTER*12 CBDT, CEDT, CABDT, CAEDT, CRSDT, CLSDT
C
#include "comcard.h"
C                                                                       
C --------------------------------------------------------------------- 
C                                                                       
C*    1.0  FILL CARD ARRAY.                                             
C          ----------------                                             
C                                                                       
      WRITE(CARD(1), 1) CBDT, CEDT
    1 FORMAT('RUN MODEL FROM ',A12,' TO ',A12)                          
                                                                        
      WRITE(CARD(4), 2) ANALPD                                          
    2 FORMAT('ANALYSIS PERIOD = ',I7)                                   
                                                                        
      WRITE(CARD(5), 3) FOREPD                                          
    3 FORMAT('FORECAST PERIOD = ',I7)                                   
                                                                        
      WRITE(CARD(6), 4) IS                                              
    4 FORMAT('WIND TIME STEP IN SECONDS = ',I7)                           
                                                                        
      WRITE(CARD(7), 5) CABDT, CAEDT
    5 FORMAT('ANALYSIS FROM ',A12,' TO ',A12)                               
                                                                        
      IF (IASS.EQ.0) THEN                                               
         CARD(8) = 'ASSIMILATION NO'                                    
      ELSE                                                              
         CARD(8) = 'ASSIMILATION YES'                                   
      END IF                                                            
                                                                        
      IF (NF.EQ.0) THEN                                                 
         CARD(9) = 'NEW FORECAST NO'                                    
      ELSE                                                              
         CARD(9) = 'NEW FORECAST YES'                                   
      END IF                                                            
                                                                        
      DO 1000 I=1,3                                                     
         IF (ISTAT(I).EQ.0) THEN                                        
            WRITE(CARD(9+I), 6) I                                       
    6       FORMAT('STATUS STORM',I1,'= UNFINISHED')                    
         ELSE                                                           
            WRITE(CARD(9+I), 7) I                                       
    7       FORMAT('STATUS STORM',I1,'= FINISHED')                      
         END IF                                                         
 1000 CONTINUE                                                          
      CARD(13)(41:52) = CRSDT
      CARD(14)(37:48) = CLSDT
C                                                                       
C*    2.0  WRITE CARD TO IDIN.                                          
C          -------------------                                          
C                                                                       
      DO 2000 I=1,JPCL
         WRITE(IDIN,'(A72)') CARD(I)                                    
 2000 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
      BLOCK DATA WAMINFOT
C
#include "comcard.h"
C
        DATA CARD /
     1'RUN MODEL FROM YYYYMMDDHHMM TO YYYYMMDDHHMM                    ',
     2'PARAM=10U/10V,                                                 ',
     3'REPRES=GG,                                                     ',
     4'ANALYSIS PERIOD = hhhhhhh                                      ',
     5'FORECAST PERIOD = hhhhhhh                                      ',
     6'WIND TIME STEP IN SECONDS = hhhhhhh                            ',
     7'ANALYSIS FROM YYYYMMDDHHMM TO YYYYMMDDHHMM                     ',
     8'ASSIMILATION [y|n]                                             ',
     9'NEW FORECAST [y|n]                                             ',
     A'STATUS STORM1= [unfinished|finished]                           ',
     1'STATUS STORM2= [unfinished|finished]                           ',
     2'STATUS STORM3= [unfinished|finished]                           ',
     3'DATE FOR OUTPUT OF BOTH RESTART FILES = YYYYMMDDHHMM           ',
     4'LAST DATE FOR SPECTRA FILE OUTPUT = YYYYMMDDHHMM               '/
      END BLOCK DATA WAMINFOT

