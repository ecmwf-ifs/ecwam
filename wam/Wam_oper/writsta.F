      SUBROUTINE WRITSTA (IDIN, CBDT, CEDT, ANALPD, FOREPD, IS,
     &                    CRSDT, CLSDT, CABDT, CAEDT, IASS, NF, ISTAT)
                                                                        
! --------------------------------------------------------------------- 

!***      *WRITSTA* - WRITES A NEW WAMINFO FILE                         

!          P. JANSSEN                                                   
!          B. HANSEN APRIL 97  : add CRSDT and CLSDT definition
!                                move COMMON CARD INTO INCLUDE FILE.

!     PURPOSE.                                                          
!     --------                                                          

!          WRITES A NEW WAMINFO FILE TO UNIT IDIN                       

!     INTERFACE.                                                        
!     ----------                                                        

!         *IDIN*     INTEGER     UNIT NUMBER.                           
!         *CBDT*     CHAR*12     BEGIN DATE TIME GROUP.                 
!         *CEDT*     CHAR*12     END DATE TIME GROUP.                   
!         *ANALPD*   INTEGER     ANALYSIS PERIOD IN SECONDS.
!         *FOREPD*   INTEGER     FORECAST PERIOD IN SECONDS.
!         *IS*       INTEGER     WIND INPUT TIME STEP IN SECONDS.
!         *CRSDT*    CHAR*12     DATE FOR RESTART FILES OUTPUT
!         *CLSDT*    CHAR*12     LAST DATE FOR SPECTRUM FILE OUTPUT
!         *CABDT*    CHAR*12     BEGIN ANALYSIS DATE TIME GROUP.            
!         *CAEDT*    CHAR*12     END ANALYSIS DATE TIME GROUP.              
!         *IASS*     INTEGER     ASSIMILATION CONDITION.                
!         *NF*       INTEGER     NEW FORECAST CONDITION.                
!         *ISTAT*    INTEGER     JOB STATUS INDICATORS.                 

!     METHOD.                                                           
!     -------                                                           

!          NEW WAMINFO FILE IS WRITTEN TO UNIT IDNIN                    

!     EXTERNALS.                                                        
!     ----------                                                        

!          NONE                                                         

!     REFERENCES.                                                       
!     -----------                                                       

!          NONE                                                         

!---------------------------------------------------------------------- 

      USE YOWCARD  , ONLY : JPCL     ,CARD

      INTEGER ANALPD, FOREPD                      
      INTEGER ISTAT(3)                                                  
      CHARACTER*12 CBDT, CEDT, CABDT, CAEDT, CRSDT, CLSDT

! --------------------------------------------------------------------- 

!*    1.0  FILL CARD ARRAY.                                             
!          ----------------                                             

      WRITE(CARD(1), 1) CBDT, CEDT
    1 FORMAT('RUN MODEL FROM ',A12,' TO ',A12)
                                                                        
      CARD(2)='PARAM=10U/10V,'

      CARD(3)='REPRES=GG,'

      WRITE(CARD(4), 2) ANALPD                                          
    2 FORMAT('ANALYSIS PERIOD = ',I7)        
                                                                        
      WRITE(CARD(5), 3) FOREPD                                          
    3 FORMAT('FORECAST PERIOD = ',I10)         
                                                                        
      WRITE(CARD(6), 4) IS                                              
    4 FORMAT('WIND TIME STEP IN SECONDS = ',I7)
                                                                        
      WRITE(CARD(7), 5) CABDT, CAEDT
    5 FORMAT('ANALYSIS FROM ',A12,' TO ',A12)
                                                                        
      IF (IASS.EQ.0) THEN                                               
        CARD(8) = 'ASSIMILATION NO'        
      ELSE                                                              
        CARD(8) = 'ASSIMILATION YES'      
      ENDIF                                                            
                                                                        
      IF (NF.EQ.0) THEN                                                 
        CARD(9) = 'NEW FORECAST NO'      
      ELSE                                                              
        CARD(9) = 'NEW FORECAST YES'    
      ENDIF                                                            
                                                                        
      DO I=1,3                                                     
        IF (ISTAT(I).EQ.0) THEN                                        
          WRITE(CARD(9+I), 6) I                                       
    6     FORMAT('STATUS STORM',I1,'= UNFINISHED')                    
        ELSE                                                           
          WRITE(CARD(9+I), 7) I                                       
    7     FORMAT('STATUS STORM',I1,'= FINISHED')                      
        ENDIF                                                         
      ENDDO
      WRITE(CARD(13),8) CRSDT
    8 FORMAT('DATE FOR OUTPUT OF BOTH RESTART FILES = ',A12)

      WRITE(CARD(14),9) CLSDT
    9 FORMAT('LAST DATE FOR SPECTRA FILE OUTPUT = ',A12)

!*    2.0  WRITE CARD TO IDIN.                                          
!          -------------------                                          

      DO I=1,JPCL
        WRITE(IDIN,'(A72)') CARD(I)                                    
      ENDDO

      RETURN                                                            
      END SUBROUTINE WRITSTA
