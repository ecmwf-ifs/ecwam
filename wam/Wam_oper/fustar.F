      SUBROUTINE FUSTAR (USMO, USA, USOI,                               
     1                  EWFG, EWOI, EWA, T, FMWA, IJS, IJL)             
                                                                        
C ----------------------------------------------------------------------
C                                                                       
C**** *FUSTAR* - ESTIMATE THE ANALYSED FRICTION VELOCITY FROM THE       
C                DURATION AND THE "MEASURED" WINDSEA ENERGY.            
C                ESTIMATE THE ANALYSED WINDSEA MEAN FREQUENCY           
C                FROM THE "MEASURED" WINDSEA ENERGY                     
C                                                                       
C     P.LIONELLO     ECMWF       FEBRUARY 1989                          
C         (FROM MODIFICATION OF A CODE BY P.JANSSEN                     
C           AND P.LIONELLO - SUMMER '87 - )                             
C                                                                       
C     PURPOSE.                                                          
C     --------                                                          
C                                                                       
C        ESTIMATE THE ANALYSED FRICTION VELOCITY FROM THE URATION AND   
C        THE "MEASURED" WINDSEA ENERGY. ESTIMATE THE ANALYSED WINDSEA   
C        MEAN FREQUENCY FROM THE "MEASURED" WINDSEA ENERGY.             
C                                                                       
C**   INTERFACE.                                                        
C     ----------                                                        
C                                                                       
C        *CALL* *FUSTAR (USMO, USA, USOI,                               
C                        EWFG, EWOI, EWA, T, FMWA, IJS, IJL)*           
C                                                                       
C        *USMO*   REAL     FIRST GUESS FRICTION VELOCITY.               
C        *USA*    REAL     ANALISED FRICTION VELOCITY.                  
C        *USOI*   REAL     MEASURED FRICTION VELOCITY.                  
C        *EWFG*   REAL     WINDSEA ENERGY (FIRST GUESS).                
C        *EWOI*   REAL     ESTIMATED WIND-SEA ENERGY FROM MEASUREMENTS. 
C        *EWA*    REAL     ANALISED WIND-SEA ENERGY.                    
C        *T*      REAL     DURATION.                                    
C        *FMWA*   REAL     ANALYSED WINDSEA MEAN FREQUENCY.             
C        *IJS*    INTEGER  FIRST INDEX IN BLOCK.                        
C        *IJL*    INTEGER  LAST  INDEX IN BLOCK.                        
C                                                                       
C      METHOD.                                                          
C      -------                                                          
C                                                                       
C       THE DURATION OF THE WIND SEA IS USED TO COMPUTE THE FRICTION    
C       VELOCITY BY NEWTONS METHOD FROM THE MODEL DURATION CURVE AND    
C       THE MEASURED ENERGY.  (FIVE ITERATIONS ARE DONE)                
C       THIS VALUE IS USED TO UPDATE THE MODEL IF IT IS                 
C       IN REASONABLE AGREEMENT WITH THE MEASUREMENT FRICTION VELOCITY. 
C       IF THERE IS NO AGREEMENT THE MEASURED WIND SPEED IS USED        
C       TO COMPUTE THE WIND SEA.                                        
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
#include "param.h"
C
#include "commpp.h"
C
#include "comtest.h"
C
C ----------------------------------------------------------------------

      REAL,DIMENSION(NINF:NSUP) :: USMO, USA, USOI, EWOI,      
     1                             EWA, T, FMWA, EWFG          
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
#include "parcons.h"
C                                                                       
      PARAMETER (EO = 955., A = 6.02*10.**(-5), B = 0.695,              
     1           AF = 0.000168, BF = -3.27)                             
C                                                                       
C       *EO,A,B* : PARAMETERS OF THE ENERGY GROWTH CURVE                
C                  ESTAR = EO*TANH(A*TSTAR**B).                         
C       *AF,BF*  : PARAMETERS IN  MEAN FREQUENCY - ENERGY RELATION      
C                  YNU(X) = (ESTAR/AF)**(1/BF).                         
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C     INLINE FUNCTION.                                                  
C     ----------------                                                  
C                                                                       
C     DIMENSIONLESS FREQUENCY AS FUNCTION OF DIMENSIONLESS ENERGY       
C                                                                       
      YNU(X) = (X/AF)**(1./BF)                                          
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    1. INITIALIZE OUTPUT FIELDS.                                      
C        -------------------------                                      
C                                                                       
 1000 CONTINUE                                                          
      NNW=0                                                             
      DO 1001 IJ=IJS,IJL                                                
         USA(IJ)  =   -99.                                              
         EWA(IJ)  =  -999.                                              
         FMWA(IJ) = -9999.                                              
 1001 CONTINUE                                                          
C                                                                       
C*    2. LOOP OVER GRID POINTS.                                         
C        ----------------------                                         
C                                                                       
 2000 CONTINUE                                                          
      DO 2001 IJ=IJS,IJL                                                
         IF (EWFG(IJ).GT.0.01) THEN                                     
C                                                                       
C*    2.1  DETERMINE FRICTION VELOCITY FROM DURATION GROWTH CURVE       
C          IF WINDSEA IS PRESENT.                                       
C           ------------------------------------------------------      
C                                                                       
            UMODEL = USMO(IJ)                                           
            USTANAL = UMODEL                                            
            EWSR = EWOI(IJ)                                             
            TX = T(IJ)                                                  
            DO 3000 ITER=1,5                                            
               TSTAR = G*TX/USTANAL                                     
               XX = A*TSTAR**B                                          
               IF (XX.GE.2.) THEN                                       
                  USTANAL = (EWSR*G**2/EO)**(0.25)                      
               ELSE                                                     
                  TT = TANH(XX)                                         
                  RKSI = EWSR*G**2/USTANAL**4/EO                        
                  USTANAL = USTANAL*(1.-(RKSI-TT)                       
     1                          /(B*XX*(1.-TT**2)-4.*TT))               
               END IF                                                   
 3000       CONTINUE                                                    
C                                                                       
C*    2.2  IF SIGNS OF (USA-USMO) AND (USOI-USMO) DIFFER, THEN          
C          THE MODEL IS WRONG ABOUT THE CORRECT RATIO BETWEEN           
C          WINDSEA AND SWELL. THEREFORE NO ASSIMILATION IS DONE.        
C          -----------------------------------------------------        
C                                                                       
            IF (USOI(IJ).GT.0.) THEN                                    
               SIGNR  = SIGN(1.,(USTANAL-UMODEL))                       
               SIGNM2 = SIGN(1.,(0.8*USOI(IJ)-UMODEL))                  
               SIGNM1 = SIGN(1.,(1.2*USOI(IJ)-UMODEL))                  
               IF (SIGNR.NE.SIGNM1.AND.SIGNR.NE.SIGNM2) THEN            
                  NNW = NNW+1                                           
               ELSE                                                     
                  USA(IJ) = USTANAL                                     
                  EWA(IJ) = EWSR                                        
                  FMWA(IJ) = YNU(EWSR*G**2/USTANAL**4)*G/USTANAL        
               END IF                                                   
            ELSE                                                        
               USA(IJ) = USTANAL                                        
               EWA(IJ) = EWSR                                           
               FMWA(IJ) = YNU(EWSR*G**2/USTANAL**4)*G/USTANAL           
            ENDIF                                                       
         ENDIF                                                          
 2001 CONTINUE                                                          
C                                                                       
      WRITE(IU06,*) ' WINDSEA / USTAR ASSIMILATION FAILED AT ',            
     1             NNW, ' POINTS'                                       
C                                                                       
      RETURN                                                            
      END                                                               
