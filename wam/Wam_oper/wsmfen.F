      SUBROUTINE WSMFEN (FSEA, EW, FM, USTT)                            
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
C**** *WSMFEN* - COMPUTES MEAN FREQUENCY AND ENERGY OF THE WINDSEA      
C****            SPECTRUM.                                              
C                                                                       
C     P.LIONELLO     ECMWF       APRIL 1990                             
C                                                                       
C     PURPOSE.                                                          
C     --------                                                          
C                                                                       
C       COMPUTES MEAN FREQUENCY AND ENERGY OF THE WINDSEA SPECTRUM.     
C                                                                       
C**   INTERFACE.                                                        
C     ----------                                                        
C                                                                       
C       *CALL* *WSMFEN (FSEA, EW, FM, USTT)*                            
C                                                                       
C        *FSEA*   REAL     WINDSEA PART OF THE WAVE SPECTRUM.           
C        *EW*     REAL     WINDSEA ENERGY .                             
C        *FM*     REAL     MEAN FREQUENCY OF THE WINDSEA.               
C        *USTT*   REAL     FRICTION VELOCITY.                           
C                                                                       
C     METHOD.                                                           
C     -------                                                           
C                                                                       
C       ENERGY AND MEAN FREQUENCY ARE FIRST COMPUTED BY INTEGRATION     
C       OF THE WINDSEA PART OF THE SPECTRUM:                            
C       THIS IMPLIES AN UNDERESTIMATION BOTH OF ENERGY AND MEAN         
C       FREQUENCY. THE UNDERESTIMATED VALUES ARE USED IN THE GROWTH     
C       CURVE PRODUCING OVERSTIMATES OF BOTH ENERGY AND MEAN FREQUENCY. 
C       THE AVERAGE OF THE TWO ESTIMATES IS TAKEN TO PROVIDE A BEST     
C       ESTIMATE.                                                       
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
#include "param.h"
C                                                                       
#include "comfred.h"
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
      DIMENSION FSEA(NANG,NFRE)                                         
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
#include "parcons.h"
c                                                                       
      PARAMETER (AF = 0.000168, BF = -3.27)                             
C                                                                       
C       *AF*   REAL   PARAMETER OF THE ENERGY - MEAN FREQUENCY RELATION.
C       *BF*   REAL   PARAMETER OF THE ENERGY - MEAN FREQUENCY RELATION.
C       *G*    REAL   ACCELLERATION OF GRAVITY.                         
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C     INLINE FUNCTIONS.                                                 
C     -----------------                                                 
C                                                                       
C     ENERGY AS FUNCTION OF MEAN FREQUENCY.                             
C                                                                       
      EN(X)=AF*X**BF                                                    
C                                                                       
C     MEAN FREQUENCY AS FUNCTION OF ENERGY.                             
C                                                                       
      YNU(X)=(X/AF)**(1./BF)                                            
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    1. INTEGRATING THE WINDSEA PART OF THE SPECTRUM.                  
C*       (UNDERESTIMATE IS IMPLIED)                                     
C        ---------------------------------------------                  
C                                                                       
C*    1.1 COMPUTATION OF THE WIND-SEA ENERGY.                           
C         -----------------------------------                           
C                                                                       
      EW = EPSMIN                                                      
      DO 1101 M=1,NFRE                                                  
         TEMP = 0.                                                      
         DO 1102 K=1,NANG                                               
            TEMP = TEMP+FSEA(K,M)                                       
 1102    CONTINUE                                                       
         EW = EW+TEMP*DFIM(M)                                           
 1101 CONTINUE                                                          
      DELT25 = FR(NFRE)/3.*DELTH                                        
      EW = EW+DELT25*TEMP                                               
C                                                                       
C*    1.2 COMPUTATION OF THE MEAN FREQUENCY.                            
C         ----------------------------------                            
C                                                                       
      FM = EPSMIN                                                     
      DO 1201 M=1,NFRE                                                  
         SPFB = 0.                                                      
         DO 1202 K=1,NANG                                               
            SPFB = SPFB+FSEA(K,M)                                       
 1202    CONTINUE                                                       
         SPFB = SPFB*DELTH                                              
         FM = FM+SPFB*FR(M)**2.*(1.1-1./1.1)/2.                         
 1201 CONTINUE                                                          
      FM = FM+SPFB*FR(NFRE)**2./3.                                      
      FM = FM/EW                                                        
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    2. ESTIMATES FROM MODEL RELATIONS.                                
C        -------------------------------                                
C                                                                       
C*    2.1 ENERGY IS DERIVED FROM UNDERESTIMATED MEAN FREQUENCY.         
C         -----------------------------------------------------         
C                                                                       
      XTEMP = FM*USTT/G                                                 
      IF (XTEMP.LE.0.) WRITE(*,*) FM,USTT,EW                            
      SPINTDI = USTT**4/G**2 * EN(FM*USTT/G)                            
C                                                                       
C*    2.2 MEAN FREQUENCY IS DERIVED FROM THE UNDERESTIMATED ENERGY.     
C         ---------------------------------------------------------     
C                                                                       
      FREQDI = G/USTT * YNU(EW*G**2/USTT**4)                            
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    3. FINAL ESTIMATE.                                                
C        ---------------                                                
C                                                                       
C*    3.1 AVERAGING THE ENERGY ESTIMATES.                               
C         -------------------------------                               
C                                                                       
      EW = (EW+SPINTDI)*0.5                                             
C                                                                       
C*    3.2 AVERAGING THE MEAN FREQUENCY ESTIMATES.                       
C         ---------------------------------------                       
C                                                                       
      FM = (FM+FREQDI)*0.5                                              
                                                                        
      RETURN                                                            
      END                                                               
