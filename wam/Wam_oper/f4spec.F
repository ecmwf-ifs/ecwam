      SUBROUTINE F4SPEC (ALFA, FM, GAMMA, SA, THETAQ, USTAR, SPINT, SPT)
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
C**** *F4SPEC* - PRODUCE A WINDSEA SPECTRUM FOR ASSIMILATION.           
C                                                                       
C      P. LIONELLO    ECMWF        FEBRUARY 1990                        
C      (DERIVED FROM JSPEC OF PETER JANSSEN  1987.)                     
C                                                                       
C**   INTERFACE.                                                        
C     ----------                                                        
C                                                                       
C       *CALL* *F4SPEC (ALFA, FM, GAMMA, SA, THETAQ, USTAR, SPINT,SPT)* 
C                                                                       
C          *ALFA*  - TOBAS' CONSTANT                                    
C          *FM*    - PEAK FREQUENCY                                     
C          *GAMMA* - OVERSHOOT PARAMETER                                
C          *SA*    - WIDTH PARAMETER                                    
C          *THETAQ*- MEAN WINDSEA DIRECTION                             
C          *USTAR* - FRICTION VELOCITY                                  
C          *SPINT* - TOTAL ENERGY                                       
C          *SPT*   - 2-D SPECTRUM                                       
C                                                                       
C      METHOD.                                                          
C      -------                                                          
C                                                                       
C      THE ONE DIMENSIONAL PARAMETRIC SPECTRUM IS EVALUATED AND         
C      DISTRIBUTED AROUND THE MEAN SEA DIRECTION.                       
C      THE TOTAL ENERGY IS COMPUTED.                                    
C                                                                       
C      EXTERNALS.                                                       
C      ----------                                                       
C                                                                       
C         *SPR*  -  SUB TO COMPUTE COS**2 SPREADING FUNCTION.           
C                                                                       
C --------------------------------------------------------------------- 
C                                                                       
#include "param.h"
      PARAMETER (ML = NFRE, KL = NANG)                                  
C                                                                       
#include "comfred.h"
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
      DIMENSION SPT(NANG,NFRE)                                          
      DIMENSION DELFR(NFRE), ET(NFRE), ST(NANG)                         
C                                                                       
C      *SPT*        REAL    2-D SPECTRUM.                               
C      *DELFR*      REAL    FREQUENCY INCREMENTS.                       
C      *ET*         REAL    1-D JONSWAP SPECTRUM.                       
C      *ST*         REAL    SPREADING FUNCTION.                         
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
#include "parcons.h"
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
C*    1. COMPUTE FREQUENCY INTERVALLS.                                  
C        -----------------------------                                  
C                                                                       
 1000 CONTINUE                                                          
      DO 1001 M=1,ML                                                    
         DELFR(M) = FR(M)*0.1                                           
 1001 CONTINUE                                                          
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
C*    2. COMPUTE TOBA SPECTRUM.                                         
C        ----------------------                                         
C                                                                       
 2000 CONTINUE                                                          
      FAK = ALFA*G*USTAR/(ZPI**3)                                       
      DO 2001 M=1,ML                                                    
         FRH = FR(M)                                                    
         EARG = (.5*((FRH-FM)/(SA*FM))**2)                              
         FJON = GAMMA**EXP(-EARG)                                       
         FMPF = 4./6.*(FM/FRH)**6                                       
         FJONH = EXP(-FMPF)                                             
         ET(M) = FAK*FJONH*FJON/FRH**4                                  
 2001 CONTINUE                                                          
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
C*  3. COMPUTATION OF SPREADING FUNCTION.                               
C      ----------------------------------                               
C                                                                       
 3000 CONTINUE                                                          
      CALL SPR (KL, NANG, THETAQ, TH, ST)                               
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
C*  4. COMPUTATION OF 2-D SPECTRUM.                                     
C      ---------------------------                                      
C                                                                       
 4000 CONTINUE                                                          
      DO 4001 M=1,ML                                                    
         DO 4002 K=1,KL                                                 
            SPT(K,M)=ET(M)*ST(K)                                        
 4002    CONTINUE                                                       
 4001 CONTINUE                                                          
C                                                                       
C --------------------------------------------------------------------- 
C                                                                       
C*  5. COMPUTE TOTAL ENERGY.                                            
C      --------------------                                             
C                                                                       
 5000 CONTINUE                                                          
C                                                                       
C*  5.1 INTEGRATE OVER DIRECTIONS.                                      
C       --------------------------                                      
C                                                                       
      KL2 = KL-2                                                        
      SPINT = 0.                                                        
      DO 5001 M=1,ML                                                    
         ET(M) = SPT(1,M)                                               
         DO 5002 K=2,KL2,2                                              
            ET(M) = ET(M)+4.*SPT(K,M)+2.*SPT(K+1,M)                     
 5002    CONTINUE                                                       
         ET(M) = (ET(M)+4.*SPT(KL,M)+SPT(1,M))*DELTH/3.                 
 5001 CONTINUE                                                          
C                                                                       
C*  5.2 INTEGRATE 1-D SPECTRUM OVER FREQUENCIES.                        
C       ----------------------------------------                        
C                                                                       
      SPINT = ET(1)*DELFR(1)                                            
      ML2 = ML-2                                                        
      DO 5003 M=2,ML2,2                                                 
         SPINT = SPINT+(4.*ET(M)*DELFR(M)+2.*ET(M+1)*DELFR(M+1))        
 5003 CONTINUE                                                          
      IF (2*(ML/2).EQ.ML) THEN                                          
         SPINT = (SPINT+4.*ET(ML)*DELFR(ML))/3.                         
      ELSE                                                              
         SPINT = (SPINT+4.*ET(ML-1)*DELFR(ML-1)                         
     1         + ET(ML)*DELFR(ML))/3.                                   
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
