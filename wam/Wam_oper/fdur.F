      SUBROUTINE FDUR (EWFG, T, USMO, IJS, IJL)                         
C                                                                       
C--------------------------------------------------------------------   
C                                                                       
C**** *FDUR* - EVALUATE THE WINDSEA DURATION FROM THE MODEL GROTH CURVE.
C                                                                       
C     P.LIONELLO     ECMWF       APRIL 1990                             
C                                                                       
C     PURPOSE.                                                          
C     --------                                                          
C                                                                       
C       EVALUATE THE WINDSEA DURATION FROM THE MODEL GROTH CURVE.       
C                                                                       
C**   INTERFACE.                                                        
C     ----------                                                        
C                                                                       
C       *CALL* *FDUR (EWFG, T, USMO, IJS, IJL)*                         
C                                                                       
C        *EWFG*   REAL     WINDSEA ENERGY (FIRST GUESS).                
C        *T*      REAL     WINDSEA DURATION (OUTPUT).                   
C        *USMO*   REAL     FIRST GUESS USTAR.                           
C        *IJS*    INTEGER  FIRST INDEX IN BLOCK.                        
C        *IJL*    INTEGER  LAST  INDEX IN BLOCK.                        
C                                                                       
C     METHOD.                                                           
C     -------                                                           
C                                                                       
C        THE DURATION IS DERIVED BY THE GROWTH CURVE                    
C                ESTAR=EO*TANH(A*TSTAR**B).                             
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
#include "param.h"
C 
C
#include "commpp.h"
C
      REAL,DIMENSION(NINF:NSUP) :: EWFG, T, USMO 
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
#include "parcons.h"
C                                                                       
      PARAMETER (EO = 955., A = 6.02*10.**(-5), B = 0.695)              
C                                                                       
C       *EO*      REAL    PARAMETER OF THE ENERGY GROWTH CURVE.         
C       *A*       REAL    PARAMETER OF THE ENERGY GROWTH CURVE.         
C       *B*       REAL    PARAMETER OF THE ENERGY GROWTH CURVE.         
C       *G*       REAL    ACCELLERATION OF GRAVITY.                     
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    1. INITIALIZE DURATION.                                           
C        -------------------                                            
C                                                                       
 1000 CONTINUE                                                          
      DO 1001 IJ=IJS,IJL                                                
         T(IJ) = -9.                                                    
 1001 CONTINUE                                                          
C                                                                       
C*    2. LOOP OVER POINTS IN BLOCK.                                     
C        --------------------------                                     
C                                                                       
 2000 CONTINUE                                                          
      DO 2001 IJ=IJS,IJL                                                
         IF (EWFG(IJ).GT.0.01) THEN                                     
C                                                                       
C        DETERMINE THE EFFECTIVE DURATION OF THE WINDSEA ACCORDING TO   
C        THE WAM MODEL GROWTH CURVE, IF THERE WINDSEA IS PRESENT.       
C                                                                       
            TSTAR = EWFG(IJ)*G**2/USMO(IJ)**4                           
            TSTAR = MIN (TSTAR/EO, .999999)                             
            TSTAR = (1.+TSTAR)/(1.-TSTAR)                               
            TSTAR = (ALOG(TSTAR)/(2.*A))**(1./B)                        
            T(IJ) = USMO(IJ)/G*TSTAR                                    
         ENDIF                                                          
 2001 CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
