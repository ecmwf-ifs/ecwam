      SUBROUTINE ANALYSE (XM, NOBS, XME, XMO, XOI, IMEAS, JMEAS, I, J,  Baae
     .                    DIST, LMAX, NDIM1, NDIM2) 
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C**** *ANALYSE* - COMPUTE FIELD AT (I,J) BY OPTIMUM INTERPOLATION.      
C                                                                       
C     PIERO LIONELLO      ECMWF     JUNE 1990                           
C                                                                       
C     PURPOSE.                                                          
C     --------                                                          
C                                                                       
C       COMPUTE THE FIELD IN THE POINT (I,J) BY OPTIMUM INTERPOLATION.  
C       ANALYSIS IS DONE FOR ONE GRID POINT.                            
C                                                                       
C**   INTERFACE.                                                        
C     ----------                                                        
C                                                                       
C       *CALL* *ANALYSE (XM, NOBS, XME, XMO, XOI, IMEAS, JMEAS, I, J)*  
C                                                                       
C         *XM*      REAL      INVERSE OF (P+M) : P=MODEL CORRELATION    
C                             MATRIX , M=MEASUREMENTS CORRELATION MATRIX
C         *NOBS*    INTEGER   NUMBER OF OBSERVATION IN THE              
C                             BOX AROUND POINT (I,J)                    
C         *XME*     REAL      FIELD FROM  MEASUREMENT                   
C         *XMO*     REAL      FIELD FROM MODEL                          
C         *XOI*     REAL      RESULT OF O.I.                            
C         *IMEAS*   INTEGER   X-INDICES OF MEASUREMENTS                 
C         *JMEAS*   INTEGER   Y-INDICES OF MEASUREMENTS                 
C         *I*       INTEGER   POINT WHERE O.I.FIELD IS PRODUCED         
C         *J*       INTEGER                                             
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
#include "param.h"
C
#include "parcons.h"
C
#include "comgrid.h"
C
#include "commap.h"
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
      DIMENSION XME(NGX,NGY), XMO(NGX,NGY), IMEAS(NDIM2), JMEAS(NDIM2)  
      DIMENSION XM(NOBS,NOBS)                                           
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
      DO L=1,NOBS-1                                                
        DO K=L+1,NOBS                                                
          XM(L,K) = XM(K,L)                                              
        ENDDO
      ENDDO
C                                                                       
      XNEW = 0.                                                         
      DO 1200 IOBS=1,NOBS                                               
         I1 = I
         I2 = IMEAS(IOBS)
         J1 = NGY - J + 1
         J2 = NGY - JMEAS(IOBS) + 1
C                  
         IF((I1.EQ.I2).AND.(J1.EQ.J2)) THEN
           DOBS = 0.
           W = 1.
         ELSE
           DELLON  = FLOAT(I1-1)*ZDELLO(J1)
     1               - FLOAT(I2-1)*ZDELLO(J2)
           COSLON  = COS(DELLON*RAD)
           DOBS2   = COSLON*COSPH(J1)*COSPH(J2) +
     1                     SINPH(J1)*SINPH(J2)
           DOBS    = ACOS(DOBS2)
           W       = EXP(-DOBS/DIST)
         ENDIF

         DO 1300 JOBS=1,NOBS                                            
            IX = IMEAS(JOBS)                                            
            JX = JMEAS(JOBS)                                            
            XNEW = XNEW + W*XM(IOBS,JOBS)*(XME(IX,JX)-XMO(IX,JX))       
 1300    CONTINUE                                                       
 1200 CONTINUE                                                          
                                                                        
      XOI = XMO(I,J) + XNEW                                             
C                                                                       
      RETURN                                                            
      END                                                               
