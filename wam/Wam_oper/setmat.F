      SUBROUTINE SETMAT (XM, NOBS, P, D, NDIM2)                         
C                                                                       
C --------------------------------------------------------------------- 
C                                                                       
C     DETERMINE MATRIX P+D                                              
C                                                                       
      DIMENSION XM(NOBS,NOBS), P(NDIM2,NDIM2), D(NDIM2,NDIM2)           
C                                                                       
      DO 1000 JOBS=1,NOBS                                               
         DO 1100 IOBS=1,NOBS                                            
            XM(IOBS,JOBS)=P(IOBS,JOBS)+D(IOBS,JOBS)                     
 1100    CONTINUE                                                       
 1000 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
