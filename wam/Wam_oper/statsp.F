      SUBROUTINE STATSP (IU06, NTOT,                                    
     1                   XMAX, YMAX, XMEAN, YMEAN, STDEVX, STDEVY,      
     2                   DMAX, BIAS, RMSE, STDEVD, SID, RV, CORR,       
     3                   SLOPES, RMSESS, RMSEUS, SLOPE1, RMSES1, RMSEU1,
     4                   SLOPE2, RINTR2, RMSES2, RMSEU2)                
                                                                        
C ----------------------------------------------------------------------
C                                                                       
C**** *STATSP*- PRINT STATISTICS.                                       
C                                                                       
C     LIANA ZAMBRESKY      GKSS/ECMWF     AUGUST 1987                   
C     H. GUNTHER           GKSS/ECMWF     AUGUST 1990                   
C                                                                       
C*    PURPOSE.                                                          
C     --------                                                          
C                                                                       
C       TO PRINT THE STATISTIC PARAMETER EVALUATED IN SUB STATSE.       
C                                                                       
C***  INTERFACE.                                                        
C     ----------                                                        
C                                                                       
C       *CALL* *STATSP (IU06, NTOT,                                     
C                      XMAX, YMAX, XMEAN, YMEAN, STDEVX, STDEVY,        
C                      DMAX, BIAS, RMSE, STDEVD, SID, RV, CORR,         
C                      SLOPES, RMSESS, RMSEUS, SLOPE1, RMSES1, RMSEU1,  
C                      SLOPE2, RINTR2, RMSES2, RMSEU2)*                 
C                                                                       
C*    VARIABLE.    TYPE.     PURPOSE.                                   
C     ---------    -------   --------                                   
C       *IU06*     INTEGER   PRINTER OUTPUT UNIT. IF EQ 0 NO OUTPUT.    
C       *NTOT*     INTEGER   NUMBER OF POINTS.                          
C       *XMAX*     REAL      MAX VALUE OF INDEPENDENT VARIABLE          
C       *YMAX*     REAL      MAX VALUE OF DEPENDENT VARIABLE            
C       *XMEAN*    REAL      MEAN VALUE OF INDEPENDENT VARIABLE         
C       *YMEAN*    REAL      MEAN VALUE OF DEPENDENT VARIABLE           
C       *STDEVX*   REAL      STANDARD DEVIATION OF INDEPENDENT VARIABLE 
C       *STDEVY*   REAL      STANDARD DEVIATION OF DEPENDENT VARIABLE   
C       *DMAX*     REAL      MAX ABSOLUTE VALUE OF DIFFERENCE           
C       *BIAS*     REAL      MEAN VALUE OF DIFFERENCES                  
C       *RMSE*     REAL      ROOT-MEAN-SQUARE ERROR                     
C       *STDEVD*   REAL      STANDARD DEVIATION OF DIFFERENCE           
C       *SID*      REAL      SCATTER INDEX                              
C       *RV*       REAL      REDUCTION OF VARIACE.                      
C       *CORR*     REAL      CORRELATION COEFFICIENT                    
C       *SLOPES*   REAL      SLOPE OF SYMMETRIC REGRESSION LINE.        
C       *RMSEUS*   REAL      UNSYSTEMATIC ROOT MEAN SQUARE ERROR        
C                            OF SYMMETRIC REGESSION LINE.               
C       *RMSESS*   REAL      SYSTEMATIC ROOT MEAN SQUARE ERROR          
C                            OF SYMMETRIC REGESSION LINE.               
C       *SLOPE1*   REAL      SLOPE OF REGRESSION LINE THROUGH ORIGIAN.  
C       *RMSEU1*   REAL      UNSYSTEMATIC ROOT MEAN SQUARE ERROR        
C                            OF 1 PARAMETER REGRESSION LINE.            
C       *RMSES1*   REAL      SYSTEMATIC ROOT MEAN SQUARE ERROR          
C                            OF 1 PARAMETER REGRESSION LINE.            
C       *SLOPE2*   REAL      SLOPE OF 2 PARAMETER LEAST SQUARES FIT.    
C       *RINTR2*   REAL      INTERCEPT OF 2 PARAMETER LEAST SQUARES FIT 
C       *RMSEU2*   REAL      UNSYSTEMATIC ROOT MEAN SQUARE ERROR        
C                            OF 2 PARAMETER REGRESSION LINE.            
C       *RMSES2*   REAL      SYSTEMATIC ROOT MEAN SQUARE ERROR          
C                            OF 2 PARAMETER REGRESSION LINE.            
C                                                                       
C     METHOD.                                                           
C     -------                                                           
C                                                                       
C       NONE.                                                           
C                                                                       
C                                                                       
C     EXTERNALS.                                                        
C     ----------                                                        
C                                                                       
C       NONE.                                                           
C                                                                       
C     REFERENCES.                                                       
C     -----------                                                       
C                                                                       
C       ZAMBRESKY,L.F.,  ECMWF TECH REPORT NUMBER 63                    
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    1. CHECK PRINTER UNIT.                                            
C        ------------------                                             
C                                                                       
 1000 CONTINUE                                                          
      IF (IU06.EQ.0) RETURN                                             
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    2. PRINT RESULTS.                                                 
C        --------------                                                 
C                                                                       
 2000 CONTINUE                                                          
      WRITE(IU06,500) NTOT                                              
 500  FORMAT(/1X,'TOTAL NUMBER OF DATA : ',I9)                          
                                                                        
      IF (NTOT.LT.2) THEN                                               
         WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++' 
         WRITE(IU06,*) ' + WARNING: NUMBER OF ENTRIES IS LESS THAN 2 +' 
         WRITE(IU06,*) ' + SOME VALUES COULD NOT BE EVALUATED.       +' 
         WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++' 
      ENDIF                                                             
                                                                        
      WRITE(IU06,501) XMAX, XMEAN, STDEVX                               
 501  FORMAT('0','INDEPENDENT VARIABLE X :',/,'       MAXIMUM VALUE',   
     1  F8.3,'  MEAN VALUE',F8.3,'  STANDARD DEVIATION',F8.3)           
                                                                        
      WRITE(IU06,502) YMAX, YMEAN, STDEVY                               
 502  FORMAT('0','DEPENDENT VARIABLE Y :',/,'       MAXIMUM VALUE',     
     1  F8.3,'  MEAN VALUE',F8.3,'  STANDARD DEVIATION',F8.3)           
                                                                        
      WRITE(IU06,503) DMAX, BIAS , STDEVD                               
 503  FORMAT( 1X,'DIFFERENCE   Y - X :',/,'  ABS. MAXIMUM VALUE',       
     1  F8.3,'        BIAS',F8.3,'  STANDARD DEVIATION',F8.3)           
                                                                        
      WRITE(IU06,504) RMSE, CORR, RV, SID                               
 504  FORMAT( 1X,' ROOT-MEAN-SQUARE ERROR ',F9.4,                       
     1           '  CORRELATION COEFFICIENT ',F9.4,/,                   
     2           '  REDUCTION OF VARIANCE  ',F9.4,                      
     3           '  SCATTER INDEX ',F9.4)                               
                                                                        
      WRITE(IU06,506)                                                   
 506  FORMAT('0','SYMMETRIC REGRESSION LINE : ')                        
      WRITE(IU06,507) SLOPES, RMSESS, RMSEUS                            
 507  FORMAT( 1X,' SLOPE', F7.3,'           ',6X,                       
     1           ' SYMMETRIC RMSE', F7.3,' UNSYMMETRIC RMSE',F7.3)      
                                                                        
      WRITE(IU06,508)                                                   
 508  FORMAT('0','1 PARAMETER REGRESSION LINE : ')                      
      WRITE(IU06,507) SLOPE1, RMSES1, RMSEU1                            
      WRITE(IU06,509)                                                   
                                                                        
 509  FORMAT('0','2 PARAMETER REGRESSION LINE : ')                      
      WRITE(IU06,510) SLOPE2, RINTR2, RMSES2, RMSEU2                    
 510  FORMAT( 1X,' SLOPE', F7.3,' INTERCEPT',F7.3,                      
     1           ' SYMMETRIC RMSE', F7.3,' UNSYMMETRIC RMSE',F7.3)      
                                                                        
      RETURN                                                            
      END                                                               
