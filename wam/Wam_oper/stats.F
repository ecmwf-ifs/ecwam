      SUBROUTINE STATS (IU06, X, Y, NDIM, NPT,                          
     1                  XMAX, YMAX, XMEAN, YMEAN, STDEVX, STDEVY,       
     2                  DMAX, BIAS, RMSE, STDEVD, SID, RV, CORR,        
     3                  SLOPES, RMSESS, RMSEUS, SLOPE1, RMSES1, RMSEU1, 
     4                  SLOPE2, RINTR2, RMSES2, RMSEU2)                 
                                                                        
C ----------------------------------------------------------------------
C                                                                       
C****  *STATS* COMPUTE STATISTICS.                                      
C                                                                       
C     LIANA ZAMBRESKY      GKSS/ECMWF     AUGUST 1987                   
C     H. GUNTHER           GKSS/ECMWF     AUGUST 1990                   
C                                                                       
C*    PURPOSE.                                                          
C     --------                                                          
C                                                                       
C       THIS SUB COMPUTES THE FOLLOWING STATISTICS                      
C                                                                       
C       1.  MAXIMUM VALUES OF X AND Y                                   
C       2.  MEAN AND STANDARD DEVIATION OF BOTH X AND Y VALUES          
C       3.  MAXIMUM VALUE OF ABSOLUTE DIFFERNCE                         
C       4.  MEAN VALUE OF DIFFERENCESM (BIAS)                           
C       5.  ROOT-MEAN-SQUARE ERROR (RMSE)                               
C       6.  STANDARD DEVIATION OF DIFFERENCES                           
C       7.  SCATTER INDEX (SI)                                          
C       8.  REDUCTION OF VARIACE.                                       
C       9.  CORRELATION COEFFICIENT                                     
C      10.  THE SYMMETRIC SLOPE AND THE UNSYMMETRIC AND SYMMETRIC RMSE  
C           OF THIS REGESSION LINE. (SLOPES, RSMEUS, RSMESS)            
C      11.  THE SLOPE OF THE REGRESSION LINE THROUGH THE ORIGAN AND     
C           THE UNSYMMETRIC AND SYMMETRIC RMSE. (SLOPE1, RSMEU1, RSMES1)
C      12.  SLOPE AND INTERCEPT OF THE 2 PARAMETER REGRESSION LINE AND  
C           THE UNSYMMETRIC AND SYMMETRIC RMSE. (SLOPE2, RINTS2, RSMEU2,
C           RSMES2)                                                     
C                                                                       
C***  INTERFACE.                                                        
C     ----------                                                        
C                                                                       
C       *CALL* *STATS (IU06, X, Y, NDIM, NPT,                           
C                      XMAX, YMAX, XMEAN, YMEAN, STDEVX, STDEVY,        
C                      DMAX, BIAS, RMSE, STDEVD, SID, RV, CORR,         
C                      SLOPES, RMSESS, RMSEUS, SLOPE1, RMSES1, RMSEU1,  
C                      SLOPE2, RINTR2, RMSES2, RMSEU2)*                 
C                                                                       
C*    VARIABLE.    TYPE.     PURPOSE.                                   
C     ---------    ------    --------                                   
C       *IU06*     INTEGER   PRINTER OUTPUT UNIT. IF EQ 0 NO OUTPUT.    
C       *X*        REAL      INDEPENDENT ARRAY                          
C       *Y*        REAL      DEPENDENT ARRAY                            
C       *NDIM*     INTEGER   ARRAY DIMENSION                            
C       *NPT*      INTEGER   NUMBER OF POINTS                           
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
C       THE COMPUTATION OF THESE STATISTICS IS STRAIGHTFORWARD.         
C                                                                       
C                                                                       
C     EXTERNALS.                                                        
C     ---------                                                         
C                                                                       
C        *STATSI* - INITIALISE SUMS FOR STATISTICS.                     
C        *STATSC* - COMPUTE    SUMS FOR STATISTICS.                     
C        *STATSE* - EVALUATE   STATISTISC FROM SUMS.                    
C        *STATSP* - PRINT      STATISTICS.                              
C                                                                       
C     REFERENCES.                                                       
C     ----------                                                        
C                                                                       
C          ZAMBRESKY,L.F.,  ECMWF TECH REPORT NUMBER 63                 
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
      DIMENSION X(NDIM),Y(NDIM)                                         
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    1. INITIALIZATION OF SUMS.                                        
C        -----------------------                                        
C                                                                       
 1000 CONTINUE                                                          
      CALL STATSI (NTOT, XMAX, YMAX, DMAX,                              
     1             SUMX, SUMY, SUMXX, SUMYY, SUMXY, SUMXY2)             
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    2. COMPUTE SUMS FOR STATISTICS.                                   
C        -----------------------------                                  
C                                                                       
 2000 CONTINUE                                                          
      CALL STATSC (X, Y, NDIM, NPT, NTOT, XMAX, YMAX, DMAX,             
     1             SUMX, SUMY, SUMXX, SUMYY, SUMXY, SUMXY2)             
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    3. EVALUATE STATISTICS.                                           
C        --------------------                                           
C                                                                       
 3000 CONTINUE                                                          
      CALL STATSE (NTOT, SUMX, SUMY, SUMXX, SUMYY, SUMXY, SUMXY2,       
     1             XMAX, YMAX, XMEAN, YMEAN, STDEVX, STDEVY,            
     2             DMAX, BIAS, RMSE, STDEVD, SID, RV, CORR,             
     3             SLOPES, RMSESS, RMSEUS, SLOPE1, RMSES1, RMSEU1,      
     4             SLOPE2, RINTR2, RMSES2, RMSEU2)                      
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    4. PRINT STATISTICS.                                              
C        -----------------                                              
C                                                                       
 4000 CONTINUE                                                          
      IF (IU06.NE.0) THEN                                               
         CALL STATSP (IU06, NTOT,                                       
     1                XMAX, YMAX, XMEAN, YMEAN, STDEVX, STDEVY,         
     2                DMAX, BIAS, RMSE, STDEVD, SID, RV, CORR,          
     3                SLOPES, RMSESS, RMSEUS, SLOPE1, RMSES1, RMSEU1,   
     4                SLOPE2, RINTR2, RMSES2, RMSEU2)                   
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
