      SUBROUTINE STATSE (NTOT, SUMX, SUMY, SUMXX, SUMYY, SUMXY, SUMXY2, 
     1                   XMAX, YMAX, XMEAN, YMEAN, STDEVX, STDEVY,      
     2                   DMAX, BIAS, RMSE, STDEVD, SID, RV, CORR,       
     3                   SLOPES, RMSESS, RMSEUS, SLOPE1, RMSES1, RMSEU1,
     4                   SLOPE2, RINTR2, RMSES2, RMSEU2)                
                                                                        
C ----------------------------------------------------------------------
C                                                                       
C**** *STATSE* - EVALUATE STATISTICS FORM CUMULATED SUMS.               
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
C       *CALL* *STATSE (NTOT, SUMX, SUMY, SUMXX, SUMYY, SUMXY, SUMXY2,  
C                      XMAX, YMAX, XMEAN, YMEAN, STDEVX, STDEVY,        
C                      DMAX, BIAS, RMSE, STDEVD, SID, RV, CORR,         
C                      SLOPES, RMSESS, RMSEUS, SLOPE1, RMSES1, RMSEU1,  
C                      SLOPE2, RINTR2, RMSES2, RMSEU2)*                 
C                                                                       
C*    VARIABLE.    TYPE.     PURPOSE.                                   
C     ---------    -------   --------                                   
C       *NTOT*     INTEGER   NUMBER OF POINTS.                          
C       *SUMX*     INTEGER   SUM OF X VALUES.                           
C       *SUMY*     INTEGER   SUM OF Y VALUES.                           
C       *SUMXX*    INTEGER   SUM OF X*X VALUES.                         
C       *SUMYY*    INTEGER   SUM OF Y*Y VALUES.                         
C       *SUMXY*    INTEGER   SUM OF X*Y VALUES.                         
C       *SUMXY2*   INTEGER   SUM OF (Y-X)*(Y-X) VALUES.                 
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
C*    1.0 INITIALISE VALUES.                                            
C         ------------------                                            
C                                                                       
 1000 CONTINUE                                                          
      XMEAN  = 0.                                                       
      YMEAN  = 0.                                                       
      STDEVX = 0.                                                       
      STDEVY = 0.                                                       
      BIAS   = 0.                                                       
      RMSE   = 0.                                                       
      STDEVD = 0.                                                       
      SID    = 0.                                                       
      RV     = 0.                                                       
      CORR   = 0.                                                       
      SLOPES = 0.                                                       
      RMSEUS = 0.                                                       
      RMSESS = 0.                                                       
      SLOPE1 = 0.                                                       
      RMSEU1 = 0.                                                       
      RMSES1 = 0.                                                       
      SLOPE2 = 0.                                                       
      RINTR2 = 0.                                                       
      RMSEU2 = 0.                                                       
      RMSES2 = 0.                                                       
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    2. ENOUGH DATA?                                                   
C        ------------                                                   
C                                                                       
 2000 CONTINUE                                                          
      IF (NTOT.LT.1) RETURN                                             
      XPT  = REAL(NTOT)                                                 
      XPTD = 1./XPT                                                     
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    3. COMPUTE MEANS, LEAST SQUARES FIT AND SCATTER INDEX.            
C        ---------------------------------------------------            
C                                                                       
 3000 CONTINUE                                                          
      XMEAN = SUMX*XPTD                                                 
      YMEAN = SUMY*XPTD                                                 
      BIAS = YMEAN-XMEAN                                                
      RMSE = SQRT(SUMXY2*XPTD)                                          
                                                                        
      IF (NTOT.GT.1) THEN                                               
         XHELP = XPT*SUMXX-SUMX*SUMX                                    
         YHELP = XPT*SUMYY-SUMY*SUMY                                    
                                                                        
         STDEVX = ABS(XHELP/(XPT*(XPT-1.)))                             
         STDEVX = SQRT(STDEVX)                                          
         STDEVY = ABS(YHELP/(XPT*(XPT-1.)))                             
         STDEVY = SQRT(STDEVY)                                          
         STDEVD = ABS((SUMXY2-XPT*BIAS**2)/(XPT-1.))                    
         STDEVD = SQRT(STDEVD)                                          
                                                                        
         IF (XMEAN.NE.0.) SID    = STDEVD/XMEAN                         
         IF (XHELP.GT.0.) RV     = 1. - XPT*SUMXY2/XHELP                
         IF (STDEVX.NE.0. .AND. STDEVY.NE.0.)                           
     1      CORR   = (SUMXY-XPT*XMEAN*YMEAN)/(STDEVX*STDEVY*(XPT-1.))   
         IF (SUMXX.NE.0.) THEN                                          
            SLOPES = SQRT(SUMYY/SUMXX)                                       
            RMSEUS = ABS((SUMYY-2*SLOPES*SUMXY+SLOPES**2*SUMXX)*XPTD)   
            RMSEUS = SQRT(RMSEUS)                                       
            RMSESS = ABS((1.-SLOPES)**2*SUMXX*XPTD)                     
            RMSESS = SQRT(RMSESS)                                       
            SLOPE1 = SUMXY/SUMXX                                        
            RMSEU1 = ABS((SUMYY-2*SLOPE1*SUMXY+SLOPE1**2*SUMXX)*XPTD)   
            RMSEU1 = SQRT(RMSEU1)                                       
            RMSES1 = ABS((1.-SLOPE1)**2*SUMXX*XPTD)                     
            RMSES1 = SQRT(RMSES1)                                       
         ENDIF                                                          
         IF (XHELP.GT.0.) THEN                                          
            SLOPE2 = (XPT*SUMXY-SUMX*SUMY)/XHELP                        
            RINTR2 = YMEAN-SLOPE2*XMEAN                                 
            RMSEU2 = ABS((SUMYY-2*SLOPE2*SUMXY+SLOPE2**2*SUMXX)*XPTD    
     1             - RINTR2**2)                                         
            RMSEU2 = SQRT(RMSEU2)                                       
            RMSES2 = ABS((1.-SLOPE2)**2*SUMXX*XPTD + RINTR2**2          
     1             - 2.*(1.-SLOPE2)*RINTR2*XMEAN)                       
            RMSES2 = SQRT(RMSES2)                                       
         ENDIF                                                          
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
