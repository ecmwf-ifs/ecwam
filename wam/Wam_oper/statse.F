      SUBROUTINE STATSE (NTOT, SUMX, SUMY, SUMXX, SUMYY, SUMXY, SUMXY2, 
     &                   XMAX, YMAX, XMEAN, YMEAN, STDEVX, STDEVY,      
     &                   DMAX, BIAS, RMSE, STDEVD, SID, RV, CORR,       
     &                   SLOPES, RMSESS, RMSEUS, SLOPE1, RMSES1, RMSEU1,
     &                   SLOPE2, RINTR2, RMSES2, RMSEU2)                
                                                                        
! ----------------------------------------------------------------------

!**** *STATSE* - EVALUATE STATISTICS FORM CUMULATED SUMS.               

!     LIANA ZAMBRESKY      GKSS/ECMWF     AUGUST 1987                   
!     H. GUNTHER           GKSS/ECMWF     AUGUST 1990                   

!*    PURPOSE.                                                          
!     --------                                                          

!       THIS SUB COMPUTES THE FOLLOWING STATISTICS                      

!       1.  MAXIMUM VALUES OF X AND Y                                   
!       2.  MEAN AND STANDARD DEVIATION OF BOTH X AND Y VALUES          
!       3.  MAXIMUM VALUE OF ABSOLUTE DIFFERNCE                         
!       4.  MEAN VALUE OF DIFFERENCESM (BIAS)                           
!       5.  ROOT-MEAN-SQUARE ERROR (RMSE)                               
!       6.  STANDARD DEVIATION OF DIFFERENCES                           
!       7.  SCATTER INDEX (SI)                                          
!       8.  REDUCTION OF VARIACE.                                       
!       9.  CORRELATION COEFFICIENT                                     
!      10.  THE SYMMETRIC SLOPE AND THE UNSYMMETRIC AND SYMMETRIC RMSE  
!           OF THIS REGESSION LINE. (SLOPES, RSMEUS, RSMESS)            
!      11.  THE SLOPE OF THE REGRESSION LINE THROUGH THE ORIGAN AND     
!           THE UNSYMMETRIC AND SYMMETRIC RMSE. (SLOPE1, RSMEU1, RSMES1)
!      12.  SLOPE AND INTERCEPT OF THE 2 PARAMETER REGRESSION LINE AND  
!           THE UNSYMMETRIC AND SYMMETRIC RMSE. (SLOPE2, RINTS2, RSMEU2,
!           RSMES2)                                                     

!***  INTERFACE.                                                        
!     ----------                                                        

!       *CALL* *STATSE (NTOT, SUMX, SUMY, SUMXX, SUMYY, SUMXY, SUMXY2,  
!                      XMAX, YMAX, XMEAN, YMEAN, STDEVX, STDEVY,        
!                      DMAX, BIAS, RMSE, STDEVD, SID, RV, CORR,         
!                      SLOPES, RMSESS, RMSEUS, SLOPE1, RMSES1, RMSEU1,  
!                      SLOPE2, RINTR2, RMSES2, RMSEU2)*                 

!*    VARIABLE.    TYPE.     PURPOSE.                                   
!     ---------    -------   --------                                   
!       *NTOT*     INTEGER   NUMBER OF POINTS.                          
!       *SUMX*     INTEGER   SUM OF X VALUES.                           
!       *SUMY*     INTEGER   SUM OF Y VALUES.                           
!       *SUMXX*    INTEGER   SUM OF X*X VALUES.                         
!       *SUMYY*    INTEGER   SUM OF Y*Y VALUES.                         
!       *SUMXY*    INTEGER   SUM OF X*Y VALUES.                         
!       *SUMXY2*   INTEGER   SUM OF (Y-X)*(Y-X) VALUES.                 
!       *XMAX*     REAL      MAX VALUE OF INDEPENDENT VARIABLE          
!       *YMAX*     REAL      MAX VALUE OF DEPENDENT VARIABLE            
!       *XMEAN*    REAL      MEAN VALUE OF INDEPENDENT VARIABLE         
!       *YMEAN*    REAL      MEAN VALUE OF DEPENDENT VARIABLE           
!       *STDEVX*   REAL      STANDARD DEVIATION OF INDEPENDENT VARIABLE 
!       *STDEVY*   REAL      STANDARD DEVIATION OF DEPENDENT VARIABLE   
!       *DMAX*     REAL      MAX ABSOLUTE VALUE OF DIFFERENCE           
!       *BIAS*     REAL      MEAN VALUE OF DIFFERENCES                  
!       *RMSE*     REAL      ROOT-MEAN-SQUARE ERROR                     
!       *STDEVD*   REAL      STANDARD DEVIATION OF DIFFERENCE           
!       *SID*      REAL      SCATTER INDEX                              
!       *RV*       REAL      REDUCTION OF VARIACE.                      
!       *CORR*     REAL      CORRELATION COEFFICIENT                    
!       *SLOPES*   REAL      SLOPE OF SYMMETRIC REGRESSION LINE.        
!       *RMSEUS*   REAL      UNSYSTEMATIC ROOT MEAN SQUARE ERROR        
!                            OF SYMMETRIC REGESSION LINE.               
!       *RMSESS*   REAL      SYSTEMATIC ROOT MEAN SQUARE ERROR          
!                            OF SYMMETRIC REGESSION LINE.               
!       *SLOPE1*   REAL      SLOPE OF REGRESSION LINE THROUGH ORIGIAN.  
!       *RMSEU1*   REAL      UNSYSTEMATIC ROOT MEAN SQUARE ERROR        
!                            OF 1 PARAMETER REGRESSION LINE.            
!       *RMSES1*   REAL      SYSTEMATIC ROOT MEAN SQUARE ERROR          
!                            OF 1 PARAMETER REGRESSION LINE.            
!       *SLOPE2*   REAL      SLOPE OF 2 PARAMETER LEAST SQUARES FIT.    
!       *RINTR2*   REAL      INTERCEPT OF 2 PARAMETER LEAST SQUARES FIT 
!       *RMSEU2*   REAL      UNSYSTEMATIC ROOT MEAN SQUARE ERROR        
!                            OF 2 PARAMETER REGRESSION LINE.            
!       *RMSES2*   REAL      SYSTEMATIC ROOT MEAN SQUARE ERROR          
!                            OF 2 PARAMETER REGRESSION LINE.            

!     METHOD.                                                           
!     -------                                                           

!       THE COMPUTATION OF THESE STATISTICS IS STRAIGHTFORWARD.         


!     EXTERNALS.                                                        
!     ----------                                                        

!       NONE.                                                           

!     REFERENCES.                                                       
!     -----------                                                       

!       ZAMBRESKY,L.F.,  ECMWF TECH REPORT NUMBER 63                    

! ----------------------------------------------------------------------

!*    1.0 INITIALISE VALUES.                                            
!         ------------------                                            

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

! ----------------------------------------------------------------------

!*    2. ENOUGH DATA?                                                   
!        ------------                                                   

      IF (NTOT.LE.0) THEN
        XMAX=0.
        YMAX=0.
      ENDIF

      IF (NTOT.LT.1) RETURN                                             
      XPT  = REAL(NTOT)                                                 
      XPTD = 1./XPT                                                     

! ----------------------------------------------------------------------

!*    3. COMPUTE MEANS, LEAST SQUARES FIT AND SCATTER INDEX.            
!        ---------------------------------------------------            

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
     &   CORR   = (SUMXY-XPT*XMEAN*YMEAN)/(STDEVX*STDEVY*(XPT-1.))   
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
     &     - RINTR2**2)                                         
          RMSEU2 = SQRT(RMSEU2)                                       
          RMSES2 = ABS((1.-SLOPE2)**2*SUMXX*XPTD + RINTR2**2          
     &     - 2.*(1.-SLOPE2)*RINTR2*XMEAN)                       
          RMSES2 = SQRT(RMSES2)                                       
        ENDIF                                                          
      ENDIF
                                                                        
      RETURN                                                            
      END SUBROUTINE STATSE




