      SUBROUTINE STATS (IU06, X, Y, NDIM, NPT,                          
     &                  XMAX, YMAX, XMEAN, YMEAN, STDEVX, STDEVY,       
     &                  DMAX, BIAS, RMSE, STDEVD, SID, RV, CORR,        
     &                  SLOPES, RMSESS, RMSEUS, SLOPE1, RMSES1, RMSEU1, 
     &                  SLOPE2, RINTR2, RMSES2, RMSEU2)                 
                                                                        
! ----------------------------------------------------------------------

!****  *STATS* COMPUTE STATISTICS.                                      

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

!       *CALL* *STATS (IU06, X, Y, NDIM, NPT,                           
!                      XMAX, YMAX, XMEAN, YMEAN, STDEVX, STDEVY,        
!                      DMAX, BIAS, RMSE, STDEVD, SID, RV, CORR,         
!                      SLOPES, RMSESS, RMSEUS, SLOPE1, RMSES1, RMSEU1,  
!                      SLOPE2, RINTR2, RMSES2, RMSEU2)*                 

!*    VARIABLE.    TYPE.     PURPOSE.                                   
!     ---------    ------    --------                                   
!       *IU06*     INTEGER   PRINTER OUTPUT UNIT. IF EQ 0 NO OUTPUT.    
!       *X*        REAL      INDEPENDENT ARRAY                          
!       *Y*        REAL      DEPENDENT ARRAY                            
!       *NDIM*     INTEGER   ARRAY DIMENSION                            
!       *NPT*      INTEGER   NUMBER OF POINTS                           
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
!     ---------                                                         

!        *STATSI* - INITIALISE SUMS FOR STATISTICS.                     
!        *STATSC* - COMPUTE    SUMS FOR STATISTICS.                     
!        *STATSE* - EVALUATE   STATISTISC FROM SUMS.                    
!        *STATSP* - PRINT      STATISTICS.                              

!     REFERENCES.                                                       
!     ----------                                                        

!          ZAMBRESKY,L.F.,  ECMWF TECH REPORT NUMBER 63                 

! ----------------------------------------------------------------------

      DIMENSION X(NDIM),Y(NDIM)                                         

! ----------------------------------------------------------------------

!*    1. INITIALIZATION OF SUMS.                                        
!        -----------------------                                        

      CALL STATSI (NTOT, XMAX, YMAX, DMAX,                              
     &             SUMX, SUMY, SUMXX, SUMYY, SUMXY, SUMXY2)             

! ----------------------------------------------------------------------

!*    2. COMPUTE SUMS FOR STATISTICS.                                   
!        -----------------------------                                  

      CALL STATSC (X, Y, NDIM, NPT, NTOT, XMAX, YMAX, DMAX,             
     &             SUMX, SUMY, SUMXX, SUMYY, SUMXY, SUMXY2)             

! ----------------------------------------------------------------------

!*    3. EVALUATE STATISTICS.                                           
!        --------------------                                           

      CALL STATSE (NTOT, SUMX, SUMY, SUMXX, SUMYY, SUMXY, SUMXY2,       
     &             XMAX, YMAX, XMEAN, YMEAN, STDEVX, STDEVY,            
     &             DMAX, BIAS, RMSE, STDEVD, SID, RV, CORR,             
     &             SLOPES, RMSESS, RMSEUS, SLOPE1, RMSES1, RMSEU1,      
     &             SLOPE2, RINTR2, RMSES2, RMSEU2)                      

! ----------------------------------------------------------------------

!*    4. PRINT STATISTICS.                                              
!        -----------------                                              

      IF (IU06.NE.0) THEN                                               
        CALL STATSP (IU06, NTOT,                                       
     &   XMAX, YMAX, XMEAN, YMEAN, STDEVX, STDEVY,         
     &   DMAX, BIAS, RMSE, STDEVD, SID, RV, CORR,          
     &   SLOPES, RMSESS, RMSEUS, SLOPE1, RMSES1, RMSEU1,   
     &   SLOPE2, RINTR2, RMSES2, RMSEU2)                   
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END SUBROUTINE STATS




