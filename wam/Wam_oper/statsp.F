      SUBROUTINE STATSP (IU06, NTOT,                                    
     &                   XMAX, YMAX, XMEAN, YMEAN, STDEVX, STDEVY,      
     &                   DMAX, BIAS, RMSE, STDEVD, SID, RV, CORR,       
     &                   SLOPES, RMSESS, RMSEUS, SLOPE1, RMSES1, RMSEU1,
     &                   SLOPE2, RINTR2, RMSES2, RMSEU2)                
                                                                        
! ----------------------------------------------------------------------

!**** *STATSP*- PRINT STATISTICS.                                       

!     LIANA ZAMBRESKY      GKSS/ECMWF     AUGUST 1987                   
!     H. GUNTHER           GKSS/ECMWF     AUGUST 1990                   

!*    PURPOSE.                                                          
!     --------                                                          

!       TO PRINT THE STATISTIC PARAMETER EVALUATED IN SUB STATSE.       

!***  INTERFACE.                                                        
!     ----------                                                        

!       *CALL* *STATSP (IU06, NTOT,                                     
!                      XMAX, YMAX, XMEAN, YMEAN, STDEVX, STDEVY,        
!                      DMAX, BIAS, RMSE, STDEVD, SID, RV, CORR,         
!                      SLOPES, RMSESS, RMSEUS, SLOPE1, RMSES1, RMSEU1,  
!                      SLOPE2, RINTR2, RMSES2, RMSEU2)*                 

!*    VARIABLE.    TYPE.     PURPOSE.                                   
!     ---------    -------   --------                                   
!       *IU06*     INTEGER   PRINTER OUTPUT UNIT. IF EQ 0 NO OUTPUT.    
!       *NTOT*     INTEGER   NUMBER OF POINTS.                          
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

!       NONE.                                                           


!     EXTERNALS.                                                        
!     ----------                                                        

!       NONE.                                                           

!     REFERENCES.                                                       
!     -----------                                                       

!       ZAMBRESKY,L.F.,  ECMWF TECH REPORT NUMBER 63                    

! ----------------------------------------------------------------------

!*    1. CHECK PRINTER UNIT.                                            
!        ------------------                                             

      IF (IU06.EQ.0) RETURN                                             

! ----------------------------------------------------------------------

!*    2. PRINT RESULTS.                                                 
!        --------------                                                 

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
     &  F8.3,'  MEAN VALUE',F8.3,'  STANDARD DEVIATION',F8.3)           
                                                                        
      WRITE(IU06,502) YMAX, YMEAN, STDEVY                               
 502  FORMAT('0','DEPENDENT VARIABLE Y :',/,'       MAXIMUM VALUE',     
     &  F8.3,'  MEAN VALUE',F8.3,'  STANDARD DEVIATION',F8.3)           
                                                                        
      WRITE(IU06,503) DMAX, BIAS , STDEVD                               
 503  FORMAT( 1X,'DIFFERENCE   Y - X :',/,'  ABS. MAXIMUM VALUE',       
     &  F8.3,'        BIAS',F8.3,'  STANDARD DEVIATION',F8.3)           
                                                                        
      WRITE(IU06,504) RMSE, CORR, RV, SID                               
 504  FORMAT( 1X,' ROOT-MEAN-SQUARE ERROR ',F9.4,                       
     &           '  CORRELATION COEFFICIENT ',F9.4,/,                   
     &           '  REDUCTION OF VARIANCE  ',F9.4,                      
     &           '  SCATTER INDEX ',F9.4)                               
                                                                        
      WRITE(IU06,506)                                                   
 506  FORMAT('0','SYMMETRIC REGRESSION LINE : ')                        
      WRITE(IU06,507) SLOPES, RMSESS, RMSEUS                            
 507  FORMAT( 1X,' SLOPE', F7.3,'           ',6X,                       
     &           ' SYMMETRIC RMSE', F7.3,' UNSYMMETRIC RMSE',F7.3)      
                                                                        
      WRITE(IU06,508)                                                   
 508  FORMAT('0','1 PARAMETER REGRESSION LINE : ')                      
      WRITE(IU06,507) SLOPE1, RMSES1, RMSEU1                            
      WRITE(IU06,509)                                                   
                                                                        
 509  FORMAT('0','2 PARAMETER REGRESSION LINE : ')                      
      WRITE(IU06,510) SLOPE2, RINTR2, RMSES2, RMSEU2                    
 510  FORMAT( 1X,' SLOPE', F7.3,' INTERCEPT',F7.3,                      
     &           ' SYMMETRIC RMSE', F7.3,' UNSYMMETRIC RMSE',F7.3)      
                                                                        
      RETURN                                                            
      END SUBROUTINE STATSP




