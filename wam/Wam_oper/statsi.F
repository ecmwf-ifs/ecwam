      SUBROUTINE STATSI (NTOT, XMAX, YMAX, DMAX,                        
     &                   SUMX, SUMY, SUMXX, SUMYY, SUMXY, SUMXY2)       
                                                                        
! ----------------------------------------------------------------------

!**** *STATSI* - INITIALIZES STATISTICS.                                

!      LIANA ZAMBRESKY        GKSS/ECMWF       AUGUST 1987              

!*    PURPOSE.                                                          
!     --------                                                          

!          PRESET SUMMATION AND MAXIMA VALUES FOR STATISTICS.           

!***  INTERFACE.                                                        
!     ----------                                                        

!     *CALL* *STATSI (NTOT, XMAX, YMAX, DMAX,                           
!                     SUMX, SUMY, SUMXX, SUMYY, SUMXY, SUMXY2)*         

!*    VARIABLE.    TYPE.     PURPOSE.                                   
!     ---------    -------   --------                                   
!       *NTOT*     INTEGER   NUMBER OF POINTS.                          
!       *XMAX*     REAL      MAX VALUE OF INDEPENDENT VARIABLE.         
!       *YMAX*     REAL      MAX VALUE OF DEPENDENT VARIABLE.           
!       *DMAX*     REAL      MAX ABSOLUTE VALUE OF DIFFERENCE.          
!       *SUMX*     REAL      SUM OF X VALUES.                           
!       *SUMY*     REAL      SUM OF Y VALUES.                           
!       *SUMXX*    REAL      SUM OF X*X VALUES.                         
!       *SUMYY*    REAL      SUM OF Y*Y VALUES.                         
!       *SUMXY*    REAL      SUM OF X*Y VALUES.                         
!       *SUMXY2*   REAL      SUM OF (Y-X)*(Y-X) VALUES.                 

!     METHOD.                                                           
!     -------                                                           

!       NONE.                                                           

!     EXTERNALS.                                                        
!     ----------                                                        

!       NONE                                                            

!     REFERENCES.                                                       
!     -----------                                                       

!       NONE.                                                           

! ----------------------------------------------------------------------

!*    1. MAXIMA AND TOTAL NUMBER.                                       
!        ------------------------                                       

      NTOT = 0                                                          
      XMAX = -999999.                                                   
      YMAX = -999999.                                                   
      DMAX = -999999.                                                   

! ----------------------------------------------------------------------

!*    2. A VARIETY OF SUMS.                                             
!        ------------------                                             

      SUMX    = 0.                                                      
      SUMY    = 0.                                                      
      SUMXX   = 0.                                                      
      SUMYY   = 0.                                                      
      SUMXY   = 0.                                                      
      SUMXY2  = 0.                                                      
                                                                        
      RETURN                                                            
      END SUBROUTINE STATSI




