      SUBROUTINE STATSC (X, Y, NDIM, NPT, NTOT, XMAX, YMAX, DMAX,       
     &                   SUMX, SUMY, SUMXX, SUMYY, SUMXY, SUMXY2)       
                                                                        
! ----------------------------------------------------------------------

!**** *STATSC* - COMPUTE SUMS FOR STATISTICS.                           

!     LIANA ZAMBRESKY        GKSS/ECMWF       AUGUST 1987               
!     H. GUNTHER             GKSS/ECMWF       AUGUST 1990               

!*    PURPOSE.                                                          
!     --------                                                          

!       THIS SUB ADDS TO THE STATISTIC SUMS THE VALUES                  
!       FROM THE INPUT FIELDS.                                          

!       1.  MAXIMUM VALUES OF X AND Y                                   
!       2.  MAXIMUM VALUE OF ABSOLUTE DIFFERNCE                         
!       3.  SUM OF X VALUES.                                            
!       4.  SUM OF Y VALUES.                                            
!       5.  SUM OF X*X VALUES.                                          
!       6.  SUM OF Y*Y VALUES.                                          
!       7.  SUM OF X*Y VALUES.                                          
!       8.  SUM OF (Y-X)*(Y-X) VALUES.                                  
!       9.  TOTAL NUMBER OF DATA PAIRS.                                 


!***  INTERFACE.                                                        
!     ----------                                                        

!       *CALL* *STATSC (X, Y, NDIM, NPT, NTOT, XMAX, YMAX, DMAX,        
!                       SUMX, SUMY, SUMXX, SUMYY, SUMXY, SUMXY2)*       

!*    VARIABLE.    TYPE.     PURPOSE.                                   
!     ---------    -------   --------                                   
!       *X*        REAL      INDEPENDENT ARRAY.                         
!       *Y*        REAL      DEPENDENT ARRAY.                           
!       *NDIM*     INTEGER   ARRAY DIMENSION.                           
!       *NPT*      INTEGER   NUMBER OF POINTS IN INPUT.                 
!       *NTOT*     INTEGER   NUMBER OF POINTS IN STATISTIC SUMS.        
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

!       THE COMPUTATION IS STRAIGHTFORWARD.                             

!     EXTERNALS.                                                        
!     ----------                                                        

!       NONE.                                                           

!     REFERENCES.                                                       
!     -----------                                                       

!       ZAMBRESKY,L.F.,  ECMWF TECH REPORT NUMBER 63                    

! ----------------------------------------------------------------------

      DIMENSION X(NDIM), Y(NDIM)                                        

! ----------------------------------------------------------------------

!*    1. MAXIMA AND TOTAL NUMBER.                                       
!        ------------------------                                       

      NTOT = NTOT+NPT                                                   
      DO I=1,NPT                                                   
        XMAX = MAX(X(I),XMAX)                                          
        YMAX = MAX(Y(I),YMAX)                                          
        DMAX = MAX(ABS(Y(I)-X(I)),DMAX)                                
      ENDDO

! ----------------------------------------------------------------------

!*    2. FORM A VARIETY OF SUMS.                                        
!        -----------------------                                        

      DO I=1,NPT                                                   
        SUMX    = SUMX+X(I)                                            
        SUMY    = SUMY+Y(I)                                            
        SUMXX   = SUMXX+X(I)*X(I)                                      
        SUMYY   = SUMYY+Y(I)*Y(I)                                      
        SUMXY   = SUMXY+X(I)*Y(I)                                      
        SUMXY2  = SUMXY2+(Y(I)-X(I))*(Y(I)-X(I))                       
      ENDDO
                                                                        
      RETURN                                                            
      END SUBROUTINE STATSC




