      SUBROUTINE STATSC (X, Y, NDIM, NPT, NTOT, XMAX, YMAX, DMAX,       
     1                   SUMX, SUMY, SUMXX, SUMYY, SUMXY, SUMXY2)       
                                                                        
C ----------------------------------------------------------------------
C                                                                       
C**** *STATSC* - COMPUTE SUMS FOR STATISTICS.                           
C                                                                       
C     LIANA ZAMBRESKY        GKSS/ECMWF       AUGUST 1987               
C     H. GUNTHER             GKSS/ECMWF       AUGUST 1990               
C                                                                       
C*    PURPOSE.                                                          
C     --------                                                          
C                                                                       
C       THIS SUB ADDS TO THE STATISTIC SUMS THE VALUES                  
C       FROM THE INPUT FIELDS.                                          
C                                                                       
C       1.  MAXIMUM VALUES OF X AND Y                                   
C       2.  MAXIMUM VALUE OF ABSOLUTE DIFFERNCE                         
C       3.  SUM OF X VALUES.                                            
C       4.  SUM OF Y VALUES.                                            
C       5.  SUM OF X*X VALUES.                                          
C       6.  SUM OF Y*Y VALUES.                                          
C       7.  SUM OF X*Y VALUES.                                          
C       8.  SUM OF (Y-X)*(Y-X) VALUES.                                  
C       9.  TOTAL NUMBER OF DATA PAIRS.                                 
C                                                                       
C                                                                       
C***  INTERFACE.                                                        
C     ----------                                                        
C                                                                       
C       *CALL* *STATSC (X, Y, NDIM, NPT, NTOT, XMAX, YMAX, DMAX,        
C                       SUMX, SUMY, SUMXX, SUMYY, SUMXY, SUMXY2)*       
C                                                                       
C*    VARIABLE.    TYPE.     PURPOSE.                                   
C     ---------    -------   --------                                   
C       *X*        REAL      INDEPENDENT ARRAY.                         
C       *Y*        REAL      DEPENDENT ARRAY.                           
C       *NDIM*     INTEGER   ARRAY DIMENSION.                           
C       *NPT*      INTEGER   NUMBER OF POINTS IN INPUT.                 
C       *NTOT*     INTEGER   NUMBER OF POINTS IN STATISTIC SUMS.        
C       *XMAX*     REAL      MAX VALUE OF INDEPENDENT VARIABLE.         
C       *YMAX*     REAL      MAX VALUE OF DEPENDENT VARIABLE.           
C       *DMAX*     REAL      MAX ABSOLUTE VALUE OF DIFFERENCE.          
C       *SUMX*     REAL      SUM OF X VALUES.                           
C       *SUMY*     REAL      SUM OF Y VALUES.                           
C       *SUMXX*    REAL      SUM OF X*X VALUES.                         
C       *SUMYY*    REAL      SUM OF Y*Y VALUES.                         
C       *SUMXY*    REAL      SUM OF X*Y VALUES.                         
C       *SUMXY2*   REAL      SUM OF (Y-X)*(Y-X) VALUES.                 
C                                                                       
C     METHOD.                                                           
C     -------                                                           
C                                                                       
C       THE COMPUTATION IS STRAIGHTFORWARD.                             
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
      DIMENSION X(NDIM), Y(NDIM)                                        
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    1. MAXIMA AND TOTAL NUMBER.                                       
C        ------------------------                                       
C                                                                       
 1000 CONTINUE                                                          
      NTOT = NTOT+NPT                                                   
      DO 1001 I=1,NPT                                                   
         XMAX = MAX(X(I),XMAX)                                          
         YMAX = MAX(Y(I),YMAX)                                          
         DMAX = MAX(ABS(Y(I)-X(I)),DMAX)                                
 1001 CONTINUE                                                          
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    2. FORM A VARIETY OF SUMS.                                        
C        -----------------------                                        
C                                                                       
 2000 CONTINUE                                                          
      DO 2001 I=1,NPT                                                   
         SUMX    = SUMX+X(I)                                            
         SUMY    = SUMY+Y(I)                                            
         SUMXX   = SUMXX+X(I)*X(I)                                      
         SUMYY   = SUMYY+Y(I)*Y(I)                                      
         SUMXY   = SUMXY+X(I)*Y(I)                                      
         SUMXY2  = SUMXY2+(Y(I)-X(I))*(Y(I)-X(I))                       
 2001 CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
