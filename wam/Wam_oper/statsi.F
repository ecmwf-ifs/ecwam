      SUBROUTINE STATSI (NTOT, XMAX, YMAX, DMAX,                        
     1                   SUMX, SUMY, SUMXX, SUMYY, SUMXY, SUMXY2)       
                                                                        
C ----------------------------------------------------------------------
C                                                                       
C**** *STATSI* - INITIALIZES STATISTICS.                                
C                                                                       
C      LIANA ZAMBRESKY        GKSS/ECMWF       AUGUST 1987              
C                                                                       
C*    PURPOSE.                                                          
C     --------                                                          
C                                                                       
C          PRESET SUMMATION AND MAXIMA VALUES FOR STATISTICS.           
C                                                                       
C***  INTERFACE.                                                        
C     ----------                                                        
C                                                                       
C     *CALL* *STATSI (NTOT, XMAX, YMAX, DMAX,                           
C                     SUMX, SUMY, SUMXX, SUMYY, SUMXY, SUMXY2)*         
C                                                                       
C*    VARIABLE.    TYPE.     PURPOSE.                                   
C     ---------    -------   --------                                   
C       *NTOT*     INTEGER   NUMBER OF POINTS.                          
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
C       NONE.                                                           
C                                                                       
C     EXTERNALS.                                                        
C     ----------                                                        
C                                                                       
C       NONE                                                            
C                                                                       
C     REFERENCES.                                                       
C     -----------                                                       
C                                                                       
C       NONE.                                                           
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    1. MAXIMA AND TOTAL NUMBER.                                       
C        ------------------------                                       
C                                                                       
 1000 CONTINUE                                                          
      NTOT = 0                                                          
      XMAX = -999999.                                                   
      YMAX = -999999.                                                   
      DMAX = -999999.                                                   
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*    2. A VARIETY OF SUMS.                                             
C        ------------------                                             
C                                                                       
 2000 CONTINUE                                                          
      SUMX    = 0.                                                      
      SUMY    = 0.                                                      
      SUMXX   = 0.                                                      
      SUMYY   = 0.                                                      
      SUMXY   = 0.                                                      
      SUMXY2  = 0.                                                      
                                                                        
      RETURN                                                            
      END                                                               
