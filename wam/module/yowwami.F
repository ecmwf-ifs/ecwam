      MODULE YOWWAMI

      IMPLICIT NONE

!     **  *WAMINFO*                                               

      CHARACTER*12 CBEGDT 
      CHARACTER*12 CENDDT 
      CHARACTER*12 CBPLTDT 
      CHARACTER*12 CEPLTDT 
      CHARACTER*12 CLSPDT 
      CHARACTER*12 CRSTDT 
      INTEGER IANALPD 
      INTEGER IFOREPD 
      INTEGER IDELWIN 
      INTEGER IASSIM 
      INTEGER NFCST 
      INTEGER ISTAT(3) 

!     VARIABLE   TYPE     PURPOSE                                  
!     --------   ----     -------                                  
!     IBEGDT    INTEGER   BEGIN DATE OF RUN (YYMMDDHHMM)           
!     IENDDT    INTEGER   END DATE OF RUN (YYMMDDHHMM)             
!     IANALPD   INTEGER   ANALYSIS PERIOD IN SECONDS
!     IFOREPD   INTEGER   FORECAST PERIOD IN SECONDS
!     IDELWIN   INTEGER   WIND TIME STEP IN SECONDS
!     CDATER    CHAR*12   DATE FOR OUTPUT OF BOTH RESTART FILES
!                         IF NOT SET IT WILL BY DEFAULT TO CENDDT 
!     CDATES    CHAR*12   LAST DATE FOR OUTPUT OF RESTART FILES
!                         IF NOT SET IT WILL BY DEFAULT TO CENDDT 
!     CBPLTDT   CHAR*12   BEGIN DATE OF ANALYSIS (YYYYMMDDHHMM)        
!     CEPLTDT   CHAR*12   END DATE OF ANALYSIS                     
!     IASSIM    INTEGER   0 = NO DATA ASSIMILATION                 
!                         1 = DATA ASSIMILATION                    
!                             RECOVERY FILES BLSPXXXX, ETC)        
!     NFCST     INTEGER   0 = FORECAST RECOVERY FILES FETCHED      
!                              1 = ANALYSIS RECOVERY FILES FETCHED      
! ----------------------------------------------------------------------

      END MODULE YOWWAMI
