C
C*    *COMMON* *TEMP* TEMPERATURE GRID
C
#if resolution==50
      PARAMETER (NXT =  720, NYT  =  325)
#elif resolution==25 && region=='m'
      PARAMETER (NXT =  561, NYT  =  289)
#else
      PARAMETER (NXT =  240, NYT  =  145)
#endif
      CHARACTER*12 CTDATE
      COMMON /TEMP/ NGXT, NGYT,
     &              TMOWEP, TMOSOP, TMOEAP, TMONOP, DELLA, DELLO,
     &              CTDATE
C
C------------------------------------------------------------------------
C
C*     VARIABLE.   TYPE.     PURPOSE.
C      ---------   -------   --------
C      *DELLA*     REAL      INCREMENT FOR LATITUDE IN TEMP-GRID (DEGREE).
C      *DELLO*     REAL      INCREMENT FOR LONGITUDE IN TEMP-GRID (DEGREE).
C      *CTDATE*    CHAR*10   DATE OF TEMPERATURE GRID
C      *NGXT*      INTEGER   NUMBER OF LONGITUDES IN T-GRID.
C      *NGYT*      INTEGER   NUMBER OF LATITUDES  IN T-GRID.
C      *TMOEAP*    REAL      MOST EASTERN LONGITUDE OF T-GRID (DEGREE).
C      *TMONOP*    REAL      MOST NORTHERN LATITUDE OF T-GRID (DEGREE).
C      *TMOSOP*    REAL      MOST SOUTHERN LATITUDE OF T-GRID (DEGREE).
C      *TMOWEP*    REAL      MOST WESTERN LONGITUDE OF T-GRID (DEGREE).
C
C------------------------------------------------------------------------
C
