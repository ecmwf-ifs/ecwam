      MODULE YOWTEMP

      IMPLICIT NONE

!*    ** *TEMP* TEMPERATURE GRID

      INTEGER, PARAMETER :: NXT=720 
      INTEGER, PARAMETER :: NYT=325 

      CHARACTER*12 CTDATE 
      INTEGER NGXT 
      INTEGER NGYT 
      REAL TMOWEP 
      REAL TMOSOP 
      REAL TMOEAP 
      REAL TMONOP 
      REAL DELLA 
      REAL DELLO 

!-----------------------------------------------------------------------

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *CTDATE*    CHAR*12   DATE OF TEMPERATURE GRID
!      *NGXT*      INTEGER   NUMBER OF LONGITUDES IN T-GRID.
!      *NGYT*      INTEGER   NUMBER OF LATITUDES  IN T-GRID.
!      *TMOEAP*    REAL      MOST EASTERN LONGITUDE OF T-GRID (DEGREE).
!      *TMONOP*    REAL      MOST NORTHERN LATITUDE OF T-GRID (DEGREE).
!      *TMOSOP*    REAL      MOST SOUTHERN LATITUDE OF T-GRID (DEGREE).
!      *TMOWEP*    REAL      MOST WESTERN LONGITUDE OF T-GRID (DEGREE).
!      *DELLA*     REAL      INCREMENT FOR LATITUDE IN TEMP-GRID (DEG).
!      *DELLO*     REAL      INCREMENT FOR LONGITUDE IN TEMP-GRID (DEG).

!-----------------------------------------------------------------------

      END MODULE YOWTEMP
