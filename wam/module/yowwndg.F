      MODULE YOWWNDG

      IMPLICIT NONE

!*    ** *WNDGRD* -  INPUT WIND GRID SPECIFICATIONS.

      REAL :: DLAM 
      REAL :: DPHI 
      REAL :: RLATS 
      REAL :: RLATN 
      REAL :: RLONL 
      REAL :: RLONR 
      INTEGER :: KCOL 
      INTEGER :: KROW 
      INTEGER :: ICODE 
      INTEGER :: ICODE_CPL 
      INTEGER :: IWPER 
      INTEGER :: ICOORD 

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *DLAM*      REAL      STEPSIZE BETWEEN LONGITUDES IN DEG.
!      *DPHI*      REAL      STEPSIZE BETWEEN LATITUDES  IN DEG.
!      *RLATS*     REAL      LATITUDE  AT (., 1) = SOUTHERN LATITUDE.
!      *RLATN*     REAL      LATITUDE  AT (.,NR) = NORTHERN LATITUDE.
!      *RLONL*     REAL      LONGITUDE AT ( 1,.) = WEST MOST LONGITUDE.
!      *RLONR*     REAL      LONGITUDE AT (NC,.) = EAST MOST LONGITUDE.
!      *KCOL*      INTEGER   NUMBER OF COLUMNES IN WIND INPUT (USED).
!      *KROW*      INTEGER   NUMBER OF ROWS     IN WIND INPUT (USED).
!      *ICODE*     INTEGER   WIND CODE FROM INPUT FILE
!                            1 = USTAR;  2 = USTRESS; 3 = U10
!      *ICODE_CPL* INTEGER   WIND CODE PASSED VIA WAVEMDL
!                            0 = NOT SET; 1 = USTAR;  2 = USTRESS; 3 = U10
!      *IWPER*     INTEGER   INDICATOR PERIODICAL GRID.
!                            0= NON-PERIODICAL;   1= PERIODICAL.
!      *ICOORD*    INTEGER   CODE FOR COORDINATE SYSTEM USED
!                            1= RECTANGULAR,EQUIDISTANT LON/LAT GRID.
!                            2= .......NOT IMPLEMENTED.
! ----------------------------------------------------------------------

      END MODULE YOWWNDG
