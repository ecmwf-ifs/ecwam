      MODULE YOWWNDG

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    ** *WNDGRD* -  INPUT WIND GRID SPECIFICATIONS.

      INTEGER(KIND=JWIM) :: ICODE 
      INTEGER(KIND=JWIM) :: ICODE_CPL 
      INTEGER(KIND=JWIM) :: IWPER 
      INTEGER(KIND=JWIM) :: ICOORD 

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
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
