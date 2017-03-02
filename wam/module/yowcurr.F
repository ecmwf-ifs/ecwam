      MODULE YOWCURR

      IMPLICIT NONE

!*    ** *CURRENT* - CURRENT FIELD.

      REAL, ALLOCATABLE :: U(:,:) 
      REAL, ALLOCATABLE :: V(:,:) 

      CHARACTER*14 :: CDTCUR
      CHARACTER*14 :: CDATECURA

      INTEGER :: IDELCUR

      LOGICAL :: LLCHKCFL
      LOGICAL :: LLCHKCFLA

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -----     --------
!      *U*         REAL      U - COMPONENT OF CURRENT (M/S).
!      *V*         REAL      V - COMPONENT OF CURRENT (M/S).
!      *CDTCUR*    CHAR*14   REFERENCE DATE FOR U AN V.
!      *CDATECURA* CHAR*14   INITIAL DATE FOR INPUT OF U AND V.
!      *IDELCUR*   INTEGER   CURRENT INPUT STEP IN SECONDS.
!      *LLCHKCFL*  LOGICAL   TRUE IF THE CFL CRITERIA HAVE TO BE CHECKED 
!      *LLCHKCFLA* LOGICAL   TRUE IF THE CFL CRITERIA HAVE TO BE CHECKED
!                            AT ALL TIME STEPS (when currents are updated). 

! ----------------------------------------------------------------------
      END MODULE YOWCURR
