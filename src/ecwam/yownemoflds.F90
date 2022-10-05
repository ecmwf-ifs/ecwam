      MODULE YOWNEMOFLDS

      USE YOWDRVTYPE  , ONLY : WAVE2OCEAN, OCEAN2WAVE

      IMPLICIT NONE

!*     ** *NEMOFLDS* NEMO FIELDS FOR COUPLED RUNS

      TYPE(WAVE2OCEAN), ALLOCATABLE :: WAM2NEMO(:,:) 
      TYPE(OCEAN2WAVE), ALLOCATABLE :: NEMO2WAM(:,:) 

      LOGICAL :: LNEMOCITHICK, LNEMOICEREST

!--------------------------------------------------------------------

!*    VARIABLE     TYPE      PURPOSE
!     --------     ----      -------
!     LNEMOCITHICK LOGICAL   SET TO TRUE IF SEA ICE THICKNESS IS 
!                            AVAILABLE FROM NEMO (E.G. LIM2 ACTIVE).
!     LNEMOICEREST LOGICAL   SET TO TRUE IF SEA ICE IS NOT RESCALED BY ICE COVER
!---------------------------------------------------------------------
      END MODULE YOWNEMOFLDS
