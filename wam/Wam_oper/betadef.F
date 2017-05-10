      REAL FUNCTION BETADEF(U10)
! ----------------------------------------------------------------------

!**** *BETADEF* - RETURNS DEFAULT VALUE OF THE CHARNOCK PARAMETER
!                 FOR GIVEN 10M NEUTRAL WIND.

!     JEAN BIDLOT

!**   INTERFACE.
!     ----------

!       *BETADEF (U10)*
!          *U10*     - 10M NEUTRAL WINDS.

!     METHOD.
!     -------

!     A QUADRATIC FIT OF THE MEAN CHARNOCK v U10 WAS DETERMINED USING
!     ONE YEAR OF MONTHLY DATA FROM THE WAVE MODEL ANALYSIS (JULY 2012 to JUNE 2013)
!     IT IS ONLY VALID UP TO U10MAX. FOR U10 VALUES LARGER THAN U10MAX,
!     THE VALUE FOR U10MAX IS ASSUMED.

!     EXTERNALS.
!     ----------

!     REFERENCES.                                                       
!     -----------  

! ----------------------------------------------------------------------
      IMPLICIT NONE

!     THESE VALUES CORRESPOND TO CY38R1 WAVE MODEL PHYSICS:
      REAL, PARAMETER :: A=0.00594770
      REAL, PARAMETER :: B=0.00054409
      REAL, PARAMETER :: C=0.00004055
      REAL, PARAMETER :: U10MAX=25.0

      REAL :: U10, Z

! ----------------------------------------------------------------------

      Z=MIN(U10,U10MAX)
      BETADEF=A+Z*(B+C*Z)

      RETURN
      END FUNCTION BETADEF
