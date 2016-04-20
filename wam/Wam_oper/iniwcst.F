      SUBROUTINE INIWCST(PRPLRADI) 

! ----------------------------------------------------------------------

!**** *INIWCST* -

!*    PURPOSE.
!     --------
!     INITIALISES CONSTANTS

!      *RPLRADI*   REAL      EARTH RADIUS REDUCTION FACTOR FOR SMALL PLANET

!**   INTERFACE.
!     ----------
!     *CALL INIWCST*

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------

      USE YOWPCONS , ONLY : PI       ,CIRC     ,ZPI      ,ZCONST   ,
     &            RAD      ,DEG      ,R        ,ZPISQRT  ,G2ZPI4   ,
     &            G

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL, INTENT(IN) :: PRPLRADI

      PI = 4.*ATAN(1.)
      ZPI= 2.*PI
      ZPISQRT = SQRT(PI)
      G2ZPI4 = G**2/ZPI**4

      ZCONST=1./(8.*PI*SQRT(2.))
      RAD = PI/180.
      DEG = 180./PI
      R = CIRC/ZPI*PRPLRADI

      RETURN
      END SUBROUTINE INIWCST
