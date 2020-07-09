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

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS , ONLY : PI       ,CIRC     ,ZPI      ,ZCONST   ,    &
     &            RAD      ,DEG      ,R        ,ZPISQRT  ,ZPI4GM2  ,    &
     &            G        ,THREEZPI

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(IN) :: PRPLRADI

      PI = 4.0_JWRB*ATAN(1.0_JWRB)
      ZPI= 2.0_JWRB*PI
      THREEZPI = 3.0_JWRB*ZPI
      ZPISQRT = SQRT(ZPI)
      ZPI4GM2 = ZPI**4/G**2

      ZCONST=1.0_JWRB/(8.0_JWRB*PI*SQRT(2.0_JWRB))
      RAD = PI/180.0_JWRB
      DEG = 180./PI
      R = CIRC/ZPI*PRPLRADI

      END SUBROUTINE INIWCST
