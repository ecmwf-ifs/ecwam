      SUBROUTINE INIWCST 

! ----------------------------------------------------------------------

!**** *INIWCST* -

!*    PURPOSE.
!     --------
!     INITIALISES CONSTANTS

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
     &            RAD      ,DEG      ,R

! ----------------------------------------------------------------------

      PI = 4.*ATAN(1.)
      ZPI= 2.*PI
      ZCONST=1./(8.*PI*SQRT(2.))
      RAD = PI/180.
      DEG = 180./PI
      R = CIRC/ZPI

      RETURN
      END SUBROUTINE INIWCST
