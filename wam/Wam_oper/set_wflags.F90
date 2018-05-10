      FUNCTION SET_WFLAGS(FLAG,N)

! ----------------------------------------------------------------------

!     J. BIDLOT     ECMWF 

!*    PURPOSE.
!     --------

!     SETS A GLOBAl FLAG FROM AN ARRAY OF LOGICAL FLAGS (FLAG)
!     DEPENDING ON THE LOGIC REQUIRED BY IFLAG


!**   INTERFACE.
!     ----------
!     *SET_WFLAGS(FLAG,N,IFLAG)*

!     *FLAG* LOGICAL  INPUT FLAGS 
!     *N*    INTEGER  SIZE OF FLAG

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
 
! ----------------------------------------------------------------------
      IMPLICIT NONE

      LOGICAL :: SET_WFLAGS

      INTEGER(KIND=JWIM) :: N, IFLAG, IFL, NTRAIN
      LOGICAL :: FLAG(N)

      SET_WFLAGS=.FALSE.

      SET_WFLAGS = FLAG(1)
      IFL=2
      DO WHILE (.NOT. SET_WFLAGS .AND. IFL.LE.N)
        SET_WFLAGS = FLAG(IFL)
        IFL=IFL+1
      ENDDO

      END FUNCTION SET_WFLAGS
