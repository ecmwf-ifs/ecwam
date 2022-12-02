      LOGICAL FUNCTION SET_WFLAGS(FLAG,N)

! ----------------------------------------------------------------------

!     J. BIDLOT     ECMWF 

!*    PURPOSE.
!     --------

!     SETS A GLOBAl FLAG FROM AN ARRAY OF LOGICAL FLAGS (FLAG)
!     DEPENDING ON THE LOGIC REQUIRED BY IFLAG
!     RETURN VALUE WILL BE TRUE, ONLY IF ANY OF THE FLAGS IN THE ARRAY ARE TRUE
!     UP TO SIZE N


!**   INTERFACE.
!     ----------
!     *SET_WFLAGS(FLAG,N)*

!     *FLAG* LOGICAL  INPUT FLAGS 
!     *N*    INTEGER  SIZE OF FLAG

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
 
! ----------------------------------------------------------------------
      IMPLICIT NONE

      !LOGICAL :: SET_WFLAGS

      LOGICAL,            INTENT(IN) :: FLAG(N)
      INTEGER(KIND=JWIM), INTENT(IN) :: N

      INTEGER(KIND=JWIM) :: IFL

      SET_WFLAGS=.FALSE.

      IF (N > 0) THEN
        SET_WFLAGS = FLAG(1)
      ENDIF
      IFL=2
      DO WHILE (.NOT. SET_WFLAGS .AND. IFL <= N)
        SET_WFLAGS = FLAG(IFL)
        IFL=IFL+1
      ENDDO

      END FUNCTION SET_WFLAGS
