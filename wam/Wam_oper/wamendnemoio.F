      SUBROUTINE WAMENDNEMOIO
!
!**** *WAMENDNEMOIO*  - End
!
!     Purpose.
!     --------
!     Call MPI finalize in case of IO servers.
!
!**   Interface.
!     ----------
!       *CALL*  *WAMENDNEMOIO*
!
!     Input:
!     -----
!
!     Output:
!     ------
!
!     Method:
!     ------
!       NEMOIO usage is controlled by environment variables
!
!     Externals:
!     ---------
!
!     Reference:
!     ---------
!
!     Author:
!     -------
!       K. Mogensen, ECMWF
!
!     Modifications.
!     --------------
!
!     -----------------------------------------------------------

      IMPLICIT NONE

      CHARACTER(LEN=128) :: CLNEMOIO ! String to store NEMOIO env.
      LOGICAL :: LNEMOIO 

#ifdef WITH_NEMO

      CALL GETENV('NEMOIOSERVER',CLNEMOIO)

! Set NEMOIO settings depending of the NEMOIO environment variables
      SELECT CASE (TRIM(CLNEMOIO))
! If NEMOIO is not set don't use NEMOIO.
      CASE ('')
         LNEMOIO    =.FALSE.
! If NEMOIO is "no" don't use NEMOIO.
      CASE ('no')
         LNEMOIO    =.FALSE.
! If NEMOIO is "yes" use NEMOIO for real.
      CASE ('yes') 
         LNEMOIO    =.TRUE.
      CASE DEFAULT
         CALL ABOR1('Invalid value of NEMOIO environment variable '//
     &        'in ininemoio.F90')
      END SELECT

      IF(LNEMOIO) CALL NEMOGCMCOUP_END_IOSERVER

#endif

      END SUBROUTINE WAMENDNEMOIO
