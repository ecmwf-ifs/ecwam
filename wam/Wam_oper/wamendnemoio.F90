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

#ifdef WITH_NEMO

      CALL NEMOGCMCOUP_END_IOSERVER

#endif

      END SUBROUTINE WAMENDNEMOIO
