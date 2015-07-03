      SUBROUTINE WAMININEMOIO(LNEMOIO,LNEMOIOSERVER)
!     
!**** *WAMININEMOIO*  - Initialize NEMO IO server (if present)
!     
!     Purpose.
!     --------
!     Get a IFS MPI communicator to use when NEMO IO server is active.
!
!**   Interface.
!     ----------
!       *CALL*  *WAMININEMOIO*
!
!     Input:
!     -----
!
!     Output:
!     ------
!       MPLUSERCOMM in MPL_MODULE is updated.
!
!     Method:
!     ------
!       NEMOIO usage is controlled by environment variables
!
!     Externals:
!     ---------
!       GETENV - Get enviroment variables
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

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE MPL_MODULE
      USE MPL_MPIF
      
      IMPLICIT NONE

      LOGICAL, INTENT(OUT) :: LNEMOIOSERVER
      LOGICAL, INTENT(INOUT) :: LNEMOIO 
      CHARACTER(LEN=128) :: CLNEMOIO              ! String to store NEMOIO env.
      INTEGER(KIND=JWIM) :: ILOCAL_COMM ! Local communicator 
      INTEGER(KIND=JWIM) :: IERR
      
#ifdef WITH_NEMO

!     Get the enviroment variables

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
         CALL ABOR1('Invalid value of NEMOIO environment variable'//
     &              'in ininemoio.F90')
      END SELECT

      IF (LNEMOIO) THEN
         CALL NEMOGCMCOUP_INIT_IOSERVER( ILOCAL_COMM, LNEMOIOSERVER )
         CALL NEMOGCMCOUP_INIT_IOSERVER_2( MPI_COMM_WORLD )
         LMPLUSERCOMM=.TRUE.
         MPLUSERCOMM=ILOCAL_COMM
      ELSE
         LNEMOIOSERVER=.FALSE.
      ENDIF
      
#else

      LNEMOIOSERVER=.FALSE.
      LNEMOIO=.FALSE.

#endif

      END
