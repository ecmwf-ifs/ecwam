!> \file yowmpimodule
!> \brief just a wrapper for different mpi modules

module yowMpiModule
#ifdef ECMWF
  use MPL_MPIF
#else
! JB to fix ??????
!  use MPI
#endif
end module
