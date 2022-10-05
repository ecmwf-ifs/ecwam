!> \file yowmpimodule
!> \brief just a wrapper for different mpi modules

module yowMpiModule
#ifdef WAM_HAVE_MPI_F08
  use mpi_f08
#else
  use MPL_MPIF
#endif
end module
