!> Has fancy data
module yowDatapool
  use yowMpiModule
  implicit none


#ifdef USE_SINGLE
  !> single precision. Enable with compiler flag -DUSE_SINGLE
  integer,parameter :: rkind = 4
#else
  !> double precision. Default real datatype
  integer,parameter :: rkind = 8
#endif


  logical, parameter :: debugPrePartition = .false.
  logical, parameter :: debugPostPartition = .true.
  logical, parameter :: debugParmetis = .false.
  !> write partition information into file fort.600.
  !> one can display partition with katerfempresenter. just open the system.dat and
  !> click on the partitions icon
  logical, parameter :: debugPartition = .false.

  !> Number of threads
  integer, save :: nTasks = 0

  !> The thread id.
  !> starts by 0. The first Thread has rank 0
  integer, save :: myrank = 0

  !> MPI Communicator.
  !> Should be MPI_COMM_WORLD. If pdlib is run into a existing MPI enviroment, comm is set to a new communicator
  integer,public,save :: comm = MPI_COMM_WORLD

  !> MPI Integer Type.
  !> Should be MPI_INTEGER
  integer,parameter :: itype = MPI_INTEGER

  !> MPI Real Type
  !> Shpuld be MPI_REAL8
#ifdef USE_SINGLE
  integer, parameter :: rtype = MPI_REAL4
#else
  integer, parameter :: rtype = MPI_REAL8  
#endif

end module yowDatapool
