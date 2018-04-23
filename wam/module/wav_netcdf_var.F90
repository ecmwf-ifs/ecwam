MODULE WAV_netcdf_var
#ifdef NETCDF_OUTPUT_WAM
      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal, MyRankLocal
# if defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
      USE coupling_var, only : WAV_COMM_WORLD
# endif
      implicit none
      integer :: idxOutput = 0
      integer :: idxFile = 1
      integer :: idxInfile = 1
      integer :: DoOutput
# ifndef MODEL_COUPLING_ATM_WAV
      integer :: NHIS = 1
      integer :: NDEFHIS = 6
# endif
      integer :: NETCDF_nbVar = 11
      real, allocatable :: NETCDF_var(:,:)
      real, allocatable :: NETCDF_var_gridded(:,:,:)
      logical :: NETCDF_initialized = .FALSE.
      integer, allocatable :: NETCDF_var_rqst(:)
      integer, allocatable :: NETCDF_var_stat(:,:)
      integer, allocatable :: NETCDF_var_type(:)
!      integer, dimension(:), pointer :: NETCDF_rqst
!      integer, dimension(:,:), pointer :: NETCDF_stat
!      integer, dimension(:), pointer :: NETCDF_type
      integer, allocatable :: ListIJS(:), ListIJL(:), ListIJ_OFFSET(:)
      LOGICAL DoNETCDF_sync
      integer MPI_COMM_NETCDF
      real*8 :: DeltaTimeNetcdfWAV = -1
      real*8 :: FileSizeTimeWAV = -1
      integer :: NDEFHIS_WAV = -1
      real*8 :: WAV_NetcdfPresTime = 0
      integer idxUcurr, idxVcurr, idxZeta
      integer idxcfl1, idxcfl2, idxcfl3
      integer NETCDF_X, NETCDF_Y
      integer MaxMNPloc
      integer, allocatable :: ListMNPloc(:), ListNP_RESloc(:)
      integer, allocatable :: ListIPLGloc(:,:)
#endif
END MODULE WAV_netcdf_var
