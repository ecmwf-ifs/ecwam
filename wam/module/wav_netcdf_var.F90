MODULE WAV_netcdf_var
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
#ifdef NETCDF_OUTPUT_WAM
      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal, MyRankLocal
      implicit none
      INTEGER(KIND=JWIM) :: idxOutput = 0
      INTEGER(KIND=JWIM) :: idxFile = 1
      INTEGER(KIND=JWIM) :: idxInfile = 1
      INTEGER(KIND=JWIM) :: DoOutput
# ifndef MODEL_COUPLING_ATM_WAV
      INTEGER(KIND=JWIM) :: NHIS = 1
      INTEGER(KIND=JWIM) :: NDEFHIS = 6
# endif
      INTEGER(KIND=JWIM) :: NETCDF_nbVar = 11
      REAL(KIND=JWRB), allocatable :: NETCDF_var(:,:)
      REAL(KIND=JWRB), allocatable :: NETCDF_var_gridded(:,:,:)
      logical :: NETCDF_initialized = .FALSE.
      INTEGER(KIND=JWIM), allocatable :: NETCDF_var_rqst(:)
      INTEGER(KIND=JWIM), allocatable :: NETCDF_var_stat(:,:)
      INTEGER(KIND=JWIM), allocatable :: NETCDF_var_type(:)
!      INTEGER(KIND=JWIM), dimension(:), pointer :: NETCDF_rqst
!      INTEGER(KIND=JWIM), dimension(:,:), pointer :: NETCDF_stat
!      INTEGER(KIND=JWIM), dimension(:), pointer :: NETCDF_type
      INTEGER(KIND=JWIM), allocatable :: ListIJS(:), ListIJL(:), ListIJ_OFFSET(:)
      LOGICAL DoNETCDF_sync
      INTEGER(KIND=JWIM) :: MPI_COMM_NETCDF
      REAL(KIND=JWRU) :: DeltaTimeNetcdfWAV = -1
      REAL(KIND=JWRU) :: FileSizeTimeWAV = -1
      INTEGER(KIND=JWIM) :: NDEFHIS_WAV = -1
      REAL(KIND=JWRU) :: WAV_NetcdfPresTime = 0
      INTEGER(KIND=JWIM) :: idxUcurr, idxVcurr, idxZeta
      INTEGER(KIND=JWIM) :: idxcfl1, idxcfl2, idxcfl3
      INTEGER(KIND=JWIM) :: NETCDF_X, NETCDF_Y
      INTEGER(KIND=JWIM) :: MaxMNPloc
      INTEGER(KIND=JWIM), allocatable :: ListMNPloc(:), ListNP_RESloc(:)
      INTEGER(KIND=JWIM), allocatable :: ListIPLGloc(:,:)
#endif
END MODULE WAV_netcdf_var
