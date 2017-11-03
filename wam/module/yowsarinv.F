      MODULE YOWSARINV

!=======================================================================
!     VARIABLES USED FOR SARINV 
!     To replace the file  sarinv.par  originally in Sar directory.
!=======================================================================

      IMPLICIT NONE

      INTEGER, PARAMETER :: jpmax1 = 3    ! number of maxima selcted
      INTEGER, PARAMETER :: jpwam  = 3000 ! maximum number of wam wave spectra per block
      INTEGER, PARAMETER :: jpnsar = 3000 ! maximum number of observed SAR spectra

      INTEGER jpksar         ! wave number grid points of observed SAR spectrum
      INTEGER jptsar         ! direction grid points of observed SAR spectrum
      INTEGER jpisar         ! number of grid points of a SAR spectrum = jpksar * jptsar
      INTEGER jpdata         ! the array size of the record containing the observed SAR spectrum
                             ! jpdata = 24 + jpksar + jptsar + jpisar + 1

      INTEGER jpspec         ! spectral length
      INTEGER jpimaf         ! length of the real image fields = jpspec - 1
      INTEGER jpim2          ! = jpimaf * 2
      INTEGER jphalf         ! = jpimaf / 2 + 1
      INTEGER jpw            ! = jpimaf * 5 / 2


!     Parameters for 2d-fft cfft2d
!
      INTEGER, PARAMETER :: KSX = 64           ! NUMBER OF WAVE NUMBERS IN SAR SPECTRUM (X-D).
      INTEGER, PARAMETER :: KSY = 64           ! NUMBER OF WAVE NUMBERS IN SAR SPECTRUM (Y-D).
      INTEGER, PARAMETER :: KSXT = 2*KSX
      INTEGER, PARAMETER :: KSYT = 2*KSY
      INTEGER, PARAMETER :: LOT = KSYT
      INTEGER, PARAMETER :: KSXLOT = 4*KSXT*LOT
      INTEGER, PARAMETER :: KSXP1 = KSX+1
      INTEGER, PARAMETER :: KSXM1 = KSX-1
      INTEGER, PARAMETER :: KSYP1 = KSY+1
      INTEGER, PARAMETER :: KSYM1 = KSY-1
      INTEGER, PARAMETER :: MGLPR = 10
      INTEGER, PARAMETER :: NFRE = 30          ! NUMBER OF FREQUENCIES IN MODEL SPECTRUM.
      INTEGER, PARAMETER :: NANG = 24          ! NUMBER OF DIRECT. IN MODEL SPECTRUM.
      INTEGER, PARAMETER :: NANG1 = NANG+1
      INTEGER, PARAMETER :: NBLO = 20          ! NUMBER OF BLOCKS IN WAVE MODEL GRID.
      INTEGER, PARAMETER :: NIBLO = 3000       ! NUMBER OF INDICES IN A BLOCK.
      INTEGER, PARAMETER :: NOR = 13

      INTEGER, PARAMETER :: ntable = 100+2*(ksxt*2)
      INTEGER, PARAMETER :: nwork = 4*ksxt*ksyt
 
      REAL, PARAMETER :: zgrav = 9.806   ! gravitation constant (!)
      REAL, PARAMETER :: xeps = 1.0E-12  ! small real number (epsilon)
!-----------------------------------------------------------------------

!
! FOR SIRCDATA ONLY
! -----------------
!
      INTEGER, PARAMETER :: jpksx = 64      ! az. dimension of input SAR spectrum
      INTEGER, PARAMETER :: jpksy = 64      ! range dimension of input SAR spectrum

      INTEGER, PARAMETER :: maxrlen = 10000 ! maximal record length (for sircdata) in bytes

      REAL, PARAMETER :: pixelx = 31.4   ! azimuthal pixel size
      REAL, PARAMETER :: pixely = 35.7   ! range pixel size


! ----------------------------------------------------------------------
      END MODULE YOWSARINV
