      MODULE YOWSARINV

!=======================================================================
!     VARIABLES USED FOR SARINV 
!     To replace the file  sarinv.par  originally in Sar directory.
!=======================================================================

      IMPLICIT NONE

      INTEGER, PARAMETER :: jpmax1 = 3    ! number of maxima selcted
     &                    , jpwam  = 3000 ! maximum number of wam wave spectra per block
     &                    , jpnsar = 3000 ! maximum number of observed SAR spectra

      INTEGER jpksar         ! wave number grid points of observed SAR spectrum
     &      , jptsar         ! direction grid points of observed SAR spectrum
     &      , jpisar         ! number of grid points of a SAR spectrum = jpksar * jptsar
     &      , jpdata         ! the array size of the record containing the observed SAR spectrum
                             ! jpdata = 24 + jpksar + jptsar + jpisar + 1

     &      , jpspec         ! spectral length
     &      , jpimaf         ! length of the real image fields = jpspec - 1
     &      , jpim2          ! = jpimaf * 2
     &      , jphalf         ! = jpimaf / 2 + 1
     &      , jpw            ! = jpimaf * 5 / 2


!     Parameters for 2d-fft cfft2d
!
      INTEGER, PARAMETER :: KSX = 64           ! NUMBER OF WAVE NUMBERS IN SAR SPECTRUM (X-D).
     &                    , KSY = 64           ! NUMBER OF WAVE NUMBERS IN SAR SPECTRUM (Y-D).
     &                    , KSXT = 2*KSX
     &                    , KSYT = 2*KSY
     &                    , LOT = KSYT
     &                    , KSXLOT = 4*KSXT*LOT
     &                    , KSXP1 = KSX+1
     &                    , KSXM1 = KSX-1
     &                    , KSYP1 = KSY+1
     &                    , KSYM1 = KSY-1
     &                    , MGLPR = 10
     &                    , NFRE = 30          ! NUMBER OF FREQUENCIES IN MODEL SPECTRUM.
     &                    , NANG = 24          ! NUMBER OF DIRECT. IN MODEL SPECTRUM.
     &                    , NANG1 = NANG+1
     &                    , NBLO = 20          ! NUMBER OF BLOCKS IN WAVE MODEL GRID.
     &                    , NIBLO = 3000       ! NUMBER OF INDICES IN A BLOCK.
     &                    , NOR = 13

     &                    , ntable = 100+2*(ksxt*2)
     &                    , nwork = 4*ksxt*ksyt
 
      REAL, PARAMETER :: zgrav = 9.806   ! gravitation constant (!)
     &                 , xeps = 1.0E-12  ! small real number (epsilon)
!-----------------------------------------------------------------------

c
c FOR SIRCDATA ONLY
c    -----------------
c
      INTEGER, PARAMETER :: jpksx = 64      ! az. dimension of input SAR spectrum
     &                    , jpksy = 64      ! range dimension of input SAR spectrum

     &                    , maxrlen = 10000 ! maximal record length (for sircdata) in bytes

      REAL, PARAMETER :: pixelx = 31.4   ! azimuthal pixel size
     &                 , pixely = 35.7   ! range pixel size


! ----------------------------------------------------------------------
      END MODULE YOWSARINV
