      MODULE YOWSARAS

      IMPLICIT NONE

      INTEGER, PARAMETER :: MPARTSW = 15

      REAL, PARAMETER :: PPSARCF = 2.
      REAL, PARAMETER :: PPSDM = 1.
      REAL, PARAMETER :: SARCORDIS = 300000.

      INTEGER :: NSPEC
      INTEGER :: NSPECW
      INTEGER :: JPSIW
      INTEGER :: MAXMPART=0

      INTEGER, ALLOCATABLE :: NW1D(:,:)
      INTEGER, ALLOCATABLE :: NIGW(:)
      INTEGER, ALLOCATABLE :: NIJW(:)
      INTEGER, ALLOCATABLE :: NTOSIW1D(:)
      INTEGER, ALLOCATABLE :: NSIW1D(:,:)
      INTEGER, ALLOCATABLE :: IGSAR(:) 
      INTEGER, ALLOCATABLE :: IJSAR(:) 

      REAL, ALLOCATABLE :: SPEC(:,:,:,:)
      REAL, ALLOCATABLE :: LONG(:,:)
      REAL, ALLOCATABLE :: LAT(:,:)
      REAL, ALLOCATABLE :: COSSARLAT(:)
      REAL, ALLOCATABLE :: SINSARLAT(:)
      REAL, ALLOCATABLE :: U10(:)
      REAL, ALLOCATABLE :: THW(:)

      REAL, ALLOCATABLE :: SPECW(:,:,:)
      REAL, ALLOCATABLE :: DIST(:,:)
      REAL, ALLOCATABLE :: UW10(:)
      REAL, ALLOCATABLE :: THWW(:)


!     VARIABLE   TYPE     PURPOSE                                  
!     --------   ----     -------                                  
!     *MPARTSW*  INTEGER  MAXIMUM NUMBER OF PARTITIONINGS (for arrays).
!     *PPSARCF   REAL     PPSARCF*SARCORDIS IS THE MAXIMUM INFLUENCE
!                         DISTANCE FOR THE SAR ASSIMILATION. 
!     *PPSDM*    REAL     RATIO SIGMA_OBS/SIGMA_FG.
!     *SARCORDIS*REAL     SAR MODEL CORRELATION DISTANCE in m 

!     *NSPEC*    INTEGER  NUMBER OF VALID OBSERVATIONS.
!     *MSPECW*   INTEGER  MAXIMUM NUMBER OF ANALYSED POINTS.
!     *JPSIW*    INTEGER  MAXIMUM  NUMBER OF SAR POINTS INFLUENCING
!                         ONE WAM POINT.
!     *MAXMPART  INTEGER  MAXIMUM NUMBER OF PARTITIONINGS USED AFTER
!                         PARTITIONING OF SAR AND MODEL SPECTRA. 
!     *NW1D*     INTEGER  INDEX IN THE 1D ARRAY OF MODEL POINTS 
!                         INFLUENCED BY THE SAR DATA POINTS FOR EACH
!                         WAM POINT ( referred to by its block indices)
!                         =0 FOR POINTS NOT INFLUENCED BY DATA POINT.
!     *NIGW*     INTEGER  BLOCK INDEX OF POINT IN 1D ARRAY.
!     *NIJW*     INTEGER  INDEX IN A BLOCK OF POINT IN 1D ARRAY.
!     *NTOSIW1D* INTEGER  NUMBER OF SAR DATA POINTS INFLUENCING A WAM
!                         POINT INDEXED IN THE 1D ARRAY.
!     *NSIW1D*   INTEGER  INDEX OF THE SAR SPECTRUM IN THE SAR DATA
!                         ARRAY IFLUENCING THE WAM POINT INDEXED IN THE
!                         1D ARRAY.
!     *IGSAR*    INTEGER  BLOCK NUMBER OF THE CLOSEST MODEL SEA POINT
!                         O will indicate that it is outside the model
!                         sea point grid. 
!     *IJSAR*    INTEGER  BLOCK INDEX OF THE CLOSEST MODEL SEA POINT
!                         O will indicate that it is outside the model
!                         sea point grid. 
!     *SPEC*     REAL     INPUT SAR SPECTRA FOR FIRST GUESS SPECTRUM = 1
!                         SAR RETRIEVED SPECTRUM = 2.
!     *LONG*     REAL     LONGITUDE OF SPEC (IN DEGREE).
!     *LAT*      REAL     LATITUDE OF SPEC (IN DEGREE).
!     *COSSARLAT*REAL     COSINE OF LAT.
!     *SINSARLAT*REAL     SINE OF LAT.
!     *U10*      REAL     FIRST GUESS U10 AT SPEC LOCATION. 
!     *THW*      REAL     FIRST GUESS WIND DIRECTION AT SPEC LOCATION. 
!     *SPECW*    REAL     WAM SPECTRA TO BE ANALYSED.
!     *DIST*     REAL     DISTANCE BETWEEN WAM AND SAR POINT
!                         INFLUENCING WAM POINT.
!     *UW10*     REAL     WIND SPEED AT XLATW,XLONW. 
!     *THWW*     REAL     WIND DIRECTION  AT XLATW,XLONW. 

! ----------------------------------------------------------------------
      END MODULE YOWSARAS
