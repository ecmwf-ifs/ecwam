      MODULE YOWGRIBHD

      IMPLICIT NONE

!*    ** *GRIB HEADERS:  COMMON (DEFAULT) GRIB HEADERS.

      INTEGER :: IMDLGRBID_G=107 !! see below the rule on how to select IMDLGRBID_G
      INTEGER :: IMDLGRBID_M=207

      INTEGER, PARAMETER :: NKSEK0=2
      INTEGER, PARAMETER :: NKSEK1=2048
      INTEGER, PARAMETER :: NKSEK2=2070
      INTEGER, PARAMETER :: NKSEK3=2
      INTEGER, PARAMETER :: NKSEK4=252

      INTEGER, PARAMETER :: NPSEK2=96
      INTEGER, PARAMETER :: NPSEK3=2

      INTEGER,ALLOCATABLE :: KSEK0(:)
      INTEGER,ALLOCATABLE :: KSEK1(:)
      INTEGER,ALLOCATABLE :: KSEK2(:)
      INTEGER,ALLOCATABLE :: KSEK3(:)
      INTEGER,ALLOCATABLE :: KSEK4(:)
 
      INTEGER :: IFS_KSEC1(NKSEK1)

      REAL, ALLOCATABLE :: PSEK2(:)
      REAL, ALLOCATABLE :: PSEK3(:)

      REAL, PARAMETER :: PPMISS=-3
      REAL, PARAMETER :: PPEPS=1.0E-10
      REAL, PARAMETER :: PPREC=0.
      REAL, PARAMETER :: PPRESOL=0.002

      INTEGER NTENCODE

      INTEGER, PARAMETER :: NGRBRESI=16
      INTEGER, PARAMETER :: NGRBRESS=9

      CHARACTER*1, PARAMETER :: HOPERI='M' 
      CHARACTER*1, PARAMETER :: HOPERS='C' 

      LOGICAL :: LGRHDIFS
      LOGICAL :: LNEWLVTP



!*    VARIABLE.   TYPE.     PURPOSE.
!     ---------   -------   --------
!     IMDLGRBID_G INTEGER   GLOBAL MODEL IDENTIFICATION FOR GRIB CODING
!                           IT CAN ALSO BE MODIFIED IN THE IMPUT NAMELIST.
!     IMDLGRBID_M INTEGER   LAW MODEL IDENTIFICATION FOR GRIB CODING
!                           IT CAN ALSO BE MODIFIED IN THE IMPUT NAMELIST.
! 
! The generating process identification numbers (model numbers) for GRIB
! headers are allocated in pre-defined ranges for ECMWF GRIB coded
! fields. The field in the GRIB code for this number is 1 octet and with
! the value 255 indicating 'missing value', only numbers 1-254 are
! available for use.

! The atmospheric model allocated numbers are in the range 121 to 203.
! The Global wave model allocated numbers are in the range 104 to 120.
! The LAW model allocated numbers are in the range 204 to 220.
 
! The numbers 221 to 254 stay reserved for the moment.

! This pre-allocation was introduced to enable some Member States, which
! use the model number to identify products, to write their software in
! such a way that model number changes did not cause them problems with
! hard-coded model numbers. 

!     NKSEK0      INTEGER   DIMENSION OF KSEK0
!     NKSEK1      INTEGER   DIMENSION OF KSEK1
!     NKSEK2      INTEGER   DIMENSION OF KSEK2
!     NKSEK3      INTEGER   DIMENSION OF KSEK3
!     NKSEK4      INTEGER   DIMENSION OF KSEK4
!     NPSEK2      INTEGER   DIMENSION OF PSEK2
!     NPSEK3      INTEGER   DIMENSION OF PSEK3
!     KSEK0       INTEGER   DEFAULT VALUES  FOR KSEC0
!     KSEK1       INTEGER   DEFAULT VALUES  FOR KSEC1
!     KSEK2       INTEGER   DEFAULT VALUES  FOR KSEC2
!     KSEK3       INTEGER   DEFAULT VALUES  FOR KSEC3
!     KSEK4       INTEGER   DEFAULT VALUES  FOR KSEC4
!     IFS_KSEC1   INTEGER   KSEC1 FOR IFS INPUT TO WAM. 
!     PSEK2       REAL      DEFAULT VALUES  FOR PSEC2
!     PSEK3       REAL      DEFAULT VALUES  FOR PSEC3
!     PPMISS      REAL      ALL SPECTRAL VALUES LESS OR EQUAL PPMISS ARE
!                           REPLACED BY THE MISSING DATA INDICATOR
!     PPEPS       REAL      SMALL NUMBER USED IN SPECTRAL PACKING OF 251
!     PPREC       REAL      REFERENCE VALUE FOR SPECTRAL PACKING OF 251
!     PPRESOL     REAL      MAXIMUN RESOLUTION POSSIBLE WHEN ENCODING 
!                           SPECTRA (PARAMETER 251).
!     NTENCODE    INTEGER   TOTAL NUMBER OF GRID POINTS FOR ENCODING
!     NGRBRESI    INTEGER   NUMBER OF BITS USED TO ENCODE INTEGRATED
!                           PARAMETERS
!     NGRBRESS    INTEGER   NUMBER OF BITS USED TO ENCODE SPECTRA
!     HOPERI      CHARACTER GRIB ENCODING ACTION FOR INTEGRATED FIELDS.
!     HOPERS      CHARACTER GRIB ENCODING ACTION FOR SPECTRA. 
!     LGRHDIFS    LOGICAL   IF TRUE THAN ISEC1 WILL USE INFORMATION 
!                           AS PROVIDED BY THE IFS.
!     LNEWLVTP    LOGICAL   IF TRUE THE NEW LEVTYPE DEFINITION WILL BE USED.

 
! ----------------------------------------------------------------------
      END MODULE YOWGRIBHD
