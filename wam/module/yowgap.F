      MODULE YOWGAP

      IMPLICIT NONE

      REAL, PARAMETER :: GAPVAL=0.1E-19

      INTEGER :: MGAP

      INTEGER, ALLOCATABLE :: LANG(:,:) 
      INTEGER, ALLOCATABLE :: HANG(:,:) 
      INTEGER, ALLOCATABLE :: LFRE(:,:) 
      INTEGER, ALLOCATABLE :: HFRE(:,:) 

      LOGICAL, ALLOCATABLE :: EQUIGAP(:,:,:)

!     VARIABLE   TYPE     PURPOSE                                  
!     --------   ----     -------                                  
!     *GAPVAL*   REAL     VALUE GIVEN TO SPECTRAL TO TAG THEM AS GAPS
!     *MGAP*     INTEGER  MAXIMUM NUMBER OF GAPS IN THE TRANSFORMED
!                         SPECTRA.
!     *LANG*     INTEGER  LOWER LEFT AND UPPER RIGHT CORNER OF A GAP.  
!     *HANG*     INTEGER  LOWER LEFT AND UPPER RIGHT CORNER OF A GAP.
!     *LFRE*     INTEGER  LOWER LEFT AND UPPER RIGHT CORNER OF A GAP.
!     *HFRE*     INTEGER  LOWER LEFT AND UPPER RIGHT CORNER OF A GAP.
!     *EQUIGAP*  LOGICAL  FLAG INDICATING WHETHER GAPS ARE CONTAINED
!                         IN ONE ANOTHER OR OVERLAPPING.
!                         EQUIVALENCE OF GAPS, WHERE
!                         EQUIGAP(ISPEC,IGAP,GAP) = .TRUE. MEANS,
!                         THAT IN SPECTRUM NUMBER ISPEC THE GAP
!                        x  NUMBERS IGAP AND JGAP ARE EQUIVALENT,
!                         ONLY USED IN SUBROUTINE MAKEFRAMES

! ----------------------------------------------------------------------
      END MODULE YOWGAP
