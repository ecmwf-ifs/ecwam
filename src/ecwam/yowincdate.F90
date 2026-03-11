! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

! ======================================================================

MODULE YOWINCDATE

! ======================================================================

!****     *YOWINCDATE* - MODULE FOR DATE TIME GROUP INCREMENT

!     J. BIDLOT   FEB 2007    RE-INRODUCING THE OLD WAY WITHOUT 
!                             THE NEED FOR ECLIB. ADDING SECONDS
!     REFACTORED  FEB 2026    GENERIC INTERFACE FOR 32 AND 64 BIT

!**   PURPOSE.
!     --------
!       PROVIDES GENERIC INTERFACE FOR UPDATING DATE TIME GROUP
!       WITH BOTH 32-BIT AND 64-BIT INTEGER SHIFTS.

!**   INTERFACE.
!     ----------

!       USE YOWINCDATE, ONLY: INCDATE
!       *CALL* *INCDATE (CDATE,ISHIFT)*
!         *CDATE*  CHAR*12 - DATE TIME GROUP (YYYYMMDDHHMM) OR
!                  CHAR*14 - DATE TIME GROUP (YYYYMMDDHHMMSS)
!         *ISHIFT* INTEGER(JWIM) or INTEGER(JWIB) - TIME INCREMENT IN SECONDS

!     METHOD.
!     -------

!       GENERIC INTERFACE DISPATCHES TO APPROPRIATE TYPED VERSION.
!       THE DATE AND TIME CAN BE SUPPLIED IN THE Y2K COMPLIANT
!       14 CHARACTER OR 12 CHARACTER FORMAT. IF THE 10 CHARACTER
!       FORMAT IS USED, AN ERROR MESSAGE WILL BE ISSUED. 

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWIB, JWRB, JWRU

      IMPLICIT NONE
      PRIVATE

      PUBLIC :: INCDATE

      INTERFACE INCDATE
        MODULE PROCEDURE INCDATE_JWIM, INCDATE_JWIB
      END INTERFACE INCDATE

      CONTAINS

! ======================================================================

      SUBROUTINE INCDATE_JWIM (CDATE, ISHIFT)

! ======================================================================

!****     *INCDATE_JWIM* - UPDATE DATE TIME GROUP (32-BIT VERSION)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(INOUT) :: CDATE
      INTEGER(KIND=JWIM), INTENT(IN) :: ISHIFT

      INTEGER(KIND=JWIM) :: IL, IRET
      INTEGER(KIND=JWIM) :: IYEAR, IMON, IDAY, IHOUR, IMIN, ISEC
      INTEGER(KIND=JWIM) :: MON(12)

      CHARACTER(LEN=80) :: CLFMT

      LOGICAL :: LLND

      DATA MON /31,28,31,30,31,30,31,31,30,31,30,31/

! ----------------------------------------------------------------------

#define ISHIFT_KIND JWIM
#include "incdate_template.inc"
#undef ISHIFT_KIND

      RETURN

      END SUBROUTINE INCDATE_JWIM

! ======================================================================

      SUBROUTINE INCDATE_JWIB (CDATE, ISHIFT)

! ======================================================================

!****     *INCDATE_JWIB* - UPDATE DATE TIME GROUP (64-BIT VERSION)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(INOUT) :: CDATE
      INTEGER(KIND=JWIB), INTENT(IN) :: ISHIFT

      INTEGER(KIND=JWIM) :: IL, IRET
      INTEGER(KIND=JWIM) :: IYEAR, IMON, IDAY, IHOUR, IMIN, ISEC
      INTEGER(KIND=JWIM) :: MON(12)

      CHARACTER(LEN=80) :: CLFMT

      LOGICAL :: LLND

      DATA MON /31,28,31,30,31,30,31,31,30,31,30,31/

! ----------------------------------------------------------------------

#define ISHIFT_KIND JWIB
#include "incdate_template.inc"
#undef ISHIFT_KIND

      RETURN

      END SUBROUTINE INCDATE_JWIB

! ======================================================================

      INTEGER(KIND=JWIM) FUNCTION MFEB_LENGTH(IYEAR)

! ======================================================================

!     LENGTH OF FEBRUARY IN DAYS FOR YEAR IYEAR (YYYY)

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IYEAR

      IF (MOD(IYEAR,400) == 0) THEN
        MFEB_LENGTH=29
      ELSEIF (MOD(IYEAR,100) == 0) THEN
        MFEB_LENGTH=28
      ELSEIF (MOD(IYEAR,4) == 0) THEN
        MFEB_LENGTH=29
      ELSE
        MFEB_LENGTH=28
      ENDIF

      END FUNCTION MFEB_LENGTH

END MODULE YOWINCDATE
