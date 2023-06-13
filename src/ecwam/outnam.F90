! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE OUTNAM (KNANG, KNFRE, KNGX, KNGY) 

! ----------------------------------------------------------------------

!*       PARAMETER NAMELIST AS IN COMMON PARAM OF THE WAM-MODEL.
!        -------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KNANG, KNFRE, KNGX, KNGY

      INTEGER(KIND=JWIM) :: iang, ifre, igx, igy


!!! it is used by Altimeter software !!!
      NAMELIST /PARWAM/ iang, ifre, igx, igy


!***  2. DETERMINE PARAMETERS NAMELIST.
!     ---------------------------------

      iang    = KNANG
      ifre    = KNFRE
      igx     = KNGX
      igy     = KNGY

!***  3. WRITE TO FILE.
!     -----------------

      OPEN(UNIT=91, FILE='PARWAM')
      WRITE (91,PARWAM)
      CLOSE (UNIT=91)

      END SUBROUTINE OUTNAM
