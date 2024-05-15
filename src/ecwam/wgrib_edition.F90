! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTEGER FUNCTION WGRIB_EDITION (IPARAMID)

! ----------------------------------------------------------------------

!**** *WGRIB_EDITION* - FUNCTION TO DETERMINE WHETHER A WAVE PARAMETER WITH IPARAMID
!                       CAN ONLY BE ENCODED IN GRIB 2.
!                       ALL OTHER PARAMETERS WILL USE WHAT IS SPECFIED WITH NGRIB_VERSION

! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWGRIBHD, ONLY : NGRIB_VERSION

! ----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JWIM), INTENT(IN) :: IPARAMID   ! grib parameter id (6 digit) 

! ----------------------------------------------------------------------

!!debile
       write(0,*) 'wgrib_edition ',IPARAMID, NGRIB_VERSION


IF ( IPARAMID >= 140131 .AND. IPARAMID <= 140134 ) THEN
  WGRIB_EDITION = 2
ELSE
  WGRIB_EDITION = NGRIB_VERSION
ENDIF

!!debile
       write(0,*) 'wgrib_edition ',WGRIB_EDITION

END FUNCTION WGRIB_EDITION
