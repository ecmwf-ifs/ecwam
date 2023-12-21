! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

PROGRAM combine_bathy_laked

! ----------------------------------------------------------------------

!**** *COMBINE_BATHY_LAKED* 

!*    PURPOSE.
!     --------

!     TO COMBINE BATHY WITH LAKE DEPTH DATA

! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWABORT , ONLY : WAM_ABORT
USE YOWGRIBINFO, ONLY : WVGETGRIDINFO
USE YOWMAP   , ONLY : NGX, NGY, IPER, IRGG, IQGAUSS,                   &
                   &  AMOWEP, AMOSOP, AMOEAP, AMONOP, XDELLA, XDELLO,  &
                   &  NLONRGG
USE YOWPCONS , ONLY : ZMISS
USE YOWGRIB  , ONLY : IGRIB_GET_VALUE, IGRIB_CLOSE_FILE, IGRIB_RELEASE, &
                    & IGRIB_SET_VALUE

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "wvopenbathy.intfb.h"

INTEGER(KIND=JWIM) :: IU06, IU07, KGRIB_HANDLE
INTEGER(KIND=JWIM) :: NUMBEROFVALUES 

REAL(KIND=JWRB), ALLOCATABLE :: VALUES(:)

LOGICAL :: LLSCANNS

! ----------------------------------------------------------------------



IU06 = 6 
IU07 = -1
KGRIB_HANDLE = -99


!     READ INPUT BATHYMETRY:
!     ----------------------

CALL WVOPENBATHY (IU06, IU07, KGRIB_HANDLE)

IF ( KGRIB_HANDLE > 0 ) THEN
  !! GRIB INPUT:
!    ----------

!    GRID INFO:
!    ---------

  CALL WVGETGRIDINFO(IU06, KGRIB_HANDLE, &
 &                   NGX, NGY, IPER, IRGG, IQGAUSS, KLONRGG, LLSCANNS, &
 &                   AMOWEP, AMOSOP, AMOEAP, AMONOP, XDELLA, XDELLO )


  CALL IGRIB_CLOSE_FILE(IU07)
  CALL IGRIB_RELEASE(KGRIB_HANDLE)

ELSE

  CALL WAM_ABORT("Adding Lake information to binary bathymetry not available",__FILENAME__,__LINE__)

ENDIF ! GRIB OR BINARY INPUT


END PROGRAM 
