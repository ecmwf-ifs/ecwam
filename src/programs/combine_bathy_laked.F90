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

INTEGER(KIND=JWIM), PARAMETER :: NPARAM=3 
INTEGER(KIND=JWIM) :: IP, LFILE
INTEGER(KIND=JWIM) :: IU06, IU07, KGRIB_HANDLE_BATHY
INTEGER(KIND=JWIM), DIMENSION(NPARAM) :: IULAKE, KGRIB_HANDLE_LAKE, IPARAMID
INTEGER(KIND=JWIM) :: NUMBEROFVALUES 

REAL(KIND=JWRB), ALLOCATABLE :: VALUES(:)

CHARACTER(LEN=80) :: FILENAME
CHARACTER(LEN=80), DIMENSION(NPARAM) :: INFILENAME

LOGICAL :: LLEXIST

LOGICAL :: LLSCANNS

! ----------------------------------------------------------------------

IU06 = 6 

IU07 = -1
KGRIB_HANDLE_BATHY = -99

IULAKE(:) = -1
KGRIB_HANDLE_LAKE(:) = -99

!  INPUT FILE:
!  -----------

! LAND SEA MASK
INFILENAME(1)='lsm'
IPARAMID(1)=172
! LAKE COVER
INFILENAME(2)='clake'
IPARAMID(2)=128026
! LAKE DEPTH
INFILENAME(3)='lakedl'
IPARAMID(3)=228007

! BATHY:
CALL WVOPENBATHY (IU06, IU07, KGRIB_HANDLE_BATHY)

DO IP = 1, NPARAM
  FILENAME = INFILENAME(IP)
  LFILE=0
  LLEXIST=.FALSE.
  IF (FILENAME /= ' ') LFILE=LEN_TRIM(FILENAME)
  INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)

  IF (.NOT. LLEXIST) THEN
    WRITE(IU06,*) '*****************************************************************'
    WRITE(IU06,*) '*                                                               *'
    WRITE(IU06,*) '*  FATAL ERROR IN                               *'
    WRITE(IU06,*) '*  =============================                                *'
    WRITE(IU06,*) '*  WAVE MODEL INPUT FILE ',FILENAME(1:LFILE), ' IS MISSING !!!!'
    WRITE(IU06,*) '*                                                               *'
    WRITE(IU06,*) '*****************************************************************'
    CALL WAM_ABORT("INPUT FILE '"//FILENAME(1:LFILE)//"' IS MISSING !!!",__FILENAME__,__LINE__)
  ENDIF

ENDDO


IF ( KGRIB_HANDLE_BATHY > 0 ) THEN
  !! GRIB INPUT:
!    ----------

!    GRID INFO:
!    ---------

  CALL WVGETGRIDINFO(IU06, KGRIB_HANDLE_BATHY, &
 &                   NGX, NGY, IPER, IRGG, IQGAUSS, NLONRGG, LLSCANNS, &
 &                   AMOWEP, AMOSOP, AMOEAP, AMONOP, XDELLA, XDELLO )





  CALL IGRIB_CLOSE_FILE(IU07)
  CALL IGRIB_RELEASE(KGRIB_HANDLE_BATHY)

ELSE

  CALL WAM_ABORT("Adding Lake information to binary bathymetry not available",__FILENAME__,__LINE__)

ENDIF ! GRIB OR BINARY INPUT


END PROGRAM 
