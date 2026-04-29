! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
SUBROUTINE WVOPENBATHY (IU06, IU07, KGRIB_HANDLE)  

!****  *WVOPENBATHY* - OPENS THE FILE CONTAINING THE MODEL BATHYMETRY (see FILENAME) AND
!                      DETERMINES WHETHER IT IS A GRIB FILE OR THE OLD BINARY FORMAT

! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWABORT , ONLY : WAM_ABORT
USE YOWGRIB  , ONLY : IGRIB_OPEN_FILE, IGRIB_NEW_FROM_FILE, IGRIB_GET_VALUE, JPGRIB_END_OF_FILE
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
IMPLICIT NONE
#include "iwam_get_unit.intfb.h"

INTEGER(KIND=JWIM), INTENT(IN)  :: IU06  ! standard out unit
INTEGER(KIND=JWIM), INTENT(OUT) :: IU07  ! open unit for input file
INTEGER(KIND=JWIM), INTENT(OUT) :: KGRIB_HANDLE  ! grib handle if input file is in grib

INTEGER(KIND=JWIM):: LFILE, KHANDLE, IRET, IERR, IEDITION
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
CHARACTER(LEN=80) :: FILENAME
LOGICAL :: LLEXIST

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WVOPENBATHY',0,ZHOOK_HANDLE)

FILENAME='wam_grid_tables'

LFILE=0
LLEXIST=.FALSE.
IF (FILENAME /= ' ') LFILE=LEN_TRIM(FILENAME)
INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)

IF (.NOT. LLEXIST .OR. LFILE == 0) THEN
  WRITE(IU06,*) '*****************************************************************'
  WRITE(IU06,*) '*                                                               *'
  WRITE(IU06,*) '*  FATAL ERROR IN SUB. WVOPENBATHY                              *'
  WRITE(IU06,*) '*  =============================                                *'
  WRITE(IU06,*) '*  WAVE MODEL INPUT FILE ',FILENAME(1:LFILE), ' IS MISSING !!!!'
  WRITE(IU06,*) '*                                                               *'
  WRITE(IU06,*) '*****************************************************************'
  CALL WAM_ABORT("WAVE MODEL INPUT FILE '"//FILENAME(1:LFILE)//"' IS MISSING !!!",__FILENAME__,__LINE__)
ENDIF

! Assume first that the input is in grib
IU07=-99
KGRIB_HANDLE=-99

CALL IGRIB_OPEN_FILE(IU07,FILENAME(1:LFILE),'r')

CALL IGRIB_NEW_FROM_FILE(IU07, KHANDLE, IRET)

IF (IRET /= JPGRIB_END_OF_FILE) THEN
  CALL IGRIB_GET_VALUE(KHANDLE,'editionNumber', IEDITION, IERR)
  IF ( IERR == 0 ) THEN
    KGRIB_HANDLE=KHANDLE
    WRITE(IU06,*) ' GRIB INPUT OF MODEL BATHYMETRY  '
    WRITE(IU06,*) ' GRIB HANDLE IU07 =  ', IU07
    WRITE(IU06,*) ' IN GRIB EDITION = ' ,IEDITION
    WRITE(IU06,*) ''
    CALL FLUSH(IU06)
  ENDIF
ENDIF


IF ( KGRIB_HANDLE < 0 ) THEN
  WRITE(IU06,*) ' BINARY INPUT OF MODEL BATHYMETRY  '
  WRITE(IU06,*) ''
  IU07 = IWAM_GET_UNIT(IU06, FILENAME(1:LFILE), 'r', 'u',0,'READ')
    CALL FLUSH(IU06)
ENDIF

IF (LHOOK) CALL DR_HOOK('WVOPENBATHY',1,ZHOOK_HANDLE)

END SUBROUTINE WVOPENBATHY
