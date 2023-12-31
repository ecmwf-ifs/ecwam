! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
SUBROUTINE WVOPENSUBBATHY (KGRIB_HANDLE)  

!****  *WVOPENSUBBATHY* - OPENS THE FILE CONTAINING THE MODEL SUBGRID BATHYMETRY INFO (see FILENAME) AND
!                         DETERMINES WHETHER IT IS A GRIB FILE OR THE OLD BINARY FORMAT

! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWABORT , ONLY : WAM_ABORT
USE YOWTEST  , ONLY : IU06
USE YOWSTAT  , ONLY : IPROPAGS
USE YOWUBUF  , ONLY : NPROPAGS
USE YOWUNIT  , ONLY : IU08
USE YOWGRIB  , ONLY : IGRIB_OPEN_FILE, IGRIB_NEW_FROM_FILE, IGRIB_GET_VALUE, JPGRIB_END_OF_FILE
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
IMPLICIT NONE
#include "iwam_get_unit.intfb.h"

INTEGER(KIND=JWIM), INTENT(OUT) :: KGRIB_HANDLE    ! grib handle if input file is in grib

INTEGER(KIND=JWIM):: LFILE, KHANDLE, IRET, IERR, IEDITION
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
CHARACTER(LEN=80) :: FILENAME
CHARACTER(LEN=1)  :: C1
LOGICAL :: LLEXIST

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WVOPENSUBBATHY',0,ZHOOK_HANDLE)

IF (IPROPAGS < 0 .OR. IPROPAGS > NPROPAGS) THEN
  WRITE(IU06,*) '***************************************'
  WRITE(IU06,*) '*                                     *'
  WRITE(IU06,*) '*  FATAL ERROR IN SUB. WVOPENSUBBATHY *'
  WRITE(IU06,*) '*  WRONG VALUE FOR IPROPAGS:          *'
  WRITE(IU06,*) '*  IPROPAGS= ',IPROPAGS
  WRITE(IU06,*) '*                                     *'
  WRITE(IU06,*) '***************************************'
  CALL FLUSH(IU06)
  CALL WAM_ABORT("Wrong value for IPROPAGS",__FILENAME__,__LINE__)
ENDIF

WRITE(C1,'(I1)') IPROPAGS 
FILENAME='wam_subgrid_'//C1

LFILE=0
LLEXIST=.FALSE.
IF (FILENAME /= ' ') LFILE=LEN_TRIM(FILENAME)
INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)

IF (.NOT. LLEXIST) THEN
  WRITE(IU06,*) '*****************************************************************'
  WRITE(IU06,*) '*                                                               *'
  WRITE(IU06,*) '*  FATAL ERROR IN SUB. WVOPENSUBBATHY                           *'
  WRITE(IU06,*) '*  =============================                                *'
  WRITE(IU06,*) '*  WAVE MODEL INPUT FILE ',FILENAME(1:LFILE), ' IS MISSING !!!!'
  WRITE(IU06,*) '*                                                               *'
  WRITE(IU06,*) '*****************************************************************'
  CALL WAM_ABORT("WAVE MODEL INPUT FILE '"//FILENAME(1:LFILE)//"' IS MISSING !!!",__FILENAME__,__LINE__)
ENDIF

! Assume first that the input is in grib
IU08(IPROPAGS)=-99
KGRIB_HANDLE=-99

CALL IGRIB_OPEN_FILE(IU08(IPROPAGS),FILENAME(1:LFILE),'r')

CALL IGRIB_NEW_FROM_FILE(IU08(IPROPAGS), KHANDLE, IRET)

IF (IRET /= JPGRIB_END_OF_FILE) THEN
  CALL IGRIB_GET_VALUE(KHANDLE,'editionNumber', IEDITION, IERR)
  IF ( IERR == 0 ) THEN
    KGRIB_HANDLE=KHANDLE
    WRITE(IU06,*) ' GRIB INPUT OF MODEL SUBGRID BATHYMETRY  '
    WRITE(IU06,*) ' GRIB HANDLE IU08(IPROPAGS) =  ', IU08(IPROPAGS)
    WRITE(IU06,*) ' IN GRIB EDITION = ' ,IEDITION
    WRITE(IU06,*) ''
    CALL FLUSH(IU06)
  ENDIF
ENDIF


IF ( KGRIB_HANDLE < 0 ) THEN
  WRITE(IU06,*) ' BINARY INPUT OF MODEL SUBGRID BATHYMETRY  '
  WRITE(IU06,*) ''
  IU08(IPROPAGS) = IWAM_GET_UNIT(IU06, FILENAME(1:LFILE), 'r', 'u',0,'READWRITE')
  CALL FLUSH(IU06)
ENDIF

IF (LHOOK) CALL DR_HOOK('WVOPENSUBBATHY',1,ZHOOK_HANDLE)

END SUBROUTINE WVOPENSUBBATHY
