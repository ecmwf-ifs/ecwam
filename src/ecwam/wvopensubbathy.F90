! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
SUBROUTINE WVOPENSUBBATHY (IREAD, NPR, FILENAME, KFILE_HANDLE, KGRIB_HANDLE)

!****  *WVOPENSUBBATHY* - OPENS THE FILE CONTAINING THE MODEL SUBGRID BATHYMETRY INFO (see FILENAME)
!                         AND DETERMINES WHETHER IT IS A GRIB FILE OR THE OLD BINARY FORMAT

! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWABORT  , ONLY : WAM_ABORT
USE YOWTEST   , ONLY : IU06
USE YOWSTAT   , ONLY : IPROPAGS
USE YOWMPP    , ONLY : IRANK
USE YOWUBUF   , ONLY : NPROPAGS
USE YOWUNIT   , ONLY : IU08
USE YOWGRIB   , ONLY : IGRIB_OPEN_FILE, IGRIB_NEW_FROM_FILE, IGRIB_GET_VALUE, JPGRIB_END_OF_FILE
USE EC_LUN    , ONLY : NULERR
USE MPL_MODULE, ONLY : MPL_BARRIER, MPL_BROADCAST
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
IMPLICIT NONE
#include "iwam_get_unit.intfb.h"

INTEGER(KIND=JWIM), INTENT(IN) :: IREAD            ! task accessing the data file
INTEGER(KIND=JWIM), INTENT(IN) :: NPR              ! total number of tasks 
CHARACTER(LEN=72),  INTENT(OUT) :: FILENAME      ! input file name associated with the grib handels 
INTEGER(KIND=JWIM), INTENT(OUT) :: KFILE_HANDLE    ! grib file handle if input file is in grib
INTEGER(KIND=JWIM), INTENT(OUT) :: KGRIB_HANDLE    ! grib handle if input file is in grib

INTEGER(KIND=JWIM):: LFILE, KHANDLE, IRET, IERR, IEDITION
INTEGER(KIND=JWIM) :: IBUF(3)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
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

IF (IRANK == IREAD) THEN
  LFILE=0
  LLEXIST=.FALSE.
  IF (FILENAME /= ' ') LFILE=LEN_TRIM(FILENAME)
  INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)

  IF (.NOT. LLEXIST .OR. LFILE == 0) THEN
    WRITE(IU06,*) '*****************************************************************'
    WRITE(IU06,*) '*                                                               *'
    WRITE(IU06,*) '*  FATAL ERROR IN SUB. WVOPENSUBBATHY                           *'
    WRITE(NULERR,*) '*  FATAL ERROR IN SUB. WVOPENSUBBATHY                           *'
    WRITE(IU06,*) '*  =============================                                *'
    WRITE(IU06,*) '*  WAVE MODEL INPUT FILE ',FILENAME(1:LFILE), ' IS MISSING !!!!'
    WRITE(NULERR,*) '*  WAVE MODEL INPUT FILE ',FILENAME(1:LFILE), ' IS MISSING !!!!'
    WRITE(IU06,*) '*                                                               *'
    WRITE(IU06,*) '*****************************************************************'
    CALL WAM_ABORT("WAVE MODEL INPUT FILE '"//FILENAME(1:LFILE)//"' IS MISSING !!!",__FILENAME__,__LINE__)
  ENDIF

  ! Assume first that the input is in grib
  KFILE_HANDLE = -99
  KGRIB_HANDLE = -99
  IU08(IPROPAGS) = -99

  CALL IGRIB_OPEN_FILE(KFILE_HANDLE,FILENAME(1:LFILE),'r')

  CALL IGRIB_NEW_FROM_FILE(KFILE_HANDLE, KHANDLE, IRET)

  IF (IRET /= JPGRIB_END_OF_FILE) THEN
    CALL IGRIB_GET_VALUE(KHANDLE,'editionNumber', IEDITION, IERR)
    IF ( IERR == 0 ) THEN
      KGRIB_HANDLE=KHANDLE
      WRITE(IU06,*) ''
      WRITE(IU06,*) ' GRIB INPUT OF MODEL SUBGRID BATHYMETRY  '
      WRITE(IU06,*) ' GRIB HANDLE KFILE_HANDLE =  ', KFILE_HANDLE
      WRITE(IU06,*) ' IN GRIB EDITION = ' ,IEDITION
      CALL FLUSH(IU06)
    ENDIF
  ENDIF


  IF ( KGRIB_HANDLE < 0 ) THEN
    WRITE(IU06,*) ' BINARY INPUT OF MODEL SUBGRID BATHYMETRY  '
    WRITE(IU06,*) ''
    IU08(IPROPAGS) = IWAM_GET_UNIT(IU06, FILENAME(1:LFILE), 'r', 'u',0,'READ')
    CALL FLUSH(IU06)
  ENDIF

ENDIF

IF ( NPR > 1 ) THEN
  IF (IRANK == IREAD) THEN
    IBUF(1) = KFILE_HANDLE
    IBUF(2) = KGRIB_HANDLE
    IBUF(3) = IU08(IPROPAGS)
  ENDIF

  CALL MPL_BROADCAST(IBUF(1:3),KROOT=IREAD, KTAG=4, CDSTRING='WVOPENSUBBATHY:')

  IF (IRANK /= IREAD) THEN
    KFILE_HANDLE = IBUF(1)
    KGRIB_HANDLE = IBUF(2)
    IU08(IPROPAGS) = IBUF(3)
  ENDIF

ENDIF

IF (LHOOK) CALL DR_HOOK('WVOPENSUBBATHY',1,ZHOOK_HANDLE)

END SUBROUTINE WVOPENSUBBATHY
