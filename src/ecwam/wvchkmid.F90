! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE WVCHKMID (IU06, KGRIB_HANDLE, FILENAME)

! ----------------------------------------------------------------------

!**** *WVCHKMID* CHECKS MODEL GRIB IDENTIFIERS

! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWABORT , ONLY : WAM_ABORT
USE YOWGRIBHD, ONLY : CDATECLIM, CEXPVERCLIM
USE YOWGRIB  , ONLY : IGRIB_GET_VALUE

USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

IMPLICIT NONE
#include "abort1.intfb.h"

INTEGER(KIND=JWIM), INTENT(IN) :: IU06 
INTEGER(KIND=JWIM), INTENT(IN) :: KGRIB_HANDLE
CHARACTER(LEN=*),   INTENT(IN) :: FILENAME


INTEGER(KIND=JWIM) :: IRET, IDATECLIM, ITIMECLIM, IDATE, ITIME 

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

CHARACTER(LEN=4)  :: CEXPVER
CHARACTER(LEN=12) :: C12
CHARACTER(LEN=14) :: CDATE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WVCHKMID',0,ZHOOK_HANDLE)

IF ( KGRIB_HANDLE > 0 ) THEN

! CHECK MODEL IDENTIFIERS

  CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'expver', C12, KRET=IRET)
  IF (IRET /= 0) THEN
    CEXPVER='****'
  ELSE
    CEXPVER=C12(1:4)
  ENDIF

  IF ( CEXPVER /= CEXPVERCLIM ) THEN
    WRITE(IU06,*) '*****************************************'
    WRITE(IU06,*) '*  FATAL ERROR(S) IN SUB. WVCHKMID      *'
    WRITE(IU06,*) '*  ===============================      *'
    WRITE(IU06,*) '* CALLED FROM ', FILENAME
    WRITE(IU06,*) '* THE PROGRAM HAS DETECTED DIFFERENT    *'
    WRITE(IU06,*) '* MODEL GRIB CEXPVERCLIM                *' 
    WRITE(IU06,*) '* MAKE SURE YOU HAVE RUN PREPROC !!!!   *'
    WRITE(IU06,*) '* CEXPVER =      ',CEXPVER
    WRITE(IU06,*) '* CEXPVERCLIM = ',CEXPVERCLIM
    WRITE(IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.     *'
    WRITE(IU06,*) '* ---------------   --------------      *'
    WRITE(IU06,*) '*****************************************'
    CALL ABORT1
  ENDIF

  CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'dataDate', IDATE)
  CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'time', ITIME)
  CDATE = CDATECLIM
  READ(CDATE(1:8),'(I8)') IDATECLIM
  READ(CDATE(9:12),'(I4)') ITIMECLIM
  IF ( IDATE /= IDATECLIM .OR. ITIME /= ITIMECLIM ) THEN
    WRITE(IU06,*) '*****************************************'
    WRITE(IU06,*) '*  FATAL ERROR(S) IN SUB. WVCHKMID      *'
    WRITE(IU06,*) '*  ===============================      *'
    WRITE(IU06,*) '* CALLED FROM ', FILENAME
    WRITE(IU06,*) '* THE PROGRAM HAS DETECTED DIFFERENT    *'
    WRITE(IU06,*) '* MODEL GRIB DATECLIM.                  *' 
    WRITE(IU06,*) '* MAKE SURE YOU HAVE RUN PREPROC !!!!   *'
    WRITE(IU06,*) '* IDATE /= IDATECLIM ',IDATE, IDATECLIM
    WRITE(IU06,*) '* ITIME /= ITIMECLIM ',ITIME, ITIMECLIM 
    WRITE(IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.     *'
    WRITE(IU06,*) '* ---------------   --------------      *'
    WRITE(IU06,*) '*****************************************'
    CALL ABORT1
  ENDIF

ELSE
  CALL WAM_ABORT("GRIB HANDLE not defined !!",__FILENAME__,__LINE__)
ENDIF

IF (LHOOK) CALL DR_HOOK('WVCHKMID',1,ZHOOK_HANDLE)

END SUBROUTINE WVCHKMID
