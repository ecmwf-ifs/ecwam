      SUBROUTINE WGRIBOUT (IU06, ITEST,
     &                     LFDB, CDFDBSF, KFDB, LFDBOPEN,
     &                     IU,
     &                     IGRIB_HANDLE, ISIZE, KGRIB_BUFR) 
! ----------------------------------------------------------------------

!****  *WGRIBOUT* OUTPUT GRIB MESSAGE. 

!       PURPOSE.
!       --------

!      OUTPUT GRIB MESSAGE TO FDB (ECMWF) OR TO FILE.

!**    INTERFACE.
!      ----------
!        *CALL* *WGRIBOUT (IU06, ITEST,
!                          LFDB, CDFDBSF, KFDB, LFDBOPEN, IU,
!                          IGRIB_HANDLE, ISIZE, KGRIB_BUFR) 
!          *IU06*    LOGFILE OUTPUT UNIT.
!          *ITEST*   TEST OUTPUT GIVEN IF ITEST GT 2.
!          *LFDB*    LOGICAL SWITCH TO ACTIVATE FDB PROCESSING
!                    .T. TO ACTIVATE, .F. TO IGNORE
!          *CDFDBSF* FDB ROOT DIRECTORY.
!          *KFDB*    FDB DATA BASE REFERENCE.
!          *LFDBOPEN*  TRUE IF DATA BASE IS OPENED.
!          *IU*      UNIT TO WRITE O FILE
!                    !!! IT WILL NEED TO BE OPEN AN
!                    !!! WITH GRIB_API SOFTWARE !!!
!          *IGRIB_HANDLE* GRIB MESSAGE HANDLE.
!          *ISIZE*   GRIB MESSAGE SIZE.
!          *KGRIB_BUFR* GRIB MESSAGE. 

!      METHOD.
!      -------

!      EXTERNALS.
!      ----------

!      REFERENCES.
!      -----------

! ----------------------------------------------------------------------

      USE GRIB_API_INTERFACE
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "abort1.intfb.h"
#include "wgrib2fdb.intfb.h"

      INTEGER, INTENT(IN) :: IGRIB_HANDLE, ISIZE
      INTEGER, INTENT(IN) :: IU06, ITEST, IU
      INTEGER, INTENT(INOUT) :: KFDB
      INTEGER, DIMENSION(ISIZE), INTENT(IN) :: KGRIB_BUFR
      INTEGER :: IERR, ITABPAR, ICLASS
      INTEGER(KIND=JPKSIZE_T) :: KBYTES

      REAL :: ZHOOK_HANDLE

      CHARACTER(LEN=*), INTENT(IN) :: CDFDBSF
      CHARACTER(LEN=4) :: CSTREAM
      CHARACTER(LEN=4) :: CEXPVER
      CHARACTER(LEN=12) :: C12

      LOGICAL, INTENT(IN) :: LFDB
      LOGICAL, INTENT(INOUT) :: LFDBOPEN

! ----------------------------------------------------------------------

#ifdef ECMWF 
      IF (LHOOK) CALL DR_HOOK('WGRIBOUT',0,ZHOOK_HANDLE)
#endif

      CALL GSTATS(1709,0)

!*    1. SAVE ENCODED RESULT
!        -------------------

      IF(ITEST.GT.0) THEN
        CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'paramId',ITABPAR)
        CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'expver',C12)
        CEXPVER=C12(1:4)
        CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'stream',C12)
        CSTREAM=C12(1:4)
        CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'type',ICLASS)
        WRITE(IU06,*)'  '
        WRITE(IU06,*)'   SUB. WGRIBOUT : PARAMETER ',ITABPAR,
     &                                 ' EXPVER=',CEXPVER,
     &                                 ' STREAM= ', CSTREAM,
     &                                 ' CLASS = ', ICLASS
        CALL FLUSH(IU06)
      ENDIF


      IF(LFDB) THEN
!*     2.1 WRITE DATA TO FDB
!          -----------------

        CALL WGRIB2FDB (IU06, ITEST,
     &                  IGRIB_HANDLE, ISIZE, KGRIB_BUFR,
     &                  CDFDBSF, KFDB, LFDBOPEN,
     &                  IERR)
        IF(ITEST.GT.0) THEN
          WRITE(IU06,*) '   SUB. WGRIBOUT: GRIB DATA WRITTEN TO FDB',
     &     '(REFERENCE=',KFDB,')'
           CALl FLUSH(IU06)
        ENDIF

        IF(IERR.NE.0)THEN
          WRITE(IU06,*) ' ------------------------'
          WRITE(IU06,*) ' ERROR ACCESSING FDB '
          WRITE(IU06,*) ' FDB ERROR CODE IS ',IERR
          WRITE(IU06,*) ' ------------------------'
          WRITE(*,*) ' ------------------------'
          WRITE(*,*) ' ERROR ACCESSING FDB '
          WRITE(*,*) ' FDB ERROR CODE IS ',IERR
          WRITE(*,*) ' ------------------------'
          CALL ABORT1
        ENDIF

      ELSE

!*      2.2 WRITE PACKED DATA TO FILE CONNECTED TO UNIT IU.
!           ----------------------------------------------

        CALL IGRIB_GET_MESSAGE_SIZE(IGRIB_HANDLE,KBYTES)
        CALL IGRIB_WRITE_BYTES(IU,KGRIB_BUFR,KBYTES)

        IF (ITEST.GT.2) THEN
            WRITE(IU06,*) '   SUB. WGRIBOUT : DATA WRITTEN TO FILE'
        ENDIF
      ENDIF

      CALL GSTATS(1709,1)

#ifdef ECMWF 
      IF (LHOOK) CALL DR_HOOK('WGRIBOUT',1,ZHOOK_HANDLE)
#endif

      RETURN
      END SUBROUTINE WGRIBOUT
