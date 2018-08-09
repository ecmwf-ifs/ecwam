SUBROUTINE WGRIBENOUT (IU06, ITEST, I1, I2, FIELD,                      &
     &                       ITABLE, IPARAM, KLEV, IK, IM,              &
     &                       CDATE, IFCST, MARSTYPE,                    &
     &                       LFDB, CDFDBSF, KFDB, LFDBOPEN, IU)

! ----------------------------------------------------------------------

!****  *WGRIBENOUT*  ENCODES WAM MODEL FIELD INTO GRIB CODE AND OUTPUT

!       J. BIDLOT    ECMWF JULY 2009: USE GRIB API 

!       PURPOSE.
!       --------
!         SUBROUTINE PACKS WAVE FIELDS INTO THE GRIB CODE
!         AND ARCHIVE INTO FDB OR INTO FILE.
!         !!!! PRESET_WGRIB_TEMPLATE NEEDS to BE CALLED BEFORE TO
!         !!!! INITIALISE THE GRIB TEMPLATES.

!**    INTERFACE.
!      ----------
!        *CALL* *WGRIBENOUT (IU06, ITEST, I1, I2, FIELD,
!                          ITABLE, IPARAM, KLEV, IK, IM,
!                          CDATE, IFCST,
!                          LFDB, CDFDBSF, KFDB, LFDBOPEN, IU,
!          *IU06*    LOGFILE OUTPUT UNIT.
!          *ITEST*   TEST OUTPUT GIVEN IF ITEST GT 2.
!          *I1*      FIRST DIMENSION OF FIELD.
!          *I2*      SECOND DIMENSION OF FIELD.
!          *FIELD*   FIELD TO BE PACKED.
!          *ITABLE*  GRIB TABLE NUMBER.
!          *IPARAM*  PARAMETER IDENTIFIER.
!          *KLEV*    REFERENCE LEVEL IN full METER
!                    (SHOULD BE 0 EXCEPT FOR 233 AND 245)
!          *IK*      DIRECTION INDEX,
!                    ONLY MEANINGFUL FOR SPECTRAL PARAMETERS.
!          *IM*      FREQUENCY INDEX,
!                    ONLY MEANINGFUL FOR SPECTRAL PARAMETERS.
!          *CDATE*   DATE AND TIME.
!          *IFCST*   FORECAST STEP IN HOURS.
!          *MARSTYPE* TYPE OF CURRENT FIELD
!           *LFDB*    LOGICAL SWITCH TO ACTIVATE FDB PROCESSING
!                    .T. TO ACTIVATE, .F. TO IGNORE
!          *IU*      LOGICAL UNIT FOR PACKED FIELD DATA.
!                    !!! IT WILL NEED TO BE OPEN AN
!                    !!! WITH GRIB_API SOFTWARE !!!

!      METHOD.
!      -------

!      EXTERNALS.
!      ----------

!      REFERENCES.
!      -----------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMPP   , ONLY : NPRECI
      USE GRIB_API_INTERFACE
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "wgribencode_model.intfb.h"
#include "wgribencode.intfb.h"
#include "wgribout.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IU06, ITEST, I1, I2
      INTEGER(KIND=JWIM), INTENT(IN) :: ITABLE, IPARAM, KLEV, IK, IM
      INTEGER(KIND=JWIM), INTENT(IN) :: IFCST
      INTEGER(KIND=JWIM), INTENT(IN) :: IU 
      INTEGER(KIND=JWIM), INTENT(INOUT) :: KFDB

      REAL(KIND=JWRB), INTENT(INOUT) :: FIELD(I1,I2)

      CHARACTER(LEN=2), INTENT(IN) :: MARSTYPE
      CHARACTER(LEN=14), INTENT(IN) :: CDATE
      CHARACTER(LEN=*), INTENT(INOUT) :: CDFDBSF

      LOGICAL, INTENT(IN) :: LFDB
      LOGICAL, INTENT(INOUT) :: LFDBOPEN

      INTEGER(KIND=JWIM) :: IGRIB_HANDLE
      INTEGER(KIND=JWIM) :: ISIZE
      INTEGER(KIND=JPKSIZE_T) :: KBYTES
      INTEGER(KIND=JWIM), ALLOCATABLE :: KGRIB_BUFR(:)

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('WGRIBENOUT',0,ZHOOK_HANDLE)

      IF(ITEST.GT.0) THEN
        WRITE(IU06,*) '   SUB. WGRIBENOUT CALLED FOR ',IPARAM
        CALL FLUSH(IU06)
      ENDIF

! ----------------------------------------------------------------------

!*    1. ENCODE RESULT
!        -------------

      CALL WGRIBENCODE_MODEL(IU06, ITEST, I1, I2, FIELD, &
     &                 ITABLE, IPARAM, KLEV, IK , IM,    &
     &                 CDATE, IFCST, MARSTYPE,           &
     &                 IGRIB_HANDLE)


!*    2. SAVE ENCODED RESULT
!        -------------------

      CALL IGRIB_GET_MESSAGE_SIZE(IGRIB_HANDLE,KBYTES)
      ISIZE=(KBYTES+NPRECI-1)/NPRECI
      ALLOCATE(KGRIB_BUFR(ISIZE))
      CALL IGRIB_GET_MESSAGE(IGRIB_HANDLE,KGRIB_BUFR)

      CALL WGRIBOUT(IU06, ITEST,                       &
     &              LFDB, CDFDBSF, KFDB, LFDBOPEN,     &
     &              IU,                                &
     &              IGRIB_HANDLE, ISIZE, KGRIB_BUFR )


!     RELEASE GRIB HANDLE
      CALL IGRIB_RELEASE(IGRIB_HANDLE)
      DEALLOCATE(KGRIB_BUFR)

      IF (LHOOK) CALL DR_HOOK('WGRIBENOUT',1,ZHOOK_HANDLE)

END SUBROUTINE WGRIBENOUT
