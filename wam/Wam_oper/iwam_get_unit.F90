FUNCTION IWAM_GET_UNIT (KUSO, CDNAME, CDACCESS, CDFORM, KRECL)

!     FIND A FREE UNIT TO ACCESS A FILE (CDNAME) AND OPEN IT.

!--------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

!--------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM) :: IWAM_GET_UNIT
      INTEGER(KIND=JWIM) :: KRECL, KUSO

      INTEGER(KIND=JWIM) :: MINUNIT, MAXUNIT, ILEN, JUME

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      LOGICAL :: LLOPENED, LLUSED

      CHARACTER(LEN=*) :: CDNAME
      CHARACTER(LEN=1) :: CDACCESS, CDFORM
      CHARACTER(LEN=20) :: CLACCESS, CLFORM, CLPOSITION

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('IWAM_GET_UNIT',0,ZHOOK_HANDLE)

      CLACCESS   = "SEQUENTIAL"
      CLFORM     = "FORMATTED"
      CLPOSITION = "ASIS"

      MINUNIT=140
      MAXUNIT=300

      ILEN=LEN(CDNAME)
      if(ILEN == 0) then
        WRITE(*,*) '!!!!!!!!! PROBLEM IN WAM IWAM_GET_UNIT !!!!!!'
        WRITE(*,*) '!!!!!!!!! LENGTH OF FILENAME=0  !!!!!!'
        CALL ABORT1
      ENDIF
      ILEN=LEN_TRIM(CDNAME)
      if(ILEN == 0) then
        WRITE(*,*) '!!!!!!!!! PROBLEM IN WAM IWAM_GET_UNIT !!!!!!'
        WRITE(*,*) '!!!!!!!!! FILENAME ALL BLANKS   !!!!!!'
        CALL ABORT1
      ENDIF
         
      INQUIRE ( FILE=CDNAME, OPENED=LLOPENED, NUMBER=JUME)

      IF (.NOT. LLOPENED) THEN
        IF ( CDACCESS .EQ. 'D' .OR. CDACCESS .EQ. 'd') CLACCESS="DIRECT"
        IF ( CDFORM   .EQ. 'U' .OR. CDFORM   .EQ. 'u') CLFORM="UNFORMATTED"
        IF ( CDACCESS .EQ. 'A' .OR. CDACCESS .EQ. 'a') CLPOSITION="APPEND" 
 
        DO JUME=MINUNIT,MAXUNIT,1
          INQUIRE ( UNIT=JUME, OPENED=LLUSED)
          IF ( .NOT. LLUSED) EXIT
        ENDDO
 
        IF(JUME.GT.MAXUNIT) THEN
          WRITE(*,*) '!!!!!!!!! PROBLEM IN WAM IWAM_GET_UNIT !!!!!!'
          WRITE(*,*) '!!!!!!!!! NO FREE UNIT FOUND BETWEEN !'
          WRITE(*,*) '!!!!!!!!! ',MINUNIT,' AND ',MAXUNIT
          WRITE(*,*) '!!!!!!!!! CALL ABORT                 !'
          CALL ABORT1
        ENDIF
        IF (KUSO .GE. 0 ) WRITE(KUSO,2004)        &
     &   '  IWAM_GET_UNIT: U=',JUME,              &
     &   ' F=',CDNAME(1:LEN_TRIM(CDNAME)),        &
     &   ' F=',CLFORM(1:LEN_TRIM(CLFORM)),        &
     &   ' A=',CLACCESS(1:LEN_TRIM(CLACCESS)),    &
     &   ' P=',CLPOSITION(1:LEN_TRIM(CLPOSITION))
 2004   FORMAT(A,I3,A,A,A,A,A,A,A,A)

        IF ( CLACCESS(1:6) .EQ. "DIRECT" ) THEN
          OPEN (UNIT=JUME, FILE=CDNAME, FORM=CLFORM, ACCESS=CLACCESS, RECL=KRECL, ACTION=READ)
        ELSE
          OPEN (UNIT=JUME, FILE=CDNAME, FORM=CLFORM, ACCESS=CLACCESS, POSITION=CLPOSITION, ACTION=READ)
        ENDIF
      ENDIF

      IWAM_GET_UNIT=JUME

      IF (LHOOK) CALL DR_HOOK('IWAM_GET_UNIT',1,ZHOOK_HANDLE)

END FUNCTION IWAM_GET_UNIT
