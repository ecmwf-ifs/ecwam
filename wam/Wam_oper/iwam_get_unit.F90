INTEGER(KIND=JWIM) FUNCTION IWAM_GET_UNIT (KUSO, CDNAME, CDACCESS, CDFORM, KRECL, CDACTION)

!     FIND A FREE UNIT TO ACCESS A FILE (CDNAME) AND OPEN IT.

!--------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

!--------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KUSO
      CHARACTER(LEN=*),  INTENT(IN) :: CDNAME, CDACTION
      CHARACTER(LEN=1), INTENT(IN) :: CDACCESS, CDFORM
      INTEGER(KIND=JWIM), INTENT(IN) :: KRECL

!!      INTEGER(KIND=JWIM) :: IWAM_GET_UNIT

      INTEGER(KIND=JWIM) :: MINUNIT, MAXUNIT, ILEN, JUME

      LOGICAL :: LLOPENED, LLUSED

      CHARACTER(LEN=20) :: CLACCESS, CLFORM, CLPOSITION

! ----------------------------------------------------------------------

      CLACCESS   = "SEQUENTIAL"
      CLFORM     = "FORMATTED"
      CLPOSITION = "ASIS"

      MINUNIT=140
      MAXUNIT=300

      ILEN=LEN(CDNAME)
      IF ( ILEN == 0) then
        WRITE(*,*) '!!!!!!!!! PROBLEM IN WAM IWAM_GET_UNIT !!!!!!'
        WRITE(*,*) '!!!!!!!!! LENGTH OF FILENAME=0  !!!!!!'
        CALL ABORT1
      ENDIF
      ILEN=LEN_TRIM(CDNAME)
      IF ( ILEN == 0) then
        WRITE(*,*) '!!!!!!!!! PROBLEM IN WAM IWAM_GET_UNIT !!!!!!'
        WRITE(*,*) '!!!!!!!!! FILENAME ALL BLANKS   !!!!!!'
        CALL ABORT1
      ENDIF
         
      INQUIRE ( FILE=CDNAME, OPENED=LLOPENED, NUMBER=JUME)

      IF (.NOT. LLOPENED) THEN
        IF ( CDACCESS == 'D' .OR. CDACCESS == 'd') CLACCESS="DIRECT"
        IF ( CDFORM   == 'U' .OR. CDFORM   == 'u') CLFORM="UNFORMATTED"
        IF ( CDACCESS == 'A' .OR. CDACCESS == 'a') CLPOSITION="APPEND" 
 
        DO JUME=MINUNIT,MAXUNIT,1
          INQUIRE ( UNIT=JUME, OPENED=LLUSED)
          IF ( .NOT. LLUSED) EXIT
        ENDDO
 
        IF (JUME > MAXUNIT) THEN
          WRITE(*,*) '!!!!!!!!! PROBLEM IN WAM IWAM_GET_UNIT !!!!!!'
          WRITE(*,*) '!!!!!!!!! NO FREE UNIT FOUND BETWEEN !'
          WRITE(*,*) '!!!!!!!!! ',MINUNIT,' AND ',MAXUNIT
          WRITE(*,*) '!!!!!!!!! CALL ABORT                 !'
          CALL ABORT1
        ENDIF
        IF (KUSO >= 0 ) WRITE(KUSO,2004)        &
     &   '  IWAM_GET_UNIT: U=',JUME,              &
     &   ' F=',CDNAME(1:LEN_TRIM(CDNAME)),        &
     &   ' F=',CLFORM(1:LEN_TRIM(CLFORM)),        &
     &   ' A=',CLACCESS(1:LEN_TRIM(CLACCESS)),    &
     &   ' P=',CLPOSITION(1:LEN_TRIM(CLPOSITION))
 2004   FORMAT(A,I3,A,A,A,A,A,A,A,A)

        IF ( CLACCESS(1:6) == "DIRECT" ) THEN
          OPEN (UNIT=JUME, FILE=CDNAME, FORM=CLFORM, ACCESS=CLACCESS, RECL=KRECL, ACTION=CDACTION)
        ELSE
          OPEN (UNIT=JUME, FILE=CDNAME, FORM=CLFORM, ACCESS=CLACCESS, POSITION=CLPOSITION, ACTION=CDACTION)
        ENDIF
      ENDIF

      IWAM_GET_UNIT=JUME

END FUNCTION IWAM_GET_UNIT
