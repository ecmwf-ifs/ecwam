C ======================================================================
C
      SUBROUTINE URAOPC (csubna, ksub, llsuvi, lltest)
C
C ======================================================================
C
C****  *URAOPC*  - SET OUTPUT CONTROL SWITCHES.
C
C     B. HANSEN           ECMWF       MARCH   1992
C
C**   PURPOSE.
C     --------
C       SETTING OF OUTPUT CONTROL LOGICALS.
C
C**   INTERFACE.
C     ----------
C       CALL  URAOPC (CSUBNA, KSUB, LASUVI, LATEST)
C
C       I/   *CSUBNA*  NAME OF CALLING SUB-ROUTINE.
C       I/   *KSUB*    SUB-ROUTINE LEVEL OF CALLLING SUB-ROUTINE.
C        /O  *LASUVI*  LOGICAL SWITCH (TRUE: SV-OUTPUT,
C                                      FALSE: NO SV_OUTPUT).
C        /O  *LATEST*  LOGICAL SWITCH (TRU: TEST-OUTPUT,
C                                      FALSE: NO TEST_OUTPUT).
C
C     EXTERNALS.
C     ----------
C       NONE.
C
C     METHOD.
C     -------
C       NONE.
C
C     REFERENCES.
C     -----------
C        NONE
C
C ----------------------------------------------------------------------
C
#include <outctrl.h>
      DATA lflow / .FALSE. /
C ----------------------------------------------------------------------
C
      INTEGER i_get_unit, iun
      SAVE iun
C
C ----------------------------------------------------------------------
C
C* 1.0   INITIALIZATION.
C  ---------------------
C
 1000 CONTINUE
      lltest = .FALSE.
      llsuvi = .FALSE.
C
      IF (lflow) THEN
        iun=i_get_unit(-1,'flow','s','f',0)
        WRITE(iun,9000)' in routine >',csubna,'< ksub is >',ksub,'<'
      ENDIF
C
      IF (ksub .LE. msuvi) THEN
        llsuvi = .TRUE.
      ELSE
        DO 1001 j=1,jpctl
          IF (csuvi(j) .EQ. csubna ) THEN
              llsuvi = .TRUE.
              GOTO 1002
          ENDIF
 1001   CONTINUE
 1002   CONTINUE
      ENDIF
C
      IF (ksub .LE. mtest) THEN
        lltest = .TRUE.
      ELSE
        DO 1003 j=1,jpctl
          IF (ctest(j) .EQ. csubna ) THEN
              lltest = .TRUE.
              GOTO 1004
          ENDIF
 1003   CONTINUE
 1004   CONTINUE
      ENDIF
      RETURN
 9000 FORMAT(A13,A7,A11,I2,A1)
      END
