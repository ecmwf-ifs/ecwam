      SUBROUTINE OUTUBUF (IU08, IU18, IFORM)

C ----------------------------------------------------------------------
C
C**** *OUTUBUF* - ROUTINE TO WRITE COMMON UBUF TO DISK.
C
C     H.GUNTHER            ECMWF       04/04/1990
C
C*    PURPOSE.
C     -------
C
C       TO WRITE OUT THE COMPUTED COMMON UBUF.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *OUTUBUF (IU08, IU18, IFORM)*
C         *IU08*   - LOGICAL UNIT FOR  UNFORMATED WRITE
C         *IU18*   - LOGICAL UNIT FOR    FORMATED WRITE
C         *IFORM*  - FORMAT OPTION  = 1  UNFORMATED WRITE
C                                   = 2  FORMATED WRITE
C                                   OTHERWISE BOTH
C     METHOD.
C     -------
C
C       FORMATED AND OR UNFORMATED WRITE TO UNITS.
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "comubuf.h"
C
C ----------------------------------------------------------------------
C
  998 FORMAT(10I8)
  999 FORMAT(5E16.7)
C
C ----------------------------------------------------------------------
C
C*    1. WRITE COMMON UBUF UNFORMATED.
C        -----------------------------
C
 1000 CONTINUE
      IF (IFORM.NE.2) WRITE (IU08) KLAT, KLON
C
C ----------------------------------------------------------------------
C
C*    2. WRITE COMMON UBUF FORMATED.
C        ---------------------------
C
 2000 CONTINUE
      IF (IFORM.NE.1) WRITE (IU18,998) KLAT, KLON
      RETURN
      END
