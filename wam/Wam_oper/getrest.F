      SUBROUTINE GETREST (IU06, IU12, IU14, IU15, IGL,
     1                    CDATE, CDATEF)
C ----------------------------------------------------------------------
C
C**** *GETREST*  GETS RESTART FILES
C
C     HEINZ GUNTHER   GKSS/ECMWF     MAY 1989
C
C*    PURPOSE.
C     --------
C
C       FETCHES RESTART FILES.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *GETREST (IU06, IU12, IU14, IU15, IGL, IFORCA)*
C
C         *IU06*    INTEGER    OUTPUT UNIT FOR PRINTER MESSAGES.
C         *IU12*    INTEGER    OUTPUT UNIT BLOCKS OF SPECTRA.
C         *IU14*    INTEGER    OUTPUT UNIT SECOND LAT OF BLOCKS.
C         *IU15*    INTEGER    OUTPUT UNIT LAST WINDFIELDS.
C         *IGL*     INTEGER    NUMBER OF BLOCKS.
C         *CDATE*  CHAR*10    DATE/TIME OF LAST FIELD IN FILES.
C         *CDATEF* CHAR*10    START DATE OF FORECAST PERIOD.
C
C     METHOD.
C     -------
C
C       NONE.
C
C     EXTERNALS.
C     ----------
C
C       *GSFILE*    - GETS OR SAVES FILES.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C      
C ---------------------------------------------------------------------
C
#include "comwami.h"
C
      CHARACTER*10 CDATE,CDATEF
C
C ----------------------------------------------------------------------
C
C*    1. ANALYSIS OR FORECATS FILES.
C        ---------------------------
C
 1000 CONTINUE
C      
      ICH = 0
C
      IF (CDATE.LE.CDATEF) THEN
C
C*    1.1 ANALYSIS FILES.
C         ---------------
C
C*    1.1.1 SPECTRA FILE.
C           -------------
C
         CALL GSFILE (IU06, IU12, ICH, CDATE, CDATEF, 'BLS', 'G')
C
C*    1.1.2 SECOND LATITUDE FILE IF MULTIBLOCK VERSION.
C           -------------------------------------------
C
         IF (IGL.GT.1) THEN
            CALL GSFILE (IU06, IU14, ICH, CDATE, CDATEF, 'SLA', 'G')
         ENDIF
C
C*    1.1.3 LAST WIND FIELDS.
C           -----------------
C
         CALL GSFILE (IU06, IU15, ICH, CDATE, CDATEF, 'LAW', 'G')
      ELSE
C
C*    1.2 FORECAST FILES.
C         ---------------
C
C*    1.2.1 SPECTRA FILE.
C           -------------
C
         CALL GSFILE (IU06, IU12, ICH, CDATE, CDATEF, 'BLS', 'G')
C
C*    1.2.2 SECOND LATITUDE FILE IF MULTIBLOCK VERSION.
C           -------------------------------------------
C
         IF (IGL.GT.1) THEN
            CALL GSFILE (IU06, IU14, ICH, CDATE, CDATEF, 'SLA', 'G')
         ENDIF
C
C*    1.2.3 LAST WIND FIELDS.
C           -----------------
C
         CALL GSFILE (IU06, IU15, ICH, CDATE, CDATEF, 'LAW', 'G')
      ENDIF
C
      RETURN
      END
