      SUBROUTINE SAVREST (IU06, IU04, IU12, IU14, IU15, IGL,
     1                    CDATE, CDATEF, CDATER )
C
C ----------------------------------------------------------------------
C
C**** *SAVREST* - SAVES RESTART FILES.
C
C     HEINZ GUNTHER      GKSS/ECMWF     MAY 1989
C     J. BIDLOT          ECMWF          NOVEMBER 1996 : ADDED  CDATER
C     B. HANSEN                         APRIL   1997 : RESTART FACILITY
C     J. BIDLOT          ECMWF          SEPTEMBER : USE OF PBIO 
C
C*    PURPOSE.
C     --------
C
C       COPIES RESTART FILES.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *SAVREST (IU06, IU04, IU12, IU14, IU15, IGL,
C                        CDATE, CDATEF,CDATER)*
C         *IU06*   INTEGER   OUTPUT UNIT FOR PRINTER MESSAGES.
C         *IU04*   INTEGER   OUTPUT UNIT WAMINFO FILE.
C         *IU12*   INTEGER   OUTPUT UNIT BLOCKS OF SPECTRA.
C         *IU14*   INTEGER   OUTPUT UNIT SECOND LAT OF BLOCKS.
C         *IU15*   INTEGER   OUTPUT UNIT LAST WINDFIELDS.
C         *IGL*    INTEGER   NUMBER OF BLOCKS.
C         *CDATE*  CHAR*10   DATE/TIME OF LAST FIELD IN FILES.
C         *CDATEF* CHAR*10   START DATE OF FORECAST PERIOD.
C         *CDATER* CHAR*10   OUTPUT DATE FOR BOTH RESTART FILES
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
C ---------------------------------------------------------------------
C
!
#include "commesspass.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
#include "comwami.h"
C
      CHARACTER*10 CDATE,CDATEF,CDATER
C
C ---------------------------------------------------------------------
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

           CALL GSFILE (IU06, IU12, ICH, CDATE, CDATEF, 'BLS', 'S')
C
C*    1.1.2 SECOND LATITUDE FILE IF MULTIBLOCK.
C           -----------------------------------
C
         IF (IGL.GT.1) THEN
            CALL GSFILE (IU06, IU14, ICH, CDATE, CDATEF, 'SLA', 'S')
         ENDIF
C
C*    1.1.3 LAST WIND FIELDS.
C           -----------------
C
         IF((CDATE .EQ. CDATER ).OR.(CDATER .EQ. '0000000000')) THEN
           CALL GSFILE (IU06, IU15, ICH, CDATE, CDATEF, 'LAW', 'S')
         ENDIF

      ELSE
C
C*    1.1 FORECAST FILES.
C         ---------------
C
C*    1.2.1 SPECTRA FILE.
C           -------------
C
           CALL GSFILE (IU06, IU12, ICH, CDATE, CDATEF, 'BLS', 'S')
C
C*    1.2.2 SECOND LATITUDE FILE IF MULTIBLOCK.
C           -----------------------------------
C
         IF (IGL.GT.1) THEN
            CALL GSFILE (IU06, IU14, ICH, CDATE, CDATEF, 'SLA', 'S')
         ENDIF
C
C*    1.2.3 LAST WIND FIELDS.
C           -----------------
C
          IF((CDATE .EQ. CDATER ).OR.(CDATER .EQ. '0000000000')) THEN
            CALL GSFILE (IU06, IU15, ICH, CDATE, CDATEF, 'LAW', 'S')
          ENDIF
      ENDIF
C
      RETURN
      END
