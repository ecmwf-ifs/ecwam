      SUBROUTINE GRFIELD (WHME, USME, NSTART, NEND)
C
C--------------------------------------------------------------------
C
C**** *GRFIELD* - READS GRIDDED MEASUREMENTS AND DISTRIBUTES THEM.
C
C     B.HANSEN       ECMWF       JANUARY 1998
C
C     PURPOSE.
C     --------
C
C       PREPARE THE DATA FOR OPTIMAL INTERPOLATION.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *GRFIELD (WHME, USME, LDFAILD)*
C
C       *WHME*   REAL    OUTPUT FIELD OF MEASURED WAVE HEIGHTS.
C       *USME*   REAL    OUTPUT FIELD OF MEASURED SURFACE WINDS.
C       *NSTART* INTEGER INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
C       *NEND*   INTEGER INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
C
C
C     METHOD.
C     -------
C
C        READ PRE PREPARED SATELLITE DATA FROM FILE AND DISTRIBUTE TO
C        OTHER PES.
C
C      EXTERNALS.
C      ----------
C
C         MPBCASTGRDFLD  - BROADCAST 2-D DATA FIELD FROM 1 PE TO ALL
C                          OTHERS.
C
C      MODIFICATONS.
C      -------------
C         NONE.
C
C-----------------------------------------------------------------------
C
#include "param.h"
C
#include "comgrid.h"
C
#include "comice.h"
C
#include "commap.h"
C
#include "commesspass.h"
C
#include "commpp.h"
C
#include "comstat.h"
C
#include "comtest.h"
C
#include "comintp.h"
C
C ----------------------------------------------------------------------
C
      LOGICAL LLFAILED, LLOPENED
      INTEGER IREAD
      CHARACTER* 3 CLFID
      CHARACTER*12 CTIME
      REAL,   DIMENSION(NGX,NGY) :: WHME, USME
      INTEGER,DIMENSION(NPROC)   :: NSTART, NEND

C
C*     VARIABLE     TYPE      PURPOSE.
C      --------     ----      --------
C
C     *WHME*        REAL      GRID WAVE HEIGHT.
C     *USME*        REAL      GRID FRICTION VELOCITY.
C     *IREAD*       INTGER    PROCESSOR WHICH WILL READ THE INPUT FILES.
C     *CLFID        CHARACTER THREE CHARACTER FILE IDENTIFIER.
C
      DATA CLFID, IREAD /'RFL', 1/
C
C ----------------------------------------------------------------------
C
C*    1. OPEN FILE AND READ.
C        -------------------
C
      IF(LMESSPASS .AND. IRANK .EQ. IREAD .OR. .NOT.LMESSPASS) THEN
        IFIELD = 0
        LLFAILED = .FALSE.
        IFAIL = 1
C
C*       FIND SPARE UNIT TO ATTACH INPUT FILE.
C
        UNIT: DO IUME=33,99,1
          INQUIRE ( UNIT=IUME, OPENED=LLOPENED)
          IF ( .NOT. LLOPENED) THEN
             WRITE(IU06,'("  GRFIELD : TAKES UNIT ", I2,
     .                    " AS ALTIMETER DATA INPUT UNIT")') IUME
             EXIT UNIT
          ENDIF
        ENDDO UNIT
C
C*       IF FOUND: SEEK INPUT FILE, READ HEADER AND FIELDS.
C*                 IFIELD MAY BE MISSING INDICATING THERE ARE NO FIELDS.
C
        IF (IUME.LT.100) THEN
          CALL CONFILE (IU06, IUME, CDTPRO, CLFID, IFAIL)
          IF(IFAIL.EQ.0) THEN
            READ (IUME,ERR=9000) CTIME, IFIELD
            READ (IUME,ERR=9000) WHME, USME
            CLOSE (IUME)
          ELSE
            LLFAILED = .TRUE.
            WRITE(IU06,'("WARNING IN GRFIELD:  W A R N I N G  !!!!!!")')
            WRITE(IU06,'("NO FILES FOUND"                            )')
            WRITE(IU06,'("WARNING IN GRFIELD:  W A R N I N G  !!!!!!")')
          ENDIF
        ELSE
C
C*       IF NOT FOUND: ISSUE WARNING MESSAGE.
C
          LLFAILED = .TRUE.
          WRITE(IU06,'("WARNING IN GRFIELD:  W A R N I N G  !!!!!!!!")')
          WRITE(IU06,'("NO SPARE UNIT TO READ SATELLITE DATA FOUND"  )')
          WRITE(IU06,'("WARNING IN GRFIELD:  W A R N I N G  !!!!!!!!")')
        ENDIF
C
C*       IF ALL ATTEMPTS OF THE NEW METHOD FAILED TRY THE OLD CALL TO
C*       GRDATA.
C
 9000   IF ( LLFAILED .OR. IFIELD .EQ. 0 ) THEN
          WRITE(IU06,'("WARNING IN GRFIELD:    W A R N I N G  !!!!!!")')
          WRITE(IU06,'("NO GRIDDED FIELDS FOUND IN FILE"             )')
          WRITE(IU06,'("TRY THE OLD CALL TO GRDATA                  ")')
          WRITE(IU06,'("WARNING IN GRFIELD:    W A R N I N G  !!!!!!")')
          IF (IFAIL .EQ. 0) CLOSE (IUME)
          CALL GRDATA (WHME, USME, NSTART, NEND) 
        ENDIF
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    2. SEND FIELDS TO ALL OTHER PES.
C        -----------------------------
C
      IF(LMESSPASS) THEN
        CALL MPBCASTGRDFLD(IREAD,KTAG,NGX,NGY,WHME)
        KTAG=KTAG+1
        CALL MPBCASTGRDFLD(IREAD,KTAG,NGX,NGY,USME)
        KTAG=KTAG+1
      ENDIF
C
C ----------------------------------------------------------------------
C
C     3. APPLY SCREENING OF DATA BASED ON MODEL FIELD
C        --------------------------------------------
C        
C        EXCLUDE POINTS WHICH ARE ON MODEL ICE
C
      IF(LMESSPASS) THEN
         DO IJ = NSTART(1), NEND(NPROC)
            IF (ICEMASK(IJ,1).EQ.0) THEN
              IX = IXLG(IJ,1)
              IY = NGY- KXLT(IJ,1) +1
              USME(IX,IY)=-3.
              WHME(IX,IY)=-3.
            ENDIF
         ENDDO
      ELSE
         DO IG = 1, IGL
            DO IJ = IJS(IG),IJL(IG)
               IF (ICEMASK(IJ,IG).EQ.0) THEN
                 IX = IXLG(IJ,IG)
                 IY = NGY- KXLT(IJ,IG) +1
                 USME(IX,IY)=-3.
                 WHME(IX,IY)=-3.
               ENDIF
            ENDDO
         ENDDO
      ENDIF

C*       OBTAIN FRICTION VELOCITIES FROM OBSERVED WINDS AND THE DRAG.
C
      FV: DO J=1,NGY
        DO I=1,NGX
          IF (USME(I,J) .GT. 0) USME(I,J) = SQRT(CDG(I,J)) * USME(I,J)
          IF (CDG(I,J) .LT. 0.000) THEN
            USME(I,J)=-3.
            WHME(I,J)=-3.
          ENDIF
        ENDDO
      ENDDO FV
C
      RETURN
      END
