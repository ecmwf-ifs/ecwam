      SUBROUTINE OUTBS (FL3, FL, IJS, IJL, IG, IGL, IU25, IU26,
     1                  THWOLD,USOLD,U10NEW,THWNEW,USNEW,TAUW,
     2                  NSTART,NEND)
C ----------------------------------------------------------------------
C
C**** *OUTBS* - MODEL OUTPUT FROM BLOCK TO FILE, PRINTER AND COMMON.
C
C     H. GUNTHER         GKSS/ECMWF         JUNE 1990
C     J. BIDLOT          ECMWF   FEBRUARY 1996 MESSAGE PASSING
C
C*    PURPOSE.
C     --------
C
C       CONTROL OUTPUT OF WAVE AND WIND FIELDS.
C
C
C**   INTERFACE.
C     ----------
C   *CALL* *OUTBS (FL3, FL, IJS, IJL, IG, IGL, IU25, IU26,
C    1                  THWOLD,USOLD,U10NEW,THWNEW,USNEW,TAUW)
C    2                  NSTART,NEND)*
C      *FL3*    - SPECTRUM.
C      *FL*     - WORKING ARRAY FOR SWELL OUTPUT.
C      *IJS*    - INDEX OF FIRST GRIDPOINT
C      *IJL*    - INDEX OF LAST GRIDPOINT
C      *IG*     - BLOCK NUMBER
C      *IGL*    - TOTAL NUMBER OF BLOCKS.
C      *IU25*   - OUTPUT UNIT FOR SPECTRA.
C      *IU26*   - OUTPUT UNIT FOR SEA AND SWELL SPECTRA.
C      *THWOLD*    INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
C      *USOLD*     INTERMEDIATE STORAGE OF MODULUS OF FRICTION
C                  VELOCITY.
C      *U10NEW*    NEW WIND SPEED IN M/S.
C                  WIND VELOCITY.
C      *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
C                  NOTATION (POINTING ANGLE OF WIND VECTOR,
C                  CLOCKWISE FROM NORTH).
C      *USNEW*     NEW FRICTION VELOCITY IN M/S.
C      *TAUW*      WAVE STRESS IN (M/S)**2
C      *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
C      *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
C
C
C     EXTERNALS.
C     ----------
C
C       *FEMEAN*    - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT.
C       *OUTERS*    - OUTPUT OF ERS-1 UWA COLOCATION SPECTRA.
C       *OUTGRID*   - SAVE BLOCKED PARAMETERS INTO GRID ARRAYS.
C       *OUTSPP*    - OUTPUT OF SPECTRA AT SELECTED POINTS.
C       *SEMEAN*    - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.
C       *SEPWISW*   - SEPARATION OF WINDSEA AND SWELL.
C       *STHQ*      - COMPUTATION OF MEAN WAVE DIRECTION AT EACH
C                     GRID POINT.
C
C     METHOD.
C     -------
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
#include "comstat.h"
C
#include "comtest.h"
C
#include "commpp.h"
C
C ----------------------------------------------------------------------
C
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
      REAL,DIMENSION(NINF-1,NSUP,NANG,NFRE) :: FL3, FL
      REAL,DIMENSION(NINF:NSUP,NBLO) :: THWOLD,USOLD,TAUW
      REAL,DIMENSION(NINF:NSUP) :: U10NEW,THWNEW,USNEW

C
C ----------------------------------------------------------------------
C
C*    1. COMPUTE MEAN PARAMETERS.
C        ------------------------
C
      CALL SEMEAN (FL3, IJS, IJL)
      CALL FEMEAN (FL3, IJS, IJL)
      CALL STHQ   (FL3, IJS, IJL)
      CALL MEANSQS (FL3, IJS, IJL, USNEW)
      IF (ITEST.GE.3) THEN
         IF (ITESTB.GE.IG) THEN
            WRITE(IU06,*) '      SUB. OUTBS: INTEGRATED',
     1                    ' PARAMETERS COMPUTED FOR OUTPUT'
         ENDIF
      ENDIF
C
C*    2. WINDSEA SWELL SEPARATION.
C        -------------------------
C
      IF (CDTINTS.EQ.CDTPRO .OR. CDTSPS.EQ.CDTPRO) THEN
         CALL SEPWISW (FL3, FL, IJS, IJL, IG,THWOLD,USOLD,USNEW)
         IF (ITEST.GE.3) THEN
            IF (ITESTB.GE.IG) THEN
               WRITE(IU06,*) '      SUB. OUTBS: SWELL /SEA',
     1                       ' SEPARATION DONE'
            ENDIF
         ENDIF
      ENDIF
C
C*    3. OUTPUT OF SPECTRA AT SELECTED GRID POINTS.
C        ------------------------------------------
C
      IF (CDTSPT.EQ.CDTPRO .OR. CDTSPS.EQ.CDTPRO) THEN
         CALL OUTSPP (FL3, FL, IG, IU25, IU26,NSTART,NEND,
     1                U10NEW,THWNEW,USNEW)
         IF (ITEST.GE.3) THEN
            IF (ITESTB.GE.IG) THEN
               WRITE(IU06,*) '      SUB. OUTBS: OUTPUT OF SPECTRA',
     1                       ' AT SELECTED POINTS DONE'
            ENDIF
         ENDIF
      ENDIF
C
C*    4. OUTPUT OF INTEGRATED PARAMETERS.
C        --------------------------------
C
      IF (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO) THEN
         CALL OUTGRID (FL3, IJS, IJL, IG,U10NEW,THWNEW,USNEW,TAUW,
     1                 NSTART,NEND)
         IF (ITEST.GE.3) THEN
            IF (ITESTB.GE.IG) THEN
               WRITE(IU06,*) '      SUB. OUTBS: OUTPUT OF ',
     1                       'INTEGRATED PARAMETERS DONE'
            ENDIF
         ENDIF
      ENDIF
C
C*    5. OUTPUT OF SPECTRA FOR ERS_1 COLLOCATION.
C        ----------------------------------------
C
      CALL OUTERS (FL3, IG, IGL, CDTPRO,NSTART,NEND,U10NEW,THWNEW,USNEW)
      IF (ITEST.GE.3) THEN
         IF (ITESTB.GE.IG) THEN
            WRITE(IU06,*) '      SUB. OUTBS: OUTPUT OF SPECTRA',
     1                    ' FOR ERS-1 COLOCATION DONE'
         ENDIF
      ENDIF

      RETURN
      END
