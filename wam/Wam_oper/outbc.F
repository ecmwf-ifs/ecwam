      SUBROUTINE OUTBC (FL3, IJS, IJL, IG, IU19, THWNEW, USNEW)

C ----------------------------------------------------------------------
C
C**** *OUTBC* - OUTPUT OF THE COARSE GRID BOUNDARY VALUES.
C
C     R. PORTZ     MPI          JANUARY 1991
!     J. BIDLOT    ECMWF        FEBRARY 1996  MESSAGE PASSING
C
C*    PURPOSE.
C     --------
C
C        WRITE THE BOUNDARY VALUE OUTPUT FILE.
C
C**   INTERFACE.
C     ----------
C
C    *CALL* *OUTBC (FL3, IJS, IJL, IG, IU19, THWNEW, USNEW)*
C      *FL3*     - BLOCK OF SPECTRA.
C      *IJS*     - INDEX OF FIRST GRIDPOINT.
C      *IJL*     - INDEX OF LAST GRIDPOINT.
C      *IG*      - BLOCK NUMBER.
C      *IU19*    - OUTPUT UNIT OF BOUNDARY VALUES.
C      *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
C                  NOTATION (POINTING ANGLE OF WIND VECTOR,
C                  CLOCKWISE FROM NORTH).
C      *USNEW*     NEW FRICTION VELOCITY IN M/S.

C
C     METHOD.
C     -------
C
C       SEQUENCIAL UNFORMATED WRITE TO UNIT.
C
C     EXTERNALS.
C     ----------
C
C       *FEMEAN*    - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT.
C       *SEMEAN*    - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.
C       *STHQ*      - COMPUTATION OF MEAN WAVE DIRECTION AT EACH
C                     GRID POINT.
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
#include "comcbou.h"
C
#include "comfred.h"
C
#include "commap.h"
C
#include "commean.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
#include "comstat.h"
C
#include "comtest.h"
C
#include "comwind.h"
C
C ----------------------------------------------------------------------
C
!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      DIMENSION FL3(NINF-1:NSUP,NANG,NFRE)

      REAL,DIMENSION(NINF:NSUP) :: THWNEW,USNEW
C
C ----------------------------------------------------------------------
C
C*    1. COMPUTE MEAN PARAMETERS.
C        ------------------------
C
 1000 CONTINUE
      CALL SEMEAN (FL3, IJS, IJL)
      CALL FEMEAN (FL3, IJS, IJL)
      CALL STHQ   (FL3, IJS, IJL)
      IF (ITEST.GE.3) THEN
         IF (ITESTB.GE.IG) THEN
            WRITE(IU06,*) '      SUB. OUTBC: INTEGRATED',
     1          ' PARAMETERS COMPUTED FOR BOUNDARY POINTS OUTPUT'
         ENDIF
      ENDIF
C
C*    2. LOOP OVER BOUNDARY POINTS.
C        --------------------------
C
 2000 CONTINUE
      DO 2001 NGOU=1,NBOUNC
         IF (IG .LT. IGARC(NGOU)) GOTO 2002
         IF (IG.EQ.IGARC(NGOU)) THEN
            IJ = IJARC(NGOU)
            IX = IXLG(IJ,IG)
            KX = KXLT(IJ,IG)
            XLON = AMOWEP + REAL(IX-1)*ZDELLO(KX)
            XLAT = AMOSOP + REAL(KX-1)*XDELLA
C
C*    2.1 WRITE INFORMATION TO FILE IU19.
C         -------------------------------
C
            XIJ   = REAL(IJ)
            XIG   = REAL(IG)
            WRITE(IU19) XLON, XLAT, CDTPRO, EMEAN(IJ),
     +                  THQ(IJ), FMEAN(IJ), USNEW(IJ), THWNEW(IJ)
            WRITE(IU19) ((FL3(IJ,K,M),K=1,NANG),M=1,NFRE)
         ENDIF
 2001 CONTINUE
 2002 CONTINUE

      RETURN
      END
