      SUBROUTINE INTPOL (F3, F1, IJS, IJL, IG)

C ----------------------------------------------------------------------
C
C**** *INTPOL* - TRANSFORMATION OF SPECTRA FROM SIGMA TO OMEGA.
C
C     S.D.HASSELMANN      MPI            1.1.91
C     H. GUNTHER          GKSS/ECMWF     1.2.91  MODIFIED FOR CYCLE_4
C
C*    PURPOSE.
C     --------
C
C       TRANSFORMATION OF A MOVING COORDINATE SYSTEM TO AN ABSOLUTE
C       COORDINATE SYSTEM.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *INTPOL (F3, F1, IJS, IJL, IG)*
C         *F3*   - BLOCK OF SPECTRA (INPUT).
C         *F1*   - BLOCK OF SPECTRA (OUTPUT).
C         *IJS*  - INDEX OF FIRST GRIDPOINT.
C         *IJL*  - INDEX OF LAST GRIDPOINT.
C         *IG*   - BLOCK NUMBER.
C
C     METHOD.
C     -------
C
C       SCATTERING TO NEIGHBOURING POINT.
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
#include "comcurr.h"
C
#include "comfred.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
#include "comshal.h"
C
#include "comstat.h"
C
C ----------------------------------------------------------------------
C
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: F3,F1
C
C ----------------------------------------------------------------------
C
#include "parcons.h"
C
      PARAMETER (CO = 1.1, ALN00 = 24.1589,
     1           FRE0 = CO-1., PI2G = ZPI/G)
C
C*     VARIABLE.   TYPE.     PURPOSE.
C      ---------   -------   --------
C        *CO*      REAL      FREQUENCY RATIO.
C        *ALN00*   REAL      1. / LOG10(CO)
C
      DIMENSION FNEF(NIBLC), NEWF(NIBLC), NEWF1(NIBLC), KNEW(NIBLC),
     1          GWP(NIBLC), GWM(NIBLC), WAVN(NIBLC)
C
C ----------------------------------------------------------------------
C
C*    0. INITIAL OUTPUT ARRAY WITH ZERO.
C        -------------------------------
C
         F1 = 0.0
C
C*    1. LOOP OVER FREQUENCIES.
C        ----------------------
C
      DO 1000 M = 1, NFRE
         IF (ISHALLO.NE.1) THEN
            DO 1201 IJ = IJS, IJL
               WAVN(IJ) = TFAK(INDEP(IJ),M)/ZPI
 1201       CONTINUE
         ELSE
            DO 1202 IJ = IJS, IJL
               WAVN(IJ) = PI2G*FR(M)*FR(M)
 1202       CONTINUE
         ENDIF
C
C*    1.3 LOOP OVER DIRECTONS.
C         --------------------
C
         DO 1301 K = 1, NANG
C
C*    1.3.1 NEW FREQUENCY AND DIRECTION AT ALL GRIDPOINTS.
C           ----------------------------------------------
C
            DO 1311 IJ = IJS, IJL
               FNEF(IJ) = FR(M)
     1            + WAVN(IJ)*(COSTH(K)*V(IJ,IG) + SINTH(K)*U(IJ,IG))
               IF (FNEF(IJ).GT.0.) THEN
                  KNEW(IJ) = K
               ELSE
                  KNEW(IJ) = MOD(K+NANG/2-1,NANG) + 1
                  FNEF(IJ) = -FNEF(IJ)
               ENDIF
 1311       CONTINUE
C
C*    1.3.2 NEW FREQUENCY BIN NUMBER AT ALL GRIDPOINTS.
C           -------------------------------------------
C
            DO 1321 IJ = IJS, IJL
               IF (FNEF(IJ).LE.FR(1)/CO) THEN
                  NEWF(IJ) = -1
               ELSE
                  NEWF(IJ) =INT(LOG10(FNEF(IJ)/FR(1))*ALN00+1.000000001)
               ENDIF
 1321       CONTINUE
C
C*    1.3.3 INTERPOLATED ENERGY DENSITIES AT ALL GRIDPOINTS.
C           ------------------------------------------------
C
            DO 1331 IJ = IJS, IJL
               FNEW = FNEF(IJ)
               NEWM = NEWF(IJ)
               IF (NEWM.LT.NFRE.AND.NEWM.GE.1) THEN
                  NEWM1 = NEWM + 1
                  GWH = DFIM(M)/(FR(NEWM1)-FR(NEWM)) *F3(IJ,K,M)
                  GWM(IJ) = GWH*(FR(NEWM1)-FNEW)/DFIM(NEWM)
                  GWP(IJ) = GWH*(FNEW-FR(NEWM))/DFIM(NEWM1)
                  NEWF1(IJ) = NEWM1
               ELSEIF (NEWM.EQ.0) THEN
                  GWH = CO*DFIM(M)/(FRE0*FR(1)) * F3(IJ,K,M)
                  GWP(IJ) = GWH*(FNEW-FR(1)/CO)/DFIM(1)
                  NEWF (IJ) = -1
                  NEWF1(IJ) = 1
               ELSEIF (NEWM.EQ.NFRE) THEN
                  GWH = DFIM(M)/(FRE0*FR(NFRE)) * F3(IJ,K,M)
                  GWM(IJ) = GWH*(CO*FR(NFRE)-FNEW)/DFIM(NFRE)
                  NEWF1(IJ) = -1
               ELSE
                  NEWF (IJ) = -1
                  NEWF1(IJ) = -1
               ENDIF
 1331       CONTINUE
C
C*    1.3.4 NEW SPECTRUM AT ALL GRIDPOINTS.
C           -------------------------------
C
            DO 1341 IJ = IJS, IJL
               NEWM  = NEWF (IJ)
               NEWM1 = NEWF1(IJ)
               KH = KNEW(IJ)
               IF (NEWM .NE.-1)
     1            F1(IJ,KH,NEWM ) = F1(IJ,KH,NEWM ) + GWM(IJ)
               IF (NEWM1.NE.-1)
     1            F1(IJ,KH,NEWM1) = F1(IJ,KH,NEWM1) + GWP(IJ)
 1341       CONTINUE
C
C*    BRANCH BACK TO 1.3 FOR NEXT DIRECTION.
C
 1301    CONTINUE
C
C*    BRANCH BACK TO 1. FOR NEXT FREQUENCY.
C
 1000 CONTINUE

      RETURN
      END
