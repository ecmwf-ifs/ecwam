      SUBROUTINE READPRE (IU07, IREFRA)

C ----------------------------------------------------------------------
C
C**** *READPRE*  READ GRID OUTPUT FROM PREPROC.
C
C     H. GUNTHER      GKSS/ECMWF     MAY 1990
C
C*    PURPOSE.
C     --------
C
C       INPUT OF PREPROC GRID OUTPUT.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *READPRE (IU07, IREFRA)*
C          *IU07 *  - INPUT UNIT OF PREPROC GRID FILE.
C          *IREFRA *- REFRACTION OPTION.
C
C     METHOD.
C     -------
C
C       UNFORMATED READ FROM UNIT.
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
#include "comcoup.h"
C
#include "comcurr.h"
C
#include "comfred.h"
C
#include "comgrid.h"
C
#include "comindn.h"
C
#include "commap.h"
C
#include "comcout.h"
C
#include "comshal.h"
C
#include "comtabl.h"
C
C ----------------------------------------------------------------------
C
C*    1. READ COMMON FREDIR (FREQUENCY DIRECTION GRID).
C        ----------------------------------------------
C
 1000 CONTINUE
      READ(IU07) FR, DFIM, GOM, C,
     1           DELTH, DELTR, TH, COSTH, SINTH
C
C*    2. READ COMMON GRIDPAR (GENERAL GRID ORGANISATION).
C        ------------------------------------------------
C
 2000 CONTINUE
      READ(IU07) DELPHI, DELLAM, NLONRGG, SINPH, COSPH,
     1           IGL, IJS, IJL2, IJLS, IJL, IJLT
C
C*    3. READ COMMON MAP (LONG. AND LAT. INDICES OF GRID POINTS).
C        --------------------------------------------------------
C
 3000 CONTINUE
      READ(IU07) IXLG, KXLT, NX, NY, IPER,
     1           AMOWEP, AMOSOP, AMOEAP, AMONOP, XDELLA, XDELLO,
     2           ZDELLO, IRGG
C
C*    4. READ COMMON INDNL (NON-LINEAR INTERACTION).
C       --------------------------------------------
C
 4000 CONTINUE
      READ(IU07) IKP, IKP1, IKM, IKM1, K1W, K2W, K11W, K21W,
     1           AF11, FKLAP, FKLAP1, FKLAM, FKLAM1,
     2           ACL1, ACL2,  CL11, CL21, DAL1, DAL2, FRH
C
C*    5. READ COMMON COUPLE.
C        -------------------
C
 5000 CONTINUE
      READ(IU07) BETAMAX, ZALP, ALPHA, XKAPPA, XNLEV
C
C*    6. READ COMMON TABLE (INDICES OF OUTPUT POINTS).
C        ---------------------------------------------
C
 6000 CONTINUE
      READ(IU07) TAUT, DELTAUW, DELU, TAUHFT, DELUST, DELALP
C
C*    7. READ COMMON COUT (INDICES OF OUTPUT POINTS).
C        --------------------------------------------
C
 7000 CONTINUE
      READ(IU07)  NGOUT, IGAR, IJAR
C
C*    8. READ COMMON SHALLO (DEPTH AND SHALLOW WATER TABLES).
C        ----------------------------------------------------
C
 8000 CONTINUE
      READ(IU07)  DEPTH, DEPTHA, DEPTHD, TCGOND, TFAK, TSIHKD
C
C*    9. READ COMMON CURRENT.
C        --------------------
C
 9000 CONTINUE
      IF (IREFRA.EQ.2) READ(IU07)  U, V

      RETURN
      END
