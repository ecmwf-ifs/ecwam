      SUBROUTINE GRADI (IG, IREFRA, DDPHI, DDLAM, DUPHI, DULAM,
     1                  DVPHI, DVLAM)

C ----------------------------------------------------------------------
C
C**** *GRADI* - CALCULATES DEPTH AND CURRENT VELOCITY GRADIENTS.
C
C     K.P. HUBBERT              AUGUST   1988
C     H. GUNTHER    ECMWF/GKSS  DECEMBER 1990  MODIFIED FOR CYCLE_4.
C
C*    PURPOSE.
C     --------
C
C       CALCULATES DEPTH AND CURRENT VELOCITY GRADIENTS OF A BLOCK.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *GRADI (IG, IREFRA, DDPHI, DDLAM, DUPHI, DULAM,
C                      DVPHI, DVLAM)*
C          *IG*     - BLOCK NUMBER.
C          *IREFRA* - REFRACTION OPTION.
C          *DDPHI*  - LATITUDE DEPTH GRADIENT.
C          *DDLAM*  - LONGITUDE DEPTH GRADIENT.
C          *DUPHI*  - LATITUDE  U-COMPONENT GRADIENT.
C          *DULAM*  - LONGITUDE U-COMPONENT GRADIENT.
C          *DVPHI*  - LATITUDE  V-COMPONENT GRADIENT.
C          *DVLAM*  - LONGITUDE V-COMPONENT GRADIENT.
C
C     METHOD.
C     ------
C
C       CENTRAL DIFFERENCING FOR DEPTH AND CURRENT VELOCITY GRADIENTS.
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
#include "commap.h"
C
#include "comgrid.h"
C
#include "comshal.h"
C
#include "comubuf.h"
C
C ----------------------------------------------------------------------
C
      DIMENSION DDPHI(NIBLD), DDLAM(NIBLD), DUPHI(NIBLC), DULAM(NIBLC),
     1          DVPHI(NIBLC), DVLAM(NIBLC)
C
C ----------------------------------------------------------------------
C
C*    1. INITIALISE.
C        -----------
C
 1000 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    2. CALCULATE DEPTH GRADIENTS.
C        --------------------------
C
 2000 CONTINUE
      DO 2001 IJ=1,NIBLD
         DDPHI(IJ) = 0.0
         DDLAM(IJ) = 0.0
 2001 CONTINUE
      DO 2002 IJ=IJS(IG),IJL(IG)
         IPP = KLAT(IJ,2)
         IPM = KLAT(IJ,1)
         IF (IPP.GT.0 .AND. IPM.GT.0) THEN
            DDPHI(IJ) = (DEPTH(IPP,IG)-DEPTH(IPM,IG))/(2.*DELPHI)
         ENDIF
         ILP = KLON(IJ,2)
         ILM = KLON(IJ,1)
         KX  = KXLT(IJ,IG)
         IF (ILP.GT.0 .AND. ILM.GT.0) THEN
            DDLAM(IJ) = (DEPTH(ILP,IG)-DEPTH(ILM,IG))/(2.*DELLAM(KX))
         ENDIF
 2002 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    3. CALCULATE CURRENT VELOCITY GRADIENTS.
C        -------------------------------------
C
 3000 CONTINUE
      IF (IREFRA.EQ.2) THEN
         DO 3001 IJ=1,NIBLC
            DUPHI(IJ) = 0.0
            DULAM(IJ) = 0.0
            DVPHI(IJ) = 0.0
            DVLAM(IJ) = 0.0
 3001    CONTINUE
         DO 3002 IJ=IJS(IG),IJL(IG)
            IPP = KLAT(IJ,2)
            IPM = KLAT(IJ,1)
            IF (IPP.GT.0 .AND. IPM.GT.0) THEN
               DUPHI(IJ) = (U(IPP,IG)-U(IPM,IG))/(2.*DELPHI)
               DVPHI(IJ) = (V(IPP,IG)-V(IPM,IG))/(2.*DELPHI)
            ENDIF
            ILP = KLON(IJ,2)
            ILM = KLON(IJ,1)
            KX  = KXLT(IJ,IG)
            IF (ILP.GT.0 .AND. ILM.GT.0) THEN
               DULAM(IJ) = (U(ILP,IG)-U(ILM,IG))/(2.*DELLAM(KX))
               DVLAM(IJ) = (V(ILP,IG)-V(ILM,IG))/(2.*DELLAM(KX))
            ENDIF
 3002    CONTINUE
      ENDIF

      RETURN
      END
