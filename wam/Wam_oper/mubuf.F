      SUBROUTINE MUBUF (IA1, IG, IU08, IU18, IFORM)

C ----------------------------------------------------------------------
C
C**** *MUBUF* - ROUTINE TO ARRANGE COMMON UBUF FOR ONE BLOCK.
C
C     H.GUNTHER            ECMWF       04/04/1990
C
C*    PURPOSE.
C     -------
C
C       TO ARRANGE NEIGHBOUR GRID POINT INDICES FOR A BLOCK
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *MUBUF (IA1, IG, IU08, IU18, IFORM)*
C          *IA1*     - TOPOGRAPHIC DATA.
C          *IG*      - BLOCK NUMBER.
C          *IU08*    - LOGICAL UNIT FOR OUTPUT OF GRID BLOCKING
C                      COMMON UBUF (UNFORMATED)
C          *IU18*    - LOGICAL UNIT FOR OUTPUT OF GRID BLOCKING
C                      COMMON UBUF (FORMATED)
C          *IFORM*   - OUTPUT FORMAT OPTION = 1 UNFORMATED
C                                           = 2 FORMATED
C                                           OTHERWISE BOTH
C
C     METHOD.
C     -------
C
C       THE INDICES OF THE NEXT POINTS ON LAT. AND LONG. ARE
C       COMPUTED. ZERO INDICATES A LAND POINT IS NEIGHBOUR.
C       THE FINAL COMMON UBUF IS WRITTEN OUT.
C
C     EXTERNALS.
C     ----------
C
C       *OUTUBUF*   - WRITE OUT COMMON UBUF.
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
#include "comgrid.h"
C
#include "commap.h"
C
#include "comubuf.h"
C
C ----------------------------------------------------------------------
C
      DIMENSION IA1(NGX, NGY)
C
C ----------------------------------------------------------------------
C
C*    1. INITIALISE ARRAYS.
C        ------------------
C
 1000 CONTINUE
      DO 1001 IJ=1,NIBLO
         KLAT(IJ,1) = 0
         KLAT(IJ,2) = 0
         KLON(IJ,1) = 0
         KLON(IJ,2) = 0
 1001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    2. COMPUTE INDICES OF NEIGHBOUR SEA POINTS.
C        ----------------------------------------
C
 2000 CONTINUE
C
C*    2.1 LONGITUDE NEIGHBOURS.
C         ---------------------
C
      DO 2100 IP = 1,IJLT(IG)
         I = IXLG(IP,IG)
         K = KXLT(IP,IG)
         IF (I.GT.1) THEN
           IF (IA1(I-1,K).GT.-990) KLON(IP,1) = IP-1
         ELSE
            IF (IPER.EQ.1 .AND. IA1(NLONRGG(K),K).GT.-990) THEN
               KLON(IP,1) = IP
               DO 2101 IH=2,NLONRGG(K)
                  IF (IA1(IH,K).GT.-990) KLON(IP,1) = KLON(IP,1)+1
 2101          CONTINUE
            ENDIF
         ENDIF
         IF (I.LT.NLONRGG(K)) THEN
           IF (IA1(I+1,K).GT.-990) KLON(IP,2) = IP+1
         ELSE
           IF (IPER.EQ.1 .AND. IA1(1,K).GT.-990) THEN
              KLON(IP,2) = IP
              DO 2102 IH=NLONRGG(K)-1,1,-1
                 IF (IA1(IH,K).GT.-990) KLON(IP,2) = KLON(IP,2)-1
 2102         CONTINUE
            ENDIF
         ENDIF
 2100 CONTINUE
C
C*    2.2 LATITUDE NEIGHBOURS.
C         --------------------
C
      DO 2200 IP = 1,IJLT(IG)
         I = IXLG(IP,IG)
         K = KXLT(IP,IG)
         IF (K.GT.1) THEN
            XMIN = FLOAT(I-1)*ZDELLO(K)/ZDELLO(K-1)
            IMIN = NINT(XMIN) + 1
            IF (IA1(IMIN,K-1).GT.-990) THEN
               DO 2201 IH = IP,1,-1
                  IF (IXLG(IH,IG).EQ.IMIN .AND.
     1                KXLT(IH,IG).EQ.K-1) GOTO 2202
 2201          CONTINUE
 2202          CONTINUE
               KLAT(IP,1) = IH
            ENDIF
         ENDIF
         IF (K.LT.NY) THEN
            XPLUS = FLOAT(I-1)*ZDELLO(K)/ZDELLO(K+1)
            IPLUS = NINT(XPLUS) + 1
            IF (IA1(IPLUS,K+1).GT.-990) THEN
               DO 2203 IH = IP,IJLT(IG)
                  IF (IXLG(IH,IG).EQ.IPLUS .AND.
     1                KXLT(IH,IG).EQ.K+1) GOTO 2204
 2203          CONTINUE
 2204          CONTINUE
               KLAT(IP,2) = IH
           ENDIF
         ENDIF
 2200 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    3. OUTPUT OF COMMON UBUF FOR ONE BLOCK.
C        ------------------------------------
C
 3000 CONTINUE
      CALL OUTUBUF (IU08, IU18, IFORM)

      RETURN
      END
