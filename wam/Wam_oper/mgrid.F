      SUBROUTINE MGRID (IA1)

C ----------------------------------------------------------------------
C
C**** *MGRID* - ROUTINE TO ARRANGE WAMODEL GRID.
C
C     H.GUNTHER            ECMWF       04/04/1990
C
C*    PURPOSE.
C     -------
C
C       TO ARRANGE WAMODEL GRID FOR A GIVEN AREA AND COMPUTE VARIOUS
C       MODEL CONSTANTS.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *MGRID (IA1)*
C          *IA1*     - TOPOGRAPHIC DATA OF PART
C
C     METHOD.
C     -------
C
C       THE NUMBER OF SEA POINTS PER LATITUDE IS COUNTED AND MODEL
C       BLOCKS OF MAXIMUM LENGTH OF NIBLO ARE CONSTRUCTED.
C
C     EXTERNALS.
C     ----------
C
C       *MBLOCK*    - SUB. TO GENERATE A BLOCK.
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
C ----------------------------------------------------------------------
C
      DIMENSION IA1(NGX,NGY)
C
C ----------------------------------------------------------------------
C
      DIMENSION IPP(NGY)
C
C          *IPP*   INTEGER   NUMBER OF SEA POINTS PER LATITUDE.
C
C ----------------------------------------------------------------------
C
C*    1. COUNT NUMBER OF SEA POINTS PER LATITUDE.
C        ----------------------------------------
C
 1000 CONTINUE
      DO 1001 K=1,NY
         IPP(K) = 0
         DO 1002 I=1,NLONRGG(K)
            IF (IA1(I,K).GT.-990) THEN
              IPP(K) = IPP(K) + 1
            ENDIF
 1002    CONTINUE
 1001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    2. MAKE BLOCKS.
C        ------------
C
 2000 CONTINUE
      IL = 0
      KA = 1
      DO 2001 K = 1,NY
         IL = IL + IPP(K)
         IF (IL.GT.NIBLO) THEN
            KE = K-1
            CALL MBLOCK (IA1, KA, KE, IPP)
            KA = KE-1
            IL = IPP(KA)+IPP(KE)+IPP(KE+1)
         ENDIF
 2001 CONTINUE
      CALL MBLOCK (IA1, KA, NY, IPP)

      RETURN
      END
