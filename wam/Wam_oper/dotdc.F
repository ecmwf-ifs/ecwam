      SUBROUTINE DOTDC (IJS, IJL, IG, IGL, ISHALLO, IREFRA, SL, FCONST)

C ----------------------------------------------------------------------
C
C**** *DOTDC* - SIGMA DOT FOR CURRENT REFRACTION.
C
C     H. GUNTHER   GKSS/ECMWF   17/02/91
C     J. BIDLOT    ECMWF FEBRUARY 1997 : MESSAGE PASSING
C
C*    PURPOSE.
C     --------
C
C       READ (FROM IU16) FOR PROPAGATION ..
C           - DEPTH AND CURRENT PART OF THETA DOT.
C           - SIGMA DOT TERM.
C       SCATTER SIGMA/ SINH (2*K*D) TABLE.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *DOTDC (IJS, IJL, IG, IGL, ISHALLO, IREFRA, SL, FCONST)*
C          *IJS*     - INDEX OF FIRST GRID POINT.
C          *IJL*     - INDEX OF LAST  GRID POINT.
C          *IG*      - BLOCK NUMBER.
C          *IGL*     - NUMBER OF BLOCKS.
C          *ISHALLO* - SHALLOW WATER OPTION.
C          *IREFRA*  - REFRACTION OPTION.
C          *SL*      - SIGMA DOT
C          *FCONST * - SIGMA/SINH(2KD)
C
C     METHOD.
C     -------
C
C       SEQUENTIAL READ FROM UNIT.
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
#include "commesspass.h"
C
#include "commpp.h"
C
#include "comrefd.h"
C
#include "comshal.h"
C
#include "comunit.h"
C
C ----------------------------------------------------------------------
C
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: SL
      REAL,DIMENSION(NIBLO,NFRE):: FCONST
C
C -----------------------------------------------------------------------
C
      CHARACTER*30 FILENAME

C*    1. READ COMMON REFDOT.
C        -------------------
C
 1000 CONTINUE
      IF (IG.EQ.1) REWIND IU16
      IF (IGL.NE.1) THEN
         IF (ISHALLO.NE.1) THEN
            IF (IREFRA.EQ.1) THEN
               READ (IU16) THDD
            ELSE
               READ (IU16) THDD, THDC, SL
            ENDIF
         ELSE
            READ (IU16) THDC, ((SL(IJ,K,NFRE),IJ=1,NIBLC),K=1,NANG)
         ENDIF
      ELSE
        IF(LMESSPASS) THEN
           FILENAME='fort.16.%p'
           CALL EXPAND_STRING(IRANK,NPROC,0,0,FILENAME,1)
           OPEN(IU16,FILE=FILENAME,STATUS='UNKNOWN')

           READ(IU16) SL
           CLOSE(IU16)

        ELSE
           IF (ISHALLO.NE.1) THEN
              IF (IREFRA.EQ.2) THEN
                 READ (IU16) SL
              ENDIF
           ELSE
              READ (IU16) ((SL(IJ,K,NFRE),IJ=1,NIBLC),K=1,NANG)
           ENDIF
         ENDIF
      ENDIF
C
C*    2. IF DEEP WATER RETURN.
C        ---------------------
C
 2000 CONTINUE
      IF (ISHALLO.EQ.1) RETURN
C
C*    3. GATHER SIGMA /SINH(2KD) FROM TABLE.
C        -----------------------------------
C
 3000 CONTINUE
      DO 3001 M=1,NFRE
         DO 3002 IJ=IJS,IJL
            FCONST(IJ,M) = TSIHKD(INDEP(IJ),M)
 3002     CONTINUE
 3001 CONTINUE

      RETURN
      END
