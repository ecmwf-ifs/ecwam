      SUBROUTINE SPLITBL (FL, IJA, IJE, IU14)

C ----------------------------------------------------------------------
C
C**** *SPLITBL*   SAVES SECOND LATITUDE OF BLOCK IG
C****             FOR BLOCK IG-1 FOR COMPUTATION OF NEXT TIME STEP.
C
C     H. GUNTHER       GKSS/ECMWF         SEPTEMBER 1989
C
C*    PURPOSE.
C     --------
C
C       INTERMEDIATE STORAGE OF WAVE INFORMATION FOR USE
C       AT NEXT TIME STEP (SEE SUB FILLBL).
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *SPLITBL (FL, IJA, IJE, IU14)*
C          *FL*      REAL ARRAY   BLOCK OF SPECTRA
C          *IJA*     INTEGER      FIRST GRID POINT TO BE SAVED.
C          *IJE*     INTEGER      LAST  GRID POINT TO BE SAVED.
C          *IU14*    INTEGER      UNIT FOR OUTPUT OF SECOND LATITUDE.
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     METHOD.
C     -------
C
C       UNFORMATED WRITE TO SPECIFIED UNIT.
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
#include "commpp.h"
C
C ----------------------------------------------------------------------
C
      REAL,DIMENSION (NINF-1:NSUP,NANG,NFRE) :: FL
C
C ----------------------------------------------------------------------
C
C*    1. OUTPUT OF SECOND LATITUDE OF BLOCK FOR LAST
C        LATITUDE OF PREVIOUS BLOCK.
C        ---------------------------------------------
C
 1000 CONTINUE
      WRITE(IU14) (((FL(IJ,K,M),IJ=IJA,IJE),K=1,NANG),M=1,NFRE)

      RETURN
      END
