      SUBROUTINE MFREDIR (ML, KL)

C ----------------------------------------------------------------------
C
C**** *MFREDIR* - ROUTINE TO COMPUTE FREQUENCY DIRECTION CONSTANTS.
C
C     H. GUNTHER   ECMWF   2/4/90.
C     B. HANSEN    ECMWF   4/11/97. ACTIVATE THE ROTATION OF SPECTRA.
C
C*    PURPOSE.
C     --------
C
C       INITIATES THE FREQUENCY AND DIRECTION CONSTANTS WHICH ARE
C       SAVED IN COMMON FREDIR.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *MFREDIR (ML, KL)*
C          *ML*  - NUMBER OF FREQUENCIES
C          *KL*  - NUMBER OF DIRECTIONS
C
C     METHOD.
C     -------
C
C       STARTING FROM THE FIRST FREQUENCY THE NEXT  ARE INCREMENTED
C       BY THE FACTOR *CO* (SEE PARAMETER STATEMENT). THE DIRECTIONS
C       ARE EQUALLY DISTRIBUTED OVER THE CIRCLE STARTING WITH ZERO.
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
#include "comfred.h"
C
#include "comtest.h"
C
C ----------------------------------------------------------------------
C
#include "parcons.h"
C
C*    *PARAMETER*  SUB. CONSTANTS.
C
      PARAMETER (CO = 1.1)
C
C*     VARIABLE.   TYPE.     PURPOSE.
C      ---------   -------   --------
C      CO          REAL      FREQUENCY RATIO.
C
C ----------------------------------------------------------------------
C
C*    1. FREQUENCY DEPENDENT CONSTANTS.
C        ------------------------------
C
 1000 CONTINUE
C
C*    1.1 COMPUTE FREQUENCIES.
C         --------------------
C
      DO 1100 M=2,ML
         FR(M) = CO*FR(M-1)
 1100 CONTINUE
C
C*    1.2 COMPUTE DEEP WATER GROUP VELOCITIES.
C        ------------------------------------
C
      DO 1200 M=1,ML
         GOM(M) = G/(4.*PI*FR(M))
 1200 CONTINUE
C
C*    1.3 COMPUTE PHASE VELOCITY IN DEEP WATER.
C         -------------------------------------
C
      DO 1301 M = 1,ML
         C(M) = G/(ZPI*FR(M))
 1301 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    2. COMPUTATION OF DIRECTIONS, BANDWIDTH, SIN AND COS.
C        --------------------------------------------------
C
 2000 CONTINUE
      DELTH = ZPI/REAL(KL)
      DELTR = DELTH*R
      DO 2001 K=1,KL
CCCC        TH(K) = REAL(K-1)*DELTH
CCCC the previous line should be used if spectra should not be rotated.
CCCC the next line should be used if rotated spectra are used
        TH(K) = REAL(K-1)*DELTH + 0.5*DELTH
        COSTH(K) = COS(TH(K))
        SINTH(K) = SIN(TH(K))
 2001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    3. COMPUTATION FREQUENCY DIRECTION AREAS
C        -------------------------------------
C
 3000 CONTINUE
      CO1 = 0.5*(CO-1.)*DELTH
      DFIM(1)= CO1*FR(1)
      DO 3001 M=2,ML-1
         DFIM(M)=CO1 * (FR(M)+FR(M-1))
 3001 CONTINUE
      DFIM(ML)=CO1*FR(ML-1)
C
C ----------------------------------------------------------------------
C
C*    4. PRINTER PROTOCOL
C         ---------------
C
 4000 CONTINUE
      WRITE (IU06,'(''1FREQUENCY AND DIRECTION GRID'')')
      WRITE (IU06,'(''0NUMBER OF FREQUENCIES IS  ML = '',I3)') ML
      WRITE (IU06,'('' NUMBER OF DIRECTIONS  IS  KL = '',I3)') KL
      WRITE (IU06,'(''0MODEL FREQUENCIES IN HERTZ:'')')
      WRITE (IU06,'(1X,13F10.5)') (FR(M),M=1,ML)
      WRITE (IU06,'(''0MODEL FREQUENCY INTERVALLS TIMES DIRECTION'',
     1              '' INTERVALL IN HERTZ*RADIENS'')')
      WRITE (IU06,'(1X,13F10.5)') (DFIM(M),M=1,ML)
      WRITE (IU06,'(''0MODEL DEEP WATER GROUPVELOCITY IN M/S:'')')
      WRITE (IU06,'(1X,13F10.5)') (GOM(M),M=1,ML)
      WRITE (IU06,'(''0MODEL DEEP WATER PHASEVELOCITY IN M/S:'')')
      WRITE (IU06,'(1X,13F10.5)') (C(M),M=1,ML)
      WRITE (IU06,'(''0MODEL DIRECTIONS IN DEGREE'',
     1              '' (CLOCKWISE FROM NORTH):'')')
      WRITE (IU06,'(1X,13F10.5)') (TH(K)*DEG,K=1,KL)

      RETURN
      END
