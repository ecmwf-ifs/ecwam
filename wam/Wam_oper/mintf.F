      SUBROUTINE MINTF

C ----------------------------------------------------------------------
C
C**** *MINTF* - MAKE INTERPOLATION TABLES FOR BOUNDARY INPUT.
C
C     R. PORTZ     MPI         15/01/1991
C     H. GUNTHER   GKSS/ECMWF  15/01/1991
C
C*    PURPOSE.
C     -------
C
C       GENERATE SPACE INTERPOLATION TABLES USED FOR BOUNDARY
C       VALUE INPUT INTO A FINE GRID MODEL.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *MINT*
C
C     METHOD.
C     -------
C
C       NONE.
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
#include "parameter_preproc.h"
C
#include "comcpbo.h"
C
#include "comfpbo.h"
C
#include "commap.h"
C
C ----------------------------------------------------------------------
C
C*    1. RATIOS OF GRID INCREMENTS.
C        --------------------------
C
 1000 CONTINUE
      IDELLA = NINT(DPHIAC/XDELLA)
      IDELLO = NINT(DLAMAC/XDELLO)
C
C*    2. SOUTHERN MOST LATITUDE OF FINE GRID.
C        ------------------------------------
C
 2000 CONTINUE
      NI    = IDELLO - 1
      PHI   = AMOSOP
C
C*    2.1 LOOP OVER COARSE GRID POINTS.
C         -----------------------------
C
      DO 2101 I = 1, NX, IDELLO
C
C*    2.2 INTERPOLATION WEIGHT FOR INTERMEDIATE POINTS.
C         ---------------------------------------------
C
         IF (I.NE.NX) THEN
            DO 2201 N = 1, NI
               BFW (I+N) = REAL(N) / REAL(IDELLO)
 2201       CONTINUE
         ENDIF
C
C*    2.3 INDICES OF COARSE GRID OUTPUT POINTS.
C         -------------------------------------
C
         XLAMDA = AMOWEP + REAL(I-1) * XDELLO
         DO 2301 M=1,NBOUNC
            IF (BLATC(M).EQ.PHI .AND. BLNGC(M).EQ.XLAMDA) THEN
               IBFL(I) = M
               IBFR(I) = M
               DO 2302 N = 1, NI
                  IF (I.NE.1) THEN
                     IBFR(I-N) = M
                  ENDIF
                  IF (I.NE.NX) THEN
                     IBFL(I+N) = M
                  ENDIF
 2302          CONTINUE
               GOTO 2101
            ENDIF
 2301    CONTINUE
 2101 CONTINUE
C
C*    3. NORTHERN MOST LATITUDE OF FINE GRID.
C        ------------------------------------
C
 3000 CONTINUE
      PHI = AMONOP
      IS  = NX + 2*(NY-2) + 1
      IE  = 2*(NX+NY-2)
C
C*    3.1 LOOP OVER COARSE GRID POINTS.
C         -----------------------------
C
      DO 3101 I = IS, IE, IDELLO
C
C*    3.2 INTERPOLATION WEIGHT FOR INTERMEDIATE POINTS.
C         ---------------------------------------------
C
         IF (I.NE.IE) THEN
            DO 3201 N = 1, NI
               BFW (I+N) = REAL(N) / REAL(IDELLO)
 3201       CONTINUE
         ENDIF
C
C*    3.3 INDICES OF COARSE GRID OUTPUT POINTS.
C         -------------------------------------
C
         XLAMDA = AMOWEP + (I-IS) * XDELLO
         DO 3301 M=1,NBOUNC
            IF (BLATC(M).EQ.PHI .AND. BLNGC(M).EQ.XLAMDA) THEN
               IBFL(I) = M
               IBFR(I) = M
               DO 3302 N = 1, NI
                  IF (I.NE.IS) THEN
                     IBFR(I-N) = M
                  ENDIF
                  IF (I.NE.IE) THEN
                     IBFL(I+N) = M
                  ENDIF
 3302          CONTINUE
               GOTO 3101
            ENDIF
 3301    CONTINUE
 3101 CONTINUE
C
C*    4. WESTERN MOST LONGITUDE OF FINE GRID.
C        ------------------------------------
C
 4000 CONTINUE
      XLAMDA = AMOWEP
      NI = IDELLA - 1
      NSTEP = 2 * IDELLA
      IE = NX + 2*(NY-2) - 1
      IS = NX + 2*NI + 1
      K = 1
C
C*    4.1 WEIGHTS AND LEFT INDICES FOR FIRST COARSE GRID SECTION.
C        --------------------------------------------------------
C
      DO 4101 N = 1, NI
         IBFL(NX-1+2*N) = IBFL(1)
         BFW (NX-1+2*N) = REAL(N) / REAL(IDELLA)
 4101 CONTINUE
C
C*    4.2 LOOP OVER COARSE GRID POINTS.
C         -----------------------------
C
      DO 4201 I = IS, IE, NSTEP
C
C*    4.3 INTERPOLATION WEIGHT FOR INTERMEDIATE POINTS.
C         ---------------------------------------------
C
         DO 4301 N = 1, NI
            BFW (I+2*N) = REAL(N) / REAL(IDELLA)
 4301    CONTINUE
C
C*    4.4 INDICES OF COARSE GRID OUTPUT POINTS.
C         -------------------------------------
C
         K = K + 1
         PHI = AMOSOP + (K-1) * DPHIAC
         DO 4401 M=1,NBOUNC
            IF (BLNGC(M).EQ.XLAMDA .AND. BLATC(M).EQ.PHI) THEN
               IBFL(I) = M
               IBFR(I) = M
               DO 4402 N = 1, NI
                  IBFL(I+2*N) = M
                  IBFR(I-2*N) = M
 4402          CONTINUE
               GOTO 4201
            ENDIF
 4401    CONTINUE
 4201 CONTINUE
C
C*    4.5 RIGHT INDICES FOR LAST COARSE GRID SECTION.
C         -------------------------------------------
C
      IF (IBFR(K) .NE. 0) THEN
         K = NX + (2* (NY-2)) + 1
         DO 4501 N = 1, NI
            IBFR(K-2*N) = IBFR(K)
 4501    CONTINUE
      ENDIF
C
C*    5. EASTERN MOST LONGITUDE OF FINE GRID.
C        ------------------------------------
C
 5000 CONTINUE
      XLAMDA = AMOEAP
      IS = NX + 2*NI + 2
      IE = NX + 2*(NY-2)
      K = 1
C
C*    5.1 WEIGHTS AND LEFT INDICES FOR FIRST COARSE GRID SECTION.
C        --------------------------------------------------------
C
      DO 5101 N = 1, NI
         IBFL(NX+2*N) = IBFL(NX)
         BFW (NX+2*N) = REAL(N) / REAL(IDELLA)
 5101 CONTINUE
C
C*    5.2 LOOP OVER COARSE GRID POINTS.
C         -----------------------------
C
      DO 5201 I = IS, IE, NSTEP
C
C*    5.3 INTERPOLATION WEIGHT FOR INTERMEDIATE POINTS.
C         ---------------------------------------------
C
         DO 5301 N = 1, NI
            BFW (I+2*N) = REAL(N) / REAL(IDELLA)
 5301    CONTINUE
C
C*    5.4 INDICES OF COARSE GRID OUTPUT POINTS.
C         -------------------------------------
C
         K = K + 1
         PHI = AMOSOP + (K-1) * DPHIAC
         DO 5401 M=1,NBOUNC
            IF (BLNGC(M).EQ.XLAMDA .AND. BLATC(M).EQ.PHI) THEN
               IBFL(I) = M
               IBFR(I) = M
               DO 5402 N = 1, NI
                  IBFL(I+2*N) = M
                  IBFR(I-2*N) = M
 5402          CONTINUE
               GOTO 5201
            ENDIF
 5401    CONTINUE
 5201 CONTINUE
C
C*    5.5 RIGHT INDICES FOR LAST COARSE GRID SECTION.
C         -------------------------------------------
C
      IF (IBFR(K) .NE. 0) THEN
         K = 2*(NX + NY - 2)
         M = NX + 2*NY - 2
         DO 5501 N = 1, NI
            IBFR(M-2*N) = IBFR(K)
 5501    CONTINUE
      ENDIF

      RETURN
      END
