      SUBROUTINE MBOXB (NBOUN, AMOWEB, AMOSOB, AMOEAB,
     1                  AMONOB, BLATB, BLNGB)

C ----------------------------------------------------------------------
C
C**** *MBOXB* - MAKE BOX OF FINE GRID IN COARSE GRID.
C
C     R. PORTZ     MPI         15/01/1991
C
C*    PURPOSE.
C     -------
C
C       COMPUTE LATITUDES AND LONGITUDES OF COARSE GRID POINTS
C       AT NEST BOUNDARY.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *MBOXB (NBOUN, AMOWEB, AMOSOB, AMOEAB, AMONOB,
C                      BLATB, BLNGB)*
C          *NBOUN*  - NUMBER OF BOUNDARY POINTS.
C          *AMOWEB* - WESTERN  LONGITUDE OF COARSE GRID.
C          *AMOSOB* - SOUTHERN LATITUDE  OF COARSE GRID.
C          *AMOEAB* - EASTERN  LONGITUDE OF COARSE GRID.
C          *AMONOB* - NORTHERN LATITUDE  OF COARSE GRID.
C          *BLATB*  - LATITUDE  OF BOUNDARY POINTS.
C          *BLNGB*  - LONGITUDE OF BOUNDARY POINTS.
C
C     METHOD.
C     -------
C
C       NONE.
C
C     EXTERNALS.
C     ----------
C
C       *ABORT1*     - TERMINATES PROCESSING.
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
#include "commap.h"
C
#include "comtest.h"
C
C ----------------------------------------------------------------------
C
      DIMENSION BLNGB(NBMAX), BLATB(NBMAX)
C
C ----------------------------------------------------------------------
C
C*    1. COMPUTED THE SQUARE BOX
C     --------------------------
C
 1000 CONTINUE
      NLNGB = NINT((AMOEAB - AMOWEB) / XDELLO) + 1
C
      NLATB = (NINT((AMONOB - AMOSOB) / XDELLA) - 1) * 2
C
      NBOUN = (NLNGB * 2) + NLATB
C
C     DIMENSION CHECK: NBMAX > NBOUN
C
      IF (NBMAX .LT. NBOUN) THEN
         WRITE(IU06,*) ' **********************************'
         WRITE(IU06,*) ' *                                *'
         WRITE(IU06,*) ' *   FATAL ERROR IN SUB. MBOX     *'
         WRITE(IU06,*) ' *   ========================     *'
         WRITE(IU06,*) ' *  NUMBER OF BOUNDARY POINTS     *'
         WRITE(IU06,*) ' *  EXCEEDS DIMENSION.            *'
         WRITE(IU06,*) ' *  DIMENSION IS        NBMAX = ', NBMAX
         WRITE(IU06,*) ' *  NUMBER OF POINTS IS NBOUN = ', NBOUN
         WRITE(IU06,*) ' *                                *'
         WRITE(IU06,*) ' **********************************'
         CALL ABORT1
      ENDIF
C
C     COMPUTES THE BOUNDARY POINTS FOR THE FIRST AND THE LAST
C     LATITUDE
C
      K = NLATB + NLNGB
      DO 1001 I = 1, NLNGB
         BLATB(I) = AMOSOB
         BLNGB(I) = AMOWEB + REAL(I-1) * XDELLO
         BLATB(K+I) = AMONOB
         BLNGB(K+I) = BLNGB(I)
 1001 CONTINUE
C
C     COMPUTED THE EAST AND THE WEST BOUNDARY POINT FOR
C     EACH LATITUDE
C
      K = 0
      DO 1002 I = 2,NLATB,2
         K = K + 1
         BLATB(NLNGB+I-1) = AMOSOB + REAL(K) * XDELLA
         BLNGB(NLNGB+I-1) = AMOWEB
         BLATB(NLNGB+I)   = BLATB(NLNGB+I-1)
         BLNGB(NLNGB+I)   = AMOEAB
 1002 CONTINUE
C
      RETURN
      END
