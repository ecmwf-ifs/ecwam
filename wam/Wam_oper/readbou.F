      SUBROUTINE READBOU (IU09, IU10, IU06)

C ----------------------------------------------------------------------
C
C**** *READBOU* - READ COARSE AND FINE GRID BOUNDARY INFORMATION FILE.
C
C     R. PORTZ     MPI          JANUARY 1991
C
C*    PURPOSE.
C     --------
C
C       READ COMMON CBOUND AND FBOUND AS WRITTEN BY PREPROC.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *READBOU (IU09, IU10, IU06)*
C          *IU09*    - INPUT  UNIT OF COMMON CBOUND.
C          *IU10*    - INPUT  UNIT OF COMMON FBOUND.
C          *IU06*    - PRINTER OUTPUT UNIT.
C
C     METHOD.
C     -------
C
C       SEQUENCIAL UNFORMATED WRITE TO UNIT.
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
#include "comcbou.h"
C
#include "comfbou.h"
C
C ----------------------------------------------------------------------
C
 1000 CONTINUE
C
      IF (IBOUNC.EQ.1) THEN
C
C     READ COMMON BOUNC
C
         READ (IU09,ERR=2000) NBOUNC
         READ (IU09,ERR=2000) (IGARC(I),I=1,NBOUNC),
     +                        (IJARC(I),I=1,NBOUNC)
C
C        CHECK THE DIMENSION
C
         IF (NMAXC.LT.NBOUNC) THEN
            WRITE(IU06,*) '*********************************'
            WRITE(IU06,*) '*  FATAL ERROR IN SUB. READBOU  *'
            WRITE(IU06,*) '*  ===========================  *'
            WRITE(IU06,*) '*                               *'
            WRITE(IU06,*) '*   NMAXC < NBOUNC              *'
            WRITE(IU06,*) '*                               *'
            WRITE(IU06,*) '*   NMAXC: ', NMAXC, ' NBOUNC: ', NBOUNC
            WRITE(IU06,*) '*                               *'
            WRITE(IU06,*) '*********************************'
            CALL ABORT1
         ENDIF
      ENDIF
C
      IF (IBOUNF.EQ.1) THEN
C
C     READ COMMON BOUNF
C
         READ (IU10,ERR=2000) NBOUNF
         READ (IU10,ERR=2000) (IGARF(I),I=1,NBOUNF),
     +                        (IJARF(I),I=1,NBOUNF),
     +                        (IBFL(I),I=1,NBOUNF),
     +                        (IBFR(I),I=1,NBOUNF),
     +                        (BFW(I),I=1,NBOUNF)
C
C        CHECK THE DIMENSION
C
         IF (NMAXF.LT.NBOUNF) THEN
            WRITE(IU06,*) '*********************************'
            WRITE(IU06,*) '*  FATAL ERROR IN SUB. READBOU  *'
            WRITE(IU06,*) '*  ===========================  *'
            WRITE(IU06,*) '*                               *'
            WRITE(IU06,*) '*   NMAXF < NBOUNF              *'
            WRITE(IU06,*) '*                               *'
            WRITE(IU06,*) '*   NMAXF: ', NMAXF, ' NBOUNF: ', NBOUNF
            WRITE(IU06,*) '*                               *'
            WRITE(IU06,*) '*********************************'
            CALL ABORT1
         ENDIF
      ENDIF
      RETURN
 2000 CONTINUE
      WRITE(IU06,*) '****************************************'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '*    FATAL ERROR IN SUB. READBOU       *'
      WRITE(IU06,*) '*    ===========================       *'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '*        READ ERROR                    *'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '*        PROGRAM ABORTS                *'
      WRITE(IU06,*) '****************************************'
      CALL ABORT1
      END
