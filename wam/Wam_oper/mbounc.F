      SUBROUTINE MBOUNC (IU09, IU19, IFORM)

C ----------------------------------------------------------------------
C
C**** *MBOUNC* - MAKE COARSE GRID BOUNDARY.
C
C     R. PORTZ     MPI         15/01/1991
C
C*    PURPOSE.
C     -------
C
C       COMPUTE ALL INFORMATION FOR COARSE GRID BOUNDARY VALUE
C       OUTPUT (COMMON CBOUND).
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *MBOUNC (IU09, IU19, IFORM)*
C          *IU09*   - LOGICAL UNIT FOR  UNFORMATED WRITE.
C          *IU19*   - LOGICAL UNIT FOR    FORMATED WRITE.
C          *IFORM*  - FORMAT OPTION  = 1  UNFORMATED WRITE.
C                                    = 2  FORMATED WRITE.
C                                      OTHERWISE BOTH.
C
C     METHOD.
C     -------
C
C       NONE.
C
C     EXTERNALS.
C     ----------
C
C       *FINDB*     - FIND BLOCK AND SEA POINT NUMBERS.
C       *MBOXB*     - FIND LAT AND LONG OF BOUNDARY POINTS.
C       *PACKI*     - PACKS AN INTEGER ARRAY.
C       *PACKC*     - PACKS A  REAL    ARRAY.
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
#include "commap.h"
C
#include "comgrid.h"
C
#include "comtest.h"
C
C ----------------------------------------------------------------------
C
  998 FORMAT(10I8)
  999 FORMAT(5E16.7)
C
C ----------------------------------------------------------------------
C
C*    1. INITIAL.
C        --------
C
 1000 CONTINUE
      DPHIAC = XDELLA
      DLAMAC = XDELLO
C
      DO 1001 I = 1,NBMAX
         IJARC(I) = 0
         IGARC(I) = 0
         BLATC(I) = 0.
         BLNGC(I) = 0.
 1001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    2. COMPUTED THE SQUARE BOX.
C        ------------------------
C
 2000 CONTINUE
      CALL MBOXB (NBOUNC, AMOWEC, AMOSOC, AMOEAC, AMONOC, BLATC, BLNGC)
C
C ----------------------------------------------------------------------
C
C*    3. SEARCH BLOCK NUMBER AND SEA POINT NUMBER.
C        -----------------------------------------
C
 3000 CONTINUE
      CALL FINDB (NBMAX, NBOUNC, BLATC, BLNGC, IGARC, IJARC)
C
C ----------------------------------------------------------------------
C
C*    4. PACKED ALL ARRAYS.
C        -------------------
C
      CALL PACKR (NBOUNC, NBOUNEW, NBMAX, IGARC, BLATC)
      CALL PACKR (NBOUNC, NBOUNEW, NBMAX, IGARC, BLNGC)
      CALL PACKI (NBOUNC, NBOUNEW, NBMAX, IGARC, IJARC)
      CALL PACKI (NBOUNC, NBOUNEW, NBMAX, IGARC, IGARC)
      NBOUNC = NBOUNEW
C
C ----------------------------------------------------------------------
C
C*    5. PRINTER PROTOCOL.
C        -----------------
C
 5000 CONTINUE
      WRITE (IU06,'(''1BOUNDARY OUTPUT POINTS OF THIS GRID:'')')
      WRITE (IU06,*) '===================================='
      WRITE (IU06,*) '  NUMBER OF BOUNDARY POINTS IS NBOUNC = ', NBOUNC
      WRITE (IU06,'(4X,''     |-----INPUT-----|'',
     1               ''-POINT INDEX--|'')')
      WRITE (IU06,'(4X,''  NO.    LAT.   LONG.  BLOCK.'',
     1              ''  POINT.'')')
      DO 5001 IO=1,NBOUNC
         WRITE (IU06,'(4X,I5,2F8.2,2I8)')
     1     IO, BLATC(IO), BLNGC(IO), IGARC(IO), IJARC(IO)
 5001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    6. WRITE COMMON CBOUND.
C        --------------------
C
 6000 CONTINUE
      IF (IFORM.NE.2) THEN
         WRITE(IU09) NBOUNC
         WRITE(IU09) (IGARC(I),I=1,NBOUNC), (IJARC(I),I=1,NBOUNC)
         WRITE(IU09) DLAMAC, DPHIAC, AMOSOC, AMONOC, AMOEAC, AMOWEC,
     1                (BLNGC(I),I=1,NBOUNC), (BLATC(I),I=1,NBOUNC)
      ENDIF
      IF (IFORM.NE.1) THEN
         WRITE(IU19,998) NBOUNC
         WRITE(IU19,998) (IGARC(I),I=1,NBOUNC), (IJARC(I),I=1,NBOUNC)
         WRITE(IU19,999) DLAMAC,DPHIAC, AMOSOC, AMONOC, AMOEAC, AMOWEC,
     1                    (BLNGC(I),I=1,NBOUNC), (BLATC(I),I=1,NBOUNC)
      ENDIF
C
C ----------------------------------------------------------------------
C
      RETURN
      END
