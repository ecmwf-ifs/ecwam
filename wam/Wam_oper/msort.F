C--------------------------------------------------------------------
C
      SUBROUTINE MSORT (IU06, NUM, KICE, IG, IJ)
C
C-----------------------------------------------------------------
C
C**** *MSORT* - SORT ICE POINTS
C
C     R. PORTZ         MPI HAMBURG       OCTOBER 1992
C
C
C     PURPOSE.
C     --------
C             SORT ICE-POINTS ACCORDING TO BLOCK AND INDEX
C             IN A BLOCK AND WRITE INFO FILE FOR WAMODEL (TAPE 8)
C
C
C**   INTERFACE.
C     ---------
C             *CALL MSORT*
C
C*     VARIABLE.   TYPE.     PURPOSE.
C      ---------   -------   --------
C      *IU06*      INTEGER   PRINTER OUTPUT UNIT.
C      *NUM*       INTEGER   NUMBER OF ICE POINTS IN
C                            MODEL GRID.
C      *NICE*      INTEGER   NUMBER OF ICE POINTS
C      *IG*        INTEGER   BLOCK NUMBER (ICE POINT)
C      *IJ*        INTEGER   NUMBER OF INDEX IN A BLOCK (ICE POINT)
C
C     METHOD.
C     -------
C             NONE
C
C    EXTERNAL.
C    ---------
C           *MYINDEX
C
C---------------------------------------------------------------------
C
C     THIS IS NOW A DUMMY ROUTINE
      RETURN
      END
