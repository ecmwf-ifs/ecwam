C ======================================================================
C
      SUBROUTINE sorti (iarray, index, n)
C
C ======================================================================
C
C**** *SORTI* - SORTS AN INTEGER ARRAY.
C
C     HEINZ GUNTHER    ECMWF/GKSS    JANUARY 1991
C
C     PURPOSE
C     -------
C
C        SORT AN INTEGER ARRAY BY AN INDEX ARRAY.
C
C     METHOD
C     ------
C
C       THE DATA IN THE INPUT ARRAY ARE SORTED AS DEFINED IN THE INPUT
C       INDEX ARRAY.
C       THIS ARRAY MAY BE GENEGATED BY SUB. SORTINI OR SORTINR.
C
C**   INTERFACE
C     ---------
C
C        *CALL* *SORTI (IARRAY, INDEX, N)
C
C          IARRAY  INTEGER   INPUT/ OUTPUT DATA ARRAY.
C          INDEX   INTEGER   INPUT INDEX ARRAY.
C          N       INTEGER   NUMBER OF DATA TO BE SORTED.
C
C     EXTERNALS
C     ---------
C
C          NONE
C
C     REFERENCES
C     ----------
C
C          NONE
C
C ----------------------------------------------------------------------
C
      DIMENSION IARRAY(N), INDEX(N), IWORK(N)
C
C ----------------------------------------------------------------------
C
C*   1. SORT DATA INTO WORK ARRAY.
C       --------------------------
C
 1000 CONTINUE
      DO 1001 I=1,N
         IWORK(I) = IARRAY(INDEX(I))
 1001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*   2. COPY BACK TO DATA ARRAY.
C       ------------------------
C
 2000 CONTINUE
      DO 2001 I=1,N
         IARRAY(I) = IWORK(I)
 2001 CONTINUE
      RETURN
      END
