C ======================================================================
C
      SUBROUTINE sortinr (array, index, n)
C
C ======================================================================
C
C**** *SORTINR* - SORTS AN INDEX ARRAY FORM A REAL DATA ARRAY.
C
C     HEINZ GUNTHER    ECMWF/GKSS    JANUARY 1991
C
C     PURPOSE
C     -------
C
C        FIND THE INDEX ARRAY WHICH SORTS THE DATA.
C
C     METHOD
C     ------
C
C       THE INPUT ARRAY IS SCANNED AND AN INDEX ARRAY IS GENERATED,
C       WHICH CONTAINS THE INDEX OF THE INPUT ARRAY TO SORT THE DATA
C       IN INCREASING ORDER. THE INPUT ARRAY IS NOT CHANGED.
C
C**   INTERFACE
C     ---------
C
C        *CALL* *SORTINR (ARRAY, INDEX, N)
C
C          ARRAY   REAL      INPUT DATA TO BE SCANNED
C          INDEX   INTEGER   OUTPUT ARRAY OF INDICES.
C          N       INTEGER   NUMBER OF DATA TO BE SCANNED
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
      DIMENSION ARRAY(N), INDEX(N)
C
C ----------------------------------------------------------------------
C
C*   1. INITIAL INDEX ARRAY.
C       -------------------
C
 1000 CONTINUE
      DO 1001 I=1,N
         INDEX(I) = I
 1001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*   2. CHECK FIRST TWO DATA.
C       ---------------------
C
 2000 CONTINUE
      IF (N.LT.2) RETURN
      IF (ARRAY(2).LT.ARRAY(1)) THEN
         INDEX(2)=1
         INDEX(1)=2
      ENDIF
      IF (N.EQ.2) RETURN
C
C ----------------------------------------------------------------------
C
C*   3. DO REST OF DATA.
C       ----------------
C
 3000 CONTINUE
      DO 3001 I=3,N
         IH = INDEX(I)
         IHM = INDEX(I-1)
         IF (ARRAY(IH).LT.ARRAY(IHM)) THEN
            DO 3002 J=I-1,2,-1
               INDEX(J+1) = INDEX(J)
               IF (ARRAY(INDEX(J-1)).LE.ARRAY(IH)) THEN
                  INDEX(J) = IH
                  GOTO 3001
                ENDIF
 3002       CONTINUE
            INDEX(2) = INDEX(1)
            INDEX(1) = IH
         ENDIF
 3001 CONTINUE

      RETURN
      END
