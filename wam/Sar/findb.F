
C=======================================================================
      SUBROUTINE FINDB (NDIM, NBOUN, BLATB, BLNGB,
     .                  amonop, amoeap, amosop, amowep,iper,
     .                  pdello, xdella, klonrgg, kgy,
     .                  nblo,niblo,igl,ijs,ijl,ixlg,kxlt,
     .                  IGARB, IJARB)
C
C**** *FINDB* - FIND BLOCK AND GRID POINT NUMBER.
C
C     R. PORTZ     MPI         15/01/1991
C
C     PURPOSE.
C     -------
C
C       FIND BLOCK AND GRID POINT NUMBER FOR A GIVEN ARRAY
C       OF LONGITUDES AND LATITUDES.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *FINDB (NDIM, NBOUN, BLATB, BLNGB, IGARB, IJARB)*
C          *NDIM*    - DIMENSION OF ARRAYS.
C          *NBOUN*   - NUMBER OF POINTS IN ARRAYS.
C          *BLATB*   - INPUT LATITUDES.
C          *BLNGB*   - INPUT LONGITUDES.
C          *IGARB*   - OUTPUT BLOCK NUMBERS.
C          *IJARB*   - OUTPUT SEA POINT NUMBERS.
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
C     MODIFICATIONS.
C     --------------
C
C     B. HANSEN  *ECMWF*     FEB 93'
C       COMMON BLOCKS DELETED .
C       NEEDED ARRAYS IMPORTED BY DUMMY ARGUMENTS.
C
C ----------------------------------------------------------------------
C
      DIMENSION
     .     IJS(NBLO), 
     .     IJL(NBLO),
     .     IXLG(NIBLO,NBLO), 
     .     KXLT(NIBLO,NBLO) 
C
      REAL PDELLO(KGY)
      INTEGER KLONRGG(KGY)
C
C ----------------------------------------------------------------------
C
      DIMENSION BLATB(NDIM), BLNGB(NDIM), IGARB(NDIM), IJARB(NDIM)
C
C ----------------------------------------------------------------------
C
C*    1. LOOP OVER INPUT LATITUDES, LONGITUDES.
C        --------------------------------------
C
 1000 CONTINUE
      DO 1001 IO = 1,NBOUN
C
C*    1.1 COMPUTE GRID MATRIX INDICES.
C         ----------------------------
C
        IOLT = NINT((BLATB(IO)-AMOSOP)/XDELLA+1.0)
        ALONG = MOD(BLNGB(IO)-AMOWEP+720.,360.)
        IF (IOLT.GE. 1 .AND. IOLT.LE. KGY) THEN
          IOLG = NINT(ALONG/PDELLO(IOLT)+1.0)
          IF (IOLG.EQ.(KLONRGG(IOLT)+1) .AND. IPER.EQ.1) IOLG = 1
Cbh
Cbh       IF ((AMOWEP+(IOLG-1)*KLONRGG(IOLT).EQ.AMOEAP+XDELLO)
Cbh    1     .AND. IPER.EQ.1) IOLG = 1
Cbh
C
C*    1.2 SEARCH BLOCK NUMBER AND SEA POINT NUMBER.
C         -----------------------------------------
C
          DO 1201 IG=1,IGL
          DO 1201 IJ = IJS(IG),IJL(IG)
            IF (IXLG(IJ,IG).EQ.IOLG .AND. KXLT(IJ,IG).EQ.IOLT) THEN
              IGARB(IO) = IG
              IJARB(IO) = IJ
              GOTO 1001
            ENDIF
 1201     CONTINUE
        ENDIF
        IGARB(IO) = 0
        IJARB(IO) = 0
 1001 CONTINUE

      RETURN
      END
