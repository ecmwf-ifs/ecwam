      SUBROUTINE outpp (idate, iuout, id1, id2, ngx, ngy, titl, const,
     1                  array, amowep, amosop, amoeap, amonop)

C ---------------------------------------------------------------
C
C**** *OUTPP* - FORMATED OUTPUT OF AN ARRAY.
C
C     H. GUNTHER       ECMWF    NOVEMBER 1989
C
C*    PURPOSE.
C     --------
C
C       FORMATED OUTPUT OF AN ARRAY.
C
C**   INTERFACE.
C     ---------
C
C       *CALL* *OUTPP (IDATE, IUOUT, ID1, ID2, NGX, NGY, TITL, CONST,
C                      ARRAY, AMOWEP, AMOSOP, AMOEAP, AMONOP)*
C          *IDATE*    INTEGER   DATE (YYMMDDHHMM).
C          *IUOUT*    INTEGER   OUTPUT UNIT.
C          *ID1*      INTEGER   FIRST DIMENSION.
C          *ID2*      INTEGER   SECOND DIMENSION.
C          *NGX*      INTEGER   FIRST DIMENSION  USED.
C          *NGY*      INTEGER   SECOND DIMENSION USED.
C          *TITL*     CHARACTER HEADER TO BE PRINTED.
C          *CONST*    REAL      SCALING FACTOR.
C          *ARRAY*    REAL      ARRAY TO BE PRINTED.
C          *AMOWEP*   REAL      MOST WESTERN LONGITUDE (DEGREE).
C          *AMOSOP*   REAL      MOST SOUTHERN LATITUDE (DEGREE).
C          *AMOEAP*   REAL      MOST EASTERN LONGITUDE (DEGREE).
C          *AMONOP*   REAL      MOST NORTHERN LATITUDE (DEGREE).
C
C     METHOD.
C     -------
C
C       A TWO DIMENSIONAL ARRAY IS PRINTED WITH A MAXIMUM OF
C       30 COLUMNS PER PAGE.
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

      CHARACTER TITL*100
      DIMENSION IARRAY(ID1,ID2),ARRAY(ID1,ID2),ILON(NGX),YLAT(NGY)
C
C ----------------------------------------------------------------------
C
C
C*    1. INITIALIZATION
C     -----------------
C
      DLAMA = (AMOEAP-AMOWEP)/FLOAT(NGX-1)
      DO 110 I=1,NGX
         ILON(I)=NINT(AMOWEP + (I-1)*DLAMA)
 110  CONTINUE

      DPHIA = (AMOSOP-AMONOP)/FLOAT(NGY-1)
      DO 120 J=1,NGY
         YLAT(J)=AMONOP + (J-1)*DPHIA
 120  CONTINUE

      NPTS=30
      NPAGE=(NGX+NPTS-1)/NPTS

      DO 100 J=1,NGY
      DO 100 I=1,NGX
         IARRAY(I,J) = NINT(CONST*ARRAY(I,J))
 100  CONTINUE

      ISTART=-NPTS+1
      IEND=ISTART+NPTS-1

      DO 200 NP=1,NPAGE
         WRITE (IUOUT,300) IDATE,TITL,NP
         ISTART=ISTART+NPTS
         IEND= MIN(IEND+NPTS,NGX)
         WRITE (IUOUT,302) (I,I=ISTART,IEND)
         WRITE (IUOUT,303) (ILON(I),I=ISTART,IEND)
         WRITE (IUOUT,304)
         DO 220 J=1,NGY
            WRITE (IUOUT,305) J,YLAT(J),(IARRAY(I,J),I=ISTART,IEND)
 220     CONTINUE
 200  CONTINUE

 300  FORMAT(1H1,6X,I10.10,2X,A100,5X,'PAGE ',I2,/)
 302  FORMAT(7X,'I=',30I4)
 303  FORMAT(5X,'LON=',30I4)
 304  FORMAT('   J LAT',/)
 305  FORMAT(1X,I2,F5.1,1X,30I4)

      RETURN
      END
