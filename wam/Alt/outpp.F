C ======================================================================
C
      SUBROUTINE OUTPP (cdate, kuout, kd1, kd2, kgx, kgy, cdtitl,
     .                  pconst, zfield, pmowep, pmosop, pmoeap, pmonop)
C
C ======================================================================
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
C       *CALL* *OUTPP (CDATE, kUOUT, kd1, kd2, kgx, kgy, cdtitl, pconst,
C                      zfield, pmowep, pmosop, pmoeap, pmonop)*
C          ALL VARIABLES ARE INPUT.
C          *CDATE*    CHARACTER DATE (YYYYMMDDHHMM).
C          *kUOUT*    INTEGER   OUTPUT UNIT.
C          *kd1*      INTEGER   FIRST DIMENSION.
C          *kd2*      INTEGER   SECOND DIMENSION.
C          *kgx*      INTEGER   FIRST DIMENSION  USED.
C          *kgy*      INTEGER   SECOND DIMENSION USED.
C          *cdtitl*   CHARACTER HEADER TO BE PRINTED.
C          *pconst*   REAL      SCALING FACTOR.
C          *zfield*   REAL      ARRAY TO BE PRINTED.
C          *pmowep*   REAL      MOST WESTERN LONGITUDE (DEGREE).
C          *pmosop*   REAL      MOST SOUTHERN LATITUDE (DEGREE).
C          *pmoeap*   REAL      MOST EASTERN LONGITUDE (DEGREE).
C          *pmonop*   REAL      MOST NORTHERN LATITUDE (DEGREE).
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
C
      CHARACTER cdtitl*100
      CHARACTER*12 CDATE
      DIMENSION iarray(kd1,kd2),zfield(kd1,kd2),ilon(kgx),zlat(kgy)
C
C ----------------------------------------------------------------------
C
C
C*    1. INITIALIZATION
C     -----------------
C
      zlonl=pmowep
      IF (pmowep.gt.pmoeap) zlonl=pmowep-360.
      zdelam = (pmoeap-zlonl)/FLOAT(kgx-1)
C
      DO 110 j110=1,kgx
         ilon(j110)=NINT(zlonl + (j110-1)*zdelam)
 110  CONTINUE

      zdphi = (pmosop-pmonop)/FLOAT(kgy-1)
      DO 120 j120 = 1,kgy
         zlat(j120)=pmonop + (j120-1)*zdphi
 120  CONTINUE

      ipts=30
      ipage=(kgx+ipts-1)/ipts

      DO 100 J100J=1,kgy
      DO 100 J100I=1,kgx
         iarray(J100I,J100J) = NINT(pconst*zfield(J100I,J100J))
 100  CONTINUE

      istart=-ipts+1
      iend=istart+ipts-1

      DO 200 j200=1,ipage
         WRITE (kuout,300) cdate,cdtitl,j200
         istart=istart+ipts
         iend= MIN(iend+ipts,kgx)
         WRITE (kuout,302) (I,I=istart,iend)
         WRITE (kuout,303) (ilon(ji),ji=istart,iend)
         WRITE (kuout,304)
         DO 220 j220=1,kgy
            WRITE (kuout,305)
     .      j220,zlat(j220),(iarray(ji,j220),ji=istart,iend)
 220     CONTINUE
 200  CONTINUE
 300  FORMAT(1H1,6X,A12,2X,A100,5X,'PAGE ',I2,/)
 302  FORMAT(7X,'I=',30I4)
 303  FORMAT(5X,'LON=',30I4)
 304  FORMAT('   J LAT',/)
 305  FORMAT(1X,I2,F5.1,1X,30I4)

      RETURN
      END
