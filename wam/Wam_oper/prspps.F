      SUBROUTINE PRSPPS (IUNIT, CDATE, XLONG, XLAT, NAME,
     1                   WAHTWS, WADIWS, WAMFWS, WAHTSW, WADISW, WAMFSW,
     2                   FR, THETA, SPEC, IDD, IDF, NANG, NFRE)

C -----------------------------------------------------------------
C
C**** *PRSPPS* -  PRINTS A SPECTRUM, AND WAVE HEIGHT, WAVE DIRECTION,
C****             WAVE MEAN FREQUENCY OF WIND WAVES AND SWELL.
C
C         H. GUNTHER       ECMWF     OCTOBER  1990
C
C*    PURPOSE.
C     --------
C
C        PRINT A WAVE MODEL SWELL SPECTRUM.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *PRSPPS (IUNIT, CDATE, XLONG, XLAT, NAME,
C                       WAHTWS, WADIWS, WAMFWS, WAHTSW, WADISW, WAMFSW,
C                       FR, THETA, SPEC, IDD, IDF, NANG, NFRE)*
C         *IUNIT*   INTEGER    OUTPUT UNIT.
C         *CDATE*   CHAR*12    DATE OF SPECTRUM (YYYYMMDDHHMM).
C         *XLONG*   REAL       LONGITUDE OF SPECTRUM (DEGREE).
C         *XLAT*    REAL       LATITUDE OF SPECTRUM (DEGREE).
C         *NAME*    CHARACTER  SIDE NAME.
C         *WAHTWS*  REAL       WINDSEA WAVE HEIGHT (METRES).
C         *WADIWS*  REAL       WINDSEA MEAN WAVE DIRECTION (DEGREE).
C         *WAMFWS*  REAL       WINDSEA MEAN WAVE FREQUENCY (HERTZ).
C         *WAHTSW*  REAL       SWELL WAVE HEIGHT (METRES).
C         *WADISW*  REAL       SWELL MEAN WAVE DIRECTION (DEGREE).
C         *WAMFSW*  REAL       SWELL MEAN WAVE FREQUENCY (HERTZ).
C         *FR*      REAL       FREQUENCY ARRAY IN HERTZ.
C         *THETA*   REAL       DIRECTION ARRAY IN RAD.
C         *SPEC*    REAL       SPECTRUM.
C         *IDD*     INTEGER    DIMENSION OF DIRECTION ARRAYS.
C         *IDF*     INTEGER    DIMENSION OF FREQUENCY ARRAYS.
C         *NANG*    INTEGER    NUMBER OF DIRECTIONS.
C         *NFRE*    INTEGER    NUMBER OF FREQUENCIES.
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
C ------------------------------------------------------------------
C
#include "parcons.h"
C
      CHARACTER NAME*20, CDATE*12
      DIMENSION FR(IDF), THETA(IDD), SPEC(IDD,IDF)
C
      PARAMETER (IPDIR = 12)
      DIMENSION ANG(NANG), ODSPEC(NFRE)
      CHARACTER FORM1*50, FORM2*50, FORM3*50
C
C       *IPDIR*   INTEGER   NUMBER OF DIRECTIONS PRINTED PER LINE.
C       *ANG*     REAL      DIRECTIONS IN DEGREE.
C       *ODSPEC*  REAL      1-D SPECTRUM.  (M*M/HERTZ)
C       *FORM *   CHAR      VARIABLE FORMATS.
C
C ------------------------------------------------------------------
C
C*    1. INITIALISE DIRECTIONS.
C        ----------------------
C
      DELTH = ZPI/REAL(NANG)
      DO 1000 K = 1,NANG
         ANG(K) = THETA(K)*DEG
 1000 CONTINUE
C
C*    2. COMPUTE 1-D SPECTRUM.
C        ---------------------
C
      DO 2000 M = 1,NFRE
         ODSPEC(M) = 0.0
         DO 2100 K = 1,NANG
             ODSPEC(M) = ODSPEC(M) + SPEC(K,M)
 2100    CONTINUE
         ODSPEC(M) = ODSPEC(M)*DELTH
 2000 CONTINUE
C
C*    3. PRINT SPECTRUM.
C        ---------------
C
      WRITE(IUNIT,4100) CDATE, NAME, XLONG, XLAT
 4100 FORMAT('1',2X,'SEA AND SWELL INFORMATION. DATE: ',A12,
     1       ' SIDE: ',A20,3X,'LONG.: ',F7.2,' LAT.: ',F6.2)
      WRITE(IUNIT,4200) WAHTWS, WADIWS, WAMFWS, WAHTSW, WADISW, WAMFSW
 4200 FORMAT(1X,'SEA: HS = ',F6.2,'  THETA = ',F6.0,'  FMEAN = ',F6.4,
     1     '   SWELL: HS = ',F6.2,'  THETA = ',F6.0,'  FMEAN = ',F6.4/)
      IPE = 0
      IF (NANG.GT.IPDIR) THEN
         FORM1 = '(1X,''DIR (DEG)'',T11,12F8.1,''  DIR (DEG)'')'
         WRITE (FORM1(21:22),'(I2)') IPDIR
         FORM2 = '(1X,F7.4,2X,12F8.3,F9.4)'
         WRITE (FORM2(13:14),'(I2)') IPDIR
         FORM3 = '(1X,''FREQ (HZ)'',T11, 96X,''  FREQ (HZ)'')'
         WRITE (FORM3(21:23),'(I3)') IPDIR*8
         DO 4000 IP=1,NANG-IPDIR,IPDIR
            IPE = IP+IPDIR-1
            WRITE(IUNIT,FORM1) (ANG(K),K=IP,IPE)
            WRITE(IUNIT,FORM3)
            WRITE(IUNIT,FORM2)
     1         (FR(M),(SPEC(K,M),K=IP,IPE),FR(M),M=1,NFRE)
            WRITE(IUNIT,FORM1) (ANG(K),K=IP,IPE)
            WRITE(IUNIT,'(1X)')
 4000    CONTINUE
      ENDIF
      IP = IPE+1
      IPE = NANG
      LEN = IPE-IP+1
      FORM1 = '(1X,''DIR (DEG)'',T11,12F8.1,''  DIR (DEG)'')'
      WRITE (FORM1(21:22),'(I2)') LEN
      FORM2 = '(1X,F7.4,2X,12F8.3,F9.4,F10.3)'
      WRITE (FORM2(13:14),'(I2)') LEN
      FORM3 = '(1X,''FREQ (HZ)'',T11, 96X,''  FREQ (HZ) 1-DSPEC'')'
      WRITE (FORM3(21:23),'(I3)') LEN*8
      WRITE (IUNIT,FORM1) (ANG(K),K=IP,IPE)
      WRITE (IUNIT,FORM3)
      WRITE (IUNIT,FORM2)
     1   (FR(M),(SPEC(K,M),K=IP,IPE),FR(M),ODSPEC(M),M=1,NFRE)
      WRITE (IUNIT,FORM1) (ANG(K),K=IP,IPE)

      RETURN
      END
