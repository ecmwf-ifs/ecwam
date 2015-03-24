      SUBROUTINE PRSPP (IUNIT, CDATE, XLONG, XLAT, NAME,
     &                  WAVEHT, WAVEDIR, WAVEMFR, WINDSPD, WINDDIR,
     &                  WIND10, FR, THETA, SPEC, IDD, IDF, NANG, NFRE)

! ----------------------------------------------------------------------

!***  *PRSPP* -  PRINTS A SPECTRUM, THE WAVE HEIGHT , WAVE DIRECTION,
!*               USTAR WIND DIRECTION, AND U10.

!         M. DE LAS HERAS  KNMI/PCM  FEBRUARY  1990

!*    PURPOSE.
!     --------

!        PRINT A WAVE MODEL SPECTRUM.

!**   INTERFACE.
!     ----------

!       *CALL* *PRSPP (IUNIT, CDATE, XLONG, XLAT, NAME,
!                      WAVEHT, WAVEDIR, WAVEMFR, WINDSPD, WINDDIR,
!                      WIND10, FR, THETA, SPEC, IDD, IDF, NANG, NFRE)*
!         *IUNIT*   INTEGER    OUTPUT UNIT.
!         *CDATE*   CHAR*14    DATE OF SPECTRUM (YYYYMMDDHHMMSS).
!         *XLONG*   REAL       LONGITUDE OF SPECTRUM (DEGREE).
!         *XLAT*    REAL       LATITUDE OF SPECTRUM (DEGREE).
!         *NAME*    CHARACTER  SIDE NAME.
!         *WAVEHT*  REAL       SIGNIFICANT WAVE HEIGHT (METRES).
!         *WAVEDIR* REAL       MEAN WAVE DIRECTION (DEGREE).
!         *WAVEMFR* REAL       MEAN WAVE FREQUENCY (HERTZ).
!         *WINDSPD* REAL       FRICTION VELOCITY (METRES/SECOND).
!         *WINDDIR* REAL       WIND DIRECTION (DEGREES).
!         *WIND10*  REAL       WIND SPEED U10 (METRES/SECOND).
!         *FR*      REAL       FREQUENCY ARRAY IN HERTZ.
!         *THETA*   REAL       DIRECTION ARRAY IN RAD.
!         *SPEC*    REAL       SPECTRUM.
!         *IDD*     INTEGER    DIMENSION OF DIRECTION ARRAYS.
!         *IDF*     INTEGER    DIMENSION OF FREQUENCY ARRAYS.
!         *NANG*    INTEGER    NUMBER OF DIRECTIONS.
!         *NFRE*    INTEGER    NUMBER OF FREQUENCIES.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ------------------------------------------------------------------


      USE YOWPCONS , ONLY : ZPI      ,DEG


      CHARACTER NAME*20, CDATE*14
      DIMENSION FR(IDF), THETA(IDD), SPEC(IDD,IDF)

      PARAMETER (IPDIR = 12)
      DIMENSION ANG(NANG), ODSPEC(NFRE)
      CHARACTER FORM1*50, FORM2*50, FORM3*50

!       *IPDIR*   INTEGER   NUMBER OF DIRECTIONS PRINTED PER LINE.
!       *ANG*     REAL      DIRECTIONS IN DEGREE.
!       *ODSPEC*  REAL      1-D SPECTRUM.  (M*M/HERTZ)
!       *FORM *   CHAR      VARIABLE FORMATS.

! ------------------------------------------------------------------

!*    1. INITIALISE DIRECTIONS.
!        ----------------------

      DELTH = ZPI/REAL(NANG)
      DO K = 1,NANG
        ANG(K) = THETA(K)*DEG
      ENDDO

!*    2. COMPUTE 1-D SPECTRUM.
!        ---------------------

      DO M = 1,NFRE
        ODSPEC(M) = 0.0
        DO K = 1,NANG
          ODSPEC(M) = ODSPEC(M) + SPEC(K,M)
        ENDDO
        ODSPEC(M) = ODSPEC(M)*DELTH
      ENDDO

!*    3. PRINT SPECTRUM.
!        ---------------

      WRITE(IUNIT,4100) CDATE, NAME, XLONG, XLAT
 4100 FORMAT('1',2X,'WAVE AND WIND INFORMATION. DATE: ',A14,
     &       ' SIDE: ',A20,3X,'LONG.: ',F9.2,' LAT.: ',F9.2)
      WRITE(IUNIT,4200) WAVEHT, WAVEDIR, WAVEMFR, WINDSPD, WIND10,
     &                  WINDDIR
 4200 FORMAT(3X,'HS = ',F9.2,'  THETA = ',F6.0,' FMEAN = ',F11.4,
     &          ' USTAR = ',F9.2,' U10 = ',F9.2,'  PHI = ',F6.0/)
      IPE = 0
      IF (NANG.GT.IPDIR) THEN
        FORM1 = '(1X,''DIR (DEG)'',T11,12F8.1,''  DIR (DEG)'')'
        WRITE (FORM1(21:22),'(I2)') IPDIR
        FORM2 = '(1X,F7.4,2X,12F8.3,F9.4)'
        WRITE (FORM2(13:14),'(I2)') IPDIR
        FORM3 = '(1X,''FREQ (HZ)'',T11, 96X,''  FREQ (HZ)'')'
        WRITE (FORM3(21:23),'(I3)') IPDIR*8
        DO IP=1,NANG-IPDIR,IPDIR
          IPE = IP+IPDIR-1
          WRITE(IUNIT,FORM1) (ANG(K),K=IP,IPE)
          WRITE(IUNIT,FORM3)
          WRITE(IUNIT,FORM2)
     &     (FR(M),(SPEC(K,M),K=IP,IPE),FR(M),M=1,NFRE)
          WRITE(IUNIT,FORM1) (ANG(K),K=IP,IPE)
          WRITE(IUNIT,'(1X)')
        ENDDO
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
     &   (FR(M),(SPEC(K,M),K=IP,IPE),FR(M),ODSPEC(M),M=1,NFRE)
      WRITE (IUNIT,FORM1) (ANG(K),K=IP,IPE)

      RETURN
      END SUBROUTINE PRSPP


