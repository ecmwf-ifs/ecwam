      SUBROUTINE OUTXT (KU06,KU1,KU2,KU3,KU4,KDEL,KASSIM,KOUT,
     1                       YA1,YA2,YA3,YA4)
C ----------------------------------------------------------------------
C
C**** *OUTXT* - WRITES USAGE OF OUTPUT UNITS.
C
C     B. HASNEN           ECMWF             APRIL 1993
C
C*    PURPOSE.
C     --------
C
C       REDUCTION OF CODE REPITION IN THE CALLING SUBROUTINE *USERIN*.
C
C
C**   INTERFACE.
C     ----------
C          *KU06*      - PRINTER OUTPUT UNIT.
C          *KU1*       - UNIT FOR UNFORMATTED DATA.
C          *KU2*       - UNIT FOR UNFORMATTED DATA OF THE DATA ASSIMILATION.
C          *KU3*       - UNIT FOR GRIB FORMATTED DATA.
C          *KU4*       - UNIT FOR GRIB FORMATTED DATA OF THE DATA ASSIM.
C          *KDEL*      - TIMESTEP OF THE OUTPUT.
C          *KASSIM*    - NO ASSIMILATION IF KASSIM = 0, DATA ASSIM IF 1.
C          *KOUT*      - NUMBER OF OUTPUT TIMES.
C          *YA1*       - NAME OF THE VARIABLE OF KU1.
C          *YA2*       - NAME OF THE VARIABLE OF KU2.
C          *YA3*       - NAME OF THE VARIABLE OF KU3.
C          *YA4*       - NAME OF THE VARIABLE OF KU4.
C
C     EXTERNALS.
C     ----------
C       NONE.
C
C     METHOD.
C     -------
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
      CHARACTER* 3 YOU(4)
      CHARACTER* 4 YA1,YA2,YA3,YA4
      CHARACTER*10 YOKDEL
      CHARACTER*26 YUL
      INTEGER KU(4)
C
C*    1. NOT MUCH TO SAY.
C        ----------------
C
      KU(1)=KU1
      KU(2)=KU2
      KU(3)=KU3
      KU(4)=KU4
      DO 10 J=1,4
        IF(KU(J).GT.111) THEN
          WRITE(YOU(J),'(A3)')KU(J)
        ELSE
          WRITE(YOU(J),'(I3)')KU(J)
        ENDIF
   10 CONTINUE
      IF ( KOUT .EQ. 0 ) THEN
        WRITE(YOKDEL,'(I10)')KDEL
        YUL=' EVERY '//YOKDEL//' SECONDS '
      ELSE
        YUL=' '
      ENDIF
      IF ( KASSIM .EQ.0 ) THEN
        WRITE(KU06,*)'   TO PRINTER AND/OR UNIT '//YA1//' = '//YOU(1)
     1             //' AND/OR'
        WRITE(KU06,*)'   GRIBBED DATA TO UNIT   '//YA2//' = '//YOU(2)
     1             //YUL
      ELSE
        WRITE(KU06,*)'   TO PRINTER AND/OR UNIT '//YA1//' = '//YOU(1)
     1             //' ('//YA3//' = '//YOU(3)//') AND/OR'
        WRITE(KU06,*)'   GRIBBED DATA TO UNIT   '//YA2//' = '//YOU(2)
     1             //' ('//YA4//' = '//YOU(4)//') '//YUL
      ENDIF
      RETURN
      END
