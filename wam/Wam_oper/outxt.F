      SUBROUTINE OUTXT (KU06,KU1,KU2,KU3,KU4,KDEL,KASSIM,KOUT,
     &                       YA1,YA2,YA3,YA4)
! ----------------------------------------------------------------------

!**** *OUTXT* - WRITES USAGE OF OUTPUT UNITS.

!     B. HANSEN           ECMWF             APRIL 1993

!*    PURPOSE.
!     --------

!       REDUCTION OF CODE REPETITION IN THE CALLING SUBROUTINE *USERIN*.


!**   INTERFACE.
!     ----------
!          *KU06*      - PRINTER OUTPUT UNIT.
!          *KU1*       - UNIT FOR UNFORMATTED DATA.
!          *KU2*       - UNIT FOR UNFORMATTED DATA OF DATA ASSIMILATION.
!          *KU3*       - UNIT FOR GRIB FORMATTED DATA.
!          *KU4*       - UNIT FOR GRIB FORMATTED DATA OF THE DATA ASSIM.
!          *KDEL*      - TIMESTEP OF THE OUTPUT.
!          *KASSIM*    - NO ASSIMILATION IF KASSIM = 0, DATA ASSIM IF 1.
!          *KOUT*      - NUMBER OF OUTPUT TIMES.
!          *YA1*       - NAME OF THE VARIABLE OF KU1.
!          *YA2*       - NAME OF THE VARIABLE OF KU2.
!          *YA3*       - NAME OF THE VARIABLE OF KU3.
!          *YA4*       - NAME OF THE VARIABLE OF KU4.

!     EXTERNALS.
!     ----------
!       NONE.

!     METHOD.
!     -------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: KU06,KU1,KU2,KU3,KU4,KDEL,KASSIM,KOUT
      CHARACTER(LEN=4), INTENT(IN) :: YA1,YA2,YA3,YA4

      INTEGER :: J
      INTEGER :: KU(4)

      CHARACTER(LEN=3) :: YOU(4)
      CHARACTER(LEN=10) :: YOKDEL
      CHARACTER(LEN=26) :: YUL

!*    1. NOT MUCH TO SAY.
!        ----------------

      KU(1)=KU1
      KU(2)=KU2
      KU(3)=KU3
      KU(4)=KU4
      DO J=1,4
        IF(KU(J).GT.111) THEN
          WRITE(YOU(J),'(A3)')KU(J)
        ELSE
          WRITE(YOU(J),'(I3)')KU(J)
        ENDIF
      ENDDO
      IF ( KOUT .EQ. 0 ) THEN
        WRITE(YOKDEL,'(I10)')KDEL
        YUL=' EVERY '//YOKDEL//' SECONDS '
      ELSE
        YUL=' '
      ENDIF
      IF ( KASSIM .EQ.0 ) THEN
        WRITE(KU06,*)'   TO PRINTER AND/OR UNIT '//YA1//' = '//YOU(1)
     &             //' AND/OR'
        WRITE(KU06,*)'   GRIBBED DATA TO UNIT   '//YA2//' = '//YOU(2)
     &             //YUL
      ELSE
        WRITE(KU06,*)'   TO PRINTER AND/OR UNIT '//YA1//' = '//YOU(1)
     &             //' ('//YA3//' = '//YOU(3)//') AND/OR'
        WRITE(KU06,*)'   GRIBBED DATA TO UNIT   '//YA2//' = '//YOU(2)
     &             //' ('//YA4//' = '//YOU(4)//') '//YUL
      ENDIF
      RETURN
      END SUBROUTINE OUTXT
