      SUBROUTINE ROTSPEC (NFRE, NANG, ML, KL, FL1, FL3, RTHET)

!----------------------------------------------------------------------

!**** *ROTSPEC* - ROUTINE TO ROTATE THE SPECTRUM.

!     EVA BAUER      MPI  HAMBURG    MAY 1990.
!     H. GUNTHER          GKSS/ECMWF JAN. 1991   MODIFIED FOR CYCLE_4

!*    PURPOSE.
!     -------

!       TO ROTATE THE SPECTRUM.

!**   INTERFACE.
!     ----------

!       *CALL* *ROTSPEC (NFRE, NANG, ML, KL, FL1, FL3, RTHET)*
!          *NFRE*  - FREQUENCY DIMENSION OF SPECTRA.
!          *NANG*  - DIRECTION DIMENSION OF SPECTRA.
!          *ML*    - NUMBER OF FREQUENCIES.
!          *KL*    - NUMBER OF DIRECTIONS.
!          *FL1*   - SPECTRUM TO BE ROTATED.
!          *FL3*   - ROTATED SPECTRUM.
!          *RTHET* - TURNING ANGLE IN RADIANS, CLOCKWISE.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ---------

!       NONE.

!     REFERENCES.
!     -----------

!       NONE.

! -------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NFRE, NANG, ML, KL
      REAL, INTENT(IN) :: RTHET
      REAL, DIMENSION(NANG,NFRE), INTENT(IN) :: FL1
      REAL, DIMENSION(NANG,NFRE), INTENT(OUT) :: FL3


      INTEGER :: K, M
      INTEGER :: INC, KC, KC1
      REAL :: ZPI, FTH, ADIF, BDIF

! --------------------------------------------------------------------

      ZPI=8.*ATAN(1.)

      FTH = MOD(RTHET+ZPI,ZPI)
      FTH = FTH * REAL(KL) / ZPI
      INC = INT(FTH)
      ADIF = FTH - INC
      BDIF = 1. - ADIF

      DO K=1,KL
        KC  = K  - INC
        IF (KC .LT. 1) KC  = KC + KL
        KC1 = KC - 1
        IF (KC1.LT. 1) KC1 = KC1+ KL

        DO M=1,ML
          FL3(K,M) = BDIF * FL1(KC,M) + ADIF * FL1(KC1,M)
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE ROTSPEC
