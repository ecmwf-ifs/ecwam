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

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: NFRE, NANG, ML, KL
      REAL(KIND=JWRB), INTENT(IN) :: RTHET
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(OUT) :: FL3


      INTEGER(KIND=JWIM) :: K, M
      INTEGER(KIND=JWIM) :: INC, KC, KC1
      REAL(KIND=JWRB) :: ZPI, FTH, ADIF, BDIF

! --------------------------------------------------------------------

      ZPI=8.0_JWRB*ATAN(1._JWRB)

      FTH = MOD(RTHET+ZPI,ZPI)
      FTH = FTH * REAL(KL) / ZPI
      INC = INT(FTH)
      ADIF = FTH - INC
      BDIF = 1.0_JWRB - ADIF

      DO K=1,KL
        KC  = K  - INC
        IF (KC .LT. 1) KC  = KC + KL
        KC1 = KC - 1
        IF (KC1.LT. 1) KC1 = KC1+ KL

        DO M=1,ML
          FL3(K,M) = BDIF * FL1(KC,M) + ADIF * FL1(KC1,M)
        ENDDO
      ENDDO

      END SUBROUTINE ROTSPEC
