      SUBROUTINE PEAKFRI (F, IJS, IJL, IPEAKF, EPEAKF)

! ----------------------------------------------------------------------

!**** *PEAKFRI* - COMPUTATION OF INDEX OF THE PEAK FREQUENCY AND
!                 THE CORRESPONDING VALUE OF THE 1-D SPECTRAM
!                 AT EACH GRID POINT.

!     JEAN BIDLOT    ECMWF           MARCH 2000 

!*    PURPOSE.
!     --------

!       COMPUTE PEAK FREQUENCY INDEX AT EACH GRID POINT.
!       IF NOT PEAK FOUND THAN EPEAKF WILL BE SET TO 0. AND
!       IPEAKF TO NFRE

!**   INTERFACE.
!     ----------

!       *CALL* *PEAKFRI (F, IJS, IJL, IPEAKF, EPEAKF)*
!          *F*     - SPECTRUM.
!          *IJS*   - INDEX OF FIRST GRIDPOINT
!          *IJL*   - INDEX OF LAST GRIDPOINT
!          *IPEAKF - INDEX OF PEAK FREQUENCY 
!          *EPEAKF* - 1-D SPECTRUM AT PEAK FREQUENCY

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWFRED  , ONLY : DELTH    ,FR
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IPEAKF(IJS:IJL) 
      INTEGER :: IJ,K,M,IJS,IJL
      REAL, DIMENSION(IJS:IJL) :: TEMP, EPEAKF
      REAL :: F(NINF-1:NSUP,NANG,NFRE)

!*    1. INITIALIZE ARRAYS
!        -----------------

      DO IJ = IJS,IJL
        EPEAKF(IJ) = 0.
        IPEAKF(IJ) = NFRE 
      ENDDO

!*    2. LOOP OVER FREQUENCIES
!        ---------------------

      DO M = 1,NFRE

!*    2.1 COMPUTE 1-D SPECTRUM (WITHOUT DELTA THETA)
!        -------------------------------------------

        DO IJ = IJS,IJL
          TEMP(IJ) = 0.
        ENDDO
        DO K = 1,NANG
          DO IJ = IJS,IJL
            TEMP(IJ) = TEMP(IJ) + F(IJ,K,M)
          ENDDO
        ENDDO

!*    2.2 DEFINE PEAK FREQUENCY
!         ---------------------

        DO IJ = IJS,IJL
          IF (EPEAKF(IJ).LT.TEMP(IJ)) THEN
            EPEAKF(IJ) = TEMP(IJ)
            IPEAKF(IJ) = M 
          ENDIF
        ENDDO

      ENDDO

!     RESCALE EPEAKF TO 1-D SPECTRUM OF WAVE ENERGY

      DO IJ = IJS,IJL
        EPEAKF(IJ) = EPEAKF(IJ)*DELTH 
      ENDDO

      RETURN
      END SUBROUTINE PEAKFRI
