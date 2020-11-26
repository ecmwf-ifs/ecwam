      SUBROUTINE FNDPRT (IJS, IJL, NPMAX,                               &
     &                   NPEAK, MIJ, NTHP, NFRP,                        &
     &                   FLLOW, LLCOSDIFF, FLNOISE,                     &
     &                   FL3, SWM,                                      &
     &                   ENE, DIR, PER)

! ----------------------------------------------------------------------

!**** *FNDPRT* -  CALCULATE PARTITION MASKS

!     D.PETTENUZZO     MAY 2011
!     MODIFIED MAY 2012 TO SUIT ECMWF CODES


!*    PURPOSE.
!     --------

!     FIND ALL THE POINTS IN THE FREQUENCY DIRECTION DOMAIN
!     WHICH BELONG TO A GIVEN PARTITION IDENTIFIED BY ITS
!     PEAK

!**   INTERFACE.
!     ----------

!       *CALL* *FNDPRT (IJS, IJL, NPMAX,
!     &                 NPEAK, MIJ, NTHP, NFRP,
!     &                 FLLOW, LLCOSDIFF, FLNOISE,
!     &                 FL3, SWM,
!     &                 ENE, DIR, PER)

!          *IJS*    - INDEX OF FIRST GRIDPOINT
!          *IJL*    - INDEX OF LAST GRIDPOINT
!          *NPMAX*  - MAXIMUM NUMBER OF PARTITIONS
!          *NPEAK*  - NUMBER OF PEAKS
!          *MIJ*    - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!          *NTHP*   - DIRECTION INDEX FOR ALL PEAKS
!          *NFRP*   - FREQUENCY INDEX FOR ALL PEAKS
!          *FLLOW*  - MINIMUM ENERGY LEVEL
!          *LLCOSDIFF* LOGICAL (COSDIFF.LT.-0.4)
!          *FLNOISE* - SPECTRAL NOISE LEVEL
!          *FL3*    - BLOCK OF SPECTRA
!          *SWM*    - ORIGINAL SWELL MASK
!                     IT MIGHT NEED TO BE ADJUSTED
!                     THIS IS POSSIBLE BECAUSE WE SUPPLY SWELL SPECTRA
!                     IN WHICH THE WIND SEA HAS BEEN REMOVED. SOME OF
!                     THE ENERGY LEAKING OUT FROM THE WINDSEA SPECTRUM
!          *EME*    - MEAN WAVE ENERGY FOR EACH WAVE SYSTEM 
!          *DIR*    - MEAN WAVE DIRECTION FOR EACH WAVE SYSTEM
!          *PER*    - MEAN PERIOD BASED ON -1 MOMENT FOR EACH WAVE SYSTEM.


!     METHOD.
!     -------

!       STEEPEST ASCENT

!     EXTERNALS.
!     ----------
!       *PARMEAN* - COMPUTE MEAN PARAMETERS.

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "parmean.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, NPMAX
      INTEGER(KIND=JWIM), INTENT(IN), DIMENSION(IJS:IJL) :: MIJ
      INTEGER(KIND=JWIM), INTENT(INOUT), DIMENSION(IJS:IJL) :: NPEAK
      INTEGER(KIND=JWIM), INTENT(IN), DIMENSION(IJS:IJL,NPMAX) :: NTHP, NFRP

      REAL(KIND=JWRB), INTENT(IN), DIMENSION(IJS:IJL) :: FLNOISE
      REAL(KIND=JWRB), INTENT(IN), DIMENSION(IJS:IJL,NANG,NFRE) :: FLLOW, FL3
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: SWM
      REAL(KIND=JWRB), INTENT(OUT), DIMENSION(IJS:IJL,0:NPMAX) :: DIR, PER, ENE

      LOGICAL, INTENT(IN), DIMENSION(IJS:IJL,NANG) :: LLCOSDIFF

      INTEGER(KIND=JWIM) :: ITHC, IFRC
      INTEGER(KIND=JWIM) :: IJ, M, K, IP, NITT
      INTEGER(KIND=JWIM) :: NANGH, KK, KKMIN, KKMAX
      INTEGER(KIND=JWIM) :: IFRL, ITHL, ITHR
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL) :: MMIN, MMAX
      INTEGER(KIND=JWIM), DIMENSION(1-NANG:2*NANG) :: KLOC

      REAL(KIND=JWRB) :: HALF_SECTOR
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE) :: W2
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE) :: W1
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE,NPMAX,IJS:IJL) :: SPEC

      LOGICAL :: LLCHANGE, LLADD
      LOGICAL :: LLADDPART
      LOGICAL, DIMENSION(IJS:IJL,NANG,NFRE) :: LLW3

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('FNDPRT',0,ZHOOK_HANDLE)

!     ASSUME THAT ONE HAS ONLY SEARCH A SECTOR IN DIRECTION (not the full 360 degrees)
!     AROUND THE PEAK
!!!   THIS IS A STRONG ASSUMPTION ON THE DIRECTIONALITY OF SWELL !!!

      HALF_SECTOR=75.0_JWRB
      NANGH=NINT((HALF_SECTOR/360.0_JWRB)*NANG)+1

!     POINTS BELOW THE MINIMUM LEVEL BELONG TO THE WINDSEA PART
!     AND CAN BE EXCLUDED FROM THE PARTITIONS OF THE SWELL SPECTRUM
      DO IJ=IJS,IJL
        MMIN(IJ)=NFRE
        MMAX(IJ)=0
      ENDDO
      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            IF(FL3(IJ,K,M).LE.FLLOW(IJ,K,M)) THEN
              W1(IJ,K,M) = 1._JWRB
              LLW3(IJ,K,M) = .FALSE.
            ELSE
              W1(IJ,K,M) = 0._JWRB
              LLW3(IJ,K,M) = .TRUE.
              MMIN(IJ)=MIN(M,MMIN(IJ))
              MMAX(IJ)=M
            ENDIF
          ENDDO
        ENDDO
      ENDDO


!     MAIN LOOP ON GRID POINTS
      DO IJ=IJS,IJL
        DO IP=1,NPEAK(IJ)
          ITHC=NTHP(IJ,IP)
          IFRC=NFRP(IJ,IP)

          KKMIN=ITHC-NANGH
          KKMAX=ITHC+NANGH
          DO KK=KKMIN,KKMAX
            KLOC(KK)=1+MOD(NANG+KK-1,NANG)
          ENDDO

!*        1.  SET UP THE W2 MAP
!             -----------------

          DO M=1,NFRE
            DO K=1,NANG
              W2(K,M) = 0.0_JWRB
            ENDDO
          ENDDO

          DO K=ITHC-1, ITHC+1
            ITHL = 1 + MOD(NANG+K-1,NANG)
            DO M=IFRC-1, IFRC+1
              IFRL = MAX(1,MIN(NFRE,M))
              IF (W1(IJ,ITHL,IFRL) .LE. 0.5_JWRB) W2(ITHL,IFRL)=0.5_JWRB
            ENDDO
          ENDDO

          IF (W1(IJ,ITHC,IFRC) .LT. 0.25_JWRB) W2(ITHC,IFRC)=1.0_JWRB

!         FIND IF MORE HIGH FREQUENCY BINS HAVE BECOME EXCLUDED
          OUT0: DO M=MMAX(IJ),MMIN(IJ),-1
            DO KK=KKMIN,KKMAX
              K=KLOC(KK)
              IF (W1(IJ,K,M) .LT. 1.0_JWRB) THEN
                MMAX(IJ)=M
                EXIT OUT0
              ENDIF
            ENDDO
          ENDDO OUT0

!*        2.  ITTERATE SEARCH
!             ---------------

          NITT=0

!     2.a Branch point

  200 CONTINUE
          NITT=NITT+1
          LLCHANGE=.FALSE.

!         2.b Determine central points

          DO M=MMIN(IJ),MIN(MIJ(IJ),MMAX(IJ))
!           by definition bins beyond M=MIJ are never extremas
!           and bins above MMAX are excluded.
            DO KK=KKMIN,KKMAX
              K = KLOC(KK)

              IF (LLW3(IJ,K,M)) THEN
                IF (W2(K,M).EQ.0.5_JWRB .AND.                           &
     &              W1(IJ,K,M).LT.0.5_JWRB    ) THEN
                  LLADD = .TRUE.
                  OUT1: DO ITHR=K-1, K+1
                    ITHL = 1 + MOD(NANG+ITHR-1,NANG)
                    DO IFRL=MAX(1,M-1), MIN(NFRE,M+1)
                      IF (W2 (ITHL,IFRL).EQ.0.0_JWRB .AND.             &
     &                    FL3(IJ,ITHL,IFRL).GT.FL3(IJ,K,M)) THEN
                        LLADD = .FALSE.
                        EXIT OUT1
                      ENDIF
                    ENDDO
                  ENDDO OUT1
                  IF (LLADD) THEN
                    W2(K,M) = 1.0_JWRB
                    LLCHANGE = .TRUE. 
                  ENDIF
                ENDIF
              ENDIF

            ENDDO
          ENDDO

!         2.c Determine peripherical points

          DO M=MMIN(IJ),MMAX(IJ)
            DO KK=KKMIN,KKMAX
              K=KLOC(KK)

              IF (LLW3(IJ,K,M) .AND. W1(IJ,K,M).LT.1.0_JWRB) THEN
                IF (W2(K,M).EQ.0.0_JWRB) THEN
                  OUT2: DO ITHR=K-1, K+1
                    ITHL = 1 + MOD(NANG+ITHR-1,NANG)
                    DO IFRL=MAX(1,M-1), MIN(NFRE,M+1)
                      IF (W2(ITHL,IFRL).EQ.1.0_JWRB) THEN
                        W2(K,M) = 0.5_JWRB
                        LLCHANGE = .TRUE.
                        EXIT OUT2
                      ENDIF
                    ENDDO
                  ENDDO OUT2
                ENDIF
              ENDIF

            ENDDO
          ENDDO

!*        2.d Branch back ?

          IF ( LLCHANGE .AND. NITT.LT.25 ) GOTO 200


!*        3   UPDATE THE OVERALL MAP

          DO M=1, NFRE
            DO K=1, NANG
              W1(IJ,K,M) = W1(IJ,K,M) + W2(K,M)
            ENDDO
          ENDDO

          DO M=1,NFRE
            DO K=1,NANG
              SPEC(K,M,IP,IJ) = FL3(IJ,K,M)*W2(K,M)
            ENDDO
          ENDDO

        ENDDO ! loop IP

!       CHECK THAT UNASSIGNED BINS ARE IN THE WIND SECTOR, OTHERWISE THEY SHOULD BE
!       ADDED AS A PARTITION IF ABOVE NOISE LEVEL
        LLADDPART = .FALSE.
        IF(NPEAK(IJ).LT.NPMAX) THEN
          DO K=1,NANG
            DO M=1,NFRE
              IF(LLCOSDIFF(IJ,K) .AND. W1(IJ,K,M).LE.0.0_JWRB .AND.     &
     &           FL3(IJ,K,M) .GT. FLNOISE(IJ) ) THEN
                W2(K,M) = 1.0_JWRB
                W1(IJ,K,M) = 1.0_JWRB
                LLADDPART = .TRUE.
              ELSE
                W2(K,M) = 0.0_JWRB
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        IF(LLADDPART) THEN
          NPEAK(IJ) = NPEAK(IJ)+1
          IP=NPEAK(IJ)
          DO M=1,NFRE
            DO K=1,NANG
              SPEC(K,M,IP,IJ) = FL3(IJ,K,M)*W2(K,M)
            ENDDO
          ENDDO
        ENDIF

      ENDDO ! loop IJ

      CALL PARMEAN(IJS, IJL, NPMAX, NPEAK,                              &
     &             SPEC,                                                &
     &             ENE, DIR, PER)


!     UPDATE SWELL MASK (IF IN THE WIND DIRECTION)
      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            IF(.NOT.LLCOSDIFF(IJ,K)) THEN
              SWM(IJ,K,M) = SWM(IJ,K,M)*MAX(W1(IJ,K,M),1.0_JWRB)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('FNDPRT',1,ZHOOK_HANDLE)

      END SUBROUTINE FNDPRT
