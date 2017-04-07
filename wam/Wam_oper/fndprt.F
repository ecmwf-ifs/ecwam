      SUBROUTINE FNDPRT (IJS, IJL, NPMAX,
     &                   NPEAK, MIJ, NTHP, NFRP,
     &                   FLLOW, LLCOSDIFF,
     &                   FL3,
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
!     &                 FLLOW, LLCOSDIFF,
!     &                 FL3,
!     &                 ENE, DIR, PER)

!          *IJS*    - INDEX OF FIRST GRIDPOINT
!          *IJL*    - INDEX OF LAST GRIDPOINT
!          *NPMAX*  - MAXIMUM NUMBER OF PARTITIONS
!          *NPEAK*  - NUMBER OF PEAKS
!          *MIJ*    - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!          *NTHP*   - DIRECTION INDEX FOR ALL PEAKS
!          *NFRP*   - FREQUENCY INDEX FOR ALL PEAKS
!          *FL3*    - BLOCK OF SPECTRA
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

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL, NPMAX
      INTEGER, INTENT(IN), DIMENSION(IJS:IJL) :: MIJ
      INTEGER, INTENT(INOUT), DIMENSION(IJS:IJL) :: NPEAK
      INTEGER, INTENT(IN), DIMENSION(IJS:IJL,NPMAX) :: NTHP, NFRP

      REAL, INTENT(IN), DIMENSION(IJS:IJL,NANG,NFRE) :: FLLOW, FL3
      REAL, INTENT(OUT), DIMENSION(IJS:IJL,0:NPMAX) :: DIR, PER, ENE

      LOGICAL, INTENT(IN), DIMENSION(IJS:IJL,NANG) :: L LCOSDIFF

      INTEGER :: ITHC, IFRC
      INTEGER :: IJ, M, K, IP, NITT
      INTEGER :: NANGH, KK, KKMIN, KKMAX
      INTEGER :: IFRL, ITHL, ITHR
      INTEGER, DIMENSION(IJS:IJL) :: MMIN, MMAX
      INTEGER, DIMENSION(1-NANG:2*NANG) :: KLOC

      REAL :: HALF_SECTOR
      REAL :: ZHOOK_HANDLE
      REAL, DIMENSION(NANG,NFRE) :: W2
      REAL, DIMENSION(IJS:IJL,NANG,NFRE) :: W1
      REAL, DIMENSION(NANG,NFRE,NPMAX,IJS:IJL) :: SPEC

      LOGICAL :: LLCHANGE, LLADD
      LOGICAL :: LLADDPART
      LOGICAL, DIMENSION(IJS:IJL,NANG,NFRE) :: LLW3

! ----------------------------------------------------------------------

!     ASSUME THAT ONE HAS ONLY SEARCH A SECTOR IN DIRECTION (not the full 360 degrees)
!     AROUND THE PEAK
!!!   THIS IS A STRONG ASSUMPTION ON THE DIRECTIONALITY OF SWELL !!!
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('FNDPRT',0,ZHOOK_HANDLE)
#endif


      HALF_SECTOR=70.
      NANGH=NINT((HALF_SECTOR/360.)*NANG)+1

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
              W1(IJ,K,M) = 1.
              LLW3(IJ,K,M) = .FALSE.
            ELSE
              W1(IJ,K,M) = 0.
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
              W2(K,M) = 0.
            ENDDO
          ENDDO

          DO K=ITHC-1, ITHC+1
            ITHL = 1 + MOD(NANG+K-1,NANG)
            DO M=IFRC-1, IFRC+1
              IFRL = MAX(1,MIN(NFRE,M))
              IF ( W1(IJ,ITHL,IFRL) .LE. 0.5 ) W2(ITHL,IFRL) = 0.5
            ENDDO
          ENDDO

          IF ( W1(IJ,ITHC,IFRC) .LT. 0.25 ) W2(ITHC,IFRC) = 1.0

!         FIND IF MORE HIGH FREQUENCY BINS HAVE BECOME EXCLUDED
          OUT0: DO M=MMAX(IJ),MMIN(IJ),-1
            DO KK=KKMIN,KKMAX
              K=KLOC(KK)
              IF ( W1(IJ,K,M).LT.1.0 ) THEN
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
              K=KLOC(KK)

              IF ( LLW3(IJ,K,M) ) THEN
                IF ( W2(K,M).EQ.0.5 .AND. W1(IJ,K,M).LT.0.5 ) THEN
                  LLADD    = .TRUE.
                  OUT1: DO ITHR=K-1, K+1
                    ITHL = 1 + MOD(NANG+ITHR-1,NANG)
                    DO IFRL=MAX(1,M-1), MIN(NFRE,M+1)
                      IF ( W2 (ITHL,IFRL).EQ.0. .AND.
     &                   SPEC(ITHL,IFRL,IP,IJ).GT.SPEC(K,M,IP,IJ) ) THEN
                            LLADD = .FALSE.
                            EXIT OUT1
                      ENDIF
                    ENDDO
                  ENDDO OUT1
                  IF ( LLADD ) THEN
                    W2(K,M) = 1.
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

              IF ( LLW3(IJ,K,M) .AND. W1(IJ,K,M).LT.1.0 ) THEN
                IF ( W2(K,M).EQ.0. ) THEN
                  OUT2: DO ITHR=K-1, K+1
                    ITHL = 1 + MOD(NANG+ITHR-1,NANG)
                    DO IFRL=MAX(1,M-1), MIN(NFRE,M+1)
                      IF ( W2(ITHL,IFRL).EQ.1. ) THEN
                        W2(K,M) = 0.5
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
!       ADDED AS A PARTITION
        LLADDPART=.FALSE.
        IF(NPEAK(IJ).LT.NPMAX) THEN
          DO K=1,NANG
            DO M=1,NFRE
              IF(LLCOSDIFF(IJ,K) .AND. W1(IJ,K,M).LE.0.0) THEN
                W2(K,M) = 1.
                LLADDPART=.TRUE.
              ELSE
                W2(K,M) = 0.
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        IF(LLADDPART) THEN
          NPEAK(IJ)=NPEAK(IJ)+1
          IP=NPEAK(IJ)
          DO M=1,NFRE
            DO K=1,NANG
              SPEC(K,M,IP,IJ) = FL3(IJ,K,M)*W2(K,M)
            ENDDO
          ENDDO
        ENDIF

      ENDDO ! loop IJ

      CALL PARMEAN(IJS, IJL, NPMAX, NPEAK,
     &             SPEC,
     &             ENE, DIR, PER)


#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('FNDPRT',1,ZHOOK_HANDLE)
#endif

      END
