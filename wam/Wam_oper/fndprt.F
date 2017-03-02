      SUBROUTINE FNDPRT ( ITHC, IFRC, MIJ, MMIN, MMAX,
     &                    SPEC, LLW3, W1, W2 )

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

!       *CALL* *FNDPRT (ITHC, IFRC, MIJ, MMIN, MMAX, SPEC, LLW3, W1, W2)*
!          *ITCH* - PEAK DISCRETE DIRECTION
!          *IFRC* - PEAK DISCRETE FREQUENCY
!          *MIJ*  - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!          *MMIN* - FREQUENCY INDEX BELOW WHICH SPECTRAL BINS ARE ALREADY EXCLUDED
!          *MMAX* - LAST FREQUENCY INDEX OF WHERE SPECTRAL BINS ARE ALREADY EXCLUDED
!          *SPEC* - BLOCK OF SPECTRA
!          *LLW3* - MAP OF BINS THAT ARE ALWAYS EXCLUDED (WINDSEA)
!          *W1*   - OVERALL MAP OF BINS 
!          *W2*   - MAP OF BINS

!     METHOD.
!     -------

!       STEEPEST ASCENT

!     EXTERNALS.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ITHC, IFRC, MIJ, MMIN
      INTEGER, INTENT(INOUT) :: MMAX
      REAL, INTENT(IN), DIMENSION(NANG,NFRE) :: SPEC
      LOGICAL, INTENT(IN), DIMENSION(NANG,NFRE) :: LLW3
      REAL, INTENT(INOUT), DIMENSION(NANG,NFRE) :: W1, W2

      INTEGER :: M, K, NITT
      INTEGER :: NANGH, KK, KKMIN, KKMAX
      INTEGER :: IFRL, ITHL, ITHR, IFRH
      INTEGER, DIMENSION(1-NANG:2*NANG) :: KLOC
      REAL :: HALF_SECTOR
      LOGICAL :: LLCHANGE, LLADD
      REAL :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

!     ASSUME THAT ONE HAS ONLY SEARCH A SECTOR IN DIRECTION (not the full 360 degrees)
!     AROUND THE PEAK
!!!   THIS IS A STRONG ASSUMPTION ON THE DIRECTIONALITY OF SWELL !!!
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('FNDPRT',0,ZHOOK_HANDLE)
#endif
      HALF_SECTOR=70.
      NANGH=NINT((HALF_SECTOR/360.)*NANG)+1
      KKMIN=ITHC-NANGH
      KKMAX=ITHC+NANGH
      DO KK=KKMIN,KKMAX
        KLOC(KK)=1+MOD(NANG+KK-1,NANG)
      ENDDO

!*    1.  SET UP THE W2 MAP
!         -----------------

      DO M=1, NFRE
        DO K=1, NANG
          W2(K,M) = 0.
        ENDDO
      ENDDO

      IFRL   = MAX ( 1 , IFRC-1 )
      IFRH   = MIN ( NFRE , IFRC+1 )

      DO K=ITHC-1, ITHC+1
        ITHL = 1 + MOD(NANG+K-1,NANG)
        DO M=IFRC-1, IFRC+1
          IFRL = MAX(1,MIN(NFRE,M))
          IF ( W1(ITHL,IFRL) .LE. 0.5 ) W2(ITHL,IFRL) = 0.5
        ENDDO
      ENDDO

      IF ( W1(ITHC,IFRC) .LT. 0.25 ) W2(ITHC,IFRC) = 1.0

!     FIND IF MORE HIGH FREQUENCY BINS HAVE BECOME EXCLUDED
      OUT0: DO M=MMAX,MMIN,-1
        DO KK=KKMIN,KKMAX
          K=KLOC(KK)
          IF ( W1(K,M).LT.1.0 ) THEN
            MMAX=M
            EXIT OUT0
          ENDIF
        ENDDO
      ENDDO OUT0

!*    2.  ITTERATE SEARCH
!         ---------------

      NITT=0

!     2.a Branch point

  200 CONTINUE
      NITT=NITT+1
      LLCHANGE=.FALSE.

!     2.b Determine central points

      DO M=MMIN,MIN(MIJ,MMAX)
!       by definition bins beyond M=MIJ are never extremas
!       and bins above MMAX are excluded.
        DO KK=KKMIN,KKMAX
          K=KLOC(KK)

          IF ( LLW3(K,M) ) THEN
            IF ( W2(K,M).EQ.0.5 .AND. W1(K,M).LT.0.5 ) THEN
              LLADD    = .TRUE.
              OUT1: DO ITHR=K-1, K+1
                ITHL = 1 + MOD(NANG+ITHR-1,NANG)
                DO IFRL=MAX(1,M-1), MIN(NFRE,M+1)
                  IF ( W2 (ITHL,IFRL).EQ.0. .AND.
     &                SPEC(ITHL,IFRL).GT.SPEC(K,M) ) THEN
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

!     2.c Determine peripherical points

      DO M=MMIN,MMAX
        DO KK=KKMIN,KKMAX
          K=KLOC(KK)

          IF ( LLW3(K,M) .AND. W1(K,M).LT.1.0 ) THEN
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

!*    2.d Branch back ?

      IF ( LLCHANGE .AND. NITT.LT.25 ) GOTO 200


!*    3   UPDATE THE OVERALL MAP

      DO M=1, NFRE
        DO K=1, NANG
          W1(K,M) = W1(K,M) + W2(K,M)
        ENDDO
      ENDDO
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('FNDPRT',1,ZHOOK_HANDLE)
#endif


      RETURN
      END
