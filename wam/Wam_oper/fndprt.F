      SUBROUTINE FNDPRT ( ITHC, IFRC, SPEC, W1, W2 )

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

!       *CALL* *FNDPRT (ITHC, IFRC, SPEC, W1, W2)*
!          *ITCH* - PEAK DISCRETE DIRECTION
!          *IFRC* - PEAK DISCRETE FREQUENCY
!          *SPEC* - BLOCK OF SPECTRA
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

! ----------------------------------------------------------------------

      INTEGER :: IFRC, ITHC, IFRL, ITHL, IFRH, ITHH
      REAL    :: SPEC(NANG,NFRE), W1(NANG,NFRE), W2(NANG,NFRE)
      LOGICAL :: CHANGE, ADD

! ----------------------------------------------------------------------


!*    1.  SET UP THE W2 MAP
!         -----------------

      DO M=1, NFRE
        DO K=1, NANG
          W2(K,M) = 0.
        ENDDO
      ENDDO

      IFRL   = MAX ( 1 , IFRC-1 )
      IFRH   = MIN ( NFRE , IFRC+1 )
      ITHL   = 1 + MOD(NANG+ITHC-2,NANG)
      ITHH   = 1 + MOD(ITHC,NANG)

      DO K=ITHC-1, ITHC+1
        ITHL = 1 + MOD(NANG+K-1,NANG)
        DO M=IFRC-1, IFRC+1
          IFRL = MAX(1,MIN(NFRE,M))
          IF ( W1(ITHL,IFRL) .LE. 0.5 ) W2(ITHL,IFRL) = 0.5
        ENDDO
      ENDDO

      IF ( W1(ITHC,IFRC) .LT. 0.25 ) W2(ITHC,IFRC) = 1.0


!*    2.  ITTERATE SEARCH
!         ---------------

      NITT   = 0

!     2.a Branch point

  200 CONTINUE
      NITT   = NITT + 1
      CHANGE = .FALSE.

!     2.b Determine central points

      DO M=1, NFRE
        DO K=1, NANG

          IF ( W2(K,M).EQ.0.5 .AND. W1(K,M).LT.0.5 ) THEN
              ADD    = .TRUE.
              DO ITHR=K-1, K+1
                ITHL = 1 + MOD(NANG+ITHR-1,NANG)
                DO IFRL=MAX(1,M-1), MIN(NFRE,M+1)
                  IF ( W2 (ITHL,IFRL).EQ.0. .AND.
     &                SPEC(ITHL,IFRL).GT.SPEC(K,M) ) ADD = .FALSE.
                ENDDO
              ENDDO
              IF ( ADD ) W2(K,M) = 1.
              CHANGE = CHANGE .OR. ADD
            ENDIF

        ENDDO
      ENDDO

!     2.c Determine periferical points

      DO M=1, NFRE
        DO K=1, NANG

          IF ( W2(K,M).EQ.0. ) THEN
              ADD    = .FALSE.
              DO ITHR=K-1, K+1
                ITHL = 1 + MOD(NANG+ITHR-1,NANG)
                DO IFRL=MAX(1,M-1), MIN(NFRE,M+1)
                  IF ( W2 (ITHL,IFRL).EQ.1. ) ADD = .TRUE.
                ENDDO
              ENDDO
              IF ( ADD ) W2(K,M) = 0.5
              CHANGE = CHANGE .OR. ADD
            ENDIF

        ENDDO
      ENDDO

!*    2.d Branch back ?

      IF ( CHANGE .AND. NITT.LT.25 ) GOTO 200


!*    3   UPDATE THE OVERALL MAP

      DO M=1, NFRE
        DO K=1, NANG
          W1(K,M) = W1(K,M) + W2(K,M)
        ENDDO
      ENDDO


      RETURN
      END
