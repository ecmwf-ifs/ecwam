      SUBROUTINE TAUHF(ML)

! ----------------------------------------------------------------------

!**** *TAUHF* - COMPUTATION OF HIGH-FREQUENCY STRESS.
!                              HIGH-FREQUENCY ENERGY FLUX.

!     PETER A.E.M. JANSSEN    KNMI      OCTOBER 90

!*    PURPOSE.
!     ---------

!       COMPUTE HIGH-FREQUENCY WAVE STRESS AND ENERGY FLUX
!       FOR BOTH ECMWF PHYSICS AND METREO FRANCE PHYSICS.

!**   INTERFACE.
!     ----------

!       *CALL* *TAUHF(ML)*
!             *ML*  NUMBER OF FREQUENCIES.

!     METHOD.
!     -------

!       SEE REFERENCE FOR WAVE STRESS CALCULATION.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : ALPHA    ,TAUWSHELTER, ITSHELT
      USE YOWPCONS , ONLY : G
      USE YOWTABL  , ONLY : IUSTAR   ,IALPHA   ,USTARM   ,TAUHFT   ,
     &            DELUST   ,DELALP   ,ALPHAMCOEF,
     &            XLEVTAILM, ILEVTAIL ,TAUHFT2 ,DELTAIL  ,
     &            PHIHFT   ,PHIHFT2

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ML

      INTEGER, PARAMETER :: JTOT=19  ! must be odd !!! 
      INTEGER :: M, L, K, I, J
      INTEGER, DIMENSION(0:IUSTAR) :: MIJ

      REAL :: ALPHAM, GM1
      REAL, DIMENSION(0:IUSTAR) :: USTAR, Z0, XLEVTAIL, TAU, PHI
      REAL, DIMENSION(0:IALPHA) :: CHARN

! ----------------------------------------------------------------------

!*    1. PRELIMINARY CALCULATIONS.
!        -------------------------

      GM1= 1./G

      DELUST = USTARM/REAL(IUSTAR)
      DO K=0,IUSTAR
        USTAR(K) = MAX(REAL(K)*DELUST,0.000001)
      ENDDO

      ALPHAM = ALPHAMCOEF*ALPHA
      DELALP = ALPHAM/REAL(IALPHA)
      DO L=0,IALPHA
        CHARN(L) = ALPHA+REAL(L)*DELALP
      ENDDO

      DELTAIL= XLEVTAILM/REAL(ILEVTAIL)

      IF (ABS(TAUWSHELTER).LE.0.0) THEN
        ITSHELT=0
      ELSE
        ITSHELT=1
      ENDIF

      CALL INIT_X0TAUHF

      IF (ITSHELT.EQ.0) THEN

!     TABLES FOR ECMWF PHYSICS:
!     -------------------------

      IF(.NOT.ALLOCATED(TAUHFT)) ALLOCATE(TAUHFT(0:IUSTAR,0:IALPHA,ML))
      IF(.NOT.ALLOCATED(PHIHFT)) ALLOCATE(PHIHFT(0:IUSTAR,0:IALPHA,ML))

!*    2. CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
!        ------------------------------------------------

      DO K=0,IUSTAR
        XLEVTAIL(K) = 0.
      ENDDO
!$OMP PARALLEL DO PRIVATE(M,K,L,MIJ,Z0,TAU,PHI)
      DO M=1,ML
        DO L=0,IALPHA
          DO K=0,IUSTAR
            MIJ(K) = M
            Z0(K) = USTAR(K)**2*CHARN(L)*GM1
          ENDDO
          CALL TAU_PHI_HF(0, IUSTAR, MIJ, USTAR, Z0, XLEVTAIL,
     &                    TAU, PHI)
          DO K=0,IUSTAR
            TAUHFT(K,L,M)=TAU(K)
            PHIHFT(K,L,M)=PHI(K)
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ELSE

!     TABLES FOR METEO FRANCE PHYSICS:
!     -------------------------------

      IF(.NOT.ALLOCATED(TAUHFT2)) 
     & ALLOCATE(TAUHFT2(0:IUSTAR,0:IALPHA,0:ILEVTAIL))
      IF(.NOT.ALLOCATED(PHIHFT2)) 
     & ALLOCATE(PHIHFT2(0:IUSTAR,0:IALPHA,0:ILEVTAIL))

!*    3. CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
!        ------------------------------------------------

      DO K=0,IUSTAR
        MIJ(K) = ML
      ENDDO
!$OMP PARALLEL DO PRIVATE(I,K,L,Z0,XLEVTAIL,TAU,PHI)
      DO I=0,ILEVTAIL
        DO L=0,IALPHA
          DO K=0,IUSTAR
            XLEVTAIL(K) = REAL(I)*DELTAIL
            Z0(K) = USTAR(K)**2*CHARN(L)*GM1
          ENDDO
          CALL TAU_PHI_HF(0, IUSTAR, MIJ, USTAR, Z0, XLEVTAIL,
     &                    TAU, PHI)
          DO K=0,IUSTAR
            TAUHFT2(K,L,I)=TAU(K)
            PHIHFT2(K,L,I)=PHI(K)
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ENDIF ! ITSHELT

      END SUBROUTINE TAUHF
