      SUBROUTINE TAUHF(ML)

! ----------------------------------------------------------------------

!**** *TAUHF* - COMPUTATION OF HIGH-FREQUENCY STRESS.

!     PETER A.E.M. JANSSEN    KNMI      OCTOBER 90

!*    PURPOSE.
!     ---------

!       COMPUTE HIGH-FREQUENCY WAVE STRESS
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

      USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,ALPHA    ,XKAPPA   ,
     &             TAUWSHELTER, ITSHELT
      USE YOWFRED  , ONLY : FR
      USE YOWPCONS , ONLY : G        ,ZPI
      USE YOWTABL  , ONLY : IUSTAR   ,IALPHA   ,USTARM   ,TAUHFT   ,
     &            DELUST   ,DELALP   ,ALPHAMCOEF,
     &            XLEVTAILM, ILEVTAIL ,TAUHFT2 ,DELTAIL
      USE YOWSTAT  , ONLY : IPHYS

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, PARAMETER :: JTOT=250

      INTEGER, INTENT(IN) :: ML
      INTEGER :: M, L, K, I, J
      REAL :: ALPHAM, CONST1, OMEGA, OMEGAC, OMEGACC
      REAL :: X0, UST, Z0, UST0, TAUW, TAUW0, XLEVTAIL
      REAL :: YC,DELY,Y,CM,ZX,ZARG,ZMU,ZLOG,ZBETA
      REAL, ALLOCATABLE :: W(:)

! ----------------------------------------------------------------------

!*    1. PRELIMINARY CALCULATIONS.
!        -------------------------

      DELUST = USTARM/REAL(IUSTAR)

      ALPHAM = ALPHAMCOEF*ALPHA
      DELALP = ALPHAM/REAL(IALPHA)

      DELTAIL= XLEVTAILM/REAL(ILEVTAIL)

      CONST1 = BETAMAX/XKAPPA**2

      ALLOCATE(W(JTOT))
      W=1.
      W(1)=0.5
      W(JTOT)=0.5

      IF (ABS(TAUWSHELTER).LE.0.0) THEN
        ITSHELT=0
      ELSE
        ITSHELT=1
      ENDIF

      IF (ITSHELT.EQ.0) THEN

!     TABLES FOR ECMWF PHYSICS:
!     -------------------------

      IF(.NOT.ALLOCATED(TAUHFT)) ALLOCATE(TAUHFT(0:IUSTAR,0:IALPHA,ML))

!$OMP PARALLEL DO PRIVATE(M,L,K,J,OMEGAC,X0,UST,Z0,OMEGACC,
!$OMP+ YC,DELY,Y,OMEGA,CM,ZX,ZARG,ZMU,ZLOG,ZBETA)
      DO M=1,ML

        OMEGAC = ZPI*FR(M)

        DO L=0,IALPHA
          DO K=0,IUSTAR
            TAUHFT(K,L,M) = 0.
          ENDDO
        ENDDO

!*    2. CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
!        ------------------------------------------------

        X0 = 0.05
        DO L=0,IALPHA
          DO K=0,IUSTAR
            UST      = MAX(REAL(K)*DELUST,0.000001)
            Z0       = UST**2*(ALPHA+REAL(L)*DELALP)/G
            OMEGACC  = MAX(OMEGAC,X0*G/UST)
            YC       = OMEGACC*SQRT(Z0/G)
            DELY     = MAX((1.-YC)/REAL(JTOT-1),0.)
            DO J=1,JTOT
              Y        = YC+REAL(J-1)*DELY
              OMEGA    = Y*SQRT(G/Z0)
              CM       = G/OMEGA
              ZX       = UST/CM +ZALP
              ZARG     = MIN(XKAPPA/ZX,20.)
              ZMU      = MIN(G*Z0/CM**2*EXP(ZARG),1.)

              ZLOG         = MIN(LOG(ZMU),0.)
              ZBETA        = CONST1*ZMU*ZLOG**4
              TAUHFT(K,L,M)= TAUHFT(K,L,M)+W(J)*ZBETA/Y*DELY
            ENDDO
          ENDDO
        ENDDO

      ENDDO
!$OMP END PARALLEL DO

      ELSE

!     TABLES FOR METEO FRANCE PHYSICS:
!     -------------------------------

      IF(.NOT.ALLOCATED(TAUHFT2)) 
     & ALLOCATE(TAUHFT2(0:IUSTAR,0:IALPHA,0:ILEVTAIL))

!*    3. CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
!        ------------------------------------------------

!$OMP PARALLEL DO PRIVATE(L,K,I,J,OMEGAC,X0,UST0,UST,Z0,OMEGACC,
!$OMP+ YC,DELY,XLEVTAIL,TAUW0,TAUW,Y,OMEGA,CM,ZX,ZARG,ZMU,ZLOG,ZBETA)
      DO L=0,IALPHA
        OMEGAC = ZPI*FR(ML)

        X0 = 0.05
        DO K=0,IUSTAR
          UST0     = MAX(REAL(K)*DELUST,0.000001)
          UST      = UST0
          Z0       = UST0**2*(ALPHA+FLOAT(L)*DELALP)/G
          OMEGACC  = MAX(OMEGAC,X0*G/UST)
          YC       = OMEGACC*SQRT(Z0/G)
          DELY     = MAX((1.-YC)/REAL(JTOT),0.)
          DO I=0,ILEVTAIL
            TAUHFT2(K,L,I)=0.
            XLEVTAIL=REAL(I)*DELTAIL
            TAUW0=UST0**2
            TAUW=TAUW0
            DO J=1,JTOT
              Y        = YC+REAL(J-1)*DELY
              OMEGA    = Y*SQRT(G/Z0)
              CM       = G/OMEGA
              ZX       = UST/CM +ZALP
              ZARG     = MIN(XKAPPA/ZX,20.)
              ZMU      = MIN(G*Z0/CM**2*EXP(ZARG),1.)

              ZLOG         = MIN(ALOG(ZMU),0.)
              ZBETA        = CONST1*ZMU*ZLOG**4
              TAUHFT2(K,L,I)=
     &           TAUHFT2(K,L,I)+W(J)*ZBETA*(UST/UST0)**2/Y*DELY
              TAUW         =TAUW-W(J)*UST**2*ZBETA*XLEVTAIL/Y*DELY
              UST          =SQRT(MAX(TAUW,0.))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ENDIF ! ITSHELT

      DEALLOCATE(W)

      RETURN
      END SUBROUTINE TAUHF
