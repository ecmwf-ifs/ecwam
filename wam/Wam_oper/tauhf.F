      SUBROUTINE TAUHF(ML)

! ----------------------------------------------------------------------

!**** *TAUHF* - COMPUTATION OF HIGH-FREQUENCY STRESS.

!     PETER A.E.M. JANSSEN    KNMI      OCTOBER 90

!*    PURPOSE.
!     ---------

!       COMPUTE HIGH-FREQUENCY WAVE STRESS

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

      USE YOWFRED  , ONLY : FR
      USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,ALPHA    ,XKAPPA
      USE YOWPCONS , ONLY : G        ,ZPI
      USE YOWTABL  , ONLY : IUSTAR   ,IALPHA   ,USTARM   ,TAUHFT   ,
     &            DELUST   ,DELALP

! ----------------------------------------------------------------------

      INTEGER, PARAMETER :: JTOT=250

      REAL, ALLOCATABLE :: W(:)

! ----------------------------------------------------------------------

!*    1. PRELIMINARY CALCULATIONS.
!        -------------------------

      ALPHAMCOEF = 40.
      ALPHAM = ALPHAMCOEF*ALPHA
      DELUST = USTARM/REAL(IUSTAR)
      DELALP = ALPHAM/REAL(IALPHA)

      CONST1 = BETAMAX/XKAPPA**2

      ALLOCATE(W(JTOT))
      W=1.
      W(1)=0.5
      W(JTOT)=0.5

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

      DEALLOCATE(W)

      RETURN
      END SUBROUTINE TAUHF
