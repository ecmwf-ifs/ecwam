      SUBROUTINE TAUHF (FRMAX)

! ----------------------------------------------------------------------

!**** *TAUHF* - COMPUTATION OF HIGH-FREQUENCY STRESS.

!     PETER A.E.M. JANSSEN    KNMI      OCTOBER 90

!*    PURPOSE.
!     ---------

!       COMPUTE HIGH-FREQUENCY WAVE STRESS

!**   INTERFACE.
!     ----------

!       *CALL* *TAUHF (FRMAX)*
!          *FRMAX - LAST MODEL FREQUENCY FR(ML).

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

      USE YOWPCONS , ONLY : G        ,ZPI
      USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,ALPHA    ,XKAPPA
      USE YOWTABL  , ONLY : IUSTAR   ,IALPHA   ,USTARM   ,TAUHFT   ,
     &            DELUST   ,DELALP

! ----------------------------------------------------------------------

      PARAMETER (JTOT = 50)

! ----------------------------------------------------------------------

!*    1. PRELIMINARY CALCULATIONS.
!        -------------------------

      ALPHAM = 10.*ALPHA
      DELUST = USTARM/REAL(IUSTAR)
      DELALP = ALPHAM/REAL(IALPHA)

      CONST1 = BETAMAX/XKAPPA**2
      OMEGAC = ZPI*FRMAX

      DO L=0,IALPHA
        DO K=0,IUSTAR
          TAUHFT(K,L) = 0.
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
          DELY     = MAX((1.-YC)/REAL(JTOT),0.)
          DO J=1,JTOT
            Y        = YC+REAL(J-1)*DELY
            OMEGA    = Y*SQRT(G/Z0)
            CM       = G/OMEGA
            ZX       = UST/CM +ZALP
            ZARG     = MIN(XKAPPA/ZX,20.)
            ZMU      = MIN(G*Z0/CM**2*EXP(ZARG),1.)

            ZLOG         = MIN(LOG(ZMU),0.)
            ZBETA        = CONST1*ZMU*ZLOG**4
            TAUHFT(K,L)  = TAUHFT(K,L)+ZBETA/Y*DELY
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE TAUHF
