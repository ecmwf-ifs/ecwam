      SUBROUTINE JONSWAP (ALPHAJ, GAMMA, SA, SB, FP, IJS, IJL, ET)

! ---------------------------------------------------------------------

!**** *JONSWAP* - ROUTINE TO COMPUTE THE 1-D JONSWAP SPECTRUM.

!     

!*    PURPOSE.
!     --------

!     TO COMPUTE THE 1-D JONSWAP SPECTRUM. 

!**   INTERFACE.
!     ----------

!       *CALL* *JONSWAP (ALPHAJ, GAMMA, SA, SB, FP, IJS, IJL, ET)*
!        *ALPHAJ*  -  OVERALL ENERGY LEVEL OF JONSWAP SPECTRUM
!        *GAMMA*   -  OVERSHOOT FACTOR.
!        *SA*      -  LEFT PEAK WIDTH.
!        *SB*      -  RIGHT PEAK WIDTH.
!        *FP*      -  PEAK FREQUENCY.
!        *IJS*     -  FIRST POINT IN BLOCK.
!        *IJL*     -  LAST  POINT IN BLOCK.
!        *ET*      -  JOSWAP SPECTRUM AT DISCRETISED FRQUENCY BINS.


!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.


!     REFERENCES.
!     -----------

!       NONE.

! ----------------------------------------------------------------------
      USE YOWFRED  , ONLY : FR
      USE YOWPARAM , ONLY : NFRE
      USE YOWPCONS , ONLY : G        ,ZPI
! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: M, IJ, IJS, IJL

      REAL :: GAMMA, SA, SB, SIGMA, G2ZPI4, G2ZPI4FRH5M
      REAL :: FRH, EARG, FJON, FMPF, FJONH

      REAL :: ALPHAJ(IJS:IJL), FP(IJS:IJL)
      REAL :: ET(IJS:IJL,NFRE)

! ----------------------------------------------------------------------

      G2ZPI4=G**2/ZPI**4

      DO M=1,NFRE
        FRH = FR(M)
        G2ZPI4FRH5M=G2ZPI4/FRH**5
        DO IJ=IJS,IJL
          IF (FP(IJ).NE.0.) THEN
            IF (FRH.GT.FP(IJ)) THEN
              SIGMA = SB
            ELSE
              SIGMA = SA
            ENDIF
            EARG = .5*((FRH-FP(IJ)) / (SIGMA*FP(IJ)))**2
            EARG = MIN(EARG,20.)
            FJON = GAMMA**EXP(-EARG)
            FMPF = 1.25*(FP(IJ)/FRH)**4
            FMPF = MIN(FMPF,20.)
            FJONH = EXP(-FMPF)
            ET(IJ,M) = ALPHAJ(IJ)*G2ZPI4FRH5M*FJONH*FJON
          ELSE
            ET(IJ,M) = 0.
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE JONSWAP
