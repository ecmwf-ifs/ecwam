      SUBROUTINE JONSWAP (ALPHAJ, ZGAMMA, SA, SB, FP, IJS, IJL, ET)

! ---------------------------------------------------------------------

!**** *JONSWAP* - ROUTINE TO COMPUTE THE 1-D JONSWAP SPECTRUM.

!     

!*    PURPOSE.
!     --------

!     TO COMPUTE THE 1-D JONSWAP SPECTRUM. 

!**   INTERFACE.
!     ----------

!       *CALL* *JONSWAP (ALPHAJ, ZGAMMA, SA, SB, FP, IJS, IJL, ET)*
!        *ALPHAJ*  -  OVERALL ENERGY LEVEL OF JONSWAP SPECTRUM
!        *ZGAMMA*   -  OVERSHOOT FACTOR.
!        *SA*      -  LEFT PEAK WIDTH.
!        *SB*      -  RIGHT PEAK WIDTH.
!        *FP*      -  PEAK FREQUENCY.
!        *IJS*     -  FIRST POINT IN BLOCK.
!        *IJL*     -  LAST  POINT IN BLOCK.
!        *ET*      -  JOSWAP SPECTRUM AT DISCRETISED FRQUENCY BINS.


!     METHOD.
!     -------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR
      USE YOWPARAM , ONLY : NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,ZPI4GM2
! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: M, IJ, IJS, IJL

      REAL(KIND=JWRB) :: ZGAMMA, SA, SB, SIGMA, G2ZPI4FRH5M
      REAL(KIND=JWRB) :: FRH, EARG, FJON, FMPF, FJONH

      REAL(KIND=JWRB) :: ALPHAJ(IJS:IJL), FP(IJS:IJL)
      REAL(KIND=JWRB) :: ET(IJS:IJL,NFRE)

! ----------------------------------------------------------------------

      DO M=1,NFRE
        FRH = FR(M)
        G2ZPI4FRH5M=1.0_JWRB/(FRH**5*ZPI4GM2)
        DO IJ=IJS,IJL
          IF (ALPHAJ(IJ) /= 0.0_JWRB .AND. FP(IJ) /= 0.0_JWRB) THEN
            IF (FRH > FP(IJ)) THEN
              SIGMA = SB
            ELSE
              SIGMA = SA
            ENDIF
            EARG = 0.5_JWRB*((FRH-FP(IJ)) / (SIGMA*FP(IJ)))**2
            EARG = MIN(EARG,50.0_JWRB)
            FJON = ZGAMMA**EXP(-EARG)
            FMPF = 1.25_JWRB*(FP(IJ)/FRH)**4
            FMPF = MIN(FMPF,50.0_JWRB)
            FJONH = EXP(-FMPF)
            ET(IJ,M) = ALPHAJ(IJ)*G2ZPI4FRH5M*FJONH*FJON
          ELSE
            ET(IJ,M) = 0.0_JWRB
          ENDIF
        ENDDO
      ENDDO

      END SUBROUTINE JONSWAP
