      SUBROUTINE JONSWAP (ALPHAJ, GAMMA, SA, SB, FP, IJS, IJL, ET)

! ----------------------------------------------------------------------

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

      REAL :: GAMMA, SA, SB
      REAL :: FRH, EARG, FJON, FMPF, FJONH, PIRHF, F, ALF 

      REAL :: SIGMA(IJS:IJL,NFRE) 
      REAL :: ALPHAJ(IJS:IJL), FP(IJS:IJL)
      REAL :: ET(IJS:IJL,NFRE)

! ----------------------------------------------------------------------
      PIRHF(F,ALF)=ALF*G**2/(ZPI**4*F**5)

! ----------------------------------------------------------------------

      DO M=1,NFRE
        FRH = FR(M)
        DO IJ=IJS,IJL
          SIGMA(IJ,M) = SA
          IF (FRH.GT.FP(IJ)) SIGMA(IJ,M) = SB
        ENDDO
      ENDDO

      DO M=1,NFRE
        FRH = FR(M)
        DO IJ=IJS,IJL
          IF (FP(IJ).NE.0.) THEN
            EARG = .5*((FRH-FP(IJ)) / (SIGMA(IJ,M)*FP(IJ)))**2
            EARG = MIN(EARG,20.)
            FJON = GAMMA**EXP(-EARG)
            FMPF = 1.25*(FP(IJ)/FRH)**4
            FMPF = MIN(FMPF,20.)
            FJONH = EXP(-FMPF)
            ET(IJ,M) = PIRHF(FRH,ALPHAJ(IJ))*FJONH*FJON
          ELSE
            ET(IJ,M) = 0.
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE JONSWAP 
