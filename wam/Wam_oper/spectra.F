      SUBROUTINE SPECTRA (IJS, IJL, IG, FL1)

! ----------------------------------------------------------------------

!**** *SPECTRA*  - COMPUTATION OF 2-D SPECTRA FOR ONE BLOCK.

!     S. HASSELMANN  - JULY 1990
!     H. GUNTHER     - DECEMBER 1990   MODIFIED FOR CYCLE_4.
!     J. BIDLOT      - FEBRUARY 1996   MESSAGE PASSING

!*    PURPOSE.
!     --------

!       INITIALISATION OF A BLOCK BY 2-D SPECTRA.

!**   INTERFACE.
!     ----------

!       *CALL* *SPECTRA (IJS, IJL, IG, FL1)*
!          *IJS*     INTEGER  FIRST POINT IN BLOCK.
!          *IJL*     INTEGER  LAST  POINT IN BLOCK.
!          *IG*      INTEGER  BLOCK NUMBER.
!          *FL1*      REAL      2-D SPECTRUM FOR EACH GRID POINT

!     METHOD.
!     -------
!       1-D JONSWAP SPECTRA AND COSINE**2 SPREADING FUNCTIONS ARE
!       COMPUTED FROM GIVEN WINDS AND PARAMETERS AT EACH GRID POINT.
!       THE 1-D SPECTRA ARE SPREAD OVER THE DIRECTIONS BY MULTIPLICATION
!       WITH THE SPREADING FUNCTION.

!     EXTERNALS.
!     ----------

!       *SPT*       - COMPUTATION OF COS**2 SPREADING FUNCTION.


!     REFERENCES.
!     -----------

!       K.HASSELMAN,D.B.ROSS,P.MUELLER AND W.SWELL
!          A PARAMETRIC WAVE PREDICTION MODEL
!          JOURNAL OF PHYSICAL OCEANOGRAPHY, VOL. 6, NO. 2, MARCH 1976

! ----------------------------------------------------------------------

      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWJONS  , ONLY : FP       ,ALPHJ    ,THES     ,GAMMA    ,
     &            SA       ,SB
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI

! ----------------------------------------------------------------------

      REAL, ALLOCATABLE :: ST(:,:)
      REAL, ALLOCATABLE :: ET(:,:)
      REAL, ALLOCATABLE :: STH(:)

!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL1

! ----------------------------------------------------------------------
      PIRHF(F,ALF)=ALF*G**2/(ZPI**4*F**5)

! ----------------------------------------------------------------------

      ALLOCATE(ST(IJS:IJL,NANG))
      ALLOCATE(ET(IJS:IJL,NFRE))
      ALLOCATE(STH(NANG))

!*    1. COMPUTE JONSWAP SPECTRUM.
!        -------------------------

      DO M=1,NFRE
        DO IJ=IJS,IJL
          IF (FP(IJ).NE.0.) THEN
            FRH = FR(M)
            SIGMA = SA
            IF (FRH.GT.FP(IJ)) SIGMA = SB
            EARG = .5*((FRH-FP(IJ)) / (SIGMA*FP(IJ)))**2
            FJON = 1
            IF (EARG.LT.99.) FJON = GAMMA**EXP(-EARG)
            FMPF = 1.25*(FP(IJ)/FRH)**4
            FJONH = 0.
            IF (FMPF.LT.99.) FJONH = EXP(-FMPF)
            ET(IJ,M) = PIRHF(FRH,ALPHJ(IJ))*FJONH*FJON
          ELSE
            ET(IJ,M) = 0.
          ENDIF
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    2. COMPUTATION OF SPREADING FUNCTION.
!        ----------------------------------

      DO IJ=IJS,IJL
        CALL SPR (NANG, NANG, THES(IJ), TH, STH)
        DO K=1,NANG
          ST(IJ,K) = STH(K)
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    3. COMPUTATION OF 2-D SPECTRUM.
!        ----------------------------

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            FL1(IJ,K,M) = ET(IJ,M)*ST(IJ,K)
          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(ST)
      DEALLOCATE(ET)
      DEALLOCATE(STH)

      RETURN
      END SUBROUTINE SPECTRA
