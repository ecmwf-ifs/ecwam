      SUBROUTINE SCOSFL (F, IJS, IJL, MM, MEANCOSFL)

! ----------------------------------------------------------------------

!**** *SCOSFL* - COMPUTATION OF THE CENTERED FOURIER COEFFICIENT FOR
!                A GIVEN FREQUENCY AT EACH GRID POINT.

!     J.-R. BIDLOT        ECMWF       MARCH 2000

!*    PURPOSE.
!     --------

!     TO COMPUTE THE CENTERED FOURIER COEFFICIENT FOR A GIVEN FREQUENCY 

!**   INTERFACE.
!     ----------

!       *CALL* *SCOSFL (F, IJS, IJL, MM, MEANCOSFL)*
!          *F*  - SPECTRUM.
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *MM*  - FREQUENCY INDEX (array)
!          *MEANCOSFL* - CENTERED FOURIER COEFFICIENT

!     METHOD.
!     -------

!       INTEGRATION OF SPECTRUM TIMES COS(THETA-MEAN THETA)
!       OVER DIRECTION.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWFRED  , ONLY : DELTH    ,TH       ,COSTH    ,SINTH
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL), INTENT(IN) :: MM
      REAL(KIND=JWRB), DIMENSION(IJS:IJL, NANG, NFRE), INTENT(IN) :: F
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: MEANCOSFL


      INTEGER(KIND=JWIM) :: IJ, K
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: MEANDIR, SI, CI

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SCOSFL',0,ZHOOK_HANDLE)

!*    1. INITIALISE ARRAYS
!        -----------------

      DO IJ=IJS,IJL
        SI(IJ) = 0.0_JWRB
        CI(IJ) = 0.0_JWRB
        MEANCOSFL(IJ) = 0.0_JWRB
      ENDDO

! ----------------------------------------------------------------------

!*    2. FIND MEAN DIRECTION
!        -------------------

      DO K=1,NANG
        DO IJ=IJS,IJL
          SI(IJ) = SI(IJ)+SINTH(K)*F(IJ,K,MM(IJ))
          CI(IJ) = CI(IJ)+COSTH(K)*F(IJ,K,MM(IJ))
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        IF (CI(IJ).EQ.0.0_JWRB .AND. SI(IJ).EQ.0.0_JWRB) THEN
          MEANDIR(IJ) = 0.0_JWRB
        ELSE
          MEANDIR(IJ) = ATAN2(SI(IJ),CI(IJ))
        ENDIF
      ENDDO

!*    3. COMPUTE FOURRIER COEFFICIENT 
!        -----------------------------

      DO K=1,NANG
        DO IJ=IJS,IJL
          MEANCOSFL(IJ)=MEANCOSFL(IJ)+COS(TH(K)-MEANDIR(IJ))*F(IJ,K,MM(IJ)) 
        ENDDO
      ENDDO
      DO IJ=IJS,IJL
        MEANCOSFL(IJ) = DELTH*MEANCOSFL(IJ)
      ENDDO

      IF (LHOOK) CALL DR_HOOK('SCOSFL',1,ZHOOK_HANDLE)

      END SUBROUTINE SCOSFL
