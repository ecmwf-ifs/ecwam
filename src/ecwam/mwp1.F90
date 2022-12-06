      SUBROUTINE MWP1 (KIJS, KIJL, F, MEANWP1)

! ----------------------------------------------------------------------

!**** *MWP1* - COMPUTATION OF MEAN PERIOD BASED ON THE FIRST MOMENT
!              FOR EACH GRID POINT.

!     J-R BIDLOT    ECMWF     MARCH 2000

!*    PURPOSE.
!     --------

!       COMPUTE MEAN PERIOD 1 AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *MWP1 (KIJS, KIJL, F, MEANWP1)*
!              *KIJS*    - INDEX OF FIRST GRIDPOINT
!              *KIJL*    - INDEX OF LAST GRIDPOINT
!              *F*       - SPECTRUM.
!              *MEANWP1* - MEAN PERIOD BASED ON THE FIRST MOMENT.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,DFIM_SIM ,DFIMFR_SIM,DELTH  ,WETAIL, WP1TAIL
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_ODD
      USE YOWPCONS , ONLY : EPSMIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), INTENT(IN) :: F(KIJS:KIJL,NANG,NFRE)
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: MEANWP1


      INTEGER(KIND=JWIM) :: IJ, K, M
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: DELT25, COEF_FR, FR1M1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TEMP, EM

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('MWP1',0,ZHOOK_HANDLE)

      DO IJ=KIJS,KIJL
        EM(IJ) = 0.0_JWRB
        MEANWP1(IJ) = 0.0_JWRB
      ENDDO

      DO M=1,NFRE_ODD
        DO IJ=KIJS,KIJL
          TEMP(IJ) = 0.0_JWRB
        ENDDO
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            TEMP(IJ) = TEMP(IJ)+F(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=KIJS,KIJL
          EM(IJ) = EM(IJ)+DFIM_SIM(M)*TEMP(IJ)
          MEANWP1(IJ) = MEANWP1(IJ)+DFIMFR_SIM(M)*TEMP(IJ)
        ENDDO
      ENDDO

!*    ADD TAIL CORRECTION TO MEAN PERIOD AND
!     AND TRANSFORM TO PERIOD.

      FR1M1 = 1.0_JWRB/FR(1)
      DELT25 = WETAIL*FR(NFRE_ODD)*DELTH
      COEF_FR = WP1TAIL*DELTH*FR(NFRE_ODD)**2

      DO IJ=KIJS,KIJL
        EM(IJ) = EM(IJ)+DELT25*TEMP(IJ)
        MEANWP1(IJ) = MEANWP1(IJ)+COEF_FR*TEMP(IJ)
        IF(EM(IJ) > 0.0_JWRB .AND. MEANWP1(IJ) > EPSMIN ) THEN
          MEANWP1(IJ) = EM(IJ)/MEANWP1(IJ)
          MEANWP1(IJ) = MIN(MEANWP1(IJ),FR1M1)
        ELSE
          MEANWP1(IJ) = 0.0_JWRB
        ENDIF
      ENDDO

      IF (LHOOK) CALL DR_HOOK('MWP1',1,ZHOOK_HANDLE)

      END SUBROUTINE MWP1