      SUBROUTINE MEANSQS_LF(NFRE_EFF, IJS, IJL, KIJS, KIJL, F, XMSS)

! ----------------------------------------------------------------------

!**** *MEANSQS_LF* - COMPUTATION OF RESOLVED MEAN SQUARE SLOPE UP TO WAVE FREQUENCY FR(NFRE_EFF)

!*    PURPOSE.
!     --------

!       COMPUTE MEAN SQUARE SLOPE AT EACH GRID POINT FOR THE RESOLVED SPECTRUM UP TO FREQUENCY FR(NFRE_EFF).
!       NFRE_EFF <= NFRE

!**   INTERFACE.
!     ----------

!       *CALL* *MEANSQS_LF (NFRE_EFF, IJS, IJL, KIJS, KIJL, F, XMSS)*
!              *XKMSS*   - WAVE NUMBER CUT OFF
!              *IJS:IJL* - 1st DIMENSION of F
!              *KIJS*    - INDEX OF FIRST GRIDPOINT
!              *KIJL*    - INDEX OF LAST GRIDPOINT
!              *F*       - SPECTRUM.
!              *XMSS*    - MEAN SQUARE SLOPE (OUTPUT).

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

      USE YOWPCONS , ONLY : G        ,GM1      ,ZPI
      USE YOWFRED  , ONLY : FR       ,ZPIFR    ,DFIM
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSHAL  , ONLY : TFAK     ,INDEP
      USE YOWSTAT  , ONLY : ISHALLO

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: NFRE_EFF
      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: XMSS 


      INTEGER(KIND=JWIM) :: IJ, M, K, KFRE

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NFRE) :: FD
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: TEMP1, TEMP2

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('MEANSQS_LF',0,ZHOOK_HANDLE)

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      KFRE=MIN(NFRE_EFF, NFRE)

      XMSS(:) = 0.0_JWRB

      IF (ISHALLO.EQ.1) THEN

!*    2.1 DEEP WATER INTEGRATION.
!         -----------------------

        DO M = 1, KFRE
          FD(M) = DFIM(M)*(ZPIFR(M))**4*GM1**2
        ENDDO

        DO M = 1, KFRE
          DO IJ = KIJS, KIJL
            TEMP2(IJ) = 0.0_JWRB
          ENDDO
          DO K = 1, NANG
            DO IJ = KIJS, KIJL
              TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
            ENDDO
          ENDDO
          DO IJ = KIJS, KIJL
            XMSS(IJ) = XMSS(IJ)+FD(M)*TEMP2(IJ)
          ENDDO
        ENDDO
!SHALLOW
      ELSE

!*    2.2 SHALLOW WATER INTEGRATION.
!         --------------------------

        DO M = 1, KFRE
          DO IJ = KIJS, KIJL
            TEMP1(IJ) = DFIM(M)*TFAK(INDEP(IJ),M)**2
            TEMP2(IJ) = 0.0_JWRB
          ENDDO
          DO K = 1, NANG
            DO IJ = KIJS, KIJL
              TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
            ENDDO
          ENDDO
          DO IJ = KIJS, KIJL
            XMSS(IJ) = XMSS(IJ)+TEMP1(IJ)*TEMP2(IJ)
          ENDDO
        ENDDO
      ENDIF
!SHALLOW

      IF (LHOOK) CALL DR_HOOK('MEANSQS_LF',1,ZHOOK_HANDLE)

      END SUBROUTINE MEANSQS_LF
