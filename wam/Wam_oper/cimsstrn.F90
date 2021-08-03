      SUBROUTINE CIMSSTRN(IJS, IJL, KIJS, KIJL, GFL, WAVNUM, DEPTH, STRN)

! ----------------------------------------------------------------------

!**** *CIMSSTRN* - COMPUTATION OF THE MEAN SQUARE WAVE STRAIN IN SEA ICE.

!     J. BIDLOT  ECMWF  JANUARY 2013.

!*    PURPOSE.
!     --------

!       COMPUTES MEAN SQUARE WAVE STRAIN AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *CIMSSTRN (IJS, IJL, KIJS, KIJL, GFL, WAVNUM, DEPTH, STRN)*
!              *IJS:IJL* - 1st DIMENSION OF GFL
!              *KIJS*    - INDEX OF FIRST GRIDPOINT
!              *KIJL*    - INDEX OF LAST GRIDPOINT
!              *GFL*     - SPECTRUM.
!              *WAVNUM*  - OPEN WATER WAVE NUMBER
!              *DEPTH*   - WATER DEPTH
!              *STRN*    - MEAN SQUARE WAVE STRAIN IN ICE (OUTPUT).

!     METHOD.
!     -------

!      !!! IT ASSUMES SO DEFAULT SETTING FOR THE MECHANICAL PROPERTIES OF
!          THE SEA ICE (SEE AKI_ICE) !!!!!!!

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWICE   , ONLY : CITHICK  ,FLMIN
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: GFL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE), INTENT(IN) :: WAVNUM

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: STRN


      INTEGER(KIND=JWIM) :: IJ,M,K,JD
      REAL(KIND=JWRB) :: AKI_ICE
      REAL(KIND=JWRB) :: F1LIM 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: XKI, E, SUME 
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('CIMSSTRN',0,ZHOOK_HANDLE)

!*    1. INITIALISE
!        ----------

      F1LIM=FLMIN/DELTH

      DO IJ=KIJS,KIJL
        STRN(IJ) = 0.0_JWRB
      ENDDO

! ----------------------------------------------------------------------

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      DO M=1,NFRE
        DO IJ=KIJS,KIJL
          XKI(IJ)=AKI_ICE(G,WAVNUM(IJ,M),DEPTH(IJ),ROWATER,CITHICK(IJ))
          E(IJ)=0.5_JWRB*CITHICK(IJ)*XKI(IJ)**3/WAVNUM(IJ,M)
        ENDDO

        DO IJ=KIJS,KIJL
          SUME(IJ) = 0.0_JWRB
        ENDDO
        DO K=1,NANG
          DO IJ=KIJS,KIJL
            SUME(IJ) = SUME(IJ)+GFL(IJ,K,M)
          ENDDO
        ENDDO

        DO IJ=KIJS,KIJL
          IF(SUME(IJ).GT.F1LIM) THEN
            STRN(IJ) = STRN(IJ)+E(IJ)**2*SUME(IJ)*DFIM(M)
          ENDIF
        ENDDO

      ENDDO

      IF (LHOOK) CALL DR_HOOK('CIMSSTRN',1,ZHOOK_HANDLE)

      END SUBROUTINE CIMSSTRN
