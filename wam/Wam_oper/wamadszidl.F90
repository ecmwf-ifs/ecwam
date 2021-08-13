      SUBROUTINE WAMADSZIDL(IJS, IJL, ADS, ZIDL)
! ----------------------------------------------------------------------

!*    PURPOSE.
!     --------

!     RE-INITIALISES ADS AND ZIDL TO THE VALUES PROVIDED BY FIELDG
!     IF COUPLED !!!!

!**   INTERFACE.
!     ----------
!     *CALL* *WAMADSZIDL*(ADS, ZIDL)
!      *IJS:IJL - FIRST DIMENSION OF ARRAYS ADS and ZIDL
!     *ADS*      AIR DENSITY IN KG/M3.
!     *ZIDL*     CONVECTIVE VELOCITy SCALE 

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWMAP   , ONLY : IFROMIJ  ,JFROMIJ
      USE YOWSTAT  , ONLY : NPROMA_WAM
      USE YOWWIND  , ONLY : FIELDG
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(INOUT):: ADS, ZIDL


      INTEGER(KIND=JWIM) :: IJ, IX, JY
      INTEGER(KIND=JWIM):: JKGLO, KIJS, KIJL, NPROMA

      REAL(KIND=JWRB):: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WAMADSZIDL',0,ZHOOK_HANDLE)

!     KEEP CORRESPONDING CONTRIBUTION 
      IF(LWCOU) THEN
          NPROMA=NPROMA_WAM
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ,IX,JY)
          DO JKGLO=IJS,IJL,NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,IJL)
            DO IJ=KIJS,KIJL
              IX = IFROMIJ(IJ)
              JY = JFROMIJ(IJ)
              ADS(IJ) = FIELDG(IX,JY)%AIRD
              ZIDL(IJ)= FIELDG(IX,JY)%ZIDL
            ENDDO
          ENDDO
!$OMP     END PARALLEL DO
      ENDIF

      IF (LHOOK) CALL DR_HOOK('WAMADSZIDL',1,ZHOOK_HANDLE)
 
      END SUBROUTINE WAMADSZIDL
