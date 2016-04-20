      SUBROUTINE WAMADSZIDL(ADS, ZIDL)
! ----------------------------------------------------------------------

!*    PURPOSE.
!     --------

!     RE-INITIALISES ADS AND ZIDL TO THE VALUES PROVIDED BY FIELDG
!     IF COUPLED !!!!

!**   INTERFACE.
!     ----------
!     *CALL* *WAMADSZIDL*(ADS, ZIDL)
!     *ADS*      AIR DENSITY IN KG/M3.
!     *ZIDL*     CONVECTIVE VELOCITy SCALE 

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWMAP   , ONLY : IFROMIJ  ,JFROMIJ
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NBLO
      USE YOWSTAT  , ONLY : NPROMA_WAM
      USE YOWSPEC  , ONLY : NSTART   ,NEND
      USE YOWWIND  , ONLY : FIELDG
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      REAL,DIMENSION(NINF:NSUP,NBLO), INTENT(INOUT) :: ADS, ZIDL 

      INTEGER :: IG, IJ, IX, JY
      INTEGER:: JKGLO,KIJS,KIJL,NPROMA

      REAL:: ZHOOK_HANDLE

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('WAMADSZIDL',0,ZHOOK_HANDLE)
#endif

!     KEEP CORRESPONDING CONTRIBUTION 
      IF(LWCOU) THEN
        DO IG=1,NBLO
          NPROMA=NPROMA_WAM
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ,IX,JY)
          DO JKGLO=NSTART(IRANK),NEND(IRANK),NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,NEND(IRANK))
            DO IJ=KIJS,KIJL
              IX = IFROMIJ(IJ,IG)
              JY = JFROMIJ(IJ,IG)
              ADS(IJ,IG) = FIELDG(IX,JY)%AIRD
              ZIDL(IJ,IG)= FIELDG(IX,JY)%ZIDL
            ENDDO
          ENDDO
!$OMP     END PARALLEL DO
        ENDDO
      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('WAMADSZIDL',1,ZHOOK_HANDLE)
#endif
      RETURN
      END SUBROUTINE WAMADSZIDL
