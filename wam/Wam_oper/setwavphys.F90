SUBROUTINE SETWAVPHYS 

! ----------------------------------------------------------------------

! SETS THE WIND INPUT PARAMETERS FOR THE PHYSIC CHOICE (IPHYS)
! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWALTAS , ONLY : EGRCRV   ,AGRCRV   ,BGRCRV   ,AFCRV    ,BFCRV,    &
     &                ESH      ,ASH      ,BSH      ,ASWKM    ,BSWKM
USE YOWCOUP  , ONLY : LLGCBZ0
USE YOWPHYS  , ONLY : BETAMAX  ,ZALP     ,ALPHA    ,  ALPHAPMAX,  &
     &                TAUWSHELTER, TAILFACTOR, TAILFACTOR_PM,     &
     &                SWELLF7  ,SSDSC2
USE YOWSTAT  , ONLY : IPHYS
USE YOWTEST  , ONLY : IU06
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"

REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SETWAVPHYS',0,ZHOOK_HANDLE)

      IF (IPHYS.EQ.0) THEN
!       JANSSSEN WIND INPUT PHYSICS:
        ALPHA   = 0.0060_JWRB
        BETAMAX = 1.30_JWRB
        ZALP    = 0.008_JWRB
        ALPHAPMAX = 0.031_JWRB
        TAUWSHELTER = 0.0_JWRB
        IF(LLGCBZ0) THEN
          TAILFACTOR = 2.0_JWRB
        ELSE 
          TAILFACTOR = 2.5_JWRB
        ENDIF
        TAILFACTOR_PM = 0.0_JWRB   ! i.e. not used

!!!     EMPIRICAL CONSTANCE FOR  SPECTRAL UPDATE FOLLOWING DATA ASSIMILATION
        EGRCRV = 1108.0_JWRB
        AGRCRV = 0.06E+6_JWRB
        BGRCRV = 9.70_JWRB
        AFCRV = 4.0E-4_JWRB
        BFCRV = -3.0_JWRB
        ESH = 1711.0_JWRB
        ASH = 8.0E-4_JWRB 
        BSH = 0.96_JWRB 
        ASWKM = 0.0981_JWRB
        BSWKM = 0.425_JWRB

      ELSE IF (IPHYS.EQ.1) THEN
!       ARDHUIN ET AL. (2010) WIND INPUT PHYSICS
        ALPHA   = 0.0065_JWRB
        ZALP    = 0.008_JWRB
        IF(LLGCBZ0) THEN
          ALPHAPMAX = 0.03_JWRB
          TAUWSHELTER = 0.0_JWRB
          BETAMAX = 1.36_JWRB
          TAILFACTOR = 2.4_JWRB
          TAILFACTOR_PM = 0.0_JWRB

          SWELLF7 = 4.14E05_JWRB

          SSDSC2 = -2.0E-5_JWRB
        ELSE 
          ALPHAPMAX = 0.031_JWRB
          TAUWSHELTER = 0.25_JWRB
          BETAMAX = 1.40_JWRB
          TAILFACTOR = 2.5_JWRB
          TAILFACTOR_PM = 3.0_JWRB

          SWELLF7 = 3.6E05_JWRB

          SSDSC2 = -2.2E-5_JWRB
        ENDIF

!!!     EMPIRICAL CONSTANCE FOR  SPECTRAL UPDATE FOLLOWING DATA ASSIMILATION
        EGRCRV = 1065.0_JWRB
        AGRCRV = 0.0655E+6_JWRB
        BGRCRV =  10.906_JWRB
        AFCRV = 2.453E-4_JWRB
        BFCRV = -3.1236_JWRB
        ESH = 1711.0_JWRB
        ASH = 8.0E-4_JWRB 
        BSH = 0.96_JWRB 
        ASWKM = 0.0981_JWRB
        BSWKM = 0.425_JWRB

      ELSE
        WRITE (IU06,*) '*************************************'
        WRITE (IU06,*) '*                                   *'
        WRITE (IU06,*) '*  ERROR IN SETWAVPHYS                 *'
        WRITE (IU06,*) '*  UKNOWN PHYSICS SELECTION :       *'
        WRITE (IU06,*) '*  IPHYS =' , IPHYS
        WRITE (IU06,*) '*                                   *'
        WRITE (IU06,*) '*************************************'
        CALL ABORT1
      ENDIF

IF (LHOOK) CALL DR_HOOK('SETWAVPHYS',1,ZHOOK_HANDLE)

END SUBROUTINE SETWAVPHYS
