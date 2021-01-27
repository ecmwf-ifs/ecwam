SUBROUTINE SETWAVPHYS 

! ----------------------------------------------------------------------

! SETS THE WIND INPUT PARAMETERS FOR THE PHYSIC CHOICE (IPHYS)
! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWALTAS , ONLY : EGRCRV   ,AGRCRV   ,BGRCRV   ,AFCRV    ,BFCRV,    &
     &                ESH      ,ASH      ,BSH      ,ASWKM    ,BSWKM
USE YOWCOUP  , ONLY : LLGCBZ0  ,LLNORMAGAM
USE YOWPHYS  , ONLY : BETAMAX  ,ZALP     ,ALPHAMIN ,ALPHA    ,ALPHAPMAX,&
     &                TAUWSHELTER, TAILFACTOR, TAILFACTOR_PM,           &
     &                DELTA_THETA_RN, RN1_RN, DTHRN_A, DTHRN_U,         &
     &                ANG_GC_A, ANG_GC_B, ANG_GC_C, ANG_GC_D, ANG_GC_E
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
        ZALP    = 0.008_JWRB
        ALPHAPMAX = 0.03_JWRB
        TAUWSHELTER = 0.0_JWRB
        TAILFACTOR = 2.5_JWRB

!       DIRECTIONALITY CORRECTION FACTORS IN THE GOWTH RATE RENORMALISATION 
        DELTA_THETA_RN = 0.75_JWRB
        RN1_RN = 0.25_JWRB
        DTHRN_A = 0.0_JWRB
        DTHRN_U = 33.0_JWRB

        IF(LLGCBZ0) THEN
          ALPHAMIN = 0.0005_JWRB
          ALPHA   = 0.0065_JWRB

          IF(LLNORMAGAM) THEN
            BETAMAX = 1.3_JWRB
          ELSE
            BETAMAX = 1.22_JWRB
          ENDIF

          ! ANGULAR ADJUSTMENT PARAMETERS FOR THE GRAVITY-CAPILLARY MODEL
          ANG_GC_A = 0.62_JWRB
           !!! currently no adjustment
          ANG_GC_B = 0.00_JWRB
          ANG_GC_C = 0.40_JWRB
          ANG_GC_D = 11.5_JWRB
          ANG_GC_E = 0.1_JWRB

        ELSE 
          ALPHAMIN = 0.0001_JWRB
          ALPHA   = 0.0065_JWRB
          BETAMAX = 1.20_JWRB
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
        ZALP    = 0.008_JWRB
        TAILFACTOR = 2.5_JWRB

!       directionality correction factors in the gowth rate renormalisation 
        DELTA_THETA_RN = 0.25_JWRB
        RN1_RN = 0.25_JWRB
        DTHRN_A = 1.0_JWRB
        DTHRN_U = 33.0_JWRB

        IF(LLGCBZ0) THEN
          ALPHAMIN = 0.0005_JWRB
          ALPHA   = 0.0065_JWRB
          ALPHAPMAX = 0.030_JWRB
          TAUWSHELTER = 0.25_JWRB
          TAILFACTOR_PM = 0.0_JWRB

          ! ANGULAR ADJUSTMENT PARAMETERS FOR THE GRAVITY-CAPILLARY MODEL
          ANG_GC_A = 0.50_JWRB
          ANG_GC_B = 0.20_JWRB
          ANG_GC_C = 0.40_JWRB
          ANG_GC_D = 11.5_JWRB
          ANG_GC_E = 0.1_JWRB

          IF(LLNORMAGAM) THEN
            BETAMAX = 1.44_JWRB
          ELSE
            BETAMAX = 1.44_JWRB
          ENDIF

        ELSE 
          ALPHAMIN = 0.0001_JWRB
          ALPHA   = 0.0065_JWRB
          ALPHAPMAX = 0.031_JWRB
          TAILFACTOR_PM = 3.0_JWRB
          IF(LLNORMAGAM) THEN
            BETAMAX = 1.40_JWRB
            TAUWSHELTER = 0.0_JWRB
          ELSE
            BETAMAX = 1.40_JWRB
            TAUWSHELTER = 0.25_JWRB
          ENDIF
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
