      SUBROUTINE CIREDUCE (IJS, IJL, CGROUP, CICOVER, CITHICK, CIWA)

! ----------------------------------------------------------------------

!**** *CIREDUCE* - COMPUTE SEA ICE REDUCTION FACTOR FOR SOURCE TERMS 
!                  AND THE SEA ICE WAVE ATTENUATION FACTORS

!           IF THERE IS NO SEA ICE INFORMATION OR
!           ALL SEA ICE COVER POINTS WILL BE MASKED
!           THEN CIWA WILL BE SET ON THE FIRST CALL. NOTHING WILL BE DONE
!           IN ALL FOLLOWING CALLS

!!!! currently also setting parametric sea ice thickness !!!!

!*    PURPOSE.
!     --------

!       CIREDUCE COMPUTES SEA ICE SOURCE TERM REDUCTION FACTOR.

!**   INTERFACE.
!     ----------

!       *CALL* *CIREDUCE (IJS, IJL, CGROUP, CICOVER, CITHICK, CIWA)

!          *IJS*     - INDEX OF FIRST GRIDPOINT
!          *IJL*     - INDEX OF LAST GRIDPOINT
!          *CGROUP*  - GROUP SPEED.
!          *CICOVER* - SEA ICE COVER.
!          *CITHICK* - SEA ICE THICKNESS. 
!          *CIWA*-     SEA ICE WAVE ATTENUATION FACTOR. 

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWICE   , ONLY : LICERUN  ,LMASKICE 
      USE YOWPARAM , ONLY : NFRE
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWSTAT  , ONLY : NPROMA_WAM

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "ciwaf.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB),DIMENSION(IJS:IJL,NFRE), INTENT(IN) :: CGROUP
      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(IN) :: CICOVER
      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(IN) :: CITHICK
      REAL(KIND=JWRB),DIMENSION(IJS:IJL,NFRE), INTENT(INOUT) :: CIWA


      INTEGER(KIND=JWIM) :: IJ, JH, M 
      INTEGER(KIND=JWIM) :: JKGLO,KIJS,KIJL,NPROMA

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      LOGICAL, SAVE :: LLFRST

      DATA LLFRST / .TRUE. /

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('CIREDUCE',0,ZHOOK_HANDLE)


        IF( .NOT. LICERUN .OR. LMASKICE ) THEN
          IF(LLFRST) THEN
            LLFRST=.FALSE.
!           NO REDUCTION, EITHER THERE IS NO SEA ICE INFORMATION OR
!           ALL SEA ICE COVER POINTS WILL BE MASKED
            NPROMA=NPROMA_WAM
            CALL GSTATS(1493,0)
!$OMP       PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ,M) 
            DO JKGLO=IJS,IJL,NPROMA
              KIJS=JKGLO
              KIJL=MIN(KIJS+NPROMA-1,IJL)
              DO M=1,NFRE
                DO IJ=KIJS,KIJL
                  CIWA(IJ,M)=1.0_JWRB
                ENDDO
              ENDDO
            ENDDO
!$OMP       END PARALLEL DO
            CALL GSTATS(1493,1)
          ENDIF
        ELSE

          CALL GSTATS(1493,0)
!         DETERMINE THE WAVE ATTENUATION FACTOR
          NPROMA=NPROMA_WAM
!$OMP     PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
          DO JKGLO=IJS,IJL,NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,IJL)
            CALL CIWAF(KIJS, KIJL, CGROUP(KIJS:KIJL,:), CICOVER(KIJS), CITHICK(KIJS), CIWA(KIJS:KIJL,:))
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1493,1)
        ENDIF

      IF (LHOOK) CALL DR_HOOK('CIREDUCE',1,ZHOOK_HANDLE)

      END SUBROUTINE CIREDUCE
