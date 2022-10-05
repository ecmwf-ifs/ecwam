SUBROUTINE CIREDUCE (CGROUP, CICOVER, CITHICK, CIWA)

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

!       *CALL* *CIREDUCE (CGROUP, CICOVER, CITHICK, CIWA)

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

      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE 
      USE YOWPARAM , ONLY : NFRE

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

#include "ciwaf.intfb.h"

      REAL(KIND=JWRB),DIMENSION(NPROMA_WAM, NFRE, NCHNK), INTENT(IN)    :: CGROUP
      REAL(KIND=JWRB),DIMENSION(NPROMA_WAM, NCHNK), INTENT(IN)          :: CICOVER
      REAL(KIND=JWRB),DIMENSION(NPROMA_WAM, NCHNK), INTENT(IN)          :: CITHICK
      REAL(KIND=JWRB),DIMENSION(NPROMA_WAM, NFRE, NCHNK), INTENT(INOUT) :: CIWA


      INTEGER(KIND=JWIM) :: IJ, M 
      INTEGER(KIND=JWIM) :: ICHNK

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      LOGICAL, SAVE :: LLFRST

      DATA LLFRST / .TRUE. /

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CIREDUCE',0,ZHOOK_HANDLE)


        IF( .NOT. LICERUN .OR. LMASKICE ) THEN

          IF (LLFRST) THEN
            LLFRST=.FALSE.
!           NO REDUCTION, EITHER THERE IS NO SEA ICE INFORMATION OR
!           ALL SEA ICE COVER POINTS WILL BE MASKED
            CALL GSTATS(1493,0)
!$OMP       PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, M, IJ) 
            DO ICHNK = 1, NCHNK
               CIWA(:, :, ICHNK) = 1.0_JWRB
            ENDDO
!$OMP       END PARALLEL DO
            CALL GSTATS(1493,1)
          ENDIF

        ELSE

          CALL GSTATS(1493,0)
!         DETERMINE THE WAVE ATTENUATION FACTOR
!$OMP     PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ICHNK)
          DO ICHNK = 1, NCHNK
            CALL CIWAF(1, NPROMA_WAM, CGROUP(:,:,ICHNK), CICOVER(:,ICHNK), CITHICK(:,ICHNK), CIWA(:,:,ICHNK))
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1493,1)
        ENDIF

IF (LHOOK) CALL DR_HOOK('CIREDUCE',1,ZHOOK_HANDLE)

END SUBROUTINE CIREDUCE
