! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE CIREDUCE_LOKI_GPU (WVPRPT, FF_NOW)

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
      USE YOWDRVTYPE ,ONLY: FREQUENCY, FORCING_FIELDS
      USE YOWSTAT, ONLY: LUPDATE_GPU_GLOBALS

! ----------------------------------------------------------------------
      IMPLICIT NONE

#include "ciwaf.intfb.h"

      TYPE(FREQUENCY), INTENT(INOUT)            :: WVPRPT
      TYPE(FORCING_FIELDS), INTENT(IN)          :: FF_NOW


      INTEGER(KIND=JWIM) :: IJ, M 
      INTEGER(KIND=JWIM) :: ICHNK

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE, ZHOOK_HANDLE_DATA_OFFLOAD

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
#ifdef OMPGPU
!$omp target teams distribute parallel do collapse(3) map(to:WVPRPT)
#else
!$acc kernels present(WVPRPT)
#endif
            DO ICHNK = 1, NCHNK
               DO M = 1, NFRE
                 DO IJ = 1, NPROMA_WAM
                   WVPRPT%CIWA(IJ,M,ICHNK) = 1.0_JWRB
                 ENDDO
               ENDDO
            ENDDO
#ifdef OMPGPU
!$omp end target teams distribute parallel do
#else
!$acc end kernels
#endif
            CALL GSTATS(1493,1)
          ENDIF

        ELSE

IF(LUPDATE_GPU_GLOBALS)THEN
IF (LHOOK) CALL DR_HOOK('DATA_OFFLOAD',0,ZHOOK_HANDLE_DATA_OFFLOAD)
!$loki update_device
IF (LHOOK) CALL DR_HOOK('DATA_OFFLOAD',1,ZHOOK_HANDLE_DATA_OFFLOAD)
ENDIF
          CALL GSTATS(1493,0)
!         DETERMINE THE WAVE ATTENUATION FACTOR

!$loki structured-data present(FF_NOW, WVPRPT)

          DO ICHNK = 1, NCHNK
            CALL CIWAF(1, NPROMA_WAM, WVPRPT%CGROUP(:,:,ICHNK), FF_NOW%CICOVER(:,ICHNK), &
&                      FF_NOW%CITHICK(:,ICHNK), WVPRPT%CIWA(:,:,ICHNK))
          ENDDO

!$loki end structured-data

          CALL GSTATS(1493,1)
        ENDIF

IF (LHOOK) CALL DR_HOOK('CIREDUCE',1,ZHOOK_HANDLE)

END SUBROUTINE CIREDUCE_LOKI_GPU
