! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE PROENVHALO (NINF, NSUP,                            &
&                      WAVNUM, CGROUP, OMOSNH2KD,            &
&                      DEPTH, DELLAM1, COSPHM1, UCUR, VCUR,   &
&                      BUFFER_EXT)

! ----------------------------------------------------------------------

!**** *PROENVHALO* - WAVE PROPGATION

!*    PURPOSE.
!     --------

!     PRODUCES ARRAYS WITH GRID POINTS VALUES AND THEIR HALO
!     FOR ENVIRONMENT VARIABLES FOR THE WAVE PROPGATION 

! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : WVPRPT_LAND
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, KIJL4CHNK, IJFROMCHNK
      USE YOWPARAM , ONLY : NFRE     , NFRE_RED
      USE YOWSHAL  , ONLY : BATHYMAX

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE YOWDRVTYPE, ONLY: ENVIRONMENT, FREQUENCY
#ifdef WAM_GPU
      USE YOWDRVTYPE, ONLY: GET_DEVICE_DATA_RDONLY
#endif

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "mpexchng.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: NINF, NSUP ! HALO EXTEND NINF to NSUP+1

      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, 3*NFRE_RED + 5), INTENT(OUT) :: BUFFER_EXT ! OMEGA / SINH(2KD)

      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NFRE, NCHNK), INTENT(IN) :: WAVNUM, CGROUP, OMOSNH2KD
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NCHNK), INTENT(IN) :: DEPTH, DELLAM1, COSPHM1, UCUR, VCUR

      INTEGER(KIND=JWIM) :: IJ, M
      INTEGER(KIND=JWIM) :: ICHNK, IJB

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE, ZHOOK_HANDLE_MPI

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PROENVHALO',0,ZHOOK_HANDLE)
!$acc data present(WAVNUM,CGROUP,OMOSNH2KD,DELLAM1,COSPHM1,DEPTH,UCUR,VCUR) &
!$acc present(BUFFER_EXT)

!!! mapping chuncks to block ONLY for actual grid points !!!!
#ifdef WAM_GPU
#ifdef OMPGPU
!$omp target teams distribute parallel do collapse(3) private(IJB)
#else
!$acc parallel loop gang vector collapse(3) private(IJB)
#endif
#else
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, IJB, M, IJ)
#endif
      DO ICHNK = 1, NCHNK
!       THE EXTRA ITERATION M = NFRE_RED+1 COPIES THE FIELDS WITHOUT A FREQUENCY INDEX
        DO M = 1, NFRE_RED+1
          DO IJ = 1, NPROMA_WAM
            IF (IJ <= KIJL4CHNK(ICHNK)) THEN
              IJB = IJFROMCHNK(1, ICHNK) + IJ - 1
              IF (M <= NFRE_RED) THEN
                BUFFER_EXT(IJB, M) = WAVNUM(IJ, M, ICHNK)
                BUFFER_EXT(IJB, M + NFRE_RED) = CGROUP(IJ, M, ICHNK)
                BUFFER_EXT(IJB, M + 2*NFRE_RED) = OMOSNH2KD(IJ, M, ICHNK)
              ELSE
                BUFFER_EXT(IJB, 3*NFRE_RED+1) = DELLAM1(IJ, ICHNK)
                BUFFER_EXT(IJB, 3*NFRE_RED+2) = COSPHM1(IJ, ICHNK)
                BUFFER_EXT(IJB, 3*NFRE_RED+3) = DEPTH(IJ, ICHNK)
                BUFFER_EXT(IJB, 3*NFRE_RED+4) = UCUR(IJ, ICHNK)
                BUFFER_EXT(IJB, 3*NFRE_RED+5) = VCUR(IJ, ICHNK)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
#ifdef WAM_GPU
#ifdef OMPGPU
!$omp end target teams distribute parallel do
#else
!$acc end parallel loop
#endif
#else
!$OMP END PARALLEL DO
#endif

#ifdef WAM_GPU
      CALL GET_DEVICE_DATA_RDONLY(WVPRPT_LAND)
#endif
      IF (LHOOK) CALL DR_HOOK('MPI_TIME',0,ZHOOK_HANDLE_MPI)
      CALL MPEXCHNG(BUFFER_EXT, 3*NFRE_RED+5, 1, 1)
      IF (LHOOK) CALL DR_HOOK('MPI_TIME',1,ZHOOK_HANDLE_MPI)

#ifdef OMPGPU
      !$omp target map(to:WVPRPT_LAND)
      !$omp parallel do
#else
      !$acc kernels present(WVPRPT_LAND)
#endif
      DO M = 1, NFRE_RED
        BUFFER_EXT(NSUP+1,M) = WVPRPT_LAND%WAVNUM(M)
        BUFFER_EXT(NSUP+1,NFRE_RED+M) = WVPRPT_LAND%CGROUP(M)
        BUFFER_EXT(NSUP+1,2*NFRE_RED+M) = WVPRPT_LAND%OMOSNH2KD(M)
      ENDDO

      BUFFER_EXT(NSUP+1,3*NFRE_RED+1) = 0.0_JWRB
      BUFFER_EXT(NSUP+1,3*NFRE_RED+2) = 0.0_JWRB 
      BUFFER_EXT(NSUP+1,3*NFRE_RED+3) = BATHYMAX
      BUFFER_EXT(NSUP+1,3*NFRE_RED+4) = 0.0_JWRB 
      BUFFER_EXT(NSUP+1,3*NFRE_RED+5) = 0.0_JWRB 
#ifdef OMPGPU
      !$omp end target
#else
      !$acc end kernels
#endif

!$acc end data

IF (LHOOK) CALL DR_HOOK('PROENVHALO',1,ZHOOK_HANDLE)

END SUBROUTINE PROENVHALO
