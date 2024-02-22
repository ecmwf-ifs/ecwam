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
&                      WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT, &
&                      DELLAM1_EXT, COSPHM1_EXT,              &
&                      DEPTH_EXT, U_EXT, V_EXT )

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

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "mpexchng.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: NINF, NSUP ! HALO EXTEND NINF to NSUP+1

      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NFRE_RED), INTENT(OUT) :: WAVNUM_EXT  ! WAVE NUMBER
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NFRE_RED), INTENT(OUT) :: CGROUP_EXT  ! GROUP VELOCITY
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NFRE_RED), INTENT(OUT) :: OMOSNH2KD_EXT ! OMEGA / SINH(2KD)
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(OUT) :: DELLAM1_EXT ! 1/DELLA
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(OUT) :: COSPHM1_EXT ! 1/COSPH
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(OUT) :: DEPTH_EXT ! WATER DEPTH
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(OUT) :: U_EXT ! U-COMPONENT OF SURFACE CURRENT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(OUT) :: V_EXT ! V-COMPONENT OF SURFACE CURRENT

      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NFRE, NCHNK), INTENT(IN) :: WAVNUM, CGROUP, OMOSNH2KD
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NCHNK), INTENT(IN) :: DEPTH, DELLAM1, COSPHM1, UCUR, VCUR

      INTEGER(KIND=JWIM) :: IJ, M
      INTEGER(KIND=JWIM) :: ICHNK, KIJS, KIJL, IJSB, IJLB

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PROENVHALO',0,ZHOOK_HANDLE)
!$acc data present(WAVNUM,CGROUP,OMOSNH2KD,DELLAM1,COSPHM1,DEPTH, UCUR,VCUR,&
!$acc WAVNUM_EXT,CGROUP_EXT,OMOSNH2KD_EXT,DELLAM1_EXT,COSPHM1_EXT,DEPTH_EXT,U_EXT,V_EXT)

!!! mapping chuncks to block ONLY for actual grid points !!!!
#ifndef _OPENACC
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, IJSB, KIJL, IJLB, M)
#endif /*_OPENACC*/
!$acc kernels loop private(ICHNK, KIJS, IJSB, KIJL, IJLB)
      DO ICHNK = 1, NCHNK
        KIJS = 1
        IJSB = IJFROMCHNK(KIJS, ICHNK)
        KIJL = KIJL4CHNK(ICHNK)
        IJLB = IJFROMCHNK(KIJL, ICHNK)

!$acc loop        
        DO M = 1, NFRE_RED
          WAVNUM_EXT(IJSB:IJLB, M) = WAVNUM(KIJS:KIJL, M,ICHNK)
          CGROUP_EXT(IJSB:IJLB, M) = CGROUP(KIJS:KIJL, M,ICHNK)
          OMOSNH2KD_EXT(IJSB:IJLB, M) = OMOSNH2KD(KIJS:KIJL, M,ICHNK)
        ENDDO

        DELLAM1_EXT(IJSB:IJLB) = DELLAM1(KIJS:KIJL,ICHNK)
        COSPHM1_EXT(IJSB:IJLB) = COSPHM1(KIJS:KIJL,ICHNK)
        DEPTH_EXT(IJSB:IJLB) = DEPTH(KIJS:KIJL,ICHNK)
        U_EXT(IJSB:IJLB) = UCUR(KIJS:KIJL,ICHNK)
        V_EXT(IJSB:IJLB) = VCUR(KIJS:KIJL,ICHNK)
      ENDDO
!$acc end kernels
#ifndef _OPENACC
!$OMP END PARALLEL DO
#endif /*_OPENACC*/

!$acc enter data copyin(WVPRPT_LAND)
!$acc enter data copyin(WVPRPT_LAND%WAVNUM,WVPRPT_LAND%CGROUP,WVPRPT_LAND%OMOSNH2KD)
!$acc data present(WVPRPT_LAND) copyin(BATHYMAX)
!!    should be combined into one single data exchange, when we start using this option.... !!!
      CALL MPEXCHNG(WAVNUM_EXT, NFRE_RED, 1, 1)
      !$acc kernels
      WAVNUM_EXT(NSUP+1,1:NFRE_RED) = WVPRPT_LAND%WAVNUM(1:NFRE_RED)
      !$acc end kernels

      CALL MPEXCHNG(CGROUP_EXT, NFRE_RED, 1, 1)
      !$acc kernels
      CGROUP_EXT(NSUP+1,1:NFRE_RED) =  WVPRPT_LAND%CGROUP(1:NFRE_RED)
      !$acc end kernels

      CALL MPEXCHNG(OMOSNH2KD_EXT, NFRE_RED, 1, 1)
      !$acc kernels
      OMOSNH2KD_EXT(NSUP+1,1:NFRE_RED) = WVPRPT_LAND%OMOSNH2KD(1:NFRE_RED)
      !$acc end kernels

      CALL MPEXCHNG(DELLAM1_EXT, 1, 1, 1)
      !$acc kernels
      DELLAM1_EXT(NSUP+1) = 0.0_JWRB 
      !$acc end kernels

      CALL MPEXCHNG(COSPHM1_EXT, 1, 1, 1)
      !$acc kernels
      COSPHM1_EXT(NSUP+1) = 0.0_JWRB 
      !$acc end kernels

      CALL MPEXCHNG(DEPTH_EXT, 1, 1, 1)
      !$acc kernels
      DEPTH_EXT(NSUP+1) = BATHYMAX
      !$acc end kernels

      CALL MPEXCHNG(U_EXT, 1, 1, 1)
      !$acc kernels
      U_EXT(NSUP+1) = 0.0_JWRB
      !$acc end kernels

      CALL MPEXCHNG(V_EXT, 1, 1, 1)
      !$acc kernels
      V_EXT(NSUP+1) = 0.0_JWRB
      !$acc end kernels
!$acc end data
!$acc exit data delete(WVPRPT_LAND%WAVNUM,WVPRPT_LAND%CGROUP,WVPRPT_LAND%OMOSNH2KD)
!$acc exit data delete(WVPRPT_LAND)
!$acc end data

IF (LHOOK) CALL DR_HOOK('PROENVHALO',1,ZHOOK_HANDLE)

END SUBROUTINE PROENVHALO
