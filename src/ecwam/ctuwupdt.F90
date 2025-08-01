! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE CTUWUPDT (IJS, IJL, NINF, NSUP,                  &
 &                   BLK2GLO,                               &
 &                   CGROUP_EXT, OMOSNH2KD_EXT,             &
 &                   COSPHM1_EXT, DEPTH_EXT, U_EXT, V_EXT )


! ----------------------------------------------------------------------

!**** *CTUWUPDT*

!*    PURPOSE.
!     --------

!     COMPUTES THE WEIGHTS FOR THE CTU ADVECTION SCHEME (IPROPAGS=2).

! -------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
USE YOWDRVTYPE  , ONLY : WVGRIDGLO

USE YOWCURR  , ONLY : LLCFLCUROFF
USE YOWFRED  , ONLY : COSTH    ,SINTH
USE YOWGRID  , ONLY : NPROMA_WAM, COSPH
USE YOWREFD  , ONLY : THDD     ,THDC     ,SDOT
USE YOWMPP   , ONLY : IRANK    ,NPROC
USE YOWPARAM , ONLY : NANG     ,NFRE_RED
USE YOWSTAT  , ONLY : IFRELFMAX, DELPRO_LF, IDELPRO, IREFRA
USE YOWTEST  , ONLY : IU06
USE YOWUBUF  , ONLY : SUMWN    ,                                            &
&                     JXO      ,JYO      ,KCR      ,KPM      ,MPM,          &
&                     WLATN    ,WLONN    ,WCORN    ,WKPMN    ,WMPMN    ,    &
&                     LLWLATN  ,LLWLONN  ,LLWCORN  ,LLWKPMN  ,LLWMPMN  ,    &
&                     KLON, KLAT, WLAT, KCOR, WCOR
USE YOWFRED  , ONLY : FR       ,DELTH, COSTH    ,SINTH
USE YOWPCONS , ONLY : ZPI


USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE OML_MOD  , ONLY : OML_GET_MAX_THREADS
      
! ----------------------------------------------------------------------
      IMPLICIT NONE

#include "abort1.intfb.h"
#include "ctuw.intfb.h"
#include "ctuwdrv.intfb.h"
#include "ctuwini.intfb.h"

INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL   ! GRID POINTS WITHIN A BLOCK
INTEGER(KIND=JWIM), INTENT(IN) :: NINF, NSUP ! GRID POINT WITH HALO EXTEND NINF:NSUP+1 
TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO  ! BLOCK TO GRID TRANSFORMATION
REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NFRE_RED), INTENT(IN) :: CGROUP_EXT  ! GROUP VELOCITY
REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NFRE_RED), INTENT(IN) :: OMOSNH2KD_EXT ! OMEGA / SINH(2KD)
REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: COSPHM1_EXT ! 1/COSPH
REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: DEPTH_EXT ! WATER DEPTH
REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: U_EXT ! U-COMPONENT OF SURFACE CURRENT
REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1), INTENT(IN) :: V_EXT ! V-COMPONENT OF SURFACE CURRENT


INTEGER(KIND=JWIM) :: IJ, K, M, J, KM1, KP1
INTEGER(KIND=JWIM) :: IC, ICL, ICR
INTEGER(KIND=JWIM) :: MSTART, MEND
INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA, MTHREADS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JWRB) :: DELPRO
REAL(KIND=JWRB), DIMENSION(NINF:NSUP,2) :: WLATM1 ! 1 - WLAT
REAL(KIND=JWRB), DIMENSION(NINF:NSUP,4) :: WCORM1 ! 1 - WCOR
REAL(KIND=JWRB), DIMENSION(NINF:NSUP,2) :: DP     ! COS PHI FACTOR

LOGICAL :: LL2NDCALL
LOGICAL, DIMENSION(IJS:IJL) :: LCFLFAIL

LOGICAL, SAVE :: LFRSTCTU
DATA LFRSTCTU /.TRUE./

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CTUWUPDT',0,ZHOOK_HANDLE)

!$acc data create(WLATM1,WCORM1,DP)

! DEFINE JXO, JYO, KCR
IF (LFRSTCTU) THEN

  IF (.NOT. ALLOCATED(MPM)) ALLOCATE(MPM(NFRE_RED,-1:1))
  !$acc kernels
  DO M=1,NFRE_RED
    MPM(M,-1)= MAX(1,M-1)
    MPM(M,0) = M
    MPM(M,1) = MIN(NFRE_RED,M+1)
  ENDDO
  !$acc end kernels

  IF (.NOT. ALLOCATED(KPM)) ALLOCATE(KPM(NANG,-1:1))
  IF (.NOT. ALLOCATED(JXO)) ALLOCATE(JXO(NANG,2))
  IF (.NOT. ALLOCATED(JYO)) ALLOCATE(JYO(NANG,2))
  IF (.NOT. ALLOCATED(KCR)) ALLOCATE(KCR(NANG,4))

!$acc update device(JXO, JYO, KCR, KPM)
 !$acc kernels
  DO K=1,NANG

    KM1 = K-1
    IF (KM1 < 1) KM1 = NANG
    KPM(K,-1)=KM1

    KPM(K,0)=K

    KP1 = K+1
    IF (KP1 > NANG) KP1 = 1
    KPM(K,1)=KP1


    IF (COSTH(K) >= 0.0_JWRB) THEN
      JYO(K,1)=1
      JYO(K,2)=2
      IF (SINTH(K) >= 0.0_JWRB) THEN
        JXO(K,1)=1
        JXO(K,2)=2
        KCR(K,1)=3
        KCR(K,2)=2
        KCR(K,3)=4
        KCR(K,4)=1
      ELSE
        JXO(K,1)=2
        JXO(K,2)=1
        KCR(K,1)=2
        KCR(K,2)=3
        KCR(K,3)=1
        KCR(K,4)=4
      ENDIF
    ELSE
      JYO(K,1)=2
      JYO(K,2)=1
      IF (SINTH(K) >= 0.0_JWRB) THEN
        JXO(K,1)=1
        JXO(K,2)=2
        KCR(K,1)=4
        KCR(K,2)=1
        KCR(K,3)=3
        KCR(K,4)=2
      ELSE
        JXO(K,1)=2
        JXO(K,2)=1
        KCR(K,1)=1
        KCR(K,2)=4
        KCR(K,3)=2
        KCR(K,4)=3
      ENDIF
    ENDIF
  ENDDO
  !$acc end kernels

  LFRSTCTU = .FALSE.

ENDIF


! THE CTU IS USED, COMPUTE THE WEIGHTS

IF (.NOT. ALLOCATED(SUMWN)) ALLOCATE(SUMWN(IJS:IJL,NANG,NFRE_RED))
IF (.NOT. ALLOCATED(WLATN)) ALLOCATE(WLATN(IJS:IJL,NANG,NFRE_RED,2,2))
IF (.NOT. ALLOCATED(WLONN)) ALLOCATE(WLONN(IJS:IJL,NANG,NFRE_RED,2))
IF (.NOT. ALLOCATED(WCORN)) ALLOCATE(WCORN(IJS:IJL,NANG,NFRE_RED,4,2))
IF (.NOT. ALLOCATED(WKPMN)) ALLOCATE(WKPMN(IJS:IJL,NANG,NFRE_RED,-1:1))

IF (IREFRA == 2 .OR. IREFRA == 3) THEN
  IF (.NOT. ALLOCATED(WMPMN)) ALLOCATE(WMPMN(IJS:IJL,NANG,NFRE_RED,-1:1))

#ifndef _OPENACC
  IF (.NOT. ALLOCATED(LLWLATN)) ALLOCATE(LLWLATN(NANG,NFRE_RED,2,2))
  IF (.NOT. ALLOCATED(LLWLONN))  ALLOCATE(LLWLONN(NANG,NFRE_RED,2))
  IF (.NOT. ALLOCATED(LLWCORN)) ALLOCATE(LLWCORN(NANG,NFRE_RED,4,2))
  IF (.NOT. ALLOCATED(LLWKPMN)) ALLOCATE(LLWKPMN(NANG,NFRE_RED,-1:1))
  IF (.NOT. ALLOCATED(LLWMPMN)) ALLOCATE(LLWMPMN(NANG,NFRE_RED,-1:1))
#endif
ENDIF



! SOME INITIALISATION FOR *CTUW*
!! NPROMA=NPROMA_WAM
   MTHREADS=1
#ifndef _OPENACC
   MTHREADS=OML_GET_MAX_THREADS()
#endif /*_OPENACC*/
   NPROMA=(IJL-IJS+1)/MTHREADS + 1

#ifdef _OPENACC
!$acc data present(KLAT,WLAT,KCOR,WCOR,WLATN,WLONN,WCORN)
#else
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO, KIJS, KIJL)
#endif /*_OPENACC*/
DO JKGLO = IJS, IJL, NPROMA
  KIJS=JKGLO
  KIJL=MIN(KIJS+NPROMA-1,IJL)
  CALL CTUWINI (KIJS, KIJL, NINF, NSUP, BLK2GLO, COSPHM1_EXT,   &
 &                  WLATM1, WCORM1, DP)
ENDDO
#ifdef _OPENACC
!$acc end data
#else
!$OMP  END PARALLEL DO
#endif /*_OPENACC*/


! COMPUTES THE WEIGHTS
! --------------------

IF (IFRELFMAX <= 0 ) THEN
  ! *PROPAGS2* IS NOT SPLIT, COMPUTE THE WEIGHTS FOR ALL WAVES 
  DELPRO = REAL(IDELPRO, JWRB)
  MSTART = 1 
  MEND   = NFRE_RED
ELSE
  ! *PROPAGS2* IS SPLIT, COMPUTE THE WEIGHTS FOR THE FAST WAVES
  DELPRO = DELPRO_LF
  MSTART = 1 
  MEND   = IFRELFMAX 
ENDIF


CALL CTUWDRV (DELPRO, MSTART, MEND,                 &
 &           IJS, IJL, NINF, NSUP,                  &
 &           BLK2GLO,                               &
 &           WLATM1, WCORM1, DP,                    &
 &           CGROUP_EXT, OMOSNH2KD_EXT,             &
 &           COSPHM1_EXT, DEPTH_EXT, U_EXT, V_EXT )


! IF *PROPAGS* IS SPLIT BETWEEN FAST AND SLOW WAVES,
! *CTUWUPDT* IS CALLED AGAIN FOR THE SLOW WAVES

IF (IFRELFMAX > 0 .AND. IFRELFMAX < NFRE_RED) THEN
  DELPRO = REAL(IDELPRO, JWRB)
  MSTART = IFRELFMAX + 1 
  MEND   = NFRE_RED

  CALL CTUWDRV (DELPRO, MSTART, MEND,                    &
 &              IJS, IJL, NINF, NSUP,                    &
 &              BLK2GLO,                                 &
 &              WLATM1, WCORM1, DP,                      &
 &              CGROUP_EXT, OMOSNH2KD_EXT,               &
 &              COSPHM1_EXT, DEPTH_EXT, U_EXT, V_EXT )

ENDIF


! FIND THE LOGICAL FLAGS THAT WILL LIMIT THE EXTEND OF THE CALCULATION IN PROPAGS2
! IN CASE REFRACTION IS USED

#ifndef _OPENACC
IF (IREFRA == 2 .OR. IREFRA == 3) THEN

  DO ICL=1,2
    DO IC=1,2
      DO K=1,NANG
        DO M=1,NFRE_RED
          LLWLATN(K,M,IC,ICL)=.FALSE.
          DO IJ=IJS,IJL
            IF (WLATN(IJ,K,M,IC,ICL) > 0.0_JWRB) THEN
              LLWLATN(K,M,IC,ICL)=.TRUE.
              EXIT
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  DO IC=1,2
    DO M=1,NFRE_RED
      DO K=1,NANG
        LLWLONN(K,M,IC)=.FALSE.
        DO IJ=IJS,IJL
          IF (WLONN(IJ,K,M,IC) > 0.0_JWRB) THEN
            LLWLONN(K,M,IC)=.TRUE.
            EXIT
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  DO ICL=1,2
    DO ICR=1,4
      DO M=1,NFRE_RED
        DO K=1,NANG
          LLWCORN(K,M,ICR,ICL)=.FALSE.
          DO IJ=IJS,IJL
            IF (WCORN(IJ,K,M,ICR,ICL) > 0.0_JWRB) THEN
              LLWCORN(K,M,ICR,ICL)=.TRUE.
              EXIT
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  DO IC=-1,1
    DO M=1,NFRE_RED
      DO K=1,NANG
        LLWKPMN(K,M,IC)=.FALSE.
        DO IJ=IJS,IJL
          IF (WKPMN(IJ,K,M,IC) > 0.0_JWRB) THEN
            LLWKPMN(K,M,IC)=.TRUE.
            EXIT
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  DO IC=-1,1
    DO M=1,NFRE_RED
      DO K=1,NANG
        LLWMPMN(K,M,IC)=.FALSE.
        DO IJ=IJS,IJL
          IF (WMPMN(IJ,K,M,IC) > 0.0_JWRB) THEN
            LLWMPMN(K,M,IC)=.TRUE.
            EXIT
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO

ENDIF
#endif

IF (ALLOCATED(THDD)) DEALLOCATE(THDD)
IF (ALLOCATED(THDC)) DEALLOCATE(THDC)
IF (ALLOCATED(SDOT)) DEALLOCATE(SDOT)

!$acc end data

IF (LHOOK) CALL DR_HOOK('CTUWUPDT',1,ZHOOK_HANDLE)

END SUBROUTINE CTUWUPDT
