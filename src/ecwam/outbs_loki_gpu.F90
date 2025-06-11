! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE OUTBS_LOKI_GPU (MIJ, FL1, XLLWS,                            &
 &                WVPRPT, WVENVI, FF_NOW, INTFLDS, NEMO2WAM,  &
 &                BOUT)
! ----------------------------------------------------------------------

!**** *OUTBS* - MODEL OUTPUT FROM BLOCK TO FILE, PRINTER AND COMMON.

!*    PURPOSE.
!     --------

!       CONTROL OUTPUT OF WAVE AND WIND FIELDS (except spectrum).

!**   INTERFACE.
!     ----------
!      *CALL*OUTBS (MIJ, FL1, XLLWS,
!                   WVPRPT, WVENVI, FF_NOW, INTFLDS, NEMO2WAM, BOUT)
!      *MIJ*     - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!      *FL1*     - INPUT SPECTRUM.
!      *XLLWS*   - WINDSEA MASK FROM INPUT SOURCE TERM
!      *WVENVI*  - WAVE ENVIRONMENT (depth, currents,...)
!      *FF_NOW*  - FORCING FIELDS
!      *INTFLDS* - INTEGRATED/DERIVED PARAMETERS
!      *NEMO2WAM*- FIELDS FRON OCEAN MODEL to WAM
!      *BOUT*    - OUTPUT PARAMETERS BUFFER 



!     EXTERNALS.
!     ----------

!       *OUTBLOCK*  - GET ALL OUTPUT PARAMETERS
!   
!     METHOD.
!     -------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : ENVIRONMENT, FREQUENCY, FORCING_FIELDS,  &
     &                         INTGT_PARAM_FIELDS, OCEAN2WAVE

      USE YOWCOUT  , ONLY : JPPFLAG  ,NIPRMOUT
      USE YOWCOUP  , ONLY : LLNORMWAMOUT
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSTAT  , ONLY : LUPDATE_GPU_GLOBALS_OUTBS

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

#include "outblock.intfb.h"

      INTEGER(KIND=JWIM), DIMENSION(NPROMA_WAM, NCHNK), INTENT(IN)          :: MIJ
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(IN) :: XLLWS
      TYPE(FREQUENCY), INTENT(IN)                                           :: WVPRPT
      TYPE(ENVIRONMENT), INTENT(IN)                                         :: WVENVI
      TYPE(FORCING_FIELDS), INTENT(IN)                                      :: FF_NOW
      TYPE(INTGT_PARAM_FIELDS), INTENT(IN)                                  :: INTFLDS
      TYPE(OCEAN2WAVE), INTENT(IN)                                          :: NEMO2WAM
      REAL(KIND=JWRB), POINTER, CONTIGUOUS, INTENT(INOUT)                   :: BOUT(:,:,:)


      INTEGER(KIND=JWIM) :: ICHNK

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE, ZHOOK_HANDLE_DATA_OFFLOAD

      LOGICAL :: LDREPROD

! ----------------------------------------------------------------------


!*    1. COMPUTE MEAN PARAMETERS.
!        ------------------------

!     COMPUTE MEAN PARAMETERS
IF (LUPDATE_GPU_GLOBALS_OUTBS) THEN
IF (LHOOK) CALL DR_HOOK('DATA_OFFLOAD',0,ZHOOK_HANDLE_DATA_OFFLOAD)
!$loki update_device
IF (LHOOK) CALL DR_HOOK('DATA_OFFLOAD',1,ZHOOK_HANDLE_DATA_OFFLOAD)
ENDIF

IF (LHOOK) CALL DR_HOOK('OUTBS',0,ZHOOK_HANDLE)
      CALL GSTATS(1502,0)
!$loki structured-data present(MIJ,WVPRPT,WVENVI,INTFLDS,FF_NOW,NEMO2WAM,BOUT)

      DO ICHNK = 1, NCHNK
        CALL OUTBLOCK(1, NPROMA_WAM, MIJ(:,ICHNK),                        &
     &                FL1(:,:,:,ICHNK), XLLWS(:,:,:,ICHNK),               &
     &                WVPRPT%WAVNUM(:,:,ICHNK), WVPRPT%CINV(:,:,ICHNK), WVPRPT%CGROUP(:,:,ICHNK), &
     &                WVENVI%DEPTH(:,ICHNK), WVENVI%UCUR(:,ICHNK), WVENVI%VCUR(:,ICHNK), &
     &                WVENVI%IODP(:,ICHNK),WVENVI%IBRMEM(:,ICHNK),                               &
     &                INTFLDS%ALTWH(:,ICHNK), INTFLDS%CALTWH(:,ICHNK), INTFLDS%RALTCOR(:,ICHNK), &
     &                INTFLDS%USTOKES(:,ICHNK), INTFLDS%VSTOKES(:,ICHNK), INTFLDS%STRNMS(:,ICHNK), &
     &                INTFLDS%TAUXD(:,ICHNK), INTFLDS%TAUYD(:,ICHNK), INTFLDS%TAUOCXD(:,ICHNK), &
     &                INTFLDS%TAUOCYD(:,ICHNK), INTFLDS%TAUOC(:,ICHNK), &
     &                INTFLDS%TAUICX(:,ICHNK), INTFLDS%TAUICY(:,ICHNK), INTFLDS%PHIOCD(:,ICHNK), &
     &                INTFLDS%PHIEPS(:,ICHNK), INTFLDS%PHIAW(:,ICHNK), &
     &                FF_NOW%AIRD(:,ICHNK), FF_NOW%WDWAVE(:,ICHNK), FF_NOW%CICOVER(:,ICHNK), &
     &                FF_NOW%WSWAVE(:,ICHNK), FF_NOW%WSTAR(:,ICHNK), &
     &                FF_NOW%UFRIC(:,ICHNK), FF_NOW%TAUW(:,ICHNK), &
     &                FF_NOW%Z0M(:,ICHNK), FF_NOW%Z0B(:,ICHNK), FF_NOW%CHRNCK(:,ICHNK), &
     &                FF_NOW%CITHICK(:,ICHNK), &
     &                NEMO2WAM%NEMOCICOVER(:,ICHNK), &
     &                NEMO2WAM%NEMOCITHICK(:, ICHNK), NEMO2WAM%NEMOUCUR(:,ICHNK), &
     &                NEMO2WAM%NEMOVCUR(:, ICHNK), &
     &                BOUT(:,:,ICHNK))
      ENDDO

!$loki end structured-data
      CALL GSTATS(1502,1)

      LUPDATE_GPU_GLOBALS_OUTBS = .FALSE.


IF (LHOOK) CALL DR_HOOK('OUTBS',1,ZHOOK_HANDLE)

END SUBROUTINE OUTBS_LOKI_GPU
