! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE WAMINTGR (CDTPRA, CDATE, CDATEWH, CDTIMP, CDTIMPNEXT,  &
 &                   BLK2GLO,                                     &
 &                   WVENVI, WVPRPT, FF_NOW, FF_NEXT, INTFLDS,    &
 &                   WAM2NEMO, MIJ, FL1, XLLWS)


! ----------------------------------------------------------------------

!**** *WAMINTGR* - 3-G WAM MODEL - TIME INTEGRATION OF WAVE FIELDS.

!*    PURPOSE.
!     --------

!       COMPUTATION OF THE 2-D FREQUENCY-DIRECTION WAVE SPECTRUM AT ALL
!       GRID POINTS FOR A GIVEN INITIAL SPECTRUM AND FORCING SURFACE
!       STRESS FIELD.

!     REFERENCE.
!     ----------

!         IFS DOCUMENTATION, part VII 

! -------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
USE YOWDRVTYPE  , ONLY : WVGRIDGLO, ENVIRONMENT, FREQUENCY, FORCING_FIELDS,  &
 &                       INTGT_PARAM_FIELDS, WAVE2OCEAN

USE YOWCOUP  , ONLY : LWNEMOCOU, NEMONTAU
USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
USE YOWPARAM , ONLY : NIBLO, NANG, NFRE
USE YOWPCONS , ONLY : EPSMIN
USE YOWSTAT  , ONLY : CDTPRO, IDELPRO, IDELT, IDELWI, LLSOURCE
USE YOWWIND  , ONLY : CDAWIFL, CDATEWO, CDATEFL

USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
      
! ----------------------------------------------------------------------

IMPLICIT NONE

#include "implsch.intfb.h"
#include "incdate.intfb.h"
#include "newwind.intfb.h"
#include "propag_wam.intfb.h"

CHARACTER(LEN=14), INTENT(IN)                                            :: CDTPRA  ! DATE FOR CALL PROPAGATION
CHARACTER(LEN=14), INTENT(INOUT)                                         :: CDATE  ! CURRENT DATE
CHARACTER(LEN=14), INTENT(INOUT)                                         :: CDATEWH ! DATE OF THE NEXT FORCING FIELDS
CHARACTER(LEN=14), INTENT(INOUT)                                         :: CDTIMP  ! START DATE OF SOURCE FUNCTION INTEGRATION
CHARACTER(LEN=14), INTENT(INOUT)                                         :: CDTIMPNEXT  ! NEXT START DATE OF SOURCE FUNCTION INTEGRATION
TYPE(WVGRIDGLO), INTENT(IN)                                              :: BLK2GLO  ! BLOCK TO GRID TRANSFORMATION
TYPE(ENVIRONMENT), INTENT(INOUT)                                         :: WVENVI !  WAVE ENVIRONMENT FIELDS
TYPE(FREQUENCY), INTENT(INOUT)                                           :: WVPRPT  ! WAVE PROPERTIES FIELDS
TYPE(FORCING_FIELDS), INTENT(INOUT)                                      :: FF_NOW  ! FORCING FIELDS AT CURRENT TIME
TYPE(FORCING_FIELDS), INTENT(IN)                                         :: FF_NEXT  !  DATA STRUCTURE WITH THE NEXT FORCING FIELDS
TYPE(INTGT_PARAM_FIELDS), INTENT(INOUT)                                  :: INTFLDS  ! INTEGRATED/DERIVED PARAMETERS
TYPE(WAVE2OCEAN), INTENT(INOUT)                                          :: WAM2NEMO  ! WAVE FIELDS PASSED TO NEMO
INTEGER(KIND=JWIM), DIMENSION(NPROMA_WAM, NCHNK), INTENT(INOUT)          :: MIJ  ! LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE
REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(INOUT) :: FL1
REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(INOUT) :: XLLWS ! TOTAL WINDSEA MASK FROM INPUT SOURCE TERM


INTEGER(KIND=JWIM) :: IJ, K, M
INTEGER(KIND=JWIM) :: ICHNK
INTEGER(KIND=JWIM) :: IDELWH

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


LOGICAL, SAVE :: LLNEWFILE

DATA LLNEWFILE / .FALSE. /

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WAMINTGR',0,ZHOOK_HANDLE)

!*     PROPAGATION TIME
!      ----------------

IF (CDATE == CDTPRA) THEN
  CALL PROPAG_WAM(BLK2GLO, WVENVI, WVPRPT, FL1)
  CDATE = CDTPRO
ENDIF


!* RETRIEVING NEW FORCING FIELDS IF NEEDED.
!  ----------------------------------------
CALL NEWWIND(CDTIMP, CDATEWH, LLNEWFILE,      &
 &           WVPRPT, FF_NOW, FF_NEXT)

! IT IS TIME TO INTEGRATE THE SOURCE TERMS
! ----------------------------------------
IF (CDATE >= CDTIMPNEXT) THEN
! COMPUTE UPDATE DUE TO SOURCE TERMS
  CALL GSTATS(1431,0)
  IF (LLSOURCE) THEN
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ICHNK)
    DO ICHNK = 1, NCHNK
      CALL IMPLSCH (1, NPROMA_WAM, FL1(:,:,:,ICHNK),   &
 &                  WVPRPT%WAVNUM(:,:,ICHNK), WVPRPT%CGROUP(:,:,ICHNK), WVPRPT%CIWA(:,:,ICHNK), &
 &                  WVPRPT%CINV(:,:,ICHNK), WVPRPT%XK2CG(:,:,ICHNK), WVPRPT%STOKFAC(:,:,ICHNK), &
 &                  WVENVI%EMAXDPT(:,ICHNK), WVENVI%INDEP(:,ICHNK), &
 &                  WVENVI%DEPTH(:,ICHNK), WVENVI%IOBND(:,ICHNK), WVENVI%IODP(:,ICHNK), &
 &                  FF_NOW%AIRD(:,ICHNK), FF_NOW%WDWAVE(:,ICHNK), FF_NOW%CICOVER(:,ICHNK), &
 &                  FF_NOW%WSWAVE(:,ICHNK), FF_NOW%WSTAR(:,ICHNK), &
 &                  FF_NOW%UFRIC(:,ICHNK), FF_NOW%TAUW(:,ICHNK), FF_NOW%TAUWDIR(:,ICHNK), &
 &                  FF_NOW%Z0M(:,ICHNK), FF_NOW%Z0B(:,ICHNK), FF_NOW%CHRNCK(:,ICHNK), &
 &                  FF_NOW%CITHICK(:,ICHNK), &
 &                  WAM2NEMO%NEMOUSTOKES(:,ICHNK), WAM2NEMO%NEMOVSTOKES(:,ICHNK), WAM2NEMO%NEMOSTRN(:,ICHNK), &
 &                  WAM2NEMO%NPHIEPS(:,ICHNK), WAM2NEMO%NTAUOC(:,ICHNK), WAM2NEMO%NSWH(:,ICHNK), &
 &                  WAM2NEMO%NMWP(:,ICHNK), WAM2NEMO%NEMOTAUX(:,ICHNK), WAM2NEMO%NEMOTAUY(:,ICHNK), &
 &                  WAM2NEMO%NEMOWSWAVE(:,ICHNK), WAM2NEMO%NEMOPHIF(:,ICHNK), &
 &                  INTFLDS%WSEMEAN(:,ICHNK), INTFLDS%WSFMEAN(:,ICHNK), &
 &                  INTFLDS%USTOKES(:,ICHNK), INTFLDS%VSTOKES(:,ICHNK), INTFLDS%STRNMS(:,ICHNK), &
 &                  INTFLDS%TAUXD(:,ICHNK), INTFLDS%TAUYD(:,ICHNK), INTFLDS%TAUOCXD(:,ICHNK), &
 &                  INTFLDS%TAUOCYD(:,ICHNK), INTFLDS%TAUOC(:,ICHNK), INTFLDS%PHIOCD(:,ICHNK), &
 &                  INTFLDS%PHIEPS(:,ICHNK), INTFLDS%PHIAW(:,ICHNK), &
 &                  MIJ(:,ICHNK), XLLWS(:,:,:,ICHNK) )

    ENDDO
!$OMP END PARALLEL DO

    IF (LWNEMOCOU) NEMONTAU = NEMONTAU + 1

  ELSE
!   NO SOURCE TERM CONTRIBUTION
!$OMP      PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK)  
    DO ICHNK = 1, NCHNK
      MIJ(:,ICHNK) = NFRE
      FL1(:,:,:,ICHNK) = MAX(FL1(:,:,:,ICHNK), EPSMIN)
      XLLWS(:,:,:,ICHNK) = 0.0_JWRB
    ENDDO
!$OMP      END PARALLEL DO
  ENDIF
  CALL GSTATS(1431,1)


!*       UPDATE FORCING FIELDS TIME COUNTER 
!        ----------------------------------
  IF (LLNEWFILE) THEN
    LLNEWFILE = .FALSE.
    IDELWH = MAX(IDELWI, IDELPRO)
    CALL INCDATE(CDAWIFL, IDELWH)
    CALL INCDATE(CDATEFL, IDELWH)
  ENDIF

  CDATEWO=CDATEWH
  CDTIMP=CDTIMPNEXT
  CALL INCDATE(CDTIMPNEXT, IDELT)   

ENDIF

IF (LHOOK) CALL DR_HOOK('WAMINTGR',1,ZHOOK_HANDLE)

END SUBROUTINE WAMINTGR
