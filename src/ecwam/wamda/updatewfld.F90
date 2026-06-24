SUBROUTINE UPDATEWFLD(KIJS, KIJL, MIJ,              &
 &                    HSOIB, HSAN, U10FG,           &
 &                    DEPTH, WSWAVE, WDWAVE, UFRIC, &
 &                    CICOVER, CITHICK, Z0M, Z0B, CHRNCK,    &
 &                    TAUW, TAUWDIR, AIRD, WSTAR, &
 &                    USTOKES, VSTOKES, STRNMS,     &
 &                    NEMOUSTOKES, NEMOVSTOKES, NEMOSTRN, &
 &                    WAVNUM, STOKFAC, CINV, XK2CG, &
 &                    FL1, XLLWS)
! ----------------------------------------------------------------------

!****  *UPDATEWFLD* - 

!     PURPOSE.
!     --------
!     UPDATES THE RELEVANT FIEDLS FOLLOWING ALITMETER DATA ASSIMILATION

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRO

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : EPSUS

      USE MPL_MODULE
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

#include "semean.intfb.h"
#include "stokestrn.intfb.h"
#include "update.intfb.h"
#include "updatewd.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL  ! GRID POINT INDEXES
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(OUT) :: MIJ  ! LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
      REAL(KIND=JWRB), DIMENSION(KIJS:kIJL), INTENT(IN) :: HSOIB  ! HS FROM OI ASSIMILATION
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: HSAN  ! ANALYSED HS AFTER SPECTRUM UPDATE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: U10FG ! FIRST GUESS U10
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: WSWAVE, WDWAVE, UFRIC, CICOVER, CITHICK
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: Z0M, Z0B, CHRNCK, TAUW, TAUWDIR, AIRD, WSTAR
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: USTOKES, VSTOKES, STRNMS
      REAL(KIND=JWRO), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: NEMOUSTOKES, NEMOVSTOKES, NEMOSTRN
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM, STOKFAC, CINV, XK2CG
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG, NFRE), INTENT(INOUT) :: FL1  ! SPECTRUM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG, NFRE), INTENT(OUT) :: XLLWS  ! WINDSEA MASK FROM INPUT SOURCE TERM


      INTEGER(KIND=JWIM) :: IJ

      REAL(KIND=JWRB) :: ZINV16
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: EMEAN, ETOIB, USA 

      LOGICAL :: LLEPSMIN

!     -------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('UPDATEWFLD',0,ZHOOK_HANDLE)


!*    3.1.1 TRANSFERS RESULTS OF OPTIMAL INT TO BLOCKS.
!           -------------------------------------------

      LLEPSMIN=.TRUE.
      ZINV16=1.0_JWRB/16.0_JWRB

      DO IJ=KIJS,KIJL
        USA(IJ) = MAX(UFRIC(IJ), EPSUS)
        U10FG(IJ) = WSWAVE(IJ)
      ENDDO

      DO IJ=KIJS,KIJL
        IF (HSOIB(IJ) > 0.0_JWRB) THEN
          ETOIB(IJ) = HSOIB(IJ)**2*ZINV16
        ELSE
          ETOIB(IJ) = -99.0_JWRB
        ENDIF
      ENDDO

!     UPDATE SPECTRA AND POTENTIALLY THE FRICTION VELOCITY
!     ----------------------------------------------------
      CALL UPDATE (KIJS, KIJL, FL1,       &
     &             WAVNUM, DEPTH,         &
     &             ETOIB, USA, WDWAVE)

!       THE WIND FIELD CURRENTLY USED IS UPDATED.
!       -----------------------------------------
      CALL UPDATEWD(KIJS, KIJL, FL1,      &
     &              WAVNUM, CINV, XK2CG,  &
     &              USA, WSWAVE, WDWAVE, UFRIC, &
     &              Z0M, Z0B, CHRNCK, TAUW, TAUWDIR, &
     &              AIRD, WSTAR, CICOVER, &
     &              MIJ, XLLWS)


!     Spectral update might not have taken all update to Hs from OI 
!     Compute the corresponding analysis Hs
      CALL SEMEAN (FL1, KIJS, KIJL, EMEAN, LLEPSMIN)
      DO IJ=KIJS,KIJL
        IF (EMEAN(IJ) > 0.0_JWRB) THEN
          HSAN(IJ)=4.0_JWRB*SQRT(EMEAN(IJ))
        ELSE
          HSAN(IJ)=0.0_JWRB
        ENDIF
      ENDDO

!     UPDATE SURFACE STOKES DRIFT and MEAN SQUARE STRAIN INTO THE SEA ICE

      CALL STOKESTRN (KIJS, KIJL, FL1, WAVNUM, STOKFAC, DEPTH, WSWAVE, WDWAVE, CICOVER, CITHICK, &
 &                    USTOKES, VSTOKES, STRNMS, NEMOUSTOKES, NEMOVSTOKES, NEMOSTRN)

IF (LHOOK) CALL DR_HOOK('UPDATEWFLD',1,ZHOOK_HANDLE)

END SUBROUTINE UPDATEWFLD
