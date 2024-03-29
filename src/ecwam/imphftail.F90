! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE IMPHFTAIL (KIJS, KIJL, MIJ, FLM, WAVNUM, XK2CG, FL1) 
! ----------------------------------------------------------------------

!**** *IMPHFTAIL* - IMPOSE A HIGH FREQUENCY TAIL TO THE SPECTRUM


!*    PURPOSE.
!     --------

!     IMPOSE A HIGH FREQUENCY TAIL TO THE SPECTRUM ABOVE FREQUENCY INDEX MIJ


!**   INTERFACE.
!     ----------

!       *CALL* *IMPHFTAIL (KIJS, KIJL, MIJ, FLM, WAVNUM, XK2CG, FL1)
!          *KIJS*    - INDEX OF FIRST GRIDPOINT
!          *KIJL*    - INDEX OF LAST GRIDPOINT
!          *MIJ*     - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!          *FLM*     - SPECTAL DENSITY MINIMUM VALUE
!          *WAVNUM*  - WAVENUMBER
!          *XK2CG*   - (WAVNUM)**2 * GROUP SPEED
!          *FL1*     - SPECTRUM (INPUT AND OUTPUT).

!     METHOD.
!     -------

!     EXTERNALS.
!     ---------

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE YOWPARAM , ONLY : NANG     ,NFRE

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM), DIMENSION(KIJL), INTENT(IN) :: MIJ
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG), INTENT(IN) :: FLM 
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: WAVNUM, XK2CG
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: FL1


      INTEGER(KIND=JWIM) :: IJ, K, M

      REAL(KIND=JWRB) :: AKM1, TFAC
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJL) :: TEMP1, TEMP2

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('IMPHFTAIL',0,ZHOOK_HANDLE)

!*    DIAGNOSTIC TAIL.
!     ----------------

      DO IJ=KIJS,KIJL
        TEMP1(IJ) = 1.0_JWRB/XK2CG(IJ,MIJ(IJ))/WAVNUM(IJ,MIJ(IJ))

        DO M=MIJ(IJ)+1,NFRE
          TEMP2(IJ) = 1.0_JWRB/XK2CG(IJ,M)/WAVNUM(IJ,M)
          TEMP2(IJ) = TEMP2(IJ)/TEMP1(IJ)

!*    MERGE TAIL INTO SPECTRA.
!     ------------------------
          DO K=1,NANG
            TFAC = FL1(IJ,K,MIJ(IJ))
            FL1(IJ,K,M) = MAX(TEMP2(IJ)*TFAC,FLM(IJ,K))
          ENDDO
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('IMPHFTAIL',1,ZHOOK_HANDLE)

      END SUBROUTINE IMPHFTAIL
