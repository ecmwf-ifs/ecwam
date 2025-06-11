! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE FRCUTINDEX_DEFAULT (KIJS, KIJL, FM, FMWS, UFRIC, CICOVER,     &
     &                       MIJ)

! ----------------------------------------------------------------------

!**** *FRCUTINDEX_DEFAULT* - RETURNS THE LAST FREQUENCY INDEX OF
!                    PROGNOSTIC PART OF SPECTRUM.

!**   INTERFACE.
!     ----------

!       *CALL* *FRCUTINDEX_DEFAULT (KIJS, KIJL, FM, FMWS, CICOVER, MIJ)
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *FM*     - MEAN FREQUENCY
!          *FMWS*   - MEAN FREQUENCY OF WINDSEA
!          *UFRIC*  - FRICTION VELOCITY IN M/S
!          *CICOVER*- CICOVER 
!          *MIJ*    - LAST FREQUENCY INDEX for imposing high frequency tail


!     METHOD.
!     -------

!*    COMPUTES LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!*    FREQUENCIES LE 2.5*MAX(FMWS,FM).


!!! be aware that if this is NOT used, for iphys=1, the cumulative dissipation has to be
!!! re-activated (see module yowphys) !!!


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,DFIM       ,FRATIO   ,FLOGSPRDM1, &
     &                ZPIFR,                                            &
     &                DELTH          ,RHOWG_DFIM ,FRIC
      USE YOWICE   , ONLY : CITHRSH_TAIL
      USE YOWPARAM , ONLY : NFRE
      USE YOWPCONS , ONLY : G        ,EPSMIN
      USE YOWPHYS  , ONLY : TAILFACTOR, TAILFACTOR_PM

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM), INTENT(OUT) :: MIJ(KIJL)
      REAL(KIND=JWRB),DIMENSION(KIJL), INTENT(IN) :: FM, FMWS, UFRIC, CICOVER


      INTEGER(KIND=JWIM) :: IJ, M

      REAL(KIND=JWRB) :: FPMH, FPPM, FM2, FPM, FPM4
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('FRCUTINDEX_DEFAULT',0,ZHOOK_HANDLE)

!*    COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!*    FREQUENCIES LE MAX(TAILFACTOR*MAX(FMNWS,FM),TAILFACTOR_PM*FPM),
!*    WHERE FPM IS THE PIERSON-MOSKOWITZ FREQUENCY BASED ON FRICTION
!*    VELOCITY. (FPM=G/(FRIC*ZPI*USTAR))
!     ------------------------------------------------------------

      FPMH = TAILFACTOR/FR(1)
      FPPM = TAILFACTOR_PM*G/(FRIC*ZPIFR(1))

      DO IJ=KIJS,KIJL
        IF (CICOVER(IJ) <= CITHRSH_TAIL) THEN
          FM2 = MAX(FMWS(IJ),FM(IJ))*FPMH
          FPM = FPPM/MAX(UFRIC(IJ),EPSMIN)
          FPM4 = MAX(FM2,FPM)
          MIJ(IJ) = NINT(LOG10(FPM4)*FLOGSPRDM1)+1
          MIJ(IJ) = MIN(MAX(1,MIJ(IJ)),NFRE)
        ELSE
          MIJ(IJ) = NFRE
        ENDIF
      ENDDO

      IF (LHOOK) CALL DR_HOOK('FRCUTINDEX_DEFAULT',1,ZHOOK_HANDLE)

      END SUBROUTINE FRCUTINDEX_DEFAULT
