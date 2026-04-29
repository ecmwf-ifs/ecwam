! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SWLDISSIP_ZBRY (KIJS, KIJL, FL1, FLD, SL,          &
     &                        WAVNUM, CGROUP,                       &
     &                        UFRIC, RAORW)
! ----------------------------------------------------------------------

!**** *SWLDISSIP_ZBRY* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.
!
!     JOSH KOUSAL & JEAN BIDLOT    ECMWF 2023
!
!*    PURPOSE.
!     --------
!     Turbulent dissipation of narrow-banded swell as described in
!     Babanin (2011, Section 7.5). 
!
!**   INTERFACE.
!     ----------

!       *CALL* *SWLDISSIP_ZBRY (KIJS, KIJL, FL1, FLD, SL,*
!                            WAVNUM, CGROUP,
!                            UFRIC, RAORW)*
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *FL1*    - SPECTRUM.
!          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*     - TOTAL SOURCE FUNCTION ARRAY
!          *WAVNUM* - WAVE NUMBER
!          *CGROUP* - GROUP SPEED
!          *UFRIC*  - FRICTION VELOCITY IN M/S.
!          *RAORW*  - RATIO AIR DENSITY TO WATER DENSITY


!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!     Babanin 2011: Cambridge Press, 295-321, 463pp.

!     ORIGIN.
!     ----------
!     Adapted from Babanin Young Donelan & Banner (ZBRY) physics 
!     as implemented as ST6 in WAVEWATCH-III
!     Implementation into ECWAM DECEMBER 2021 by J. Kousal 


! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR      , TH     ,ZPIFR   ,FRATIO    ,DELTH, DFIM,&
      &                     SIG     , DDEN
      USE YOWPCONS , ONLY : G        ,ZPI
      USE YOWPARAM , ONLY : NANG    ,NFRE
      USE YOWPHYS  , ONLY : LLSWL6CSTB1, ZSWL6B1

      USE YOMHOOK  , ONLY : LHOOK   ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(IN) :: FL1
REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: FLD, SL
REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: WAVNUM, CGROUP
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: UFRIC, RAORW 

INTEGER(KIND=JWIM) :: IJ, M, I, J, M2, K2, K, NANGD
INTEGER(KIND=JWIM), DIMENSION(KIJL) :: MPEAK

REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: ABAND, KMAX, ANAR, BN, DDIS
REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE) :: KK
REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE) :: A
REAL(KIND=JWRB), DIMENSION(KIJL)      :: B1, SUMDIR_IJ

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SWLDISSIP_ZBRY',0,ZHOOK_HANDLE)

DO M = 1, NFRE
  DO K = 1, NANG                    ! Apply to all directions
    DO IJ = KIJS,KIJL
      A(IJ,K,M) = FL1(IJ,K,M) * CGROUP(IJ,M) / ( ZPI * SIG(M) ) ! ACTION DENSITY SPECTRUM
    END DO
  END DO
END DO

!/ 0) --- Initialize parameters -------------------------------------- /
DO M = 1, NFRE
  DO IJ = KIJS,KIJL
    ABAND(IJ,M) = 0.0_JWRB
    DDIS(IJ,M)  = 0.0_JWRB
  END DO
END DO
DO M = 1, NFRE
  DO K = 1, NANG
    DO IJ = KIJS,KIJL
      ABAND(IJ,M) = ABAND(IJ,M) + A(IJ,K,M)
    END DO
  END DO
END DO

!/ 1) --- Choose calculation of steepness a*k ------------------------ /
!/        Replace the measure of steepness with the spectral
!         saturation after Banner et al. (2002) ---------------------- /
DO M = 1,NFRE
  DO IJ = KIJS,KIJL
    KMAX(IJ,M) = 0.0_JWRB
  END DO
  DO K = 1,NANG
    DO IJ = KIJS,KIJL
      KK(IJ,K,M) = A(IJ,K,M)
      KMAX(IJ,M) = MAX(KMAX(IJ,M), KK(IJ,K,M))
    END DO
  END DO
END DO

DO M = 1,NFRE
  DO K = 1,NANG
    DO IJ = KIJS,KIJL
      IF (KMAX(IJ,M).LT.1.0E-34_JWRB) THEN
        KK(IJ,K,M) = 1.0_JWRB
      ELSE
        KK(IJ,K,M) = KK(IJ,K,M)/KMAX(IJ,M)
      END IF
    END DO
  END DO
END DO

DO M = 1,NFRE
  DO IJ = KIJS,KIJL
    SUMDIR_IJ(IJ) = 0.0_JWRB
  END DO
  DO K = 1,NANG
    DO IJ = KIJS,KIJL
      SUMDIR_IJ(IJ) = SUMDIR_IJ(IJ) + KK(IJ,K,M)
    END DO
  END DO
  DO IJ = KIJS,KIJL
    ANAR(IJ,M) = 1.0_JWRB/( SUMDIR_IJ(IJ) * DELTH )
    BN(IJ,M) = ANAR(IJ,M) * ( ABAND(IJ,M) * SIG(M) * DELTH ) * WAVNUM(IJ,M)**3
  END DO
END DO

!
IF (.NOT.LLSWL6CSTB1) THEN
!/    --- A constant value for B1 attenuates swell too strong in the
!/        western central Pacific (i.e. cross swell less than 1.0m).
!/        Workaround is to scale B1 with steepness a*kp, where kp is
!/        the peak wavenumber. ZSWL6B1 remains a scaling constant, but
!/        with different magnitude.  --------------------------------- /
  DO IJ = KIJS,KIJL
    MPEAK(IJ) = MAXLOC(ABAND(IJ,:),1)         ! Index for peak
    SUMDIR_IJ(IJ) = 0.0_JWRB
  END DO
  DO K = 1,NFRE
    DO IJ = KIJS,KIJL
      SUMDIR_IJ(IJ) = SUMDIR_IJ(IJ) + ABAND(IJ,K)*DDEN(K)/CGROUP(IJ,K)
    END DO
  END DO
  DO IJ = KIJS,KIJL
    B1(IJ) = ZSWL6B1*(2.0_JWRB*SQRT(SUMDIR_IJ(IJ))*WAVNUM(IJ,MPEAK(IJ)))
  END DO
END IF
!
!/ 2) --- Calculate the derivative term only (in units of 1/s) ------- /
DO M = 1,NFRE
  DO IJ = KIJS,KIJL
    IF (ABAND(IJ,M) .GT. 1.0E-30_JWRB) THEN
      DDIS(IJ,M) = -(2.0_JWRB/3.0_JWRB) * B1(IJ) * SIG(M) * SQRT(BN(IJ,M))
    END IF
  END DO
END DO
!
!/ 3) --- Apply dissipation term of derivative to all directions ----- /
DO M = 1,NFRE
  DO K = 1, NANG
    DO IJ = KIJS,KIJL
      SL(IJ,K,M)  = SL(IJ,K,M)  + DDIS(IJ,M)*FL1(IJ,K,M)
      FLD(IJ,K,M) = FLD(IJ,K,M) + DDIS(IJ,M)
    END DO
  END DO
END DO


IF (LHOOK) CALL DR_HOOK('SWLDISSIP_ZBRY',1,ZHOOK_HANDLE)

END SUBROUTINE SWLDISSIP_ZBRY
