! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SDISSIP_ZBRY (KIJS, KIJL, FL1, FLD, SL,          &
     &                        WSWAVE, WAVNUM, CGROUP,             &
     &                        UFRIC, RAORW)
! ----------------------------------------------------------------------

!**** *SDISSIP_ZBRY* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.
!
!     JOSH KOUSAL & JEAN BIDLOT    ECMWF 2023
!
!*    PURPOSE.
!     --------
!      Observation-based source term for dissipation after Babanin et al.
!      (2010) following the implementation by Rogers et al. (2012). The
!      dissipation function Sds accommodates an inherent breaking term T1
!      and an additional cumulative term T2 at all frequencies above the
!      peak. The forced dissipation term T2 is an integral that grows
!      toward higher frequencies and dominates at smaller scales
!      (Babanin et al. 2010).
!
!**   INTERFACE.
!     ----------

!       *CALL* *SDISSIP_ZBRY (KIJS, KIJL, FL1, FLD, SL,*
!                            WSWAVE, WAVNUM, CGROUP,
!                            UFRIC, RAORW)*
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *FL1*    - SPECTRUM.
!          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*     - TOTAL SOURCE FUNCTION ARRAY
!          *WSWAVE* - WIND SPEED IN M/S.
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

!       Babanin et al. 2010: JPO 40(4), 667-683
!       Rogers  et al. 2012: JTECH 29(9) 1329-1346

!     ORIGIN.
!     ----------
!     Adapted from Babanin Young Donelan & Banner (ZBRY) physics 
!     as implemented as ST6 in WAVEWATCH-III
!     Implementation into ECWAM DECEMBER 2021 by J. Kousal 


! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR      , TH     ,ZPIFR   ,FRATIO    ,DELTH, DFIM,&
      &                     SIG     , DF
      USE YOWPCONS , ONLY : G        ,ZPI
      USE YOWPARAM , ONLY : NANG    ,NFRE
      USE YOWSTAT  , ONLY : LLLOWWINDS
      USE YOWPHYS  , ONLY : LLSDS6ET, ISDS6P1, ISDS6P2, ZSDS6A1, ZSDS6A2

      USE YOMHOOK  , ONLY : LHOOK   ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(IN) :: FL1
REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: FLD, SL
REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: WAVNUM, CGROUP
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: WSWAVE, UFRIC, RAORW 

INTEGER(KIND=JWIM) :: IJ, K, M, I, J

REAL(KIND=JWRB), DIMENSION(NFRE) :: FREQ ! frequencies [Hz]
REAL(KIND=JWRB), DIMENSION(NFRE) :: DFII ! frequency bandwiths [Hz]
REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: ANAR ! directional narrowness
REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: EDENS    ! spectral density E(f)
REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: ETDENS ! threshold spec. density ET(f)
REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: EXDENS ! excess spectral density EX(f)
REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: NEXDENS! normalised excess spec.dens.
REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: T1  ! inherent breaking term
REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: T2  ! forced dissipation term
REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: T12 ! =T1+T2 or combined dissipation
REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: ADF ! temporary variable
REAL(KIND=JWRB) :: BNT ! empirical constant for wave breaking probability
REAL(KIND=JWRB) :: XFAC ! temporary variableis
REAL(KIND=JWRB), DIMENSION(KIJL) :: EDENSMAX ! temporary variable
REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE) :: A
REAL(KIND=JWRB), DIMENSION(KIJL)           :: CUMADF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SDISSIP_ZBRY',0,ZHOOK_HANDLE)

DO M = 1, NFRE
   DO K = 1, NANG                    ! Apply to all directions
      DO IJ = KIJS,KIJL
         A(IJ,K,M) = FL1(IJ,K,M) * CGROUP(IJ,M) / ( ZPI * SIG(M) ) ! ACTION DENSITY SPECTRUM
      END DO
   END DO
END DO

!/ 0) --- Initialize essential parameters ---------------------------- /
FREQ(:) = FR(1:NFRE)
BNT     = 0.035_JWRB**2
DO M = 1, NFRE
   DO IJ = KIJS,KIJL
      ANAR(IJ,M)    = 1.0_JWRB
      T1(IJ,M)      = 0.0_JWRB
      T2(IJ,M)      = 0.0_JWRB
      NEXDENS(IJ,M) = 0.0_JWRB
   END DO
END DO
!
!/ 1) --- Calculate threshold spectral density, spectral density, and
!/        the level of exceedence EXDENS(f) -------------------------- /
   DO M = 1, NFRE
      DO IJ = KIJS,KIJL
         ETDENS(IJ,M) = ( ZPI * BNT ) / ( ANAR(IJ,M) * CGROUP(IJ,M) * WAVNUM(IJ,M)**3 )
         EDENS(IJ,M) = 0.0_JWRB
      END DO
      DO K = 1, NANG
         DO IJ = KIJS,KIJL
            EDENS(IJ,M) = EDENS(IJ,M) + FL1(IJ,K,M)
         END DO
      END DO
      DO IJ = KIJS,KIJL
         EDENS(IJ,M) = EDENS(IJ,M) * DELTH   ! E(f)
         EXDENS(IJ,M) = MAX(0.0_JWRB,EDENS(IJ,M)-ETDENS(IJ,M))
      END DO
   END DO
!
!/    --- normalise by a generic spectral density -------------------- /
IF (LLSDS6ET) THEN
   DO M = 1,NFRE
      DO IJ = KIJS,KIJL
         NEXDENS(IJ,M) = EXDENS(IJ,M) / ETDENS(IJ,M)    ! normalise by threshold spectral density
      END DO
   END DO
ELSE                            ! normalise by spectral density
   DO IJ = KIJS,KIJL
      EDENSMAX(IJ) = MAXVAL(EDENS(IJ,1:NFRE))*1.0E-5_JWRB
   END DO
   DO M = 1,NFRE
      DO IJ = KIJS,KIJL
         IF (EDENS(IJ,M) .GT. EDENSMAX(IJ)) THEN
            NEXDENS(IJ,M) = EXDENS(IJ,M) / EDENS(IJ,M)
         END IF
      END DO
   END DO
END IF
!
!/ 2) --- Calculate inherent breaking component T1 ------------------- /
   DO M = 1,NFRE
      DO IJ = KIJS,KIJL
         T1(IJ,M) = ZSDS6A1 * ANAR(IJ,M) * FREQ(M) * (NEXDENS(IJ,M)**ISDS6P1)
      END DO
   END DO
!
!/ 3) --- Calculate T2, the dissipation of waves induced by
!/        the breaking of longer waves T2 ---------------------------- /
   DO M = 1,NFRE
      DO IJ = KIJS,KIJL
         ADF(IJ,M) = ANAR(IJ,M) * (NEXDENS(IJ,M)**ISDS6P2)
      END DO
   END DO

XFAC   = (1.0_JWRB-1.0_JWRB/FRATIO)/(FRATIO-1.0_JWRB/FRATIO)
DO M = 1,NFRE
   DFII(M) = DF(M)
   IF (M .GT. 1 .AND. M .LT. NFRE) THEN
      DFII(M) = DFII(M) * XFAC
   END IF
END DO

DO IJ = KIJS,KIJL
   CUMADF(IJ) = 0.0_JWRB
END DO
DO M = 1,NFRE
   DO IJ = KIJS,KIJL
      CUMADF(IJ) = CUMADF(IJ) + ADF(IJ,M)*DFII(M)
      T2(IJ,M) = ZSDS6A2 * CUMADF(IJ)
   END DO
END DO

!/ 4) --- Sum up dissipation terms and apply to all directions ------- /
DO M = 1,NFRE
   DO IJ = KIJS,KIJL
      T12(IJ,M) = -1.0_JWRB * ( MAX(0.0_JWRB,T1(IJ,M))+MAX(0.0_JWRB,T2(IJ,M)) )
   END DO
END DO

IF (LLLOWWINDS) THEN
   DO M = 1,NFRE
      DO K = 1, NANG
         DO IJ = KIJS,KIJL
            IF ( WSWAVE(IJ)>=5._JWRB) THEN
               ! no dissipation for winds<5m/s (following Muhammad Yasrab's work)
               SL(IJ,K,M)  = SL(IJ,K,M)  + T12(IJ,M)*FL1(IJ,K,M)
               FLD(IJ,K,M) = FLD(IJ,K,M) + T12(IJ,M)
            END IF
         END DO
      END DO
   END DO
ELSE
   DO M = 1,NFRE
      DO K = 1, NANG
         DO IJ = KIJS,KIJL
            SL(IJ,K,M)  = SL(IJ,K,M)  + T12(IJ,M)*FL1(IJ,K,M)
            FLD(IJ,K,M) = FLD(IJ,K,M) + T12(IJ,M)
         END DO
      END DO
   END DO
END IF

IF (LHOOK) CALL DR_HOOK('SDISSIP_ZBRY',1,ZHOOK_HANDLE)

END SUBROUTINE SDISSIP_ZBRY
