! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SDISSIP_ZBRY (KIJS, KIJL, FL1, FLD, SL,          &
     &                        WSWAVE, WAVNUM, CGROUP, XK2CG,             &
     &                        UFRIC, COSWDIF, RAORW)
! ----------------------------------------------------------------------

!**** *SDISSIP_ZBRY* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.

!     LOTFI AOUF       METEO FRANCE 2013
!     FABRICE ARDHUIN  IFREMER  2013


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

!       *CALL* *SDISSIP_ZBRY (KIJS, KIJL, FL1, FLD,SL,*
!                            WAVNUM, CGROUP, XK2CG,
!                            UFRIC, COSWDIF, RAORW)*
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *FL1*    - SPECTRUM.
!          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*     - TOTAL SOURCE FUNCTION ARRAY
!          *WAVNUM* - WAVE NUMBER
!          *CGROUP* - GROUP SPEED
!          *XK2CG*  - (WAVE NUMBER)**2 * GROUP SPEED
!          *UFRIC*  - FRICTION VELOCITY IN M/S.
!          *RAORW*  - RATIO AIR DENSITY TO WATER DENSITY
!          *COSWDIF*-  COS(TH(K)-WDWAVE(IJ))


!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!     IRANGE

!     REFERENCE.
!     ----------

!       Babanin et al. 2010: JPO 40(4), 667-683
!       Rogers  et al. 2012: JTECH 29(9) 1329-1346

!     ORIGIN.
!     ----------
!     Adapted from Babanin Young Donelan & Banner (ZBRY) physics 
!     as implemented as ST6 in WAVEWATCH-III
!     WW3 module:       W3SRC6MD
!     WW3 subroutine:   W3SDS6
!     Implementation into ECWAM DECEMBER 2021 by J. Kousal 


! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR      , TH     ,ZPIFR   ,FRATIO    ,DELTH, DFIM
      USE YOWPCONS , ONLY : G        ,ZPI
      USE YOWPARAM , ONLY : NANG    ,NFRE
      USE YOWPHYS  , ONLY : SDSBR   ,ISDSDTH ,ISB     ,IPSAT    ,      &
&                  SSDSC2  , SSDSC4, SSDSC6,  MICHE, SSDSC3, SSDSBRF1, &
&                  BRKPBCOEF ,SSDSC5, NSDSNTH,                         &
&                  INDICESSAT, SATWEIGHTS
      USE YOWSTAT , ONLY : IPHYS2_LOWWINDS

      USE YOMHOOK  , ONLY : LHOOK   ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "irange.intfb.h"      

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: FLD, SL
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: WAVNUM, CGROUP, XK2CG 
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: WSWAVE, UFRIC, RAORW 
      REAL(KIND=JWRB), DIMENSION(KIJL, NANG), INTENT(IN) :: COSWDIF 

      INTEGER(KIND=JWIM) :: IJ, K, M, I, J, M2, K2, NANGD
      INTEGER(KIND=JWIM) :: NSPEC !num. of freqs, dirs, spec. bins
      INTEGER(KIND=JWIM), DIMENSION(NANG) :: ITHN
      INTEGER(KIND=JWIM), DIMENSION(NFRE) :: IKN

      REAL(KIND=JWRB), PARAMETER :: SDS6A1  = 4.75E-6_JWRB  ! ST6 PARAM
      REAL(KIND=JWRB), PARAMETER :: SDS6A2  = 7.00E-5_JWRB  ! ST6 PARAM
      INTEGER(KIND=JWIM), PARAMETER :: SDS6P1  = 4     ! ST6 PARAM
      INTEGER(KIND=JWIM), PARAMETER :: SDS6P2  = 4     ! ST6 PARAM
      LOGICAL, PARAMETER   :: SDS6ET  = .TRUE.         ! ST6 PARAM

      REAL(KIND=JWRB), DIMENSION(NFRE) :: FREQ ! frequencies [Hz]
      REAL(KIND=JWRB), DIMENSION(NFRE) :: SIG  ! frequencies [RAD]
      REAL(KIND=JWRB), DIMENSION(NFRE) :: DFII ! frequency bandwiths [Hz]
      REAL(KIND=JWRB), DIMENSION(NFRE) :: ANAR ! directional narrowness
      REAL(KIND=JWRB), DIMENSION(NFRE) :: EDENS    ! spectral density E(f)
      REAL(KIND=JWRB), DIMENSION(NFRE) :: ETDENS ! threshold spec. density ET(f)
      REAL(KIND=JWRB), DIMENSION(NFRE) :: EXDENS ! excess spectral density EX(f)
      REAL(KIND=JWRB), DIMENSION(NFRE) :: NEXDENS! normalised excess spec.dens.
      REAL(KIND=JWRB), DIMENSION(NFRE) :: T1  ! inherent breaking term
      REAL(KIND=JWRB), DIMENSION(NFRE) :: T2  ! forced dissipation term
      REAL(KIND=JWRB), DIMENSION(NFRE) :: T12 ! =T1+T2 or combined dissipation
      REAL(KIND=JWRB), DIMENSION(NFRE) :: ADF ! temporary variable
      REAL(KIND=JWRB), DIMENSION(NFRE) :: DF  ! FREQUENCY INTERVALS
      REAL(KIND=JWRB) :: BNT ! empirical constant for wave breaking probability
      REAL(KIND=JWRB) :: XFAC, EDENSMAX ! temporary variableis

      REAL(KIND=JWRB), DIMENSION(NANG*NFRE) :: S, D, A
      REAL(KIND=JWRB), DIMENSION(NANG*NFRE) :: SIG2, CG2
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE) :: DDS

      
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: XK, CGG_WAM
      REAL(KIND=JWRB), DIMENSION(NFRE)      :: SIGP2


      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      
! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SDISSIP_ZBRY',0,ZHOOK_HANDLE)

      NSPEC = NANG * NFRE   ! NUMBER OF SPECTRAL BINS

      DO M = 1,NFRE
        SIG(M)   = ZPI*FR(M)
        SIGP2(M) = SIG(M)**2
      END DO

!     COMPUTE FREQUENCY INTERVALLS (borrowed from Wam_others/f4spec.F)
      DO M = 1,NFRE
        DF(M) = FR(M)*( FRATIO - 1.0_JWRB )
      ENDDO

      IKN    = IRANGE(1,NSPEC,NANG)   ! Index vector for elements of 1 ... NFRE
!                                    ! such that e.g. SIG(1:NFRE) = SIG2(IKN).
      DO K = 1, NANG                    ! Apply to all directions 
         SIG2   (IKN+(K-1)) = SIG
      END DO


      ! LOOP OVER LOCATIONS
      DO IJ = KIJS,KIJL

        DO K = 1, NANG                    ! Apply to all directions
           CG2    (IKN+(K-1)) = CGROUP(IJ,:)
        END DO

        A = RESHAPE( FL1(IJ,:,:) , (/NSPEC/)) * CG2 / ( ZPI * SIG2 )! ACTION DENSITY SPECTRUM
        ! WAM E(f,theta) to WW3 A(k,theta) conversion factor: CG2 / ( ZPI *SIG2 ) 

!/ 0) --- Initialize essential parameters ---------------------------- /
        IKN     = IRANGE(1,NSPEC,NANG)    ! Index vector for elements of 1,
!                                      ! 2,..., NFRE such that for example
!                                      ! SIG(1:NFRE) = SIG2(IKN).
        FREQ    = FR(1:NFRE)
        ANAR    = 1.0_JWRB
        BNT     = 0.035_JWRB**2
        T1      = 0.0_JWRB
        T2      = 0.0_JWRB
        NEXDENS = 0.0_JWRB
!
!/ 1) --- Calculate threshold spectral density, spectral density, and
!/        the level of exceedence EXDENS(f) -------------------------- /
!        ETDENS  = ( ZPI * BNT ) / ( ANAR * CGG(IJ,:) * WN(IJ,:)**3 )
        ETDENS  = ( ZPI * BNT ) / ( ANAR * CGROUP(IJ,:) * WAVNUM(IJ,:)**3 )
        !EDENS   = SUM(FL1(IJ,:,:),1) * ZPI * SIG * DELTH / CGG(IJ,:)  !E(f)
        EDENS   = SUM(FL1(IJ,:,:),1) * DELTH   !E(f)
        EXDENS  = MAX(0.0_JWRB,EDENS-ETDENS)
!
!/    --- normalise by a generic spectral density -------------------- /
        IF (SDS6ET) THEN                ! ww3_grid.inp: &SDS6 SDSET = T or F
           NEXDENS = EXDENS / ETDENS    ! normalise by threshold spectral density
        ELSE                            ! normalise by spectral density
           EDENSMAX = MAXVAL(EDENS)*1.0E-5_JWRB
           IF (ALL(EDENS .GT. EDENSMAX)) THEN
              NEXDENS = EXDENS / EDENS
           ELSE
              DO M = 1,NFRE
                 IF (EDENS(M) .GT. EDENSMAX) NEXDENS(M) = EXDENS(M) / EDENS(M)
              END DO
           END IF
        END IF
!
!/ 2) --- Calculate inherent breaking component T1 ------------------- /
        T1 = SDS6A1 * ANAR * FREQ * (NEXDENS**SDS6P1)
!
!/ 3) --- Calculate T2, the dissipation of waves induced by
!/        the breaking of longer waves T2 ---------------------------- /
        ADF    = ANAR * (NEXDENS**SDS6P2)
        XFAC   = (1.0_JWRB-1.0_JWRB/FRATIO)/(FRATIO-1.0_JWRB/FRATIO)
        DO M = 1,NFRE
           DFII(M) = DF(M) ! bug fix (spotted by Heinz): brought init into loc loop
!          IF (M .GT. 1) DFII(M) = DFII(M) * XFAC
           IF (M .GT. 1 .AND. M .LT. NFRE) DFII(M) = DFII(M) * XFAC
           T2(M) = SDS6A2 * SUM( ADF(1:M)*DFII(1:M) )
        END DO

!/ 4) --- Sum up dissipation terms and apply to all directions ------- /
        T12 = -1.0_JWRB * ( MAX(0.0_JWRB,T1)+MAX(0.0_JWRB,T2) )
        DO K = 1, NANG
           D(IKN+(K-1)) = T12
        END DO
!
        !S = D * A
!
!/ 5) --- Diagnostic output (switch !/T6) ---------------------------- /
!/T6     CALL STME21 ( TIME , IDTIME )
!/T6     WRITE (NDST,270) 'T1*E',IDTIME(1:19),(T1*EDENS)
!/T6     WRITE (NDST,270) 'T2*E',IDTIME(1:19),(T2*EDENS)
!/T6     WRITE (NDST,271) SUM(SUM(RESHAPE(S,(/ NANG,NFRE /)),1)*DDEN/CG)
!
!/T6     270 FORMAT (' TEST W3SDS6 : ',A,'(',A,')',':',70E11.3)
!/T6     271 FORMAT (' TEST W3SDS6 : Total SDS  =',E13.5)

        IF (.NOT. (IPHYS2_LOWWINDS .AND. WSWAVE(IJ)<=5._JWRB)) THEN
         ! no dissipation for U10<5m/s (following Muhammad Yasrab's work)
         ! i.e. don't update SL and FLD for low winds
         DDS = RESHAPE(D,(/NANG,NFRE/))
         DO M = 1,NFRE
            DO K = 1, NANG
               SL(IJ,K,M)  = SL(IJ,K,M)  + DDS(K,M)*FL1(IJ,K,M)
               FLD(IJ,K,M) = FLD(IJ,K,M) + DDS(K,M)
            END DO
         END DO
        END IF

      END DO
      ! END LOOP OVER LOC

      IF (LHOOK) CALL DR_HOOK('SDISSIP_ZBRY',1,ZHOOK_HANDLE)

      END SUBROUTINE SDISSIP_ZBRY
