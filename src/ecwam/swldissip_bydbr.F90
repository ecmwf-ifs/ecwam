! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SWLDISSIP_BYDBR (KIJS, KIJL, FL1, FLD, SL,          &
     &                        WAVNUM, CGROUP, XK2CG,             &
     &                        UFRIC, COSWDIF, RAORW)
! ----------------------------------------------------------------------

!**** *SWLDISSIP_BYDBR* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.

!     LOTFI AOUF       METEO FRANCE 2013
!     FABRICE ARDHUIN  IFREMER  2013


!*    PURPOSE.
!     --------
!     Turbulent dissipation of narrow-banded swell as described in
!     Babanin (2011, Section 7.5). 
!
!**   INTERFACE.
!     ----------

!       *CALL* *SWLDISSIP_BYDBR (KIJS, KIJL, FL1, FLD,SL,*
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

!     Babanin 2011: Cambridge Press, 295-321, 463pp.

!     ORIGIN.
!     ----------
!     Adapted from Babanin Young Donelan & Banner (BYDB) physics 
!     as implemented as ST6 in WAVEWATCH-III
!     WW3 module:       W3SWLDMD
!     WW3 subroutine:   W3SWL6
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

      USE YOMHOOK  , ONLY : LHOOK   ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "irange.intfb.h"      

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: FLD, SL
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: WAVNUM, CGROUP, XK2CG 
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: UFRIC, RAORW 
      REAL(KIND=JWRB), DIMENSION(KIJL, NANG), INTENT(IN) :: COSWDIF 

      INTEGER(KIND=JWIM) :: IJ, M, I, J, M2, K2, K, NANGD
      INTEGER(KIND=JWIM) :: NSPEC !num. of freqs, dirs, spec. bins
      INTEGER(KIND=JWIM), DIMENSION(NANG) :: KKD
      INTEGER(KIND=JWIM), DIMENSION(NANG) :: ITHN
      INTEGER(KIND=JWIM), DIMENSION(NFRE) :: IKN

      REAL(KIND=JWRB), PARAMETER :: SWL6B1    = 0.0041_JWRB  ! ST6 PARAM
      LOGICAL, PARAMETER         :: SWL6CSTB1 = .FALSE.      ! ST6 PARAM

      REAL(KIND=JWRB), DIMENSION(NFRE) :: ABAND, KMAX, ANAR, BN, AORB, DDIS
      REAL(KIND=JWRB), DIMENSION(NFRE) :: SIG, DDEN
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE) :: KK
      REAL(KIND=JWRB), DIMENSION(NANG*NFRE) :: S, D, A
      REAL(KIND=JWRB), DIMENSION(NANG*NFRE) :: SIG2, CG2
      REAL(KIND=JWRB)                       :: B1
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE) :: DSWL

      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: XK, CGG_WAM
      REAL(KIND=JWRB), DIMENSION(NFRE)      :: SIGP2
      REAL(KIND=JWRB), DIMENSION(NFRE) :: DF  ! FREQUENCY INTERVALS

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      
! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SWLDISSIP_BYDBR',0,ZHOOK_HANDLE)

      NSPEC = NANG * NFRE   ! NUMBER OF SPECTRAL BINS

      DO M = 1,NFRE
        SIG(M)   = ZPI*FR(M)
        SIGP2(M) = SIG(M)**2
        DDEN(M) = ZPI*DFIM(M)*SIG(M)
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

!/ 0) --- Initialize parameters -------------------------------------- /
        IKN   = IRANGE(1,NSPEC,NANG)       ! Index vector for array access, e.g.
                                          ! in form of WN(1:NFRE) == WN2(IKN).
        ABAND = SUM(RESHAPE(A,(/ NANG,NFRE /)),1) ! action density as function of wavenumber
        DDIS  = 0.0_JWRB
        D     = 0.0_JWRB
        B1    = SWL6B1                         ! empirical constant from NAMELIST

!/ 1) --- Choose calculation of steepness a*k ------------------------ /
!/        Replace the measure of steepness with the spectral
!         saturation after Banner et al. (2002) ---------------------- /
        KK    = RESHAPE(A,(/ NANG,NFRE /))
        KMAX  = MAXVAL(KK,1)
        DO M = 1,NFRE
           IF (KMAX(M).LT.1.0E-34_JWRB) THEN
              KK(1:NANG,M) = 1.0_JWRB
           ELSE
              KK(1:NANG,M) = KK(1:NANG,M)/KMAX(M)
           END IF
        END DO
        ANAR  = 1.0_JWRB/( SUM(KK,1) * DELTH )
!        BN    = ANAR * ( ABAND * SIG * DELTH ) * WN(IJ,:)**3
        BN    = ANAR * ( ABAND * SIG * DELTH ) * WAVNUM(IJ,:)**3

!
        IF (.NOT.SWL6CSTB1) THEN
!
!/    --- A constant value for B1 attenuates swell too strong in the
!/        western central Pacific (i.e. cross swell less than 1.0m).
!/        Workaround is to scale B1 with steepness a*kp, where kp is
!/        the peak wavenumber. SWL6B1 remains a scaling constant, but
!/        with different magnitude.  --------------------------------- /
            M    = MAXLOC(ABAND,1)         ! Index for peak
!           EMEAN = SUM(ABAND * DDEN / CG)  ! Total sea surface variance
!            B1    = SWL6B1*(2.0_JWRB*SQRT(SUM(ABAND*DDEN/CGG(IJ,:)))*&
!            &       WN(IJ,M))
            B1    = SWL6B1*(2.0_JWRB*SQRT(SUM(ABAND*DDEN/CGROUP(IJ,:)))*&
            &       WAVNUM(IJ,M))

!
        END IF
!
!/ 2) --- Calculate the derivative term only (in units of 1/s) ------- /
        DO M = 1,NFRE
           IF (ABAND(M) .GT. 1.0E-30_JWRB) THEN
              DDIS(M) = -(2.0_JWRB/3.0_JWRB) * B1 * SIG(M) * SQRT(BN(M))
           END IF
        END DO
!
!/ 3) --- Apply dissipation term of derivative to all directions ----- /
        DO K = 1, NANG
           D(IKN+(K-1)) = DDIS
        END DO
!
        !S = D * A
!
!       WRITE(*,*) ' B1       =',B1
!       WRITE(*,*) ' DDIS_tot =',SUM(DDIS*ABAND*DDEN/CG)
!       WRITE(*,*) ' EDENS_tot=',sum(aband*dden/cg)
!       WRITE(*,*) ' EDENS_tot=',sum(aband*sig*dth*dsii/cg)
!       WRITE(*,*) ' '
!       WRITE(*,*) ' SWL6_tot =',sum(SUM(RESHAPE(S,(/ NANG,NFRE /)),1)*DDEN/CG)

        DSWL = RESHAPE(D,(/NANG,NFRE/))
        DO M = 1,NFRE
          DO K = 1, NANG
            SL(IJ,K,M)  = SL(IJ,K,M)  + DSWL(K,M)*FL1(IJ,K,M)
            FLD(IJ,K,M) = FLD(IJ,K,M) + DSWL(K,M)
          END DO
        END DO

      END DO
      ! END LOOP OVER LOC

      IF (LHOOK) CALL DR_HOOK('SWLDISSIP_BYDBR',1,ZHOOK_HANDLE)

      END SUBROUTINE SWLDISSIP_BYDBR
