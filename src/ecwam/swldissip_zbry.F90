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
     &                        UFRIC, COSWDIF, RAORW)
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

!       *CALL* *SWLDISSIP_ZBRY (KIJS, KIJL, FL1, FLD,SL,*
!                            WAVNUM, CGROUP, 
!                            UFRIC, COSWDIF, RAORW)*
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *FL1*    - SPECTRUM.
!          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*     - TOTAL SOURCE FUNCTION ARRAY
!          *WAVNUM* - WAVE NUMBER
!          *CGROUP* - GROUP SPEED
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
!     Adapted from Babanin Young Donelan & Banner (ZBRY) physics 
!     as implemented as ST6 in WAVEWATCH-III
!     WW3 module:       W3SWLDMD
!     WW3 subroutine:   W3SWL6
!     Implementation into ECWAM DECEMBER 2021 by J. Kousal 


! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR      , TH     ,ZPIFR   ,FRATIO    ,DELTH, DFIM,
      &                     SIG     , DDEN
      USE YOWPCONS , ONLY : G        ,ZPI
      USE YOWPARAM , ONLY : NANG    ,NFRE
      USE YOWPHYS  , ONLY : LLSWL6CSTB1, ZSWL6B1

      USE YOMHOOK  , ONLY : LHOOK   ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "irange.intfb.h"      

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: FLD, SL
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: WAVNUM, CGROUP
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: UFRIC, RAORW 
      REAL(KIND=JWRB), DIMENSION(KIJL, NANG), INTENT(IN) :: COSWDIF 

      INTEGER(KIND=JWIM) :: IJ, M, I, J, M2, K2, K, NANGD
      INTEGER(KIND=JWIM) :: NSPEC !num. of freqs, dirs, spec. bins
      INTEGER(KIND=JWIM), DIMENSION(NANG) :: ITHN
      INTEGER(KIND=JWIM), DIMENSION(NFRE) :: IKN

      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: ABAND, KMAX, ANAR, BN, DDIS
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE) :: KK
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG*NFRE) :: S, D, A, CG2
      REAL(KIND=JWRB), DIMENSION(KIJL)      :: B1
      REAL(KIND=JWRB), DIMENSION(NANG*NFRE) :: SIG2 
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE) :: DSWL

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      
! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SWLDISSIP_ZBRY',0,ZHOOK_HANDLE)

      NSPEC = NANG * NFRE   ! NUMBER OF SPECTRAL BINS

      IKN    = IRANGE(1,NSPEC,NANG)   ! Index vector for elements of 1 ... NFRE
!                                     ! such that e.g. SIG(1:NFRE) = SIG2(IKN).
      DO K = 1, NANG                    ! Apply to all directions 
         SIG2   (IKN+(K-1)) = SIG
      END DO
      
      DO K = 1, NANG
        DO IJ = KIJS,KIJL
           CG2    (IJ,IKN+(K-1)) = CGROUP(IJ,:)
        END DO
      END DO

      DO IJ = KIJS,KIJL
        A(IJ,:) = RESHAPE( FL1(IJ,:,:) , (/NSPEC/)) * CG2(IJ,:) / ( ZPI * SIG2 )! ACTION DENSITY SPECTRUM
        ! WAM E(f,theta) to WW3 A(k,theta) conversion factor: CG2 / ( ZPI *SIG2 ) 
      END DO
      
      !/ 0) --- Initialize parameters -------------------------------------- /
      IKN   = IRANGE(1,NSPEC,NANG)       ! Index vector for array access, e.g.
                                         ! in form of WN(1:NFRE) == WN2(IKN).

      DO IJ = KIJS,KIJL
        ABAND(IJ,:) = SUM(RESHAPE(A(IJ,:),(/ NANG,NFRE /)),1) ! action density as function of wavenumber
        DDIS(IJ,:)  = 0.0_JWRB
        D(IJ,:)     = 0.0_JWRB
      END DO

!/ 1) --- Choose calculation of steepness a*k ------------------------ /
!/        Replace the measure of steepness with the spectral
!         saturation after Banner et al. (2002) ---------------------- /
      DO IJ = KIJS,KIJL
        KK(IJ,:,:)    = RESHAPE(A(IJ,:),(/ NANG,NFRE /))
        KMAX(IJ,:)  = MAXVAL(KK(IJ,:,:),1)
      END DO

      DO M = 1,NFRE
        DO IJ = KIJS,KIJL
           IF (KMAX(IJ,M).LT.1.0E-34_JWRB) THEN
              KK(IJ,1:NANG,M) = 1.0_JWRB
           ELSE
              KK(IJ,1:NANG,M) = KK(IJ,1:NANG,M)/KMAX(IJ,M)
           END IF
        END DO
      END DO

      DO IJ = KIJS,KIJL
        ANAR(IJ,:)  = 1.0_JWRB/( SUM(KK(IJ,:,:),1) * DELTH )
        BN(IJ,:)    = ANAR(IJ,:) * ( ABAND(IJ,:) * SIG * DELTH ) * WAVNUM(IJ,:)**3
      END DO

!
      IF (.NOT.LLSWL6CSTB1) THEN
!/    --- A constant value for B1 attenuates swell too strong in the
!/        western central Pacific (i.e. cross swell less than 1.0m).
!/        Workaround is to scale B1 with steepness a*kp, where kp is
!/        the peak wavenumber. ZSWL6B1 remains a scaling constant, but
!/        with different magnitude.  --------------------------------- /
        DO IJ = KIJS,KIJL
            M      = MAXLOC(ABAND(IJ,:),1)         ! Index for peak
            B1(IJ) = ZSWL6B1*(2.0_JWRB*SQRT(SUM(ABAND(IJ,:)*DDEN/CGROUP(IJ,:)))*&
            &        WAVNUM(IJ,M))
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
      DO K = 1, NANG
        DO IJ = KIJS,KIJL
          D(IJ,IKN+(K-1)) = DDIS(IJ,:)
        END DO
      END DO
!
      DO IJ = KIJS,KIJL
        DSWL(IJ,:,:) = RESHAPE(D(IJ,:),(/NANG,NFRE/))
      END DO   
 
      DO M = 1,NFRE
        DO K = 1, NANG
          DO IJ = KIJS,KIJL
             SL(IJ,K,M)  = SL(IJ,K,M)  + DSWL(IJ,K,M)*FL1(IJ,K,M)
             FLD(IJ,K,M) = FLD(IJ,K,M) + DSWL(IJ,K,M)
          END DO  
        END DO
      END DO


      IF (LHOOK) CALL DR_HOOK('SWLDISSIP_ZBRY',1,ZHOOK_HANDLE)

      END SUBROUTINE SWLDISSIP_ZBRY
