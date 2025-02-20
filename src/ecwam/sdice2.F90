! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SDICE2 (KIJS, KIJL, FL1, FLD, SL, SLICE,     &
     &                   WAVNUM, CGROUP,                      &
     &                   CICV)
! ----------------------------------------------------------------------

!**** *SDICE2* - COMPUTE SEA ICE WAVE ATTENUATION FACTORS DUE TO ICE FLOES
!                BOTTOM FRICTION (CAME FROM CIWABR)


!     JEAN BIDLOT       ECMWF ~2012
!     JOSH KOUSAL       ECMWF 2023


!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!       *CALL* *SDICE2 (KIJS, KIJL, FL1, FLD,SL,*
!                       WAVNUM, CGROUP,
!                       CICV)*
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *FL1*    - SPECTRUM.
!          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*     - TOTAL SOURCE FUNCTION ARRAY
!          *SLICE*  - TOTAL SOURCE FUNCTION ARRAY, ICE
!          *WAVNUM* - WAVE NUMBER
!          *CGROUP* - GROUP SPEED
!          *CICV*   - SEA ICE COVER

!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!     KOHOUT A., M. MEYLAN, D PLEW, 2011: ANNALS OF GLACIOLOGY, 2011. 
!     M.J. Doble, J.-R. Bidlot / Ocean Modelling 70 (2013), 166-173

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : DFIM
      USE YOWICE   , ONLY : CDICWA  ,ZALPFACB
      USE YOWPARAM , ONLY : NANG    ,NFRE
      USE YOWPCONS , ONLY : EPSMIN  
      USE YOWSTAT  , ONLY : IDELT

      USE YOMHOOK  , ONLY : LHOOK   ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: FLD, SL
      REAL(KIND=JWRB), DIMENSION(KIJS,NANG,NFRE), INTENT(OUT) ::        SLICE
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: WAVNUM, CGROUP
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: CICV

      REAL(KIND=JWRB), DIMENSION(NFRE)    :: XK2

      INTEGER(KIND=JWIM) :: IJ, K, M
      REAL(KIND=JWRB)    :: EWH
      REAL(KIND=JWRB)    :: ALP              !! ALP=SPATIAL ATTENUATION RATE OF ENERGY
      REAL(KIND=JWRB)    :: FLDICE
      REAL(KIND=JWRB)    :: DELT, DELTM, XIMP, DELT5, GTEMP1
      
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SDICE2',0,ZHOOK_HANDLE)

      DELT = IDELT
      DELTM = 1.0_JWRB/DELT
      XIMP = 1.0_JWRB
      DELT5 = XIMP*DELT

      DO M = 1,NFRE
         DO K = 1,NANG
            DO IJ = KIJS,KIJL

               EWH            = 4.0_JWRB*SQRT(MAX(EPSMIN,FL1(IJ,K,M)*DFIM(M)))
               XK2(M)         = WAVNUM(IJ,M)**2
               ALP            = CDICWA*XK2(M)*EWH*ZALPFACB

!              apply the source term
               FLDICE         = -ALP(IJ,M)   * CGROUP(IJ,M)   
               SLICE(IJ,K,M)  =  FL1(IJ,K,M) * FLDICE
               SL(IJ,K,M)     =  SL(IJ,K,M)  + CICV(IJ)*SLICE(IJ,K,M)
               FLD(IJ,K,M)    =  FLD(IJ,K,M) + CICV(IJ)*FLDICE
               
!              to be used for wave radiative stress calculation
               GTEMP1         =  MAX((1.0_JWRB-DELT5*FLDICE),1.0_JWRB)    
               SLICE(IJ,K,M)  =  SLICE(IJ,K,M)/GTEMP1

            END DO
         END DO
      END DO
      
      IF (LHOOK) CALL DR_HOOK('SDICE2',1,ZHOOK_HANDLE)

      END SUBROUTINE SDICE2
