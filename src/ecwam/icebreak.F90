! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE ICEBREAK (KIJS, KIJL, EMEAN, AKMEAN,                   &
     &                     CITH, IBR_MEM, ALPFAC)
! ----------------------------------------------------------------------

!**** *ICEBREAK* - COMPUTATION OF BREAK UP OF SEA ICE BY WAVES


!     JOSH KOUSAL       ECMWF 2023


!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!       *CALL* *ICEBREAK (KIJS, KIJL, EMEAN, AKMEAN, 
!                         CITH, IBR_MEM)*
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *EMEAN*  - MEAN ENERGY DENSITY 
!          *AKMEAN* - MEAN WAVE NUMBER  BASED ON sqrt(1/k)*F INTGRATION  (TODO: XKMEAN? OR ?)
!          *CITH*   - SEA ICE THICKNESS
!          *IBR_MEM*- ICEBREAK MEMORY                      
!          *ALPFAC* - FACTOR TO REDUCE ATTENUATION IN ICE. EQUAL TO 1 IN SOLID ICE

!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!     Kousal, J., Voermans, J. J., Liu, Q., Heil, P., & Babanin, A. V.
!     (2022). A two-part model for wave-sea ice interaction: Attenuation
!     and break-up. JGRO. DOI:10.1029/2022JC018571

!     Voermans, J. J., Rabault, J., Filchuk, K., Ryzhov, I., Heil, P., 
!     Marchenko, A., et al. (2020). Experimental evidence for a
!     universal threshold characterizing wave-induced sea ice break-up.
!     Cryosphere DOI:10.5194/tc-14-4265-2020

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS , ONLY : ZPI

      USE YOWTEST  , ONLY : IU06

      USE YOMHOOK  , ONLY : LHOOK   ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: CITH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: EMEAN, AKMEAN

      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: IBR_MEM
      REAL(KIND=JWRB)   , DIMENSION(KIJS:KIJL), INTENT(INOUT) :: ALPFAC

      INTEGER(KIND=JWIM) :: IJ, K, M
      REAL(KIND=JWRB)    :: IBR_CONST2, IBR_CONST3, IBR_CONST4
      REAL(KIND=JWRB)    :: IBR_CONST5, IBR
      REAL(KIND=JWRB)    :: HS, LAMBDASQ
      
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('ICEBREAK',0,ZHOOK_HANDLE)

      IBR_CONST2 = 6E+9_JWRB     ! Young's modulus          (Y=6GPa)
      !IBR_CONST2 = 0.2E+9_JWRB     ! custom value (strong ice bound)
      !IBR_CONST2 = 9E+9_JWRB     ! custom value (weak ice bound)

      IBR_CONST3 = 0.55E+6_JWRB ! flexural strength of ice (sig=0.55MPa)
      !IBR_CONST3 = 0.7E+6_JWRB ! custom value (strong ice bound)
      !IBR_CONST3 = 0.1E+6_JWRB ! custom value (weak ice bound)           

      IBR_CONST4 = 0.014_JWRB   ! non-dimensional threshold as 
                                ! in Voermans et al. 2020  (I_br=0.014)
      IBR_CONST5 = 0.01_JWRB    ! FACTOR TO DECREASE ATTENUATION BY WHEN ICE IS BROKEN

      DO IJ = KIJS,KIJL

        HS       =  4._JWRB*SQRT(EMEAN(IJ)) !TODO:check
        LAMBDASQ =  ZPI / (AKMEAN(IJ)**2)   !TODO:check

        ! ICE BREAK UP PARAMETER
        IBR      =  HS*CITH(IJ)*IBR_CONST2 / (2._JWRB*IBR_CONST3*LAMBDASQ)

        ! BREAK ICE IF IBR EXCEEDS THRESHOLD
        IF (IBR >= IBR_CONST4) THEN
          IBR_MEM(IJ) = 1
          ALPFAC(IJ)  = IBR_CONST5 
        ENDIF

      ENDDO

      
      IF (LHOOK) CALL DR_HOOK('ICEBREAK',1,ZHOOK_HANDLE)

      END SUBROUTINE ICEBREAK
