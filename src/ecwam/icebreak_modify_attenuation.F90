! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE ICEBREAK_MODIFY_ATTENUATION (KIJS, KIJL, IBRMEM, ALPFAC)
! ----------------------------------------------------------------------

!**** *ICEBREAK_MODIFY_ATTENUATION* - MODIFY ATTENUATION DUE TO ICE BREAKUP 


!     JOSH KOUSAL       ECMWF 2023


!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!       *CALL* *ICEBREAK_MODIFY_ATTENUATION (KIJS, KIJL, IBRMEM, ALPFAC)*
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *IBRMEM* - ICEBREAK MEMORY                      
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

      USE YOWICE   , ONLY : ZALPFACX, ZIBRW_THRSH

      USE YOWTEST  , ONLY : IU06

      USE YOMHOOK  , ONLY : LHOOK   ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB)   , DIMENSION(KIJL), INTENT(IN) :: IBRMEM
      REAL(KIND=JWRB)   , DIMENSION(KIJL), INTENT(INOUT) :: ALPFAC

      INTEGER(KIND=JWIM) :: IJ
      REAL(KIND=JWRB)    :: ZALPFACXINV
      
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('ICEBREAK_MODIFY_ATTENUATION',0,ZHOOK_HANDLE)

      ZALPFACXINV = 1.0_JWRB/ZALPFACX    ! FACTOR TO DECREASE ATTENUATION BY WHEN ICE IS BROKEN

      DO IJ = KIJS,KIJL

        ! USE IBRMEM TO CHANGE ATTENUATION FACTOR ALPFAC
        ! OTHERWISE ALPFAC KEEPS VALUE SET IN IMPLSCH
 
        IF  (IBRMEM(IJ) <= ZIBRW_THRSH) THEN
          ALPFAC(IJ)  = ZALPFACXINV
        ENDIF

      ENDDO

      IF (LHOOK) CALL DR_HOOK('ICEBREAK_MODIFY_ATTENUATION',1,ZHOOK_HANDLE)

      END SUBROUTINE ICEBREAK_MODIFY_ATTENUATION
