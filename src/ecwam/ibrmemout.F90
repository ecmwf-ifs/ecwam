! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE IBRMEMOUT(KIJS, KIJL, IBRMEM, CICV, IBRMEMMSK)

! ----------------------------------------------------------------------

!**** *IBRMEMOUT* - MASK NON ICE POINTS IBRMEM.

!     J. BIDLOT  ECMWF  JANUARY 2013.

!*    PURPOSE.
!     --------

!

!**   INTERFACE.
!     ----------

!       *CALL* *IBRMEMOUT (KIJS, KIJL, IBRMEM, CICV, IBRMEMMSK)*
!              *KIJS*    - INDEX OF FIRST GRIDPOINT
!              *KIJL*    - INDEX OF LAST GRIDPOINT
!              *IBRMEM*  - MODEL IBRMEM
!              *CICV*    - SEA ICE COVER
!              *IBRMEMMSK* - IBRMEM MASKED WHERE NO SEA ICE PRESENT (OUTPUT)

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS , ONLY : ZMISS

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(IN)  :: IBRMEM
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN)     :: CICV
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(OUT) :: IBRMEMMSK


      INTEGER(KIND=JWIM) :: IJ

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('IBRMEMOUT',0,ZHOOK_HANDLE)

      DO IJ=KIJS,KIJL
        IBRMEMMSK(IJ) = IBRMEM(IJ)
      ENDDO

      DO IJ = KIJS,KIJL
        IF .NOT. (CICV(IJ) > 0.0_JWRB) IBRMEMMSK(IJ) = 0 ! 2=SOLID,1=BROKEN,0=MISSING
      ENDDO        

      IF (LHOOK) CALL DR_HOOK('IBRMEMOUT',1,ZHOOK_HANDLE)
      

      END SUBROUTINE IBRMEMOUT
