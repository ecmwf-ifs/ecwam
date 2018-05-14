      SUBROUTINE OUTINT

! ----------------------------------------------------------------------

!**** *OUTINT* - OUTPUT OF INTEGRATED FIELDS WITHOUT THE I/O SERVER

!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!       *CALL* *OUTINT

!     EXTERNALS.
!     ----------
!       *OUTGRID*
!       *WGRIBENOUT*

!     METHOD.
!     -------

!       NONE.

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
  
      REAL(KIND=JWRB) :: ZHOOK_HANDLE


! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('OUTINT',0,ZHOOK_HANDLE)

      write(*,*) ' hello world !'

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('OUTINT',1,ZHOOK_HANDLE)

      END SUBROUTINE OUTINT
