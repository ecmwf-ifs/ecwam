      SUBROUTINE OUTWINT(IJS, IJL, NBOUT, BOUT)

! ----------------------------------------------------------------------

!**** *OUTWINT* -

!*    PURPOSE.
!     --------


!**   INTERFACE.
!     ----------

!     METHOD.
!     -------

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : JPPFLAG ,ITOBOUT
      USE YOWTEST  , ONLY : IU06
      USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, NBOUT
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NBOUT), INTENT(IN) :: BOUT

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTWINT',0,ZHOOK_HANDLE)
#endif

!!!! if I/O server do something



#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTWINT',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE OUTWINT
