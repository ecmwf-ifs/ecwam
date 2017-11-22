      SUBROUTINE SDEPTHLIM(IJS, IJL, EMAXDPT, FL)
! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF   NOVEMBER 2017

!*    PURPOSE.
!     --------
!     LIMITS THE SPECTRAL VARIANCE SUCH THAT THE TOTAL VARIANCE
!     DOES NOT EXCEED THE MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH

!**   INTERFACE.
!     ----------
!     *CALL* *SDEPTHLIM((IJS, IJL, EMAXDPT, FL)
!          *IJS* - INDEX OF FIRST GRIDPOINT.
!          *IJL* - INDEX OF LAST GRIDPOINT.
!          *EMAXDPT - MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH
!          *FL*  - SPECTRUM.


!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!     NONE

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : EPSMIN
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: EMAXDPT
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: FL

      INTEGER(KIND=JWIM) :: IJ, K, M
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: EM
      LOGICAL :: LLEPSMIN

! ----------------------------------------------------------------------

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SDEPTHLIM',0,ZHOOK_HANDLE)
#endif

      LLEPSMIN=.TRUE.
      CALL SEMEAN (FL, IJS, IJL, EM, LLEPSMIN)

      DO IJ=IJS,IJL
        EM(IJ)=MIN(EMAXDPT(IJ)/EM(IJ),1.0_JWRB)
      ENDDO

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            FL(IJ,K,M) = MAX(FL(IJ,K,M)*EM(IJ),EPSMIN) 
          ENDDO
        ENDDO
      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SDEPTHLIM',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE SDEPTHLIM
