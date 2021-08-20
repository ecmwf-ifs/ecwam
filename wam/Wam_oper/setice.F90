!-----------------------------------------------------------------------

      SUBROUTINE SETICE (KIJS, KIJL, FL1, CICOVER, WSWAVE, WDWAVE)

!-----------------------------------------------------------------------

!**** *SETICE* ROUTINE TO SET SPECTRA ON ICE TO NOISE LEVEL.

!     R.PORTZ      MPI         OKT.1992
!     J. BIDLOT    ECMWF       JUNE 1996  MESSAGE PASSING

!     PURPOSE.
!     -------

!          *SETICE* SET ICE SPECTRA (FL1) TO NOISE LEVEL 

!**   INTERFACE.
!     ----------

!         *CALL* *SETICE(KIJS, KIJL, FL1, CICOVER, WSWAVE, WDWAVE)*
!          *KIJS*    - LOCAL INDEX OF FIRST GRIDPOINT
!          *KIJL*    - LOCAL  INDEX OF LAST GRIDPOINT
!          *FL1*     - SPECTRA
!          *CICOVER* - SEA ICE COVER
!          *WSWAVE*  - WIND SPEED.
!          *WDWAVE*  - WIND DIRECTION (RADIANS).


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : TH
      USE YOWICE   , ONLY : FLMIN    ,CITHRSH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : EPSMIN

      USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG, NFRE), INTENT(INOUT) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: CICOVER, WSWAVE, WDWAVE


      INTEGER(KIND=JWIM) :: IJ, M, K

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: CIREDUC, TEMP, ICEFREE 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG) :: SPRD
! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SETICE',0,ZHOOK_HANDLE)

!*    1. SET SPECTRA TO NOISE LEVEL OVER ICE POINTS.
!     ----------------------------------------------

      DO K=1,NANG
        DO IJ = KIJS, KIJL
          SPRD(IJ,K)=MAX(0.0_JWRB,COS(TH(K)-WDWAVE(IJ)))**2
        ENDDO
      ENDDO

      DO IJ = KIJS,KIJL
        IF (CICOVER(IJ) > CITHRSH) THEN
          CIREDUC(IJ)=MAX(EPSMIN,(1.0_JWRB-CICOVER(IJ)))
          ICEFREE(IJ)=0.0_JWRB
        ELSE
          CIREDUC(IJ)=0.0_JWRB
          ICEFREE(IJ)=1.0_JWRB
        ENDIF
      ENDDO

      DO IJ = KIJS,KIJL
        TEMP(IJ)=CIREDUC(IJ)*FLMIN
      ENDDO
      DO M = 1, NFRE
        DO K = 1, NANG
          DO IJ = KIJS,KIJL
            FL1(IJ,K,M)=FL1(IJ,K,M)*ICEFREE(IJ)+TEMP(IJ)*SPRD(IJ,K)
          ENDDO
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('SETICE',1,ZHOOK_HANDLE)

      END SUBROUTINE SETICE
