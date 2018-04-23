!-----------------------------------------------------------------------

      SUBROUTINE SETICE (FL3, IJS, IJL, CICVR, U10NEW, THWNEW)

!-----------------------------------------------------------------------

!**** *SETICE* ROUTINE TO SET SPECTRA ON ICE TO NOISE LEVEL.

!     R.PORTZ      MPI         OKT.1992
!     J. BIDLOT    ECMWF       JUNE 1996  MESSAGE PASSING

!     PURPOSE.
!     -------

!          *SETICE* SET ICE SPECTRA (FL3) TO NOISE LEVEL 

!**   INTERFACE.
!     ----------

!         *CALL* *SETICE*

!     METHOD.
!     -------

!          NONE.


!     EXTERNALS.
!     ----------

!          .NONE.

!     REFERENCE.
!     ----------

!          NONE.

! ----------------------------------------------------------------------

      USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK
      USE YOWFRED  , ONLY : TH
      USE YOWICE   , ONLY : FLMINFR  ,CITHRSH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : EPSMIN
      USE YOWTABL  , ONLY : JUMAX    ,DELU

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL
      REAL, DIMENSION(IJS:IJL), INTENT(IN) :: CICVR, U10NEW, THWNEW
      REAL, DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: FL3

      INTEGER :: IJ, M, K
      INTEGER, DIMENSION(IJS:IJL) :: JU 

      REAL :: XJ
      REAL :: ZHOOK_HANDLE
      REAL, DIMENSION(IJS:IJL) :: CIREDUC, TEMP, ICEFREE 
      REAL, DIMENSION(IJS:IJL,NANG) :: SPRD
! ----------------------------------------------------------------------

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SETICE',0,ZHOOK_HANDLE)
#endif

!*    1. SET SPECTRA TO NOISE LEVEL OVER ICE POINTS.
!     ----------------------------------------------

      DO K=1,NANG
        DO IJ=IJS,IJL
          SPRD(IJ,K)=MAX(0.,COS(TH(K)-THWNEW(IJ)))**2
        ENDDO
      ENDDO

      DO IJ = IJS,IJL
        XJ=U10NEW(IJ)/DELU
        JU(IJ)=MIN(JUMAX, MAX(NINT(XJ),1))
        IF(CICVR(IJ).GT.CITHRSH) THEN
          CIREDUC(IJ)=MAX(EPSMIN,(1.-CICVR(IJ)))
          ICEFREE(IJ)=0.0
        ELSE
          CIREDUC(IJ)=0.0
          ICEFREE(IJ)=1.0
        ENDIF
      ENDDO

      DO M = 1, NFRE
        DO IJ = IJS,IJL
          TEMP(IJ)=CIREDUC(IJ)*FLMINFR(JU(IJ),M)
        ENDDO
        DO K = 1, NANG
          DO IJ = IJS,IJL
            FL3(IJ,K,M)=FL3(IJ,K,M)*ICEFREE(IJ)+TEMP(IJ)*SPRD(IJ,K)
          ENDDO
        ENDDO
      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SETICE',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE SETICE
