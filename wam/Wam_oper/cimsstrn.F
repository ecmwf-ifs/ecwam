      SUBROUTINE CIMSSTRN(F, IJS, IJL, STRN)

! ----------------------------------------------------------------------

!**** *CIMSSTRN* - COMPUTATION OF THE MEAN SQUARE WAVE STRAIN IN SEA ICE.

!     J. BIDLOT  ECMWF  JANUARY 2013.

!*    PURPOSE.
!     --------

!       COMPUTES MEAN SQUARE WAVE STRAIN AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *CIMSSTRN (F, IJS, IJL, STRN)*
!              *F*   - SPECTRUM.
!              *IJS* - INDEX OF FIRST GRIDPOINT
!              *IJL* - INDEX OF LAST GRIDPOINT
!              *STRN*- MEAN SQUARE WAVE STRAIN IN ICE (OUTPUT).

!     METHOD.
!     -------

!      !!! IT ASSUMES SO DEFAULT SETTING FOR THE MECHANICAL PROPERTIES OF
!          THE SEA ICE (SEE AKI_ICE) !!!!!!!

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM,JWRB,JWRU
      USE YOWICE   , ONLY : CITHICK  ,FLMIN
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER
      USE YOWSHAL  , ONLY : TFAK     ,INDEP    ,DEPTHA,   DEPTHD
      USE YOWSTAT  , ONLY : ISHALLO
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK


! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS,IJL
      INTEGER(KIND=JWIM) :: IJ,M,K,IG,JD
      REAL(KIND=JWRB) :: AKI_ICE
      REAL(KIND=JWRB) :: F1LIM 
      REAL(KIND=JWRB), INTENT(IN) :: F(IJS:IJL,NANG,NFRE)
      REAL(KIND=JWRU), DIMENSION(IJS:IJL), INTENT(OUT) :: STRN ! Has to be REAL*8
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: DPTH, XK, XKI, E, SUME 
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('CIMSSTRN',0,ZHOOK_HANDLE)
#endif

!*    1. INITIALISE
!        ----------

      IG=1
      F1LIM=FLMIN/DELTH

      DO IJ=IJS,IJL
        STRN(IJ) = 0.0
      ENDDO

! ----------------------------------------------------------------------

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      DO M=1,NFRE
        IF (ISHALLO.EQ.1) THEN
          DO IJ=IJS,IJL
            DPTH(IJ)=999.
            XK(IJ)=(ZPI*FR(M))**2/G
          ENDDO
        ELSE
          DO IJ=IJS,IJL
            JD=INDEP(IJ)
            DPTH(IJ)=DEPTHA*DEPTHD**(JD-1)
            XK(IJ)=TFAK(JD,M)
          ENDDO
        ENDIF

        DO IJ=IJS,IJL
          XKI(IJ)=AKI_ICE(G,XK(IJ),DPTH(IJ),ROWATER,CITHICK(IJ,IG))
          E(IJ)=0.5*CITHICK(IJ,IG)*XKI(IJ)**3/XK(IJ)
        ENDDO

        DO IJ=IJS,IJL
          SUME(IJ) = 0.
        ENDDO
        DO K=1,NANG
          DO IJ=IJS,IJL
            SUME(IJ) = SUME(IJ)+F(IJ,K,M)
          ENDDO
        ENDDO

        DO IJ=IJS,IJL
          IF(SUME(IJ).GT.F1LIM) THEN
            STRN(IJ) = STRN(IJ)+E(IJ)**2*SUME(IJ)*DFIM(M)
          ENDIF
        ENDDO

      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('CIMSSTRN',1,ZHOOK_HANDLE)
#endif
      RETURN
      END SUBROUTINE CIMSSTRN
