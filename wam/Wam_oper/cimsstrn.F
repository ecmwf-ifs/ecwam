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

      USE YOWICE   , ONLY : CITHICK
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER
      USE YOWSHAL  , ONLY : TFAK     ,INDEP    ,DEPTHA,   DEPTHD
      USE YOWSTAT  , ONLY : ISHALLO
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK


! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IJ,M,K,IJS,IJL,IG,JD
      REAL :: AKI_ICE
      REAL :: F(IJS:IJL,NANG,NFRE)
      REAL, DIMENSION(IJS:IJL) :: STRN
      REAL, DIMENSION(IJS:IJL) :: DPTH, XK, XKI, E, F1DDF 
      REAL :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('CIMSSTRN',0,ZHOOK_HANDLE)

!*    1. INITIALISE
!        ----------

      IG=1

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
          F1DDF(IJ) = 0.
        ENDDO
        DO K=1,NANG
          DO IJ=IJS,IJL
            F1DDF(IJ) = F1DDF(IJ)+F(IJ,K,M)*DFIM(M)
          ENDDO
        ENDDO

        DO IJ=IJS,IJL
          STRN(IJ) = STRN(IJ)+E(IJ)**2*F1DDF(IJ)
        ENDDO

      ENDDO

      IF (LHOOK) CALL DR_HOOK('CIMSSTRN',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE CIMSSTRN
