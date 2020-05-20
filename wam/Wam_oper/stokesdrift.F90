      SUBROUTINE STOKESDRIFT(F3, IJS, IJL, U10, THW, CICVR, USTOKES, VSTOKES)
 
!
!***  *STOKESDRIFT*   DETERMINES THE STOKES DRIFT
!
!     PETER JANSSEN MARCH 2009
!
!     PURPOSE.
!     --------
!
!              DETERMINATION OF STOKES DRIFT VECTOR
!
!     INTERFACE.
!     ----------
!              *CALL*  *STOKESDRIFT(F3,IJS,IJL,U10,THW,CICVR,USTOKES,VSTOKES)*
!
!                       INPUT:
!                            *F3*    - 2-D SPECTRUM
!                            *IJS*   - FIRST GRIDPOINT              
!                            *IJL*   - LAST GRIDPOINT              
!                            Auxilliary fields to specify Stokes when model sea ice cover the blocking threshold
!                            as 0.016*U10, aligned in the wind direction
!                            *U10*    - WIND SPEED IN M/S.
!                            *THW*    - WIND DIRECTION IN RADIANS.
!                            *CICVR*  - SEA ICE COVER.
!
!                       OUTPUT: 
!                            *USTOKES*   - U-COMPONENT STOKES DRIFT
!                            *VSTOKES*   - V-COMPONENT STOKES DRIFT
!
!     METHOD.
!     -------
!              DETERMINE U- AND V-COMPONENT OF STOKES DRIFT FOLLOWING
!              K.E. KENYON, J.G.R., 74, 6991-6994
!
!     EXTERNALS.
!     ----------
!              NONE
!
!
!-----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPCONS , ONLY : G        ,ZPI
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,TH       ,    &
     &                      DFIM_SIM ,FRATIO   ,COSTH    ,SINTH
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE ,CITHRSH
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_ODD
      USE YOWSHAL  , ONLY : TFAC_ST  ,INDEP
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
       
! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: U10, THW, CICVR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F3
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: USTOKES, VSTOKES


      INTEGER(KIND=JWIM) :: IJ, M, K

      REAL(KIND=JWRB), PARAMETER :: STMAX=1.5_JWRB   ! maximum magnitude (this is for safety when coupled)
      REAL(KIND=JWRB) :: CONST, FAC, FAC1, FAC2, FAC3
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: STFAC

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('STOKESDRIFT',0,ZHOOK_HANDLE)

 
!***  1. DETERMINE STOKE DRIFT VECTOR.
!     --------------------------------

      CONST = 2.0_JWRB*DELTH*ZPI**3/G*FR(NFRE_ODD)**4

!***  1.1 PERFORM INTEGRATION.
!     ------------------------
 
      DO IJ = IJS,IJL
         USTOKES(IJ) = 0.0_JWRB
         VSTOKES(IJ) = 0.0_JWRB
      ENDDO

      DO M=1,NFRE_ODD
         DO IJ = IJS,IJL
           STFAC(IJ) = TFAC_ST(INDEP(IJ),M)*DFIM_SIM(M)
         ENDDO
         DO K=1,NANG
            DO IJ = IJS,IJL
               FAC3 = STFAC(IJ)*F3(IJ,K,M)
               USTOKES(IJ) = USTOKES(IJ)+FAC3*SINTH(K)
               VSTOKES(IJ) = VSTOKES(IJ)+FAC3*COSTH(K)
            ENDDO
         ENDDO
      ENDDO
 
!***  1.2 ADD CONTRIBUTION OF UNRESOLVED WAVES.
!     -----------------------------------------
 
      DO K=1,NANG
         FAC1 = CONST*SINTH(K)
         FAC2 = CONST*COSTH(K)
         DO IJ = IJS,IJL
            USTOKES(IJ) = USTOKES(IJ)+FAC1*F3(IJ,K,NFRE_ODD)
            VSTOKES(IJ) = VSTOKES(IJ)+FAC2*F3(IJ,K,NFRE_ODD)
         ENDDO
      ENDDO


!***  1.3 Sea Ice exception
!     ---------------------
      IF (LICERUN .AND. LMASKICE) THEN
       DO IJ=IJS,IJL
         IF(CICVR(IJ) .GT. CITHRSH) THEN
           USTOKES(IJ) = 0.016_JWRB*U10(IJ)*SIN(THW(IJ))*(1.0_JWRB - CICVR(IJ))
           VSTOKES(IJ) = 0.016_JWRB*U10(IJ)*COS(THW(IJ))*(1.0_JWRB - CICVR(IJ))
         ENDIF
       ENDDO
     ENDIF

!***  1.4 Protection
!     --------------

      DO IJ = IJS,IJL
         USTOKES(IJ) = MIN(MAX(USTOKES(IJ),-STMAX),STMAX)
         VSTOKES(IJ) = MIN(MAX(VSTOKES(IJ),-STMAX),STMAX)
      ENDDO

      IF (LHOOK) CALL DR_HOOK('STOKESDRIFT',1,ZHOOK_HANDLE)

      END SUBROUTINE STOKESDRIFT
