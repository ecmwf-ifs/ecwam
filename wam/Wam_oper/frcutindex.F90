      SUBROUTINE FRCUTINDEX (IJS, IJL, FM, FMWS, USNEW, CICVR,          &
     &                       MIJ, MIJFLX, RHOWGDFTH)

! ----------------------------------------------------------------------

!**** *FRCUTINDEX* - RETURNS THE LAST FREQUENCY INDEX OF
!                    PROGNOSTIC PART OF SPECTRUM.

!**   INTERFACE.
!     ----------

!       *CALL* *FRCUTINDEX (IJS, IJL, FM, FMWS, CICVR, MIJ, MIJFLX, RHOWGDFTH)
!          *IJS*    - INDEX OF FIRST GRIDPOINT
!          *IJL*    - INDEX OF LAST GRIDPOINT
!          *FM*     - MEAN FREQUENCY
!          *FMWS*   - MEAN FREQUENCY OF WINDSEA
!          *USNEW*  - FRICTION VELOCITY IN M/S
!          *CICVR*  - CICVR 
!          *MIJ*    - LAST FREQUENCY INDEX for imposing high frequency tail
!          *MIJFLX* - LAST FREQUENCY INDEX for computing the fluxes
!          *RHOWGDFTH - WATER DENSITY * G * DF * DTHETA
!                       FOR TRAPEZOIDAL INTEGRATION BETWEEN FR(1) and FR(MIJ) 
!                       !!!!!!!!  RHOWGDFTH=0 FOR FR > FR(MIJ)


!     METHOD.
!     -------

!*    COMPUTES LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!*    FREQUENCIES LE 2.5*MAX(FMWS,FM).


!!! be aware that if this is not used, for iphys=1, the cumulative dissipation has to be
!!! re-activated (see module yowphys) !!!


!     EXTERNALS.
!     ---------

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------

!     EXTERNALS.
!     ---------

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,DFIM       ,FRATIO   ,FLOGSPRDM1, &
     &                ZPIFR,                                            &
     &                DELTH          ,RHOWG_DFIM ,FRIC
      USE YOWICE   , ONLY : CITHRSH_TAIL
      USE YOWPARAM , ONLY : NFRE
      USE YOWPCONS , ONLY : G        ,EPSMIN
      USE YOWPHYS  , ONLY : TAILFACTOR, TAILFACTOR_PM
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM), INTENT(OUT) :: MIJ(IJS:IJL)
      INTEGER(KIND=JWIM), INTENT(OUT) :: MIJFLX(IJS:IJL)

      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(IN) :: FM, FMWS, USNEW, CICVR
      REAL(KIND=JWRB),DIMENSION(IJS:IJL,NFRE), INTENT(OUT) :: RHOWGDFTH 

      INTEGER(KIND=JWIM) :: IJ, M

      REAL(KIND=JWRB) :: FPMH, FPPM, FM2, FPM, FPM4
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('FRCUTINDEX',0,ZHOOK_HANDLE)

!*    COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!*    FREQUENCIES LE MAX(TAILFACTOR*MAX(FMNWS,FM),TAILFACTOR_PM*FPM),
!*    WHERE FPM IS THE PIERSON-MOSKOWITZ FREQUENCY BASED ON FRICTION
!*    VELOCITY. (FPM=G/(FRIC*ZPI*USTAR))
!     ------------------------------------------------------------

      FPMH = TAILFACTOR/FR(1)
      FPPM = TAILFACTOR_PM*G/(FRIC*ZPIFR(1))

      MIJFLX(:) = NFRE
      DO IJ=IJS,IJL
        IF (CICVR(IJ).LE.CITHRSH_TAIL) THEN
          FM2 = MAX(FMWS(IJ),FM(IJ))*FPMH
          FPM = FPPM/MAX(USNEW(IJ),EPSMIN)
          FPM4 = MAX(FM2,FPM)
          MIJ(IJ) = NINT(LOG10(FPM4)*FLOGSPRDM1)+1
          MIJ(IJ) = MIN(MAX(1,MIJ(IJ)),NFRE)
        ELSE
          MIJ(IJ) = NFRE
        ENDIF
      ENDDO

!     SET RHOWGDFTH
      DO IJ=IJS,IJL
        DO M=1,MIJFLX(IJ)
          RHOWGDFTH(IJ,M) = RHOWG_DFIM(M)
        ENDDO
        IF(MIJFLX(IJ).NE.NFRE) RHOWGDFTH(IJ,MIJFLX(IJ))=0.5_JWRB*RHOWGDFTH(IJ,MIJFLX(IJ))
        DO M=MIJFLX(IJ)+1,NFRE
          RHOWGDFTH(IJ,M) = 0.0_JWRB
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('FRCUTINDEX',1,ZHOOK_HANDLE)

      END SUBROUTINE FRCUTINDEX
