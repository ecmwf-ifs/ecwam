      SUBROUTINE FRCUTINDEX (IJS, IJL, FM, FMWS, USNEW, CICVR,          &
     &                       MIJ, MIJFLX, RHOWGDFTH)

! ----------------------------------------------------------------------

!**** *FRCUTINDEX* - RETURNS THE LAST FREQUENCY INDEX OF
!                    PROGNOSTIC PART OF SPECTRUM.

!**   INTERFACE.
!     ----------

!       *CALL* *FRCUTINDEX (IJS, IJL, FM, FMWS, CICVR, MIJ)
!          *IJS*    - INDEX OF FIRST GRIDPOINT
!          *IJL*    - INDEX OF LAST GRIDPOINT
!          *FM*     - MEAN FREQUENCY
!          *FMWS*   - MEAN FREQUENCY OF WINDSEA
!          *USNEW*  - FRICTION VELOCITY IN M/S
!          *CICVR*  - CICVR 
!          *MIJ*    - LAST FREQUENCY INDEX used to impose the high frequency tail
!          *MIJFLX* - LAST FREQUENCY INDEX used to compute the fluxes
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
      USE YOWPHYS  , ONLY : TAILFACTOR, TAILFACTOR_PM, TAILFACTOR_FLX
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM), INTENT(OUT) :: MIJ(IJS:IJL)
!!! test a different MIJ for all fluxes calculation
      INTEGER(KIND=JWIM), INTENT(OUT) :: MIJFLX(IJS:IJL)

      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(IN) :: FM, FMWS, USNEW, CICVR
      REAL(KIND=JWRB),DIMENSION(IJS:IJL,NFRE), INTENT(OUT) :: RHOWGDFTH 

      INTEGER(KIND=JWIM) :: IJ, M

      REAL(KIND=JWRB) :: FPMH, FPPM, FM2, FPM, FPM4, FPMF 
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
      FPMF = TAILFACTOR_FLX/FR(1)

      DO IJ=IJS,IJL
        IF (CICVR(IJ).LE.CITHRSH_TAIL) THEN
          FM2 = MAX(FMWS(IJ),FM(IJ))*FPMH
          FPM = FPPM/MAX(USNEW(IJ),EPSMIN)
          FPM4 = MAX(FM2,FPM)
          MIJ(IJ) = NINT(LOG10(FPM4)*FLOGSPRDM1)+1
          MIJ(IJ) = MIN(MAX(1,MIJ(IJ)),NFRE)
!!! 
          FM2 = MAX(FMWS(IJ),FM(IJ))*FPMF
          FPM4 = MAX(FM2,FPM)
          MIJFLX(IJ) = NINT(LOG10(FPM4)*FLOGSPRDM1)+1
          MIJFLX(IJ) = MIN(MAX(1,MIJFLX(IJ)),NFRE)
        ELSE
          MIJ(IJ) = NFRE
          MIJFLX(IJ) = NFRE
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
