      SUBROUTINE FRCUTINDEX (IJS, IJL, FM, FMWS, USNEW, CICVR,
     &                       MIJ, RHOWGDFTH)

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
!          *MIJ*    - LAST FREQUENCY INDEX
!          *RHOWGDFTH - WATER DENSITY * G * DF * DTHETA
!                       FOR TRAPEZOIDAL INTEGRATION BETWEEN FR(1) and FR(MIJ) 
!                       !!!!!!!!  RHOWGDFTH=0 FOR FR > FR(MIJ)


!     METHOD.
!     -------

!*    COMPUTES LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!*    FREQUENCIES LE 2.5*MAX(FMWS,FM).

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

      USE YOWCOUP  , ONLY : TAILFACTOR, TAILFACTOR_PM
      USE YOWFRED  , ONLY : FR       ,DFIM       ,FRATIO   ,FLOGSPRDM1,
     &                DELTH          ,DFIM_END_U ,FRIC
      USE YOWICE   , ONLY : CITHRSH_TAIL
      USE YOWPARAM , ONLY : NFRE
      USE YOWPCONS , ONLY : G        ,ROWATER    ,ZPI      ,EPSMIN
      USE YOWSTAT  , ONLY : IPHYS
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL
      INTEGER, INTENT(OUT) :: MIJ(IJS:IJL)

      REAL,DIMENSION(IJS:IJL), INTENT(IN) :: FM, FMWS, USNEW, CICVR
      REAL,DIMENSION(IJS:IJL,NFRE), INTENT(OUT) :: RHOWGDFTH 

      INTEGER :: IJ, M

      REAL :: FPMH, FPPM, FM2, FPM, FPM4
      REAL :: XLOGDFRTH
      REAL :: ROG_XLOGDFRTH, HROG_XLOGDFRTH
      REAL :: ZHOOK_HANDLE

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('FRCUTINDEX',0,ZHOOK_HANDLE)
#endif
! ----------------------------------------------------------------------

!*    COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!*    FREQUENCIES LE MAX(TAILFACTOR*MAX(FMNWS,FM),TAILFACTOR_PM*FPM),
!*    WHERE FPM IS THE PIERSON-MOSKOWITZ FREQUENCY BASED ON FRICTION
!*    VELOCITY. (FPM=G/(FRIC*ZPI*USTAR))
!     ------------------------------------------------------------

      XLOGDFRTH=LOG(FRATIO)*DELTH
      ROG_XLOGDFRTH = ROWATER*G*XLOGDFRTH
      HROG_XLOGDFRTH = 0.5*ROG_XLOGDFRTH
      FPMH = TAILFACTOR/FR(1)
      FPPM = TAILFACTOR_PM*G/(FRIC*ZPI*FR(1))

      IF(IPHYS.EQ.0) THEN
        DO IJ=IJS,IJL
          IF (CICVR(IJ).LE.CITHRSH_TAIL) THEN
            FM2 = MAX(FMWS(IJ),FM(IJ))*FPMH
            FPM = FPPM/MAX(USNEW(IJ),EPSMIN)
            FPM4 = MAX(FM2,FPM)
            MIJ(IJ) = INT(LOG10(FPM4)*FLOGSPRDM1)+1
            MIJ(IJ) = MIN(MAX(1,MIJ(IJ)),NFRE)
          ELSE
            MIJ(IJ) = NFRE
          ENDIF
        ENDDO
      ELSE
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
      ENDIF

!     COMPUTE RHOWGDFTH
      DO IJ=IJS,IJL
        RHOWGDFTH(IJ,1) = HROG_XLOGDFRTH*FR(1)
        DO M=2,MIJ(IJ)-1
          RHOWGDFTH(IJ,M) = ROG_XLOGDFRTH*FR(M)
        ENDDO
        RHOWGDFTH(IJ,MIJ(IJ)) = HROG_XLOGDFRTH*FR(MIJ(IJ))
        DO M=MIJ(IJ)+1,NFRE
          RHOWGDFTH(IJ,M) = 0.
        ENDDO
      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('FRCUTINDEX',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE FRCUTINDEX
