      SUBROUTINE OUTSETWMASK (IJS, IJL, IODP, CICVR, BOUT)
! ----------------------------------------------------------------------

!**** *OUTSETWMASK* -


!*    PURPOSE.
!     --------

!    IMPOSE THE ICE MASK AND LAND MASK ON OUTPUT PARAMETERS


!**   INTERFACE.
!     ----------

!        *CALL* *OUTSETWMASK (IJS, IJL, IODP, CICVR, BOUT)
!         *IJS*    - INDEX OF FIRST LOCAL GRIDPOINT.
!         *IJL*    - INDEX OF LAST LOCAL GRIDPOINT.
!         *CICVR*  - SEA ICE COVER FIELD.
!         *IODP*   - LAND MASK IF IODP(IJ)=0
!         *BOUT*   - OUTPUT PARAMETERS

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : JPPFLAG  ,IPRMINFO, NIPRMOUT, ITOBOUT
      USE YOWICE   , ONLY : LICERUN  ,CITHRSH
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSTAT  , ONLY : LLSOURCE
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL), INTENT(IN) :: IODP
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: CICVR 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NIPRMOUT), INTENT(INOUT) :: BOUT

      INTEGER(KIND=JWIM) :: IJ, ITG, IR

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('OUTSETWMASK',0,ZHOOK_HANDLE)

      DO IR=1,JPPFLAG
        ITG=ITOBOUT(IR)
        IF(ITG.GT.0) THEN
          IF (LICERUN .AND. LLSOURCE .AND. IPRMINFO(IR,4).EQ.1) THEN
!         SEA ICE MASK IS APPLIED
            DO IJ = IJS,IJL
              IF (CICVR(IJ).GT.CITHRSH) BOUT(IJ,ITG) = ZMISS
            ENDDO
          ENDIF

          IF (IPRMINFO(IR,5).EQ.1) THEN
!           SEA MASK IS APPLIED
            DO IJ = IJS,IJL
              BOUT(IJ,ITG) = BOUT(IJ,ITG)*IODP(IJ) + (1-IODP(IJ))*ZMISS
            ENDDO
          ENDIF

        ENDIF
      ENDDO

      IF (LHOOK) CALL DR_HOOK('OUTSETWMASK',1,ZHOOK_HANDLE)

      END SUBROUTINE OUTSETWMASK
