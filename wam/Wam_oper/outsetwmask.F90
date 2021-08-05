      SUBROUTINE OUTSETWMASK (KIJS, KIJL, IODP, CICVR, BOUT)
! ----------------------------------------------------------------------

!**** *OUTSETWMASK* -


!*    PURPOSE.
!     --------

!    IMPOSE THE ICE MASK AND LAND MASK ON OUTPUT PARAMETERS


!**   INTERFACE.
!     ----------

!        *CALL* *OUTSETWMASK (KIJS, KIJL, IODP, CICVR, BOUT)
!         *KIJS*   - INDEX OF FIRST LOCAL GRIDPOINT.
!         *KIJL*   - INDEX OF LAST LOCAL GRIDPOINT.
!         *CICVR   - SEA ICE COVER FIELD.
!         *IODP*   - LAND MASK IF IODP(IJ)=0
!         *BOUT*   - OUTPUT PARAMETERS

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : JPPFLAG  ,IPRMINFO, NIPRMOUT, ITOBOUT
      USE YOWICE   , ONLY : LICERUN  ,CITHRSH
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSTAT  , ONLY : LLSOURCE

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(IN) :: IODP
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: CICVR 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NIPRMOUT), INTENT(INOUT) :: BOUT

      INTEGER(KIND=JWIM) :: IJ, ITG, IR

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

!----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('OUTSETWMASK',0,ZHOOK_HANDLE)

      DO IR=1,JPPFLAG
        ITG=ITOBOUT(IR)
        IF (ITG > 0) THEN
          IF (LICERUN .AND. LLSOURCE .AND. IPRMINFO(IR,4) == 1) THEN
!         SEA ICE MASK IS APPLIED
            DO IJ = KIJS,KIJL
              IF (CICVR(IJ).GT.CITHRSH) BOUT(IJ,ITG) = ZMISS
            ENDDO
          ENDIF

          IF (IPRMINFO(IR,5) == 1) THEN
!           SEA MASK IS APPLIED
            DO IJ = KIJS,KIJL
              BOUT(IJ,ITG) = BOUT(IJ,ITG)*IODP(IJ) + (1-IODP(IJ))*ZMISS
            ENDDO
          ENDIF

        ENDIF
      ENDDO

      IF (LHOOK) CALL DR_HOOK('OUTSETWMASK',1,ZHOOK_HANDLE)

      END SUBROUTINE OUTSETWMASK
