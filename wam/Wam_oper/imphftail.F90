      SUBROUTINE IMPHFTAIL (IJS, IJL, KIJS, KIJL, MIJ, FLM, GFL) 
! ----------------------------------------------------------------------

!**** *IMPHFTAIL* - IMPOSE A HIGH FREQUENCY TAIL TO THE SPECTRUM


!*    PURPOSE.
!     --------

!     IMPOSE A HIGH FREQUENCY TAIL TO THE SPECTRUM ABOVE FREQUENCY INDEX MIJ


!**   INTERFACE.
!     ----------

!       *CALL* *IMPHFTAIL (IJS, IJL, KIJS, KIJL, MIJ, FLM, GFL)
!          *IJS:IJL* - 1st DIMENSION OF GFL
!          *KIJS*    - INDEX OF FIRST GRIDPOINT
!          *KIJL*    - INDEX OF LAST GRIDPOINT
!          *MIJ*     - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!          *FLM*     - SPECTAL DENSITY MINIMUM VALUE
!          *GFL*     - SPECTRUM (INPUT AND OUTPUT).

!     METHOD.
!     -------

!     EXTERNALS.
!     ---------

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FRM5 
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSHAL  , ONLY : TCGOND   ,TFAK     ,INDEP
      USE YOWSTAT  , ONLY : ISHALLO

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, KIJS, KIJL
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(IN) :: MIJ
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG), INTENT(IN) :: FLM 

      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: GFL


      INTEGER(KIND=JWIM) :: IJ, K, M

      REAL(KIND=JWRB) :: AKM1, TFAC
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE) :: TEMP2

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('IMPHFTAIL',0,ZHOOK_HANDLE)

!*    DIAGNOSTIC TAIL.
!     ----------------

      IF(ISHALLO.EQ.1) THEN
        DO IJ=KIJS,KIJL
          DO M=MIJ(IJ),NFRE
            TEMP2(IJ,M) = FRM5(M)
          ENDDO
        ENDDO
      ELSE
        DO IJ=KIJS,KIJL
          DO M=MIJ(IJ),NFRE
            AKM1 = 1._JWRB/TFAK(INDEP(IJ),M)
            TEMP2(IJ,M) = AKM1**3/TCGOND(INDEP(IJ),M)
          ENDDO
        ENDDO
      ENDIF

      DO IJ=KIJS,KIJL
        DO M=MIJ(IJ)+1,NFRE
          TEMP2(IJ,M) = TEMP2(IJ,M)/TEMP2(IJ,MIJ(IJ))
        ENDDO
      ENDDO

!*    MERGE TAIL INTO SPECTRA.
!     ------------------------
      DO K=1,NANG
        DO IJ=KIJS,KIJL
          TFAC = GFL(IJ,K,MIJ(IJ))
          DO M=MIJ(IJ)+1,NFRE
            GFL(IJ,K,M) = MAX(TEMP2(IJ,M)*TFAC,FLM(IJ,K))
          ENDDO
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('IMPHFTAIL',1,ZHOOK_HANDLE)

      END SUBROUTINE IMPHFTAIL
