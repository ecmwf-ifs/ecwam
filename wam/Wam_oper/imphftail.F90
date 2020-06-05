      SUBROUTINE IMPHFTAIL (IJS, IJL, MIJ, FLM, FL1) 
! ----------------------------------------------------------------------

!**** *IMPHFTAIL* - IMPOSE A HIGH FREQUENCY TAIL TO THE SPECTRUM


!*    PURPOSE.
!     --------

!     IMPOSE A HIGH FREQUENCY TAIL TO THE SPECTRUM ABOVE FREQUENCY INDEX MIJ


!**   INTERFACE.
!     ----------

!       *CALL* *IMPHFTAIL (IJS, IJL, MIJ, FLM, FL1)
!          *IJS*    - INDEX OF FIRST GRIDPOINT
!          *IJL*    - INDEX OF LAST GRIDPOINT
!          *MIJ*    - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!          *FLM*    - SPECTAL DENSITY MINIMUM VALUE
!          *FL1*    - FREQUENCY SPECTRUM(INPUT AND OUTPUT).

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

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM), INTENT(IN) :: MIJ(IJS:IJL)
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG), INTENT(IN) :: FLM 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(INOUT) :: FL1

      INTEGER(KIND=JWIM) :: IJ, K, M

      REAL(KIND=JWRB) :: AKM1, TFAC
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NFRE) :: TEMP2

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('IMPHFTAIL',0,ZHOOK_HANDLE)

!*    DIAGNOSTIC TAIL.
!     ----------------

      IF(ISHALLO.EQ.1) THEN
        DO IJ=IJS,IJL
          DO M=MIJ(IJ),NFRE
            TEMP2(IJ,M) = FRM5(M)
          ENDDO
        ENDDO
      ELSE
        DO IJ=IJS,IJL
          DO M=MIJ(IJ),NFRE
            AKM1 = 1._JWRB/TFAK(INDEP(IJ),M)
            TEMP2(IJ,M) = AKM1**3/TCGOND(INDEP(IJ),M)
          ENDDO
        ENDDO
      ENDIF

      DO IJ=IJS,IJL
        DO M=MIJ(IJ)+1,NFRE
          TEMP2(IJ,M) = TEMP2(IJ,M)/TEMP2(IJ,MIJ(IJ))
        ENDDO
      ENDDO

!*    MERGE TAIL INTO SPECTRA.
!     ------------------------
      DO K=1,NANG
        DO IJ=IJS,IJL
          TFAC = FL1(IJ,K,MIJ(IJ))
          DO M=MIJ(IJ)+1,NFRE
            FL1(IJ,K,M) = MAX(TEMP2(IJ,M)*TFAC,FLM(IJ,K))
          ENDDO
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('IMPHFTAIL',1,ZHOOK_HANDLE)

      END SUBROUTINE IMPHFTAIL
