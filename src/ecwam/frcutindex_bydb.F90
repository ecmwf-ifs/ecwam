      SUBROUTINE FRCUTINDEX_BYDB (KIJS, KIJL, FM, UFRIC, CICOVER,     &
            &                       MIJ)

! ----------------------------------------------------------------------

!**** *FRCUTINDEX_BYDB* - RETURNS THE LAST FREQUENCY INDEX OF
!                        PROGNOSTIC PART OF SPECTRUM.

!**   INTERFACE.
!     ----------

!       *CALL* *FRCUTINDEX_BYDB (KIJS, KIJL, FM, UFRIC, CICOVER,MIJ)
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *FM*     - MEAN FREQUENCY
!          *UFRIC*  - FRICTION VELOCITY IN M/S
!          *CICOVER*- CICOVER 
!          *MIJ*    - LAST FREQUENCY INDEX for imposing high frequency tail



!     METHOD.
!     -------

!*    COMPUTES LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM 
!     ACCORDING TO BYDB   

!     EXTERNALS.
!     ---------

!     REFERENCE.
!     ----------

!     ORIGIN.
!     ----------
!     Adapted from Babanin Young Donelan & Banner (BYDB) physics 
!     as implemented as ST6 in WAVEWATCH-III
!     WW3 module:       W3SRCEMD
!     WW3 subroutine:   
!     Implementation into ECWAM DECEMBER 2021 by J. Kousal 

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : TAILFACTOR, TAILFACTOR_PM
      USE YOWFRED  , ONLY : FR       ,DFIM       ,FRATIO   ,FLOGSPRDM1, &
     &                DELTH          ,RHOWG_DFIM ,FRIC
      USE YOWICE   , ONLY : CITHRSH_TAIL
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,EPSMIN, EPSUS
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM), INTENT(OUT) :: MIJ(KIJL)

      REAL(KIND=JWRB),DIMENSION(KIJL), INTENT(IN) :: FM, UFRIC

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      INTEGER(KIND=JWIM)          :: IJ, NK, NKH, NKH1, M
      REAL(KIND=JWRB), PARAMETER  :: SIN6FC = 6.0_JWRB
      REAL(KIND=JWRB)             :: FXFM, FXPM, FACTI1, FACTI2 ! constants
      REAL(KIND=JWRB)             :: FHIGH ! Cut-off frequency in integration (rad/s)
      REAL(KIND=JWRB)             :: SIGNK ! LAST FREQUENCY [RAD]
      REAL(KIND=JWRB)             :: USTM1


! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('FRCUTINDEX_BYDB',0,ZHOOK_HANDLE)

      NK     = NFRE
      FXFM   = SIN6FC
      FXFM   = FXFM * ZPI
      FXPM   = 4.0_JWRB             !TODO: 4.0_JWRB is the factor for the tail (is this right)
      FXPM   = FXPM * G / 28.0_JWRB !TODO: should this be FRIC?
      SIGNK  = ZPI*FR(NFRE)

      DO IJ=KIJS,KIJL
        IF (CICOVER(IJ) <= CITHRSH_TAIL) THEN            

            USTM1     = 1.0_JWRB/MAX(UFRIC(IJ),EPSUS)                       ! Protect the code

            IF (FXFM .LE. 0) THEN
            FHIGH = SIGNK ! LAST FREQ i.e. let tail evolve freely
            ELSE
            FHIGH = MAX (FXFM * FM(IJ), FXPM * USTM1 )
            ENDIF


            FACTI1 = 1.0_JWRB / LOG(FRATIO)
            FACTI2 = 1.0_JWRB - LOG(ZPI*FR(1)) * FACTI1

            NKH    = MIN ( NK , INT(FACTI2+FACTI1*LOG(MAX(1.0E-7_JWRB,FHIGH))) )
            NKH1   = MIN ( NK , NKH+1 )


            IF (FXFM .LE. 0) THEN
                  FHIGH = SIGNK
            ELSE
                  FHIGH = MIN ( SIGNK, MAX(FXFM * FM(IJ), FXPM * USTM1) )
            ENDIF
            NKH    = MAX ( 2 , MIN ( NKH1 ,                           &
                  INT ( FACTI2 + FACTI1*LOG(MAX(1.0E-7_JWRB,FHIGH)) ) ) )

            MIJ(IJ) = NKH
        ELSE
            MIJ(IJ) = NFRE
        ENDIF
      END DO

      IF (LHOOK) CALL DR_HOOK('FRCUTINDEX_BYDB',1,ZHOOK_HANDLE)

      END SUBROUTINE FRCUTINDEX_BYDB
