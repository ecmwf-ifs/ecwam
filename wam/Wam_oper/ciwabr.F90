      SUBROUTINE CIWABR (IJS, IJL, KIJS, KIJL, CICOVER, GFL, CIWAB)

! ----------------------------------------------------------------------

!**** *CIWABR* - COMPUTE SEA ICE WAVE ATTENUATION FACTORS DUE TO ICE FLOES
!                BOTTOM FRICTION.

!*    PURPOSE.
!     --------

!       CIWABR COMPUTES SEA ICE WAVE ATTENUATION FACTORS DUE TO ICE FLOES
!              BOTTOM FRICTION.

!**   INTERFACE.
!     ----------

!       *CALL* *CIWABR (IJS, IJL, KIJS,KIJL,CICOVER,GFL,CIWAB)

!
!          *IJS:IJL   - 1st DIMENSION OF GFL
!          *KIJS*     - INDEX OF FIRST POINT.
!          *KIJL*     - INDEX OF LAST POINT.
!          *CICOVER*  -SEA ICE COVER.
!          *GFL*      -ENERGY SPECTRUM. 
!          *CIWAB*    -SEA ICE WAVE ATTENUATION FACTOR DUE TO ICE FLOE BOTTOM FRICTION 

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCES.                                                       
!     -----------  

!     KOHOUT A., M. MEYLAN, D PLEW, 2011: ANNALS OF GLACIOLOGY, 2011. 


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,DFIM     , DELTH     ,GOM
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE , CDICWA
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,ZPI4GM2    ,EPSMIN
      USE YOWSHAL  , ONLY : TCGOND  ,TFAK      ,INDEP
      USE YOWSTAT  , ONLY : IDELT   ,ISHALLO 
      USE YOWTEST  , ONLY : IU06    ,ITEST

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, KIJS, KIJL 
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(IN) :: CICOVER
      REAL(KIND=JWRB),DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: GFL
      REAL(KIND=JWRB),DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(OUT) :: CIWAB

      INTEGER(KIND=JWIM) :: K, M, IJ
      REAL(KIND=JWRB) :: EWH 
      REAL(KIND=JWRB) :: X, ALP
      REAL(KIND=JWRB),DIMENSION(NFRE) :: XK2 
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('CIWABR',0,ZHOOK_HANDLE)

      IF( .NOT. LICERUN .OR. LMASKICE ) THEN

        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=KIJS,KIJL
              CIWAB(IJ,K,M)=1.0_JWRB
            ENDDO
          ENDDO
        ENDDO

      ELSE

        IF (ISHALLO.NE.1) THEN
          DO M=1,NFRE
            DO K=1,NANG
              DO IJ=KIJS,KIJL
                EWH=4.0_JWRB*SQRT(MAX(EPSMIN,GFL(IJ,K,M)*DFIM(M)))
                XK2(M)=TFAK(INDEP(IJ),M)**2
                ALP=CDICWA*XK2(M)*EWH
                X=ALP*TCGOND(INDEP(IJ),M)*IDELT
                CIWAB(IJ,K,M)=1.0_JWRB-CICOVER(IJ)*(1.0_JWRB-EXP(-MIN(X,50.0_JWRB)))
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO M=1,NFRE
            XK2(M)=ZPI4GM2*FR(M)**4
          ENDDO
          DO M=1,NFRE
            DO K=1,NANG
              DO IJ=KIJS,KIJL
                EWH=4.0_JWRB*SQRT(MAX(EPSMIN,GFL(IJ,K,M)*DFIM(M)))
                ALP=CDICWA*XK2(M)*EWH
                X=ALP*GOM(M)*IDELT
                CIWAB(IJ,K,M)=1.0_JWRB-CICOVER(IJ)*(1.0_JWRB-EXP(-MIN(X,50.0_JWRB)))
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      IF (LHOOK) CALL DR_HOOK('CIWABR',1,ZHOOK_HANDLE)

      END SUBROUTINE CIWABR
