      SUBROUTINE CIWABR (IJS,IJL,CICOVER,FL,CIWAB)

! ----------------------------------------------------------------------

!**** *CIWABR* - COMPUTE SEA ICE WAVE ATTENUATION FACTORS DUE TO ICE FLOES
!                BOTTOM FRICTION.

!*    PURPOSE.
!     --------

!       CIWABR COMPUTES SEA ICE WAVE ATTENUATION FACTORS DUE TO ICE FLOES
!              BOTTOM FRICTION.

!**   INTERFACE.
!     ----------

!       *CALL* *CIWABR (IJS,IJL,CICOVER,FL,CIWAB)

!
!          *IJS*     - INDEX OF FIRST POINT.
!          *IJL*     - INDEX OF LAST POINT.
!          *CICOVER*  -SEA ICE COVER.
!          *FL*       -ENERGY SPECTRUM. 
!          *CIWAB*    -SEA ICE WAVE ATTENUATION FACTOR DUE TO ICE FLOE BOTTOM FRICTION 

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCES.                                                       
!     -----------  

!     KOHOUT A., M. MEYLAN, D PLEW, 2011: ANNALS OF GLACIOLOGY, 2011. 


! ----------------------------------------------------------------------

      USE YOWFRED  , ONLY : FR       ,DFIM     , DELTH     ,GOM
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE , CDICWA
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,EPSMIN
      USE YOWSHAL  , ONLY : TCGOND  ,TFAK      ,INDEP
      USE YOWSTAT  , ONLY : IDELT   ,ISHALLO 
      USE YOWTEST  , ONLY : IU06    ,ITEST

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: K, M, IJ, IJS, IJL 
      REAL :: ZPI4G2 
      REAL :: EWH 
      REAL :: X, ALP
      REAL,DIMENSION(NFRE) :: XK2 
      REAL,DIMENSION(IJS:IJL) :: CICOVER
      REAL,DIMENSION(IJS:IJL,NANG,NFRE) :: CIWAB
      REAL,DIMENSION(IJS:IJL,NANG,NFRE) :: FL
      REAL :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('CIWABR',0,ZHOOK_HANDLE)

      IF( .NOT. LICERUN .OR. LMASKICE ) THEN

        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=IJS,IJL
              CIWAB(IJ,K,M)=1.
            ENDDO
          ENDDO
        ENDDO

      ELSE

        IF (ISHALLO.NE.1) THEN
          DO M=1,NFRE
            DO K=1,NANG
              DO IJ=IJS,IJL
                EWH=4.*SQRT(MAX(EPSMIN,FL(IJ,K,M)*DFIM(M)))
                XK2(M)=TFAK(INDEP(IJ),M)**2
                ALP=CDICWA*XK2(M)*EWH
                X=ALP*TCGOND(INDEP(IJ),M)*IDELT
                CIWAB(IJ,K,M)=1.-CICOVER(IJ)*(1.-EXP(-MIN(X,50.)))
              ENDDO
            ENDDO
          ENDDO
        ELSE
          ZPI4G2=ZPI**4/G**2
          DO M=1,NFRE
            XK2(M)=ZPI4G2*FR(M)**4
          ENDDO
          DO M=1,NFRE
            DO K=1,NANG
              DO IJ=IJS,IJL
                EWH=4.*SQRT(MAX(EPSMIN,FL(IJ,K,M)*DFIM(M)))
                ALP=CDICWA*XK2(M)*EWH
                X=ALP*GOM(M)*IDELT
                CIWAB(IJ,K,M)=1.-CICOVER(IJ)*(1.-EXP(-MIN(X,50.)))
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      IF (LHOOK) CALL DR_HOOK('CIWABR',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE CIWABR
