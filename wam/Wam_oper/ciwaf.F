      SUBROUTINE CIWAF (IJS,IJL,CICOVER,CITHICK,CIWA)

! ----------------------------------------------------------------------

!**** *CIWAF* - COMPUTE SEA ICE WAVE ATTENUATION FACTORS DUE TO ICE FLOES
!               SCATTERING.

!*    PURPOSE.
!     --------

!       CIWAF COMPUTES SEA ICE WAVE ATTENUATION FACTORS

!**   INTERFACE.
!     ----------

!       *CALL* *CIWAF (IJS,IJL,CICOVER,CITHICK,CIWA)

!
!          *IJS*     - INDEX OF FIRST POINT.
!          *IJL*     - INDEX OF LAST POINT.
!          *CICOVER*  -SEA ICE COVER.
!          *CITHICK*  -SEA ICE THICKNESS. 
!          *CIWA*     -SEA ICE WAVE ATTENUATION FACTOR. 

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCES.                                                       
!     -----------  

!     KOHOUT AND MEYLAN, 2008: JGR, 113, doi:10.1029/2007JC004424


! ----------------------------------------------------------------------

      USE YOWFRED  , ONLY : FR       ,GOM
      USE YOWICE   , ONLY : NICT   ,NICH     ,TICMIN   ,HICMIN   ,
     &              DTIC   ,DHIC   ,CIDEAC   ,CIBLOCK
      USE YOWPARAM , ONLY : NFRE
      USE YOWPCONS , ONLY : EPSMIN
      USE YOWSHAL  , ONLY : TCGOND   ,INDEP 
      USE YOWSTAT  , ONLY : IDELT   ,ISHALLO 
      USE YOWTEST  , ONLY : IU06    ,ITEST

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL 
      REAL,DIMENSION(IJS:IJL), INTENT(IN) :: CICOVER
      REAL,DIMENSION(IJS:IJL), INTENT(INOUT) :: CITHICK
      REAL,DIMENSION(IJS:IJL,NFRE), INTENT(OUT) :: CIWA

      INTEGER :: ICM, I, MAXICM
      INTEGER :: M, IJ
      INTEGER :: IT, IT1, IH, IH1 

      REAL :: CIDMIN, CIDMAX, CIDMEAN, CIFRGL, CIFRGMT
      REAL :: SD, SN 
      REAL :: TW, X
      REAL :: A, B, C
      REAL :: CIDEAC_INT, WT, WT1, WH, WH1 
      REAL :: ZHOOK_HANDLE
      REAL,DIMENSION(IJS:IJL) :: DINV
      REAL,DIMENSION(IJS:IJL,NFRE) :: ALP 

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('CIWAF',0,ZHOOK_HANDLE)
#endif

!     following  Dumont et al. (2011), eqn (13):
      ! sea ice fragility
      CIFRGL=0.955
      ! sea ice foes minimum size (m)
      CIDMIN=20.
      ! sea ice fragmentation number (>2)
      CIFRGMT=2

      A=200.
      C=300.
!     the maximum number of fragmentation step is set
!     for the limit value for sea ice = 0 (CIDMAX=A)
!     to insure that CIDMEAN will be increasing with increasing sea ice cover
      MAXICM=INT(LOG(A/CIDMIN)/LOG(CIFRGMT))
      DO IJ=IJS,IJL
        IF(CITHICK(IJ).GT.0.0) THEN
          ! sea ice foes maxmimum size (m)
!!!       testing making it function of sea ice cover
          CIDMAX=A+C*CICOVER(IJ)
          ICM=MIN(INT(LOG(CIDMAX/CIDMIN)/LOG(CIFRGMT)),MAXICM)
          SN=0.0
          SD=0.0
           DO I=0,ICM
            X=(CIFRGMT**2*CIFRGL)**I
            SN=SN+X*CIDMAX/CIFRGMT**I
            SD=SD+X
          ENDDO
        CIDMEAN=SN/SD
        DINV(IJ)=1./CIDMEAN
        ELSE
          DINV(IJ)=CIDMIN
        ENDIF
      ENDDO


      DO M=1,NFRE
        TW=1./FR(M)
        IT=FLOOR((TW-TICMIN)/DTIC+1)
        IT=MAX(1,MIN(IT,NICT))
        IT1=IT+1
        IT1=MAX(1,MIN(IT1,NICT))
        WT1=MAX(MIN(1.,(TW-(TICMIN+(IT-1)*DTIC))/DTIC),0.)
        WT=1.-WT1
        DO IJ=IJS,IJL
          IF(CITHICK(IJ).GT.0.0) THEN
            IH=FLOOR((CITHICK(IJ)-HICMIN)/DHIC+1)
            IH=MAX(1,MIN(IH,NICH))
            IH1=IH+1
            IH1=MAX(1,MIN(IH1,NICH))
            WH1=MAX(MIN(1.,(CITHICK(IJ)-(HICMIN+(IH-1)*DHIC))/DHIC),0.)
            WH=1.-WH1
            CIDEAC_INT=WT*(WH*CIDEAC(IT,IH)+ WH1*CIDEAC(IT,IH1)) +
     &                WT1*(WH*CIDEAC(IT1,IH)+WH1*CIDEAC(IT1,IH1))
!!!            ALP(IJ,M)=CICOVER(IJ)*CIDEAC_INT*DINV(IJ)
            ALP(IJ,M)=CICOVER(IJ)*EXP(CIDEAC_INT)*DINV(IJ)
          ELSE
            ALP(IJ,M)=0.0
          ENDIF
        ENDDO
      ENDDO

      IF (ISHALLO.NE.1) THEN
        DO M=1,NFRE
          DO IJ=IJS,IJL
            X=ALP(IJ,M)*TCGOND(INDEP(IJ),M)*IDELT
            IF(X.LT.EPSMIN) THEN
              CIWA(IJ,M)=1.0
            ELSE IF(CICOVER(IJ).GT.CIBLOCK) THEN
              CIWA(IJ,M)=0.0
            ELSE
              CIWA(IJ,M)=EXP(-MIN(X,50.))
            ENDIF
          ENDDO
        ENDDO
      ELSE
        DO M=1,NFRE
          DO IJ=IJS,IJL
            X=ALP(IJ,M)*GOM(M)*IDELT
            IF(X.LT.EPSMIN) THEN
              CIWA(IJ,M)=1.0
            ELSE IF(CICOVER(IJ).GT.CIBLOCK) THEN
              CIWA(IJ,M)=0.0
            ELSE
              CIWA(IJ,M)=EXP(-MIN(X,50.))
            ENDIF
          ENDDO
        ENDDO
      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('CIWAF',1,ZHOOK_HANDLE)
#endif
      RETURN
      END SUBROUTINE CIWAF
