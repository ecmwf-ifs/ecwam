      SUBROUTINE IMPLSCH (FL3, FL, IJS, IJL, IG,
     &                    THWOLD, USOLD,
     &                    TAUW, Z0OLD,
     &                    ROAIRO, ZIDLOLD,
     &                    CICVR, CIWA,
     &                    U10NEW, THWNEW, USNEW,
     &                    Z0NEW, ROAIRN, ZIDLNEW,
     &                    SL, FCONST)

! ----------------------------------------------------------------------

!**** *IMPLSCH* - IMPLICIT SCHEME FOR TIME INTEGRATION OF SOURCE
!****             FUNCTIONS.

!     S.D.HASSELMANN.  MPI
!     H. GUENTHER AND L. ZAMBRESKY  OPTIMIZATION PERFORMED.
!     H. GUENTHER      GKSS/ECMWF   OCTOBER 1989  NEW WIND FIELD
!                                                 INTERFACE AND
!                                                 TIME COUNTING
!     P.A.E.M. JANSSEN KNMI         AUGUST  1990  COUPLED MODEL
!     H. GUENTHER      GKSS/ECMWF   JUNE    1991  NEW SEPARATION OF
!                                                  DIAG- AND PROGNOSTIC
!                                                  PART OF SPECTRUM.
!     P.A.E.M. JANSSEN ECMWF        FEBRUARY 1995  ADD MINIMUM VALUE
!                                                  (FLMIN).
!     J. BIDLOT ECMWF               FEBRUARY 1996 MESSAGE PASSING
!     J. BIDLOT ECMWF               FEBRUARY 1997 MESSAGE PASSING
!     J. BIDLOT ECMWF               FEBRUARY 2001 MODIFY CALLING ORDER
!                                   BETWEEN SINPUT-STRESSO-AIRSEA
!     S. ABDALLA       ECMWF        OCTOBER 2001
!                                   INCLUSION OF AIR DENSITY AND Zi/L.

!*    PURPOSE.
!     --------

!       THE IMPLICIT SCHEME ENABLES THE USE OF A TIMESTEP WHICH IS
!       LARGE COMPARED WITH THE CHARACTERISTIC DYNAMIC TIME SCALE.
!       THE SCHEME IS REQUIRED FOR THE HIGH FREQUENCIES WHICH
!       RAPIDLY ADJUST TO A QUASI-EQUILIBRIUM.

!**   INTERFACE.
!     ----------

!       *CALL* *IMPLSCH (FL3, FL, IJS, IJL, IG,
!    1                    THWOLD,USOLD,TAUW,Z0OLD,
!    &                    ROAIRO, ZIDLOLD, 
!    &                    CICVR, CIWA,
!    2                    U10NEW,THWNEW,USNEW,Z0NEW,ROAIRN,ZIDLNEW,
!    3                    SL,FCONST)
!          *FL3*    - FREQUENCY SPECTRUM(INPUT AND OUTPUT).
!          *FL*     - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *IJS*    - INDEX OF FIRST GRIDPOINT
!          *IJL*    - INDEX OF LAST GRIDPOINT
!          *IG*     - BLOCK NUMBER
!      *U10NEW*    NEW WIND SPEED IN M/S.
!      *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!      *THWOLD*    INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
!                  WIND VELOCITY.
!      *USNEW*     NEW FRICTION VELOCITY IN M/S.
!      *USOLD*     INTERMEDIATE STORAGE OF MODULUS OF FRICTION
!                  VELOCITY.
!      *Z0NEW*     ROUGHNESS LENGTH IN M.
!      *Z0OLD*     INTERMEDIATE STORAGE OF ROUGHNESS LENGTH IN
!                  M.
!      *TAUW*      WAVE STRESS IN (M/S)**2
!      *ROAIRN*    AIR DENSITY IN KG/M3.
!      *ROAIRO*    INTERMEDIATE STORAGE OF AIR DENSITY.
!      *ZIDLNEW*   Zi/L (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!      *ZIDLOLD*   INTERMEDIATE STORAGE OF Zi/L.
!      *CICVR*     SEA ICE COVER.
!      *CIWA*      SEA ICE WAVE ATTENUATION.
!      *SL*        REAL      TOTAL SOURCE FUNCTION ARRAY.
!      *FCONST*    REAL      = 1 FOR PROGNOSTIC FREQUENCY BANDS.
!                            = 0 FOR DIAGNOSTIC FREQUENCY BANDS.

!     METHOD.
!     -------

!       THE SPECTRUM AT TIME (TN+1) IS COMPUTED AS
!       FN+1=FN+DELT*(SN+SN+1)/2., WHERE SN IS THE TOTAL SOURCE
!       FUNCTION AT TIME TN, SN+1=SN+(DS/DF)*DF - ONLY THE DIAGONAL
!       TERMS OF THE FUNCTIONAL MATRIX DS/DF ARE YOWPUTED, THE
!       NONDIAGONAL TERMS ARE NEGLIGIBLE.
!       THE ROUTINE IS CALLED AFTER PROPAGATION FOR TIME PERIOD
!       BETWEEN TWO PROPAGATION CALLS - ARRAY FL3 CONTAINS THE
!       SPECTRUM AND FL IS USED AS AN INTERMEDIATE STORAGE FOR THE
!       DIAGONAL TERM OF THE FUNCTIONAL MATRIX.

!     EXTERNALS.
!     ---------

!       *AIRSEA*    - SURFACE LAYER STRESS AND ROUGHNESS LENGTH.
!       *CREWFN*    - CREATE A WIND FILE NAME.
!       *FEMEAN*    - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT.
!       *INCDATE*   - UPDATE DATE TIME GROUP.
!SHALLOW
!       *SBOTTOM*   - COMPUTES BOTTOM DISSIPATION SOURCE TERM AND
!                     LINEAR CONTRIBUTION TO FUNCTIONAL MATRIX.
!SHALLOW
!       *SDISSIP*   - COMPUTATION OF DISSIPATION SOURCE FUNCTION
!                     AND LINEAR CONTRIBUTION OF DISSIPATION TO
!                     FUNCTIONAL MATRIX IN IMPLICIT SCHEME.
!       *SEMEAN*    - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.
!       *SINPUT*    - COMPUTATION OF INPUT SOURCE FUNCTION, AND
!                     LINEAR CONTRIBUTION OF INPUT SOURCE FUNCTION
!                     TO FUNCTIONAL MATRIX IN IMPLICIT SCHEME.
!       *SNONLIN*   - COMPUTATION OF NONLINEAR TRANSFER RATE AND
!                     DIAGONAL LINEAR CONTRIBUTION OF NONLINEAR SOURCE
!                     FUNCTION TO  FUNCTIONAL MATRIX.
!       *STRESSO*   - COMPUTATION NORMALISED WAVE STRESS.
!           !!!!!!! MAKE SURE THAT SINPUT IS CALLED FIRST, STRESSO
!           !!!!!!! NEXT, AND THEN THE REST OF THE SOURCE FUNCTIONS.
!       *FRCUTINDEX*

!     REFERENCE.
!     ----------

!       S. HASSELMANN AND K. HASSELMANN, "A GLOBAL WAVE MODEL",
!       30/6/85 (UNPUBLISHED NOTE)

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWFLUX, LWNEMOCOU
      USE YOWCOUT  , ONLY : LWFLUXOUT 
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,FR5      ,
     &            FRM5     ,COFRM4   ,WETAIL   ,TH       ,
     &            C        ,FRIC     ,FLOGSPRDM1
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWICE   , ONLY : FLMINFR  ,LCIWABR
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,PHIEPS   ,PHIAW    ,
     &                      TAUOC    ,NEMEAN   ,NFMEAN
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,PI       ,EPSMIN
      USE YOWSHAL  , ONLY : DEPTH    ,TCGOND   ,TFAK     ,INDEP 
      USE YOWSTAT  , ONLY : IDELT    ,ISHALLO
      USE YOWTABL  , ONLY : JUMAX    ,DELU
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IJS, IJL, IG
!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL,DIMENSION(IJS:IJL,NANG,NFRE) :: FL,FL3,SL
      REAL,DIMENSION(IJS:IJL,NFRE) :: FCONST, CIWA
      REAL,DIMENSION(IJS:IJL) :: THWOLD,USOLD,Z0OLD,TAUW,
     &                           ROAIRO,ZIDLOLD,CICVR
      REAL,DIMENSION(IJS:IJL) :: U10NEW,THWNEW,USNEW,Z0NEW,
     &                           ROAIRN,ZIDLNEW

! ----------------------------------------------------------------------

      INTEGER :: IJ, K, M
      INTEGER :: ILEV
      INTEGER :: MIJ(IJS:IJL)
      INTEGER :: JU(IJS:IJL)

      REAL :: XJ, DELT, XIMP, DELT5, AKM1, AK2VGM1
      REAL :: GTEMP1, GTEMP2, FLHAB, FLLOWEST
      REAL :: TAU, XN, PHIDIAG 
      REAL :: ZHOOK_HANDLE
      REAL :: DELFL(NFRE)
      REAL, DIMENSION(IJS:IJL) :: EMEANWS, FMEANWS, USFM, GADIAG 
      REAL, DIMENSION(IJS:IJL) :: F1MEAN, AKMEAN, XKMEAN 
      REAL, DIMENSION(IJS:IJL) :: TAUWLF,TAUWD,PHIAWDIAG,PHIAWUNR,PHIOC,
     &                            PHIWA 
      REAL, DIMENSION(IJS:IJL,NANG) :: SPRD
      REAL, DIMENSION(IJS:IJL,NFRE) :: TEMP, TEMP2
      REAL, DIMENSION(IJS:IJL,NANG,NFRE) :: CIWAB 
      REAL, DIMENSION(IJS:IJL,NANG,NFRE) :: XLLWS

      LOGICAL :: LCFLX

      IF (LHOOK) CALL DR_HOOK('IMPLSCH',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

!*    1. INITIALISATION.
!        ---------------

      LCFLX=LWFLUX.OR.LWFLUXOUT.OR.LWNEMOCOU
! ----------------------------------------------------------------------

!*    2. COMPUTATION OF IMPLICIT INTEGRATION.
!        ------------------------------------

!         INTEGRATION IS DONE FROM CDATE UNTIL CDTPRO FOR A BLOCK
!         OF LATITUDES BETWEEN PROPAGATION CALLS.


!*    2.2 COMPUTE MEAN PARAMETERS.
!        ------------------------

      CALL FKMEAN(FL3, IJS, IJL, EMEAN(IJS), FMEAN(IJS),
     &            F1MEAN, AKMEAN, XKMEAN)

      DO K=1,NANG
        DO IJ=IJS,IJL
          SPRD(IJ,K)=MAX(0.,COS(TH(K)-THWNEW(IJ)))**2
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        XJ=U10NEW(IJ)/DELU
        JU(IJ)=MIN(JUMAX, MAX(NINT(XJ),1))
      ENDDO

!     COMPUTE DAMPING COEFFCIENT DUE TO FRICTION ON BOTTOM OF THE SEA ICE.
!!! testing sea ice attenuation (might need to restrict usage when needed)
      IF(LCIWABR) THEN
        CALL CIWABR(IJS, IJL, CICVR, FL3, CIWAB) 
      ELSE
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=IJS,IJL
              CIWAB(IJ,K,M) = 1. 
            ENDDO
          ENDDO
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    2.3 COMPUTATION OF SOURCE FUNCTIONS.
!         --------------------------------

!*      2.3.1 INITIALISE SOURCE FUNCTION AND DERIVATIVE ARRAY.
!             ------------------------------------------------

!!            FL AND SL ARE INITIALISED IN SINPUT

      ILEV=1
      CALL AIRSEA (U10NEW(IJS), TAUW(IJS), USNEW(IJS), Z0NEW(IJS),
     &   IJS, IJL, ILEV)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: AIRSEA CALLED BEFORE DO LOOP'
        CALL FLUSH (IU06)
      ENDIF

!*    2.3.2 ADD SOURCE FUNCTIONS AND WAVE STRESS.
!           -------------------------------------

      CALL SINPUT (FL3, FL, IJS, IJL, THWNEW, USNEW, Z0NEW,
     &             ROAIRN, ZIDLNEW, SL, XLLWS)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: SINPUT CALLED'
        CALL FLUSH (IU06)
      ENDIF

!     MEAN FREQUENCY CHARACTERISTIC FOR WIND SEA
      CALL FEMEANWS(FL3,IJS,IJL,EMEANWS,FMEANWS,XLLWS)

!     COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
      CALL FRCUTINDEX(IJS, IJL, FMEAN(IJS), FMEANWS, CICVR, MIJ)

      CALL STRESSO (FL3, IJS, IJL, THWNEW, USNEW, Z0NEW,
     &              ROAIRN, TAUW, TAUWLF, PHIWA, 
     &              PHIAWDIAG, PHIAWUNR, SL, 
     &              MIJ, .FALSE.)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: STRESSO CALLED'
        CALL FLUSH (IU06)
      ENDIF

      CALL AIRSEA (U10NEW(IJS), TAUW(IJS), USNEW(IJS), Z0NEW(IJS),
     &             IJS, IJL, ILEV)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: AIRSEA CALLED'
        CALL FLUSH (IU06)
      ENDIF

!     2.3.3 ADD THE OTHER SOURCE TERMS.
!           ---------------------------
      CALL SNONLIN (FL3, FL, IJS, IJL, IG, SL, AKMEAN(IJS))
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: SNONLIN CALLED'
        CALL FLUSH (IU06)
      ENDIF
      CALL SDISSIP (FL3 ,FL, IJS, IJL, IG, SL, F1MEAN, XKMEAN,
     &              PHIOC, TAUWD, MIJ, .FALSE.)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: SDISSIP CALLED'
        CALL FLUSH (IU06)
      ENDIF
!SHALLOW
      IF(ISHALLO.NE.1) CALL SBOTTOM (FL3, FL, IJS, IJL, IG, SL)
!SHALLOW

! ----------------------------------------------------------------------

!*    2.4 COMPUTATION OF NEW SPECTRA.
!         ---------------------------

!       INCREASE OF SPECTRUM IN A TIME STEP IS LIMITED TO A FINITE
!       FRACTION OF A TYPICAL F**(-4) EQUILIBRIUM SPECTRUM.

      DELT = IDELT
      XIMP = 1.0
      DELT5 = XIMP*DELT

      DO M=1,NFRE
        DELFL(M) = COFRM4(M)*DELT
      ENDDO
      DO IJ=IJS,IJL
        USFM(IJ) = USNEW(IJ)*MAX(FMEANWS(IJ),FMEAN(IJ))
      ENDDO

      DO M=1,NFRE
        DO IJ=IJS,IJL
          TEMP(IJ,M) = USFM(IJ)*DELFL(M)
        ENDDO
      ENDDO

      IF(ISHALLO.EQ.1) THEN
        DO M=1,NFRE
          DO IJ=IJS,IJL
            TEMP2(IJ,M) = FRM5(M)
          ENDDO
        ENDDO
      ELSE
        DO M=1,NFRE
          DO IJ=IJS,IJL
            AKM1 = 1./TFAK(INDEP(IJ),M)
            AK2VGM1 = AKM1**2/TCGOND(INDEP(IJ),M)
            TEMP2(IJ,M) = AKM1*AK2VGM1
          ENDDO
        ENDDO
      ENDIF

      DO K=1,NANG
        DO M=1,NFRE
          DO IJ=IJS,IJL
            GTEMP1 = MAX((1.-DELT5*FL(IJ,K,M)),1.)
            GTEMP2 = DELT*SL(IJ,K,M)/GTEMP1
            FLHAB = ABS(GTEMP2)
            FLHAB = MIN(FLHAB,TEMP(IJ,M))
            FL3(IJ,K,M) = FL3(IJ,K,M) + SIGN(FLHAB,GTEMP2)
            FLLOWEST = FLMINFR(JU(IJ),M)*SPRD(IJ,K)
            FL3(IJ,K,M)=
     &        MAX(CIWA(IJ,M)*CIWAB(IJ,K,M)*FL3(IJ,K,M),FLLOWEST)
          ENDDO
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    2.5 REPLACE DIAGNOSTIC PART OF SPECTRA BY A F**(-5) TAIL.
!         -----------------------------------------------------


!*    2.5.1 COMPUTE MEAN PARAMETERS.
!           ------------------------

      CALL FKMEAN(FL3, IJS, IJL, EMEAN(IJS), FMEAN(IJS),
     &            F1MEAN, AKMEAN, XKMEAN)

!     MEAN FREQUENCY CHARACTERISTIC FOR WIND SEA
      CALL FEMEANWS(FL3,IJS,IJL,EMEANWS,FMEANWS,XLLWS)

!*    2.5.3 COMPUTE TAIL ENERGY RATIOS.
!           ---------------------------

!     COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
      CALL FRCUTINDEX(IJS, IJL, FMEAN(IJS), FMEANWS, CICVR, MIJ)

      DO IJ=IJS,IJL
        GADIAG(IJ) = 1./TEMP2(IJ,MIJ(IJ))
      ENDDO


!*    2.5.4 MERGE TAIL INTO SPECTRA.
!           ------------------------

      DO IJ=IJS,IJL
        DO M=1,MIJ(IJ)
          FCONST(IJ,M) = 1.
          TEMP(IJ,M) = 0.
        ENDDO
        DO M=MIJ(IJ)+1,NFRE
          FCONST(IJ,M) = 0.
          TEMP(IJ,M) = TEMP2(IJ,M)*GADIAG(IJ)
        ENDDO
      ENDDO

      DO K=1,NANG
        DO IJ=IJS,IJL
          GADIAG(IJ) = FL3(IJ,K,MIJ(IJ))
        ENDDO
        DO M=1,NFRE
          DO IJ=IJS,IJL
            FLLOWEST = FLMINFR(JU(IJ),M)*SPRD(IJ,K)
            FL3(IJ,K,M) = GADIAG(IJ)*TEMP(IJ,M)
     &       + MAX(FL3(IJ,K,M),FLLOWEST)*FCONST(IJ,M)
          ENDDO
        ENDDO
      ENDDO


!*    2.5.5 REEVALUATE WIND INPUT SOURCE TERM, AND WAVE DISSIPATION.
!           -------------------------------------------------------

      CALL SINPUT (FL3, FL, IJS, IJL, THWNEW, USNEW, Z0NEW,
     &             ROAIRN, ZIDLNEW, SL, XLLWS)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: SINPUT CALLED AT THE END'
        CALL FLUSH (IU06)
      ENDIF
      CALL STRESSO (FL3, IJS, IJL, THWNEW, USNEW, Z0NEW,
     &              ROAIRN, TAUW, TAUWLF, PHIWA, 
     &              PHIAWDIAG, PHIAWUNR, SL, MIJ, 
     &              LCFLX)

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: STRESSO CALLED AT THE END'
        CALL FLUSH (IU06)
      ENDIF

      CALL AIRSEA (U10NEW(IJS), TAUW(IJS), USNEW(IJS), Z0NEW(IJS),
     &             IJS, IJL, ILEV)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: AIRSEA CALLED AT THE END'
        CALL FLUSH (IU06)
      ENDIF

      CALL SDISSIP (FL3 ,FL, IJS, IJL, IG, SL, F1MEAN, XKMEAN,
     &              PHIOC, TAUWD, MIJ, LCFLX)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: SDISSIP CALLED AT THE END'
        CALL FLUSH (IU06)
      ENDIF

!*    2.5.6 DETERMINE FLUXES FROM AIR TO WAVE AND FROM WAVE TO OCEAN.
!           -------------------------------------------------------

      IF(LCFLX) THEN
        DO IJ=IJS,IJL
          TAU       = ROAIRN(IJ)*USNEW(IJ)**2
          XN        = ROAIRN(IJ)*USNEW(IJ)**3

          PHIDIAG    = PHIAWDIAG(IJ)+PHIAWUNR(IJ)
          PHIEPS(IJ) = (PHIOC(IJ)-PHIDIAG)/XN 
          PHIAW(IJ)  = (PHIWA(IJ)+PHIAWUNR(IJ))/XN
          TAUOC(IJ)  = (TAU-TAUWLF(IJ)-TAUWD(IJ))/TAU
          NEMEAN(IJ) = EMEAN(IJ)
          NFMEAN(IJ) = FMEAN(IJ)
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    2.6 SAVE WINDS INTO INTERMEDIATE STORAGE.
!         -------------------------------------

      DO IJ=IJS,IJL
        USOLD(IJ) = USNEW(IJ)
        Z0OLD(IJ) = Z0NEW(IJ)
        ROAIRO(IJ) = ROAIRN(IJ)
        ZIDLOLD(IJ) = ZIDLNEW(IJ)
      ENDDO

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('IMPLSCH',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE IMPLSCH
