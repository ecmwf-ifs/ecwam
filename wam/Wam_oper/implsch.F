      SUBROUTINE IMPLSCH (FL3, FL, IJS, IJL, IG,
     &                    THWOLD,USOLD,TAUW,Z0OLD,
     &                    U10NEW,THWNEW,USNEW,Z0NEW,SL,FCONST)

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

!*    PURPOSE.
!     --------

!       THE IMPLICIT SCHEME ENABLES THE USE OF A TIMESTEP WHICH IS
!       LARGE YOWPARED WITH THE CHARACTERISTIC DYNAMIC TIME SCALE.
!       THE SCHEME IS REQUIRED FOR THE HIGH FREQUENCIES WHICH
!       RAPIDLY ADJUST TO A QUASI-EQUILIBRIUM.

!**   INTERFACE.
!     ----------

!       *CALL* *IMPLSCH (FL3, FL, IJS, IJL, IG,
!    1                    THWOLD,USOLD,TAUW,Z0OLD,
!    2                    U10NEW,THWNEW,USNEW,Z0NEW,SL,FCONST)
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

!     REFERENCE.
!     ----------

!       S. HASSELMANN AND K. HASSELMANN, "A GLOBAL WAVE MODEL",
!       30/6/85 (UNPUBLISHED NOTE)

! ----------------------------------------------------------------------

      USE YOWPCONS , ONLY : G        ,ZPI      ,EPSMIN
      USE YOWCOUP  , ONLY : LWCOU
      USE YOWFRED  , ONLY : FR
      USE YOWICE   , ONLY : FLMIN
      USE YOWMEAN  , ONLY : FMEAN
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO
      USE YOWSTAT  , ONLY : IDELT    ,ISHALLO
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------

!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL,FL3,SL
      REAL,DIMENSION(NINF:NSUP,NFRE) :: FCONST 
      REAL,DIMENSION(NINF:NSUP,NBLO) :: THWOLD,USOLD,Z0OLD,TAUW
      REAL,DIMENSION(NINF:NSUP) :: U10NEW,THWNEW,USNEW,Z0NEW

! ----------------------------------------------------------------------

      PARAMETER (GZPI28 = G/28./ZPI)

      DIMENSION MIJ(IJS:IJL), MFMF(IJS:IJL), GADIAG(IJS:IJL),
     &          TEMP(IJS:IJL,NFRE), DELFL(NFRE)

! ----------------------------------------------------------------------

!*    1. INITIALISATION.
!        ---------------

! ----------------------------------------------------------------------

!*    2. COMPUTATION OF IMPLICIT INTEGRATION.
!        ------------------------------------

!         INTEGRATION IS DONE FROM CDATE UNTIL CDTPRO FOR A BLOCK
!         OF LATITUDES BETWEEN PROPAGATION CALLS.


!*    2.2 COMPUTE MEAN PARAMETERS.
!         ------------------------

      CALL SEMEAN(FL3, IJS, IJL)
      CALL FEMEAN(FL3, IJS, IJL)

! ----------------------------------------------------------------------

!*    2.3 COMPUTATION OF SOURCE FUNCTIONS.
!         --------------------------------
!!!! we will limit ourselves to 1 iteration for now !!!
      NCALL=1

      DO ICALL=1,NCALL

!*      2.3.1 INITIALISE SOURCE FUNCTION AND DERIVATIVE ARRAY.
!             ------------------------------------------------

!!            FL and SL are initialised in sinput

!*      2.3.2 ADD SOURCE FUNCTIONS AND WAVE STRESS.
!             -------------------------------------
        ILEV=1
        CALL AIRSEA (U10NEW(IJS), TAUW(IJS,IG), USNEW(IJS), Z0NEW(IJS),
     &   IJS, IJL, ILEV)
        IF (ITEST.GE.2) THEN
          WRITE(IU06,*) '   SUB. IMPLSCH: AIRSEA CALLED'
          CALL FLUSH (IU06)
        ENDIF

        CALL SINPUT (FL3, FL, IJS, IJL, IG,THWOLD,USOLD,Z0OLD,
     &   THWNEW,USNEW,Z0NEW,SL,ICALL,NCALL)
        IF (ITEST.GE.2) THEN
          WRITE(IU06,*) '   SUB. IMPLSCH: SINPUT CALLED'
          CALL FLUSH (IU06)
        ENDIF

        CALL STRESSO (FL3,IJS,IJL,IG,THWNEW,USNEW,Z0NEW,TAUW,SL)
        IF (ITEST.GE.2) THEN
          WRITE(IU06,*) '   SUB. IMPLSCH: STRESSO CALLED'
          CALL FLUSH (IU06)
        ENDIF

      ENDDO

!     2.3.3 ADD THE OTHER SOURCE TERMS.
!           ---------------------------
      CALL SNONLIN (FL3, FL, IJS, IJL, IG, SL)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: SNOLIN CALLED'
        CALL FLUSH (IU06)
      ENDIF
      CALL SDISSIP (FL3 ,FL, IJS, IJL, SL)
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
        DELFL(M) = 5.E-07*G*FR(M)**(-4.)*DELT
        DO K=1,NANG
          DO IJ=IJS,IJL
            GTEMP1 = MAX((1.-DELT5*FL(IJ,K,M)),1.)
            GTEMP2 = DELT*SL(IJ,K,M)/GTEMP1
            FLHAB = ABS(GTEMP2)
            FAC   = USNEW(IJ)*FMEAN(IJ)*DELFL(M)
            FLHAB = MIN(FLHAB,FAC)
            FL3(IJ,K,M) = FL3(IJ,K,M) + SIGN(FLHAB,GTEMP2)
            FL3(IJ,K,M) = MAX(FL3(IJ,K,M),FLMIN)
          ENDDO
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    2.5 REPLACE DIAGNOSTIC PART OF SPECTRA BY A F**(-5) TAIL.
!         -----------------------------------------------------


!*    2.5.1 COMPUTE MEAN PARAMETERS.
!           ------------------------

      CALL SEMEAN(FL3, IJS, IJL)
      CALL FEMEAN(FL3, IJS, IJL)

!*    2.5.2 COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!*          FREQUENCIES LE 2.5*FMEAN.
!           ------------------------------------------------------------

      FPMH = 2.5/FR(1)
      DO IJ=IJS,IJL
        FPM4 = FMEAN(IJ)*FPMH
        MIJ(IJ) = ALOG10(FPM4)*24.1589+1.
        MIJ(IJ) = MIN(MIJ(IJ),NFRE)
      ENDDO

!*    2.5.3 COMPUTE TAIL ENERGY RATIOS.
!           ---------------------------

      DO M=1,NFRE
        DELFL(M) = (1./FR(M))**5.
      ENDDO
      DO IJ=IJS,IJL
        GADIAG(IJ) = FR(MIJ(IJ))**5.
      ENDDO

!*    2.5.4 MERGE TAIL INTO SPECTRA.
!           ------------------------

      DO M=1,NFRE
        DO IJ=IJS,IJL
          FCONST(IJ,M) = 0.
          TEMP(IJ,M) = GADIAG(IJ)*DELFL(M)
        ENDDO
      ENDDO
      DO IJ=IJS,IJL
        J = MIJ(IJ)
        DO M=1,J
          FCONST(IJ,M) = 1.
          TEMP(IJ,M) = 0.
        ENDDO
      ENDDO

      DO K=1,NANG
        DO IJ=IJS,IJL
          GADIAG(IJ) = FL3(IJ,K,MIJ(IJ))
        ENDDO
        DO M=1,NFRE
          DO IJ=IJS,IJL
            FL3(IJ,K,M) = GADIAG(IJ)*TEMP(IJ,M)
     &       + FL3(IJ,K,M)*FCONST(IJ,M)
          ENDDO
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    2.6 SAVE WINDS INTO INTERMEDIATE STORAGE.
!         -------------------------------------

      DO IJ=IJS,IJL
        USOLD(IJ,IG) = USNEW(IJ)
        Z0OLD(IJ,IG) = Z0NEW(IJ)
      ENDDO

! ----------------------------------------------------------------------

      RETURN
      END SUBROUTINE IMPLSCH
