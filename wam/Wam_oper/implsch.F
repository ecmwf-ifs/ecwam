      SUBROUTINE IMPLSCH (FL3, FL, IJS, IJL, IG,
     1                    THWOLD,USOLD,TAUW,Z0OLD,
     2                    U10NEW,THWNEW,USNEW,Z0NEW,SL,FCONST)

C ----------------------------------------------------------------------
C
C**** *IMPLSCH* - IMPLICIT SCHEME FOR TIME INTEGRATION OF SOURCE
C****             FUNCTIONS.
C
C     S.D.HASSELMANN.  MPI
C     H. GUENTHER AND L. ZAMBRESKY  OPTIMIZATION PERFORMED.
C     H. GUENTHER      GKSS/ECMWF   OCTOBER 1989  NEW WIND FIELD
C                                                 INTERFACE AND
C                                                 TIME COUNTING
C     P.A.E.M. JANSSEN KNMI         AUGUST  1990  COUPLED MODEL
C     H. GUENTHER      GKSS/ECMWF   JUNE    1991  NEW SEPARATION OF
C                                                  DIAG- AND PROGNOSTIC
C                                                  PART OF SPECTRUM.
C     P.A.E.M. JANSSEN ECMWF        FEBRUARY 1995  ADD MINIMUM VALUE (FMIN).
C     J. BIDLOT ECMWF               FEBRUARY 1996 MESSAGE PASSING
C     J. BIDLOT ECMWF               FEBRUARY 1997 MESSAGE PASSING
C
C*    PURPOSE.
C     --------
C
C       THE IMPLICIT SCHEME ENABLES THE USE OF A TIMESTEP WHICH IS
C       LARGE COMPARED WITH THE CHARACTERISTIC DYNAMIC TIME SCALE.
C       THE SCHEME IS REQUIRED FOR THE HIGH FREQUENCIES WHICH
C       RAPIDLY ADJUST TO A QUASI-EQUILIBRIUM.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *IMPLSCH (FL3, FL, IJS, IJL, IG,
C    1                    THWOLD,USOLD,TAUW,Z0OLD,
C    2                    U10NEW,THWNEW,USNEW,Z0NEW,SL,FCONST)
C          *FL3*    - FREQUENCY SPECTRUM(INPUT AND OUTPUT).
C          *FL*     - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
C          *IJS*    - INDEX OF FIRST GRIDPOINT
C          *IJL*    - INDEX OF LAST GRIDPOINT
C          *IG*     - BLOCK NUMBER
C      *U10NEW*    NEW WIND SPEED IN M/S.
C      *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
C                  NOTATION (POINTING ANGLE OF WIND VECTOR,
C                  CLOCKWISE FROM NORTH).
C      *THWOLD*    INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
C                  WIND VELOCITY.
C      *USNEW*     NEW FRICTION VELOCITY IN M/S.
C      *USOLD*     INTERMEDIATE STORAGE OF MODULUS OF FRICTION
C                  VELOCITY.
C      *Z0NEW*     ROUGHNESS LENGTH IN M.
C      *Z0OLD*     INTERMEDIATE STORAGE OF ROUGHNESS LENGTH IN
C                  M.
C      *TAUW*      WAVE STRESS IN (M/S)**2
C      *SL*        REAL      TOTAL SOURCE FUNCTION ARRAY.
C      *FCONST*    REAL      = 1 FOR PROGNOSTIC FREQUENCY BANDS.
C                            = 0 FOR DIAGNOSTIC FREQUENCY BANDS.
C
C     METHOD.
C     -------
C
C       THE SPECTRUM AT TIME (TN+1) IS COMPUTED AS
C       FN+1=FN+DELT*(SN+SN+1)/2., WHERE SN IS THE TOTAL SOURCE
C       FUNCTION AT TIME TN, SN+1=SN+(DS/DF)*DF - ONLY THE DIAGONAL
C       TERMS OF THE FUNCTIONAL MATRIX DS/DF ARE COMPUTED, THE
C       NONDIAGONAL TERMS ARE NEGLIGIBLE.
C       THE ROUTINE IS CALLED AFTER PROPAGATION FOR TIME PERIOD
C       BETWEEN TWO PROPAGATION CALLS - ARRAY FL3 CONTAINS THE
C       SPECTRUM AND FL IS USED AS AN INTERMEDIATE STORAGE FOR THE
C       DIAGONAL TERM OF THE FUNCTIONAL MATRIX.
C
C     EXTERNALS.
C     ---------
C
C       *AIRSEA*    - SURFACE LAYER STRESS AND ROUGHNESS LENGTH.
C       *CREWFN*    - CREATE A WIND FILE NAME.
C       *FEMEAN*    - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT.
C       *INCDATE*   - UPDATE DATE TIME GROUP.
CSHALLOW
C       *SBOTTOM*   - COMPUTES BOTTOM DISSIPATION SOURCE TERM AND
C                     LINEAR CONTRIBUTION TO FUNCTIONAL MATRIX.
CSHALLOW
C       *SDISSIP*   - COMPUTATION OF DISSIPATION SOURCE FUNCTION
C                     AND LINEAR CONTRIBUTION OF DISSIPATION TO
C                     FUNCTIONAL MATRIX IN IMPLICIT SCHEME.
C       *SEMEAN*    - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.
C       *SINPUT*    - COMPUTATION OF INPUT SOURCE FUNCTION, AND
C                     LINEAR CONTRIBUTION OF INPUT SOURCE FUNCTION
C                     TO FUNCTIONAL MATRIX IN IMPLICIT SCHEME.
C       *SNONLIN*   - COMPUTATION OF NONLINEAR TRANSFER RATE AND
C                     DIAGONAL LINEAR CONTRIBUTION OF NONLINEAR SOURCE
C                     FUNCTION TO  FUNCTIONAL MATRIX.
C       *STRESSO*   - COMPUTATION NORMALISED WAVE STRESS.
C           !!!!!!! MAKE SURE THAT SINPUT IS CALLED FIRST, STRESSO
C           !!!!!!! NEXT, AND THEN THE REST OF THE SOURCE FUNCTIONS.
C
C     REFERENCE.
C     ----------
C
C       S. HASSELMANN AND K. HASSELMANN, "A GLOBAL WAVE MODEL",
C       30/6/85 (UNPUBLISHED NOTE)
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "comfred.h"
C
#include "commean.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
#include "comstat.h"
C
#include "comtest.h"
C
#include "comunit.h"
C
#include "comwind.h"
C
C ----------------------------------------------------------------------
C
C     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL,FL3,SL
      REAL,DIMENSION(NINF:NSUP,NFRE) :: FCONST 
      REAL,DIMENSION(NINF:NSUP,NBLO) :: THWOLD,USOLD,Z0OLD,TAUW
      REAL,DIMENSION(NINF:NSUP) :: U10NEW,THWNEW,USNEW,Z0NEW
C
C ----------------------------------------------------------------------
C
#include "parcons.h"
C
      PARAMETER (GZPI28 = G/28./ZPI)
C
      LOGICAL NEWREAD, NEWFILE
      DIMENSION MIJ(IJS:IJL), MFMF(IJS:IJL), GADIAG(IJS:IJL),
     1          TEMP(IJS:IJL,NFRE), DELFL(NFRE)
C
C ----------------------------------------------------------------------
C
C*    1. INITIALISATION.
C        ---------------
C
 1000 CONTINUE
C 
      FMIN  = 0.0001   ! THIS MAKES HS=3cm WHEN U10 = 0.0
C
C ----------------------------------------------------------------------
C
C*    2. COMPUTATION OF IMPLICIT INTEGRATION.
C        ------------------------------------
C
C         INTEGRATION IS DONE FROM CDATE UNTIL CDTPRO FOR A BLOCK
C         OF LATITUDES BETWEEN PROPAGATION CALLS.
C
 2000 CONTINUE
C
C*    2.1 CALCULATE ROUGHNESS LENGTH AND FRICTION VELOCITIES.
C         ---------------------------------------------------
C
         CALL AIRSEA (U10NEW(IJS), TAUW(IJS,IG), USNEW(IJS), Z0NEW(IJS),
     1                IJS, IJL)
C
C ----------------------------------------------------------------------
C
C*    2.2 COMPUTE MEAN PARAMETERS.
C         ------------------------
C
 2200 CONTINUE
         CALL SEMEAN(FL3, IJS, IJL)
         CALL FEMEAN(FL3, IJS, IJL)
C
C ----------------------------------------------------------------------
C
C*    2.3 COMPUTATION OF SOURCE FUNCTIONS.
C         --------------------------------
C
 2300 CONTINUE
C
C*    2.3.1 INITIALISE SOURCE FUNCTION AND DERIVATIVE ARRAY.
C           ------------------------------------------------
C
         DO M=1,NFRE
           DO K=1,NANG
             DO IJ=IJS,IJL
               SL(IJ,K,M) = 0.
               FL(IJ,K,M) = 0.
             ENDDO
           ENDDO
         ENDDO
C
C*    2.3.2 ADD SOURCE FUNCTIONS AND WAVE STRESS.
C           -------------------------------------
C
         CALL SINPUT (FL3, FL, IJS, IJL, IG,THWOLD,USOLD,Z0OLD,
     1                THWNEW,USNEW,Z0NEW,SL)
         CALL STRESSO (FL3,IJS,IJL,IG,THWNEW,USNEW,Z0NEW,TAUW,SL)
         CALL SNONLIN (FL3, FL, IJS, IJL, IG, SL)
         CALL SDISSIP (FL3 ,FL, IJS, IJL, SL)
CSHALLOW
         IF(ISHALLO.NE.1) CALL SBOTTOM (FL3, FL, IJS, IJL, IG, SL)
CSHALLOW
C
C ----------------------------------------------------------------------
C
C
C*    2.4 COMPUTATION OF NEW SPECTRA.
C         ---------------------------
C
C       INCREASE OF SPECTRUM IN A TIME STEP IS LIMITED TO A FINITE
C       FRACTION OF A TYPICAL F**(-4) EQUILIBRIUM SPECTRUM.
C
 2400 CONTINUE
         DELT = IDELT
         XIMP = 1.0
         DELT5 = XIMP*DELT
         DO 2401 M=1,NFRE
            DELFL(M) = 5.E-07*G*FR(M)**(-4.)*DELT
            DO 2402 K=1,NANG
               DO 2403 IJ=IJS,IJL
                  GTEMP1 = MAX((1.-DELT5*FL(IJ,K,M)),1.)
                  GTEMP2 = DELT*SL(IJ,K,M)/GTEMP1
                  FLHAB = ABS(GTEMP2)
                  FAC   = USNEW(IJ)*FMEAN(IJ)*DELFL(M)
                  FLHAB = MIN(FLHAB,FAC)
                  FL3(IJ,K,M) = FL3(IJ,K,M) + SIGN(FLHAB,GTEMP2)
                  FL3(IJ,K,M) = MAX(FL3(IJ,K,M),FMIN)
 2403          CONTINUE
 2402       CONTINUE
 2401    CONTINUE
C
C ----------------------------------------------------------------------
C
C*    2.5 REPLACE DIAGNOSTIC PART OF SPECTRA BY A F**(-5) TAIL.
C         -----------------------------------------------------
C
 2500 CONTINUE
C
C*    2.5.1 COMPUTE MEAN PARAMETERS.
C           ------------------------
C
         CALL SEMEAN(FL3, IJS, IJL)
         CALL FEMEAN(FL3, IJS, IJL)
C
C*    2.5.2 COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
C*          FREQUENCIES LE MAX(4*F(PM) , 2.5*FMEAN).
C           ------------------------------------------------------------
C
         FPMH = 2.5/FR(1)
         FPM = 4.*GZPI28/FR(1)
         DO 2521 IJ=IJS,IJL
            FPM4 = FPM/(USNEW(IJ)+EPSMIN)
            MIJ(IJ) = ALOG10(FPM4)*24.1589+2.
            FPM4 = FMEAN(IJ)*FPMH
            MFMF(IJ) = ALOG10(FPM4)*24.1589+1.
 2521    CONTINUE
         DO 2522 IJ=IJS,IJL
            MIJ(IJ) = MAX(MFMF(IJ),MIJ(IJ))
            MIJ(IJ) = MIN(MIJ(IJ),NFRE)
 2522    CONTINUE
C
C*    2.5.3 COMPUTE TAIL ENERGY RATIOS.
C           ---------------------------
C
         DO 2531 M=1,NFRE
            DELFL(M) = (1./FR(M))**5.
 2531    CONTINUE
         DO 2532 IJ=IJS,IJL
            GADIAG(IJ) = FR(MIJ(IJ))**5.
 2532    CONTINUE
C
C*    2.5.4 MERGE TAIL INTO SPECTRA.
C           ------------------------
C
         DO 2541 M=1,NFRE
            DO 2542 IJ=IJS,IJL
               FCONST(IJ,M) = 0.
               TEMP(IJ,M) = GADIAG(IJ)*DELFL(M)
 2542       CONTINUE
 2541    CONTINUE
         DO 2543 IJ=IJS,IJL
            J = MIJ(IJ)
            DO 2544 M=1,J
               FCONST(IJ,M) = 1.
               TEMP(IJ,M) = 0.
 2544       CONTINUE
 2543    CONTINUE
C
         DO 2545 K=1,NANG
            DO 2546 IJ=IJS,IJL
               GADIAG(IJ) = FL3(IJ,K,MIJ(IJ))
 2546       CONTINUE
            DO 2547 M=1,NFRE
               DO 2548 IJ=IJS,IJL
                   FL3(IJ,K,M) = GADIAG(IJ)*TEMP(IJ,M)
     1                         + FL3(IJ,K,M)*FCONST(IJ,M)
 2548          CONTINUE
 2547       CONTINUE
 2545    CONTINUE
C
C ----------------------------------------------------------------------
C
C*    2.6 SAVE WINDS INTO INTERMEDIATE STORAGE.
C         -------------------------------------
C
 2600 CONTINUE
         DO 2601 IJ=IJS,IJL
            USOLD(IJ,IG) = USNEW(IJ)
            Z0OLD(IJ,IG) = Z0NEW(IJ)
 2601    CONTINUE
C
C ----------------------------------------------------------------------
C
C
      RETURN
      END
