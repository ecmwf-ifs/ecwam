      SUBROUTINE INTSPEC (NFRE, NANG, ML, KL, FR, DEL12, DEL1L,
     1                    F1, FMEAN1, EMEAN1, THETM1,
     2                    F2, FMEAN2, EMEAN2, THETM2,
     3                    FL, FMEAN,  EMEAN,  THETM )

C ----------------------------------------------------------------------
C
C**** *INTSPEC* -  INTERPOLATION OF SPECTRA.
C
C     SUSANNE HASSELMANN  MPI        JUNE 1990.
C     H. GUNTHER          GKSS/ECMWF JAN. 1991   MODIFIED FOR CYCLE_4
C
C*    PURPOSE.
C     --------
C
C       INTERPOLATION OF SPECTRA.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *INTSPEC (NFRE, NANG, ML, KL, FR, DEL12, DEL1L,
C                        F1, FMEAN1, EMEAN1, THETM1,
C                        F2, FMEAN2, EMEAN2, THETM2,
C                        FL, FMEAN,  EMEAN,  THETM )*
C         *NFRE*   - FREQUENCY DIMENSION OF SPECTRA.
C         *NANG*   - DIRECTION DIMENSION OF SPECTRA.
C         *ML*     - NUMBER OF FREQUENCIES.
C         *KL*     - NUMBER OF DIRECTIONS.
C         *FR*     - FREQUENCY ARRAY.
C         *DEL12*  - DISTANCE SPECTRUM 2 - SPECTRUM 1.
C         *DEL1L*  - DISTANCE SPECTRUM L - SPECTRUM 1.
C         *F1*     - INPUT SPECTRUM 1.
C         *FMEAN1* - INPUT MEAN FREQUENCY OF F1.
C         *EMEAN1* - INPUT MEAN ENERGY OF F1.
C         *THETM1* - INPUT MEAN DIRECTION OF F1.
C         *F2*     - INPUT SPECTRUM 2.
C         *FMEAN2* - INPUT MEAN FREQUENCY OF F2.
C         *EMEAN2* - INPUT MEAN ENERGY OF F2.
C         *THETM2* - INPUT MEAN DIRECTION OF F2.
C         *FL*     - INTEPOLATED SPECTRUM.
C         *FMEAN*  - INTEPOLATED MEAN FREQUENCY.
C         *EMEAN*  - INTEPOLATED MEAN ENERGY.
C         *THETM*  - INTEPOLATED MEAN DIRECTION.
C
C     METHOD.
C     -------
C
C       ROTATE SPECTRA ACCORDING TO MEAN OF MEAN ANGLES, TRANSFORM
C       FREQUENCIES ACCORDING TO MEAN OF MEAN FREQUENCIES ,ADJUST ENERGY
C       ACCORDCING TO MEAN OF TOTAL ENERGY AND INTERPOLATE RESULTING
C       SPECTRA.
C
C     EXTERNALS.
C     ----------
C
C       *ROTSPEC*   - ROTATES SPECTRUM.
C       *STRSPEC*   - TRANSFORM FREQUENCIES.
C
C     REFERENCES.
C     -----------
C
C       K.HASSELMANN, 1990,
C          INTERPOLATION OF WAVE SPECTRA. WAM NOTE 6/6/90.
C
C ----------------------------------------------------------------------
C
      DIMENSION F1(NANG,NFRE), F2(NANG,NFRE), FL(NANG, NFRE), FR(NFRE)
C
C ----------------------------------------------------------------------
C
#include "parcons.h"
C
      DIMENSION F3(NANG,NFRE), F4(NANG,NFRE)
C
C ----------------------------------------------------------------------
C
C*    1. INTERPOLATION WEIGHTS.
C        ----------------------
C
 1000 CONTINUE
C
      GW1 = (DEL12-DEL1L)/DEL12
      GW2 = DEL1L/DEL12
C
C ----------------------------------------------------------------------
C
C*    2. INTERPOLATE MEAN VALUES.
C        ------------------------
C
 2000 CONTINUE
C
      IF (EMEAN1.EQ.0.) THEN
C
C*    2.1 ENERGY OF SPECTRUM 1 IS ZERO.
C         -----------------------------
C
         EMEAN = GW2*EMEAN2
         FMEAN = FMEAN2
         THETM = THETM2
         DO 1101 M=1,ML
            DO 1102 K=1,KL
               FL(K,M) = GW2*F2(K,M)
 1102       CONTINUE
 1101    CONTINUE

         RETURN

      ELSE IF (EMEAN2.EQ.0.) THEN
C
C*    2.2 ENERGY OF SPECTRUM 2 IS ZERO.
C         -----------------------------
C
         EMEAN = GW1*EMEAN1
         FMEAN = FMEAN1
         THETM = THETM1
         DO 1201 M=1,ML
            DO 1202 K=1,KL
               FL(K,M) = GW1*F1(K,M)
 1202       CONTINUE
 1201    CONTINUE

         RETURN

      ENDIF
C
C*    2.3 ENERGY BOTH SPECTRUM IS GT ZERO.
C         -------------------------------
C
      EMEAN = GW1*EMEAN1+GW2*EMEAN2
      FMEAN = GW1*FMEAN1+GW2*FMEAN2
      CM1 = COS(THETM1)
      CM2 = COS(THETM2)
      SM1 = SIN(THETM1)
      SM2 = SIN(THETM2)
      CM  = GW1*CM1 + GW2*CM2
      SM  = GW1*SM1 + GW2*SM2
      THETM = ATAN2(SM,CM)
      THETM = MOD(THETM+ZPI,ZPI)
C
C ----------------------------------------------------------------------
C
C*    3. SPECTRUM 1.
C        -----------
C
 3000 CONTINUE
C
C*    3.1 ROTATE.
C         -------
C
      RTHET1 = THETM - THETM1
      CALL ROTSPEC (NFRE, NANG, ML, KL, F1, F3, RTHET1)
 160  CONTINUE
C
C*    3.2 STRETCH.
C         --------
C
      GAMMA = FMEAN1/FMEAN
      CALL STRSPEC (NFRE, NANG, ML, KL, FR, F3, GAMMA)
C
C*    3.3 ADJUST ENERGY.
C         --------------
C
      EMEANH = EMEAN/EMEAN1
      DO 3301 M=1,ML
         DO 3302 K=1,KL
            F3(K,M) = F3(K,M)*EMEANH
 3302    CONTINUE
 3301 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    4. SPECTRUM 2.
C        -----------
C
 4000 CONTINUE
C
C*    4.1 ROTATE.
C         -------
C
      RTHET2 = THETM-THETM2
      CALL ROTSPEC (NFRE, NANG, ML, KL, F2, F4, RTHET2)
C
C*    4.2 STRETCH.
C         --------
C
      GAMMA = FMEAN2/FMEAN
      CALL STRSPEC (NFRE, NANG, ML, KL, FR, F4, GAMMA)
C
C*    4.3 ADJUST ENERGY.
C         --------------
C
      EMEANH = EMEAN/EMEAN2
      DO 4301 M=1,ML
         DO 4302 K=1,KL
            F4(K,M) = F4(K,M)*EMEANH
 4302    CONTINUE
 4301 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    5. LINEAR INTERPOLATION OF NEW SPECTRA.
C        ------------------------------------
C
 5000 CONTINUE
C
      DO 5001 M=1,ML
         DO 5002 K=1,KL
            FL(K,M) = GW1*F3(K,M) + GW2*F4(K,M)
 5002    CONTINUE
 5001 CONTINUE

      RETURN
      END
