! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

      SUBROUTINE LFACTOR(S, CINV, U10, USTAR, UPROXY, USDIR, ROAIRN, &
     &                   LFACT, LREDUCE, TAUWX, TAUWY, TAU)
 
! ----------------------------------------------------------------------------
!
!  1. Purpose :
!
!      Numerical approximation for the reduction factor LFACTOR(f) to
!      reduce energy in the high-frequency part of the resolved part
!      of the spectrum to meet the constraint on total stress (TAU).
!      The constraint is TAU <= TAU_TOT (TAU_TOT = TAU_WAV + TAU_VIS),
!      thus the wind input is reduced to match our constraint.
!
!  2. Method :
!
!     1) If required, extend resolved part of the spectrum to 10Hz using
!        an approximation for the spectral slope at the high frequency
!        limit: Sin(F) prop. F**(-2) and for E(F) prop. F**(-5).
!     2) Calculate stresses:
!        total stress:     TAU_TOT  = DAIR * USTAR**2
!        viscous stress:   TAU_VIS  = DAIR * Cv * U10**2
!        viscous stress (x,y-components):
!                          TAUV_X   = TAU_VIS * COS(USDIR)
!                          TAUV_Y   = TAU_VIS * SIN(USDIR)
!        wave supported stress (x,y-components):    /10Hz
!                          TAUW_X,Y = GRAV * DWAT * | [SinX,Y(F)]/C(F) dF
!                                                   /
!        total stress (input):   TAU  = SQRT( (TAUW_X + TAUV_X)**2
!                                     + (TAUW_Y + TAUV_Y)**2 )
!     3) If TAU does not meet our constraint reduce the wind input
!        using reduction factor:
!                          LFACT(F) = MIN(1,exp((1-U/C(F))*RTAU))
!        Then alter RTAU and repeat 3) until our constraint is matched.
!
! ----------------------------------------------------------------------------
!
!**   INTERFACE.
!     ----------

!     *CALL* *LFACTOR(S, CINV, U10, USTAR, USDIR, ROAIRN, SIG, DSII, &
!     &                   LFACT, TAUWX, TAUWY, TAU)
!            *S*   - NEG. WIND INPUT ENERGY DENSITY SPECTRUM.
!          *CINV*  - INVERSE PHASE SPEED CALC. IN INPUT ROUTINE
!         *UABS*   - 10M WIND SPEED
!        *USTAR*   - NEW FRICTION VELOCITY IN M/S.
!        *USDIR*   - WIND DIRECTION
!       *ROAIRN*   - AIR DENSITY IN KG/M3
!          *LFACT* - CORRECTION FACTOR
!          *TAUNWX, TAUNWY* - NEGATIVE WAVE NORMAL STRESS COMPONENTS

!     EXTERNALS.
!     ----------
!     TAUWINDS

!     ORIGIN.
!     ----------
!     Adapted from Babanin Young Donelan & Banner (ZBRY) physics 
!     as implemented as ST6 in WAVEWATCH-III 
!     Implementation into ECWAM DECEMBER 2021 by J. Kousal 

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,TH       ,DFIM     ,COSTH  ,SINTH,&
&                           FRATIO   ,DELTH    ,FRIC, SIG,DSII   ,SIGM1,&
&                           DF       ,NFRE_EXT ,DSII_EXT ,SIG_EXT
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER  ,EPSMIN, GM1
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "tauwinds.intfb.h"
 
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(IN)  :: S ! Sin(sigma) in [m2/rad-Hz]
      REAL(KIND=JWRB), DIMENSION(NFRE),      INTENT(IN)  :: CINV
      REAL(KIND=JWRB),                       INTENT(IN)  :: U10, USTAR, UPROXY, USDIR, ROAIRN

      REAL(KIND=JWRB), DIMENSION(NFRE),      INTENT(OUT) :: LFACT
      LOGICAL,                               INTENT(OUT) :: LREDUCE
      REAL(KIND=JWRB),                       INTENT(OUT) :: TAUWX, TAUWY, TAU

      INTEGER(KIND=JWIM), PARAMETER :: ITERMAX = 80 ! Max. no. iterations
                                                    ! to find numerical LFACT soln

      REAL(KIND=JWRB), DIMENSION(NFRE_EXT) :: LF_EXT, CINV_EXT
      REAL(KIND=JWRB), DIMENSION(NFRE_EXT) :: SDENS_EXT, SDENSX_EXT, SDENSY_EXT
      REAL(KIND=JWRB), DIMENSION(NFRE_EXT) :: UCINV_EXT

      REAL(KIND=JWRB)      :: TAU_TOT, TAU_VIS, TAU_WAV
      REAL(KIND=JWRB)      :: TAUVX, TAUVY, TAUX, TAUY   
      REAL(KIND=JWRB)      :: RTAU, DRTAU, ERR   
      LOGICAL              :: OVERSHOT

      INTEGER(KIND=JWIM) :: IK, K, M, SIGN_NEW, SIGN_OLD

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('LFACTOR',0,ZHOOK_HANDLE)

!/ 1) --- Either extrapolate arrays up to 10Hz or use discrete spectral
!         grid per se. Limit the constraint to the positive part of the
!         wind input only. ---------------------------------------------- /
      IF (NFRE .LT. NFRE_EXT) THEN
            CINV_EXT(1:NFRE)          = CINV
            SDENSX_EXT(1:NFRE)        = 0.0_JWRB
            SDENSY_EXT(1:NFRE)        = 0.0_JWRB
            SDENS_EXT(1:NFRE)         = 0.0_JWRB
            DO K = 1, NANG
                  DO M = 1, NFRE
                        SDENS_EXT(M) = SDENS_EXT(M) + S(K,M)
                        SDENSX_EXT(M) = SDENSX_EXT(M) + MAX(0.0_JWRB,S(K,M))*COSTH(K)
                        SDENSY_EXT(M) = SDENSY_EXT(M) + MAX(0.0_JWRB,S(K,M))*SINTH(K)
                  END DO
            END DO
            SDENS_EXT(1:NFRE)         = SDENS_EXT(1:NFRE) * DELTH
            SDENSX_EXT(1:NFRE)        = SDENSX_EXT(1:NFRE) * DELTH
            SDENSY_EXT(1:NFRE)        = SDENSY_EXT(1:NFRE) * DELTH
!        --- Spectral slope for S_IN(F) is proportional to F**(-2) ------ /
            CINV_EXT(NFRE+1:NFRE_EXT)   = SIG_EXT(NFRE+1:NFRE_EXT)*GM1 ! 1/c=σ/g
            SDENS_EXT(NFRE+1:NFRE_EXT)  = SDENS_EXT(NFRE)  * (SIG_EXT(NFRE)/SIG_EXT(NFRE+1:NFRE_EXT))**2
            SDENSX_EXT(NFRE+1:NFRE_EXT) = SDENSX_EXT(NFRE) * (SIG_EXT(NFRE)/SIG_EXT(NFRE+1:NFRE_EXT))**2
            SDENSY_EXT(NFRE+1:NFRE_EXT) = SDENSY_EXT(NFRE) * (SIG_EXT(NFRE)/SIG_EXT(NFRE+1:NFRE_EXT))**2
      ELSE
            CINV_EXT         = CINV
            SDENSX_EXT(1:NFRE) = 0.0_JWRB
            SDENSY_EXT(1:NFRE) = 0.0_JWRB
            SDENS_EXT(1:NFRE)  = 0.0_JWRB
            DO K = 1, NANG
                  DO M = 1, NFRE
                        SDENS_EXT(M) = SDENS_EXT(M) + S(K,M)
                        SDENSX_EXT(M) = SDENSX_EXT(M) + MAX(0.0_JWRB,S(K,M))*COSTH(K)
                        SDENSY_EXT(M) = SDENSY_EXT(M) + MAX(0.0_JWRB,S(K,M))*SINTH(K)
                  END DO
            END DO
            SDENS_EXT(1:NFRE) = SDENS_EXT(1:NFRE) * DELTH
            SDENSX_EXT(1:NFRE) = SDENSX_EXT(1:NFRE) * DELTH
            SDENSY_EXT(1:NFRE) = SDENSY_EXT(1:NFRE) * DELTH
      END IF
!
!/ 2) --- Stress calculation ----------------------------------------- /
!     --- The total stress ------------------------------------------- /
      TAU_TOT  = USTAR**2 * ROAIRN
!
!     --- The viscous stress and check that it does not exceed
!         the total stress. ------------------------------------------ /
      TAU_VIS  = MAX(0.0_JWRB, -5.0E-5_JWRB*U10 + 1.1E-3_JWRB) * U10**2 * ROAIRN
      TAU_VIS  = MIN(0.95_JWRB * TAU_TOT, TAU_VIS)
!
      TAUVX    = TAU_VIS * COS(USDIR)
      TAUVY    = TAU_VIS * SIN(USDIR)
!
!     --- The wave supported stress. --------------------------------- /
      TAUWX    = TAUWINDS(SDENSX_EXT,CINV_EXT,DSII_EXT)   ! normal stress (x-component)
      TAUWY    = TAUWINDS(SDENSY_EXT,CINV_EXT,DSII_EXT)   ! normal stress (y-component)
      TAU_WAV  = SQRT(TAUWX**2 + TAUWY**2)        ! normal stress (magnitude)
!
      TAUX     = TAUVX + TAUWX                        ! total stress (x-component)
      TAUY     = TAUVY + TAUWY                        ! total stress (y-component)
      TAU      = SQRT(TAUX**2  + TAUY**2)                 ! total stress (magnitude)
      ERR      = (TAU-TAU_TOT)/TAU_TOT                    ! initial error
!
!/ 3) --- Find reduced Sin(f) = L(f)*Sin(f) to satisfy our constraint
!/        TAU <= TAU_TOT --------------------------------------------- /
      LF_EXT = 1.0_JWRB
      LREDUCE = .FALSE.
      IK = 0
!
      IF (TAU .GT. TAU_TOT) THEN

         OVERSHOT    = .FALSE.
         RTAU        = ERR / 90.0_JWRB
         DRTAU       = 2.0_JWRB

         SIGN_NEW    = INT(SIGN(1.0_JWRB,ERR))

         UCINV_EXT   = 1.0_JWRB - (UPROXY * CINV_EXT)

         DO IK=1,ITERMAX

            LF_EXT   = MIN(1.0_JWRB, EXP(UCINV_EXT * RTAU) )
            TAUWX    = TAUWINDS(SDENSX_EXT*LF_EXT,CINV_EXT,DSII_EXT)
            TAUWY    = TAUWINDS(SDENSY_EXT*LF_EXT,CINV_EXT,DSII_EXT)
            TAU_WAV  = SQRT(TAUWX**2 + TAUWY**2)
            TAUX     = TAUVX + TAUWX
            TAUY     = TAUVY + TAUWY
            TAU      = SQRT(TAUX**2 + TAUY**2)
            ERR      = (TAU-TAU_TOT) / TAU_TOT
            SIGN_OLD = SIGN_NEW
            SIGN_NEW = INT(SIGN(1.0_JWRB, ERR))

!        --- Slow down DRTAU when overshot. -------------------------- /

            IF (SIGN_NEW .NE. SIGN_OLD) OVERSHOT = .TRUE.
            IF (OVERSHOT) DRTAU = MAX(0.5_JWRB*(1.0_JWRB+DRTAU),1.00010_JWRB)

            RTAU = RTAU * (DRTAU**SIGN_NEW)            

            IF (ABS(ERR) .LT. 1.54E-4_JWRB) EXIT
            
         END DO

      END IF

      LFACT(1:NFRE) = LF_EXT(1:NFRE)
      LREDUCE = ANY(LFACT(1:NFRE) .LT. 1.0_JWRB)

      IF (LHOOK) CALL DR_HOOK('LFACTOR',1,ZHOOK_HANDLE)

      END SUBROUTINE LFACTOR
