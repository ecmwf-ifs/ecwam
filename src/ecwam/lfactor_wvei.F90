! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

      SUBROUTINE LFACTORXY(S, CINV, U10, USTAR, UPROXY, USDIR, ROAIRN, &
     &                   LFACT, TAUWX, TAUWY, TAU)
 
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
!     IRANGE

!     ORIGIN.
!     ----------
!     Adapted from Babanin Young Donelan & Banner (ZBRY) physics 
!     as implemented as ST6 in WAVEWATCH-III 
!     WW3 module:       W3SRC6MD    
!     WW3 subroutine:   LFACTOR
!     Implementation into ECWAM DECEMBER 2021 by J. Kousal 

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,TH       ,DFIM     ,COSTH  ,SINTH,&
&                           FRATIO   ,DELTH    ,FRIC, SIG,DSII   ,SIGM1,&
&                           DF       ,NFRE_EXT ,DSII_EXT ,SIG_EXT
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER  ,EPSMIN, GM1
      USE YOWPHYS  , ONLY : CDFAC    ,FRQMAX
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "irange.intfb.h"
#include "tauwinds.intfb.h"
#include "wvei.intfb.h"
#include "abort1.intfb.h"
 
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(IN)  :: S ! Sin(sigma) in omega
      REAL(KIND=JWRB), DIMENSION(NFRE),      INTENT(IN)  :: CINV
      REAL(KIND=JWRB),                       INTENT(IN)  :: U10, USTAR, UPROXY, USDIR, ROAIRN

      REAL(KIND=JWRB), DIMENSION(NFRE),      INTENT(OUT) :: LFACT
      REAL(KIND=JWRB),                       INTENT(OUT) :: TAUWX, TAUWY, TAU

      INTEGER(KIND=JWIM), PARAMETER :: ITERMAX = 80 ! Max. no. iterations
                                                    ! to find numerical LFACT soln

      REAL(KIND=JWRB), DIMENSION(NANG*NFRE) :: ECOS2, ESIN2
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE) :: SX, SY
      REAL(KIND=JWRB), DIMENSION(NFRE)  :: SDENSX_LF, SDENSY_LF
      REAL(KIND=JWRB), DIMENSION(NFRE)  :: ZA_SX, ZA_SY
      REAL(KIND=JWRB), DIMENSION(NFRE)  :: UCINV
      REAL(KIND=JWRB), DIMENSION(NFRE)  :: LF

      REAL(KIND=JWRB)               :: SDENSX_HF, SDENSY_HF, SDENSX_MF_HF, SDENSY_MF_HF
      REAL(KIND=JWRB)               :: SDENSX_MF, SDENSY_MF
      REAL(KIND=JWRB)               :: TAUWX_LF, TAUWY_LF, TAUWX_MF_HF, TAUWY_MF_HF
      REAL(KIND=JWRB)               :: TAUWX_HF, TAUWY_HF  
      REAL(KIND=JWRB)               :: ZA_EXP, ZWVEI, FRQMID

      REAL(KIND=JWRB)      :: TAU_TOT, TAU_VIS, TAU_WAV
      REAL(KIND=JWRB)      :: TAUVX, TAUVY, TAUX, TAUY

      REAL(KIND=JWRB)      :: TAU_NND, TAU_INIT(2)

      REAL(KIND=JWRB)      :: RTAU, DRTAU, ERR
      LOGICAL              :: OVERSHOT
      LOGICAL              :: LLFRQHF

      INTEGER(KIND=JWIM) :: IK, ITH, M, SIGN_NEW, SIGN_OLD
      INTEGER(KIND=JWIM) :: NK, NTH, NSPEC !num. of freqs, dirs, spec. bins
      INTEGER(KIND=JWIM), DIMENSION(NANG) :: ITHN
      INTEGER(KIND=JWIM), DIMENSION(NFRE) :: IKN

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('LFACTOR',0,ZHOOK_HANDLE)

      NTH   = NANG  ! NUMBER OF DIRS , SAME AS KL
      NK    = NFRE  ! NUMBER OF FREQS, SAME AS ML
      NSPEC = NK * NTH   ! NUMBER OF SPECTRAL BINS

      FRQMID = MIN(FRQMAX,G/(ZPI*UPROXY))  ! dynamic cutoff frequency, capped at FRQMAX
      
      ! Determines if we need special treatment for some of the high frequencies to account for LFAC
      IF (FRQMID>FR(NFRE)) THEN
         LLFRQHF = .TRUE.
      ELSE
         LLFRQHF = .FALSE.
      END IF

      ITHN   = IRANGE(1,NTH,1)    ! Index vector 1:NTH
      DO IK = 1, NK
         ECOS2 (ITHN+(IK-1)*NTH) = COSTH
         ESIN2 (ITHN+(IK-1)*NTH) = SINTH
      END DO

      SX = MAX(0.0_JWRB,S)*RESHAPE(ECOS2,(/NTH,NK/))  
      SY = MAX(0.0_JWRB,S)*RESHAPE(ESIN2,(/NTH,NK/))  

!/ ----------------------------------------------------------------- /
!/ ----------------------------------------------------------------- /
!/ Part I) --- low/high frequency contributions of TAU ------------- /

      !/ 0) --- split integral into low/high frequency contributions ------------- /
      !
      !
      !     Th=2pi,f=inf                   Th=2pi,f=FR(NFRE)         Th=2pi,f=inf  
      !          / /                         / /                      / /
      !          | | S(f,Th)/c df dTh  =     | | S(f,Th)/c df dTh +   | | S(f,Th)/c df dTh
      !         / /                         / /                      / /
      !     Th=0,f=0                     Th=0,f=0                 Th=0,f=FR(NFRE)
      !
      !
      !                                =     LF_contribution      +  HF_contribution
      !
      !
      !/ 1) --- low frequency contributions to the integral ---------------------- /
      !         -- Direct summation over available freq. bins up to FR(NFRE)  

      SDENSX_LF = SUM(SX,1) * DELTH
      SDENSY_LF = SUM(SY,1) * DELTH  

      TAUWX_LF = TAUWINDS(SDENSX_LF,CINV,DSII)   ! x-component
      TAUWY_LF = TAUWINDS(SDENSY_LF,CINV,DSII)   ! y-component

      !/ 2) --- high frequency contributions to the integral --------------------- /
      !         -- Assume spectral slope for S_IN(F) is proportional to F**(-2), then 
      !            integral collapses into easy analytic solution
      !
      !          
      !   Th=2pi,ω=ZPI*FRQMID  
      !     / /
      !     | | S(f,Th)/c df dTh = (LOG(ZPI*FRQMID) - LOG(SIG(NFRE))) * SIG(NFRE)**2 * DELTH * SUM(S(:,NFRE)) * ZPI * GM1
      !    / /
      !  Th=0,ω=SIG(NFRE)
      !
      !
      ! Determine value of spectrum at NFRE (i.e. at highest frequency). 
      !  - Note, direction dimension must remain

      ZA_SX       = SX(:,NFRE)
      ZA_SY       = SY(:,NFRE)

      ! mid + high frequency contributions
      SDENSX_MF_HF   = ZPI * GM1 * SIG(NFRE)**2 * (LOG(ZPI*FRQMAX) - LOG(SIG(NFRE))) * DELTH * SUM(ZA_SX)
      SDENSY_MF_HF   = ZPI * GM1 * SIG(NFRE)**2 * (LOG(ZPI*FRQMAX) - LOG(SIG(NFRE))) * DELTH * SUM(ZA_SY)

      TAUWX_MF_HF   = G * ROWATER * ( SDENSX_MF_HF )
      TAUWY_MF_HF   = G * ROWATER * ( SDENSY_MF_HF )

      !/ 3) --- summate low + mid + high frequency contributions to the integral ------- /

      TAUWX = TAUWX_LF + TAUWX_MF_HF
      TAUWY = TAUWY_LF + TAUWY_MF_HF
      ! WRITE (*,*) '   '
      ! WRITE (*,*) '  TAUWX_LF    ',TAUWX_LF
      ! WRITE (*,*) '  TAUWX_MF_HF ',TAUWX_MF_HF

!/ ----------------------------------------------------------------- /
!/ ----------------------------------------------------------------- /
!/ Part II) --- USTAR based TAU calculation ------------- /      
!
!/ 2) --- Stress calculation ----------------------------------------- /
!     --- The total stress ------------------------------------------- /

      TAU_TOT  = USTAR**2 * ROAIRN

!     --- The viscous stress and check that it does not exceed
!         the total stress. ------------------------------------------ /

      TAU_VIS  = MAX(0.0_JWRB, -5.0E-5_JWRB*U10 + 1.1E-3_JWRB) * U10**2 * ROAIRN
!       TAU_VIS  = MIN(0.9 * TAU_TOT, TAU_VIS)
      TAU_VIS  = MIN(0.95_JWRB * TAU_TOT, TAU_VIS)

      TAUVX    = TAU_VIS * COS(USDIR)
      TAUVY    = TAU_VIS * SIN(USDIR)

!     --- The wave supported stress (using elements calculated in Part I). -- /
!
      TAU_WAV  = SQRT(TAUWX**2 + TAUWY**2)        ! normal stress (magnitude)
      TAU_INIT = (/TAUWX,TAUWY/)                  ! unadjusted normal stress components

      TAUX     = TAUVX + TAUWX                        ! total stress (x-component)
      TAUY     = TAUVY + TAUWY                        ! total stress (y-component)
      TAU      = SQRT(TAUX**2  + TAUY**2)                 ! total stress (magnitude)
      ERR      = (TAU-TAU_TOT)/TAU_TOT                    ! initial error

!/ 3) --- Find reduced Sin(f) = L(f)*Sin(f) to satisfy our constraint
!/        TAU <= TAU_TOT --------------------------------------------- /

      LF = 1.0_JWRB
      IK = 0

      IF (LLFRQHF) THEN
         ! mid frequency contributions (from FR(NFRE) to FRQMID)
         SDENSX_MF   = (LOG(ZPI*FRQMID) - LOG(SIG(NFRE))) * SIG(NFRE)**2 * DELTH * SUM(ZA_SX) * ZPI * GM1
         SDENSY_MF   = (LOG(ZPI*FRQMID) - LOG(SIG(NFRE))) * SIG(NFRE)**2 * DELTH * SUM(ZA_SY) * ZPI * GM1
      ELSE
         ! mid + high frequency contributions (from FR(NFRE) to FRQMAX)
         SDENSX_MF_HF   = (LOG(ZPI*FRQMAX) - LOG(SIG(NFRE))) * SIG(NFRE)**2 * DELTH * SUM(ZA_SX) * ZPI * GM1
         SDENSY_MF_HF   = (LOG(ZPI*FRQMAX) - LOG(SIG(NFRE))) * SIG(NFRE)**2 * DELTH * SUM(ZA_SY) * ZPI * GM1
      END IF
      

      IF (TAU .GT. TAU_TOT) THEN
         OVERSHOT    = .FALSE.
         RTAU        = ERR / 90.0_JWRB
         DRTAU       = 2.0_JWRB
         SIGN_NEW    = INT(SIGN(1.0_JWRB,ERR)) 

         UCINV       = 1.0_JWRB - (UPROXY * CINV)
         DO IK=1,ITERMAX
            
            ! LFAC capped at 1 (i.e. always a reduction of Sin)
            LF       = MIN(1.0_JWRB, EXP(UCINV * RTAU) )

            !/ 1) --- low frequency contributions as above ---------------------- /
            TAUWX_LF = TAUWINDS(SDENSX_LF * LF,CINV,DSII)   ! x-component
            TAUWY_LF = TAUWINDS(SDENSY_LF * LF,CINV,DSII)   ! y-component

            !/ 2) --- high frequency contributions as above, but modification to include LFAC ------ /
            !
            !          
            !   Th=2pi,ω=ZPI*FRQMAX  
            !     / /
            !     | | LFAC*S(f,Th)/c df dTh = (-EXP(RTAU)*EI(ZA_EXP*SIG(NFRE))) * SIG(NFRE)**2 * DELTH * SUM(S(:,NFRE)) * ZPI * GM1
            !    / /
            !  Th=0,ω=ZPI*FRQMID
            !
            !       where,  LFAC       =  EXP(RTAU*( 1 - (U * 1/c)))
            !               ZA_EXP     = -ZPI*RTAU*U/GRAV
            !
            !       then , use EI for exponential integral solution
            !

            IF (LLFRQHF) THEN
                  ZA_EXP      = -ZPI*RTAU*UPROXY*GM1
                  ZWVEI       = -EXP(RTAU)*WVEI(ZA_EXP*(ZPI*FRQMID))
                  SDENSX_HF   = ZWVEI * (ZPI*FRQMID)**2 * ZPI * GM1 * DELTH * SUM(ZA_SX) 
                  SDENSY_HF   = ZWVEI * (ZPI*FRQMID)**2 * ZPI * GM1 * DELTH * SUM(ZA_SY)
                  TAUWX_MF_HF = G * ROWATER * ( SDENSX_MF + SDENSX_HF)
                  TAUWY_MF_HF = G * ROWATER * ( SDENSY_MF + SDENSY_HF)
            ELSE
                  TAUWX_MF_HF = G * ROWATER * ( SDENSX_MF_HF )
                  TAUWY_MF_HF = G * ROWATER * ( SDENSY_MF_HF )
            END IF

            TAUWX    = TAUWX_LF + TAUWX_MF_HF
            TAUWY    = TAUWY_LF + TAUWY_MF_HF
            ! WRITE (*,*) '   '
            ! WRITE (*,*) '  TAUWX_LF    ',TAUWX_LF
            ! WRITE (*,*) '  TAUWX_MF_HF ',TAUWX_MF_HF

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

            ! CALL ABORT1

            IF (ABS(ERR) .LT. 1.54E-4_JWRB) EXIT

         END DO

      END IF

      LFACT = LF

      IF (LHOOK) CALL DR_HOOK('LFACTOR',1,ZHOOK_HANDLE)

      END SUBROUTINE LFACTORXY
