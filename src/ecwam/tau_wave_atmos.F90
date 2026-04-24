! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

        SUBROUTINE TAU_WAVE_ATMOS(S, CINV, TAUNWX, TAUNWY )

! ----------------------------------------------------------------------
!
!  1. Purpose :
!
!     Calculated the stress for the negative part of the input term,
!     that is the stress from the waves to the atmosphere. Relevant
!     in the case of opposing winds.
!
!  2. Method :
!     1) If required, extend resolved part of the spectrum to 10Hz using
!        an approximation for the spectral slope at the high frequency
!        limit: Sin(F) prop. F**(-2) and for E(F) prop. F**(-5).
!     2) Calculate stresses:
!        stress components (x,y):      /10Hz
!            TAUNW_X,Y = GRAV * DWAT * | [SinX,Y(F)]/C(F) dF
!                                      /
!
! ----------------------------------------------------------------------------
!
!**   INTERFACE.
!     ----------

!     *CALL* *TAU_WAVE_ATMOS(S, CINV, TAUNWX, TAUNWY )

!            *S*  - NEG. WIND INPUT ENERGY DENSITY SPECTRUM.
!          *CINV* - INVERSE PHASE SPEED CALC. IN INPUT ROUTINE
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
        &                      FRATIO  ,DELTH    ,SIG      , DSII
        USE YOWMPP   , ONLY : NINF     ,NSUP
        USE YOWPARAM , ONLY : NANG     ,NFRE
        USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER  ,EPSMIN, GM1
        USE YOWPHYS  , ONLY : FRQMAX
        USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

        IMPLICIT NONE
#include "tauwinds.intfb.h"

        REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(IN)  :: S ! Sin(sigma) in [m2/rad-Hz]
        REAL(KIND=JWRB), DIMENSION(NFRE),      INTENT(IN)  :: CINV
        REAL(KIND=JWRB),                       INTENT(OUT) :: TAUNWX, TAUNWY 


        REAL(KIND=JWRB), DIMENSION(NANG,NFRE) :: SX, SY

        REAL(KIND=JWRB), DIMENSION(NFRE)  :: SDENSX_LF, SDENSY_LF
        REAL(KIND=JWRB), DIMENSION(NFRE)  :: ZA_SX, ZA_SY

        REAL(KIND=JWRB)               :: SDENSX_HF, SDENSY_HF
        REAL(KIND=JWRB)               :: TAUNWX_LF, TAUNWY_LF
        REAL(KIND=JWRB)               :: TAUNWX_HF, TAUNWY_HF
        
        REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('TAU_WAVE_ATMOS',0,ZHOOK_HANDLE)

      SX = ABS(MIN(0.0_JWRB,S)) * SPREAD(COSTH,2,NFRE)
      SY = ABS(MIN(0.0_JWRB,S)) * SPREAD(SINTH,2,NFRE)

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
      
      TAUNWX_LF = TAUWINDS(SDENSX_LF,CINV,DSII)   ! x-component
      TAUNWY_LF = TAUWINDS(SDENSY_LF,CINV,DSII)   ! y-component

      !/ 2) --- high frequency contributions to the integral --------------------- /
      !         -- Assume spectral slope for S_IN(F) is proportional to F**(-2), then 
      !            integral collapses into easy analytic solution
      !
      !          
      !   Th=2pi,ω=inf  
      !     / /
      !     | | S(f,Th)/c df dTh = (LOG(inf) - LOG(SIG(NFRE))) * SIG(NFRE)**2 * DELTH * SUM(S(:,NFRE)) * ZPI * GM1
      !    / /
      !  Th=0,ω=SIG(NFRE)
      !
      !  Instead, use FRQMAX extension (e.g. 10Hz)          
      !
      !   Th=2pi,ω=ZPI*FRQMAX  
      !     / /
      !     | | S(f,Th)/c df dTh = (LOG(ZPI*FRQMAX) - LOG(SIG(NFRE))) * SIG(NFRE)**2 * DELTH * SUM(S(:,NFRE)) * ZPI * GM1
      !    / /
      !  Th=0,ω=SIG(NFRE)
      !
      !
      ! Determine value of spectrum at NFRE (i.e. at highest frequency). 
      !  - Note, direction dimension must remain

      ZA_SX       = SX(:,NFRE)
      ZA_SY       = SY(:,NFRE)

      SDENSX_HF   = (LOG(ZPI*FRQMAX) - LOG(SIG(NFRE))) * SIG(NFRE)**2 * DELTH * SUM(ZA_SX) * ZPI * GM1
      SDENSY_HF   = (LOG(ZPI*FRQMAX) - LOG(SIG(NFRE))) * SIG(NFRE)**2 * DELTH * SUM(ZA_SY) * ZPI * GM1

      TAUNWX_HF   = G * ROWATER * ( SDENSX_HF )
      TAUNWY_HF   = G * ROWATER * ( SDENSY_HF )

      !/ 3) --- summate low + high frequency contributions to the integral ------- /
      
      TAUNWX = TAUNWX_LF + TAUNWX_HF
      TAUNWY = TAUNWY_LF + TAUNWY_HF

      IF (LHOOK) CALL DR_HOOK('TAU_WAVE_ATMOS',1,ZHOOK_HANDLE)

      END SUBROUTINE TAU_WAVE_ATMOS
