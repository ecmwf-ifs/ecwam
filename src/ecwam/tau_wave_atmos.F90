! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

        SUBROUTINE TAU_WAVE_ATMOS(S, CINV, SIG, DSII, TAUNWX, TAUNWY )

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

!     *CALL* *TAU_WAVE_ATMOS(S, CINV, SIG, DSII, TAUNWX, TAUNWY )

!            *S*  - NEG. WIND INPUT ENERGY DENSITY SPECTRUM.
!          *CINV* - INVERSE PHASE SPEED CALC. IN INPUT ROUTINE
!          *SIG*  - FREQ (RAD)
!          *DSII* - ZPI*DF
!          *TAUNWX, TAUNWY* - NEGATIVE WAVE NORMAL STRESS COMPONENTS

!     EXTERNALS.
!     ----------
!     TAUWINDS
!     IRANGE

!     ORIGIN.
!     ----------
!     Adapted from Babanin Young Donelan & Banner (BYDB) physics 
!     as implemented as ST6 in WAVEWATCH-III 
!     WW3 module:       W3SRC6MD    
!     WW3 subroutine:   TAU_WAVE_ATMOS
!     Implementation into ECWAM DECEMBER 2021 by J. Kousal 

! ----------------------------------------------------------------------
        USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

        USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,TAUWSHELTER, XKAPPA, RNU      ,RNUM
        USE YOWFRED  , ONLY : FR       ,TH       ,DFIM     ,COSTH  ,SINTH,&
        &                      FRATIO   ,DELTH
        USE YOWMPP   , ONLY : NINF     ,NSUP
        USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO
        USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER  ,EPSMIN
        USE YOWSHAL  , ONLY : TFAK     ,INDEP
        USE YOWSTAT  , ONLY : ISHALLO
        USE YOWTABL  , ONLY : IAB      ,SWELLFT
        USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

        IMPLICIT NONE
#include "irange.intfb.h"
#include "tauwinds.intfb.h"

        REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(IN)  :: S
        REAL(KIND=JWRB), DIMENSION(NFRE),      INTENT(IN)  :: CINV, SIG, DSII
        REAL(KIND=JWRB),                       INTENT(OUT) :: TAUNWX, TAUNWY 

        REAL(KIND=JWRB), PARAMETER :: FRQMAX  = 10.0_JWRB ! Upper freq. limit to extrapolate to.

        REAL(KIND=JWRB), DIMENSION(NANG*NFRE) :: ECOS2, ESIN2
        REAL, ALLOCATABLE :: IK10Hz(:), SIG10Hz(:), CINV10Hz(:)
        REAL, ALLOCATABLE :: SDENSX10Hz(:), SDENSY10Hz(:)
        REAL, ALLOCATABLE :: DSII10Hz(:), UCINV10Hz(:)

        INTEGER(KIND=JWIM) :: IK, ITH, NK10Hz
        INTEGER(KIND=JWIM) :: NK, NTH, NSPEC !num. of freqs, dirs, spec. bins
        INTEGER(KIND=JWIM), DIMENSION(NANG) :: ITHN
        INTEGER(KIND=JWIM), DIMENSION(NFRE) :: IKN

        REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

        IF (LHOOK) CALL DR_HOOK('TAU_WAVE_ATMOS',0,ZHOOK_HANDLE)

        NTH   = NANG  ! NUMBER OF DIRS , SAME AS KL
        NK    = NFRE  ! NUMBER OF FREQS, SAME AS ML
        NSPEC = NK * NTH   ! NUMBER OF SPECTRAL BINS

!/ 0) --- Find the number of frequencies required to extend arrays
!/        up to f=10Hz and allocate arrays --------------------------- /
        NK10Hz = CEILING(ALOG(FRQMAX/(SIG(1)/ZPI))/ALOG(FRATIO))+1
        NK10Hz = MAX(NK,NK10Hz)
!
        ALLOCATE(IK10Hz(NK10Hz))
        IK10Hz = REAL( IRANGE(1,NK10Hz,1) )
!
        ALLOCATE(SIG10Hz(NK10Hz))
        ALLOCATE(CINV10Hz(NK10Hz))
        ALLOCATE(DSII10Hz(NK10Hz))
        ALLOCATE(SDENSX10Hz(NK10Hz))
        ALLOCATE(SDENSY10Hz(NK10Hz))
        ALLOCATE(UCINV10Hz(NK10Hz))
!
        ITHN   = IRANGE(1,NTH,1)    ! Index vector 1:NTH
        DO IK = 1, NK
                ECOS2 (ITHN+(IK-1)*NTH) = COSTH
                ESIN2 (ITHN+(IK-1)*NTH) = SINTH
        END DO

!
!/ 1) --- Either extrapolate arrays up to 10Hz or use discrete spectral
!         grid per se. Limit the constraint to the positive part of the
!         wind input only. ---------------------------------------------- /
        IF (NK .LT. NK10Hz) THEN
                SDENSX10Hz(1:NK)        = SUM(ABS(MIN(0.0_JWRB,S))*&
                                        &     RESHAPE(ECOS2,(/NTH,NK/)),1) * DELTH
                SDENSY10Hz(1:NK)        = SUM(ABS(MIN(0.0_JWRB,S))*&
                                        &     RESHAPE(ESIN2,(/NTH,NK/)),1) * DELTH
                SIG10Hz                 = SIG(1)*FRATIO**(IK10Hz-1.0_JWRB)
                CINV10Hz(1:NK)          = CINV
                CINV10Hz(NK+1:NK10Hz)   = SIG10Hz(NK+1:NK10Hz)*0.101978_JWRB
                DSII10Hz                = 0.5_JWRB * SIG10Hz * (FRATIO-1.0_JWRB/FRATIO)
!        The first and last frequency bin:
                DSII10Hz(1)             = 0.5_JWRB * SIG10Hz(1) * (FRATIO-1.0_JWRB)
                DSII10Hz(NK10Hz)        = 0.5_JWRB * SIG10Hz(NK10Hz) * &
                                        &            (FRATIO-1.0_JWRB) / FRATIO
!
!        --- Spectral slope for S_IN(F) is proportional to F**(-2) ------ /
                SDENSX10Hz(NK+1:NK10Hz) = SDENSX10Hz(NK) * & 
                                        & (SIG10Hz(NK)/SIG10Hz(NK+1:NK10Hz))**2
                SDENSY10hz(NK+1:NK10Hz) = SDENSY10Hz(NK) * &
                                        & (SIG10Hz(NK)/SIG10Hz(NK+1:NK10Hz))**2
        ELSE
                SIG10Hz          = SIG
                CINV10Hz         = CINV
                DSII10Hz         = DSII
                SDENSX10Hz(1:NK) = SUM(ABS(MIN(0.0_JWRB,S))*&
                                & RESHAPE(ECOS2,(/NTH,NK/)),1) * DELTH
                SDENSY10Hz(1:NK) = SUM(ABS(MIN(0.0_JWRB,S))*&
                                & RESHAPE(ESIN2,(/NTH,NK/)),1) * DELTH
        END IF
!
!/ 2) --- Stress calculation ----------------------------------------- /
!     --- The wave supported stress (waves to atmosphere) ------------ /
        TAUNWX = TAUWINDS(SDENSX10Hz,CINV10Hz,DSII10Hz)   ! x-component
        TAUNWY = TAUWINDS(SDENSY10Hz,CINV10Hz,DSII10Hz)   ! y-component


        IF (LHOOK) CALL DR_HOOK('TAU_WAVE_ATMOS',1,ZHOOK_HANDLE)

        END SUBROUTINE TAU_WAVE_ATMOS
