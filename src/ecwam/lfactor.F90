! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

      SUBROUTINE LFACTOR(S, CINV, U10, USTAR, USDIR, ROAIRN, SIG, DSII, &
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
!          *SIG*   - FREQ (RAD)
!          *DSII*  - ZPI*DF
!          *LFACT* - CORRECTION FACTOR
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
!     WW3 subroutine:   LFACTOR
!     Implementation into ECWAM DECEMBER 2021 by J. Kousal 

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,TH       ,DFIM     ,COSTH  ,SINTH,&
&                      FRATIO   ,DELTH    ,FRIC
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER  ,EPSMIN
      USE YOWPHYS  , ONLY : BETAMAX  ,ZALP     ,TAUWSHELTER, XKAPPA, RNU ,RNUM, CDFAC
      USE YOWSTAT  , ONLY : ISHALLO
      USE YOWTABL  , ONLY : IAB      ,SWELLFT
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "irange.intfb.h"
#include "tauwinds.intfb.h"
 
      REAL(KIND=JWRB), DIMENSION(NANG,NFRE), INTENT(IN)  :: S
      REAL(KIND=JWRB), DIMENSION(NFRE),      INTENT(IN)  :: CINV, SIG, DSII
      REAL(KIND=JWRB),                       INTENT(IN)  :: U10, USTAR, USDIR, ROAIRN

      REAL(KIND=JWRB), DIMENSION(NFRE),      INTENT(OUT) :: LFACT
      REAL(KIND=JWRB),                       INTENT(OUT) :: TAUWX, TAUWY, TAU

      REAL(KIND=JWRB), PARAMETER :: FRQMAX  = 10.0_JWRB ! Upper freq. limit to extrap. to 
      INTEGER(KIND=JWIM), PARAMETER :: ITERMAX = 80 ! Max. no. iterations
                                                ! to find numerical LFACT soln

      REAL(KIND=JWRB), DIMENSION(NANG*NFRE) :: ECOS2, ESIN2
      REAL, ALLOCATABLE :: IK10Hz(:), LF10Hz(:), SIG10Hz(:), CINV10Hz(:)
      REAL, ALLOCATABLE :: SDENS10Hz(:), SDENSX10Hz(:), SDENSY10Hz(:)
      REAL, ALLOCATABLE :: DSII10Hz(:), UCINV10Hz(:)

      REAL(KIND=JWRB)      :: TAU_TOT, TAU_VIS, TAU_WAV
      REAL(KIND=JWRB)      :: TAUVX, TAUVY, TAUX, TAUY   
      REAL(KIND=JWRB)      :: TAU_NND, TAU_INIT(2)
      REAL(KIND=JWRB)      :: UPROXY, RTAU, DRTAU, ERR   
      LOGICAL              :: OVERSHOT
      CHARACTER(LEN=23)    :: IDTIME

      INTEGER(KIND=JWIM) :: IK, ITH, NK10Hz, M, SIGN_NEW, SIGN_OLD
      INTEGER(KIND=JWIM) :: NK, NTH, NSPEC !num. of freqs, dirs, spec. bins
      INTEGER(KIND=JWIM), DIMENSION(NANG) :: ITHN
      INTEGER(KIND=JWIM), DIMENSION(NFRE) :: IKN

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('LFACTOR',0,ZHOOK_HANDLE)

      NTH   = NANG  ! NUMBER OF DIRS , SAME AS KL
      NK    = NFRE  ! NUMBER OF FREQS, SAME AS ML
      NSPEC = NK * NTH   ! NUMBER OF SPECTRAL BINS

!/ 0) --- Find the number of frequencies required to extend arrays
!/        up to f=10Hz and allocate arrays --------------------------- /
      NK10Hz = CEILING(LOG(FRQMAX/(SIG(1)/ZPI))/LOG(FRATIO))+1
      NK10Hz = MAX(NK,NK10Hz)
!
      ALLOCATE(IK10Hz(NK10Hz))
      IK10Hz = REAL( IRANGE(1,NK10Hz,1) )
!
      ALLOCATE(SIG10Hz(NK10Hz))
      ALLOCATE(CINV10Hz(NK10Hz))
      ALLOCATE(DSII10Hz(NK10Hz))
      ALLOCATE(LF10Hz(NK10Hz))
      ALLOCATE(SDENS10Hz(NK10Hz))
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
            SDENS10Hz(1:NK)         = SUM(S,1) * DELTH
            SDENSX10Hz(1:NK)        = SUM(MAX(0.0_JWRB,S)*RESHAPE(ECOS2,(/NTH,NK/)),1) * DELTH
            SDENSY10Hz(1:NK)        = SUM(MAX(0.0_JWRB,S)*RESHAPE(ESIN2,(/NTH,NK/)),1) * DELTH
            SIG10Hz                 = SIG(1)*FRATIO**(IK10Hz-1.0_JWRB)
            CINV10Hz(1:NK)          = CINV
            CINV10Hz(NK+1:NK10Hz)   = SIG10Hz(NK+1:NK10Hz)*0.101978_JWRB ! 1/c=σ/g
            DSII10Hz                = 0.5_JWRB * SIG10Hz * (FRATIO-1.0_JWRB/FRATIO)
!        The first and last frequency bin:
            DSII10Hz(1)             = 0.5_JWRB * SIG10Hz(1) * (FRATIO-1.0_JWRB)
            DSII10Hz(NK10Hz)        = 0.5_JWRB * SIG10Hz(NK10Hz) * (FRATIO-1.0_JWRB) / FRATIO
!
!        --- Spectral slope for S_IN(F) is proportional to F**(-2) ------ /
            SDENS10Hz(NK+1:NK10Hz)  = SDENS10Hz(NK)  * (SIG10Hz(NK)/SIG10Hz(NK+1:NK10Hz))**2
            SDENSX10Hz(NK+1:NK10Hz) = SDENSX10Hz(NK) * (SIG10Hz(NK)/SIG10Hz(NK+1:NK10Hz))**2
            SDENSY10hz(NK+1:NK10Hz) = SDENSY10Hz(NK) * (SIG10Hz(NK)/SIG10Hz(NK+1:NK10Hz))**2
      ELSE
            SIG10Hz          = SIG
            CINV10Hz         = CINV
            DSII10Hz         = DSII
            SDENS10Hz(1:NK)  = SUM(S,1) * DELTH
            SDENSX10Hz(1:NK) = SUM(MAX(0.0_JWRB,S)*RESHAPE(ECOS2,(/NTH,NK/)),1) * DELTH
            SDENSY10Hz(1:NK) = SUM(MAX(0.0_JWRB,S)*RESHAPE(ESIN2,(/NTH,NK/)),1) * DELTH
      END IF
!
!/ 2) --- Stress calculation ----------------------------------------- /
!     --- The total stress ------------------------------------------- /
      TAU_TOT  = USTAR**2 * ROAIRN
!
!     --- The viscous stress and check that it does not exceed
!         the total stress. ------------------------------------------ /
      TAU_VIS  = MAX(0.0_JWRB, -5.0E-5_JWRB*U10 + 1.1E-3_JWRB) * U10**2 * ROAIRN
!       TAU_VIS  = MIN(0.9 * TAU_TOT, TAU_VIS)
      TAU_VIS  = MIN(0.95_JWRB * TAU_TOT, TAU_VIS)
!
      TAUVX    = TAU_VIS * COS(USDIR)
      TAUVY    = TAU_VIS * SIN(USDIR)
!
!     --- The wave supported stress. --------------------------------- /
      TAUWX    = TAUWINDS(SDENSX10Hz,CINV10Hz,DSII10Hz)   ! normal stress (x-component)
      TAUWY    = TAUWINDS(SDENSY10Hz,CINV10Hz,DSII10Hz)   ! normal stress (y-component)
      TAU_NND  = TAUWINDS(SDENS10Hz, CINV10Hz,DSII10Hz)   ! normal stress (non-directional)
      TAU_WAV  = SQRT(TAUWX**2 + TAUWY**2)        ! normal stress (magnitude)
      TAU_INIT = (/TAUWX,TAUWY/)                  ! unadjusted normal stress components
!
      TAUX     = TAUVX + TAUWX                        ! total stress (x-component)
      TAUY     = TAUVY + TAUWY                        ! total stress (y-component)
      TAU      = SQRT(TAUX**2  + TAUY**2)                 ! total stress (magnitude)
      ERR      = (TAU-TAU_TOT)/TAU_TOT                    ! initial error
!
!/ 3) --- Find reduced Sin(f) = L(f)*Sin(f) to satisfy our constraint
!/        TAU <= TAU_TOT --------------------------------------------- /
      !CALL STME21 ( TIME , IDTIME )
      LF10Hz = 1.0_JWRB
      IK = 0
!
      IF (TAU .GT. TAU_TOT) THEN
      OVERSHOT    = .FALSE.
      RTAU        = ERR / 90.0_JWRB
      DRTAU       = 2.0_JWRB
      SIGN_NEW    = INT(SIGN(1.0_JWRB,ERR))

      UPROXY      = FRIC * CDFAC * USTAR
      UCINV10Hz   = 1.0_JWRB - (UPROXY * CINV10Hz)
!
!/T6          WRITE (NDST,270) IDTIME, U10
!/T6          WRITE (NDST,271)
      DO IK=1,ITERMAX
            LF10Hz   = MIN(1.0_JWRB, EXP(UCINV10Hz * RTAU) )
!
            TAU_NND  = TAUWINDS(SDENS10Hz *LF10Hz,CINV10Hz,DSII10Hz)
            TAUWX    = TAUWINDS(SDENSX10Hz*LF10Hz,CINV10Hz,DSII10Hz)
            TAUWY    = TAUWINDS(SDENSY10Hz*LF10Hz,CINV10Hz,DSII10Hz)
            TAU_WAV  = SQRT(TAUWX**2 + TAUWY**2)
!
            TAUX     = TAUVX + TAUWX
            TAUY     = TAUVY + TAUWY
            TAU      = SQRT(TAUX**2 + TAUY**2)
            ERR      = (TAU-TAU_TOT) / TAU_TOT
!
            SIGN_OLD = SIGN_NEW
            SIGN_NEW = INT(SIGN(1.0_JWRB, ERR))
!/T6             WRITE (NDST,272) IK, RTAU, DRTAU, TAU, TAU_TOT, ERR, &
!/T6                              TAUWX, TAUWY, TAUVX, TAUVY, TAU_NND
!
!        --- Slow down DRTAU when overshot. -------------------------- /
            IF (SIGN_NEW .NE. SIGN_OLD) OVERSHOT = .TRUE.
            IF (OVERSHOT) DRTAU = MAX(0.5_JWRB*(1.0_JWRB+DRTAU),1.00010_JWRB)
!
            RTAU = RTAU * (DRTAU**SIGN_NEW)
!
            IF (ABS(ERR) .LT. 1.54E-4_JWRB) EXIT
      END DO
!
!           IF (IK .GE. ITERMAX) WRITE (NDST,280) IDTIME(1:19), U10, TAU, &
!                          TAU_TOT, ERR, TAUWX, TAUWY, TAUVX, TAUVY,TAU_NND
      END IF
!
      LFACT(1:NK) = LF10Hz(1:NK)
!
!/T6      WRITE (NDST,273) 'Sin ', IDTIME(1:19), SDENS10Hz*TPI
!/T6      WRITE (NDST,273) 'SinR', IDTIME(1:19), SDENS10Hz*LF10Hz*TPI
!/T6      WRITE (NDST,274) 'Sin   ', SUM(SDENS10Hz(1:NK)*DSII)
!/T6      WRITE (NDST,274) 'SinR  ', SUM(SDENS10Hz(1:NK)*LF10Hz(1:NK)*DSII)
!/T6      WRITE (NDST,274) 'SinR/C', TAUWINDS(SDENS10Hz(1:NK)*LFACT,CINV,DSII)
!
!/T6      270 FORMAT (' TEST W3SIN6 : LFACTOR SUBROUTINE CALCULATING FOR ', &
!/T6                  A,'  U10=',F5.1                                       )
!/T6      271 FORMAT (' TEST W3SIN6 : IK    RTAU     DRTAU   TAU   TAU_TOT' &
!/T6                  '    ERR    TAUW_X TAUW_Y TAUV_X TAUV_Y  TAU1D'       )
!/T6      272 FORMAT (' TEST W3SIN6 : ',I2,2F9.5,2F8.5,E10.2,4F7.4,F7.3     )
!/T6      273 FORMAT (' TEST W3SIN6 : ',A,'(',A,'):', 70E11.3               )
!        274 FORMAT (' TEST W3SIN6 : Total ',A,' =', E13.5                  )
!        280 FORMAT (' WARNING LFACTOR (TIME,U10,TAU,TAU_TOT,ERR,TAUW_XY,'  &
!                  'TAUV_XY,TAU_SCALAR): ',A,F6.1,2F7.4,E10.3,4F7.4,F7.3  )
!
      DEALLOCATE(IK10Hz,SIG10Hz,CINV10Hz,DSII10Hz,LF10Hz)
      DEALLOCATE(SDENS10Hz,SDENSX10Hz,SDENSY10Hz,UCINV10Hz)


      IF (LHOOK) CALL DR_HOOK('LFACTOR',1,ZHOOK_HANDLE)

      END SUBROUTINE LFACTOR
