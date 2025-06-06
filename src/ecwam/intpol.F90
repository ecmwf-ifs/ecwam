! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE INTPOL (KIJS, KIJL, FLR, FLA, WAVNUM, UCUR, VCUR, IRA)

! ----------------------------------------------------------------------

!**** *INTPOL* - TRANSFORMATION OF SPECTRA FROM SIGMA TO OMEGA.

!     S.D.HASSELMANN      MPI            1.1.91
!     H. GUNTHER          GKSS/ECMWF     1.2.91  MODIFIED FOR CYCLE_4

!*    PURPOSE.
!     --------

!       TRANSFORMATION OF A MOVING COORDINATE SYSTEM TO AN ABSOLUTE
!       COORDINATE SYSTEM.

!**   INTERFACE.
!     ----------

!       *CALL* *INTPOL (KIJS, KIJL, FLR, FLA, WAVNUM, UCUR, VCUR, IRA)*
!         *KIJS*  - INDEX OF FIRST GRIDPOINT.
!         *KIJL*  - INDEX OF LAST GRIDPOINT.
!         *FLR*   - SPECTRA (INPUT)
!         *FLA*   - SPECTRA (OUTPUT)
!         *WAVNUM*- WAVE NUMBER.
!         *UCUR*  - U-COMPONENT OF THE SURFACE CURRENT
!         *VCUR*  - V-COMPONENT OF THE SURFACE CURRENT
!         *IRA*   - = 1 TRANSFORMATION FROM MOVING TO ABSOLUTE COORD.
!                   =-1 TRANSFORMATION FROM ABSOLUTE TO MOVING COORD.

!     METHOD.
!     -------

!       SCATTERING TO NEIGHBOURING POINT.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCURR  , ONLY : CURRENT_MAX
      USE YOWFRED  , ONLY : FR       ,DFIM     ,COSTH    ,SINTH     ,   &
     &              DELTH  ,FRATIO   ,FLOGSPRDM1, FR5
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,EPSMIN
      USE YOWTEST  , ONLY : IU06

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(IN) :: FLR
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(OUT) :: FLA
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: WAVNUM 
      REAL(KIND=JWRB), DIMENSION (KIJL), INTENT(IN) :: UCUR, VCUR
      INTEGER(KIND=JWIM), INTENT(IN) :: IRA


      INTEGER(KIND=JWIM) :: IJ, M, K
      INTEGER(KIND=JWIM) :: NFRE_MAX, NEWM, NEWM1, KH
      INTEGER(KIND=JWIM), DIMENSION(KIJL) :: NEWF, NEWFLA, KNEW

      REAL(KIND=JWRB) :: FRE0, CDF, ZPI2GM, COEF, FMAX, FREQ, DFREQTH, FR5OFREQ5
      REAL(KIND=JWRB) :: FNEW, GWH
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NFRE) :: DFTH
      REAL(KIND=JWRB), DIMENSION(KIJL) :: OLDFL, WAVN
      REAL(KIND=JWRB), DIMENSION(KIJL) :: FNEF, GWP, GWM

      LOGICAL, DIMENSION(KIJL) :: LICE2SEA
! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('INTPOL',0,ZHOOK_HANDLE)

!*    0. INITIAL OUTPUT ARRAY WITH ZERO.
!        -------------------------------

      FRE0 = FRATIO-1.0_JWRB
      ZPI2GM = ZPI**2/G

      COEF = IRA/ZPI

!!    MAXIMUM EXTEND OF THE HIGH FRQUENCY TAIL THAT CAN BE REMAPPED ONTO THE FREQUENCY DISCRETISATION
      FMAX = FR(NFRE) + (ZPI/G)*FR(NFRE)**2 * CURRENT_MAX
      NFRE_MAX=FLOOR(LOG10(FMAX/FR(1))*FLOGSPRDM1)+1

!!!??? I believe that in order to be energy conserving, since
!!! DFIM(1) is only the part above FR(1)
!!! purpose of trapezoidal rule, one has to take DFTH(1) which 
!!! define on both side of FR(1). Similarly for FR(NFRE)
!!! Redefine the centered frequency-direction intervals
      CDF= 0.5_JWRB*(FRATIO-1.0_JWRB/FRATIO)*DELTH
      DO M=1,NFRE
        DFTH(M)= FR(M)*CDF
      ENDDO

      IF (ABS(IRA) /= 1) THEN
        WRITE (IU06,*) '**************************************'
        WRITE (IU06,*) '* SUBROUTINE INTPOL:                  '
        WRITE (IU06,*) '* IRA MUST BE = 1 OR -1               '
        WRITE (IU06,*) '* IRA = ',IRA 
        WRITE (IU06,*) '*                                    *'
        WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.  *'
        WRITE (IU06,*) '*                                    *'
        WRITE (IU06,*) '**************************************'
        CALL ABORT1
      ENDIF

      DO IJ = KIJS, KIJL
        LICE2SEA(IJ)=.TRUE.
      ENDDO

      DO K=1,NANG
        DO M=1,NFRE
          DO IJ = KIJS, KIJL
             IF (FLR(IJ,K,M) > EPSMIN) LICE2SEA(IJ) = .FALSE. 
          ENDDO
        ENDDO
      ENDDO


      DO M = 1, NFRE
        DO K = 1, NANG
          DO IJ = KIJS, KIJL
            FLA(IJ,K,M) = 0.0_JWRB
          ENDDO
        ENDDO
      ENDDO

!*    1. LOOP OVER FREQUENCIES.
!        ----------------------

      DO M = 1, NFRE_MAX

        IF ( M <= NFRE ) THEN
          FREQ = FR(M)
          DFREQTH = DFTH(M)
          DO IJ = KIJS, KIJL
            WAVN(IJ) = WAVNUM(IJ,M)
          ENDDO
        ELSE
          FREQ = FR(NFRE)*FRATIO**(M-NFRE)
          DFREQTH = FREQ*CDF
          DO IJ = KIJS, KIJL
            WAVN(IJ) = ZPI2GM * FREQ**2
          ENDDO
        ENDIF

        FR5OFREQ5 = FR5(NFRE)/(FREQ**5)

!*    1.3 LOOP OVER DIRECTONS.
!         --------------------

        DO K = 1, NANG

!*    1.3.1 NEW FREQUENCY AND DIRECTION AT ALL GRIDPOINTS.
!           ----------------------------------------------

          DO IJ = KIJS, KIJL
            FNEF(IJ) = FREQ + COEF*WAVN(IJ)*(COSTH(K)*VCUR(IJ) + SINTH(K)*UCUR(IJ))
            IF (FNEF(IJ) > 0.0_JWRB) THEN
              KNEW(IJ) = K
            ELSE
              KNEW(IJ) = MOD(K+NANG/2-1,NANG) + 1
              FNEF(IJ) = -FNEF(IJ)
            ENDIF
          ENDDO

!*    1.3.2 NEW FREQUENCY BIN NUMBER AT ALL GRIDPOINTS.
!           -------------------------------------------

          DO IJ = KIJS, KIJL
            IF (FNEF(IJ) <= FR(1)/FRATIO) THEN
              NEWF(IJ)= -1
            ELSE
              NEWF(IJ)=FLOOR(LOG10(FNEF(IJ)/FR(1))*FLOGSPRDM1)+1
            ENDIF
          ENDDO

!*    1.3.3 INTERPOLATED ENERGY DENSITIES AT ALL GRIDPOINTS.
!           ------------------------------------------------
          IF ( M <= NFRE ) THEN
            DO IJ = KIJS, KIJL
              IF (LICE2SEA(IJ)) THEN
                OLDFL(IJ)=0.0_JWRB
              ELSE
                OLDFL(IJ)=FLR(IJ,K,M)
              ENDIF
            ENDDO
          ELSE
            DO IJ = KIJS, KIJL
              IF (LICE2SEA(IJ)) THEN
                OLDFL(IJ)=0.0_JWRB
              ELSE
                !! Extend with a f**-5 tail
                OLDFL(IJ)=FLR(IJ,K,NFRE)*FR5OFREQ5
              ENDIF
            ENDDO
          ENDIF

          DO IJ = KIJS, KIJL
            FNEW = FNEF(IJ)
            NEWM = NEWF(IJ)
            IF (NEWM < NFRE .AND. NEWM >= 1) THEN
              !!!! DFTH is only defined between 1 and NFRE
              NEWM1 = NEWM + 1
              GWH = DFREQTH/(FR(NEWM1)-FR(NEWM)) *OLDFL(IJ)
              GWM(IJ) = GWH*(FR(NEWM1)-FNEW)/DFTH(NEWM)
              GWP(IJ) = GWH*(FNEW-FR(NEWM))/DFTH(NEWM1)
              NEWFLA(IJ) = NEWM1
            ELSEIF (NEWM == 0) THEN
              GWH = FRATIO*DFREQTH/(FRE0*FR(1)) * OLDFL(IJ)
              GWP(IJ) = GWH*(FNEW-FR(1)/FRATIO)/DFTH(1)
              NEWF (IJ) = -1
              NEWFLA(IJ) = 1
            ELSEIF (NEWM == NFRE) THEN
              GWH = DFREQTH/(FRE0*FR(NFRE)) * OLDFL(IJ)
              GWM(IJ) = GWH*(FRATIO*FR(NFRE)-FNEW)/DFTH(NFRE)
              NEWFLA(IJ) = -1
            ELSE
              NEWF (IJ) = -1
              NEWFLA(IJ) = -1
            ENDIF
          ENDDO

!*    1.3.4 NEW SPECTRUM AT ALL GRIDPOINTS.
!           -------------------------------

          DO IJ = KIJS, KIJL
            NEWM  = NEWF (IJ)
            NEWM1 = NEWFLA(IJ)
            KH = KNEW(IJ)
            IF (NEWM /= -1) FLA(IJ,KH,NEWM ) = FLA(IJ,KH,NEWM ) + GWM(IJ)
            IF (NEWM1 /= -1)FLA(IJ,KH,NEWM1) = FLA(IJ,KH,NEWM1) + GWP(IJ)
          ENDDO

!*    BRANCH BACK TO 1.3 FOR NEXT DIRECTION.

        ENDDO

!*    BRANCH BACK TO 1. FOR NEXT FREQUENCY.

      ENDDO

      DO M = 1, NFRE
        DO K = 1, NANG
          DO IJ = KIJS, KIJL
            FLA(IJ,K,M) = MAX(FLA(IJ,K,M),EPSMIN)
          ENDDO
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('INTPOL',1,ZHOOK_HANDLE)

      END SUBROUTINE INTPOL
