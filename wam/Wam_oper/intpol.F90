      SUBROUTINE INTPOL (GFL, FLA, IJS, IJL, KIJS, KIJL, IRA)

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

!       *CALL* *INTPOL (GFL, FLA, IJS, IJL, KIJS, KIJL, IRA)*
!         *GFL*   - SPECTRA (INPUT) (1st dimension IJS:IJL).
!         *FLA*   - SPECTRA (OUTPUT) (1st DIMENSION KIJS:KIJL).
!         *KIJS*  - INDEX OF FIRST GRIDPOINT.
!         *KIJL*  - INDEX OF LAST GRIDPOINT.
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

      USE YOWCURR  , ONLY : U        ,V
      USE YOWFRED  , ONLY : FR       ,DFIM     ,COSTH    ,SINTH     ,   &
     &              DELTH  ,FRATIO   ,FLOGSPRDM1
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,EPSMIN
      USE YOWSHAL  , ONLY : TFAK     ,INDEP
      USE YOWSTAT  , ONLY : ISHALLO
      USE YOWTEST  , ONLY : IU06

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, KIJS, KIJL, IRA

      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: GFL

      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(OUT) :: FLA


      INTEGER(KIND=JWIM) :: IJ, M, K
      INTEGER(KIND=JWIM) :: NEWM, NEWM1, KH 
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL) :: NEWF, NEWFLA, KNEW

      REAL(KIND=JWRB) :: PI2G, FRE0, CDF 
      REAL(KIND=JWRB) :: FNEW, GWH
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: DFTH(NFRE)
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: OLDFL 
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL) :: FNEF, GWP, GWM, WAVN

      LOGICAL :: LICE2SEA(KIJS:KIJL)
! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('INTPOL',0,ZHOOK_HANDLE)

!*    0. INITIAL OUTPUT ARRAY WITH ZERO.
!        -------------------------------

      PI2G = ZPI/G 
      FRE0 = FRATIO-1.0_JWRB

!!!??? I believe that in order to be energy conserving, since
!!! DFIM(1) is only the part above FR(1)
!!! purpose of trapezoidal rule, one has to take DFTH(1) which 
!!! define on both side of FR(1). Similarly for FR(NFRE)
!!! Redefine the centered frequency-direction intervals
      CDF= 0.5_JWRB*(FRATIO-1.0_JWRB/FRATIO)*DELTH
      DO M=1,NFRE
        DFTH(M)= FR(M)*CDF
      ENDDO

      IF (ABS(IRA).NE.1) THEN
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
             IF (GFL(IJ,K,M) .GT. EPSMIN) LICE2SEA(IJ) = .FALSE. 
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

      DO M = 1, NFRE
        IF (ISHALLO.NE.1) THEN
          DO IJ = KIJS, KIJL
            WAVN(IJ) = TFAK(INDEP(IJ),M)/ZPI
          ENDDO
        ELSE
          DO IJ = KIJS, KIJL
            WAVN(IJ) = PI2G*FR(M)*FR(M)
          ENDDO
        ENDIF

!*    1.3 LOOP OVER DIRECTONS.
!         --------------------

        DO K = 1, NANG

!*    1.3.1 NEW FREQUENCY AND DIRECTION AT ALL GRIDPOINTS.
!           ----------------------------------------------

          DO IJ = KIJS, KIJL
            FNEF(IJ) = FR(M) + IRA*WAVN(IJ)*(COSTH(K)*V(IJ) + SINTH(K)*U(IJ))
            IF (FNEF(IJ).GT.0.0_JWRB) THEN
              KNEW(IJ) = K
            ELSE
              KNEW(IJ) = MOD(K+NANG/2-1,NANG) + 1
              FNEF(IJ) = -FNEF(IJ)
            ENDIF
          ENDDO

!*    1.3.2 NEW FREQUENCY BIN NUMBER AT ALL GRIDPOINTS.
!           -------------------------------------------

          DO IJ = KIJS, KIJL
            IF (FNEF(IJ).LE.FR(1)/FRATIO) THEN
              NEWF(IJ)= -1
            ELSE
              NEWF(IJ)=FLOOR(LOG10(FNEF(IJ)/FR(1))*FLOGSPRDM1)+1
            ENDIF
          ENDDO

!*    1.3.3 INTERPOLATED ENERGY DENSITIES AT ALL GRIDPOINTS.
!           ------------------------------------------------
          DO IJ = KIJS, KIJL
            IF (LICE2SEA(IJ)) THEN
              OLDFL(IJ)=0.0_JWRB
            ELSE
              OLDFL(IJ)=GFL(IJ,K,M)
            ENDIF
          ENDDO

          DO IJ = KIJS, KIJL
            FNEW = FNEF(IJ)
            NEWM = NEWF(IJ)
            IF (NEWM.LT.NFRE.AND.NEWM.GE.1) THEN
              NEWM1 = NEWM + 1
              GWH = DFTH(M)/(FR(NEWM1)-FR(NEWM)) *OLDFL(IJ)
              GWM(IJ) = GWH*(FR(NEWM1)-FNEW)/DFTH(NEWM)
              GWP(IJ) = GWH*(FNEW-FR(NEWM))/DFTH(NEWM1)
              NEWFLA(IJ) = NEWM1
            ELSEIF (NEWM.EQ.0) THEN
              GWH = FRATIO*DFTH(M)/(FRE0*FR(1)) * OLDFL(IJ)
              GWP(IJ) = GWH*(FNEW-FR(1)/FRATIO)/DFTH(1)
              NEWF (IJ) = -1
              NEWFLA(IJ) = 1
            ELSEIF (NEWM.EQ.NFRE) THEN
              GWH = DFTH(M)/(FRE0*FR(NFRE)) * OLDFL(IJ)
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
            IF (NEWM .NE.-1)                                            &
     &       FLA(IJ,KH,NEWM ) = FLA(IJ,KH,NEWM ) + GWM(IJ)
            IF (NEWM1.NE.-1)                                            &
     &       FLA(IJ,KH,NEWM1) = FLA(IJ,KH,NEWM1) + GWP(IJ)
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
