      SUBROUTINE INTPOL (F3, F1, IJS, IJL, IG, IRA)

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

!       *CALL* *INTPOL (F3, F1, IJS, IJL, IG, IRA)*
!         *F3*   - SPECTRA (INPUT).
!         *F1*   - SPECTRA (OUTPUT).
!         *IJS*  - INDEX OF FIRST GRIDPOINT.
!         *IJL*  - INDEX OF LAST GRIDPOINT.
!         *IG*   - BLOCK NUMBER.
!         *IRA*  - = 1 TRANSFORMATION FROM MOVING TO ABSOLUTE COORD.
!                  =-1 TRANSFORMATION FROM ABSOLUTE TO MOVING COORD.

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

      USE YOWCURR  , ONLY : U        ,V
      USE YOWFRED  , ONLY : FR       ,DFIM     ,COSTH    ,SINTH     ,
     &              DELTH  ,FRATIO   ,FLOGSPRDM1
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWICE   , ONLY : FLMIN
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,EPSMIN
      USE YOWSHAL  , ONLY : TFAK     ,INDEP
      USE YOWSTAT  , ONLY : ISHALLO
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: IJ, IJS, IJL, IG, IRA, M, K
      INTEGER :: NEWM, NEWM1, KH 
      INTEGER, DIMENSION(IJS:IJL) :: NEWF,NEWF1,KNEW

      REAL :: PI2G, FRE0, CDF 
      REAL :: FNEW, GWH
      REAL :: ZHOOK_HANDLE
      REAL :: DFTH(NFRE)
      REAL, DIMENSION(IJS:IJL) :: FNEF,GWP,GWM,WAVN
      REAL, DIMENSION(NINF-1:NSUP,NANG,NFRE) :: F3, F1

      LOGICAL :: LICE2SEA(IJS:IJL)
! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('INTPOL',0,ZHOOK_HANDLE)


!*    0. INITIAL OUTPUT ARRAY WITH ZERO.
!        -------------------------------

      PI2G = ZPI/G 
      FRE0 = FRATIO-1.

!!!??? I believe that in order to be energy conserving, since
!!! DFIM(1) is only the part above FR(1)
!!! purpose of trapezoidal rule, one has to take DFTH(1) which 
!!! define on both side of FR(1). Similarly for FR(NFRE)
!!! Redefine the centered frequency-direction intervals
      CDF= 0.5*(FRATIO-1./FRATIO)*DELTH
      DO M=1,NFRE
        DFTH(M)= FR(M)*CDF
      ENDDO

      IF(ABS(IRA).NE.1) THEN
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

      DO IJ = IJS, IJL
        LICE2SEA(IJ)=.TRUE.
      ENDDO

      DO K=1,NANG
        DO M=1,NFRE
          DO IJ = IJS, IJL
             IF(F3(IJ,K,M) .GT. EPSMIN) LICE2SEA(IJ) = .FALSE. 
          ENDDO
        ENDDO
      ENDDO

      DO M = 1, NFRE
        DO K = 1, NANG
          DO IJ = IJS, IJL
            IF(LICE2SEA(IJ)) THEN
              F3(IJ,K,M) = 0.0
            ENDIF
          ENDDO
        ENDDO
      ENDDO


      DO M = 1, NFRE
        DO K = 1, NANG
          DO IJ = IJS, IJL
            F1(IJ,K,M) = 0.0
          ENDDO
        ENDDO
      ENDDO

!*    1. LOOP OVER FREQUENCIES.
!        ----------------------

      DO M = 1, NFRE
        IF (ISHALLO.NE.1) THEN
          DO IJ = IJS, IJL
            WAVN(IJ) = TFAK(INDEP(IJ),M)/ZPI
          ENDDO
        ELSE
          DO IJ = IJS, IJL
            WAVN(IJ) = PI2G*FR(M)*FR(M)
          ENDDO
        ENDIF

!*    1.3 LOOP OVER DIRECTONS.
!         --------------------

        DO K = 1, NANG

!*    1.3.1 NEW FREQUENCY AND DIRECTION AT ALL GRIDPOINTS.
!           ----------------------------------------------

          DO IJ = IJS, IJL
            FNEF(IJ) = FR(M)
     &       + IRA*WAVN(IJ)*(COSTH(K)*V(IJ,IG) + SINTH(K)*U(IJ,IG))
            IF (FNEF(IJ).GT.0.) THEN
              KNEW(IJ) = K
            ELSE
              KNEW(IJ) = MOD(K+NANG/2-1,NANG) + 1
              FNEF(IJ) = -FNEF(IJ)
            ENDIF
          ENDDO

!*    1.3.2 NEW FREQUENCY BIN NUMBER AT ALL GRIDPOINTS.
!           -------------------------------------------

          DO IJ = IJS, IJL
            IF (FNEF(IJ).LE.FR(1)/FRATIO) THEN
              NEWF(IJ)= -1
            ELSE
              NEWF(IJ)=FLOOR(LOG10(FNEF(IJ)/FR(1))*FLOGSPRDM1)+1
            ENDIF
          ENDDO

!*    1.3.3 INTERPOLATED ENERGY DENSITIES AT ALL GRIDPOINTS.
!           ------------------------------------------------

          DO IJ = IJS, IJL
            FNEW = FNEF(IJ)
            NEWM = NEWF(IJ)
            IF (NEWM.LT.NFRE.AND.NEWM.GE.1) THEN
              NEWM1 = NEWM + 1
              GWH = DFTH(M)/(FR(NEWM1)-FR(NEWM)) *F3(IJ,K,M)
              GWM(IJ) = GWH*(FR(NEWM1)-FNEW)/DFTH(NEWM)
              GWP(IJ) = GWH*(FNEW-FR(NEWM))/DFTH(NEWM1)
              NEWF1(IJ) = NEWM1
            ELSEIF (NEWM.EQ.0) THEN
              GWH = FRATIO*DFTH(M)/(FRE0*FR(1)) * F3(IJ,K,M)
              GWP(IJ) = GWH*(FNEW-FR(1)/FRATIO)/DFTH(1)
              NEWF (IJ) = -1
              NEWF1(IJ) = 1
            ELSEIF (NEWM.EQ.NFRE) THEN
              GWH = DFTH(M)/(FRE0*FR(NFRE)) * F3(IJ,K,M)
              GWM(IJ) = GWH*(FRATIO*FR(NFRE)-FNEW)/DFTH(NFRE)
              NEWF1(IJ) = -1
            ELSE
              NEWF (IJ) = -1
              NEWF1(IJ) = -1
            ENDIF
          ENDDO

!*    1.3.4 NEW SPECTRUM AT ALL GRIDPOINTS.
!           -------------------------------

          DO IJ = IJS, IJL
            NEWM  = NEWF (IJ)
            NEWM1 = NEWF1(IJ)
            KH = KNEW(IJ)
            IF (NEWM .NE.-1)
     &       F1(IJ,KH,NEWM ) = F1(IJ,KH,NEWM ) + GWM(IJ)
            IF (NEWM1.NE.-1)
     &       F1(IJ,KH,NEWM1) = F1(IJ,KH,NEWM1) + GWP(IJ)
          ENDDO

!*    BRANCH BACK TO 1.3 FOR NEXT DIRECTION.

        ENDDO

!*    BRANCH BACK TO 1. FOR NEXT FREQUENCY.

      ENDDO

      DO M = 1, NFRE
        DO K = 1, NANG
          DO IJ = IJS, IJL
            F1(IJ,K,M) = MAX(F1(IJ,K,M),EPSMIN)
          ENDDO
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('INTPOL',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE INTPOL
