      SUBROUTINE INTPOL (F3, F1, IJS, IJL, IG)

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

!       *CALL* *INTPOL (F3, F1, IJS, IJL, IG)*
!         *F3*   - BLOCK OF SPECTRA (INPUT).
!         *F1*   - BLOCK OF SPECTRA (OUTPUT).
!         *IJS*  - INDEX OF FIRST GRIDPOINT.
!         *IJL*  - INDEX OF LAST GRIDPOINT.
!         *IG*   - BLOCK NUMBER.

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
     &         FRATIO
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NIBLC
      USE YOWPCONS , ONLY : G        ,ZPI
      USE YOWSHAL  , ONLY : TFAK     ,INDEP
      USE YOWSTAT  , ONLY : ISHALLO

! ----------------------------------------------------------------------

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: F3,F1

      DIMENSION FNEF(NIBLC), NEWF(NIBLC), NEWF1(NIBLC), KNEW(NIBLC),
     &          GWP(NIBLC), GWM(NIBLC), WAVN(NIBLC)

! ----------------------------------------------------------------------

!*    0. INITIAL OUTPUT ARRAY WITH ZERO.
!        -------------------------------

      PI2G = ZPI/G 
      FRE0 = FRATIO-1.
      ALN00 = 1. / LOG10(FRATIO)
      F1 = 0.0

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
     &       + WAVN(IJ)*(COSTH(K)*V(IJ,IG) + SINTH(K)*U(IJ,IG))
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
              NEWF(IJ) = -1
            ELSE
              NEWF(IJ) =INT(LOG10(FNEF(IJ)/FR(1))*ALN00+1.000000001)
            ENDIF
          ENDDO

!*    1.3.3 INTERPOLATED ENERGY DENSITIES AT ALL GRIDPOINTS.
!           ------------------------------------------------

          DO IJ = IJS, IJL
            FNEW = FNEF(IJ)
            NEWM = NEWF(IJ)
            IF (NEWM.LT.NFRE.AND.NEWM.GE.1) THEN
              NEWM1 = NEWM + 1
              GWH = DFIM(M)/(FR(NEWM1)-FR(NEWM)) *F3(IJ,K,M)
              GWM(IJ) = GWH*(FR(NEWM1)-FNEW)/DFIM(NEWM)
              GWP(IJ) = GWH*(FNEW-FR(NEWM))/DFIM(NEWM1)
              NEWF1(IJ) = NEWM1
            ELSEIF (NEWM.EQ.0) THEN
              GWH = FRATIO*DFIM(M)/(FRE0*FR(1)) * F3(IJ,K,M)
              GWP(IJ) = GWH*(FNEW-FR(1)/FRATIO)/DFIM(1)
              NEWF (IJ) = -1
              NEWF1(IJ) = 1
            ELSEIF (NEWM.EQ.NFRE) THEN
              GWH = DFIM(M)/(FRE0*FR(NFRE)) * F3(IJ,K,M)
              GWM(IJ) = GWH*(FRATIO*FR(NFRE)-FNEW)/DFIM(NFRE)
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

      RETURN
      END SUBROUTINE INTPOL
