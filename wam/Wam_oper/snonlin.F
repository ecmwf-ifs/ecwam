      SUBROUTINE SNONLIN (F, FL, IJS, IJL, IG, SL, AKMEAN)

! ----------------------------------------------------------------------

!**** *SNONLIN* - COMPUTATION OF NONLINEAR TRANSFER RATE AND ITS
!****             FUNCTIONAL DERIVATIVE (DIAGONAL TERMS ONLY) AND
!****             ADDITION TO CORRESPONDING NET EXPRESSIONS.

!     S.D. HASSELMANN.  MPI

!     G. KOMEN, P. JANSSEN   KNMI             MODIFIED TO SHALLOW WATER
!     H. GUENTHER, L. ZAMBRESKY               OPTIMIZED
!     H. GUENTHER       GKSS/ECMWF  JUNE 1991 INTERACTIONS BETWEEN DIAG-
!                                             AND PROGNOSTIC PART.
!     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
!     P. JANSSEN  ECMWF  JUNE 2005       IMPROVED SCALING IN SHALLOW
!                                        WATER
!     J. BIDLOT   ECMWF  AUGUST 2006     KEEP THE OLD FORMULATION
!                                        UNDER A SWITCH (ISNONLIN = 0 for OLD
!                                                                 = 1 for NEW
!                                        BE AWARE THAT THE OLD FORMULATION
!                                        REQUIRES THE MEAN WAVE NUMBER AKMEAN.

!*    PURPOSE.
!     --------

!       SEE ABOVE.

!**   INTERFACE.
!     ----------

!       *CALL* *SNONLIN (F, FL, IJS, IJL, IG, SL)*
!          *F*   - SPECTRUM.
!          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *IG*  - BLOCK NUMBER.
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWINDN  , ONLY : IKP      ,IKP1     ,IKM      ,IKM1     ,
     &            K1W      ,K2W      ,K11W     ,K21W     ,AF11     ,
     &            FKLAP    ,FKLAP1   ,FKLAM    ,FKLAM1   ,ACL1     ,
     &            ACL2     ,CL11     ,CL21     ,DAL1     ,DAL2     ,
     &            FRH      ,ENH
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWSHAL  , ONLY : DEPTH    ,TFAK,    INDEP 
      USE YOWSTAT  , ONLY : ISHALLO  ,ISNONLIN
! ----------------------------------------------------------------------

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: F,FL,SL

      DIMENSION FTEMP(IJS:IJL), AD(IJS:IJL), DELAD(IJS:IJL), 
     &          DELAP(IJS:IJL), DELAM(IJS:IJL)

      REAL :: AKMEAN(IJS:IJL)
      REAL :: ENHFR(IJS:IJL)

! ----------------------------------------------------------------------

!*    1. SHALLOW WATER INITIALISATION (ONLY FOR THE OLD FORMULATION).
!        ----------------------------- SEE INITMDL FOR THE NEW ONE

      IF (ISHALLO.NE.1) THEN
        IF (ISNONLIN.EQ.0) THEN
          DO IJ=IJS,IJL
            ENHFR(IJ) = 0.75*DEPTH(IJ,IG)*AKMEAN(IJ)
            ENHFR(IJ) = MAX(ENHFR(IJ),.5)
            ENHFR(IJ) = 1.+(5.5/ENHFR(IJ))*(1.-.833*ENHFR(IJ)) *
     &            EXP(-1.25*ENHFR(IJ))
          ENDDO
          DO MC=1,NFRE+4
            DO IJ=IJS,IJL
              ENH(IJ,MC,IG) = ENHFR(IJ) 
            ENDDO
          ENDDO
        ENDIF
      ENDIF


!*    2. FREQUENCY LOOP.
!        ---------------

      DO MC=1,NFRE+4
        MP  = IKP (MC)
        MP1 = IKP1(MC)
        MM  = IKM (MC)
        MM1 = IKM1(MC)
        FFACP  = 1.
        FFACP1 = 1.
        FFACM1 = 1.
        FTAIL  = 1.
        IC  = MC
        IP  = MP
        IP1 = MP1
        IM  = MM
        IM1 = MM1
        IF (IP1.GT.NFRE) THEN
          FFACP1 = FRH(IP1-NFRE+1)
          IP1 = NFRE
          IF (IP .GT.NFRE) THEN
            FFACP  = FRH(IP -NFRE+1)
            IP  = NFRE
            IF (IC .GT.NFRE) THEN
              FTAIL  = FRH(IC -NFRE+1)
              IC  = NFRE
              IF (IM1.GT.NFRE) THEN
                FFACM1 = FRH(IM1-NFRE+1)
                IM1 = NFRE
              ENDIF
            ENDIF
          ENDIF
        ENDIF
        FKLAMP  = FKLAP(MC)
        FKLAMP1 = FKLAP1(MC)
        GW2 = FKLAMP1*FFACP*DAL1
        GW1 = GW2*CL11
        GW2 = GW2*ACL1
        GW4 = FKLAMP*FFACP1*DAL1
        GW3 = GW4*CL11
        GW4 = GW4*ACL1
        FKLAMPA = FKLAMP*CL11
        FKLAMPB = FKLAMP*ACL1
        FKLAMP2 = FKLAMP1*ACL1
        FKLAMP1 = FKLAMP1*CL11
        FKLAPA2 = FKLAMPA**2
        FKLAPB2 = FKLAMPB**2
        FKLAP12 = FKLAMP1**2
        FKLAP22 = FKLAMP2**2

        FKLAMM  = FKLAM(MC)
        FKLAMM1 = FKLAM1(MC)
        GW6 = FKLAMM1*DAL2
        GW5 = GW6*CL21
        GW6 = GW6*ACL2
        GW8 = FKLAMM*FFACM1*DAL2
        GW7 = GW8*CL21
        GW8 = GW8*ACL2
        FKLAMMA = FKLAMM*CL21
        FKLAMMB = FKLAMM*ACL2
        FKLAMM2 = FKLAMM1*ACL2
        FKLAMM1 = FKLAMM1*CL21
        FKLAMA2 = FKLAMMA**2
        FKLAMB2 = FKLAMMB**2
        FKLAM12 = FKLAMM1**2
        FKLAM22 = FKLAMM2**2

        IF (ISHALLO.EQ.1) THEN
          DO IJ=IJS,IJL
            FTEMP(IJ) = AF11(MC)
          ENDDO
        ELSE
!SHALLOW
          DO IJ=IJS,IJL
            FTEMP(IJ) = AF11(MC)*ENH(IJ,MC,IG)
          ENDDO
!SHALLOW
        ENDIF

!*    2.1 LOOP FOR ANGULAR SYMMETRY.
!         -------------------------

        DO KH=1,2

!*    2.1.1   ANGULAR LOOP.
!             -------------

          DO K=1,NANG
            K1  = K1W (K,KH)
            K2  = K2W (K,KH)
            K11 = K11W(K,KH)
            K21 = K21W(K,KH)

!*    2.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND
!*            DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
!             ----------------------------------------------

            IF (MC.GT.4) THEN
              DO IJ=IJS,IJL
                SAP = GW1*F(IJ,K1 ,IP ) + GW2*F(IJ,K11,IP )
     &              + GW3*F(IJ,K1 ,IP1) + GW4*F(IJ,K11,IP1)
                SAM = GW5*F(IJ,K2 ,IM ) + GW6*F(IJ,K21,IM )
     &              + GW7*F(IJ,K2 ,IM1) + GW8*F(IJ,K21,IM1)
                FIJ = F(IJ,K  ,IC )*FTAIL
                FAD1 = FIJ*(SAP+SAM)
                FAD2 = FAD1-2.*SAP*SAM
                FAD1 = FAD1+FAD2
                FCEN = FTEMP(IJ)*FIJ
                AD(IJ) = FAD2*FCEN
                DELAD(IJ) = FAD1*FTEMP(IJ)
                DELAP(IJ) = (FIJ-2.*SAM)*DAL1*FCEN
                DELAM(IJ) = (FIJ-2.*SAP)*DAL2*FCEN
              ENDDO

              DO IJ=IJS,IJL
                SL(IJ,K2 ,MM ) = SL(IJ,K2 ,MM ) + AD(IJ)*FKLAMM1
                SL(IJ,K21,MM ) = SL(IJ,K21,MM ) + AD(IJ)*FKLAMM2

                FL(IJ,K2 ,MM ) = FL(IJ,K2 ,MM ) + DELAM(IJ)*FKLAM12
                FL(IJ,K21,MM ) = FL(IJ,K21,MM ) + DELAM(IJ)*FKLAM22
              ENDDO

              IF (MM1.LE.NFRE) THEN
                DO IJ=IJS,IJL
                  SL(IJ,K2 ,MM1) = SL(IJ,K2 ,MM1) + AD(IJ)*FKLAMMA
                  SL(IJ,K21,MM1) = SL(IJ,K21,MM1) + AD(IJ)*FKLAMMB

                  FL(IJ,K2 ,MM1) = FL(IJ,K2 ,MM1) + DELAM(IJ)*FKLAMA2
                  FL(IJ,K21,MM1) = FL(IJ,K21,MM1) + DELAM(IJ)*FKLAMB2
                ENDDO

                IF (MC .LE.NFRE) THEN
                  DO IJ=IJS,IJL
                    SL(IJ,K  ,MC ) = SL(IJ,K  ,MC ) - 2.*AD(IJ)

                    FL(IJ,K  ,MC ) = FL(IJ,K  ,MC ) - 2.*DELAD(IJ)
                  ENDDO

                  IF (MP .LE.NFRE) THEN
                    DO IJ=IJS,IJL
                      SL(IJ,K1 ,MP ) = SL(IJ,K1 ,MP ) + AD(IJ)*FKLAMP1
                      SL(IJ,K11,MP ) = SL(IJ,K11,MP ) + AD(IJ)*FKLAMP2

                      FL(IJ,K1 ,MP ) = FL(IJ,K1 ,MP )
     &                               + DELAP(IJ)*FKLAP12
                      FL(IJ,K11,MP ) = FL(IJ,K11,MP )
     &                               + DELAP(IJ)*FKLAP22
                    ENDDO

                    IF (MP1.LE.NFRE) THEN
                      DO IJ=IJS,IJL
                        SL(IJ,K1 ,MP1) = SL(IJ,K1 ,MP1)
     &                                 + AD(IJ)*FKLAMPA
                        SL(IJ,K11,MP1) = SL(IJ,K11,MP1)
     &                                 + AD(IJ)*FKLAMPB

                        FL(IJ,K1 ,MP1) = FL(IJ,K1 ,MP1)
     &                                 + DELAP(IJ)*FKLAPA2
                        FL(IJ,K11,MP1) = FL(IJ,K11,MP1)
     &                                 + DELAP(IJ)*FKLAPB2
                      ENDDO
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF

            ELSE

              DO IJ=IJS,IJL
                SAP = GW1*F(IJ,K1 ,IP ) + GW2*F(IJ,K11,IP )
     &              + GW3*F(IJ,K1 ,IP1) + GW4*F(IJ,K11,IP1)
                FIJ = F(IJ,K,IM)
                FAD2 = FIJ*SAP
                FAD1 = 2.*FAD2
                FCEN = FTEMP(IJ)*FIJ
                AD(IJ) = FAD2*FCEN
                DELAD(IJ) = FAD1*FTEMP(IJ)
                DELAP(IJ) = FIJ*DAL1*FCEN
              ENDDO

              DO IJ=IJS,IJL
                SL(IJ,K  ,MC ) = SL(IJ,K  ,MC ) - 2.*AD(IJ)
                SL(IJ,K1 ,MP ) = SL(IJ,K1 ,MP ) + AD(IJ)*FKLAMP1
                SL(IJ,K11,MP ) = SL(IJ,K11,MP ) + AD(IJ)*FKLAMP2
                SL(IJ,K1 ,MP1) = SL(IJ,K1 ,MP1) + AD(IJ)*FKLAMPA
                SL(IJ,K11,MP1) = SL(IJ,K11,MP1) + AD(IJ)*FKLAMPB

                FL(IJ,K  ,MC ) = FL(IJ,K  ,MC ) - 2.*DELAD(IJ)
                FL(IJ,K1 ,MP ) = FL(IJ,K1 ,MP ) + DELAP(IJ)*FKLAP12
                FL(IJ,K11,MP ) = FL(IJ,K11,MP ) + DELAP(IJ)*FKLAP22
                FL(IJ,K1 ,MP1) = FL(IJ,K1 ,MP1) + DELAP(IJ)*FKLAPA2
                FL(IJ,K11,MP1) = FL(IJ,K11,MP1) + DELAP(IJ)*FKLAPB2
              ENDDO
            ENDIF

!*    BRANCH BACK TO 2.1.1 FOR NEXT DIRECTION.

          ENDDO

!*    BRANCH BACK TO 2.1 FOR MIRROR INTERACTIONS.

        ENDDO

!*    BRANCH BACK TO 2. FOR NEXT FREQUENCY.

      ENDDO

      RETURN
      END SUBROUTINE SNONLIN
