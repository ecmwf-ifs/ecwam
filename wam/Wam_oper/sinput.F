      SUBROUTINE SINPUT (F, FL, IJS, IJL, THWNEW, USNEW, Z0NEW,
     &                   ROAIRN,ZIDLNEW, SL)
! ----------------------------------------------------------------------

!**** *SINPUT* - COMPUTATION OF INPUT SOURCE FUNCTION.

!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990

!     OPTIMIZED BY : H. GUENTHER

!     MODIFIED BY : 
!       J-R BIDLOT NOVEMBER 1995
!       J-R BIDLOT FEBRUARY 1996-97
!       J-R BIDLOT FEBRUARY 1999 : INTRODUCE ICALL AND NCALL
!       P.A.E.M. JANSSEN MAY 2000 : INTRODUCE GUSTINESS
!       J-R BIDLOT FEBRUARY 2001 : MAKE IT FULLY IMPLICIT BY ONLY
!                                  USING NEW STRESS AND ROUGHNESS. 
!       S. ABDALLA OCTOBER 2001:  INTRODUCTION OF VARIABLE AIR
!                                 DENSITY AND STABILITY-DEPENDENT 
!                                 WIND GUSTINESS

!*    PURPOSE.
!     ---------

!       COMPUTE INPUT SOURCE FUNCTION AND STORE ADDITIVELY INTO NET
!       SOURCE FUNCTION ARRAY, ALSO COMPUTE FUNCTIONAL DERIVATIVE OF
!       INPUT SOURCE FUNCTION.
!
!       GUSTINESS IS INTRODUCED FOLL0WING THE APPROACH OF JANSSEN(1986),
!       USING A GAUSS-HERMITE APPROXIMATION SUGGESTED BY MILES(1997).
!       IN THE PRESENT VERSION ONLY TWO HERMITE POLYNOMIALS ARE UTILISED
!       IN THE EVALUATION OF THE PROBABILITY INTEGRAL. EXPLICITELY ONE THEN
!       FINDS:
!
!             <GAMMA(X)> = 0.5*( GAMMA(X(1+SIG)) + GAMMA(X(1-SIG)) )
!
!       WHERE X IS THE FRICTION VELOCITY AND SIG IS THE RELATIVE GUSTINESS
!       LEVEL.

!**   INTERFACE.
!     ----------

!     *CALL* *SINPUT (F, FL, IJS, IJL, THWNEW, USNEW, Z0NEW,
!    &                   ROAIRN,ZIDLNEW, SL)
!            *F* - SPECTRUM.
!           *FL* - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
!          *IJS* - INDEX OF FIRST GRIDPOINT.
!          *IJL* - INDEX OF LAST GRIDPOINT.
!       *THWNEW* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!        *USNEW* - NEW FRICTION VELOCITY IN M/S.
!        *Z0NEW* - ROUGHNESS LENGTH IN M.
!       *ROAIRN* - AIR DENSITY IN KG/M3
!      *ZIDLNEW* - Zi/L  USED FOR GUSTINESS.
!                  (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!           *SL* - TOTAL SOURCE FUNCTION ARRAY.

!     METHOD.
!     -------

!       SEE REFERENCE.

!     EXTERNALS.
!     ----------

!       NONE.

!     MODIFICATIONS
!     -------------

!     - REMOVAL OF CALL TO CRAY SPECIFIC FUNCTIONS EXPHF AND ALOGHF
!       BY THEIR STANDARD FORTRAN EQUIVALENT EXP and ALOGHF
!     - MODIFIED TO MAKE INTEGRATION SCHEME FULLY IMPLICIT
!     - INTRODUCTION OF VARIABLE AIR DENSITY
!     - INTRODUCTION OF WIND GUSTINESS

!     REFERENCE.
!     ----------

!       P. JANSSEN, J.P.O., 1989.
!       P. JANSSEN, J.P.O., 1991

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,XKAPPA
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO
      USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER   ,YEPS
      USE YOWSHAL  , ONLY : TFAK     ,INDEP
      USE YOWSTAT  , ONLY : ISHALLO

! ----------------------------------------------------------------------

!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: F,FL,SL
      REAL,DIMENSION(NINF:NSUP) :: THWNEW, USNEW, Z0NEW
     &                           , ROAIRN,ZIDLNEW

! ----------------------------------------------------------------------

      REAL, DIMENSION(IJS:IJL,NANG) :: TEMP1, UFAC2
      REAL, DIMENSION(IJS:IJL) :: UCN1, UCN2, ZCN, CM
      REAL, DIMENSION(IJS:IJL) :: SIG_N
      REAL, DIMENSION(IJS:IJL) :: CNSN
      REAL, DIMENSION(NFRE) :: FAC, CONST
      PARAMETER (A = 0.8/1000. , B = 0.08/1000.)

! ----------------------------------------------------------------------

      ONETHIRD = 1./3.
      BG_GUST  = 0.        ! No background gustiness (so 12. is not used)

!*    1. PRECALCULATED ANGULAR DEPENDENCE.
!        ---------------------------------

      DO K=1,NANG
        DO IJ=IJS,IJL
          TEMP1(IJ,K) = COS(TH(K)-THWNEW(IJ))
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
!       In the following U10 is estimated assuming everything is based on u*
        U10 = USNEW(IJ)/XKAPPA*LOG(10./Z0NEW(IJ))
        C_D = A+B*U10
        DC_DDU = B
        SIG_CONV = 1. + 0.5*U10/C_D*DC_DDU
        SIG_N (IJ) = SIG_CONV*(BG_GUST-0.5*ZIDLNEW(IJ))**ONETHIRD
     &                 *USNEW(IJ)/U10
      ENDDO

! ----------------------------------------------------------------------

!*    2. LOOP OVER FREQUENCIES.
!        ----------------------

      CONST1  = BETAMAX/XKAPPA**2 /ROWATER
      DO M=1,NFRE
        FAC(M) = ZPI*FR(M)
        CONST(M)=FAC(M)*CONST1
      ENDDO

      DO M=1,NFRE

!*      INVERSE OF PHASE VELOCITIES.
!       ----------------------------

        IF (ISHALLO.EQ.1) THEN
          DO IJ=IJS,IJL
            CM(IJ) = FAC(M)/G
          ENDDO
        ELSE
          DO IJ=IJS,IJL
            CM(IJ) = TFAK(INDEP(IJ),M)/FAC(M)
          ENDDO
        ENDIF

!*      PRECALCULATE FREQUENCY DEPENDENCE.
!       ----------------------------------

        DO IJ=IJS,IJL
          UCN1(IJ) = USNEW(IJ)*(1.+SIG_N(IJ))*CM(IJ) + ZALP
          UCN2(IJ) = USNEW(IJ)*(1.-SIG_N(IJ))*CM(IJ) + ZALP
          ZCN(IJ) = LOG(G*Z0NEW(IJ)*CM(IJ)**2)
          CNSN(IJ) = CONST(M) * ROAIRN(IJ)
        ENDDO

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------

        DO K=1,NANG
          DO IJ=IJS,IJL
            IF (TEMP1(IJ,K).GT.0.01) THEN
              X    = TEMP1(IJ,K)*UCN1(IJ)
              ZLOG = ZCN(IJ) + XKAPPA/X
              IF (ZLOG.LT.0.) THEN
                ZLOG2X=ZLOG*ZLOG*X
                UFAC2(IJ,K) = CNSN(IJ)*EXP(ZLOG)*ZLOG2X*ZLOG2X
              ELSE
                UFAC2(IJ,K) = 0.
              ENDIF

              X    = TEMP1(IJ,K)*UCN2(IJ)
              ZLOG = ZCN(IJ) + XKAPPA/X
              IF (ZLOG.LT.0.) THEN
                ZLOG2X=ZLOG*ZLOG*X
                UFAC2(IJ,K) = UFAC2(IJ,K)+
     &                        CNSN(IJ)*EXP(ZLOG)*ZLOG2X*ZLOG2X
              ENDIF
            ELSE
              UFAC2(IJ,K) = 0.
            ENDIF
          ENDDO
        ENDDO

!*    2.2 ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
!         ------------------------------------------------

        DO K=1,NANG
          DO IJ=IJS,IJL
            TEMPORARY  = 0.5*UFAC2(IJ,K)
!!!      see UPDATEWD FL must be determined before SL !!!
            FL(IJ,K,M) = TEMPORARY 
            SL(IJ,K,M) = TEMPORARY*F(IJ,K,M)
          ENDDO
        ENDDO

      ENDDO

      RETURN
      END SUBROUTINE SINPUT
