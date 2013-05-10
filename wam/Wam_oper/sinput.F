      SUBROUTINE SINPUT (F, FL, IJS, IJL, THWNEW, USNEW, Z0NEW,
     &                   ROAIRN, ZIDLNEW, SL, XLLWS)
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
!       P.A.E.M. JANSSEN OCTOBER 2008: INTRODUCE DAMPING WHEN WAVES ARE 
!                                      RUNNING FASTER THAN THE WIND.
!       J-R BIDLOT JANUARY 2013: SHALLOW WATER FORMULATION.

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
!    &                   ROAIRN,ZIDLNEW, SL, XLLWS)
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
!         *XLLWS*- 1 WHERE SINPUT IS POSITIVE


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
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER   ,YEPS
      USE YOWSHAL  , ONLY : TFAK     ,INDEP
      USE YOWSTAT  , ONLY : ISHALLO  ,IDAMPING
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      REAL,DIMENSION(IJS:IJL,NANG,NFRE) :: F, FL, SL
      REAL,DIMENSION(IJS:IJL) :: THWNEW, USNEW, Z0NEW, ROAIRN, ZIDLNEW

! ----------------------------------------------------------------------
      REAL, PARAMETER :: A = 0.8/1000.
      REAL, PARAMETER :: B = 0.08/1000.

      INTEGER :: IJ,IJS,IJL,IG,K,M

      REAL :: ONETHIRD, BG_GUST, CONST1, CONST3, XKAPPAD
      REAL :: U10, C_D, DC_DDU, SIG_CONV
      REAL :: RWINV
      REAL :: X1,X2,ZLOG1,ZLOG2,ZLOG2X,XV1,XV2,ZBETA1,ZBETA2
      REAL ZHOOK_HANDLE
      REAL, DIMENSION(NFRE) :: FAC, CONST
      REAL, DIMENSION(IJS:IJL) :: UCN1, UCN2, ZCN, CM, USP, USM
      REAL, DIMENSION(IJS:IJL) :: SH, XK
      REAL, DIMENSION(IJS:IJL) :: SIG_N, XV1D, XV2D
      REAL, DIMENSION(IJS:IJL) :: CNSN
      REAL, DIMENSION(IJS:IJL) :: EPSIL 
      REAL, DIMENSION(IJS:IJL) :: UCN1D,UCN2D
      REAL, DIMENSION(IJS:IJL,NANG) :: TEMP1, UFAC2
      REAL, DIMENSION(IJS:IJL,NANG) :: TEMPD
      REAL, DIMENSION(IJS:IJL,NANG,NFRE) :: XLLWS

      LOGICAL :: L1,L2
      LOGICAL, DIMENSION(IJS:IJL,NANG) :: LZ

      IF (LHOOK) CALL DR_HOOK('SINPUT',0,ZHOOK_HANDLE)

! ----------------------------------------------------------------------

      ONETHIRD = 1./3.
      BG_GUST  = 0.        ! NO BACKGROUND GUSTINESS (S0 12. IS NOT USED)
      CONST1   = BETAMAX/XKAPPA**2 
      CONST3   = 2.*XKAPPA/CONST1  ! SEE IDAMPING
      XKAPPAD  = 1.D0/XKAPPA
      RWINV = 1.0/ROWATER

      CONST3 = IDAMPING*CONST3

!*    1. PRECALCULATED ANGULAR DEPENDENCE.
!        ---------------------------------

      DO K=1,NANG
        DO IJ=IJS,IJL
          TEMP1(IJ,K) = COS(TH(K)-THWNEW(IJ))
          IF(TEMP1(IJ,K) .GT. 0.01) THEN
            LZ(IJ,K) = .TRUE.
            TEMPD(IJ,K) = 1.D0/TEMP1(IJ,K)
          ELSE
            LZ(IJ,K) = .FALSE.
            TEMPD(IJ,K) = 1.D0
          ENDIF
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
!
!       IN THE FOLLOWING U10 IS ESTIMATED ASSUMING EVERYTHING IS
!       BASED ON U*
!
        U10 = USNEW(IJ)/XKAPPA*LOG(10./Z0NEW(IJ))
        C_D = A+B*U10
        DC_DDU = B
        SIG_CONV = 1. + 0.5*U10/C_D*DC_DDU
        SIG_N(IJ) = MIN(0.5, SIG_CONV * 
     &                        (BG_GUST*USNEW(IJ)**3+
     &                         0.5*XKAPPA*ZIDLNEW(IJ)**3)**ONETHIRD
     &                           /U10
     &                 )
      ENDDO

      DO IJ=IJS,IJL
        USP(IJ) = USNEW(IJ)*(1.+SIG_N(IJ))
        USM(IJ) = USNEW(IJ)*(1.-SIG_N(IJ))
      ENDDO

      DO IJ=IJS,IJL
          EPSIL(IJ) = ROAIRN(IJ)*RWINV
      ENDDO
! ----------------------------------------------------------------------

!*    2. LOOP OVER FREQUENCIES.
!        ----------------------


      DO M=1,NFRE

        FAC(M) = ZPI*FR(M)
        CONST(M)=FAC(M)*CONST1

!*      INVERSE OF PHASE VELOCITIES.
!       ----------------------------

        IF (ISHALLO.EQ.1) THEN
          DO IJ=IJS,IJL
            XK(IJ) = FAC(M)**2/G
            CM(IJ) = FAC(M)/G
            SH(IJ) = 1.0
          ENDDO
        ELSE
          DO IJ=IJS,IJL
            XK(IJ) = TFAK(INDEP(IJ),M)
            CM(IJ) = XK(IJ)/FAC(M)
            SH(IJ) = FAC(M)**2/(G*XK(IJ)) 
          ENDDO
        ENDIF

!*      PRECALCULATE FREQUENCY DEPENDENCE.
!       ----------------------------------

        DO IJ=IJS,IJL
          UCN1(IJ) = USP(IJ)*CM(IJ) + ZALP
          UCN2(IJ) = USM(IJ)*CM(IJ) + ZALP

          UCN1D(IJ) = 1.D0/ UCN1(IJ)
          UCN2D(IJ) = 1.D0/ UCN2(IJ)

          ZCN(IJ)  = LOG(XK(IJ)*Z0NEW(IJ))
          CNSN(IJ) = CONST(M)*SH(IJ)*EPSIL(IJ)

          XV1      = -USP(IJ)*XKAPPAD*ZCN(IJ)*CM(IJ)
          XV2      = -USM(IJ)*XKAPPAD*ZCN(IJ)*CM(IJ)

          XV1D(IJ) = 1.D0/XV1
          XV2D(IJ) = 1.D0/XV2
        ENDDO

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------

        DO K=1,NANG
          DO IJ=IJS,IJL
            ZBETA1 = CONST3*(TEMP1(IJ,K)-XV1D(IJ))*UCN1(IJ)**2
            ZBETA2 = CONST3*(TEMP1(IJ,K)-XV2D(IJ))*UCN2(IJ)**2
            IF (LZ(IJ,K)) THEN
              ZLOG1 = ZCN(IJ) + XKAPPA*TEMPD(IJ,K)*UCN1D(IJ)
              ZLOG2 = ZCN(IJ) + XKAPPA*TEMPD(IJ,K)*UCN2D(IJ)
              L1    = ZLOG1.LT.0.
              L2    = ZLOG2.LT.0.
              IF (ZLOG1.LT.0.) THEN
                X1=TEMP1(IJ,K)*UCN1(IJ)
                ZLOG2X=ZLOG1*ZLOG1*X1
                UFAC2(IJ,K) = EXP(ZLOG1)*ZLOG2X*ZLOG2X+ZBETA1
                XLLWS(IJ,K,M)= 1.
              ELSE
                UFAC2(IJ,K) = ZBETA1
                XLLWS(IJ,K,M)= 0.
              ENDIF
              IF (ZLOG2.LT.0.) THEN
                X2=TEMP1(IJ,K)*UCN2(IJ)
                ZLOG2X=ZLOG2*ZLOG2*X2
                UFAC2(IJ,K) = UFAC2(IJ,K)+
     &                        EXP(ZLOG2)*ZLOG2X*ZLOG2X+ZBETA2
                XLLWS(IJ,K,M)= 1.
              ELSE
                UFAC2(IJ,K) = UFAC2(IJ,K)+ZBETA2
              ENDIF
            ELSE
              UFAC2(IJ,K) = ZBETA1+ZBETA2
              XLLWS(IJ,K,M)= 0.
            ENDIF
          ENDDO
        ENDDO

!*    2.2 ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
!         ------------------------------------------------

        DO K=1,NANG
          DO IJ=IJS,IJL
            FL(IJ,K,M) = 0.5*CNSN(IJ)*UFAC2(IJ,K)
            SL(IJ,K,M) = FL(IJ,K,M)*F(IJ,K,M)
          ENDDO
        ENDDO

      ENDDO

      IF (LHOOK) CALL DR_HOOK('SINPUT',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE SINPUT
