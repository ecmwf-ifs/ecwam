      SUBROUTINE SINPUT (F, FL, IJS, IJL, THWNEW, USNEW, Z0NEW,
     &                   ROAIRN, WSTAR, SL, XLLWS)
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
!    &                   ROAIRN, WSTAR, SL, XLLWS)
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
!        *WSTAR* - FREE CONVECTION VELOCITY SCALE (M/S).
!           *SL* - TOTAL SOURCE FUNCTION ARRAY.
!         *XLLWS*- 1 WHERE SINPUT IS POSITIVE


!     METHOD.
!     -------

!       SEE REFERENCE.

!     EXTERNALS.
!     ----------

!       WSIGSTAR. 

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
      USE YOWTEST  , ONLY : IU06
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL
      INTEGER, PARAMETER :: NSIN=2 
      INTEGER :: IJ, IG, K, M
      INTEGER :: ISIN

      REAL,DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F
      REAL,DIMENSION(IJS:IJL,NANG,NFRE), INTENT(OUT) :: FL, SL, XLLWS
      REAL,DIMENSION(IJS:IJL), INTENT(IN) :: THWNEW, USNEW, Z0NEW
      REAL,DIMENSION(IJS:IJL), INTENT(IN) :: ROAIRN, WSTAR

      REAL :: CONST1, CONST3, XKAPPAD
      REAL :: RWINV

      REAL :: X, ZLOG, ZLOG2X, ZBETA

      REAL :: ZHOOK_HANDLE
      REAL, DIMENSION(NSIN) :: WSIN
      REAL, DIMENSION(NFRE) :: FAC, CONST
      REAL, DIMENSION(IJS:IJL) :: CM
      REAL, DIMENSION(IJS:IJL) :: SH, XK
      REAL, DIMENSION(IJS:IJL) :: SIG_N
      REAL, DIMENSION(IJS:IJL) :: CNSN
      REAL, DIMENSION(IJS:IJL) :: EPSIL 
      REAL, DIMENSION(IJS:IJL,NSIN) :: SIGDEV, US, Z0, UCN, ZCN 
      REAL, DIMENSION(IJS:IJL,NSIN) :: XVD, UCND, CONST3_UCN2 
      REAL, DIMENSION(IJS:IJL,NANG) :: TEMP1, UFAC2
      REAL, DIMENSION(IJS:IJL,NANG) :: TEMPD

      LOGICAL, DIMENSION(IJS:IJL,NANG) :: LZ

! ----------------------------------------------------------------------

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SINPUT',0,ZHOOK_HANDLE)
#endif

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
            TEMPD(IJ,K) = XKAPPA/TEMP1(IJ,K)
          ELSE
            LZ(IJ,K) = .FALSE.
            TEMPD(IJ,K) = XKAPPA 
          ENDIF
        ENDDO
      ENDDO


!     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
      IF(NSIN.GT.1) CALL WSIGSTAR (IJS, IJL, USNEW, Z0NEW, WSTAR, SIG_N)


!     DEFINE WHERE SINPUT WILL BE EVALUATED IN RELATIVE TERM WRT USTAR
!     DEFINE ALSO THE RELATIVE WEIGHT OF EACH.

      IF(NSIN.EQ.1) THEN
        WSIN(1) = 1.0 
        DO IJ=IJS,IJL
          SIGDEV(IJ,1) = 1.0
        ENDDO
      ELSE IF (NSIN.EQ.2) THEN
        WSIN(1) = 0.5 
        WSIN(2) = 0.5 
        DO IJ=IJS,IJL
          SIGDEV(IJ,1) = 1.-SIG_N(IJ)
          SIGDEV(IJ,2) = 1.+SIG_N(IJ)
        ENDDO
      ELSE
         WRITE (IU06,*) '**************************************'
         WRITE (IU06,*) '*    FATAL ERROR                     *'
         WRITE (IU06,*) '*    ===========                     *'
         WRITE (IU06,*) '* IN SINPUT: NSIN > 2                *'
         WRITE (IU06,*) '* NSIN = ', NSIN
         WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.  *'
         WRITE (IU06,*) '*                                    *'
         WRITE (IU06,*) '**************************************'
         CALL ABORT1
      ENDIF


      IF(NSIN.EQ.1) THEN
        DO IJ=IJS,IJL
          US(IJ,1) = USNEW(IJ)
          Z0(IJ,1) = Z0NEW(IJ)
        ENDDO
      ELSE
        DO ISIN=1,NSIN
          DO IJ=IJS,IJL
            US(IJ,ISIN) = USNEW(IJ)*SIGDEV(IJ,ISIN)
            Z0(IJ,ISIN) = Z0NEW(IJ)
          ENDDO
        ENDDO
      ENDIF

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
          CNSN(IJ) = CONST(M)*SH(IJ)*EPSIL(IJ)
        ENDDO

        DO ISIN=1,NSIN
          DO IJ=IJS,IJL
            UCN(IJ,ISIN) = US(IJ,ISIN)*CM(IJ) + ZALP
            CONST3_UCN2(IJ,ISIN) = CONST3*UCN(IJ,ISIN)**2
            UCND(IJ,ISIN) = 1.D0/ UCN(IJ,ISIN)
            ZCN(IJ,ISIN)  = LOG(XK(IJ)*Z0(IJ,ISIN))
            XVD(IJ,ISIN) =
     &          1.D0/(-US(IJ,ISIN)*XKAPPAD*ZCN(IJ,ISIN)*CM(IJ))
          ENDDO
        ENDDO

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------

        DO K=1,NANG

          DO IJ=IJS,IJL
            XLLWS(IJ,K,M)= 0.
          ENDDO

          DO ISIN=1,1
            DO IJ=IJS,IJL
              ZBETA = CONST3_UCN2(IJ,ISIN)*(TEMP1(IJ,K)-XVD(IJ,ISIN))
              UFAC2(IJ,K) = WSIN(ISIN)*ZBETA
            ENDDO
          ENDDO
          DO ISIN=2,NSIN
            DO IJ=IJS,IJL
              ZBETA = CONST3_UCN2(IJ,ISIN)*(TEMP1(IJ,K)-XVD(IJ,ISIN))
              UFAC2(IJ,K) = UFAC2(IJ,K)+WSIN(ISIN)*ZBETA
            ENDDO
          ENDDO

          DO ISIN=1,NSIN
            DO IJ=IJS,IJL
              IF (LZ(IJ,K)) THEN
                ZLOG = ZCN(IJ,ISIN) + TEMPD(IJ,K)*UCND(IJ,ISIN)
                IF (ZLOG.LT.0.) THEN
                  X=TEMP1(IJ,K)*UCN(IJ,ISIN)
                  ZLOG2X=ZLOG*ZLOG*X
                  UFAC2(IJ,K) = UFAC2(IJ,K)+
     &                          WSIN(ISIN)*EXP(ZLOG)*ZLOG2X*ZLOG2X
                  XLLWS(IJ,K,M)= 1.
                ENDIF
              ENDIF
            ENDDO
          ENDDO

        ENDDO

!*    2.2 ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
!         ------------------------------------------------

        DO K=1,NANG
          DO IJ=IJS,IJL
            FL(IJ,K,M) = CNSN(IJ)*UFAC2(IJ,K)
            SL(IJ,K,M) = FL(IJ,K,M)*F(IJ,K,M)
          ENDDO
        ENDDO

      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SINPUT',1,ZHOOK_HANDLE)
#endif

      RETURN
      END SUBROUTINE SINPUT
