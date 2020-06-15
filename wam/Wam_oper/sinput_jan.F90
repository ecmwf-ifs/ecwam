      SUBROUTINE SINPUT_JAN (NGST, F, FL, IJS, IJL, THWNEW, USNEW, Z0NEW,     &
     &                       ROAIRN, WSTAR, SL, SPOS, XLLWS)
! ----------------------------------------------------------------------

!**** *SINPUT_JAN* - COMPUTATION OF INPUT SOURCE FUNCTION.

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

!     *CALL* *SINPUT_JAN (NGST, F, FL, IJS, IJL, THWNEW, USNEW, Z0NEW,
!    &                   ROAIRN, WSTAR, SL, SPOS, XLLWS)
!         *NGST* - IF = 1 THEN NO GUSTINESS PARAMETERISATION
!                - IF = 2 THEN GUSTINESS PARAMETERISATION
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
!         *SPOS* - ONLY POSITIVE PART OF INPUT SOURCE FUNCTION ARRAY.
!        *XLLWS* - 1 WHERE SINPUT IS POSITIVE.


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

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER   ,YEPS
      USE YOWPHYS  , ONLY : ZALP     ,XKAPPA,  BETAMAXOXKAPPA2
      USE YOWSHAL  , ONLY : TFAK     ,INDEP
      USE YOWSTAT  , ONLY : ISHALLO  ,IDAMPING
      USE YOWTEST  , ONLY : IU06
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "wsigstar.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: NGST 
      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL

      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: F
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: THWNEW, USNEW, Z0NEW
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: ROAIRN, WSTAR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(OUT) :: FL, SL, SPOS, XLLWS


      INTEGER(KIND=JWIM) :: IJ, IG, K, M
      INTEGER(KIND=JWIM) :: IGST

      REAL(KIND=JWRB) :: CONST1, CONST3, XKAPPAD
      REAL(KIND=JWRB) :: RWINV
      REAL(KIND=JWRB) :: X, ZLOG, ZLOG2X, ZBETA
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NGST) :: WSIN
      REAL(KIND=JWRB), DIMENSION(NFRE) :: FAC, CONST
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: CM
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: SH, XK
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: SIG_N
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: CNSN
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: EPSIL 
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NGST) :: SIGDEV,US,Z0,UCN,ZCN
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NGST) :: XVD, UCND, CONST3_UCN2
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG) :: TEMP1, UFAC1, UFAC2
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG) :: TEMPD

      LOGICAL, DIMENSION(IJS:IJL,NANG) :: LZ

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SINPUT_JAN',0,ZHOOK_HANDLE)

      CONST1   = BETAMAXOXKAPPA2 
      CONST3   = 2.0_JWRB*XKAPPA/CONST1  ! SEE IDAMPING
      XKAPPAD  = 1.E0_JWRB/XKAPPA
      RWINV = 1.0_JWRB/ROWATER

      CONST3 = IDAMPING*CONST3


!*    1. PRECALCULATED ANGULAR DEPENDENCE.
!        ---------------------------------

      DO K=1,NANG
        DO IJ=IJS,IJL
          TEMP1(IJ,K) = COS(TH(K)-THWNEW(IJ))
          IF(TEMP1(IJ,K) .GT. 0.01_JWRB) THEN
            LZ(IJ,K) = .TRUE.
            TEMPD(IJ,K) = XKAPPA/TEMP1(IJ,K)
          ELSE
            LZ(IJ,K) = .FALSE.
            TEMPD(IJ,K) = XKAPPA 
          ENDIF
        ENDDO
      ENDDO


!     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
      IF(NGST.GT.1) CALL WSIGSTAR (IJS, IJL, USNEW, Z0NEW, WSTAR, SIG_N)


!     DEFINE WHERE SINPUT WILL BE EVALUATED IN RELATIVE TERM WRT USTAR
!     DEFINE ALSO THE RELATIVE WEIGHT OF EACH.

      IF(NGST.EQ.1) THEN
        WSIN(1) = 1.0_JWRB 
        DO IJ=IJS,IJL
          SIGDEV(IJ,1) = 1.0_JWRB
        ENDDO
      ELSE IF (NGST.EQ.2) THEN
        WSIN(1) = 0.5_JWRB 
        WSIN(2) = 0.5_JWRB 
        DO IJ=IJS,IJL
          SIGDEV(IJ,1) = 1.0_JWRB-SIG_N(IJ)
          SIGDEV(IJ,2) = 1.0_JWRB+SIG_N(IJ)
        ENDDO
      ELSE
         WRITE (IU06,*) '**************************************'
         WRITE (IU06,*) '*    FATAL ERROR                     *'
         WRITE (IU06,*) '*    ===========                     *'
         WRITE (IU06,*) '* IN SINPUT_JAN: NGST > 2            *'
         WRITE (IU06,*) '* NGST = ', NGST
         WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.  *'
         WRITE (IU06,*) '*                                    *'
         WRITE (IU06,*) '**************************************'
         CALL ABORT1
      ENDIF


      IF(NGST.EQ.1) THEN
        DO IJ=IJS,IJL
          US(IJ,1) = USNEW(IJ)
          Z0(IJ,1) = Z0NEW(IJ)
        ENDDO
      ELSE
        DO IGST=1,NGST
          DO IJ=IJS,IJL
            US(IJ,IGST) = USNEW(IJ)*SIGDEV(IJ,IGST)
            Z0(IJ,IGST) = Z0NEW(IJ)
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
            SH(IJ) = 1.0_JWRB
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

        DO IGST=1,NGST
          DO IJ=IJS,IJL
            UCN(IJ,IGST) = US(IJ,IGST)*CM(IJ) + ZALP
            CONST3_UCN2(IJ,IGST) = CONST3*UCN(IJ,IGST)**2
            UCND(IJ,IGST) = 1.0_JWRB/ UCN(IJ,IGST)
            ZCN(IJ,IGST)  = LOG(XK(IJ)*Z0(IJ,IGST))
            XVD(IJ,IGST) = 1.0_JWRB/(-US(IJ,IGST)*XKAPPAD*ZCN(IJ,IGST)*CM(IJ))
          ENDDO
        ENDDO

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------

        DO K=1,NANG

!         WIND INPUT:
          DO IJ=IJS,IJL
            XLLWS(IJ,K,M)= 0.0_JWRB
          ENDDO

          DO IJ=IJS,IJL
            UFAC1(IJ,K) = 0.0_JWRB
          ENDDO
          DO IGST=1,NGST
            DO IJ=IJS,IJL
              IF (LZ(IJ,K)) THEN
                ZLOG = ZCN(IJ,IGST) + TEMPD(IJ,K)*UCND(IJ,IGST)
                IF (ZLOG.LT.0.0_JWRB) THEN
                  X=TEMP1(IJ,K)*UCN(IJ,IGST)
                  ZLOG2X=ZLOG*ZLOG*X
                  UFAC1(IJ,K) = UFAC1(IJ,K)+WSIN(IGST)*EXP(ZLOG)*ZLOG2X*ZLOG2X
                  XLLWS(IJ,K,M)= 1.0_JWRB
                ENDIF
              ENDIF
            ENDDO
          ENDDO

!         SWELL DAMPING:
          DO IGST=1,1
            DO IJ=IJS,IJL
              ZBETA = CONST3_UCN2(IJ,IGST)*(TEMP1(IJ,K)-XVD(IJ,IGST))
              UFAC2(IJ,K) = WSIN(IGST)*ZBETA
            ENDDO
          ENDDO
          DO IGST=2,NGST
            DO IJ=IJS,IJL
              ZBETA = CONST3_UCN2(IJ,IGST)*(TEMP1(IJ,K)-XVD(IJ,IGST))
              UFAC2(IJ,K) = UFAC2(IJ,K)+WSIN(IGST)*ZBETA
            ENDDO
          ENDDO

        ENDDO

!*    2.2 ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
!         ------------------------------------------------

        DO K=1,NANG
          DO IJ=IJS,IJL
            SPOS(IJ,K,M) = CNSN(IJ)*UFAC1(IJ,K)
            FL(IJ,K,M) = SPOS(IJ,K,M)+CNSN(IJ)*UFAC2(IJ,K)
            SPOS(IJ,K,M) = SPOS(IJ,K,M)*F(IJ,K,M)
            SL(IJ,K,M) = FL(IJ,K,M)*F(IJ,K,M)
          ENDDO
        ENDDO

      ENDDO

      IF (LHOOK) CALL DR_HOOK('SINPUT_JAN',1,ZHOOK_HANDLE)

      END SUBROUTINE SINPUT_JAN
