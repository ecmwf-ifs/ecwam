      SUBROUTINE SINPUT (F, FL, IJS, IJL, IG,THWOLD,USOLD,Z0OLD,
     &                   THWNEW,USNEW,Z0NEW,SL)
! ----------------------------------------------------------------------

!**** *SINPUT* - COMPUTATION OF INPUT SOURCE FUNCTION.

!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990

!     OPTIMIZED BY : H. GUENTHER

!     MODIFIED BY : J-R BIDLOT NOVEMBER 1995
!                   J-R BIDLOT FEBRUARY 1996-97

!*    PURPOSE.
!     ---------

!       COMPUTE INPUT SOURCE FUNCTION AND STORE ADDITIVELY INTO NET
!       SOURCE FUNCTION ARRAY, ALSO COMPUTE FUNCTIONAL DERIVATIVE OF
!       INPUT SOURCE FUNCTION.

!**   INTERFACE.
!     ----------

!     *CALL* *SINPUT (F, FL, IJS, IJL, IG,THWOLD,USOLD,Z0OLD,
!    1                   THWNEW,USNEW,Z0NEW,SL)
!          *F*   - SPECTRUM.
!          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
!          *IJS* - INDEX OF FIRST GRIDPOINT.
!          *IJL* - INDEX OF LAST GRIDPOINT.
!          *IG*  - BLOCK NUMBER.
!      *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!      *THWOLD*    INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
!                  WIND VELOCITY.
!      *USNEW*     NEW FRICTION VELOCITY IN M/S.
!      *USOLD*     INTERMEDIATE STORAGE OF MODULUS OF FRICTION
!                  VELOCITY.
!      *Z0NEW*     ROUGHNESS LENGTH IN M.
!      *Z0OLD*     INTERMEDIATE STORAGE OF ROUGHNESS LENGTH IN
!                  M.
!      *SL*        REAL      TOTAL SOURCE FUNCTION ARRAY.

!     METHOD.
!     -------

!       SEE REFERENCE.

!     EXTERNALS.
!     ----------

!       NONE.

!     MODIFICATION
!     ------------

!     REMOVAL OF CALL TO CRAY SPECIFIC FUNCTIONS EXPHF AND ALOGHF
!     BY THEIR STANDARD FORTRAN EQUIVALENT EXP and ALOGHF

!     REFERENCE.
!     ----------

!       P. JANSSEN, J.P.O., 1989.
!       P. JANSSEN, J.P.O., 1991

! ----------------------------------------------------------------------


      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO
      USE YOWPCONS , ONLY : G        ,ZPI      ,YEPS
      USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,XKAPPA
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWSHAL  , ONLY : TFAK     ,INDEP
      USE YOWSTAT  , ONLY : ISHALLO



! ----------------------------------------------------------------------

!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: F,FL,SL
      REAL,DIMENSION(NINF:NSUP,NBLO) :: THWOLD,USOLD,Z0OLD
      REAL,DIMENSION(NINF:NSUP) :: THWNEW,USNEW,Z0NEW


! ----------------------------------------------------------------------

      DIMENSION TEMP(IJS:IJL,NANG), TEMP1(IJS:IJL,NANG),
     &          UCO(IJS:IJL), ZCO(IJS:IJL), UCN(IJS:IJL), ZCN(IJS:IJL),
     &          UFAC1(IJS:IJL), UFAC2(IJS:IJL), CM(IJS:IJL)

! ----------------------------------------------------------------------

!*    1. PRECALCULATED ANGULAR DEPENDENCE.
!        ---------------------------------

      DO K=1,NANG
        TKD=TH(K)
        DO IJ=IJS,IJL
          TEMP (IJ,K) = COS(TKD-THWOLD(IJ,IG))
          TEMP1(IJ,K) = COS(TKD-THWNEW(IJ))
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    2. LOOP OVER FREQUENCIES.
!        ----------------------

      CONST1  = YEPS*BETAMAX/XKAPPA**2
      DO M=1,NFRE
        FAC = ZPI*FR(M)
        CONST=FAC*CONST1

!*      INVERSE OF PHASE VELOCITIES.
!       ----------------------------

        IF (ISHALLO.EQ.1) THEN
          DO IJ=IJS,IJL
            CM(IJ) = FAC/G
          ENDDO
        ELSE
          DO IJ=IJS,IJL
            CM(IJ) = TFAK(INDEP(IJ),M)/FAC
          ENDDO
        ENDIF

!*      PRECALCULATE FREQUENCY DEPENDENCE.
!       ----------------------------------

        DO IJ=IJS,IJL
          UCO(IJ) = USOLD(IJ,IG)*CM(IJ) + ZALP
          ZCO(IJ) = LOG(G*Z0OLD(IJ,IG)*CM(IJ)**2)
          UCN(IJ) = USNEW(IJ)*CM(IJ) + ZALP
          ZCN(IJ) = LOG(G*Z0NEW(IJ)*CM(IJ)**2)
        ENDDO

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------

        DO K=1,NANG
          DO IJ=IJS,IJL
            UFAC1(IJ) = 0.
            UFAC2(IJ) = 0.
          ENDDO
          DO IJ=IJS,IJL
            IF (TEMP(IJ,K).GT.0.01) THEN
              X    = TEMP(IJ,K)*UCO(IJ)
              ZARG = XKAPPA/X
              ZLOG = ZCO(IJ) + ZARG
              IF (ZLOG.LT.0.) THEN
                ZLOG2X=ZLOG*ZLOG*X
                UFAC1(IJ) = CONST*EXP(ZLOG)*ZLOG2X*ZLOG2X
              ENDIF
            ENDIF
          ENDDO

          DO IJ=IJS,IJL
            IF (TEMP1(IJ,K).GT.0.01) THEN
              X    = TEMP1(IJ,K)*UCN(IJ)
              ZARG = XKAPPA/X
              ZLOG = ZCN(IJ) + ZARG
              IF (ZLOG.LT.0.) THEN
                ZLOG2X=ZLOG*ZLOG*X
                UFAC2(IJ) = CONST*EXP(ZLOG)*ZLOG2X*ZLOG2X
              ENDIF
            ENDIF
          ENDDO

!*    2.2 ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
!         ------------------------------------------------

          DO IJ=IJS,IJL
            SL(IJ,K,M) = 0.5*(UFAC1(IJ)+UFAC2(IJ))*F(IJ,K,M)
            FL(IJ,K,M) = UFAC2(IJ)
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE SINPUT
