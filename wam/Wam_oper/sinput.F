      SUBROUTINE SINPUT (F, FL, IJS, IJL, IG,THWOLD,USOLD,Z0OLD,
     &                   THWNEW,USNEW,Z0NEW,SL,ICALL,NCALL)
! ----------------------------------------------------------------------

!**** *SINPUT* - COMPUTATION OF INPUT SOURCE FUNCTION.

!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990

!     OPTIMIZED BY : H. GUENTHER

!     MODIFIED BY : J-R BIDLOT NOVEMBER 1995
!                   J-R BIDLOT FEBRUARY 1996-97
!                   J-R BIDLOT FEBRUARY 1999 : INTRODUCE ICALL AND NCALL

!*    PURPOSE.
!     ---------

!       COMPUTE INPUT SOURCE FUNCTION AND STORE ADDITIVELY INTO NET
!       SOURCE FUNCTION ARRAY, ALSO COMPUTE FUNCTIONAL DERIVATIVE OF
!       INPUT SOURCE FUNCTION. 	SPECIAL FEATURES WERE ADDED TO ALLOW
!       FOR MORE THAN ONE CALL OF SINPUT PER TIME STEP IN ORDER TO
!       LIMIT MEMORY USAGE (SEE ICALL).

!**   INTERFACE.
!     ----------

!     *CALL* *SINPUT (F, FL, IJS, IJL, IG,THWOLD,USOLD,Z0OLD,
!    1                   THWNEW,USNEW,Z0NEW,SL,ICALL,NCALL)
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
!      *ICALL*     CALLING INDEX 1<= ICALL <= NCALL.
!      *NCALL*     TOTAL NUMBER OF CALLS PER TIME STEP.

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

      USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,XKAPPA
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO
      USE YOWPCONS , ONLY : G        ,ZPI      ,YEPS
      USE YOWSHAL  , ONLY : TFAK     ,INDEP
      USE YOWSTAT  , ONLY : ISHALLO

! ----------------------------------------------------------------------

!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: F,FL,SL
      REAL,DIMENSION(NINF:NSUP,NBLO) :: THWOLD,USOLD,Z0OLD
      REAL,DIMENSION(NINF:NSUP) :: THWNEW,USNEW,Z0NEW

! ----------------------------------------------------------------------

      REAL, DIMENSION(IJS:IJL,NANG) :: TEMP, TEMP1, UFAC2
      REAL, DIMENSION(IJS:IJL) :: UCO, ZCO, UCN, ZCN, CM
      REAL, DIMENSION(NFRE) :: FAC, CONST

! ----------------------------------------------------------------------

!*    1. PRECALCULATED ANGULAR DEPENDENCE.
!        ---------------------------------

!!! temp and temp1 will be reused further down and will no contain
!!! the cosine of directional difference
      IF(ICALL.EQ.1) THEN
        DO K=1,NANG
          DO IJ=IJS,IJL
            TEMP (IJ,K) = COS(TH(K)-THWOLD(IJ,IG))
            TEMP1(IJ,K) = COS(TH(K)-THWNEW(IJ))
          ENDDO
        ENDDO
      ELSE
        DO K=1,NANG
          DO IJ=IJS,IJL
            TEMP1(IJ,K) = COS(TH(K)-THWNEW(IJ))
          ENDDO
        ENDDO

      ENDIF

! ----------------------------------------------------------------------

!*    2. LOOP OVER FREQUENCIES.
!        ----------------------

      CONST1  = YEPS*BETAMAX/XKAPPA**2
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

        IF(ICALL.EQ.1) THEN
          DO IJ=IJS,IJL
            UCO(IJ) = USOLD(IJ,IG)*CM(IJ) + ZALP
            ZCO(IJ) = LOG(G*Z0OLD(IJ,IG)*CM(IJ)**2)
          ENDDO
        ENDIF

        DO IJ=IJS,IJL
          UCN(IJ) = USNEW(IJ)*CM(IJ) + ZALP
          ZCN(IJ) = LOG(G*Z0NEW(IJ)*CM(IJ)**2)
        ENDDO

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------

        IF(ICALL.EQ.1) THEN
          DO K=1,NANG
            DO IJ=IJS,IJL
              IF (TEMP(IJ,K).GT.0.01) THEN
                X  = TEMP(IJ,K)*UCO(IJ)
                ZLOG = ZCO(IJ) + XKAPPA/X
                IF (ZLOG.LT.0.) THEN
                  ZLOG2X=ZLOG*ZLOG*X
                  FL(IJ,K,M) = CONST(M)*EXP(ZLOG)*ZLOG2X*ZLOG2X
                ELSE
                  FL(IJ,K,M) = 0.
                ENDIF
              ELSE
                FL(IJ,K,M) = 0.
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        DO K=1,NANG
          DO IJ=IJS,IJL
            IF (TEMP1(IJ,K).GT.0.01) THEN
              X    = TEMP1(IJ,K)*UCN(IJ)
              ZLOG = ZCN(IJ) + XKAPPA/X
              IF (ZLOG.LT.0.) THEN
                ZLOG2X=ZLOG*ZLOG*X
                UFAC2(IJ,K) = CONST(M)*EXP(ZLOG)*ZLOG2X*ZLOG2X
              ELSE
                UFAC2(IJ,K) = 0.
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
            SL(IJ,K,M) = 0.5*(FL(IJ,K,M)+UFAC2(IJ,K))*F(IJ,K,M)
          ENDDO
        ENDDO

        IF(ICALL.EQ.NCALL) THEN
          DO K=1,NANG
            DO IJ=IJS,IJL
              FL(IJ,K,M) = UFAC2(IJ,K) 
            ENDDO
          ENDDO
        ENDIF

      ENDDO

      RETURN
      END SUBROUTINE SINPUT
