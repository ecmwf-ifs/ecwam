      SUBROUTINE SINPUT (F, FL, IJS, IJL, THWNEW, USNEW, Z0NEW, SL)
! ----------------------------------------------------------------------

!**** *SINPUT* - COMPUTATION OF INPUT SOURCE FUNCTION.

!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990

!     OPTIMIZED BY : H. GUENTHER

!     MODIFIED BY : J-R BIDLOT NOVEMBER 1995
!                   J-R BIDLOT FEBRUARY 1996-97
!                   J-R BIDLOT FEBRUARY 2001 : MAKE IT FULLY IMPLICIT
!                                              BY ONLY USING NEW
!                                              STRESS AND ROUGHNESS. 

!*    PURPOSE.
!     ---------

!       COMPUTE INPUT SOURCE FUNCTION AND STORE ADDITIVELY INTO NET
!       SOURCE FUNCTION ARRAY, ALSO COMPUTE FUNCTIONAL DERIVATIVE OF
!       INPUT SOURCE FUNCTION. 	SPECIAL FEATURES WERE ADDED TO ALLOW
!       FOR MORE THAN ONE CALL OF SINPUT PER TIME STEP IN ORDER TO
!       LIMIT MEMORY USAGE (SEE ICALL).

!**   INTERFACE.
!     ----------

!     *CALL* *SINPUT (F, FL, IJS, IJL, THWNEW, USNEW, Z0NEW, SL)
!          *F*   - SPECTRUM.
!          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
!          *IJS* - INDEX OF FIRST GRIDPOINT.
!          *IJL* - INDEX OF LAST GRIDPOINT.
!      *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!      *USNEW*     NEW FRICTION VELOCITY IN M/S.
!      *Z0NEW*     ROUGHNESS LENGTH IN M.
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
      REAL,DIMENSION(NINF:NSUP) :: THWNEW,USNEW,Z0NEW

! ----------------------------------------------------------------------

      REAL, DIMENSION(IJS:IJL,NANG) :: TEMP1, UFAC2
      REAL, DIMENSION(IJS:IJL) :: UCN, ZCN, CM
      REAL, DIMENSION(NFRE) :: FAC, CONST

! ----------------------------------------------------------------------

!*    1. PRECALCULATED ANGULAR DEPENDENCE.
!        ---------------------------------

      DO K=1,NANG
        DO IJ=IJS,IJL
          TEMP1(IJ,K) = COS(TH(K)-THWNEW(IJ))
        ENDDO
      ENDDO

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

        DO IJ=IJS,IJL
          UCN(IJ) = USNEW(IJ)*CM(IJ) + ZALP
          ZCN(IJ) = LOG(G*Z0NEW(IJ)*CM(IJ)**2)
        ENDDO

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------

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
            SL(IJ,K,M) = UFAC2(IJ,K)*F(IJ,K,M)
            FL(IJ,K,M) = UFAC2(IJ,K) 
          ENDDO
        ENDDO

      ENDDO

      RETURN
      END SUBROUTINE SINPUT
