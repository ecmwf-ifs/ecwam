      SUBROUTINE STRESSO (F,IJS, IJL, IG,THWNEW,USNEW,Z0NEW,TAUW,SL)

! ----------------------------------------------------------------------

!**** *STRESSO* - COMPUTATION OF WAVE STRESS.

!     H. GUNTHER      GKSS/ECMWF  NOVEMBER  1989 CODE MOVED FROM SINPUT.
!     P.A.E.M. JANSSEN      KNMI  AUGUST    1990
!     J. BIDLOT             ECMWF FEBRUARY  1996-97

!*    PURPOSE.
!     --------

!       COMPUTE NORMALIZED WAVE STRESS FROM INPUT SOURCE FUNCTION

!**   INTERFACE.
!     ----------

!       *CALL* *STRESSO (F,IJS, IJL, IG,THWNEW,USNEW,Z0NEW,TAUW,SL)*
!          *F*   - WAVE SPECTRUM.
!          *IJS* - INDEX OF FIRST GRIDPOINT.
!          *IJL* - INDEX OF LAST GRIDPOINT.
!          *IG*  - ACTUAL BLOCK NUMBER.
!      *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!      *USNEW*     NEW FRICTION VELOCITY IN M/S.
!      *Z0NEW*     ROUGHNESS LENGTH IN M.
!      *TAUW*      WAVE STRESS IN (M/S)**2
!      *SL*        REAL      TOTAL SOURCE FUNCTION ARRAY.

!     METHOD.
!     -------

!       THE INPUT SOURCE FUNCTION IS INTEGRATED OVER FREQUENCY
!       AND DIRECTIONS.
!       BECAUSE ARRAY *SL* IS USED, ONLY THE INPUT SOURCE
!       HAS TO BE STORED IN *SL* (CALL FIRST SINPUT, THEN
!       STRESSO, AND THEN THE REST OF THE SOURCE FUNCTIONS)

!     EXTERNALS.
!     -----------

!       NONE.

!     REFERENCE.
!     ----------

!       R SNYDER ET AL,1981.
!       G. KOMEN, S. HASSELMANN AND K. HASSELMANN, JPO, 1984.
!       P. JANSSEN, JPO, 1985

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : ALPHA
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,TH       ,
     &            COSTH    ,SINTH
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO
      USE YOWPCONS , ONLY : G        ,ZPI      ,YINVEPS
      USE YOWTABL  , ONLY : IUSTAR   ,IALPHA   ,TAUHFT   ,DELUST   ,
     &            DELALP

! ----------------------------------------------------------------------

!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: F,SL
      REAL,DIMENSION(NINF:NSUP,NBLO) :: TAUW
      REAL,DIMENSION(NINF:NSUP) :: THWNEW,USNEW,Z0NEW

! ----------------------------------------------------------------------

      DIMENSION CONSTF(NFRE), TAUHF(IJS:IJL), TEMP(IJS:IJL),
     &          XSTRESS(IJS:IJL), YSTRESS(IJS:IJL)
! ----------------------------------------------------------------------

!*    1. PRECOMPUTE FREQUENCY SCALING.
!        -----------------------------

      CONST0  = DELTH*(ZPI)**4*FR(NFRE)**5/G**2

      DO M=1,NFRE
        CONSTF(M) =ZPI*YINVEPS*FR(M)*DFIM(M)
      ENDDO

!*    2. COMPUTE WAVE STRESS OF ACTUEL BLOCK.
!        ------------------------------------

!*    2.1 PRESET STRESS ARRAYS.
!         ---------------------

      DO IJ=IJS,IJL
        XSTRESS(IJ)=0.
        YSTRESS(IJ)=0.
      ENDDO

!*    2.2 INTEGRATE INPUT SOURCE FUNCTION OVER FREQUENCY AND DIRECTIONS.
!         --------------------------------------------------------------

      DO K=1,NANG
        DO M=1,NFRE
          CONST1=CONSTF(M)*SINTH(K)
          CONST2=CONSTF(M)*COSTH(K)
          DO IJ=IJS,IJL
            XSTRESS(IJ)=XSTRESS(IJ)+SL(IJ,K,M)*CONST1
            YSTRESS(IJ)=YSTRESS(IJ)+SL(IJ,K,M)*CONST2
          ENDDO
        ENDDO
      ENDDO

!*    2.3 CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
!     ----------------------------------------------------

      DO IJ=IJS,IJL
        TEMP(IJ) = 0.
      ENDDO

      DO K=1,NANG
        TKD = TH(K)
        DO IJ=IJS,IJL
          COSW     = MAX(COS(TKD-THWNEW(IJ)),0.)
          TEMP(IJ) = TEMP(IJ)+F(IJ,K,NFRE)*COSW**3
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        UST   = MAX(USNEW(IJ),0.000001)
        XI    = UST / DELUST
        XI    = MIN(REAL(IUSTAR),XI)
        I     = MIN (IUSTAR-1, INT(XI))
        I     = MAX (0, I)
        DELI1 = MIN (1. ,XI-REAL(I))
        DELI2   = 1. - DELI1

        XJ    = (G*Z0NEW(IJ)/UST**2-ALPHA) / DELALP
        XJ    = MIN(REAL(IALPHA),XJ)
        J     = MIN (IALPHA-1, INT(XJ))
        J     = MAX (0, J)
        DELJ1 = MIN (1. ,XJ-REAL(J))
        DELJ2   = 1. - DELJ1

        TAU1 = (TAUHFT(I,J  )*DELI2 + TAUHFT(I+1,J  )*DELI1)*DELJ2
     &   + (TAUHFT(I,J+1)*DELI2 + TAUHFT(I+1,J+1)*DELI1)*DELJ1

        TAUHF(IJ) = CONST0*TEMP(IJ)*UST**2*TAU1
      ENDDO

      DO IJ=IJS,IJL
        XSTRESS(IJ) = XSTRESS(IJ)+TAUHF(IJ)*SIN(THWNEW(IJ))
        YSTRESS(IJ) = YSTRESS(IJ)+TAUHF(IJ)*COS(THWNEW(IJ))
        TAUW(IJ,IG) = SQRT(XSTRESS(IJ)**2+YSTRESS(IJ)**2)
      ENDDO

      RETURN
      END SUBROUTINE STRESSO
