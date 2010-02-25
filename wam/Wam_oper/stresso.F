      SUBROUTINE STRESSO (F,IJS, IJL, THWNEW,USNEW,Z0NEW,
     &                    ROAIRN, TAUW, SL, MIJ)

! ----------------------------------------------------------------------

!**** *STRESSO* - COMPUTATION OF WAVE STRESS.

!     H. GUNTHER      GKSS/ECMWF  NOVEMBER  1989 CODE MOVED FROM SINPUT.
!     P.A.E.M. JANSSEN      KNMI  AUGUST    1990
!     J. BIDLOT             ECMWF FEBRUARY  1996-97
!     S. ABDALLA            ECMWF OCTOBER   2001 INTRODUCTION OF VARIABLE
!                                                AIR DENSITY
!     J. BIDLOT             ECMWF           2007  ADD MIJ

!*    PURPOSE.
!     --------

!       COMPUTE NORMALIZED WAVE STRESS FROM INPUT SOURCE FUNCTION

!**   INTERFACE.
!     ----------

!       *CALL* *STRESSO (F,IJS,IJL,THWNEW,USNEW,Z0NEW,ROAIRN,TAUW,SL,MIJ)*
!         *F*       - WAVE SPECTRUM.
!         *IJS*     - INDEX OF FIRST GRIDPOINT.
!         *IJL*     - INDEX OF LAST GRIDPOINT.
!         *THWNEW*  - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                     NOTATION (POINTING ANGLE OF WIND VECTOR,
!                     CLOCKWISE FROM NORTH).
!         *USNEW*   - NEW FRICTION VELOCITY IN M/S.
!         *Z0NEW*   - ROUGHNESS LENGTH IN M.
!         *ROAIRN*  - AIR DENSITY IN KG/M3.
!         *TAUW*    - WAVE STRESS IN (M/S)**2
!         *SL*      - TOTAL SOURCE FUNCTION ARRAY.
!         *MIJ*     - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.

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
     &            COSTH    ,SINTH    ,FR5
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER  ,YINVEPS
      USE YOWTABL  , ONLY : IUSTAR   ,IALPHA   ,EPS1     ,TAUHFT   ,
     &            DELUST   ,DELALP

! ----------------------------------------------------------------------

!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      INTEGER :: MIJ(IJS:IJL)

      REAL,DIMENSION(IJS:IJL,NANG,NFRE) :: F,SL
      REAL,DIMENSION(IJS:IJL) :: THWNEW, USNEW, Z0NEW, ROAIRN, TAUW
      REAL :: CONST0(IJS:IJL)

! ----------------------------------------------------------------------

      REAL :: CONST, DFIMLOC, CNST
      REAL :: ZPIROFR(NFRE)
      REAL :: CONSTF(IJS:IJL,NFRE)
      REAL, DIMENSION(IJS:IJL) :: TAUHF, TEMP, XSTRESS, YSTRESS, UST2
! ----------------------------------------------------------------------

!*    1. PRECOMPUTE FREQUENCY SCALING.
!        -----------------------------

      CONST = DELTH*(ZPI)**4/G**2

      DO IJ=IJS,IJL
        CONST0(IJ)  = CONST*FR5(MIJ(IJ))
      ENDDO

      DO M=1,NFRE
        ZPIROFR(M) =ZPI*ROWATER*FR(M)
      ENDDO

      DO IJ=IJS,IJL
        DO M=1,MIJ(IJ)-1
          CONSTF(IJ,M) =ZPIROFR(M)*DFIM(M)
        ENDDO
        DFIMLOC=0.5*DELTH*(FR(MIJ(IJ))-FR(MIJ(IJ)-1))
        CONSTF(IJ,MIJ(IJ)) =ZPIROFR(MIJ(IJ))*DFIMLOC
      ENDDO

!*    2. COMPUTE WAVE STRESS OF ACTUEL BLOCK.
!        ------------------------------------

!*    2.2 INTEGRATE INPUT SOURCE FUNCTION OVER FREQUENCY AND DIRECTIONS.
!         --------------------------------------------------------------

      DO IJ=IJS,IJL
        XSTRESS(IJ)=0.
        YSTRESS(IJ)=0.
        DO M=1,MIJ(IJ)
          DO K=1,NANG
            CNST=SL(IJ,K,M)*CONSTF(IJ,M)
            XSTRESS(IJ)=XSTRESS(IJ)+CNST*SINTH(K)
            YSTRESS(IJ)=YSTRESS(IJ)+CNST*COSTH(K)
          ENDDO
        ENDDO
        XSTRESS(IJ)=XSTRESS(IJ)/MAX(ROAIRN(IJ),1.)
        YSTRESS(IJ)=YSTRESS(IJ)/MAX(ROAIRN(IJ),1.)
      ENDDO


!*    2.3 CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
!     ----------------------------------------------------

      K=1
      DO IJ=IJS,IJL
        COSW     = MAX(COS(TH(K)-THWNEW(IJ)),0.)
        TEMP(IJ) = F(IJ,K,MIJ(IJ))*COSW**3
      ENDDO

      DO K=2,NANG
        DO IJ=IJS,IJL
          COSW     = MAX(COS(TH(K)-THWNEW(IJ)),0.)
          TEMP(IJ) = TEMP(IJ)+F(IJ,K,MIJ(IJ))*COSW**3
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        UST   = MAX(USNEW(IJ),0.000001)
        UST2(IJ) = UST**2
        XI    = UST / DELUST
        XI    = MIN(REAL(IUSTAR),XI)
        I     = MIN (IUSTAR-1, INT(XI))
        I     = MAX (0, I)
        DELI1 = MIN (1. ,XI-REAL(I))
        DELI2   = 1. - DELI1

        XJ    = (G*Z0NEW(IJ)/UST2(IJ)-ALPHA) / DELALP
        XJ    = MIN(REAL(IALPHA),XJ)
        J     = MIN (IALPHA-1, INT(XJ))
        J     = MAX (0, J)
        DELJ1 = MAX(MIN (1. ,XJ-REAL(J)),0.)
        DELJ2   = 1. - DELJ1

        TAU1 = ( TAUHFT(I  ,J  ,MIJ(IJ))*DELI2 +
     &           TAUHFT(I+1,J  ,MIJ(IJ))*DELI1 )*DELJ2
     &       + ( TAUHFT(I  ,J+1,MIJ(IJ))*DELI2 +
     &           TAUHFT(I+1,J+1,MIJ(IJ))*DELI1 )*DELJ1

        TAUHF(IJ) = CONST0(IJ)*TEMP(IJ)*UST2(IJ)*TAU1
      ENDDO

      DO IJ=IJS,IJL
        XSTRESS(IJ) = XSTRESS(IJ)+TAUHF(IJ)*SIN(THWNEW(IJ))
        YSTRESS(IJ) = YSTRESS(IJ)+TAUHF(IJ)*COS(THWNEW(IJ))
        TAUW(IJ) = SQRT(XSTRESS(IJ)**2+YSTRESS(IJ)**2)

        TAUW(IJ) = MIN(TAUW(IJ),UST2(IJ)-EPS1)
        TAUW(IJ) = MAX(TAUW(IJ),0.)
      ENDDO

      RETURN
      END SUBROUTINE STRESSO
