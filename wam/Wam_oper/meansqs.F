      SUBROUTINE MEANSQS(F, IJS, IJL, USNEW)

! ----------------------------------------------------------------------

!**** *MEANSQS* - COMPUTATION OF MEAN SQUARE SLOPE.

!     P.A.E.M. JANSSEN
!     J. BIDLOT  ECMWF  FEBRUARY 1996  MESSAGE PASSING


!*    PURPOSE.
!     --------

!       COMPUTE MEAN SQUARE SLOPE AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *MEANSQS (F, IJS, IJL, USNEW)*
!              *F*   - SPECTRUM.
!              *IJS* - INDEX OF FIRST GRIDPOINT
!              *IJL* - INDEX OF LAST GRIDPOINT
!              *USNEW*     NEW FRICTION VELOCITY IN M/S.

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

      USE YOWPCONS , ONLY : G        ,ZPI      ,EPSMIN
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH
      USE YOWMEAN  , ONLY : FMEAN    ,SMEAN
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSHAL  , ONLY : TFAK     ,INDEP
      USE YOWSTAT  , ONLY : ISHALLO

! ----------------------------------------------------------------------

!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: F

      REAL,DIMENSION(NINF:NSUP) :: USNEW

! ----------------------------------------------------------------------

      PARAMETER (XLAMBDAC=0.0628,SURFT=0.000075)

      DIMENSION TEMP1(IJS:IJL), TEMP2(IJS:IJL)

! ----------------------------------------------------------------------

!*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
!        ------------------------------------------------

      DO IJ=IJS,IJL
        SMEAN(IJ) = EPSMIN
      ENDDO

      FS     = FR(NFRE)      
      XKC    = ZPI/XLAMBDAC
      FC     = SQRT(G*XKC+SURFT*XKC**3)/ZPI
      CONST2 = ZPI**4*FS**5/G**2

! ----------------------------------------------------------------------

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      IF (ISHALLO.EQ.1) THEN

!*    2.1 DEEP WATER INTEGRATION.
!         -----------------------

        DO M=1,NFRE
          FD = DFIM(M)*(ZPI*FR(M))**4/G**2
          DO IJ=IJS,IJL
            TEMP2(IJ) = 0.
          ENDDO
          DO K=1,NANG
            DO IJ=IJS,IJL
              TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
            ENDDO
          ENDDO
          DO IJ=IJS,IJL
            SMEAN(IJ) = SMEAN(IJ)+FD*TEMP2(IJ)
          ENDDO
        ENDDO
!SHALLOW
      ELSE

!*    2.2 SHALLOW WATER INTEGRATION.
!         --------------------------

        DO M=1,NFRE
          DO IJ=IJS,IJL
            TEMP1(IJ) = DFIM(M)*TFAK(INDEP(IJ),M)**2
            TEMP2(IJ) = 0.
          ENDDO
          DO K=1,NANG
            DO IJ=IJS,IJL
              TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
            ENDDO
          ENDDO
          DO IJ=IJS,IJL
            SMEAN(IJ) = SMEAN(IJ)+TEMP1(IJ)*TEMP2(IJ)
          ENDDO
        ENDDO
      ENDIF
!SHALLOW

!*    3. ADD TAIL CORRECTION TO MEAN SQUARE SLOPE.
!        ------------------------------------------

      DO IJ=IJS,IJL      
        CONST1    = LOG(FC/FS)  
        CP        = G/(ZPI*FMEAN(IJ))
        XI        = CP/USNEW(IJ)
        ALPHAP    = MAX(0.21/XI,0.0040)
!        CONST3    = MAX(CONST2*DELTH*TEMP2(IJ),ALPHAP)
!        CONST3    = CONST2*DELTH*TEMP2(IJ)+ALPHAP
        CONST3    = CONST2*DELTH*TEMP2(IJ)
        SMEAN(IJ) = SMEAN(IJ)+CONST1*CONST3
      ENDDO

      RETURN
      END SUBROUTINE MEANSQS
