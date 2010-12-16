      SUBROUTINE MEANSQS(F, IJS, IJL, USNEW, FM, SM)

! ----------------------------------------------------------------------

!**** *MEANSQS* - COMPUTATION OF MEAN SQUARE SLOPE.

!     P.A.E.M. JANSSEN
!     J. BIDLOT  ECMWF  FEBRUARY 1996  MESSAGE PASSING


!*    PURPOSE.
!     --------

!       COMPUTE MEAN SQUARE SLOPE AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *MEANSQS (F, IJS, IJL, USNEW, FM, SM)*
!              *F*   - SPECTRUM.
!              *IJS* - INDEX OF FIRST GRIDPOINT
!              *IJL* - INDEX OF LAST GRIDPOINT
!              *USNEW*     NEW FRICTION VELOCITY IN M/S (INPUT).
!              *FM*  - MEAN WAVE FREQUENCY (INPUT).
!              *SM*  - MEAN SQUARE SLOPE (OUTPUT).

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
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSHAL  , ONLY : TFAK     ,INDEP
      USE YOWSTAT  , ONLY : ISHALLO

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL, PARAMETER :: XLAMBDAC=0.0628
      REAL, PARAMETER :: SURFT=0.000075

      INTEGER :: IJ,M,K,IJS,IJL
      REAL :: FS, XKC, FC, CONST1, CONST2, CP, XI, ALPHAP, CONST3 
      REAL :: F(IJS:IJL,NANG,NFRE), FD(NFRE)
      REAL,DIMENSION(IJS:IJL) :: USNEW
      REAL, DIMENSION(IJS:IJL) ::  FM, SM, TEMP1, TEMP2


! ----------------------------------------------------------------------


!*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
!        ------------------------------------------------

      DO IJ=IJS,IJL
        SM(IJ) = EPSMIN
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
          FD(M) = DFIM(M)*(ZPI*FR(M))**4/G**2
        ENDDO

        DO M=1,NFRE
          DO IJ=IJS,IJL
            TEMP2(IJ) = 0.
          ENDDO
          DO K=1,NANG
            DO IJ=IJS,IJL
              TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
            ENDDO
          ENDDO
          DO IJ=IJS,IJL
            SM(IJ) = SM(IJ)+FD(M)*TEMP2(IJ)
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
            SM(IJ) = SM(IJ)+TEMP1(IJ)*TEMP2(IJ)
          ENDDO
        ENDDO
      ENDIF
!SHALLOW

!*    3. ADD TAIL CORRECTION TO MEAN SQUARE SLOPE.
!        ------------------------------------------

      DO IJ=IJS,IJL      
        CONST1    = 0.*LOG(FC/FS)  
        CP        = G/(ZPI*FM(IJ))
        XI        = CP/USNEW(IJ)
        ALPHAP    = MAX(0.21/XI,0.0040)
        CONST3    = CONST2*DELTH*TEMP2(IJ)
        SM(IJ) = SM(IJ)+CONST1*CONST3
      ENDDO

      RETURN
      END SUBROUTINE MEANSQS
