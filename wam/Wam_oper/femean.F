      SUBROUTINE FEMEAN (F, IJS, IJL)

! ----------------------------------------------------------------------

!**** *FEMEAN* - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT.

!     S.D. HASSELMANN
!     MODIFIED : P.JANSSEN (INTEGRATION OF F**-4 TAIL)
!     OPTIMIZED BY : L. ZAMBRESKY AND H. GUENTHER


!*    PURPOSE.
!     --------

!       COMPUTE MEAN FREQUENCY AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *FEMEAN (F, IJS, IJL)*
!              *F*   - SPECTRUM.
!              *IJS* - INDEX OF FIRST GRIDPOINT
!              *IJL* - INDEX OF LAST GRIDPOINT

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


      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,EPSMIN
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,FRTAIL
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,AKMEAN
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWSTAT  , ONLY : ISHALLO
      USE YOWSHAL  , ONLY : TFAK     ,INDEP



! ----------------------------------------------------------------------

      DIMENSION F(NINF-1:NSUP,NANG,NFRE)

! ----------------------------------------------------------------------

      DIMENSION TEMP1(IJS:IJL), TEMP2(IJS:IJL)

! ----------------------------------------------------------------------

!*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
!        ------------------------------------------------

      DO IJ=IJS,IJL
        FMEAN(IJ) = EPSMIN
      ENDDO
      DELT25 = FRTAIL*DELTH

! ----------------------------------------------------------------------

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      IF (ISHALLO.EQ.1) THEN

!*    2.1 DEEP WATER INTEGRATION.
!         -----------------------

        DO M=1,NFRE
          FD = DFIM(M)/FR(M)
          DO IJ=IJS,IJL
            TEMP2(IJ) = 0.
          ENDDO
          DO K=1,NANG
            DO IJ=IJS,IJL
              TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
            ENDDO
          ENDDO
          DO IJ=IJS,IJL
            FMEAN(IJ) = FMEAN(IJ)+FD*TEMP2(IJ)
          ENDDO
        ENDDO
!SHALLOW
      ELSE

!*    2.2 SHALLOW WATER INTEGRATION.
!         --------------------------

        DO IJ=IJS,IJL
          AKMEAN(IJ) = EPSMIN
        ENDDO
        DO M=1,NFRE
          FD=DFIM(M)/FR(M)
          DO IJ=IJS,IJL
            TEMP1(IJ) = DFIM(M)/SQRT(TFAK(INDEP(IJ),M))
            TEMP2(IJ) = 0.
          ENDDO
          DO K=1,NANG
            DO IJ=IJS,IJL
              TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
            ENDDO
          ENDDO
          DO IJ=IJS,IJL
            FMEAN(IJ) = FMEAN(IJ)+FD*TEMP2(IJ)
            AKMEAN(IJ) = AKMEAN(IJ)+TEMP1(IJ)*TEMP2(IJ)
          ENDDO
        ENDDO

!        ADD TAIL TO MEAN WAVENUMBER AND NORMALIZE WITH TOTAL ENERGY.

        DEL2 = DELT25*SQRT(G)/ZPI
        DO IJ=IJS,IJL
          AKMEAN(IJ) = AKMEAN(IJ)+DEL2*TEMP2(IJ)
          AKMEAN(IJ) = (EMEAN(IJ)/AKMEAN(IJ))**2
        ENDDO
      ENDIF
!SHALLOW

!*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
!*       NORMALIZE WITH TOTAL ENERGY.
!        ------------------------------------------

      DO IJ=IJS,IJL
        FMEAN(IJ) = FMEAN(IJ)+DELT25*TEMP2(IJ)
        FMEAN(IJ) = EMEAN(IJ)/FMEAN(IJ)
      ENDDO

      RETURN
      END SUBROUTINE FEMEAN
