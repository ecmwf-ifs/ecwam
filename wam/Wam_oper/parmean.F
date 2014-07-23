      SUBROUTINE PARMEAN (F, NPK, EM, THQ, MEANWP)

! ----------------------------------------------------------------------

!**** *PARMEAN* - COMPUTATION OF TOTAL ENERGY, MEAN DIRECTION,
!                  AND MEAN PERIOD BASED ON THE FIRST MOMENT

!      J-R BIDLOT    ECMWF     MARCH 2000
!      D PETTENUZZO  MAY 2012 MERGED 3 ROUTINES FOR PAR COMPUTATION

!*    PURPOSE.
!     --------

!       COMPUTE TOTAL ENERGY, MEAN DIRECTION AND MEAN PERIOD.

!**   INTERFACE.
!     ----------

!       *CALL* *PARMEAN (F,NPK,EM,THQ,MEANWP)*
!              *F*       - SPECTRUM.
!              *NPK*     - NUMBER OF PEAKS 
!              *EM*      - MEAN WAVE ENERGY FOR THE WAVE SYSTEM 
!              *THQ*     - MEAN WAVE DIRECTION
!              *MEANWP*  - MEAN PERIOD BASED ON -1 MOMENT.

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

      USE YOWFRED  ,ONLY : FR , DFIM , DFIMOFR, COSTH, SINTH
      USE YOWPARAM ,ONLY : NANG, NFRE
      USE YOWPCONS ,ONLY : EPSMIN, ZPI

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IPK, NPK, K, M
      REAL    :: F1D(NFRE)
      REAL, DIMENSION(NPK) :: TEMP, EMEAN, MEANWP, THQ, SI, CI, EM, FM
      REAL    :: F(NANG,NFRE,NPK)

! ----------------------------------------------------------------------


      DO IPK=1,NPK

        EM(IPK) = EPSMIN
        FM(IPK) = EPSMIN
        DO M=1,NFRE
          F1D(M) = F(1,M,IPK)
          DO K=2,NANG
            F1D(M) = F1D(M)+F(K,M,IPK)
          ENDDO
          EM(IPK) = EM(IPK)+F1D(M)*DFIM(M)
          FM(IPK) = FM(IPK)+F1D(M)*DFIMOFR(M)
        ENDDO

        IF(EM(IPK).GT.EPSMIN) THEN
          MEANWP(IPK) = FM(IPK)/EM(IPK)
        ELSE
          MEANWP(IPK) = 0.
        ENDIF

      ENDDO

! ----------------------------------------------------------------------

!*    MEAN WAVE DIRECTION
!     -------------------

      DO IPK=1,NPK
        SI(IPK) = 0.
        CI(IPK) = 0.

        DO K=1,NANG
          TEMP(IPK) = F(K,1,IPK)*DFIM(1)
          DO M=2,NFRE
            TEMP(IPK) = TEMP(IPK)+F(K,M,IPK)*DFIM(M)
          ENDDO
          SI(IPK) = SI(IPK)+SINTH(K)*TEMP(IPK)
          CI(IPK) = CI(IPK)+COSTH(K)*TEMP(IPK)
        ENDDO


        IF (CI(IPK).EQ.0.) CI(IPK) = EPSMIN
        THQ(IPK) = ATAN2(SI(IPK),CI(IPK))
        IF (THQ(IPK).LT.0.) THQ(IPK) = THQ(IPK) + ZPI
      ENDDO

! ----------------------------------------------------------------------

      RETURN
      END SUBROUTINE PARMEAN
