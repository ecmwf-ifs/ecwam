      SUBROUTINE PARMEAN (IJS, IJL, NPMAX, NPK,
     &                    F,
     &                    EM, THQ, MEANWP)

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

!       *CALL* *PARMEAN (IJS, IJL, NPMAX, NPK,
!                        F,
!                        EM, THQ, MEANWP)
!              *IJS*     - INDEX OF FIRST GRIDPOINT
!              *IJL*     - INDEX OF LAST GRIDPOINT
!              *NPMAX*   - MAXIMUM NUMBER OF PARTITIONS
!              *NPK*     - NUMBER OF PEAKS 
!              *F*       - SPECTRUM.
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
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL, NPMAX
      INTEGER, INTENT(IN), DIMENSION(IJS:IJL) :: NPK
      REAL, INTENT(IN) :: F(NANG,NFRE,NPMAX,IJS:IJL)
      REAL, INTENT(OUT), DIMENSION(IJS:IJL,0:NPMAX) ::  EM, THQ, MEANWP

      INTEGER :: IPK, NPK, K, M
      REAL :: ZHOOK_HANDLE
      REAL    :: F1D(NFRE)
      REAL, DIMENSION(NPK) :: TEMP, EMEAN, SI, CI, FM

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('PARMEAN',0,ZHOOK_HANDLE)
#endif

!     ENE must be initialised and index 0 must stay equal to 0
      DO IP=0,NPMAX
        DO IJ=IJS,IJL
          ENE(IJ,IP)=0.
          DIR(IJ,IP)=0.
          PER(IJ,IP)=0.
        ENDDO
      ENDDO

?????????? to finish
        IF(NPEAK(IJ).GT.0) THEN
        ENDIF


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

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('PARMEAN',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE PARMEAN
