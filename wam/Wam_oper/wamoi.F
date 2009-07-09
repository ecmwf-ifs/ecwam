      SUBROUTINE WAMOI (NOBS, IOBS4IJ, W, IJS, IJL, NDIM2,
     &                  DIST, XMO, XOI, IG) 

!**** *WAMOI* - CARRIES OUT OPTIMUM INTERPOLATION FOR EACH GRID POINT.

!     J. BIDLOT      ECMWF      JULY 2007 SPLIT FROM *OIFIELD* 

!     PURPOSE.                                                          
!     --------                                                          
!       TO CARRY OUT OPTIMUM INTERPOLATION FOR A GRID POINT.

!**   INTERFACE.                                                        
!     ----------                                                        

!       *CALL* *WAMOI ()*      
!         *NOBS*    INTEGER  NUMBER OF OBSERVATIONS INFLUENCING GRID POINT. 
!         *IOBS4IJ* INTEGER  INDEX TO ALL OBSERVATIONS INFLUENCING GRID POINT.
!         *W*       REAL     ERROR CORRELATION OF ALL OBSERVATIONS
!                            INFLUENCING POINT IJ.
!         *IJS*     INTEGER  INDEX OF FIRST GRIDPOINT.
!         *IJL*     INTEGER  INDEX OF LAST GRIDPOINT.
!         *NDIM2*   INTEGER  MAX DIMENSION OF IOBS4IJ AND W.
!         *DIST*    REAL     CORRELATION DISTANCE.
!         *XMO*     REAL     MODEL FIRST GUESS.
!         *XOI*     REAL     O.I. ON OUTPUT.
!         *IG*      INTEGER  OBSOLETE BLOCK INDEX.

! ----------------------------------------------------------------------

      USE YOWALTAS , ONLY : IJALT    , XLONOBS  , SIGRATIO2, DIFFALTFG
      USE YOWGRID  , ONLY : SINPH    ,COSPH
      USE YOWMAP   , ONLY : IXLG     ,KXLT
      USE YOWPCONS , ONLY : RAD

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IJS, IJL, IJ
      INTEGER :: NDIM2 
      INTEGER :: NOBS(IJS:IJL)
      INTEGER :: IG 
      INTEGER :: IOBS, JOBS, KOBS
      INTEGER :: I1, I2, K1, K2
      INTEGER :: IOBS4IJ(IJS:IJL,NDIM2)

      REAL :: DIST 
      REAL :: XMO(IJS:IJL)
      REAL :: XOI(IJS:IJL)
      REAL :: XNEW
      REAL :: DELLON, COSLON, DOBS, DOBS2
      REAL :: COND
      REAL :: DIFF 
      REAL :: W(IJS:IJL,NDIM2)

      REAL, ALLOCATABLE, DIMENSION(:) :: V
      REAL, ALLOCATABLE, DIMENSION(:,:) :: XM, D, P 

! ----------------------------------------------------------------------

!     LOOP ON EACH GRID POINT


      DO IJ=IJS,IJL

      IF (NOBS(IJ).GT.0) THEN

        XNEW = 0.

        ALLOCATE(D(NOBS(IJ),NOBS(IJ)))
        ALLOCATE(P(NOBS(IJ),NOBS(IJ)))
        ALLOCATE(XM(NOBS(IJ),NOBS(IJ)))

!       DETERMINE OBSERVATION CORRELATION MATRIX  D

        DO IOBS=1,NOBS(IJ)
          DO JOBS=IOBS+1,NOBS(IJ)
            D(IOBS,JOBS)=0.
            D(JOBS,IOBS)=0.
          ENDDO
          D(IOBS,IOBS) = SIGRATIO2(IOBS4IJ(IJ,IOBS))
        ENDDO

!       DETERMINE MODEL CORRELATION MATRIX AT OBSERVATION
!       POINTS.
                                                                        
        DO IOBS=1,NOBS(IJ)
          I1 = IXLG(IJALT(IOBS4IJ(IJ,IOBS),1),IG) 
          K1 = KXLT(IJALT(IOBS4IJ(IJ,IOBS),1),IG) 

          DO JOBS=IOBS+1,NOBS(IJ)

            I2 = IXLG(IJALT(IOBS4IJ(IJ,JOBS),1),IG)
            K2 = KXLT(IJALT(IOBS4IJ(IJ,JOBS),1),IG) 

            IF(IJALT(IOBS4IJ(IJ,IOBS),1).EQ.IJALT(IOBS4IJ(IJ,JOBS),1))
     &      THEN
              DOBS = 0.
            ELSE
              DELLON = XLONOBS(IOBS4IJ(IJ,IOBS))
     &         - XLONOBS(IOBS4IJ(IJ,JOBS))
              COSLON = COS(DELLON*RAD)
              DOBS2  = COSLON*COSPH(K1)*COSPH(K2) +
     &         SINPH(K1)*SINPH(K2)
              DOBS = ACOS(MAX(MIN(DOBS2,1.),-1.))
            ENDIF

            P(IOBS,JOBS) = EXP(-DOBS/DIST)
            P(JOBS,IOBS) = P(IOBS,JOBS)
          ENDDO
          P(IOBS,IOBS) = 1.
        ENDDO

!       SUM BOTH MATRICES
        DO IOBS=1,NOBS(IJ)
          DO JOBS=1,NOBS(IJ)
            XM(IOBS,JOBS)=D(IOBS,JOBS)+P(IOBS,JOBS)
          ENDDO
        ENDDO

        DEALLOCATE(P)
        DEALLOCATE(D)            

!       INVERSE MATRIX D
        COND=0.

        CALL WAM_SYMINV(XM,NOBS(IJ),NOBS(IJ),COND)

        DO IOBS=1,NOBS(IJ)-1
          DO JOBS=IOBS+1,NOBS(IJ)
           XM(IOBS,JOBS) = XM(JOBS,IOBS)
          ENDDO
        ENDDO


!       PRODUCE FIELD AT POINT I,J ACCORDING TO OPTIMUM INTERPOLATION

        DO IOBS=1,NOBS(IJ)
          DIFF = 0.
          DO JOBS=1,NOBS(IJ)
            KOBS=IOBS4IJ(IJ,JOBS)
            DIFF=DIFF+XM(IOBS,JOBS)*DIFFALTFG(KOBS)
          ENDDO
          XNEW=XNEW+W(IJ,IOBS)*DIFF
        ENDDO

        DEALLOCATE(XM)

!       UPDATE FIRST GUESS
        XOI(IJ) = XMO(IJ) + XNEW

      ENDIF
      ENDDO

      RETURN                                                            
      END SUBROUTINE WAMOI
