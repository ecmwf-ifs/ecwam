SUBROUTINE WAMOI (NOBS, IOBS4IJ, W, KIJS, KIJL, NDIM2, &
 &                BLK2GLO,                             &
 &                DIST, XMO, XOI) 

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
!         *KIJS*    INTEGER  INDEX OF FIRST GRIDPOINT.
!         *KIJL*    INTEGER  INDEX OF LAST GRIDPOINT.
!         *NDIM2*   INTEGER  MAX DIMENSION OF IOBS4IJ AND W.
!         *BLK2GLO*          BLOCK TO GRID TRANSFORMATION
!         *DIST*    REAL     CORRELATION DISTANCE.
!         *XMO*     REAL     MODEL FIRST GUESS.
!         *XOI*     REAL     O.I. ON OUTPUT.

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO

      USE YOWALTAS , ONLY : IJALT    , XLONOBS  , SIGRATIO2, DIFFALTFG, &
     &                 ALTUNDATA
      USE YOWGRID  , ONLY : SINPH    ,COSPH
      USE YOWPARAM , ONLY : LLUNSTR
      USE YOWPCONS , ONLY : RAD
      USE YOWSPHERE, ONLY : SPHERICAL_COORDINATE_DISTANCE

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

#include "wam_syminv.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: NOBS(KIJS:KIJL)
      INTEGER(KIND=JWIM), INTENT(IN) :: IOBS4IJ(KIJS:KIJL,NDIM2)
      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM), INTENT(IN) :: NDIM2
      TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO
      REAL(KIND=JWRB), INTENT(IN) :: W(KIJS:KIJL,NDIM2)
      REAL(KIND=JWRB), INTENT(IN) :: DIST(KIJS:KIJL)
      REAL(KIND=JWRB), INTENT(IN) :: XMO(KIJS:KIJL)
      REAL(KIND=JWRB), INTENT(OUT) :: XOI(KIJS:KIJL)


      INTEGER(KIND=JWIM) :: IJ, I1, I2, K1, K2
      INTEGER(KIND=JWIM) :: IOBS, JOBS, KOBS
      INTEGER(KIND=JWIM) :: IOBS_LOC, JOBS_LOC

      REAL(KIND=JWRB) :: XNEW
      REAL(KIND=JWRB) :: DELLON, COSLON, DOBS, DOBS2
      REAL(KIND=JWRB) :: COND
      REAL(KIND=JWRB) :: DIFF 
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:) :: V
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:) :: XM, D, P 

      REAL(KIND=JWRU)  :: XLONID, XLATID, XLONJD, XLATJD, DISTD

! ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WAMOI',0,ZHOOK_HANDLE)

!     LOOP ON EACH GRID POINT

      DO IJ=KIJS,KIJL
      IF (NOBS(IJ) > 0) THEN

        XNEW = 0._JWRB

        ALLOCATE(D(NOBS(IJ),NOBS(IJ)))
        ALLOCATE(P(NOBS(IJ),NOBS(IJ)))
        ALLOCATE(XM(NOBS(IJ),NOBS(IJ)))

!       DETERMINE OBSERVATION CORRELATION MATRIX  D

        DO IOBS=1,NOBS(IJ)
          IOBS_LOC = IOBS4IJ(IJ,IOBS)
          DO JOBS=IOBS+1,NOBS(IJ)
!!            JOBS_LOC = IOBS4IJ(IJ,JOBS)
            D(IOBS,JOBS)=0._JWRB
            D(JOBS,IOBS)=0._JWRB
          ENDDO
          D(IOBS,IOBS) = SIGRATIO2(IOBS_LOC)
        ENDDO

!       DETERMINE MODEL CORRELATION MATRIX AT OBSERVATION
!       POINTS.
 
        IF (LLUNSTR) THEN
          DO IOBS=1,NOBS(IJ)
            IOBS_LOC = IOBS4IJ(IJ,IOBS)
            XLATID = ALTUNDATA(IOBS_LOC,1)
            XLONID = ALTUNDATA(IOBS_LOC,2)
            DO JOBS=IOBS+1,NOBS(IJ)
              JOBS_LOC = IOBS4IJ(IJ,JOBS)
              IF (IJALT(IOBS_LOC,1) == IJALT(JOBS_LOC,1)) THEN
                DOBS = 0._JWRB
              ELSE
                XLATJD = ALTUNDATA(JOBS_LOC,1)
                XLONJD = ALTUNDATA(JOBS_LOC,2)
                CALL SPHERICAL_COORDINATE_DISTANCE(XLONID,XLONJD,XLATID,XLATJD,DISTD)
                DOBS = REAL(DISTD,JWRB)
              ENDIF

              P(IOBS,JOBS) = EXP(-DOBS/DIST(IJ))
              P(JOBS,IOBS) = P(IOBS,JOBS)
            ENDDO
            P(IOBS,IOBS) = 1._JWRB
          ENDDO
        ELSE
          DO IOBS=1,NOBS(IJ)
            IOBS_LOC = IOBS4IJ(IJ,IOBS)
            I1 = BLK2GLO%IXLG(IJALT(IOBS_LOC,1))
            K1 = BLK2GLO%KXLT(IJALT(IOBS_LOC,1))
            DO JOBS=IOBS+1,NOBS(IJ)
              JOBS_LOC = IOBS4IJ(IJ,JOBS)
              I2 = BLK2GLO%IXLG(IJALT(JOBS_LOC,1))
              K2 = BLK2GLO%KXLT(IJALT(JOBS_LOC,1))

              IF (IJALT(IOBS_LOC,1) == IJALT(JOBS_LOC,1)) THEN
                DOBS = 0._JWRB
              ELSE
                DELLON = XLONOBS(IOBS_LOC)- XLONOBS(JOBS_LOC)
                COSLON = COS(DELLON*RAD)
                DOBS2  = COSLON*COSPH(K1)*COSPH(K2) +  SINPH(K1)*SINPH(K2)
                DOBS = ACOS(MAX(MIN(DOBS2,1.),-1.0_JWRB))
              ENDIF

              P(IOBS,JOBS) = EXP(-DOBS/DIST(IJ))
              P(JOBS,IOBS) = P(IOBS,JOBS)
            ENDDO
            P(IOBS,IOBS) = 1._JWRB
          ENDDO
        ENDIF

!       SUM BOTH MATRICES
        DO IOBS=1,NOBS(IJ)
          DO JOBS=1,NOBS(IJ)
            XM(IOBS,JOBS)=D(IOBS,JOBS)+P(IOBS,JOBS)
          ENDDO
        ENDDO

        DEALLOCATE(P)
        DEALLOCATE(D)            

!       INVERSE MATRIX D
        COND=0._JWRB

        CALL WAM_SYMINV(XM,NOBS(IJ),NOBS(IJ),COND)

        DO IOBS=1,NOBS(IJ)-1
          DO JOBS=IOBS+1,NOBS(IJ)
           XM(IOBS,JOBS) = XM(JOBS,IOBS)
          ENDDO
        ENDDO

!       PRODUCE FIELD AT POINT I,J ACCORDING TO OPTIMUM INTERPOLATION

        DO IOBS=1,NOBS(IJ)
          DIFF = 0._JWRB
          DO JOBS=1,NOBS(IJ)
            KOBS=IOBS4IJ(IJ,JOBS)
            DIFF=DIFF+XM(IOBS,JOBS)*DIFFALTFG(KOBS)
          ENDDO
          XNEW=XNEW+W(IJ,IOBS)*DIFF
        ENDDO

        DEALLOCATE(XM)

!       UPDATE FIRST GUESS
        XOI(IJ) = XMO(IJ) + XNEW

      ELSE
!     NEGATIVE VALUES WILL REMAIN WHERE O.I. WILL PRODUCE NO RESULTS  
        XOI(IJ) =  -1._JWRB
      ENDIF
      ENDDO

IF (LHOOK) CALL DR_HOOK('WAMOI',1,ZHOOK_HANDLE)

END SUBROUTINE WAMOI
