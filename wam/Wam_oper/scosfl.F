      SUBROUTINE SCOSFL (F, IJS, IJL, MM, MEANCOSFL)

! ----------------------------------------------------------------------

!**** *SCOSFL* - COMPUTATION OF THE CENTERED FOURIER COEFFICIENT FOR
!                A GIVEN FREQUENCY AT EACH GRID POINT.

!     J.-R. BIDLOT        ECMWF       MARCH 2000

!*    PURPOSE.
!     --------

!     TO COMPUTE THE CENTERED FOURIER COEFFICIENT FOR A GIVEN FREQUENCY 

!**   INTERFACE.
!     ----------

!       *CALL* *SCOSFL (F, IJS, IJL, MM, MEANCOSFL)*
!          *F*  - SPECTRUM.
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *MM*  - FREQUENCY INDEX (array)
!          *MEANCOSFL* - CENTERED FOURIER COEFFICIENT

!     METHOD.
!     -------

!       INTEGRATION OF SPECTRUM TIMES COS(THETA-MEAN THETA)
!       OVER DIRECTION.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : ZPI      ,EPSMIN
      USE YOWFRED  , ONLY : DELTH    ,TH       ,COSTH    ,SINTH
      USE YOWMPP   , ONLY : NINF     ,NSUP

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IJ,K,IJS,IJL
      INTEGER,DIMENSION(IJS:IJL) :: MM
      REAL,DIMENSION(IJS:IJL) :: MEANDIR, MEANCOSFL
      REAL :: SI(IJS:IJL), CI(IJS:IJL), F(NINF-1:NSUP,NANG,NFRE)

! ----------------------------------------------------------------------

!*    1. INITIALISE ARRAYS
!        -----------------

      DO IJ=IJS,IJL
        SI(IJ) = 0.
        CI(IJ) = 0.
        MEANCOSFL(IJ) = 0.
      ENDDO

! ----------------------------------------------------------------------

!*    2. FIND MEAN DIRECTION
!        -------------------

      DO K=1,NANG
        DO IJ=IJS,IJL
          SI(IJ) = SI(IJ)+SINTH(K)*F(IJ,K,MM(IJ))
          CI(IJ) = CI(IJ)+COSTH(K)*F(IJ,K,MM(IJ))
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        IF (CI(IJ).EQ.0. .AND. SI(IJ).EQ.0.) THEN
          MEANDIR(IJ) = 0.
        ELSE
          MEANDIR(IJ) = ATAN2(SI(IJ),CI(IJ))
        ENDIF
      ENDDO

!*    3. COMPUTE FOURRIER COEFFICIENT 
!        -----------------------------

      DO K=1,NANG
        DO IJ=IJS,IJL
          MEANCOSFL(IJ)=MEANCOSFL(IJ)+
     &                  COS(TH(K)-MEANDIR(IJ))*F(IJ,K,MM(IJ))
        ENDDO
      ENDDO
      DO IJ=IJS,IJL
        MEANCOSFL(IJ) = DELTH*MEANCOSFL(IJ)
      ENDDO

      RETURN
      END SUBROUTINE SCOSFL
