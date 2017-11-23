      SUBROUTINE GAPINTERPOL( 
!              DIMENSIONS
     &             NSPEC, ISPECS ,ISPECE, MGAP ,NGAP, NANG, NFRE,
!              INPUT ( SPEC ALSO OUTPUT)
     &             SPEC ,
     &             FX , FY , LANG, HANG , LFRE , HFRE ,
     &             IU06 )

!-------------------------------------------------------------------

!     PURPOSE:
!     --------

!     INTERPOLATE GAPS.

!      AUTHOR.
!      -------
!      S.HASSELMANN, MPI HAMBURG, 1993.
!      J. WASZKEWITZ,MPI HAMBURG, 1993.

!     method.
!     -------
!     compute the 6 coefficients of a 2d parabola fit by minimizing
!     a cost function over all solutions gained at spectral points
!     with values > 0.
!     derivatives of the cost function with respect to the parabola
!     coefficients deliver a set of 6 equations with 6 unknowns
!     which are solved with nag routine f04jgf.

!-------------------------------------------------------------------
!     MODULE :
!     --------

      USE YOWGAP   , ONLY : GAPVAL

!-------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"

!     INTERFACE:
!     ----------


      INTEGER :: NGAP(ISPECS:ISPECE)
      INTEGER NSPEC, ISPECS, ISPECE, MGAP, NANG, NFRE
!                 DIMENSIONS.
      REAL SPEC(0:NSPEC,NANG,NFRE)
!                 2D SPECTRUM.
      REAL FX(NANG,NFRE), 
     &     FY(NANG,NFRE)
!                 POLAR COORDINATES.
      INTEGER LANG(ISPECS:ISPECE,MGAP), 
     &        HANG(ISPECS:ISPECE,MGAP),
     &        LFRE(ISPECS:ISPECE,MGAP),
     &        HFRE(ISPECS:ISPECE,MGAP)
!                 LOWER LEFT AND UPPER RIGHT CORNER OF GAP.
      INTEGER IU06
!     STANDARD OUTPUT UNIT

!     LOCAL VARIABLES:
!     ----------------

      INTEGER ISPEC, IGAP, IANG, IFRE, I12, I, J
!                 LOOP INDICES.
      INTEGER LANG0(2),HANG0(2),LFRE0,HFRE0
!                 CORNERS OF FRAMES FOR TEMPORARILY USE
      REAL X, Y
!                 FREQUENCY VALUES IN CARTHESIAN COORDINATES
      REAL Z_REAL
      DOUBLE PRECISION DL_DOUBLE
!                 USED TO CHOOSE BETWEEN SGELSS OR DGELSS (LAPACK)
      REAL A(6,6), B(6), T(6), S(6), WORK(1000), RCOND
      INTEGER IRANK, INFO
!                 INFO CONTAINS AN ERROR FLAG AFTER CALLING THE
!                 LAPACK ROUTINE SGELSS

!-------------------------------------------------------------------

!     1.1 CHECK POSITION OF CORNERS WITH RESPECT TO 0 DEGREE.
!     ------------------------------------------------------

!     LOOP OVER SPECTRA.
      DO ISPEC = ISPECS,ISPECE 
        DO IGAP = 1,NGAP(ISPEC)
!     LOOP OVER GAPS.
          IF( LANG(ISPEC,IGAP).EQ.1.AND.HANG(ISPEC,IGAP).EQ.NANG.OR.
     &        LANG(ISPEC,IGAP)-HANG(ISPEC,IGAP).EQ.2.OR.
     &        LANG(ISPEC,IGAP)-HANG(ISPEC,IGAP)+NANG.EQ.2 )THEN
            LANG0(1) = 1
            HANG0(1) = NANG
            LANG0(2) = 1
            HANG0(2) = 0
          ELSE
            IF( LANG(ISPEC,IGAP).GT.HANG(ISPEC,IGAP) )THEN
              LANG0(1) = 1
              HANG0(1) = HANG(ISPEC,IGAP) + 1
              LANG0(2) = LANG(ISPEC,IGAP) - 1
              HANG0(2) = NANG
            ELSE
              LANG0(1) = MAX( LANG(ISPEC,IGAP)-1 , 1 )
              HANG0(1) = MIN( HANG(ISPEC,IGAP)+1 , NANG )
              IF( LANG(ISPEC,IGAP).EQ.1 )THEN
                LANG0(2) = NANG
                HANG0(2) = NANG
              ELSE IF( HANG(ISPEC,IGAP).EQ.NANG )THEN
                LANG0(2) = 1
                HANG0(2) = 1
              ELSE
                LANG0(2) = 1
                HANG0(2) = 0
              END IF
            END IF
          END IF
          LFRE0 = MAX(1,LFRE(ISPEC,IGAP)-1)
          HFRE0 = MIN(NFRE,HFRE(ISPEC,IGAP)+1)

!     1.2 COMPUTE COEFFICIENT MATRIX FOR LINEAR SYSTEM OF 6 
!         EQUATIONS.
!     ---------------------------------------------------------

          DO I = 1,6
            DO J = 1,6
              A(I,J) = 0.
            END DO
            B(I) = 0.
          END DO
          DO I12 = 1,2
            DO IANG = LANG0(I12),HANG0(I12)
              DO IFRE = LFRE0,HFRE0
                IF( SPEC(ISPEC,IANG,IFRE).NE.0. )THEN
                  X = FX( IANG , IFRE )
                  Y = FY( IANG , IFRE )
                  T(1) = 1.
                  T(2) = X
                  T(3) = Y
                  T(4) = X*X
                  T(5) = Y*Y
                  T(6) = X*Y
                  DO I = 1,6
                    DO J= 1,6
                      A(I,J) = A(I,J) + T(I) * T(J)
                    END DO
                    B(I) = B(I) + T(I) * SPEC(ISPEC,IANG,IFRE)
                  END DO
                END IF
              END DO
            END DO
          END DO

!      1.3 SOLVE SYSTEM OF EQUATIONS WITH RESPECT TO COEFFICIENTS 
!          FOR 2D PARABOLA FIT.
!          (TRY TO BE SMART AND DETECT COMPILER-FLAG PROMOTION OF REALS.
!           THIS ASSUMES THAT DOUBLE PRECISION IS NOT ALSO PROMOTED.)
!      ------------------------------------------------------------


          INFO  = 0
          RCOND = 5.E-4
          IF (KIND(Z_REAL) == KIND(DL_DOUBLE)) THEN
            CALL DGELSS (6,6,1,A,6,B,6,S,RCOND,IRANK,WORK,SIZE(WORK),
     &                   INFO)
          ELSE
            CALL SGELSS (6,6,1,A,6,B,6,S,RCOND,IRANK,WORK,SIZE(WORK),
     &                   INFO)
          ENDIF

          IF( INFO < 0 )THEN
            WRITE(IU06,*) 'ERROR! SUBROUTINE GAPINTERPOL:'
            WRITE(IU06,*) 'SGELSS RETURNED EXIT CODE INFO=',INFO
            WRITE(IU06,*) 'INFO<0 MEANS INFO''TH ARGUMENT IS INVALID'
            WRITE(IU06,*) 'SEE DOCUMENTATION FOR LAPACK LIBRARY'
            WRITE(IU06,*) 'FOR MORE INFORMATION.'
            CALL ABORT1
          END IF

          IF( INFO.EQ.0 )THEN
            B(1) = MAX(B(1),0.)
            DO I12 = 1,2
              DO IANG = LANG0(I12),HANG0(I12)
                DO IFRE = LFRE0,HFRE0
                  IF( SPEC(ISPEC,IANG,IFRE).EQ.GAPVAL )THEN
                    X = FX( IANG , IFRE )
                    Y = FY( IANG , IFRE )
                    SPEC(ISPEC,IANG,IFRE) = MAX( 0. ,
     &                B(1)+B(2)*X+B(3)*Y+B(4)*X*X+B(5)*Y*Y+B(6)*X*Y)
                  END IF
                END DO
              END DO
            END DO
          END IF

!       END OF LOOPS OVER NUMBER OF GAPS.

        END DO
!    
!     END OF LOOP OVER NUMBER OF SPECTRA.

      END DO

      RETURN

      END SUBROUTINE GAPINTERPOL
