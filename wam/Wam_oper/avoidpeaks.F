      SUBROUTINE AVOIDPEAKS( 
!                   DIMENSIONS
     &                 NSPEC, MGAP ,NGAP, MPART, NANG, NFRE,
!                   INPUT 
     &                 SPEC ,
     &                 NPEAK , PEAKANG , PEAKFRE ,
!                   ARRAYS TO BE ADJUSTED
     &                 LANG, HANG , LFRE , HFRE )

!-------------------------------------------------------------------

!     PURPOSE:
!     --------

!     THE FRAMES AROUND THE GAPS ARE REDUCED TO AVOID PEAKS 
!     OTHERWISE THE PARABOLA FIT CANNOT BE COMPUTED.

!     METHOD:
!     -------

!     IF A PEAK LIES INSIDE A FRAME, THE FRAME IS SPLIT INTO TWO
!     FRAMES TO AVOID PEAKS.

!     AUTHOR.
!     -------

!     S.HASSELMANN MPI HAMNBURG. 1993

!-------------------------------------------------------------------

      IMPLICIT NONE

!     INTERFACE:
!     ----------
      INTEGER NSPEC, MGAP, NGAP(NSPEC), MPART, NANG, NFRE
!                 DIMENSIONS.
      REAL SPEC(0:NSPEC,NANG,NFRE)
!                 2D SPECTRUM.
      INTEGER LANG(NSPEC,MGAP), 
     &        HANG(NSPEC,MGAP),
     &        LFRE(NSPEC,MGAP),
     &        HFRE(NSPEC,MGAP)
!                 THE LEFT LOWER AND RIGHT UPPER CORNER OF THE FRAME
!                 IN ANG AND FRE, WILL BE ADJUSTED.
      INTEGER NPEAK(NSPEC)
!                 NUMBER OF PEAKS
      INTEGER PEAKANG(NSPEC,MPART),
     &        PEAKFRE(NSPEC,MPART)
!                 THE DIRECTIONAL AND FREQUENCY INDICES OF THE PEAK
 
!     LOCAL  VARIABLES:
!     -----------------

      INTEGER ISPEC, IGAP, IANG, IFRE
!                  LOOP INDICES.
      INTEGER IANGP1, IANGM1, IFREP1, IFREM1
!                  LOOP INDICES + OR - 1
      REAL SPEC0
!                  TEMPORARY STORAGE
!-------------------------------------------------------------------

!     IF PEAK FOUND IN FRAME BOUNDARY MOVED TO PEAK.
!     ----------------------------------------------

      DO ISPEC=1,NSPEC
        DO IGAP=1,NGAP(ISPEC)
          DO IANG=LANG(ISPEC,IGAP),HANG(ISPEC,IGAP)
            DO IFRE=LFRE(ISPEC,IGAP),HFRE(ISPEC,IGAP)
              IANGP1 = IANG + 1
              IF( IANGP1.GT.NANG ) IANGP1 = 1
              IANGM1 = IANG - 1
              IF( IANGM1.EQ.0 ) IANGM1 = NANG
              IFREP1 = IFRE + 1
              IF( IFREP1.GT.NFRE ) IFREP1 = NFRE
              IFREM1 = IFRE - 1
              IF( IFREM1.EQ.0 ) IFREM1 = 1
              SPEC0 = SPEC(ISPEC,IANG,IFRE)
              IF( SPEC0.GT.0..AND.
     &            SPEC0.GE.SPEC(ISPEC,IANG  ,IFREM1).AND.
     &            SPEC0.GE.SPEC(ISPEC,IANG  ,IFREP1).AND.
     &            SPEC0.GE.SPEC(ISPEC,IANGP1,IFRE  ).AND.
     &            SPEC0.GE.SPEC(ISPEC,IANGM1,IFRE  ) )THEN
                IF((IANG-LANG(ISPEC,IGAP)).LT.
     &            (HANG(ISPEC,IGAP)-IANG)) THEN
                  LANG(ISPEC,IGAP) =IANG
                ELSE
                  HANG(ISPEC,IGAP) =IANG
                END IF
              END IF
            END DO
          END DO
        END DO
      END DO

      RETURN

      END SUBROUTINE AVOIDPEAKS
