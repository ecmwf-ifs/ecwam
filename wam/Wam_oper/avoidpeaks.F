      SUBROUTINE AVOIDPEAKS( 
!                   DIMENSIONS
     &                 NSPEC,ISPECS,ISPECE,MGAP,NGAP,MPART,NANG,NFRE,
!                   INPUT 
     &                 SPEC ,
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
      INTEGER :: NGAP(ISPECS:ISPECE)
      INTEGER NSPEC, ISPECS, ISPECE, MGAP, MPART, NANG, NFRE
!                 DIMENSIONS.
      REAL SPEC(0:NSPEC,NANG,NFRE)
!                 2D SPECTRUM.
      INTEGER LANG(ISPECS:ISPECE,MGAP), 
     &        HANG(ISPECS:ISPECE,MGAP),
     &        LFRE(ISPECS:ISPECE,MGAP),
     &        HFRE(ISPECS:ISPECE,MGAP)
!                 THE LEFT LOWER AND RIGHT UPPER CORNER OF THE FRAME
!                 IN ANG AND FRE, WILL BE ADJUSTED.
 
!     LOCAL  VARIABLES:
!     -----------------

      INTEGER ISPEC, IGAP, IANG, IFRE
!                  LOOP INDICES.
      INTEGER IANGP1, IANGM1, IFREP1, IFREM1
!                  LOOP INDICES + OR - 1
!     LOWER AND UPPER LIMIT ON FIRST DIMENSION POINTINF TO NUMBER OF 
!     SPECTRA.
      REAL SPEC0
!                  TEMPORARY STORAGE
!-------------------------------------------------------------------

!     IF PEAK FOUND IN FRAME BOUNDARY MOVED TO PEAK.
!     ----------------------------------------------

      DO ISPEC=ISPECS ,ISPECE
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
