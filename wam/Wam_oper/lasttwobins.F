      SUBROUTINE LASTTWOBINS(
!                  DIMENSIONS
     &               NSPEC , MPART , NPART , NANG , NFRE ,
!                  INPUT/OUTPUT
     &               SPEC , PART , PARTINFO , PEAKANG , PEAKFRE ,FX,FY,
     &               MEANE, SPREAD ,
!                  WORK SPACE
     &               SUMX , SUMY , SUMXX , SUMYY ,
!                  FREQ TABLE.
     &               FR , CONTROL, IU06 )

!-----------------------------------------------------------------------

      IMPLICIT NONE

!     PURPOSE:
!     --------

!     IF PEAK LIES IN LAST TWO FREQUENCY BINS AND IS NOT WINDSEA,
!     COMBINE WITH NEAREST NOT WINDSEA PEAK
!!    IF THE NUMBER OF FREQUENCY IF 30 THEN DO THE SAME FOR THE LAST 5.

!     EXERNALS:
!     ---------

!     DOCOMBINE - TO DO THE COMBINING OF TWO PARTITIONS TO ONE

!     INTERFACE:
!     ----------

      INTEGER NSPEC, MPART, NPART(NSPEC), NANG, NFRE
!                  ARRAY DIMENSIONS.
      REAL SPEC(0:NSPEC,NANG,NFRE)
!                  2D SPECTRA.
      INTEGER PART(NSPEC,NANG,NFRE),
     &        PARTINFO(NSPEC,MPART),
     &        PEAKANG(NSPEC,NANG*NFRE/2),
     &        PEAKFRE(NSPEC,NANG*NFRE/2)
!      PART(..) GIVES NUMBER OF PARTITIONING FOR EACH SPECTRAL BIN.
!      PARTINFO =1 -> WIND SEA, =2 -> MIXED WIND SEA-SWELL, =3 -> SWELL.
!      INDEX OF PEAK DIRECTION AND PEAK FREQUENCY OF PARTITIONINGS.
      REAL FX(NSPEC,MPART), FY(NSPEC,MPART),
     &     MEANE(NSPEC,MPART),
     &     SPREAD(NSPEC,MPART),
!            MEAN WAVE NUMBERS, MEAN ENERGY AND SPREAD OF PARTITIONINGS.
     &     SUMX(NSPEC,MPART), SUMY(NSPEC,MPART),
     &     SUMXX(NSPEC,MPART), SUMYY(NSPEC,MPART)
!                  WORK SPACE.
      REAL FR(NFRE)
!                  TABLE OF FREQUENCIES.
      LOGICAL CONTROL
      INTEGER IU06

!     VARIABLES:
!     ----------

      INTEGER ISPEC, IPART, JPART, KPART, KANG, KFRE, NBINEXC
!                   LOOP INDEXES.
      REAL F0, FX0, FY0
!                   WORK SPACE.
      REAL DIST, DIST0
!                   SQUARES OF DISTANCES

!-----------------------------------------------------------------------

      IF ( NFRE . LT. 30 ) THEN
        NBINEXC=2
      ELSE
        NBINEXC=5
      ENDIF

!     1. FIND NEAREST NEIGHBOUR.
!     ---------------------------
      DO ISPEC = 1,NSPEC
        DO IPART = 1,NPART(ISPEC)
          KPART = -1
          DO WHILE( PEAKFRE(ISPEC,IPART).GT.NFRE-NBINEXC
     &        .AND. PARTINFO(ISPEC,IPART).EQ.0
     &        .AND. IPART.LE.NPART(ISPEC) .AND. KPART.NE.0 )
!           COMBINE WITH NEAREST NOT WINDSEA PEAK
!           FIRST: FIND NEAREST NOT WINDSEA PEAK TO IPEAK
            KPART = 0
            DIST = 0
            DO JPART = 1,NPART(ISPEC)
              IF( JPART.NE.IPART .AND. PARTINFO(ISPEC,JPART).EQ.0 
     &            .AND. PEAKFRE(ISPEC,JPART).LE.NFRE-NBINEXC )THEN
                DIST0 = (FX(ISPEC,IPART)-FX(ISPEC,JPART))**2
     &              + (FY(ISPEC,IPART)-FY(ISPEC,JPART))**2
                IF( KPART.EQ.0 .OR. DIST0.LT.DIST )THEN
                  KPART = JPART
                  DIST = DIST0
                END IF
              END IF
            END DO

!      2. COMBINE WAVE SYSTEMS.
!      ------------------------

            IF( KPART.NE.0 )THEN
              KANG = PEAKANG(ISPEC,KPART)
              KFRE = PEAKFRE(ISPEC,KPART)
              F0 = FR(PEAKFRE(ISPEC,IPART))
              FX0 = FX(ISPEC,KPART)
              FY0 = FY(ISPEC,KPART)
              CALL DOCOMBINE( IPART , KPART , ISPEC ,
     &            NSPEC , MPART , NPART , NANG , NFRE ,
     &            SPEC , PART , PARTINFO ,
     &            PEAKANG , PEAKFRE , FX, FY,
     &            MEANE , SPREAD ,
     &            SUMX , SUMY , SUMXX , SUMYY ,
     &            FR , CONTROL, IU06 )
!             THE PEAK IN THE LAST TWO BINS IS UNBELIEVABLE,
!             SO THE COMBINED PART IS SAID TO HAVE 
!             THE PEAK WHERE THE OTHER PART HAS ITS PEAK
!             (THE PART-INDEX OF THE COMBINED PEAK IS MIN(IPART,KPART))
              PEAKANG(ISPEC,MIN(IPART,KPART)) = KANG
              PEAKFRE(ISPEC,MIN(IPART,KPART)) = KFRE
              FX(ISPEC,MIN(IPART,KPART)) = FX0
              FY(ISPEC,MIN(IPART,KPART)) = FY0
            END IF
          END DO
        END DO
      END DO

      RETURN

      END SUBROUTINE LASTTWOBINS
