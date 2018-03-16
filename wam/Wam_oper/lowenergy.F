      SUBROUTINE LOWENERGY( 
!                 DIMENSIONS
     &              NSPEC , MPART , NPART , NANG , NFRE ,
!                 INPUT/OUTPUT
     &              SPEC , PART , PARTINFO , PEAKANG , PEAKFRE , FX, FY,
     &              LOWEST, MEANE, SPREAD , 
!                 WORK SPACE.
     &              SUMX , SUMY , SUMXX , SUMYY ,
!                 INPUT TABLE
     &              FR , CONTROL, IU06 )

!---------------------------------------------------------------------

!     PURPOSE:
!     --------

!     IF MEAN ENERGY OF PARTITIONING IS < LOWEST
!     COMBINE WITH NEAREST NOT WINDSEA PEAK

!     EXERNALS:
!     ---------

!     DOCOMBINE - COMBINES TWO PARTITIONINGS.

!---------------------------------------------------------------------

      IMPLICIT NONE
#include "docombine.intfb.h"

!     INTERFACE:
!     ----------

      INTEGER NSPEC, MPART, NPART(NSPEC), NANG , NFRE
!                  ARRAY DIMENSIONS.
      REAL SPEC(0:NSPEC,NANG,NFRE)
!                  2D SPECTRA.
      INTEGER PART(NSPEC,NANG,NFRE),
     &        PARTINFO(NSPEC,MPART),
     &        PEAKANG(NSPEC,NANG*NFRE/2),
     &        PEAKFRE(NSPEC,NANG*NFRE/2)
!       PART(..) GIVES NUMBER OF PARTITIONING FOR EACH SPECTRAL BIN.
!       partinfo =1 -> wind sea, =2 -> mixed wind sea-swell =0 -> swell
!       , 3=old windsea.
!       INDEX OF PEAK DIRECTION AND PEAK FREQUENCY OF PARTITIONINGS.
      REAL FX(NSPEC,MPART), FY(NSPEC,MPART),
     &     MEANE(NSPEC,MPART),
     &     SPREAD(NSPEC,MPART),
!             MEAN WAVE NUMBERS, MEAN ENERGY AND SPREAD OF PARTITIONINGS.
     &     SUMX(NSPEC,MPART), SUMY(NSPEC,MPART),
     &     SUMXX(NSPEC,MPART), SUMYY(NSPEC,MPART)
!                  WORK SPACE.
      REAL FR(NFRE)
!                  TABLE OF FREQUENCIES.
      LOGICAL CONTROL
      INTEGER IU06

!     local variables:
!     ----------------

      INTEGER ISPEC, IPART, JPART, KPART
!                   LOOP INDEXES.
      REAL DIST, DIST0
!                   FOR SQUARES OF DISTANCES
      REAL LOWEST
!                   THE LOWEST WAVE HEIGHT WHICH IS NOT COMBINED

!-----------------------------------------------------------------------

      DO ISPEC = 1,NSPEC
        DO IPART = 1,NPART(ISPEC)
          KPART = -1
          DO WHILE( MEANE(ISPEC,IPART).LT.(LOWEST/4)**2
     &        .AND. PARTINFO(ISPEC,IPART).EQ.0
     &        .AND. IPART.LE.NPART(ISPEC) .AND. KPART.NE.0 )
!           COMBINE WITH NEAREST NOT WINDSEA PEAK
!           FIRST: FIND NEAREST NOT WINDSEA PEAK TO IPEAK
            KPART = 0
            DIST = 100.
            DO JPART = 1,NPART(ISPEC)
               IF( JPART.NE.IPART)THEN
                DIST0 = (FX(ISPEC,IPART)-FX(ISPEC,JPART))**2
     &              + (FY(ISPEC,IPART)-FY(ISPEC,JPART))**2
                IF( KPART.EQ.0 .OR. DIST0.LT.DIST )THEN
                  KPART = JPART
                  DIST = DIST0
                END IF
              END IF
            END DO
!           SECOND: COMBINE 
            IF( KPART.NE.0 )THEN
              CALL DOCOMBINE( 
!                    DIMENSIONS
     &                 IPART , KPART , ISPEC ,
     &                 NSPEC , MPART , NPART , NANG , NFRE ,
!                    INPUT/OUTPUT
     &                 SPEC , PART , PARTINFO ,
     &                 PEAKANG , PEAKFRE , FX, FY,
     &                 MEANE , SPREAD , 
!                    WORK SPACE 
     &                 SUMX , SUMY , SUMXX , SUMYY ,
!                    INPUT TABLE
     &                 FR , CONTROL, IU06 )
            END IF
          END DO
        END DO
      END DO

      RETURN

      END SUBROUTINE LOWENERGY
