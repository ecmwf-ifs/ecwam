      SUBROUTINE HALFSPREAD( 
!                  DIMENSIONS
     &               NSPEC , MPART , NPART , NANG , NFRE ,
!                  INPUT/OUTPUT
     &               SPEC,PART,PARTINFO,PEAKANG,PEAKFRE ,FX,FY,
     &               MEANE , SPREAD ,
!                  WORK SPACE
     &               SUMX , SUMY , SUMXX , SUMYY ,
!                  FREQU. TABLE
     &               FR , CONTROL, IU06 )

!-----------------------------------------------------------------------

!     PURPOSE:
!     --------

!     COMBINES PEAKS IF DISTANCE BETWEEN THEM IS LESS THAN HALF
!     THEIR SPREAD.

!     EXERNALS:
!     ---------

!     DOCOMBINE - TO DO THE COMBINING OF TWO PARTITIONS TO ONE

!-----------------------------------------------------------------------

      IMPLICIT NONE

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
!     PART(..) GIVES NUMBER OF PARTITIONING FOR EACH SPECTRAL BIN.
!     partinfo =1 -> wind sea, =2 -> mixed wind sea-swell, =0 -> swell
!      , 3=old windsea.
!      INDEX OF PEAK DIRECTION AND PEAK FREQUENCY OF PARTITIONINGS.
      REAL FX(NSPEC,MPART), FY(NSPEC,MPART),
     &     MEANE(NSPEC,MPART),
     &     SPREAD(NSPEC,MPART),
!                  mean wave numbers, mean energy and
!                  spread of partitionings.
     &     SUMX(NSPEC,MPART), SUMY(NSPEC,MPART),
     &     SUMXX(NSPEC,MPART), SUMYY(NSPEC,MPART)
!                  WORK SPACE.
      REAL FR(NFRE)
!                  TABLE OF FREQUENCIES.
      LOGICAL CONTROL
      INTEGER IU06

!     local variables:
!     ----------------

      INTEGER ISPEC, IPART, JPART
!                  LOOP INDEXES.

!-----------------------------------------------------------------------

      DO ISPEC = 1,NSPEC
        IPART = 1
        DO WHILE(IPART.LE.NPART(ISPEC)-1.AND.
     &                                  PARTINFO(ISPEC,IPART).EQ.0)
          JPART  = IPART + 1
          DO WHILE( IPART.NE.0 .AND. JPART.LE.NPART(ISPEC).AND.
     &                                  PARTINFO(ISPEC,JPART).EQ.0)
            IF( (FX(ISPEC,IPART)-FX(ISPEC,JPART))**2
     &          +(FY(ISPEC,IPART)-FY(ISPEC,JPART))**2
     &          .LT. 0.5*MIN(SPREAD(ISPEC,IPART),SPREAD(ISPEC,JPART))
     &          )THEN
              CALL DOCOMBINE( 
!                     DIMENSIONS
     &                  IPART , JPART , ISPEC ,
     &                  NSPEC , MPART , NPART , NANG , NFRE ,
!                     INPUT/OUTPUT
     &                  SPEC , PART , PARTINFO ,
     &                  PEAKANG , PEAKFRE , FX, FY,
     &                  MEANE , SPREAD ,
     &                  SUMX , SUMY , SUMXX , SUMYY ,
     &                  FR , CONTROL, IU06 )
              IF( CONTROL ) WRITE(IU06,*) ' REASON: ',
     &            'DISTANCE IS LESS THAN HALF OF MINIMAL SPREAD'
              IPART = 0
            ELSE
              JPART = JPART + 1
            END IF
          END DO
          IPART = IPART + 1
        END DO
      END DO

      RETURN

      END SUBROUTINE HALFSPREAD
