      SUBROUTINE WINDSEA( 
!                  DIMENSIONS
     &               NSPEC , MPART , NPART , NANG , NFRE ,
!                  INPUT
     &               SPEC , PART , 
!                  OUTPUT
     &               PARTINFO , 
!                  INPUT
     &               PEAKANG , PEAKFRE , FX, FY,
     &               MEANE ,  MEANANG, MEANFRE ,
     &               SPREAD , 
!                  WORK SPACE
     &               SUMX , SUMY , SUMXX , SUMYY ,
!                  INPUT TABLES
     &               U10 , THW ,TH0 , CONTROL, IU06 )

!-----------------------------------------------------------------------

!     PURPOSE:
!     --------

!     DECIDES WHICH ARE WINDSEA (=1), MIXED (=2), AND 
!     SWELL (=0) , OLD WINDSEA (=3) SYSTEMS,
!     STORES THIS INFORMATION IN ARRAY PARTINFO AND
!     COMBINES ALL WINDSEA PEAKS.
!     (AS A SIDE EFFECT THE ARRAYS FX AND FY ARE INITIALIZED)

!     EXERNALS:
!     ---------

!     DOCOMBINE - TO DO THE COMBINING OF TWO PARTITIONS TO ONE

!-----------------------------------------------------------------------

!     MODULE :
!     --------

      USE YOWFRED  , ONLY : FR       ,DELTH
      USE YOWPCONS , ONLY : G        ,ZPI

!-----------------------------------------------------------------------

      IMPLICIT NONE

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
!           PART(..) GIVES NUMBER OF PARTITIONING FOR EACH SPECTRAL BIN.
!           PARTINFO =1 -> WIND SEA, =2 -> MIXED WIND SEA-SWELL, 
!           =0 -> swell , 3=old windsea.
!           INDEX OF PEAK DIRECTION AND PEAK FREQUENCY OF PARTITIONINGS.

      REAL FX(NSPEC,MPART), FY(NSPEC,MPART),
     &     MEANE(NSPEC,MPART),MEANANG(NSPEC,MPART),
     &     MEANFRE(NSPEC,MPART),
     &     SPREAD(NSPEC,MPART),
!            MEAN WAVE NUMBERS, MEAN ENERGY AND SPREAD OF PARTITIONINGS.
     &     SUMX(NSPEC,MPART), SUMY(NSPEC,MPART),
     &     SUMXX(NSPEC,MPART), SUMYY(NSPEC,MPART)
!                  WORK SPACE.
      REAL U10(NSPEC) , THW(NSPEC)
!                  U10, WIND DIRECTION.
      REAL TH0
!                      FIRST DIRECTION IN SPECTRUM.
      LOGICAL CONTROL
      INTEGER IU06

!     LOCAL VARIABLES:
!     ----------------

      INTEGER  ISPEC, IPART, JPART
!                   LOOP INDEXES.
      REAL TWEG, THWQ, CM, CT, THQP, CTE, CMP, CTP,
     &     SQF, THQPP, THQPM, CTEP, CTEM, SQSPREAD, CTETWO

      LOGICAL LWINDSEA

!---------------------------------------------------------------------

!     1. SET ARRAY PARTINFO TO ZERO.
!     ------------------------------

      DO ISPEC = 1,NSPEC
        DO IPART = 1,NPART(ISPEC)
          PARTINFO(ISPEC,IPART) = 0
        END DO
      END DO

!     2. FIND ALL WINDSEA AND MIXED PEAKS.
!     ------------------------------------

      TWEG = 1.3 * ZPI / G

      DO ISPEC = 1,NSPEC
        THWQ = THW(ISPEC)
        LWINDSEA=.FALSE.
        DO IPART = 1,NPART(ISPEC)
          THQP = TH0 + DELTH * (PEAKANG(ISPEC,IPART)-1)
          CM = TWEG * FR(PEAKFRE(ISPEC,IPART))
          CT = CM * U10(ISPEC)
          CTE = CT * COS(THQP-THWQ) - 1.
            CTETWO = (2. / 1.3) * CT * COS(THQP-THWQ) - 1.
          IF( CTE.GT.0 )THEN
            LWINDSEA=.TRUE.
            PARTINFO(ISPEC,IPART) = 1
            ELSE IF ( CTETWO.GT.0 ) THEN
               LWINDSEA=.TRUE.
               PARTINFO(ISPEC,IPART) = 3
          ENDIF
        END DO
        IF(.NOT.LWINDSEA) THEN
          DO IPART = 1,NPART(ISPEC)
            SQSPREAD = SQRT(SPREAD(ISPEC,IPART)/2.)
            THQP = MEANANG(ISPEC,IPART)
            CM = TWEG * MEANFRE(ISPEC,IPART)
            CMP = CM + TWEG * SQSPREAD
            CTP = CMP * U10(ISPEC)
            SQF = SQSPREAD/(MEANFRE(ISPEC,IPART))
            THQPP = THQP + SQF
            THQPM = THQP - SQF
            CTEP = CTP * COS(THQPP-THWQ) - 1.
            CTEM = CTP * COS(THQPM-THWQ) - 1.
            IF( CTEP.GT.0. .OR. CTEM.GT.0. )THEN
              PARTINFO(ISPEC,IPART) = 2
            END IF
          END DO
        END IF
      END DO

!     3.COMBINE ALL WINDSEA PEAKS.
!     ----------------------------

      DO ISPEC = 1,NSPEC
        IPART = 0
        DO JPART = 1,NPART(ISPEC)
          IF( PARTINFO(ISPEC,JPART).GT.0 )THEN
            IF( IPART.EQ.0 )THEN
              IPART = JPART
            ELSE
              DO WHILE( PARTINFO(ISPEC,JPART).gt.0
     &            .AND.JPART.LE.NPART(ISPEC) )
                CALL DOCOMBINE( 
!                      DIMENSIONS
     &                   IPART , JPART , ISPEC ,
     &                   NSPEC , MPART , NPART , NANG , NFRE ,
!                      INPUT/OUTPUT
     &                   SPEC , PART , PARTINFO ,
     &                   PEAKANG , PEAKFRE , FX, FY,
     &                   MEANE , SPREAD , 
!                      WORK SPACE
     &                   SUMX , SUMY , SUMXX , SUMYY ,
!                      TABLES.
     &                   FR , CONTROL, IU06 )
                IF( CONTROL )
     &              WRITE(IU06,*) ' REASON: BOTH ARE WINDSEA PEAKS'
              END DO
            END IF
          END IF
        END DO
      END DO

      RETURN

      END SUBROUTINE WINDSEA
