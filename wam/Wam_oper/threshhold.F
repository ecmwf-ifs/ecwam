      SUBROUTINE THRESHHOLD( 
!                  DIMENSIONS
     &               NSPEC , MPART , NPART , NANG , NFRE ,
!                  INPUT/OUTPUT
     &               SPEC,PART,PARTINFO,PEAKANG,PEAKFRE,FX,FY,
     &               MEANE,SPREAD,
!                  WORK SPACE
     &               SUMX,SUMY,SUMXX,SUMYY,
!                  INPUT TABLES
     &               FR , LEVEL , CONTROL, IU06 )

!-----------------------------------------------------------------------

!     PURPOSE:
!     --------

!     COMBINES PARTS IF MINIMUM ENERGY BETWEEN PEAKS IS LARGER THAN
!     >LEVEL<  TIMES THE LOWEST ENERGY OF THE PEAKS. 

!     METHOD:
!     -------

!     FOR EVERY SPECTRUM AND EVERY PAIR OF PEAKS "WALK" THE MOST 
!     DIRECT WAY IN THE GRID FROM ONE PEAK TO THE OTHER. 

!     EXERNALS:
!     ---------

!     DOCOMBINE - COMBINING OF TWO PARTITIONS.

!-----------------------------------------------------------------------

      IMPLICIT NONE

!     INTERFACE:
!     ----------

      INTEGER NSPEC, MPART, NPART(NSPEC), NANG, NFRE
!                  ARRAY DIMENSIONS.
      REAL SPEC(0:NSPEC,NANG,NFRE)
!                  2D SPECTRUM.
      INTEGER PART(NSPEC,NANG,NFRE),
     &        PARTINFO(NSPEC,MPART),
     &        PEAKANG(NSPEC,NANG*NFRE/2),
     &        PEAKFRE(NSPEC,NANG*NFRE/2)
!      PART(..) GIVES NUMBER OF PARTITIONING FOR EACH SPECTRAL BIN.
!      partinfo =1 -> wind sea, =2 -> mixed wind sea-swell =0 -> swell
!       , 3=old windsea.
!      INDEX OF PEAK DIRECTION AND PEAK FREQUENCY OF PARTITIONINGS.

      REAL FX(NSPEC,MPART), FY(NSPEC,MPART),
     &     MEANE(NSPEC,MPART),
     &     SPREAD(NSPEC,MPART),
!                  MEAN WAVE NUMBERS, MEAN ENERGY AND SPREAD OF 
!                  PARTITIONINGS.
     &     SUMX(NSPEC,MPART), SUMY(NSPEC,MPART),
     &     SUMXX(NSPEC,MPART), SUMYY(NSPEC,MPART)
!                  WORK SPACE.
      REAL FR(NFRE), LEVEL
!                  TABLE OF FREQUENCIES, THRESHHOLD LEVEL.
      LOGICAL CONTROL
      INTEGER IU06

!     LOCAL VARIABLES:
!     ----------------

      INTEGER ISPEC, IPART, JPART, IANG, IFRE
!                   LOOP INDEXES. 
      LOGICAL COMBINE
      INTEGER DISTANG, DISTFRE
!                   ABSOLUT VALUE OF DISTANCE IN GRIDPOINTS BETWEEN
!                   BOTH PEAKS IN ANGLE OR FREQUENCY RSP.
      INTEGER DIST
!                   THE DISTANCE IN GRIDPOINTS BETWEEN TWO PEAKS. 
      INTEGER STEPANG, STEPFRE
!                   DIRECTION FOR "WALKING" IN THE GRID FROM ONE PEAK TO
!                   THE OTHER,CONTAINS -1 OR 1.
      INTEGER IDIST
!                   COUNTS STEPS WHILE "WALKING" IN THE GRID FROM ONE 
!                   PEAK TO THE OTHER.
      INTEGER SUMANG, SUMFRE
!                  SUMS TO DECIDE IF TO MAKE A ANG-STEP OR A FRE-STEP
!                  while "walking" in the grid from one peak to
!                  the other
      INTEGER DIST2ANG, DIST2FRE
!                  TWICE THE VALUE DISTANG OR DISTFRE RSP. 
      REAL PEAKENERGY
      INTEGER TEMP1, TEMP2, IPA, JPA
!                  FOR TEMPORARILY USE

!---------------------------------------------------------------------

      DO ISPEC = 1,NSPEC
        IPART = 1
        DO WHILE( IPART.LE.NPART(ISPEC)-1 .AND.
     &                                PARTINFO(ISPEC,IPART).EQ.0)
          JPART = IPART + 1
          DO WHILE( IPART.NE.0 .AND. JPART.LE.NPART(ISPEC) .AND.
     &                                PARTINFO(ISPEC,JPART).EQ.0)
!           "WALKING" IN THE GRID FROM PEAK I TO PEAK J
            DISTANG = PEAKANG(ISPEC,JPART) - PEAKANG(ISPEC,IPART)
            IF( DISTANG.LT.0 ) DISTANG = DISTANG + NANG
            IF( DISTANG.LT.NANG/2 )THEN
!             "WALK" FROM "LEFT" TO "RIGHT"
              STEPANG = 1
              COMBINE = .TRUE.
            ELSE IF( DISTANG.GT.NANG/2 )THEN
!             "WALK" FROM "RIGHT" TO "LEFT"
              DISTANG = NANG - DISTANG
              STEPANG = -1
              COMBINE = .TRUE.
            ELSE
!             DO NOT COMBINE
              STEPANG = 0
              COMBINE = .FALSE.
            ENDIF
            DISTFRE = PEAKFRE(ISPEC,JPART)-PEAKFRE(ISPEC,IPART)
            STEPFRE = SIGN(1,DISTFRE)
            DISTFRE = ABS(DISTFRE)
            DIST = DISTANG  + DISTFRE
            PEAKENERGY = MIN( 
     &          SPEC(ISPEC,PEAKANG(ISPEC,IPART),PEAKFRE(ISPEC,IPART)),
     &          SPEC(ISPEC,PEAKANG(ISPEC,JPART),PEAKFRE(ISPEC,JPART)))
            IANG = PEAKANG(ISPEC,IPART)
            IFRE = PEAKFRE(ISPEC,IPART)

            IPA=PART(ISPEC,IANG,IFRE)
            JPA=PART(ISPEC,PEAKANG(ISPEC,JPART),PEAKFRE(ISPEC,JPART))

            SUMANG = ABS(DISTANG)
            SUMFRE = ABS(DISTFRE)
            DIST2ANG = 2 * DISTANG
            DIST2FRE = 2 * DISTFRE
            IDIST = 1
            DO WHILE( COMBINE .AND. IDIST.LT.DIST )
              IF( SUMANG.GT.SUMFRE )THEN
                IANG = IANG + STEPANG
                IF( IANG.LE.0 )THEN
                  IANG = IANG + NANG
                ELSE IF( IANG.GT.NANG )THEN
                  IANG = IANG - NANG
                END IF
                SUMFRE = SUMFRE + DIST2FRE
                IDIST = IDIST + 1
                IF( SPEC(ISPEC,IANG,IFRE).LT.LEVEL*PEAKENERGY )
     &              COMBINE = .FALSE.
              ELSE IF( SUMFRE.GT.SUMANG )THEN
                IFRE = IFRE + STEPFRE
                SUMANG = SUMANG + DIST2ANG
                IDIST = IDIST + 1
                IF( SPEC(ISPEC,IANG,IFRE).LT.LEVEL*PEAKENERGY )
     &              COMBINE = .FALSE.
              ELSE
                TEMP1 = IANG
                IANG = IANG + STEPANG
                IF( IANG.LE.0 )THEN
                  IANG = IANG + NANG
                ELSE IF( IANG.GT.NANG )THEN
                  IANG = IANG - NANG
                END IF
                IF( SPEC(ISPEC,IANG,IFRE).LT.LEVEL*PEAKENERGY )
     &              COMBINE = .FALSE.
                TEMP2 = IANG
                IANG = TEMP1
                IFRE = IFRE + STEPFRE
                IF( SPEC(ISPEC,IANG,IFRE).LT.LEVEL*PEAKENERGY )
     &              COMBINE = .FALSE.
                IANG = TEMP2
                SUMANG = SUMANG + DIST2ANG
                SUMFRE = SUMFRE + DIST2FRE
                IDIST = IDIST + 2
              END IF
                IF(PART(ISPEC,IANG,IFRE).NE.IPA.AND.
     &             PART(ISPEC,IANG,IFRE).NE.JPA) 
     &              COMBINE = .FALSE.

            END DO
            IF( COMBINE )THEN
              CALL DOCOMBINE( IPART , JPART , ISPEC ,
     &                   NSPEC , MPART , NPART , NANG , NFRE ,
     &            SPEC , PART , PARTINFO ,
     &            PEAKANG , PEAKFRE , FX, FY,
     &            MEANE , SPREAD ,
     &            SUMX , SUMY , SUMXX , SUMYY ,
     &            FR , CONTROL, IU06 )
              IF( CONTROL ) WRITE(IU06,*) ' REASON: ',
     &            'MINIMUM ENERGY BETWEEN THEM IS HIGHER THAN ',
     &            LEVEL,' OF THE MINIMUM PEAK ENERGY'
              IPART = 0
            ELSE
              JPART = JPART + 1
            END IF
          END DO
          IPART = IPART + 1
        END DO
      END DO

      RETURN

      END SUBROUTINE THRESHHOLD
