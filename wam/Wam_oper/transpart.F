      SUBROUTINE TRANSPART( 
!                    DIMENSIONS
     &                 NSPEC , MPART , NPARTW ,NANG , NFRE ,
!                    INPUT/OUTPUT
     &                 SPECW , 
!                    WORK SPACE
     &                 SPECT , 
!                    INPUT
     &                 PARTW ,CORRELTAW ,MEANWE,MEANWANG,MEANWFRE,
     &                 MEANSE , MEANSANG , MEANSFRE ,
     &                 FRATIO,
!                    OUTPUT
     &                 ADJUST , INTROTATE , FRACROTATE , 
     &                 STRETCH , INTLOGSTRETCH ,
     &                 COPO , COPODIFF  )

!-------------------------------------------------------------------

!     PURPOSE:
!     --------

!     TRANSFORMS PARTITIONINGS OF WAM FIRST GUESS SPECTRA TO 
!     MATCH THE MEAN VALUES TO THE PARTITIONINGS OF THE SAR 
!     RETRIEVED SPECTRA

!     S. HASSELMANN, MPI HAMBURG, 1993.
!     J. WASZKEWITZ, MPI HAMBURG, 1993.

!-------------------------------------------------------------------

!     MODULE :
!     --------

      USE YOWFRED  , ONLY : DELTH
      USE YOWGAP   , ONLY : GAPVAL 

!-------------------------------------------------------------------

      IMPLICIT NONE

!     INTERFACE:
!     ----------

      INTEGER NSPEC, MPART, NPARTW(NSPEC), NANG, NFRE
!                   DIMENSIONS.
      REAL SPECW(0:NSPEC,NANG,NFRE), 
     &     SPECT(NSPEC,NANG,NFRE)
!                   WAM SPECTRM BEFORE AND AFTER CORRECTION 
!                   AND WORK SPACE.
      INTEGER PARTW(NSPEC,NANG,NFRE)
!                   ASSIGNMENT OF A PARTITIONING NUMBER TO EACH 
!                   DIRECTIONAL-FREQUENCY BIN.
      INTEGER CORRELTAW(NSPEC,MPART)
!                   CROSSASSIGNMENT OF WAM PARTITIONINGS TO SAR
!                   PARTITIONINGS.
      REAL MEANWE(NSPEC,MPART),
     &     MEANWANG(NSPEC,MPART),
     &     MEANWFRE(NSPEC,MPART),
     &     MEANSE(NSPEC,MPART),
     &     MEANSANG(NSPEC,MPART), 
     &     MEANSFRE(NSPEC,MPART)
!                   MEAN VALUES FOR WAM AND SAR PARTITIONINGS.
      REAL ADJUST(NSPEC,MPART)
!                   FACTOR FOR ENERGY ADJUSTMENT.
      INTEGER INTROTATE(NSPEC,MPART)
!                   ROTAION PARAMETER(NUMBER OF GRID POINTS)
      REAL FRACROTATE(NSPEC,MPART)
!                   FRACTIONAL ROTATION PARAMETER.
      REAL STRETCH(NSPEC,MPART)
!                   FREQUENCY TRANSFORMATION PARAMETER.
      INTEGER INTLOGSTRETCH(NSPEC,MPART)
!                   FREQUENCY TRANSFORMATION PARAMETER(GRID)
      REAL COPO(NFRE), COPODIFF(NFRE)
      REAL FRATIO 

!     LOCAL VARIABLES:
!     ----------------

      INTEGER ISPEC, IPARTW, IPARTS, IANG, IFRE
!                    LOOP INDICES.
      INTEGER JANG, JANGP1, JFRE, JFREP1
!                    GRIDPOINT INDICES AFTER INTERPOLATED
!                    JANGP1 = "JANG PLUS 1" , JFREP1 = "JFRE PLUS 1"
      REAL AANG, BANG, AFRE, BFRE
!                    FOR INTERPOLATION
      REAL ROTATE
!                    THE NUMBER OF GRIDPOINTS THE PARTITIONING 
!                    MUST BE ROTATED
      REAL LOGCO, ENERGY, STRETCHFAC

!-------------------------------------------------------------------

!     INITIALIZATION:
!     ---------------

      LOGCO = LOG(FRATIO)
      COPO(1) = FRATIO 
      DO IFRE = 2,NFRE
        COPO(IFRE) = COPO(IFRE-1) * FRATIO 
        COPODIFF(IFRE) = COPO(IFRE) - COPO(IFRE-1)
      END DO

      DO ISPEC = 1,NSPEC
        INTROTATE(ISPEC,1)=0
        DO IPARTW = 1,NPARTW(ISPEC)
          IPARTS = CORRELTAW(ISPEC,IPARTW)
          IF( IPARTS.EQ.0 )THEN
            ADJUST(ISPEC,IPARTW) = 1.
            INTROTATE(ISPEC,IPARTW) = 0
            FRACROTATE(ISPEC,IPARTW) = 0.
            STRETCH(ISPEC,IPARTW) = 1.
            INTLOGSTRETCH(ISPEC,IPARTW) = 0
          ELSE
            ADJUST(ISPEC,IPARTW)
     &        = MEANSE(ISPEC,IPARTS) / MEANWE(ISPEC,IPARTW)
            ROTATE
     &        = ( MEANSANG(ISPEC,IPARTS) - MEANWANG(ISPEC,IPARTW) )
            ROTATE=ROTATE/ DELTH
            INTROTATE(ISPEC,IPARTW) = INT(ROTATE)
            FRACROTATE(ISPEC,IPARTW) = ROTATE - INT(ROTATE)
            STRETCH(ISPEC,IPARTW)
     &        = MEANSFRE(ISPEC,IPARTS) / MEANWFRE(ISPEC,IPARTW)
            INTLOGSTRETCH(ISPEC,IPARTW)
     &        = INT( LOG(STRETCH(ISPEC,IPARTW))/LOGCO+10000 ) - 10000
          END IF
        END DO
      END DO

!     ROTATE, STRETCH AND ADJUST THE WAM FIRST GUESS PARTITIONINGS
!     IN ONE STEP AND ADD THEM TO ARRAY SPECT:
!     ------------------------------------------------------------

      DO JANG = 1,NANG
        DO JFRE = 1,NFRE
          DO ISPEC = 1,NSPEC
            SPECW(ISPEC,JANG,JFRE) = MAX(SPECW(ISPEC,JANG,JFRE),
     &                                   0.1e-15)
            SPECT(ISPEC,JANG,JFRE) = GAPVAL
          END DO
        END DO
      END DO

      DO IANG = 1,NANG
        DO IFRE = 1,NFRE
          DO ISPEC = 1,NSPEC
            IPARTW = PARTW(ISPEC,IANG,IFRE)           
            JANG = IANG + INTROTATE(ISPEC,IPARTW)
            IF( JANG.GT.NANG ) JANG = JANG - NANG
            IF( JANG.LT.1 ) JANG = JANG + NANG
            IF((INTROTATE(ISPEC,IPARTW)+
     &          FRACROTATE(ISPEC,IPARTW)).GT.0) THEN
             JANGP1 = JANG + 1
            ELSE
             JANGP1 = JANG - 1
            END IF
            IF( JANGP1.GT.NANG ) JANGP1 = JANGP1 - NANG
            IF( JANGP1.LT.1 ) JANGP1 = JANGP1 + NANG
            JFRE = IFRE + INTLOGSTRETCH(ISPEC,IPARTW)
            IF(JFRE.LT.1) THEN
              JFRE=1
              STRETCHFAC=1.
            ELSE
              STRETCHFAC=STRETCH(ISPEC,IPARTW)
            ENDIF
            JFREP1 = JFRE + 1
            IF( JFRE.GE.1 .AND. JFREP1.LE.NFRE )THEN
              BANG = ABS(FRACROTATE(ISPEC,IPARTW))
              AANG = 1 - BANG
              AFRE = ( COPO(JFREP1)
     &            - COPO(IFRE) * STRETCHFAC )
     &            / COPODIFF(JFREP1)    
              BFRE = 1 - AFRE
              ENERGY = SPECW(ISPEC,IANG,IFRE) * ADJUST(ISPEC,IPARTW)
              SPECT(ISPEC,JANG  ,JFRE  )
     &            = SPECT(ISPEC,JANG  ,JFRE  ) +
     &                   ENERGY * AANG * AFRE 
              SPECT(ISPEC,JANG  ,JFREP1)
     &            = SPECT(ISPEC,JANG  ,JFREP1) +
     &                   ENERGY * AANG * BFRE 
              SPECT(ISPEC,JANGP1,JFRE  ) 
     &            = SPECT(ISPEC,JANGP1,JFRE  ) +
     &                   ENERGY * BANG * AFRE 
              SPECT(ISPEC,JANGP1,JFREP1) 
     &            =  SPECT(ISPEC,JANGP1,JFREP1) +
     &                   ENERGY * BANG * BFRE 
            END IF 
          END DO
        END DO
      END DO

!     TRANSFER SPECT TO SPECW
!     -----------------------

      DO IANG = 1,NANG
        DO IFRE = 1,NFRE
          DO ISPEC = 1,NSPEC
            SPECW(ISPEC,IANG,IFRE) = SPECT(ISPEC,IANG,IFRE)
          END DO
        END DO
      END DO

      RETURN

      END SUBROUTINE TRANSPART
