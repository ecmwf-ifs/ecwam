      SUBROUTINE TRANSPART( 
!                    DIMENSIONS
     &                 NSPEC , MPART , NPARTW ,NANG , NFRE ,
!                    INPUT/OUTPUT
     &                 SPECW , NW1D, 
!                    INPUT
     &                 PARTW ,CORRELTAW ,MEANWE,MEANWANG,MEANWFRE,
     &                 MEANSE , MEANSANG , MEANSFRE ,
     &                 FRATIO,
!                    OUTPUT
     &                 ADJUST , INTROTATE , FRACROTATE , 
     &                 STRETCH , INTLOGSTRETCH )

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
      USE YOWGRID  , ONLY : IJS      ,IJL
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWSPEC,   ONLY : FL1
      USE YOWSTAT  , ONLY : NWAM_BLKS

!-------------------------------------------------------------------

      IMPLICIT NONE

!     INTERFACE:
!     ----------

      INTEGER NSPEC, MPART, NPARTW(NSPEC), NANG, NFRE
!                   DIMENSIONS.
      REAL :: SPECW(0:NSPEC,NANG,NFRE) 
!                   WAM SPECTRM BEFORE AND AFTER CORRECTION 

      INTEGER :: NW1D(NINF-1:NSUP,1)

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

      REAL FRATIO 

!     LOCAL VARIABLES:
!     ----------------

      INTEGER ISPEC, IPARTW, IPARTS, IANG, IFRE, JP, JBL
      INTEGER NPROMA, JKGLO, KIJS, KIJL
!                    LOOP INDICES.
      INTEGER JANG, JANGP1, JFRE, JFREP1
!                    GRIDPOINT INDICES AFTER INTERPOLATED
!                    JANGP1 = "JANG PLUS 1" , JFREP1 = "JFRE PLUS 1"
      REAL AANG, BANG, AFRE, BFRE
!                    FOR INTERPOLATION
      REAL ROTATE
!                    THE NUMBER OF GRIDPOINTS THE PARTITIONING 
!                    MUST BE ROTATED
      REAL COPO(NFRE+1), COPODIFF(NFRE+1)
!                   POWERS OF FRATIO AND DIFFERENCES OF POWERS OF 
!                   FRATIO, ONLY USED IN SUBROUTINE TRANSFORM.
      REAL LOGCO, ENERGY, STRETCHFAC

      REAL :: SPECT(NSPEC,NANG,NFRE+1)
!                   WORK SPACE.

!-------------------------------------------------------------------

!     INITIALIZATION:
!     ---------------

      LOGCO = LOG(FRATIO)
      COPO(1) = FRATIO 
      DO IFRE = 2,NFRE+1
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
     &        = (MEANSE(ISPEC,IPARTS) * MEANWFRE(ISPEC,IPARTW)) /
     &          (MEANWE(ISPEC,IPARTW) * MEANSFRE(ISPEC,IPARTS))
!!!!! ??????
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

!     restore the unsmoothed model spectra
!!!  will not work if multi block !!!

      JBL=1
      NPROMA=(IJL(JBL)-IJS(JBL)+1)/NWAM_BLKS+1
!$OMP     PARALLEL DO SCHEDULE(STATIC) 
!$OMP+    PRIVATE(JKGLO,KIJS,KIJL,JP,JANG,JFRE)
      DO JKGLO = IJS(JBL),IJL(JBL),NPROMA
        KIJS=JKGLO
        KIJL=MIN(KIJS+NPROMA-1,IJL(JBL))
        DO JP = KIJS,KIJL
          IF(NW1D(JP,JBL).GT.0) THEN
            DO JANG = 1,NANG
              DO JFRE = 1,NFRE
                SPECW(NW1D(JP,JBL),JANG,JFRE)=FL1(JP,JANG,JFRE)
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO
!$OMP     END PARALLEL DO


      DO JANG = 1,NANG
        DO JFRE = 1,NFRE+1
          DO ISPEC = 1,NSPEC
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
!!!!              JFRE=1
!!! it will artificially put energy in the lowest bins otherwise
              STRETCHFAC=1.
            ELSE
              STRETCHFAC=STRETCH(ISPEC,IPARTW)
            ENDIF
            JFREP1 = JFRE + 1

            IF( JFRE.GE.1 .AND. JFRE.LE.NFRE )THEN
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

!     ADD UNASSIGNED SAR PARTITIONS
!     -----------------------------

!!    no longer in use

      RETURN

      END SUBROUTINE TRANSPART
