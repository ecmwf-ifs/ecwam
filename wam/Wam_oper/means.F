      SUBROUTINE MEANS(
!                   DIMENSIONS
     &                NSPEC , MPART , NPART , NANG, NFRE ,
!                   INPUT
     &                SPEC , PART ,
!                   OUTPUT
     &                MEANE , MEANANG , MEANFRE ,
!                   WORK SPACE.
     &                SUM0 , SUMX0 , SUMY0 )

!-----------------------------------------------------------------------

!     PURPOSE:
!     --------

!     COMPUTES MEAN DIRECTION AND MEAN FREQUENCY OF 2D SPECTRA.
!     FOR THE CALCULATION OF THE SPECTRAL PARTITIONING.

!-----------------------------------------------------------------------

!     MODULE :
!     --------

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,COSTH    ,
     &           SINTH     ,FRTAIL
      USE YOWPCONS , ONLY : EPSMIN

!-----------------------------------------------------------------------

      IMPLICIT NONE

!     INTERFACE:
!     ----------

      INTEGER NSPEC, MPART, NPART(NSPEC), NANG, NFRE
!                 ARRAY DIMENSIONS.
      REAL SPEC(0:NSPEC,NANG,NFRE)
!                 2D SPECTRA.
      INTEGER PART(NSPEC,NANG,NFRE)
!          PART(..) GIVES NUMBER OF PARTITIONING FOR EACH SPECTRAL BIN.
      REAL MEANE(NSPEC,MPART),
     &     MEANANG(NSPEC,MPART), MEANFRE(NSPEC,MPART)
!                 INTEGRATED VALUES OF PARTITIONINGS.
      REAL SUM0(NSPEC,MPART),
     &     SUMX0(NSPEC,MPART),SUMY0(NSPEC,MPART)
!                 FOR TEMPORARY SUMS

!     LOCAL VARIABLES:
!     ----------------

      INTEGER ISPEC, IPART, IANG, IFRE
!                   LOOP INDEXES.
      REAL :: FACTOR(NFRE)
      REAL FACTOR_TAIL

!-----------------------------------------------------------------------

!     1. COMPUTE MEAN DIRECTIONS.
!     --------------------------

      DO IPART = 1,MPART
        DO ISPEC = 1,NSPEC
          SUMX0(ISPEC,IPART) = 0
          SUMY0(ISPEC,IPART) = 0
        END DO
      END DO

      DO IANG = 1,NANG
        DO IPART = 1,MPART
          DO ISPEC = 1,NSPEC
            SUM0(ISPEC,IPART) = 0
          END DO
        END DO
        DO IFRE=1,NFRE
          DO ISPEC=1,NSPEC
            SUM0(ISPEC,PART(ISPEC,IANG,IFRE))
     &          = SUM0(ISPEC,PART(ISPEC,IANG,IFRE))
     &          + SPEC(ISPEC,IANG,IFRE) * DFIM(IFRE)
          END DO
        END DO
        DO IPART = 1,MPART
          DO ISPEC = 1,NSPEC
            SUMX0(ISPEC,IPART) = SUMX0(ISPEC,IPART) 
     &          + COSTH(IANG) * SUM0(ISPEC,IPART)
            SUMY0(ISPEC,IPART) = SUMY0(ISPEC,IPART) 
     &          + SINTH(IANG) * SUM0(ISPEC,IPART)
          END DO
        END DO
      END DO

      DO ISPEC = 1,NSPEC
        DO IPART = 1,NPART(ISPEC)
          MEANANG(ISPEC,IPART)
     &        = ATAN2(SUMY0(ISPEC,IPART),SUMX0(ISPEC,IPART))
        END DO
      END DO

!----------------------------------------------------------------------

!    2. COMPUTE MEAN FREQUENCY.
!    --------------------------

      DO IPART = 1,MPART
        DO ISPEC = 1,NSPEC
          MEANFRE(ISPEC,IPART) = EPSMIN
        END DO
      END DO

      DO IFRE = 1,NFRE
        FACTOR(IFRE) = DFIM(IFRE) / FR(IFRE)
      END DO

      DO IFRE = 1,NFRE
        DO IPART = 1,MPART
          DO ISPEC = 1,NSPEC
            SUM0(ISPEC,IPART) = 0
          END DO
        END DO
        DO IANG = 1,NANG
          DO ISPEC = 1,NSPEC
            SUM0(ISPEC,PART(ISPEC,IANG,IFRE))
     &          = SUM0(ISPEC,PART(ISPEC,IANG,IFRE))
     &          + SPEC(ISPEC,IANG,IFRE)
          END DO
        END DO
        DO IPART = 1,MPART
          DO ISPEC = 1,NSPEC
            MEANFRE(ISPEC,IPART) = MEANFRE(ISPEC,IPART)
     &          + FACTOR(IFRE)  * SUM0(ISPEC,IPART)
          END DO
        END DO
      END DO

!     ADD TAIL CONTRIBUTION.

      FACTOR_TAIL = FRTAIL*DELTH

      DO IPART = 1,MPART
        DO ISPEC = 1,NSPEC
          MEANFRE(ISPEC,IPART) = MEANFRE(ISPEC,IPART)
     &        + FACTOR_TAIL * SUM0(ISPEC,IPART)
          MEANFRE(ISPEC,IPART)
     &        = MEANE(ISPEC,IPART) / MEANFRE(ISPEC,IPART)
        END DO
      END DO

      RETURN

      END SUBROUTINE MEANS
