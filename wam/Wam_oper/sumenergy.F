      SUBROUTINE SUMENERGY( 
!                  DIMENSIONS
     &               NSPEC , MPART , NPART , NANG , NFRE ,
!                  INPUT
     &               SPEC , PART ,
!                  OUTPUT
     &               MEANE , SPREAD , 
!                  WORK SPACE
     &               SUMX , SUMY , SUMXX , SUMYY ,
     &               SUM0 , SUMX0 , SUMY0 , SUMXX0 , SUMYY0 )

!-----------------------------------------------------------------------

!     PURPOSE:
!     --------

!     COMPUTE MEAN ENERGY AND SQUARES OF SPREAD FOR ALL PARTITIONINGS.

!-----------------------------------------------------------------------

!     MODULES:
!     --------

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,COSTH    ,
     &           SINTH     ,WETAIL
      USE YOWPCONS , ONLY : EPSMIN

!-----------------------------------------------------------------------

      IMPLICIT NONE

!     INTERFACE:
!     ----------

      INTEGER NSPEC, MPART, NPART(NSPEC), NANG, NFRE
!                    ARRAY DIMENSIONS.
      REAL SPEC(0:NSPEC,NANG,NFRE)
!                    2D SPECTRA.
      INTEGER PART(NSPEC,NANG,NFRE)
!                    PART(..) GIVES NUMBER OF PARTITIONING FOR EACH 
!                    SPECTRAL BIN.
      REAL MEANE(NSPEC,MPART),
     &     SPREAD(NSPEC,MPART),
!                    MEAN ENERGY, SPECTRAL SPREAD.
     &    SUMX(NSPEC,MPART), SUMY(NSPEC,MPART),
     &    SUMXX(NSPEC,MPART), SUMYY(NSPEC,MPART)
!                    WORK SPACE.
      REAL SUM0(NSPEC,MPART),
     &     SUMX0(NSPEC,MPART),SUMXX0(NSPEC,MPART),
     &     SUMY0(NSPEC,MPART),SUMYY0(NSPEC,MPART)
!                    TEMPORARY SPACE.

!     LOCAL VARIABLES:
!     ----------------

      INTEGER ISPEC, IPART, IANG, IFRE
!                  LOOP INDEXES.
      REAL FX, FY
!                  WAVE NUMBERS.
      REAL FACTOR 

!-----------------------------------------------------------------------

!     COMPUTE SUMS AND SQUARE SUMS
!     ---------------------------

      DO IPART = 1,MPART
        DO ISPEC = 1,NSPEC
          MEANE(ISPEC,IPART) = EPSMIN 
          SUMX(ISPEC,IPART) = 0.
          SUMY(ISPEC,IPART) = 0.
          SUMXX(ISPEC,IPART) = 0.
          SUMYY(ISPEC,IPART) = 0.
        END DO
      END DO

      DO IFRE = 1,NFRE

        DO IPART = 1,MPART
          DO ISPEC = 1,NSPEC
            SUM0(ISPEC,IPART) = 0.
            SUMX0(ISPEC,IPART) = 0. 
            SUMY0(ISPEC,IPART) = 0. 
            SUMXX0(ISPEC,IPART) = 0. 
            SUMYY0(ISPEC,IPART) = 0. 
          END DO
        END DO

        DO IANG = 1,NANG
          FX = FR(IFRE) * COSTH(IANG)
          FY = FR(IFRE) * SINTH(IANG)
          DO ISPEC = 1,NSPEC
            SUM0(ISPEC,PART(ISPEC,IANG,IFRE))
     &          = SUM0(ISPEC,PART(ISPEC,IANG,IFRE))
     &          + SPEC(ISPEC,IANG,IFRE)
            SUMX0(ISPEC,PART(ISPEC,IANG,IFRE))
     &          = SUMX0(ISPEC,PART(ISPEC,IANG,IFRE))
     &          + FX * SPEC(ISPEC,IANG,IFRE)
            SUMXX0(ISPEC,PART(ISPEC,IANG,IFRE))
     &          = SUMXX0(ISPEC,PART(ISPEC,IANG,IFRE))
     &          + FX**2 * SPEC(ISPEC,IANG,IFRE)
            SUMY0(ISPEC,PART(ISPEC,IANG,IFRE))
     &          = SUMY0(ISPEC,PART(ISPEC,IANG,IFRE))
     &          + FY * SPEC(ISPEC,IANG,IFRE)
            SUMYY0(ISPEC,PART(ISPEC,IANG,IFRE))
     &          = SUMYY0(ISPEC,PART(ISPEC,IANG,IFRE))
     &          + FY**2 * SPEC(ISPEC,IANG,IFRE)
          END DO
        END DO

        DO IPART = 1,MPART
          DO ISPEC = 1,NSPEC
            MEANE(ISPEC,IPART) = MEANE(ISPEC,IPART)
     &          + DFIM(IFRE) * SUM0(ISPEC,IPART)
            SUMX(ISPEC,IPART) = SUMX(ISPEC,IPART)
     &          + DFIM(IFRE) * SUMX0(ISPEC,IPART)
            SUMY(ISPEC,IPART) = SUMY(ISPEC,IPART)
     &          + DFIM(IFRE) * SUMY0(ISPEC,IPART)
            SUMXX(ISPEC,IPART) = SUMXX(ISPEC,IPART)
     &          + DFIM(IFRE) * SUMXX0(ISPEC,IPART)
            SUMYY(ISPEC,IPART) = SUMYY(ISPEC,IPART)
     &          + DFIM(IFRE) * SUMYY0(ISPEC,IPART)
          END DO
        END DO
      END DO

!     ADD TAIL ENERGY
!     ---------------

      FACTOR = WETAIL*FR(NFRE)*DELTH

      DO IPART = 1,MPART
        DO ISPEC = 1,NSPEC
          MEANE(ISPEC,IPART) = MEANE(ISPEC,IPART)
     &        + FACTOR * SUM0(ISPEC,IPART)
        END DO
      END DO

!     COMPUTE SPREAD
!     --------------

      DO IPART = 1,MPART
        DO ISPEC = 1,NSPEC
          SPREAD(ISPEC,IPART) = MAX
     &        ( SUMXX(ISPEC,IPART)/MEANE(ISPEC,IPART)
     &        + SUMYY(ISPEC,IPART)/MEANE(ISPEC,IPART)
     &        - (SUMX(ISPEC,IPART)/MEANE(ISPEC,IPART))**2
     &        - (SUMY(ISPEC,IPART)/MEANE(ISPEC,IPART))**2
     &        , 0. )
        END DO
      END DO

      RETURN

      END SUBROUTINE SUMENERGY
