      SUBROUTINE TRANSMEANS(
!                     DIMENSIONS
     &                  NSPEC ,  NANG , NFRE ,
!                     INPUT
     &                  SPECW ,
!                     OUTPUT
     &                  MEANTE , MEANTANG , MEANTFRE ,
!                     WORK SPACE
     &                  SUM0 , SUMX0 , SUMY0 )

!-------------------------------------------------------------------

!     PURPOSE:
!     --------

!     MEAN VALUES OF TRANSFORMED SPECTRA
!     COMPUTES MEAN ENERGY,MEAN ANGLE AND MEAN FREQUENCY

!-------------------------------------------------------------------
!     MODULE :
!     --------

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,COSTH    ,
     &           SINTH     ,FRTAIL   ,WETAIL
      USE YOWPCONS , ONLY : EPSMIN

!-------------------------------------------------------------------

      IMPLICIT NONE

!     INTERFACE:
!     ----------

      INTEGER NSPEC, NANG, NFRE
!                     DIMENSIONS.
      REAL SPECW(0:NSPEC,NANG,NFRE)
!                     2D SPECTRUM.
      REAL MEANTE(NSPEC),
     &     MEANTANG(NSPEC), 
     &     MEANTFRE(NSPEC)
!                     MEAN PARAMETERS.
      REAL SUM0(NSPEC), SUMX0(NSPEC), SUMY0(NSPEC)
!                     WORK SPACE FOR TEMPORARY SUMS.   
!     local variables:
!     ----------------

      INTEGER ISPEC, IANG, IFRE
!                     LOOP VARIABLES.
      REAL FACTOR

!-------------------------------------------------------------------
!     INITIALIZE:
!     -----------

!     COMPUTE MEAN ENERGIES:
!     ----------------------

      DO ISPEC = 1,NSPEC
        MEANTE(ISPEC) = EPSMIN 
      END DO
      DO IFRE = 1,NFRE
        DO ISPEC = 1,NSPEC
          SUM0(ISPEC) = 0.
        END DO
        DO IANG = 1,NANG
          DO ISPEC = 1,NSPEC
            SUM0(ISPEC) = SUM0(ISPEC) + SPECW(ISPEC,IANG,IFRE)
          END DO
        END DO
        DO ISPEC = 1,NSPEC
          MEANTE(ISPEC) = MEANTE(ISPEC) + DFIM(IFRE) * SUM0(ISPEC)
        END DO
      END DO

!     ADD TAIL CONTRIBUTION
      FACTOR = WETAIL*FR(NFRE)*DELTH
      DO ISPEC = 1,NSPEC
        MEANTE(ISPEC) = MEANTE(ISPEC) + FACTOR * SUM0(ISPEC)
      END DO

!     COMPUTE MEAN ANGLES:
!     --------------------

      DO ISPEC = 1,NSPEC
        SUMX0(ISPEC) = 0.
        SUMY0(ISPEC) = 0.
      END DO

      DO IANG = 1,NANG
        DO ISPEC = 1,NSPEC
          SUM0(ISPEC) = 0.
        END DO
        DO IFRE=1,NFRE
          DO ISPEC=1,NSPEC
            SUM0(ISPEC)
     &          = SUM0(ISPEC) + SPECW(ISPEC,IANG,IFRE) * DFIM(IFRE)
          END DO
        END DO
        DO ISPEC = 1,NSPEC
          SUMX0(ISPEC) = SUMX0(ISPEC) + COSTH(IANG) * SUM0(ISPEC)
          SUMY0(ISPEC) = SUMY0(ISPEC) + SINTH(IANG) * SUM0(ISPEC)
        END DO
      END DO

      DO ISPEC = 1,NSPEC
        IF( SUMY0(ISPEC).NE.0..OR.SUMX0(ISPEC).NE.0. )THEN
          MEANTANG(ISPEC) = ATAN2(SUMY0(ISPEC),SUMX0(ISPEC))
        ELSE
          MEANTANG(ISPEC) = 0.
        END IF
      END DO

!     COMPUTE MEAN FREQUENCIES:
!     -------------------------

      DO ISPEC = 1,NSPEC
        MEANTFRE(ISPEC) = EPSMIN
      END DO

      DO IFRE = 1,NFRE
        DO ISPEC = 1,NSPEC
          SUM0(ISPEC) = 0.
        END DO
        DO IANG = 1,NANG
          DO ISPEC = 1,NSPEC
            SUM0(ISPEC) = SUM0(ISPEC) + SPECW(ISPEC,IANG,IFRE)
          END DO
        END DO
        FACTOR = DFIM(IFRE) / FR(IFRE)
        DO ISPEC = 1,NSPEC
          MEANTFRE(ISPEC) = MEANTFRE(ISPEC) + FACTOR * SUM0(ISPEC)
        END DO
      END DO

      FACTOR = FRTAIL*DELTH
      DO ISPEC = 1,NSPEC
        MEANTFRE(ISPEC) = MEANTFRE(ISPEC) + FACTOR * SUM0(ISPEC)
        IF( MEANTFRE(ISPEC).NE.0. )THEN
          MEANTFRE(ISPEC) = MEANTE(ISPEC) / MEANTFRE(ISPEC)
        ELSE
          MEANTFRE(ISPEC) = 0.
        END IF
      END DO

      RETURN

      END SUBROUTINE TRANSMEANS
