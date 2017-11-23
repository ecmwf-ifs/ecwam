      SUBROUTINE TUSTREAS(
!           DIMENSIONS
     &            NSPEC, MPART, NPARTW, NANG, NFRE,
!           INPUT (FOR SPECTRUM ALSO OUTPUT)
     &            SPECW, NW1D, PARTW, CORRELTAW,
     &            MEANWE, MEANWANG, MEANWFRE,
     &            MEANSE, MEANSANG, MEANSFRE,
!           OUTPUT
     &            MEANTE, MEANTANG, MEANTFRE,
!           INPUT PARAMETERS 
     &            FRATIO , FR , COSTH , SINTH, CONTROL,
     &            IU06)

!-------------------------------------------------------------------
!     MODULE :
!     --------

      USE YOWMPP   , ONLY : NINF     ,NSUP

!-------------------------------------------------------------------

      IMPLICIT NONE
#include "fillgaps.intfb.h"
#include "transmeans.intfb.h"
#include "transpart.intfb.h"

!     PURPOSE:
!     --------

!     FOR DATA ASSIMILATION WITHOUT SUPERIMPOSING EXTRA SAR WAVE 
!     SYSTEMS ONTO FIRST GUESS SPECTRA.

!     TRANSFORMS PARTITIONINGS OF WAM FIRST GUESS SPECTRA TO MEAN 
!     PARAMETERS OF INVERTED SAR SPECTRA.

!     TO MAKE THIS ROUTINE AS FAST AS POSSIBLE, THE LOOP OVER 
!     ALL SPECTRA ARE VECTORIZED INNER LOOP.

!     METHOD:
!     -------

!     THE WAM FIRST GUESS PARTITIONINGS ARE ROTATED,
!     STRETCHED TO MATCH THE MEAN ENERGY, MEAN DIRECTION,
!     AND MEAN FREQUENCY OF THE PARTITIONINGS OF THE 
!     INVERTED SAR SPECTRA. THE PARTITIONINGS ARE SUPER IMPOSED 
!     TO DERIVE ONE COMBINED WAVE SPECTRUM. 
!     GAPS BETWEEN TWO PARTITIONINGS ARE 
!     INTERPOLATED BY A 2D PARABOLIC INTERPOLATION. FOR OVERLAPPING 
!     PARTITIONINGS THE MAXIMAL VALUE IS TAKEN.

!     LOOPS OVER ALL SPECTRA ARE VECTORIZED.
  
!     BEFORE THIS ROUTINE IS CALLED THE CROSS CORRELATION OF 
!     WAM FIRSRT GUESS AND INVERTED SAR PARTITIONINGS
!     HAS TO BE CARRIED OUT (CORRELTAW)!

!     
!     PARAMETERS FOR DIMENSIONS:
!     --------------------------

!     AUTHOR:
!     -------

!     SUSANNE HASSELMANN, 1993 MPI FUER METEOROLOGIE HAMBURG
!     JUERGEN WASZKEWTIZ, 1993 MPI FUER METEOROLOGIE HAMBURG

!     INTERFACE:
!     ----------

      INTEGER NSPEC
!                     (INPUT) THE ACTUAL NUMBER OF SPECTRA.
      INTEGER MPART
!                     (INPUT) THE MAXIMAL NUMBER (=DIMENSION)
!                     OF PARTITIONINGS.
      INTEGER NPARTW(NSPEC)
!                     (INPUT) ARRAY OF THE ACTUAL NUMBERS OF 
!                     PARTITIONINGS OF THE ISPEC(TH) WAM FIRST GUESS
!                     SPECTRUM.
      INTEGER NANG
!                     (INPUT) THE ACTUAL NUMBER OF SPECTRAL
!                     DIRECTION AND THE DIMENSION.
      INTEGER NFRE
!                     THE ACTUAL NUMBER OF SPECTRAL FREQUENCY BINS
!                      AND THE DIMENSION.
      REAL SPECW(0:NSPEC,NANG,NFRE)
!                     (INPUT,OUTPUT) AT INPUT: WAM FIRST GUESS 
!                     SPECTRA.
!                     AT OUTPUT: CORRECTED WAM SPECTRA
      INTEGER :: NW1D(NINF-1:NSUP,1)

      INTEGER PARTW(NSPEC,NANG,NFRE)
!                     (INPUT) PARTITIONINGS OF WAM FIRST GUESS 
!                     SPECTRA, WHERE PARTW(ISPEC,IANG,IFRE) IS 
!                     THE NUMBER OF THE PARTITION OF THE WAM 
!                     FIRST GUESS SPECTRUM ISPEC
!                     AT DIRECTION IANG AND FREQUENCY IFRE
      INTEGER CORRELTAW(NSPEC,MPART)
!                     (INPUT) CROSS ASSIGNMENT OF WAM AND SAR 
!                     PARTITIONINGS, WHERE CORRELTAW(ISPEC,IPARTW) 
!                     MEANS THE PARTITIONING IPARTW
!                     OF THE WAM FIRST GUESS SPECTRUM ISPEC
!                     IS ASSIGNED TO PARTITIONING CORRELTAW(IPARTW) 
!                     OF THE INVERTED SAR SPECTRUM
!                     (A ZERO MEANS, THAT THERE IS NO CROSS 
!                     ASSIGNMENT.
      REAL MEANWE(NSPEC,MPART)
!                     (INPUT) MEAN ENERGIES OF PARTITIONINGS OF WAM 
!                     FIRST GUESS SPECTRA, WHERE MEANWE(ISPEC,IPART)
!                     IS THE ENERGY OF SPECTRUM ISPEC 
!                     PARTITIONING IPART
      REAL MEANWANG(NSPEC,MPART)
!                     (INPUT) MEAN DIRECTIONS IN RADIANCE OF 
!                     PARTITIONINGS OF WAM FIRST GUESS SPECTRA.
      REAL MEANWFRE(NSPEC,MPART)
!                     (INPUT) MEAN FREQUENCIES OF PARTITIONINGS
!                     OF WAM FIRST GUESS SPECTRA.
      REAL MEANSE(NSPEC,MPART)
!                     (INPUT) MEAN ENERGIES OF PARTITIONINGS
!                     OF INVERTED SAR SPECTRA, 
      REAL MEANSANG(NSPEC,MPART)
!                     (INPUT) MEAN DIRECTIONS IN RADIANCE OF 
!                     PARTITIONINGS OF INVERTED SAR SPECTRA.
      REAL MEANSFRE(NSPEC,MPART)
!                     (INPUT) MEAN FREQUENCIES OF PARTITIONINGS
!                     OF INVERTED SAR SPECTRAS.
      REAL MEANTE(NSPEC)
!                     (OUTPUT) ARRAY OF TOTAL MEAN ENERGIES OF 
!                     TRANSFORMED WAM FIRST GUESS SPECTRA
      REAL MEANTANG(NSPEC)
!                     (OUTPUT) ARRAY OF MEAN DIRECTIONS IN RADIAN
!                     OF TRANSFORMED WAM FIRST GUESS SPECTRA
      REAL MEANTFRE(NSPEC)
!                     (OUTPUT) MEAN FREQUENCIES OF TRANSFORMED
!                     WAM FIRST GUESS SPECTRA.
      REAL FRATIO 
!                     (INPUT) RATIO BETWEEN FREQUENCIES 
!                     (FRATIO=FRE(2)/FRE(1))
      REAL FR(NFRE)
!                       FREQUENCIES, IS READ IN SUBROUTINE 
!                       READPRES AND USED IN MEANST AND FILLGAPS
      REAL COSTH(NANG),
     &     SINTH(NANG)
!                       COSINE AND SINE, READ IN SUBROUTINE 
!                       READPRES AND USED ONLY IN SUBROUTINE 
!                       MEANST.
      LOGICAL CONTROL
!                     (INPUT) CONTROLS  OUTPUT MESSAGES ON UNIT IU06 ?
!                     .TRUE. = YES , .FALSE. = NO
      INTEGER IU06
!                     USER OUTPUT UNIT.

!     EXTERNALS:
!     ----------

!     AVOIDPEAKS     - REDUCES SIZE OF FRAMES AROUND GAPS TO 
!                      AVOID PEAKS INSIDE A FRAME.
!     F04JGE         - FIND SOLUTION OF LINEAR LEAST SQUARES 
!                      PROBLEM (NAG)
!     FILLGAPS       - FILL GAPS BETWEEN PARTITIONINGS OF CORRECTED 
!                      SPECTRA WITH A 2D PARABOLIC INTERPOLATION.
!     GAPINTERPOL    - 2D PARABOLIC INTERPOLATION.
!     MAKEFRAMES     - BUILDS FRAMES AROUND GAPS BETWEEN WAVE 
!                      SYSTEMS.
!     TRANSPART      - TRANSFORMS SPECTRAL PARTITIONINGS.
!     TRANSMEANS     - COMPUTES MEANS OF TRANSFORMED TOTAL SPECTRA

!     SUBROUTINE TREE:

!     TUSTREAS
!       |
!       TRANSPART
!       |
!       FILLGAPS
!       | |
!       | MAKEFRAMES
!       | |
!       | AVOIDPEAKS
!       | |
!       | GAPINTERPOL
!       |   |
!       |   F04JGE
!       |
!       TRANSMEANS
!      
!       
!     VARIABLES:
!     ----------

      REAL ADJUST(NSPEC,MPART)
!                       ENERGY ADJUSTMENT PARAMETER.
      INTEGER INTROTATE(NSPEC,MPART)
!                       DIRECTIONAL GRIDPOINT ADJUSTMENT PARAMETER.
!                       (GRID POINTS), ONLY USED IN SUBROUTINE 
!                       TRANSFORM.
      REAL FRACROTATE(NSPEC,MPART)
!                       DIRECTIONAL GRIDPOINT ADJUSTMENT PARAMETER.
!                       (FRACTIONAL PART), ONLY USED IN SUBROUTINE 
!                       TRANSFORM
      REAL STRETCH(NSPEC,MPART)
!                       FREQUENCY ADJUSTMENT PARAMETER.
      INTEGER INTLOGSTRETCH(NSPEC,MPART)
!                       INTLOGSTRETCH = INT(LOG(STRETCH)/LOGCO), 
!                       ONLY USED IN SUBROUTINE TRANSFORM
      REAL FX(NANG,NFRE),
     &     FY(NANG,NFRE)
!                       TABLE OF FREQUENCY-DIRECTIONAL 
!                       POLAR COORDINATES.
!                       USED IN SUBROUTINE GAPINTERPOL AND TRANSMEANS
      REAL SUM0(NSPEC), SUMX0(NSPEC), SUMY0(NSPEC)
!                       TEMPORARILY STORAGE, ONLY USED IN 
!                       SUBROUTINE TRANSMEANS.
      INTEGER IANG, IFRE

!-------------------------------------------------------------------

      IF( CONTROL ) WRITE(IU06,*) 'CALL OF SUBROUTINE TUSTREAS'

!-------------------------------------------------------------------

!     2. INITIALIZATION.
!     ------------------

      DO IANG = 1,NANG
        DO IFRE = 1,NFRE
          FX(IANG,IFRE) = FR(IFRE) * COSTH(IANG)
          FY(IANG,IFRE) = FR(IFRE) * SINTH(IANG)
        END DO
      END DO

!-------------------------------------------------------------------

!     3. CORRECT PARTITIONINGS OF MODEL FIRST GUESS SPECTRA 
!        ACCORDING TO MEAN VALUES OF PARTITIONINGS OF SAR 
!        INVERTED SPECTRA.
!     ------------------------------------------------------

      IF( CONTROL ) THEN
        WRITE(IU06,*) 'CALLING SUBROUTINE TRANSPART'
        CALL FLUSH(IU06)
      ENDIF
      CALL TRANSPART( 
!                    DIMENSIONS
     &                 NSPEC , MPART , NPARTW , NANG , NFRE,
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

!     4. INTERPOLATE GAPS BETWEEN SUPERIMPOSED PARTITIONINGS BY 
!        2D PARABOLIC INTERPOLATION.
!     ---------------------------------------------------------

      IF( CONTROL ) THEN
        WRITE(IU06,*) 'CALLING SUBROUTINE FILLGAPS'
        CALL FLUSH(IU06)
      ENDIF
      CALL FILLGAPS( 
!                  DIMENSIONS
     &               NSPEC , MPART , NANG , NFRE ,
!                  INPUT (SPECW ALSO OUTPUT)
     &               SPECW , 
     &               FX , FY , FR , CONTROL , IU06 )

!-------------------------------------------------------------------

!     5. COMPUTE MEAN PARAMETERS OF TOTAL CORRECTED SPECTRA.
!     ------------------------------------------------------

      IF( CONTROL ) THEN
        WRITE(IU06,*) 'CALLING SUBROUTINE TRANSMEANS'
        CALL FLUSH(IU06)
      ENDIF
      CALL TRANSMEANS(
!                  DIMENSIONS
     &                   NSPEC ,  NANG , NFRE ,
!                  INPUT 
     &                   SPECW ,
!                  OUTPUT
     &                   MEANTE , MEANTANG , MEANTFRE ,
!                  WORK ARRAYS
     &                   SUM0 , SUMX0 , SUMY0 )

      IF( CONTROL ) THEN
        WRITE(IU06,*) 'END OF SUBROUTINE TUSTREAS'
        CALL FLUSH(IU06)
      ENDIF

      RETURN

      END SUBROUTINE TUSTREAS
