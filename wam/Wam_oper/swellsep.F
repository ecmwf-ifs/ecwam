      SUBROUTINE SWELLSEP( 
!             DIMENSIONS 
     &          NSPEC , MPART , NANG , NFRE ,
!             INPUT  
     &          SPEC , US , THW , TH0 , FR , COSTH , SINTH,
     &          LEVEL , CONTROL, LOWEST,
!             OUTPUT
     &          NPART , PART , PARTINFO , MAXMPART,
     &          MEANE , MEANANG , MEANFRE, PKFRE)


!-------------------------------------------------------------------

!     PURPOSE:
!     --------

!     PARTITIONING OF AN ARRAY OF 2-D WAVE SPECTRA INTO WINDSEA AND
!     SWELL PARTITIONS.

!     METHODS:
!     --------

!     MOUNTAINEERING CONCEPT: THE HUMP DOMAIN OF A SPECTRAL PEAK 
!     CONSISTS OF ALL SPECTRAL POINTS WHOSE PATHS OF STEEPEST ASCENT
!     LEADS TO THAT PEAK. (SUBROUTINE PARTITIONING)

!     PARTITIONINGS ARE COMBINED TO ONE IF(SUBROUTINE COMBINEPEAKS):
!     - THE PEAKS ARE ONLY ONE GRID POINT APART 
!             (SUBROUTINE ONESTEPAPART).
!     - THE PEAK LIES IN LAST TWO FREQUENCY BINS 
!             (SUBROUTINE LASTTWOBINS).
!     - THE DISTANCE BETWEEN PEAKS IS LESS THAN HALF THE 
!              SPREAD OF THE PARTITIONING (SUBROUTINE HALFSPREAD).
       
!     - THE MINIMUM ENERGY BETWEEN PEAKS IS >LEVEL * THE 
!       MINIMUM PEAK ENERGY (SUBROUTINE THRESHHOLD).
!     - THE MEAN ENERGY IS < LOWEST.
!     ALL WINDSEA PEAKS ARE COMBINED TO ONE.

!     REFERENCES:
!     -----------

!     S. HASSELMANN, K. HASSELMANN, C. BRUENING, 1993:
!     EXTRACTION OF WAVE SPECTRA FROM SAR IMAGE SPECTRA. IN:
!     DYNAMICS AND MODELING OF OCEAN WAVES. SECTION V.4.3.
!     EDITOR G. KOMEN. KLUWER & CO, NETHERLANDS.

!     MEMORANDUM FOR WAM GROUP BY KLAUS HASSELMANN 30/7/92.

!     AUTHOR:
!     -------

!     SUSANNE HASSELMANN, 1992
!     JUERGEN WASZKEWITZ, 1993
!     JEAN BIDLOT, 1999

!     EXTERNALS:
!     ----------

!     COMBINEPEAKS      -  COMBINE PEAKS.
!     DOCOMBINE         -  DO THE COMBINING OF TWO PARTITIONS TO 
!                          ONE..
!     HALFSPREAD        -  COMBINE PEAKS IF DISTANCE BETWEEN PEAKS 
!                          IS LESS THAN HALF OF SPREADS.
!     LASTTWOBINS       -  COMBINE PEAKS IN LAST TWO BINS WITH 
!                          CLOSEST PEAK.
!     LOWENERGY         -  COMBINE PARTITIONINGS WITH TOO LOW ENERGY
!     MEANS             -  COMPUTE MEAN ANGLE AND MEAN FREQUENCY.
!     ONESTEPAPART      -  COMBINE PEAKS WHICH ARE ONLY ONE GTID 
!                          POINT APART.
!     PARTITIONING      -  PARTITION SPECTRA.
!     SUMENERGY         -  COMPUTE ENERGY AND SPREAD.
!     THRESHHOLD        -  COMBINE PEAKS AFTER THRESSHOLD CHECK
!     WINDSEA           -  WINDSEA - SWELL SEPARATION.

!     SUBROUTINE TREE:

!     SWELLSEP
!       |
!       PARTITIONING
!       |
!       COMBINEPEAKS
!       | |
!       | ONESTEPAPART
!       | |
!       | SUMENERGY
!       | |
!       | WINDSEA
!       | | |
!       | | DOCOMBINE
!       | |
!       | LASTTWOBINS
!       | | |
!       | | DOCOMBINE
!       | |
!       | HALFSPREAD
!       | | |
!       | | DOCOMBINE
!       | |
!       | THRESHHOLD
!       | | |
!       | | DOCOMBINE
!       | |
!       | LOWENERGY
!       |   |
!       |   DOCOMBINE
!       |  
!       |
!       MEANS

!-----------------------------------------------------------------------

!     MODULES
!     -------

      USE YOWICE   , ONLY : FLMIN
      USE YOWTEST  , ONLY : IU06

!-----------------------------------------------------------------------

!     INTERFACE:
!     ----------

      IMPLICIT NONE

      INTEGER NSPEC  ! NUMBER OF SPECTRA TO BE PARTITIONED
      INTEGER MPART  ! MAXIMUM NUMBER OF PARTITIONING.
      INTEGER NPART(NSPEC)
!                      NUMBER OF PARTITIONINGS,
!                      WHERE NPART(ISPEC) IS THE NUMBER OF 
!                      PARTITIONINGS OF SPECTRUM ISPEC.
      INTEGER NANG   ! NUMBER OF SPECTRAL DIRECTIONS
      INTEGER NFRE   ! NUMBER OF SPECTRAL FREQUENCIES
      REAL SPEC(0:NSPEC,NANG,NFRE) ! SPECTRA
      INTEGER PART(NSPEC,NANG,NFRE)
!                      PARTITIONINGS, WHERE PART(ISPEC,IANG,IFRE) IS
!                      THE NUMBER OF THE PARTITIONING OF THE 
!                      SPECTRUM ISPEC AT DIRECTION IANG AND
!                      FREQUENCY IFRE
      INTEGER PARTINFO(NSPEC,MPART)
!                      INFORMATION ABOUT PARTITIONING:
!                      0=SWELL , 1=WINDSEA, 2=MIXED
      INTEGER :: MAXMPART
!     MAXIMUM NUMBER OF PARTITIONINGS USED

      REAL LOWEST
!                      IF HS IS LOWER THAN LOWEST WAVE SYSTEM IS
!                      COMBINED WITH CLOSEST SYSTEM.
      REAL MEANE(NSPEC,MPART)
!                      MEAN ENERGIES, WHERE MEANE(ISPEC,IPART)
!                      IS THE ENERGY OF SPECTRUM NUMBER ISPEC
!                      PARTITIONING IPART
      REAL MEANANG(NSPEC,MPART)
!                      MEAN DIRECTIONS IN RADIAN, LIKE ABOVE
      REAL MEANFRE(NSPEC,MPART)
!                      MEAN FREQUENCIES IN HZ, LIKE ABOVE
      REAL PKFRE(NSPEC,MPART)
!                      PEAK FREQUENCIES IN HZ, LIKE ABOVE
      REAL US(NSPEC)  ! USTAR
      REAL THW(NSPEC) ! WIND DIRECTION IN DEGREE
      REAL TH0        ! FIRST DIRECTION IN SPECTRUM.
      REAL LEVEL
!                      THE LEVEL FOR THE THRESHHOLD CHECK,
!                      SEE EXPLANATION OF SUBROUTINE THRESHHOLD
      LOGICAL CONTROL
!                      STANDARD OUTPUT CONTROL MESSAGES ON UNIT IU06 
!                      .TRUE. = YES , .FALSE. = NO


!     PARAMETER:
!     ----------

!     VARIABLES:
!     ----------

      INTEGER PEAKANG( NSPEC , NANG*NFRE/2 )
!                      CONTAINS DIRECTIONAL INDICES OF PEAKS.
      INTEGER PEAKFRE( NSPEC , NANG*NFRE/2 )
!                      CONTAINS FREQUENCY INDICES OF PEAKS.
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FX, FY, SPREAD
!            MEAN WAVE NUMBERS, SQUARES OF THE SPREAD
!
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SUMX, SUMXX, SUMY, SUMYY 
!                      WORK SPACE FOR SUMMATION WHILE 
!                      COMPUTING MEANE AND SPREAD.
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SUM0, SUMX0, SUMXX0, SUMY0,
     &                                     SUMYY0 
!                      WORK SPACE ONLY USED IN SUBROUTINE SUMENERGY.
      REAL FR(NFRE) ! FREQUENCY ARRAY.
      REAL COSTH(NANG) , SINTH(NANG) ! COSIN AND SIN TABLE

!-------------------------------------------------------------------

      IF( CONTROL ) WRITE(IU06,*) 'CALL OF SUBROUTINE SWELLSEP'

!     1. PARTITION SPECTRA.
!     ---------------------
      IF( CONTROL ) WRITE(IU06,*) 'CALLING SUBROUTINE PARTITIONING'
      CALL PARTITIONING( 
     &               NSPEC , MPART , NPART , NANG , NFRE ,
     &               SPEC , FLMIN, IU06,
     &               PART , PEAKANG , PEAKFRE, MAXMPART )

!     2. COMBINE PEAKS.
!     -----------------

      ALLOCATE(FX(NSPEC,MPART))
      ALLOCATE(FY(NSPEC,MPART))
      ALLOCATE(SPREAD(NSPEC,MPART))
      ALLOCATE(SUMX(NSPEC,MPART))
      ALLOCATE(SUMXX(NSPEC,MPART))
      ALLOCATE(SUMY(NSPEC,MPART))
      ALLOCATE(SUMYY(NSPEC,MPART))
      ALLOCATE(SUM0(NSPEC,MPART))
      ALLOCATE(SUMX0(NSPEC,MPART))
      ALLOCATE(SUMXX0(NSPEC,MPART))
      ALLOCATE(SUMY0(NSPEC,MPART))
      ALLOCATE(SUMYY0(NSPEC,MPART))

      IF( CONTROL ) WRITE(IU06,*) 'CALLING SUBROUTINE COMBINEPEAKS'
      CALL COMBINEPEAKS(
!                  DIMENSIONS
     &               NSPEC ,MPART , NPART , NANG , NFRE ,
!                  INPUT
     &               SPEC ,
!                  INPUT AND OUTPUT 
     &               PART , PARTINFO , PEAKANG, PEAKFRE, FX, FY,
     &               LOWEST,MEANE , SPREAD , 
     &               MEANANG, MEANFRE , PKFRE ,
!                  WORK SPACE
     &               SUMX , SUMY , 
     &               SUMXX , SUMYY ,
     &               SUM0 , SUMX0 , SUMY0 , SUMXX0 , SUMYY0 ,
!                  INPUT OF PRECOMPUTED TABLES.
     &               US , THW , TH0 , FR , COSTH , SINTH ,
     &               LEVEL , CONTROL, IU06 )

!     3. COMPUTE INTEGRATED VALUES HS,MEAN FREQUENCY,MEAN DIRECTION.
!     --------------------------------------------------------------

      IF( CONTROL ) WRITE(IU06,*) 'CALLING SUBROUTINE MEANS'
      CALL MEANS(
!                   DIMENSIONS
     &                NSPEC , MPART , NPART , NANG , NFRE ,
!                   INPUT
     &                SPEC , PART ,
!                   OUTPUT
     &                MEANE , MEANANG , MEANFRE , PKFRE ,
!                   WORKSPACE
     &                SUM0 , SUMX0 , SUMY0 )

!     4. DEALLOCATE 
!     -------------

      DEALLOCATE(FX)
      DEALLOCATE(FY)
      DEALLOCATE(SPREAD)
      DEALLOCATE(SUMX)
      DEALLOCATE(SUMXX)
      DEALLOCATE(SUMY)
      DEALLOCATE(SUMYY)
      DEALLOCATE(SUM0)
      DEALLOCATE(SUMX0)
      DEALLOCATE(SUMXX0)
      DEALLOCATE(SUMY0)
      DEALLOCATE(SUMYY0)

      IF( CONTROL ) WRITE(IU06,*) 'ENDING SUBROUTINE SWELLSEP'
      IF( CONTROL ) CALL FLUSH(IU06)

      RETURN

      END SUBROUTINE SWELLSEP
