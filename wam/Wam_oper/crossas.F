      SUBROUTINE CROSSAS(
!             DIMENSIONS
     &          NSPEC , MPART , NPARTW, NPARTS ,
!             INPUT
     &          DISMIN,
     &          MEANWANG, MEANWFRE ,
     &          MEANSANG, MEANSFRE ,
     &          CONTROL,
!             OUTPUT
     &          CORRELTAB)

!-----------------------------------------------------------------------

!     PURPOSE:
!     --------

!     USED IN DATA ASSIMILATION SCHEME.

!     CROSSASSIGN WAM AND SAR PARTITIONINGS.

!     TO MAKE THIS ROUTINE AS FAST AS POSSIBLE, THE LOOP OVER ALL
!     SPECTRA IS A VECTORIZED INNER LOOP.
     
!     AUTHOR:
!     -------

!     SUSANNE HASSELMANN 1993 AT MPI FUER METEOROLOGIE HAMBURG
!     JUERGEN WASZKEWTIZ 1993 AT MPI FUER METEOROLOGIE HAMBURG

!     PARAMETERS:
!     -----------

!     EXTERNALS:
!     ----------
!     NONE.

!-----------------------------------------------------------------------

!     MODULES:
!     --------

      USE YOWPCONS , ONLY : G        ,ZPI
      USE YOWTEST  , ONLY : IU06

!     INTERFACE:
!     ----------

      IMPLICIT NONE

      INTEGER NSPEC
!                     NUMBER OF SPECTRA.
      INTEGER MPART
!                     MAXIMAL NUMBER (=DIMENSION SIZE)
!                     OF PARTITIONINGS
      INTEGER NPARTW(NSPEC)
!                     NUMBER OF PARTITIONINGS
!                     OF THE WAM FIRST GUESS SPECTRA
      INTEGER NPARTS(NSPEC)
!                     NUMBER OF 
!                     PARTITIONINGS OF THE INVERTED SAR SPECTRA
      REAL DISMIN
!                     SAR AND WAM WAVE SYSTEMS ARE CROSSASSIGNED IF 
!                     THEIR DISTANCE IS LESS THAN DISMIN.
      INTEGER CORRELTAB(NSPEC,MPART)
!                     CORRELATION TABLE, WHERE CORRELTAB(ISPEC,IPARTW)
!                     MEANS, THAT PARTITIONING IPARTW OF WAM FIRST
!                     GUESS SPECTRUM ISPEC IS CORRELATED WITH 
!                     PARTITIONIG CORRELTAB(ISPEC,IPARTW)
!                     OF INVERTED SAR SPECTRUM ISPEC
!                     (A ZERO MEANS, NO CROSSASSIGNMENT)
      INTEGER FLAG(MPART)
!                     >0 = CORRELTAB
!                     <0 NO CROSSASSIGNMENT TO A WAM-PARTITIONING.
      REAL MEANWANG(NSPEC,MPART)
!                     MEAN DIRECTIONS IN RADIANCE
!                     OF PARTITIONINGS OF WAM FIRST GUESS SPECTRA.
      REAL MEANWFRE(NSPEC,MPART)
!                     MEAN FREQUENCIES OF PARTITIONINGS
!                     OF WAM FIRST GUESS SPECTRA.
      REAL MEANSANG(NSPEC,MPART)
!                     MEAN DIRECTIONS IN RADIANCE
!                     OF PARTITIONINGS OF INVERTED SAR SPECTRA.
      REAL MEANSFRE(NSPEC,MPART)
!                     MEAN FREQUENCIES OF PARTITIONINGS
!                     OF INVERTED SAR SPECTRA.
      LOGICAL CONTROL
!                     (INPUT) WRITE CONTROL MESSAGES ON UNIT 6 ?
!                     .TRUE. = YES , .FALSE. = NO

!     VARIABLES:
!     ----------

      INTEGER ISPEC, IPARTW, IPARTS
!                      LOOP INDEXES.
      INTEGER NPARTWMAX, NPARTSMAX
!                      MAXIMUM VALUE OF NPARTW/NPARTS OVER ALL 
!                      SPECTRA
      REAL, ALLOCATABLE :: DIST(:,:,:)
!                      THE SQUARE OF THE DISTANCE BETWEEN TWO SYSTEMS
!                      IN K-SPACE.
      REAL MINDIST
!                      MINIMAL DISTANCE FOR CROSS ASSIGNMENT.
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FWX, FWY
!                      MEAN WAVE NUMBERS OF WAM PARTITIONINGS.
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FSX, FSY
!                      MEAN WAVE NUMBERS OF SAR PARTITIONINGS.
      REAL AKABS
      INTEGER IPARWH,IPARWH1,IPARSH,IPSH
      REAL MINDIS(MPART),MINDISH

!---------------------------------------------------------------------

      IF( CONTROL ) WRITE(IU06,*) 'CALL OF SUBROUTINE CORREL'

      NPARTWMAX = 0
      NPARTSMAX = 0
      DO ISPEC = 1,NSPEC
        NPARTWMAX = MAX(NPARTWMAX,NPARTW(ISPEC))
        NPARTSMAX = MAX(NPARTSMAX,NPARTS(ISPEC))
      END DO
      ALLOCATE(FWX(NSPEC,NPARTWMAX))
      ALLOCATE(FWY(NSPEC,NPARTWMAX))
      ALLOCATE(FSX(NSPEC,NPARTSMAX))
      ALLOCATE(FSY(NSPEC,NPARTSMAX))
      ALLOCATE(DIST(NSPEC,NPARTWMAX,NPARTSMAX))

!      COMPUTE WAVE NUMBER ARRAYS.

      DO ISPEC = 1,NSPEC
        DO IPARTW=1,NPARTWMAX
          FWX(ISPEC,IPARTW)=0.1E-10
          FWY(ISPEC,IPARTW)=0.1E-10
        END DO
      END DO
      DO ISPEC = 1,NSPEC
        DO IPARTS=1,NPARTSMAX
          FSX(ISPEC,IPARTs)=0.1E-10
          FSY(ISPEC,IPARTs)=0.1E-10
        END DO
      END DO
      DO ISPEC = 1,NSPEC
        DO IPARTW = 1,NPARTW(ISPEC)
          AKABS=(ZPI*MEANWFRE(ISPEC,IPARTW))**2/G
          FWX(ISPEC,IPARTW)
     &         = AKABS * COS( MEANWANG(ISPEC,IPARTW) )
          FWY(ISPEC,IPARTW)
     &         = AKABS * SIN( MEANWANG(ISPEC,IPARTW) )
        END DO
      END DO
      DO ISPEC = 1,NSPEC
        DO IPARTS = 1,NPARTS(ISPEC)
          AKABS=(ZPI*MEANSFRE(ISPEC,IPARTS))**2/G
          FSX(ISPEC,IPARTS)
     &         = AKABS * COS( MEANSANG(ISPEC,IPARTS) )
          FSY(ISPEC,IPARTS)
     &         = AKABS * SIN( MEANSANG(ISPEC,IPARTS) )
        END DO
      END DO

!     COMPUTE DISTANCE BETWEEN SAR AND WAM PEAKS.

      DO IPARTS = 1,NPARTSMAX
        DO IPARTW = 1,NPARTWMAX
          DO ISPEC = 1,NSPEC
            AKABS=(0.5*(FWX(ISPEC,IPARTW)**2+FSX(ISPEC,IPARTS)**2+
     &           FWY(ISPEC,IPARTW)**2+FSY(ISPEC,IPARTS)**2))
            AKABS=MAX(AKABS,0.1E-10)
            DIST(ISPEC,IPARTW,IPARTS)
     &          = (( FWX(ISPEC,IPARTW) - FSX(ISPEC,IPARTS) ) ** 2
     &          + ( FWY(ISPEC,IPARTW) - FSY(ISPEC,IPARTS) ) ** 2)/
     &         AKABS
          END DO
        END DO
      END DO

!     COMPUTE MINIMUM DISTANCE FOR EACH WAM TO SAR PEAKS AND 
!     CROSS ASSIGN WAM AND SAR WAVE SYSTEMS.

      DO ISPEC = 1,NSPEC
        DO IPARTS=1,NPARTS(ISPEC)
         FLAG (IPARTS)=-1
         MINDIS(IPARTS)=100.
        END DO
        DO IPARTW = 1,NPARTW(ISPEC)
          CORRELTAB(ISPEC,IPARTW)=0
          MINDIST=100.
          DO IPARTS = 1,NPARTS(ISPEC)
            IF( DIST(ISPEC,IPARTW,IPARTS).LT.DISMIN )THEN
               IF(DIST(ISPEC,IPARTW,IPARTS).LT.MINDIST) THEN
                IF( FLAG (IPARTS).LT.0 ) THEN
                   IF(CORRELTAB(ISPEC,IPARTW).GT.0) THEN
                    FLAG (CORRELTAB(ISPEC,IPARTW))=-1
                    MINDIS(CORRELTAB(ISPEC,IPARTW))=100.
                   END IF
                   MINDIST=DIST(ISPEC,IPARTW,IPARTS)
                   CORRELTAB(ISPEC,IPARTW) = IPARTS
                   FLAG(IPARTS)=ipartw
                   MINDIS(IPARTS)=MINDIST
                 ELSEIF(DIST(ISPEC,IPARTW,IPARTS).LT.MINDIS(IPARTS))
     &                         THEN
                   MINDIST=DIST(ISPEC,IPARTW,IPARTS)
                   CORRELTAB(ISPEC,IPARTW) = IPARTS

                   IF ( FLAG(IPARTS) .NE. IPARTW .AND.
     &                 FLAG (IPARTS) .GT. 0) THEN
                    CORRELTAB(ISPEC, FLAG(IPARTS)) = 0
                    MINDISH=DISMIN
                    IPARWH= FLAG(IPARTS)
                    IPSH=0
                    IPARWH1=-1

                    DO IPARSH=1,NPARTS(ISPEC)
                     IF ( FLAG(IPARSH) .LT. 0 .AND.
     &                   DIST(ISPEC,IPARWH,IPARSH) .LT.  MINDISH .AND.
     &                   IPARSH .NE. IPARTS) THEN
                      MINDISH=DIST(ISPEC,IPARWH,IPARSH)
                      IPSH=IPARSH
                      IPARWH1=IPARWH
                     END IF
                    END DO
                    IF(IPSH.GT.0) THEN
                     CORRELTAB(ISPEC,IPARWH1) = IPSH
                     FLAG(IPSH)=IPARWH1
                     MINDIS(IPSH)=MINDISH
                    END IF
                   END IF
                   FLAG(IPARTS)=IPARTW
!                 ELSE
!                   CORRELTAB(ISPEC,IPARTW) = 0
                 END IF
               END IF
            ELSEIF(MINDIST.EQ.100.) THEN
                CORRELTAB(ISPEC,IPARTW) = 0
            END IF
           END DO
        END DO
       END DO

      DEALLOCATE(FWX)
      DEALLOCATE(FWY)
      DEALLOCATE(FSX)
      DEALLOCATE(FSY)
      DEALLOCATE(DIST)

      IF( CONTROL ) WRITE(IU06,*) 'ENDING SUBROUTINE CROSSAS'

      RETURN

      END SUBROUTINE CROSSAS
