      SUBROUTINE DEV( 
!             DIMENSIONS
     &          NSPEC , MPART,
     &          NPARTW, NPARTS , 
!             INPUT
     &          NCORRELTAB,
     &          XMEANEW , XMEANANGW, XMEANFREW ,
     &          XMEANES , XMEANANGS, XMEANFRES ,
     &          CONTROL,
!             OUTPUT
     &          ED,XKD,YKD) 

!-------------------------------------------------------------------

!     PURPOSE:
!     --------
!     COMPUTE DIFFERENCES BETWEEN WAM AND SAR MEAN ENERGY AND 
!     MEAN WAVE NUMBERS FOR OPTIMAL INTERPOLATION.
!     SAR-WAM
!     
!     AUTHOR:
!     -------

!     S.HASSELMANN MPI FUER METEOROLOGIE, HAMBURG, GERMANY.
!     PIERO LIONELLO (UNIV. OF PADUA) HAMBURG 11/93

!     PARAMETERS:
!     -----------
!     EXTERNALS:
!     ----------
!**********************************************************************

!     MODULES:
!     --------

      USE YOWPCONS , ONLY : G        ,ZPI
      USE YOWTEST  , ONLY : IU06

!     INTERFACE:
!     ----------

      IMPLICIT NONE

      INTEGER
     & NSPEC,
!        NUMBER OF SAR DATA POINTS.
     & MPART,
!         MAXIMAL NUMBER OF PARTITIONINGS
     & NPARTW(NSPEC),NPARTS(NSPEC)
!         NUMBER OF PARTITIONINGSIN WAM OR SAR SPECTRA
      REAL 
     & XMEANEW(NSPEC,MPART),XMEANES(NSPEC,MPART),
!         MEAN ENERGIES OF PARTITIONINGS
     & XMEANANGW(NSPEC,MPART),XMEANANGS(NSPEC,MPART),
!         MEAN DIRECTIONS 
     & XMEANFREW(NSPEC,MPART),XMEANFRES(NSPEC,MPART)
!         MEAN FREQUENCIES OF PARTITIONINGS
      INTEGER 
     & NCORRELTAB(NSPEC,MPART)
!         INDEX OF PARTITIONING MPART IN WAM SPECTRUM NSPEC
!         CROSSASSIGNED TO PARTITIONING NCORRELTAB OF SAR 
!         SPECTRUM NSPEC.

      LOGICAL CONTROL
!         (INPUT) WRITE CONTROL MESSAGES ON UNIT 6 ?
!         .TRUE. = YES , .FALSE. = NO
      REAL   
     & ED(NSPEC,0:MPART),XKD(NSPEC,0:MPART),
     & YKD(NSPEC,0:MPART)
!        DIFFERENCES IN ENERGY AND MEAN WAVE NUMBERS KX,KY 
!        OF WAM AND SAR WAVE SYSTEMS AT DATA POINTS FOR SAR 
!        PSRTITIONING INDEX MPART.

!     VARIABLES:
!     ----------

      INTEGER JJPARTW, JJPARTS, JSPEC, IDSAR
!             COUNTERS
      REAL 
     & FWX(MPART,NSPEC), FWY(MPART,NSPEC),
!        MEAN WAVE NUMBERS OF WAM PARTITIONINGS.
     & FSX(MPART,NSPEC), FSY(MPART,NSPEC)
!        MEAN WAVE NUMBERS OF SAR PARTITIONINGS.
      REAL AKABS
!        INTERMEDIATE STORAGE.

!-------------------------------------------------------------------

      IF( CONTROL ) WRITE(IU06,*) 'CALL OF SUBROUTINE DEV'

!     1. COMPUTE WAVE NUMBER ARRAYS.
!     ------------------------------

      DO JSPEC = 1,NSPEC
        DO JJPARTW = 1,NPARTW(JSPEC)
          AKABS=(ZPI*XMEANFREW(JSPEC,JJPARTW))**2/G
          FWX(JJPARTW,JSPEC)
     &         = AKABS * COS( XMEANANGW(JSPEC,JJPARTW) )
          FWY(JJPARTW,JSPEC)
     &         = AKABS * SIN( XMEANANGW(JSPEC,JJPARTW) )
        END DO
        DO JJPARTS = 1,NPARTS(JSPEC)
          AKABS=(ZPI*XMEANFRES(JSPEC,JJPARTS))**2/G
          FSX(JJPARTS,JSPEC)
     &         = AKABS * COS( XMEANANGS(JSPEC,JJPARTS) )
          FSY(JJPARTS,JSPEC)
     &         = AKABS * SIN( XMEANANGS(JSPEC,JJPARTS) )
        END DO
      END DO

!-------------------------------------------------------------------

!     2. COMPUTES DEVIATIONS.
!     -----------------------
      DO JSPEC=1,NSPEC 
        DO JJPARTS = 0,MPART
          ED (JSPEC,JJPARTS) = 0.
          XKD(JSPEC,JJPARTS) = 0.
          YKD(JSPEC,JJPARTS) = 0.
        END DO
      END DO

      DO JSPEC=1,NSPEC 
        DO JJPARTW = 1,NPARTW(JSPEC)
          IDSAR=NCORRELTAB(JSPEC,JJPARTW)
          IF (IDSAR.NE.0) THEN
            ED(JSPEC,IDSAR)=
     &       XMEANES(JSPEC,IDSAR)-XMEANEW(JSPEC,JJPARTW)
            XKD(JSPEC,IDSAR)=
     &       FSX(IDSAR,JSPEC)-FWX(JJPARTW,JSPEC)
            YKD(JSPEC,IDSAR)=
     &       FSY(IDSAR,JSPEC)-FWY(JJPARTW,JSPEC)
          END IF
        END DO
      END DO

      IF( CONTROL ) WRITE(IU06,*) 'ENDING SUBROUTINE DEV'

      RETURN
      END SUBROUTINE DEV
