      SUBROUTINE DEV( 
!             DIMENSIONS
     &          NSPEC , MPART,
     &          NPARTW, NPARTS , 
!             INPUT
     &          NCORRELTAB, XKLOW, XKUP, XKCUTOFF,
     &          XMEANEW , XMEANANGW, XMEANFREW ,
     &          PKFREW, PARTINFOW,
     &          XMEANES , XMEANANGS, XMEANFRES ,
     &          PKFRES, PARTINFOS,
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
      REAL :: PKFREW(NSPEC,MPART),PKFRES(NSPEC,MPART)
!         PEAK FREQUENCIES OF PARTITIONINGS
      INTEGER, DIMENSION (NSPEC,MPART) :: PARTINFOW, PARTINFOS
!                   INFORMATION FOR WINDSEA-SWELL:
!                   WINDSEA (=1), MIXED (=2), OR SWELL (=0)
      REAL :: XKLOW(NSPEC), XKUP(NSPEC)
!     LOWER AND UPPER SAR WAVE NUMBER. 
      REAL :: XKCUTOFF(NSPEC)
!     CUT OFF WAVE LENGHT FROM THE SAR INVERSION.
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
      REAL :: COEF, AKA
      REAL 
     & FWX(MPART,NSPEC), FWY(MPART,NSPEC),
!        MEAN WAVE NUMBERS OF WAM PARTITIONINGS.
     & FSX(MPART,NSPEC), FSY(MPART,NSPEC)
!        MEAN WAVE NUMBERS OF SAR PARTITIONINGS.
      REAL AKAW(NSPEC,MPART), AKAS(NSPEC,MPART)
!        INTERMEDIATE STORAGE.

!-------------------------------------------------------------------

      IF( CONTROL ) WRITE(IU06,*) 'CALL OF SUBROUTINE DEV'

!     1. COMPUTE WAVE NUMBER ARRAYS.
!     ------------------------------

      COEF = ZPI**2/G
      DO JSPEC = 1,NSPEC
        DO JJPARTW = 1,NPARTW(JSPEC)
          AKA=COEF * XMEANFREW(JSPEC,JJPARTW)**2
          FWX(JJPARTW,JSPEC)
     &         = AKA * COS( XMEANANGW(JSPEC,JJPARTW) )
          FWY(JJPARTW,JSPEC)
     &         = AKA * SIN( XMEANANGW(JSPEC,JJPARTW) )

          AKAW(JSPEC,JJPARTW)=COEF * PKFREW(JSPEC,JJPARTW)**2
        END DO
        DO JJPARTS = 1,NPARTS(JSPEC)
          AKA=COEF * XMEANFRES(JSPEC,JJPARTS)**2
          FSX(JJPARTS,JSPEC)
     &         = AKA * COS( XMEANANGS(JSPEC,JJPARTS) )
          FSY(JJPARTS,JSPEC)
     &         = AKA * SIN( XMEANANGS(JSPEC,JJPARTS) )

          AKAS(JSPEC,JJPARTS)=COEF * PKFRES(JSPEC,JJPARTS)**2
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

!!!! take only innovation if swell in both systems
!!!! and peak wave number less than cutoff wave number in SAR systems
!!!! and peak wave number less than upper SAR wave number 
!!!! and peak wave number greater than lower SAR wave number 

            IF( PARTINFOW(JSPEC,JJPARTW).EQ.0          .AND.
     &          PARTINFOS(JSPEC,IDSAR).EQ.0            .AND.
     &          AKAS(JSPEC,IDSAR).LT.XKCUTOFF(JSPEC)   .AND.
     &          AKAS(JSPEC,IDSAR).LT.XKUP(JSPEC)       .AND.
     &          AKAS(JSPEC,IDSAR).GT.XKLOW(JSPEC)   ) THEN 
              ED(JSPEC,IDSAR)=
     &          XMEANES(JSPEC,IDSAR)-XMEANEW(JSPEC,JJPARTW)
              XKD(JSPEC,IDSAR)=
     &          FSX(IDSAR,JSPEC)-FWX(JJPARTW,JSPEC)
              YKD(JSPEC,IDSAR)=
     &          FSY(IDSAR,JSPEC)-FWY(JJPARTW,JSPEC)
            ENDIF
          END IF
        END DO
      END DO

      IF( CONTROL ) WRITE(IU06,*) 'ENDING SUBROUTINE DEV'

      RETURN
      END SUBROUTINE DEV
