      SUBROUTINE COMBINEPEAKS(
!                  DIMENSIONS
     &               NSPEC ,MPART , NPART , NANG , NFRE ,
!                  INPUT
     &               SPEC ,
!                  INPUT AND OUTPUT
     &               PART , PARTINFO , PEAKANG , PEAKFRE , FX, FY,
     &               LOWEST,MEANE , SPREAD ,  
     &               MEANANG, MEANFRE , PKFRE ,
!                  WORK SPACE
     &               SUMX , SUMY ,SUMXX , SUMYY ,
     &               SUM0 , SUMX0 , SUMY0 , SUMXX0 , SUMYY0 ,
!                  INPUT OF GENERAL PARAMETERS.
     &               US , THW , TH0 , FR , COSTH , SINTH ,
     &               LEVEL , CONTROL, IU06 )

!-------------------------------------------------------------------

!     PURPOSE:
!     --------

!     COMBINE PARTITIONINGS IF
!     - THE PEAKS ARE ONLY ONE GRID POINT APART 
!         (SUBROUTINE ONESTEPAPART)
!     - THE PEAK LIES IN LAST TWO FREQUENCY BINS
!         (SUBROUTINE LASTTWOBINS)
!     - THE DISTANCE BETWEEN PEAKS IS LESS THAN HALF THE SPREAD OF THE
!       PARTITIONING (SUBROUTINE HALFSPREAD)
!     - THE MINIMUM ENERGY BETWEEN PEAKS IS HIGHER THAN LEVEL * MINIMUM
!       PEAK ENERGY (SUBROUTINE THRESHHOLD)
!     ALL WINDSEA PEAKS ARE COMBINED TO ONE.

!---------------------------------------------------------------------

      IMPLICIT NONE

!     INTERFACE:
!     ----------

      INTEGER NSPEC, MPART, NPART(NSPEC), NANG, NFRE
!          ARRAY DIMENSIONS.
      REAL SPEC(0:NSPEC,NANG,NFRE)
!          2D INPUT SPECTRA TO BE PARTITIONED.
      INTEGER PART(NSPEC,NANG,NFRE),
     &  PARTINFO(NSPEC,MPART),
     &  PEAKANG(NSPEC,NANG*NFRE/2),
     &  PEAKFRE(NSPEC,NANG*NFRE/2)
!          PART(..) GIVES NUMBER OF PARTITIONING FOR EACH SPECTRAL BIN.
!          partinfo =1 -> wind sea, =2 -> mixed wind sea-swell,
!                   =0 -> swell, , 3=old windsea
!                  index of peak direction and peak frequency of
!                  partitionings.
      REAL LOWEST
!                  WAVE SYSTEM OF WAVE HEIGHT LESS THAN LOWEST IS 
!                  COMBINED WITH NEXT SYSTEM.
      REAL FX(NSPEC,MPART), FY(NSPEC,MPART),
     &  MEANE(NSPEC,MPART),
     &  SPREAD(NSPEC,MPART),
!          MEAN WAVE NUMBERS, MEAN ENERGY AND SPREAD OF PARTITIONINGS.
     &  SUMX(NSPEC,MPART), SUMY(NSPEC,MPART),
     &  SUMXX(NSPEC,MPART), SUMYY(NSPEC,MPART)
!          WORK SPACE.
      REAL SUM0(NSPEC,MPART),
     &  SUMX0(NSPEC,MPART),SUMXX0(NSPEC,MPART),
     &  SUMY0(NSPEC,MPART),SUMYY0(NSPEC,MPART)
!          FOR TEMPORARY STORAGE.
      REAL US(NSPEC), THW(NSPEC),
     &  FR(NFRE), COSTH(NANG), SINTH(NANG)
!          USTAR, WIND DIRECTION,FREQUENCIES,COS AND SIN TABLES 
!          FREQU./DIRECTIONAL INCREMENT.
      REAL TH0
!                      FIRST DIRECTION IN SPECTRUM.
      REAL LEVEL
!          combine systems if minimum energy between peaks
!          is higher than
!          LEVEL * MINIMUM PEAK ENERGY (SUBROUTINE THRESHHOLD)
       
      REAL MEANANG(NSPEC,MPART)
!                      MEAN DIRECTIONS IN RADIANCE, LIKE ABOVE
      REAL MEANFRE(NSPEC,MPART)
!                      MEAN FREQUENCIES IN HZ, LIKE ABOVE
      REAL PKFRE(NSPEC,MPART)
!                      PEAK FREQUENCIES IN HZ, LIKE ABOVE

      INTEGER ISPEC,IPART
      LOGICAL CONTROL
!          CONTROLS PRINT OUTPUT.
      INTEGER IU06
!         PRINT OUTPUT UNIT.
      INTEGER :: MAXNPARTRE

!-----------------------------------------------------------------------

!     1. COMBINE PARTITIONINGS IF PEAKS ARE ONLY ONE GRID POINT APART 
!     -----------------------------------------------------------------

      IF( CONTROL ) WRITE(IU06,*) '  CALLING SUBROUTINE ONESTEPAPART'
      CALL ONESTEPAPART( 
!               DIMENSIONS
     &            NSPEC , MPART , NPART , NANG , NFRE ,
!               INPUT/OUTPUT
     &            SPEC , PART , PEAKANG , PEAKFRE )


!     2. INITIALIZE ARRAYS FX AND FY
!     -----------------------------------

      DO ISPEC = 1,NSPEC
        DO IPART = 1,NPART(ISPEC)
          FX(ISPEC,IPART) = FR(PEAKFRE(ISPEC,IPART))
     &        * COSTH(PEAKANG(ISPEC,IPART))    
          FY(ISPEC,IPART) = FR(PEAKFRE(ISPEC,IPART))
     &        * SINTH(PEAKANG(ISPEC,IPART))    
        END DO
      END DO

!     3. COMPUTE MEAN ENERGY OF PARTITIONINGS.
!     ----------------------------------------

      IF( CONTROL ) WRITE(IU06,*) '  CALLING SUBROUTINE SUMENERGY'
      CALL SUMENERGY( 
!               DIMENSIONS
     &            NSPEC , MPART , NPART , NANG , NFRE ,
!               INPUT
     &            SPEC , PART ,
!               OUTPUT
     &            MEANE , SPREAD , SUMX , SUMY , SUMXX , SUMYY ,
!               WORK SPACE
     &            SUM0 , SUMX0 , SUMY0 , SUMXX0 , SUMYY0 )

      IF( CONTROL ) WRITE(IU06,*) '  CALLING SUBROUTINE MEANS'
      CALL MEANS(
!               DIMENSIONS
     &            NSPEC , MPART , NPART , NANG , NFRE ,
!                   INPUT
     &            SPEC , PART ,
!               OUTPUT
     &            MEANE , MEANANG, MEANFRE , PKFRE ,
!               WORKSPACE
     &            SUM0 , SUMX0 , SUMY0 )

!     4. ALL WINDSEA PEAKS ARE COMBINED TO ONE.
!     -----------------------------------------

      IF( CONTROL ) WRITE(IU06,*) '  CALLING SUBROUTINE WINDSEA'
      CALL WINDSEA( 
!               DIMENSIONS
     &            NSPEC , MPART , NPART , NANG , NFRE ,
!               INPUT/OUTPUT
     &            SPEC , PART ,
     &            PARTINFO , PEAKANG , PEAKFRE , FX, FY,
     &            MEANE ,  MEANANG, MEANFRE ,
     &            SPREAD ,
     &            SUMX , SUMY , SUMXX , SUMYY ,
     &            US , THW , TH0 , CONTROL, IU06 )

!     5. COMBINE IF THE PEAK LIES IN LAST TWO FREQUENCY BINS.
!     -------------------------------------------------------

      IF( CONTROL ) WRITE(IU06,*) '  CALLING SUBROUTINE LASTTWOBINS'
      CALL LASTTWOBINS( 
!                DIMENSIONS
     &             NSPEC , MPART , NPART , NANG , NFRE ,
!                INPUT/OUTPUT
     &             SPEC , PART , PARTINFO , PEAKANG , PEAKFRE , FX, FY,
     &             MEANE , SPREAD , SUMX , SUMY , SUMXX , SUMYY ,
     &             FR , CONTROL, IU06 )

!     6. COMBINE IF DISTANCE BETWEEN PEAKS IS LESS THAN HALF 
!         THE SPREAD OF THE PARTITIONINGS. 
!     ---------------------------------------------------------

      IF( CONTROL ) WRITE(IU06,*) '  CALLING SUBROUTINE HALFSPREAD'
      CALL HALFSPREAD( 
!                 DIMENSIONS
     &              NSPEC , MPART , NPART , NANG , NFRE ,
!                 INPUT/OUTPUT
     &              SPEC , PART , PARTINFO , PEAKANG , PEAKFRE , FX, FY,
     &              MEANE , SPREAD , 
!                 WORK SPACE
     &             SUMX , SUMY , SUMXX , SUMYY ,
!                 PRECOMPUTED TABLES FOR FREQUENCIES.
     &              FR , CONTROL, IU06 )

!     7. COMBINE IF MINIMUM ENERGY BETWEEN PEAKS IS HIGHER
!        THAN LEVEL * MINIMUM PEAK ENERGY.
!     -----------------------------------------------------

      IF( CONTROL ) WRITE(IU06,*) '  CALLING SUBROUTINE THRESHHOLD'
      CALL THRESHHOLD( 
!                 DIMENSIONS
     &              NSPEC , MPART , NPART , NANG , NFRE ,
!                 INPUT/OUTPUT
     &              SPEC , PART , PARTINFO , PEAKANG , PEAKFRE , FX, FY,
     &              MEANE , SPREAD ,
     &              SUMX , SUMY , SUMXX , SUMYY , 
!                 PRECOMPUTED TABLES FOR FREQUENCIES,THRESHHOLD LEVEL..
     &              FR , LEVEL , CONTROL, IU06 )

!     6. LOW ENERGY PARTITIONINGS ARE COMBINED WITH CLOSEST PARTITIONING
!     ------------------------------------------------------------------

      IF( CONTROL ) WRITE(IU06,*) '  CALLING SUBROUTINE LOWENERGY'
      CALL LOWENERGY( 
!                 DIMENSIONS.
     &              NSPEC , MPART , NPART , NANG , NFRE ,
!                 INPUT/OUTPUT
     &              SPEC , PART , PARTINFO , PEAKANG , PEAKFRE , FX, FY,
     &              LOWEST, MEANE, SPREAD ,
!                 WORK SPACE
     &              SUMX , SUMY , SUMXX , SUMYY,
!                 PRECOMPUTED TABLES FOR FREQUENCIES.
     &              FR , CONTROL, IU06 )

!!! temporary diagnostic
      MAXNPARTRE=0
      DO ISPEC = 1,NSPEC
        MAXNPARTRE=MAX(MAXNPARTRE,NPART(ISPEC))
      ENDDO
      WRITE(IU06,*) '  AT THE END OF SUBROUTINE COMBINEPEAKS'
      WRITE(IU06,*) '  THE MAXIMUN NUMBER OF PARTITIONS WAS REDUCED TO'
      WRITE(IU06,*) '  MAXNPARTRE = ',MAXNPARTRE

      RETURN

      END SUBROUTINE COMBINEPEAKS
