      SUBROUTINE FILLGAPS( 
!                  DIMENSIONS
     &                NSPEC , MPART , NANG , NFRE,
!                  INPUT (FOR SPEC ALSO OUTPUT)
     &                SPEC ,
     &                FX , FY , FR , CONTROL , IU06 )

!-------------------------------------------------------------------

!     PURPOSE:
!     --------

!     FILL GAPS (=HOLES) IN THE TRANSFORMED SPECTRA 
!     (BUT SOME ZEROS COULD SURVIVE THIS ROUTINE)

!     METHOD:
!     -------

!     ROUTINE BUILDS FRAMES AROUND GAPS,
!     FRAMES ARE REDUCED TO AVOID PEAKS, 
!     AND THEN AN INTERPOLATION IS MADE FOR FILLING THE GAPS.

!     EXTERNALS:
!     ----------

!     MAKEFRAMES  - MAKES THE FRAME.
!     AVOIDPEAKS  - REDUCES SIZE OF FRAMES TO AVOID PEAKS.
!     GAPINTERPOL - INTERPOLATES TO FILL THE GAPS.

!     AUTHOR.
!     -------
!     S.HASSELMANN, MPI HAMBURG, 1993.
!     J.WASZKEWITZ, MPI HAMBURG, 1993.

!-------------------------------------------------------------------

!     MODULE :
!     --------

      USE YOWICE   , ONLY : FLMIN
      USE YOWGAP   , ONLY : MGAP     ,LANG     ,HANG     ,LFRE     ,
     &            HFRE     ,EQUIGAP 

!-------------------------------------------------------------------

      IMPLICIT NONE

!     INTERFACE:
!     ----------

      INTEGER NSPEC, MPART, NANG, NFRE
!               ARRAY DIMENSIONS.
      REAL SPEC(0:NSPEC,NANG,NFRE)
!               2D SPECTRA.
      REAL FX(NANG,NFRE), 
     &     FY(NANG,NFRE), 
     &     FR(NFRE)
!                POLAR COORDINATES AND FREQUENCIES.
      LOGICAL CONTROL
      INTEGER IU06

!     LOCAL VARIABLES:
!     ----------------

      INTEGER NGAP(NSPEC)
!                       THE NUMBER OF GAPS
      INTEGER NPEAK(NSPEC)
!               NUMBER PF PARTITIONINGS IN A SPECTRUM.
      INTEGER PEAKANG(NSPEC,MPART), 
     &        PEAKFRE(NSPEC,MPART)
!               FREQ-DIRECT. INDEX OF PEAKS
      INTEGER ISPEC, IGAP, IANG, IFRE
!                LOOP INDICES

!-------------------------------------------------------------------
!     
!!!! the default value of MGAP was found to be 2*MPART but I don't
!!!  see the reason. Is there a better values ?? 

      MGAP=2*MPART

      ALLOCATE(LANG(NSPEC,MGAP)) 
      ALLOCATE(HANG(NSPEC,MGAP))
      ALLOCATE(LFRE(NSPEC,MGAP))
      ALLOCATE(HFRE(NSPEC,MGAP))
      ALLOCATE(EQUIGAP(NSPEC,MGAP,MGAP))

!     1.1 BUILD FRAMES AROUND GAPS.
!     ----------------------------

      CALL MAKEFRAMES( 
!          DIMENSIONS
     &             NSPEC , NGAP , NANG , NFRE ,
!          INPUT
     &             SPEC )

      IF( CONTROL )THEN
        DO ISPEC = 1,NSPEC
          WRITE(IU06,*) 'SPECTRA NUMBER ',ISPEC,':'
          WRITE(IU06,*) NGAP(ISPEC),' GAPS before avoidpeaks '
          IF( NGAP(ISPEC).GT.0 ) WRITE(IU06,*) 'CORNERS OF FRAMES:'
          DO IGAP = 1,NGAP(ISPEC)
            WRITE(IU06,'(''  FREQUENCIES = '',E8.2,'' , '',E8.2,
     &                ''  ANGLES 
     &= '',I3,'' , '',I3)')
     &          FR(LFRE(ISPEC,IGAP)),FR(HFRE(ISPEC,IGAP)),
     &          NINT((LANG(ISPEC,IGAP)-1.)/NANG*360.),
     &          NINT((HANG(ISPEC,IGAP)-1.)/NANG*360.)
          END DO
        END DO
      END IF

!     1.2 REDUCE SIZE OF FRAMES IF PEAKS ARE CONTAINED.
!     -------------------------------------------------

      CALL AVOIDPEAKS( 
!          DIMENSIONS
     &             NSPEC, MGAP ,NGAP, MPART, NANG, NFRE,
!          INPUT
     &             SPEC ,
     &             NPEAK , PEAKANG , PEAKFRE ,
!          ARRAYS TO BE ADJUSTED
     &             LANG, HANG , LFRE , HFRE )

      IF( CONTROL )THEN
        DO ISPEC = 1,NSPEC
          WRITE(IU06,*) 'SPECTRA NUMBER ',ISPEC,':'
          WRITE(IU06,*) NGAP(ISPEC),' GAPS'
          IF( NGAP(ISPEC).GT.0 ) WRITE(IU06,*) 'CORNERS OF FRAMES:'
          DO IGAP = 1,NGAP(ISPEC)
            WRITE(IU06,'(''  FREQUENCIES = '',E8.2,'' , '',E8.2,
     &                ''  ANGLES 
     &= '',I3,'' , '',I3)')
     &          FR(LFRE(ISPEC,IGAP)),FR(HFRE(ISPEC,IGAP)),
     &          NINT((LANG(ISPEC,IGAP)-1.)/NANG*360.),
     &          NINT((HANG(ISPEC,IGAP)-1.)/NANG*360.)
          END DO
        END DO
      END IF

!-------------------------------------------------------------------

!     2. INTERPOLATE FRAMES WITH 2D PARABOLA FIT.
!     -------------------------------------------


      CALL GAPINTERPOL( 
!          DIMENSIONS
     &             NSPEC , MGAP , NGAP , NANG , NFRE,
!          INPUT ( SPEC ALSO OUTPUT)
     &             SPEC ,
     &             FX , FY , LANG, HANG , LFRE , HFRE ,
     &             IU06 )

!     3. IMPOSE MINIMUM SPECTRAL VALUE
!     --------------------------------

      DO IANG = 1,NANG
        DO IFRE = 1,NFRE
          DO ISPEC = 1,NSPEC
            SPEC(ISPEC,IANG,IFRE) = MAX(SPEC(ISPEC,IANG,IFRE),FLMIN)
          END DO
        END DO
      END DO


!     4. DEALLOCATE ARRAYS
!     --------------------

      IF(ALLOCATED(LANG)) DEALLOCATE(LANG) 
      IF(ALLOCATED(HANG)) DEALLOCATE(HANG)
      IF(ALLOCATED(LFRE)) DEALLOCATE(LFRE)
      IF(ALLOCATED(HFRE)) DEALLOCATE(HFRE)
      IF(ALLOCATED(EQUIGAP)) DEALLOCATE(EQUIGAP)     

      RETURN

      END SUBROUTINE FILLGAPS
