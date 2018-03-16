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

      USE YOWGAP   , ONLY : MGAP     ,LANG     ,HANG     ,LFRE     ,
     &            HFRE     ,EQUIGAP 

!-------------------------------------------------------------------

      IMPLICIT NONE
#include "avoidpeaks.intfb.h"
#include "gapinterpol.intfb.h"
#include "makeframes.intfb.h"

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

      INTEGER, ALLOCATABLE :: NGAP(:)
!                       THE NUMBER OF GAPS
      INTEGER ISPEC, IGAP, IANG, IFRE
!                LOOP INDICES
      INTEGER :: NSPECSPLIT, ILOOP, NLOOP, ISPECS, ISPECE 

!-------------------------------------------------------------------

      MGAP=MPART+1

!     INTRODUCE NLOOP TO LIMIT THE MEMORY USED BY FILLGAPS WITHOUT
!     DISTROYING TOO MUCH THE VECTORISATION
      NSPECSPLIT=2000
      IF(NSPEC.LT.INT(1.1*NSPECSPLIT)) NSPECSPLIT=NSPEC
      NLOOP=((NSPEC-1)/NSPECSPLIT)+1

      IF( CONTROL ) THEN
         WRITE(IU06,*) 'NUMBER OF LOOPS IN FILLGAPS = ',NLOOP
         CALL FLUSH(IU06)
      ENDIF

      DO ILOOP = 1,NLOOP

        ISPECS=(ILOOP-1)*NSPECSPLIT+1
        ISPECE=MIN(NSPEC,ILOOP*NSPECSPLIT)
 
        ALLOCATE(NGAP(ISPECS:ISPECE)) 
        ALLOCATE(LANG(ISPECS:ISPECE,MGAP)) 
        ALLOCATE(HANG(ISPECS:ISPECE,MGAP))
        ALLOCATE(LFRE(ISPECS:ISPECE,MGAP))
        ALLOCATE(HFRE(ISPECS:ISPECE,MGAP))
        ALLOCATE(EQUIGAP(ISPECS:ISPECE,MGAP,MGAP))

!       1.1 BUILD FRAMES AROUND GAPS.
!       ----------------------------

        IF( CONTROL ) THEN
           WRITE(IU06,*) 'CALLING SUBROUTINE MAKEFRAMES, ILOOP= ',ILOOP
           CALL FLUSH(IU06)
        ENDIF
        CALL MAKEFRAMES( 
!            DIMENSIONS
     &               NSPEC ,ISPECS ,ISPECE , NGAP , NANG , NFRE ,
!            INPUT
     &               SPEC )

!       1.2 REDUCE SIZE OF FRAMES IF PEAKS ARE CONTAINED.
!       -------------------------------------------------

        IF( CONTROL ) THEN
           WRITE(IU06,*) 'CALLING SUBROUTINE AVOIDPEAKS, ILOOP= ',ILOOP
          CALL FLUSH(IU06)
        ENDIF
        CALL AVOIDPEAKS( 
!            DIMENSIONS
     &               NSPEC,ISPECS ,ISPECE , MGAP ,NGAP, MPART,NANG,NFRE,
!            INPUT
     &               SPEC ,
!            ARRAYS TO BE ADJUSTED
     &               LANG, HANG , LFRE , HFRE )

!-------------------------------------------------------------------

!       2. INTERPOLATE FRAMES WITH 2D PARABOLA FIT.
!       -------------------------------------------


        IF( CONTROL ) THEN
           WRITE(IU06,*) 'CALLING SUBROUTINE GAPINTERPOL, ILOOP= ',ILOOP
           CALL FLUSH(IU06)
        ENDIF
        CALL GAPINTERPOL( 
!            DIMENSIONS
     &               NSPEC , ISPECS ,ISPECE, MGAP , NGAP , NANG , NFRE,
!            INPUT ( SPEC ALSO OUTPUT)
     &               SPEC ,
     &               FX , FY , LANG, HANG , LFRE , HFRE ,
     &               IU06 )

!       3. IMPOSE MINIMUM SPECTRAL VALUE
!       --------------------------------

        DO IANG = 1,NANG
          DO IFRE = 1,NFRE
            DO ISPEC = ISPECS ,ISPECE 
              SPEC(ISPEC,IANG,IFRE) = MAX(SPEC(ISPEC,IANG,IFRE),0.)
            END DO
          END DO
        END DO


!       4. DEALLOCATE ARRAYS
!       --------------------

        IF(ALLOCATED(NGAP)) DEALLOCATE(NGAP) 
        IF(ALLOCATED(LANG)) DEALLOCATE(LANG) 
        IF(ALLOCATED(HANG)) DEALLOCATE(HANG)
        IF(ALLOCATED(LFRE)) DEALLOCATE(LFRE)
        IF(ALLOCATED(HFRE)) DEALLOCATE(HFRE)
        IF(ALLOCATED(EQUIGAP)) DEALLOCATE(EQUIGAP)     

      ENDDO

      RETURN

      END SUBROUTINE FILLGAPS
