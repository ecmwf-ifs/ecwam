      SUBROUTINE PARTITIONING( 
!                  DIMENSIONS
     &               NSPEC , MPART , NPART , NANG , NFRE ,
!                  INPUT
     &               SPEC , FLMIN, IU06, 
!                  OUTPUT
     &               PART , EQUIPART , PEAKANG , PEAKFRE, MAXMPART )

!-----------------------------------------------------------------------

!     PURPOSE:
!     --------

!     PARTITIONS SPECTRA.

!     SUSANNE HASSELMANN, MPI HAMBURG, 1992.
!     JUERGEN WASZKEWITZ,MPI HAMBURG, 1993 (VECTORIZATION)

!     METHOD:
!     -------

!     THE SPECTRAL PARTITIONING IS ACHIEVED BY FIRST RUNNING THROUGH
!     THE ARRAY OF ALL SPECTRAL GRID POINTS, DETERMINING FOR EACH
!     GRID POINT ITS HIGHEST NEIGHBOUR. EACH POINT IS ASSIGNED A PART
!     NUMBER. THE FOLLOWING CASES OCCUR:
!     - THE GRID POINT AND ITS HIGHEST NEIGHBOUR HAS NO PART NUMBER:
!       THEN BOTH ARE ASSIGNED THE SAME NEW PART NUMBER.
!     - ONE OF THEM HAS ALREADY A PART NUMBER, THE OTHER DOES NOT.
!       THEN THE OTHER IS ASSIGNED THE SAME PART NUMBER.
!     - BOTH, THE GRID POINT AND THE HIGHEST NEIGHBOUR, HAVE PART 
!       NUMBERS (WHICH ARE DIFFERENT). THEN IN A SPECIAL LIST
!       PART NUMBERS ARE SAVED, THAT THE PART NUMBER OF THE GRID POINT
!       IS EQUIVALENT TO THE PART NUMBER OF THE HIGHEST NEIGHBOUR. 
!     AFTER ASSIGNING THE PART NUMBERS TO EACH GRID POINT, ALL
!     EQUIVALENT PART NUMBERS ARE CHANGED TO ONE UNIQUE PART NUMBER.

!     TO MAKE THIS ROUTINE AS FAST AS POSSIBLE, THE LOOP OVER ALL
!     SPECTRA IS A VECTORIZED INNER LOOP!

!     
!-----------------------------------------------------------------------

      IMPLICIT NONE

!     INTERFACE:
!     ----------
!     
      INTEGER NSPEC, MPART, NPART(NSPEC), NANG, NFRE
!                  ARRAY DIMENSIONS.
      REAL SPEC(0:NSPEC,NANG,NFRE)
!                  2D SPECTRA. 
      REAL :: FLMIN
!             MINIMUM VALUES FOR SPECTRAL COMPONENTS
      INTEGER :: IU06
!             STANDARD OUTPUT UNIT
      INTEGER PART(NSPEC,NANG,NFRE),
     &        PEAKANG(NSPEC,NANG*NFRE/2),
     &        PEAKFRE(NSPEC,NANG*NFRE/2)
!       PART(..) GIVES NUMBER OF PARTITIONING FOR EACH SPECTRAL BIN.
!       partinfo =1 -> wind sea, =2 -> mixed wind sea-swell =0 -> swell
!        , 3=old windsea.
!       INDEX OF PEAK DIRECTION AND PEAK FREQUENCY OF PARTITIONINGS.
      INTEGER EQUIPART( NSPEC , NANG*NFRE/2 )
!                  CONTAINS EQUIVALENT PART NUMBERS
      INTEGER :: MAXMPART
!     MAXIMUM NUMBER OF PARTITIONINGS USED

!     local variables:
!     ----------------
      INTEGER ISPEC, IPART, JPART , IANG, IFRE, MAXJPART
!                   LOOP INDEXES.
      INTEGER IANGP1, IANGM1
!                  "IANG PLUS 1" , "IANG MINUS 1"
      INTEGER HANG, HFRE
!                   ANGLE AND FREQUENCY INDEX OF HIGHEST NEIGHBOUR
      REAL HSPEC
!                   VALUE OF HIGHEST NEIGHBOUR
      REAL HSPEC0
!                   FOR TEMPORARY USE

!-----------------------------------------------------------------------

!     1.INTITIALIZE.
!     --------------

      DO ISPEC = 1,NSPEC
        NPART(ISPEC) = 0
        DO IANG = 1,NANG
          DO IFRE = 1,NFRE
            PART(ISPEC,IANG,IFRE) = 0
          END DO
        END DO
        DO IPART = 1,NANG*NFRE/2
          EQUIPART(ISPEC,IPART) = 0
          PEAKANG(ISPEC,IPART) = 0
          PEAKFRE(ISPEC,IPART) = 0
        END DO
      END DO

!     2. PARTITION SPECTRA.
!     ---------------------

      DO IANG = 1,NANG
        DO IFRE = 1,NFRE
          DO ISPEC = 1,NSPEC

            IF( SPEC(ISPEC,IANG,IFRE).LE.FLMIN )THEN
!             in WAM the minimum values of spectral components is set
!             to FLMIN 
              PART(ISPEC,IANG,IFRE) = 1

            ELSE
!             FIND HIGHEST NEIGHBOUR

              HANG = IANG
              HFRE = IFRE
              HSPEC = SPEC(ISPEC,HANG,HFRE)
              IF( IANG.NE.NANG )THEN
                IANGP1 = IANG + 1
              ELSE
                IANGP1 = 1
              ENDIF
              HSPEC0 = SPEC(ISPEC,IANGP1,IFRE)
              IF( HSPEC0.GE.HSPEC )THEN
                HANG = IANGP1
                HFRE = IFRE
                HSPEC = HSPEC0
              END IF
              IF( IANG.NE.1 )THEN
                IANGM1 = IANG - 1
              ELSE
                IANGM1 = NANG
              ENDIF
              HSPEC0 = SPEC(ISPEC,IANGM1,IFRE)
              IF( HSPEC0.GT.HSPEC )THEN
                HANG = IANGM1
                HFRE = IFRE
                HSPEC = HSPEC0
              END IF
              IF( IFRE.NE.NFRE )THEN
                HSPEC0 = SPEC(ISPEC,IANG,IFRE+1)
                IF( HSPEC0.GE.HSPEC )THEN
                  HANG = IANG
                  HFRE = IFRE+1
                  HSPEC = HSPEC0
                END IF
              END IF
              IF( IFRE.NE.1 )THEN
                HSPEC0 = SPEC(ISPEC,IANG,IFRE-1)
                IF( HSPEC0.GT.HSPEC )THEN
                  HANG = IANG
                  HFRE = IFRE-1
                END IF
              END IF

!             SET PART NUMBER TO POINT OR NEIGHBPOUR

              IF( IANG.EQ.HANG .AND. IFRE.EQ.HFRE )THEN
!               POINT IS PEAK
                IF( PART(ISPEC,IANG,IFRE).EQ.0 )THEN
                  NPART(ISPEC) = NPART(ISPEC) + 1
                  PART(ISPEC,IANG,IFRE) = NPART(ISPEC)
                END IF
                PEAKANG(ISPEC,PART(ISPEC,IANG,IFRE)) = IANG
                PEAKFRE(ISPEC,PART(ISPEC,IANG,IFRE)) = IFRE
              ELSE
!               POINT IS NO PEAK
                IF( PART(ISPEC,IANG,IFRE).NE.0 )THEN
                  IF( PART(ISPEC,HANG,HFRE).NE.0 )THEN
                   EQUIPART(ISPEC,PART(ISPEC,IANG,IFRE))
     &                = PART(ISPEC,HANG,HFRE)
                  ELSE
                    PART(ISPEC,HANG,HFRE) = PART(ISPEC,IANG,IFRE)
                  END IF
                ELSE
                  IF( PART(ISPEC,HANG,HFRE).NE.0 )THEN
                    PART(ISPEC,IANG,IFRE) = PART(ISPEC,HANG,HFRE)
                  ELSE
                    NPART(ISPEC) = NPART(ISPEC) + 1
                    PART(ISPEC,IANG,IFRE) = NPART(ISPEC)
                    PART(ISPEC,HANG,HFRE) = NPART(ISPEC)
                  ENDIF
                END IF
              END IF
            END IF
          END DO
        END DO
      END DO


!    3. ASSIGN EQUIVALENT PART NUMBERS TO ONE UNIQUE NUMBER.
!     ------------------------------------------------------

!     FIRST: ASSIGN TO EVERY PART NUMBER THIS UNIQUE NUMBER

      DO ISPEC = 1,NSPEC
        JPART = 0
        DO IPART = 1,NPART(ISPEC)
          IF( EQUIPART(ISPEC,IPART).EQ.0 )THEN
            JPART = JPART + 1
            EQUIPART(ISPEC,IPART) = -JPART
            PEAKANG(ISPEC,JPART) = PEAKANG(ISPEC,IPART)
            PEAKFRE(ISPEC,JPART) = PEAKFRE(ISPEC,IPART)
          END IF
        END DO
        NPART(ISPEC) = JPART
        IF(NPART(ISPEC).EQ.0) THEN
          NPART(ISPEC)=1
          PEAKANG(ISPEC,NPART(ISPEC))=1
          PEAKFRE(ISPEC,NPART(ISPEC))=1
        ENDIF

        MAXJPART = MAX(MAXJPART,NPART(ISPEC))

!       NOW NPART(ISPEC) CONTAINS THE NUMBER OF PARTITIONS
!       AND THE LIST EQUIPART( , ) ENDS WITH THE FIRST 0.

        IPART = 1
        DO WHILE( EQUIPART(ISPEC,IPART).NE.0.AND.IPART.LE.NANG*NFRE/2 ) 
          JPART = IPART
          DO WHILE( JPART.GE.0 )
            JPART = EQUIPART(ISPEC,JPART)
          END DO
          EQUIPART(ISPEC,IPART) = JPART
          IPART = IPART + 1
        END DO
      END DO

      MAXJPART=0
      DO ISPEC = 1,NSPEC
        MAXJPART = MAX(MAXJPART,NPART(ISPEC))
      END DO
      MAXMPART = MAX(MAXMPART,MAXJPART)

      WRITE(IU06,*) ' '
      WRITE(IU06,*) '  SUBROUTING PARTITIONING :'
      WRITE(IU06,*) '  THE MAXIMUM NUMBER OF PARTITIONING IS: ',
     &                 MAXJPART
      WRITE(IU06,*) '  FOR A TOTAL NUMBER OF PARTITIONING OF: ',
     &                 MAXMPART
      WRITE(IU06,*) ' '
      CALL FLUSH(IU06)


!!! we need to find a more flexible scheme which will not crashed !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF( MAXJPART.GT.MPART )THEN
        WRITE(IU06,*) '************************************************'
        WRITE(IU06,*) 'ERROR! SUBROUTINE PARTITIONING:'
        WRITE(IU06,*) 'DIMENSION MPART IS TOO SMALL'
        WRITE(IU06,*) 'FOR THE PARTITIONING OF SPECTRA '
        WRITE(IU06,*) 'IT MUST BE MPART >= ',MAXJPART
        WRITE(IU06,*) 'CHANGE PARAMETER MPARTSW'
        WRITE(IU06,*) 'IN MODULE YOWSARAS'
        WRITE(IU06,*) '************************************************'
        CALL ABORT1 
      ENDIF


!     THEN: ASSIGN TO EVERY POINT THIS UNIQUE PART NUMBER

      DO IANG = 1,NANG
        DO IFRE = 1,NFRE
          DO ISPEC =1,NSPEC
            PART(ISPEC,IANG,IFRE) 
     &         = MAX( 1 , -EQUIPART(ISPEC,PART(ISPEC,IANG,IFRE)) )
          END DO
        END DO
      END DO
 
      RETURN

      END SUBROUTINE PARTITIONING
