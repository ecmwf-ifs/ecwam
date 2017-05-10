      SUBROUTINE ONESTEPAPART(
!                  DIMENSIONS
     &               NSPEC , MPART , NPART , NANG , NFRE ,
!                  INPUT/OUTPUT
     &               SPEC , PART , PEAKANG , PEAKFRE )

!-----------------------------------------------------------------------

!     PURPOSE:
!     --------

!     COMBINE PEAKS WHICH ARE ONLY ONE INDEX APART

!-----------------------------------------------------------------------

      IMPLICIT NONE

!     INTERFACE:
!     ----------
!     
      INTEGER NSPEC, MPART, NPART(NSPEC), NANG, NFRE
!                  DIMENSIONS.
      REAL SPEC(0:NSPEC,NANG,NFRE)
!                  2D SPECTRA.
      INTEGER PART(NSPEC,NANG,NFRE),
     &        PEAKANG(NSPEC,NANG*NFRE/2),
     &        PEAKFRE(NSPEC,NANG*NFRE/2)
!       PART(..) GIVES NUMBER OF PARTITIONING FOR EACH SPECTRAL BIN.
!       partinfo =1 -> wind sea, =2 -> mixed wind sea-swell =0 -> swell
!        , 3=old windsea.
!       INDEX OF PEAK DIRECTION AND PEAK FREQUENCY OF PARTITIONINGS.

!     LOCAL VARIABLES:
!     ----------------
      INTEGER ISPEC, IPART, JPART, IANG, IFRE
!                  LOOP INDEXES.

!-----------------------------------------------------------------------

      DO ISPEC = 1,NSPEC
        DO IPART = 1,NPART(ISPEC)-1
          DO JPART = IPART+1,NPART(ISPEC)
            DO WHILE( 
     &          ((ABS(PEAKANG(ISPEC,IPART)-PEAKANG(ISPEC,JPART)).EQ.1
     &          .OR.ABS(PEAKANG(ISPEC,IPART)-PEAKANG(ISPEC,JPART))
     &          .EQ.NANG-1)
     &          .AND.PEAKFRE(ISPEC,IPART).EQ.PEAKFRE(ISPEC,JPART)
     &          .OR.PEAKANG(ISPEC,IPART).EQ.PEAKANG(ISPEC,JPART)
     &          .AND.ABS(PEAKFRE(ISPEC,IPART)-PEAKFRE(ISPEC,JPART))
     &          .EQ.1)
     &          .AND. JPART.LE.NPART(ISPEC))
!             CHANGE JPART TO IPART AND NPART(ISPEC) TO JPART
              DO IANG = 1,NANG
                DO IFRE = 1,NFRE
                  IF( PART(ISPEC,IANG,IFRE).EQ.JPART )THEN
                    PART(ISPEC,IANG,IFRE) = IPART
                  ELSE IF( PART(ISPEC,IANG,IFRE).EQ.NPART(ISPEC) )THEN
                    PART(ISPEC,IANG,IFRE) = JPART
                  END IF
                END DO
              END DO
              PEAKANG(ISPEC,JPART) = PEAKANG(ISPEC,NPART(ISPEC))
              PEAKFRE(ISPEC,JPART) = PEAKFRE(ISPEC,NPART(ISPEC))
              NPART(ISPEC) = NPART(ISPEC) - 1
            END DO
          END DO
        END DO
      END DO

      RETURN

      END SUBROUTINE ONESTEPAPART
