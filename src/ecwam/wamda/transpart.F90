      SUBROUTINE TRANSPART(                                             &
     &                 IJSG, IJLG, FL1,                                 &
!                    DIMENSIONS
     &                 NSPEC , MPART , NPARTW ,NANG , NFRE ,            &
!                    INPUT/OUTPUT
     &                 SPECW , NW1D,                                    &
!                    INPUT
     &                 PARTW ,CORRELTAW ,MEANWE,MEANWANG,MEANWFRE,      &
     &                 MEANSE , MEANSANG , MEANSFRE ,                   &
     &                 FRATIO,                                          &
!                    OUTPUT
     &                 ADJUST , INTROTATE , FRACROTATE ,                &
     &                 STRETCH , INTLOGSTRETCH )

!-------------------------------------------------------------------

!     PURPOSE:
!     --------

!     TRANSFORMS PARTITIONINGS OF WAM FIRST GUESS SPECTRA TO 
!     MATCH THE MEAN VALUES TO THE PARTITIONINGS OF THE SAR 
!     RETRIEVED SPECTRA

!     S. HASSELMANN, MPI HAMBURG, 1993.
!     J. WASZKEWITZ, MPI HAMBURG, 1993.

!-------------------------------------------------------------------

!     MODULE :
!     --------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : DELTH
      USE YOWGAP   , ONLY : GAPVAL 
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, ICHNKFROMIJ, IPRMFROMIJ

!-------------------------------------------------------------------

      IMPLICIT NONE

!     INTERFACE:
!     ----------

      INTEGER(KIND=JWIM), INTENT(IN) :: IJSG, IJLG
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(IN) :: FL1

      INTEGER(KIND=JWIM) :: NSPEC, MPART, NPARTW(NSPEC), NANG, NFRE
!                   DIMENSIONS.
      REAL(KIND=JWRB) :: SPECW(0:NSPEC,NANG,NFRE) 
!                   WAM SPECTRM BEFORE AND AFTER CORRECTION 

      INTEGER(KIND=JWIM) :: NW1D(IJSG:IJLG)

      INTEGER(KIND=JWIM) :: PARTW(NSPEC,NANG,NFRE)
!                   ASSIGNMENT OF A PARTITIONING NUMBER TO EACH 
!                   DIRECTIONAL-FREQUENCY BIN.
      INTEGER(KIND=JWIM) :: CORRELTAW(NSPEC,MPART)
!                   CROSSASSIGNMENT OF WAM PARTITIONINGS TO SAR
!                   PARTITIONINGS.
      REAL(KIND=JWRB) :: MEANWE(NSPEC,MPART), &
     &     MEANWANG(NSPEC,MPART),             &
     &     MEANWFRE(NSPEC,MPART),             &
     &     MEANSE(NSPEC,MPART),               &
     &     MEANSANG(NSPEC,MPART),             &
     &     MEANSFRE(NSPEC,MPART)
!                   MEAN VALUES FOR WAM AND SAR PARTITIONINGS.
      REAL(KIND=JWRB) :: ADJUST(NSPEC,MPART)
!                   FACTOR FOR ENERGY ADJUSTMENT.
      INTEGER(KIND=JWIM) :: INTROTATE(NSPEC,MPART)
!                   ROTAION PARAMETER(NUMBER OF GRID POINTS)
      REAL(KIND=JWRB) :: FRACROTATE(NSPEC,MPART)
!                   FRACTIONAL ROTATION PARAMETER.
      REAL(KIND=JWRB) :: STRETCH(NSPEC,MPART)
!                   FREQUENCY TRANSFORMATION PARAMETER.
      INTEGER(KIND=JWIM) :: INTLOGSTRETCH(NSPEC,MPART)
!                   FREQUENCY TRANSFORMATION PARAMETER(GRID)

      REAL(KIND=JWRB) :: FRATIO 

!     LOCAL VARIABLES:
!     ----------------

      INTEGER(KIND=JWIM) :: ISPEC, IPARTW, IPARTS, IANG, IFRE, JP
      INTEGER(KIND=JWIM) :: IP, IK
      INTEGER(KIND=JWIM) :: NPROMA, JKGLO, KIJS, KIJL
!                    LOOP INDICES.
      INTEGER(KIND=JWIM) :: JANG, JANGP1, JFRE, JFREP1
!                    GRIDPOINT INDICES AFTER INTERPOLATED
!                    JANGP1 = "JANG PLUS 1" , JFREP1 = "JFRE PLUS 1"
      REAL(KIND=JWRB) :: AANG, BANG, AFRE, BFRE
!                    FOR INTERPOLATION
      REAL(KIND=JWRB) :: ROTATE
!                    THE NUMBER OF GRIDPOINTS THE PARTITIONING 
!                    MUST BE ROTATED
      REAL(KIND=JWRB) :: COPO(NFRE+1), COPODIFF(NFRE+1)
!                   POWERS OF FRATIO AND DIFFERENCES OF POWERS OF 
!                   FRATIO, ONLY USED IN SUBROUTINE TRANSFORM.
      REAL(KIND=JWRB) :: LOGCO, ENERGY, STRETCHFAC

      REAL(KIND=JWRB) :: SPECT(NSPEC,NANG,NFRE+1)
!                   WORK SPACE.

!-------------------------------------------------------------------

!     INITIALIZATION:
!     ---------------

      LOGCO = LOG(FRATIO)
      COPO(1) = FRATIO 
      DO IFRE = 2,NFRE+1
        COPO(IFRE) = COPO(IFRE-1) * FRATIO 
        COPODIFF(IFRE) = COPO(IFRE) - COPO(IFRE-1)
      END DO

      DO ISPEC = 1,NSPEC
        INTROTATE(ISPEC,1)=0
        DO IPARTW = 1,NPARTW(ISPEC)
          IPARTS = CORRELTAW(ISPEC,IPARTW)
          IF ( IPARTS == 0 )THEN
            ADJUST(ISPEC,IPARTW) = 1._JWRB
            INTROTATE(ISPEC,IPARTW) = 0
            FRACROTATE(ISPEC,IPARTW) = 0._JWRB
            STRETCH(ISPEC,IPARTW) = 1._JWRB
            INTLOGSTRETCH(ISPEC,IPARTW) = 0
          ELSE
            ADJUST(ISPEC,IPARTW)                                   &
     &        = (MEANSE(ISPEC,IPARTS) * MEANWFRE(ISPEC,IPARTW)) /  &
     &          (MEANWE(ISPEC,IPARTW) * MEANSFRE(ISPEC,IPARTS))
!!!!! ??????
            ROTATE = ( MEANSANG(ISPEC,IPARTS) - MEANWANG(ISPEC,IPARTW) )
            ROTATE=ROTATE/ DELTH
            INTROTATE(ISPEC,IPARTW) = INT(ROTATE)
            FRACROTATE(ISPEC,IPARTW) = ROTATE - INT(ROTATE)
            STRETCH(ISPEC,IPARTW) = MEANSFRE(ISPEC,IPARTS) / MEANWFRE(ISPEC,IPARTW)
            INTLOGSTRETCH(ISPEC,IPARTW) = INT( LOG(STRETCH(ISPEC,IPARTW))/LOGCO+10000 ) - 10000
          END IF
        END DO
      END DO

!     ROTATE, STRETCH AND ADJUST THE WAM FIRST GUESS PARTITIONINGS
!     IN ONE STEP AND ADD THEM TO ARRAY SPECT:
!     ------------------------------------------------------------

!     restore the unsmoothed model spectra
!!!  will not work if multi block !!!

      NPROMA=NPROMA_WAM
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO, KIJS, KIJL, JP, IK, IP, JANG, JFRE)
      DO JKGLO = IJSG, IJLG, NPROMA
        KIJS=JKGLO
        KIJL=MIN(KIJS+NPROMA-1, IJLG)
        DO JP = KIJS,KIJL
          IF (NW1D(JP) > 0) THEN
            IK = ICHNKFROMIJ(JP)
            IP = IPRMFROMIJ(JP)
            DO JANG = 1,NANG
              DO JFRE = 1,NFRE
                SPECW(NW1D(JP),JANG,JFRE) = FL1(IP, JANG, JFRE, IK)
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO


      DO JANG = 1,NANG
        DO JFRE = 1,NFRE+1
          DO ISPEC = 1,NSPEC
            SPECT(ISPEC,JANG,JFRE) = GAPVAL
          END DO
        END DO
      END DO

      DO IANG = 1,NANG
        DO IFRE = 1,NFRE
          DO ISPEC = 1,NSPEC
            IPARTW = PARTW(ISPEC,IANG,IFRE)           

            JANG = IANG + INTROTATE(ISPEC,IPARTW)
            IF ( JANG > NANG ) JANG = JANG - NANG
            IF ( JANG < 1 ) JANG = JANG + NANG
            IF ((INTROTATE(ISPEC,IPARTW)+FRACROTATE(ISPEC,IPARTW)) > 0) THEN
             JANGP1 = JANG + 1
            ELSE
             JANGP1 = JANG - 1
            END IF
            IF ( JANGP1 > NANG ) JANGP1 = JANGP1 - NANG
            IF ( JANGP1 < 1 ) JANGP1 = JANGP1 + NANG

            JFRE = IFRE + INTLOGSTRETCH(ISPEC,IPARTW)
            IF (JFRE < 1) THEN
!!!!              JFRE=1
!!! it will artificially put energy in the lowest bins otherwise
              STRETCHFAC=1.
            ELSE
              STRETCHFAC=STRETCH(ISPEC,IPARTW)
            ENDIF
            JFREP1 = JFRE + 1

            IF ( JFRE >= 1 .AND. JFRE <= NFRE )THEN
              BANG = ABS(FRACROTATE(ISPEC,IPARTW))
              AANG = 1 - BANG
              AFRE = ( COPO(JFREP1) - COPO(IFRE) * STRETCHFAC ) / COPODIFF(JFREP1)
              BFRE = 1 - AFRE
              ENERGY = SPECW(ISPEC,IANG,IFRE) * ADJUST(ISPEC,IPARTW)
              SPECT(ISPEC,JANG  ,JFRE  ) = SPECT(ISPEC,JANG  ,JFRE  ) + ENERGY * AANG * AFRE 
              SPECT(ISPEC,JANG  ,JFREP1) = SPECT(ISPEC,JANG  ,JFREP1) + ENERGY * AANG * BFRE
              SPECT(ISPEC,JANGP1,JFRE  ) = SPECT(ISPEC,JANGP1,JFRE  ) + ENERGY * BANG * AFRE
              SPECT(ISPEC,JANGP1,JFREP1) = SPECT(ISPEC,JANGP1,JFREP1) + ENERGY * BANG * BFRE 
            END IF 
          END DO
        END DO
      END DO

!     TRANSFER SPECT TO SPECW
!     -----------------------

      DO IANG = 1,NANG
        DO IFRE = 1,NFRE
          DO ISPEC = 1,NSPEC
            SPECW(ISPEC,IANG,IFRE) = SPECT(ISPEC,IANG,IFRE)
          END DO
        END DO
      END DO

!     ADD UNASSIGNED SAR PARTITIONS
!     -----------------------------

!!    no longer in use

      END SUBROUTINE TRANSPART
