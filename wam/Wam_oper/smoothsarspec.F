       SUBROUTINE SMOOTHSARSPEC(NSPEC,NANG,NFRE,SPEC) 
 
!---------------------------------------------------------------------

!     PURPOSE
!     -------
!     TO SMOOTH SPECTRA THAT ARE USED IN THE PARTITIONING OF THE SAR
!     ASSIMILATION.

!     AUTHOR
!     ------
!     J. BIDLOT  ECMWF APRIL 2001.

!     EXTERNALS
!     ---------

!**********************************************************************

      USE YOWFRED  , ONLY : DFIM
      USE YOWPCONS , ONLY : EPSMIN

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: NSPEC,NANG,NFRE,ISPEC,K,M,KM,KP

      REAL :: C16
      REAL :: SPEC(0:NSPEC,NANG,NFRE)
      REAL :: WORK(NSPEC,NANG,NFRE)
      REAL, DIMENSION(NSPEC) ::  EM, EW, TEMP

!     COMPUTE THE ENERGY OF THE ORIGINAL SPECTRA

      DO ISPEC = 1,NSPEC 
        EM(ISPEC) = EPSMIN 
      ENDDO

      DO M=1,NFRE
        DO ISPEC = 1,NSPEC
          TEMP(ISPEC) = 0.
        ENDDO
        DO K=1,NANG
          DO ISPEC = 1,NSPEC
            TEMP(ISPEC) = TEMP(ISPEC)+SPEC(ISPEC,K,M)
          ENDDO
        ENDDO
        DO ISPEC = 1,NSPEC 
          EM(ISPEC) = EM(ISPEC)+DFIM(M)*TEMP(ISPEC)
        ENDDO
      ENDDO


!     SMOOTH SPECTRA

      DO K=1,NANG
        KM  = K - 1
        IF (KM .LT. 1) KM = NANG
        KP = K + 1
        IF (KP.GT. NANG) KP = 1

        DO M=2,NFRE-1
          DO ISPEC = 1,NSPEC
            WORK(ISPEC,K,M) = SPEC(ISPEC,KM,M-1) +
     &                     2.*SPEC(ISPEC,KM,M)   + 
     &                        SPEC(ISPEC,KM,M+1) + 
     &                     2.*SPEC(ISPEC,K,M-1)  + 
     &                     4.*SPEC(ISPEC,K,M)    + 
     &                     2.*SPEC(ISPEC,K,M+1)  + 
     &                        SPEC(ISPEC,KP,M-1) + 
     &                     2.*SPEC(ISPEC,KP,M)   + 
     &                        SPEC(ISPEC,KP,M+1)
          ENDDO
        ENDDO
      ENDDO

      C16=1./16

      DO K=1,NANG
        KM  = K - 1
        IF (KM .LT. 1) KM = NANG
        KP = K + 1
        IF (KP.GT. NANG) KP = 1

        DO M=2,NFRE-1
          DO ISPEC = 1,NSPEC
            SPEC(ISPEC,K,M) = C16 * WORK(ISPEC,K,M)
          ENDDO
        ENDDO
      ENDDO

!     COMPUTE THE ENERGY OF THE SMOOTHED SPECTRA

      DO ISPEC = 1,NSPEC 
        EW(ISPEC) = EPSMIN 
      ENDDO

      DO M=1,NFRE
        DO ISPEC = 1,NSPEC
          TEMP(ISPEC) = 0.
        ENDDO
        DO K=1,NANG
          DO ISPEC = 1,NSPEC
            TEMP(ISPEC) = TEMP(ISPEC)+SPEC(ISPEC,K,M)
          ENDDO
        ENDDO
        DO ISPEC = 1,NSPEC 
          EW(ISPEC) = EW(ISPEC)+DFIM(M)*TEMP(ISPEC)
        ENDDO
      ENDDO

!     ENSURE THE ENERGY IS CONSERVED BETWEEN ORIGINAL AND SMOOTHED
!     SPECTRA.

      DO M=1,NFRE
        DO K=1,NANG
          DO ISPEC = 1,NSPEC
            SPEC(ISPEC,K,M) = (EM(ISPEC)/EW(ISPEC))*SPEC(ISPEC,K,M)
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE SMOOTHSARSPEC 
