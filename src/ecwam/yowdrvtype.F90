! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWDRVTYPE

!     DERIVED TYPES DEFINITION:

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU, JWRO

      IMPLICIT NONE

!*    **  *VARIABLES DEPENDENT ON GLOBAL GRID POINTS ON STRUCTURED GRID
!!! (ideally this be removed, but currently still need for the propagation halo, and the output of global fields
!!!  but also for the wave DA. It should be possible to localise this...)

      TYPE WVGRIDGLO
        INTEGER(KIND=JWIM), DIMENSION(:), ALLOCATABLE :: IXLG ! WEST-EAST INDEX FOR A GIVEN IJ
        INTEGER(KIND=JWIM), DIMENSION(:), ALLOCATABLE :: KXLT ! NORTH-SOUTH INDEX FOR A GIVEN IJ
        CONTAINS
           PROCEDURE :: ALLOC=>WVGRIDGLO_ALLOC
           PROCEDURE :: DEALLOC=>WVGRIDGLO_DEALLOC
      END TYPE WVGRIDGLO 


!*    **  *VARIABLES DEPENDENT ON LOCAL GRID POINTS ON STRUCTURED GRID

      TYPE WVGRIDLOC
        INTEGER(KIND=JWIM), DIMENSION(:,:), ALLOCATABLE :: IFROMIJ ! WEST-EAST INDEX FOR A GIVEN IJ (LOCAL VERSION OF IXLG)
        INTEGER(KIND=JWIM), DIMENSION(:,:), ALLOCATABLE :: KFROMIJ ! SOUTH-NORTH INDEX FOR A GIVEN IJ (LOCAL VERSION OF KXLT)
        INTEGER(KIND=JWIM), DIMENSION(:,:), ALLOCATABLE :: JFROMIJ ! NORTH-SOUTH INDEX FOR A GIVEN IJ (LOCAL VERSION OF NGY-KXLT+1)
        CONTAINS
           PROCEDURE :: ALLOC=>WVGRIDLOC_ALLOC
           PROCEDURE :: DEALLOC=>WVGRIDLOC_DEALLOC
      END TYPE WVGRIDLOC 


!*    **  *VARIABLES DEPENDENT ON EARTH GEOMETRY, DEPTH or CURRENTS (but not frequency or direction
!          i.e. the environment in which the waves evolve 

      TYPE ENVIRONMENT
        INTEGER(KIND=JWIM), DIMENSION(:,:), ALLOCATABLE :: INDEP   ! DEPTH INDEX FOR A BLOCK.
        INTEGER(KIND=JWIM), DIMENSION(:,:), ALLOCATABLE :: IODP    ! 0 OVER LAND, 1 OVER SEA.
        INTEGER(KIND=JWIM), DIMENSION(:,:), ALLOCATABLE :: IOBND   ! 0 OVER BOUNDARY POINTS, 1 OTHERWISE.
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: DELLAM1    ! 1./DELLAM AT BLOCK POINTS.
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: COSPHM1    ! 1./COSPH AT BLOCK POINTS.
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: DEPTH      ! WATER DEPTH IN METRES.
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: EMAXDPT    ! MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH
                                                                   ! EMAXDPT=0.0625*(GAM_B_J*DEPTH)**2
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: UCUR       ! U-COMPONENT OF SURFACE CURRENT (m/s)
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: VCUR       ! V-COMPONENT OF SURFACE CURRENT (m/s)
        CONTAINS
           PROCEDURE :: ALLOC=>ENVIRONMENT_ALLOC
           PROCEDURE :: DEALLOC=>ENVIRONMENT_DEALLOC
      END TYPE ENVIRONMENT


!*    **  *VARIABLES DEPENDENT ON FREQUENCY (but not direction)
!          i.e. wave properties dependent on water depth or sea ice

      TYPE FREQUENCY
        REAL(KIND=JWRB), DIMENSION(:,:,:), ALLOCATABLE :: WAVNUM     ! WAVE NUMBER
        REAL(KIND=JWRB), DIMENSION(:,:,:), ALLOCATABLE :: CINV       ! RECIPROCAL OF THE PHASE VELOCITY (1/c)
        REAL(KIND=JWRB), DIMENSION(:,:,:), ALLOCATABLE :: CGROUP     ! GROUP SPEED
        REAL(KIND=JWRB), DIMENSION(:,:,:), ALLOCATABLE :: XK2CG      ! (WAVE NUMBER)**2 * GROUP SPEED
        REAL(KIND=JWRB), DIMENSION(:,:,:), ALLOCATABLE :: OMOSNH2KD  ! OMEGA / SINH(2KD)
        REAL(KIND=JWRB), DIMENSION(:,:,:), ALLOCATABLE :: STOKFAC    ! FACTOR TO COMPUTE SURFACE STOKES DRIFT FROM SPECTRUM 2*G*K**2/(OMEGA*TANH(2KD))
        REAL(KIND=JWRB), DIMENSION(:,:,:), ALLOCATABLE :: CIWA       ! SEA ICE WAVE ATTENUATION
        CONTAINS
           PROCEDURE :: ALLOC=>FREQUENCY_ALLOC
           PROCEDURE :: DEALLOC=>FREQUENCY_DEALLOC
      END TYPE FREQUENCY

      TYPE FREQUENCY_LAND
        REAL(KIND=JWRB), DIMENSION(:), ALLOCATABLE :: WAVNUM     ! WAVE NUMBER
        REAL(KIND=JWRB), DIMENSION(:), ALLOCATABLE :: CINV       ! RECIPROCAL OF THE PHASE VELOCITY (1/c)
        REAL(KIND=JWRB), DIMENSION(:), ALLOCATABLE :: CGROUP     ! GROUP SPEED
        REAL(KIND=JWRB), DIMENSION(:), ALLOCATABLE :: XK2CG      ! (WAVE NUMBER)**2 * GROUP SPEED
        REAL(KIND=JWRB), DIMENSION(:), ALLOCATABLE :: OMOSNH2KD  ! OMEGA / SINH(2KD)
        REAL(KIND=JWRB), DIMENSION(:), ALLOCATABLE :: STOKFAC    ! FACTOR TO COMPUTE SURFACE STOKES DRIFT FROM SPECTRUM 2*G*K**2/(OMEGA*TANH(2KD))
        REAL(KIND=JWRB), DIMENSION(:), ALLOCATABLE :: CIWA       ! SEA ICE WAVE ATTENUATION
        CONTAINS
           PROCEDURE :: ALLOC=>FREQUENCY_LAND_ALLOC
           PROCEDURE :: DEALLOC=>FREQUENCY_LAND_DEALLOC
      END TYPE FREQUENCY_LAND

!*    **  *VARIABLES USED FOR FORCING INPUT AND COMPUTATIONS.
!          See *INIT_FIELDG* for the initialisation

      TYPE FORCING_FIELDS 
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: UWND    ! U COMPONENT ON WAVE MODEL GRID of
                                                                ! 10m wind, or friction velocity or surfact stress
                                                                ! See ICODE and ICODE_CPL
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: VWND    ! V COMPONENT ON WAVE MODEL GRID of
                                                                ! 10m wind, or friction velocity or surfact stress
                                                                ! See ICODE and ICODE_CPL
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: AIRD    ! AIR DENSITY ON WAVE MODEL GRID
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: WSTAR   ! CONVECTIVE VELOCITY ON WAVE MODEL GRID
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: CICOVER ! SEA ICE FRACTION ON WAVE MODEL GRID
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: CITHICK ! SEA ICE THICKNESS ON WAVE MODEL GRID
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: LKFR    ! LAKE FRACTION ON WAVE MODEL GRID
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: USTRA   ! U-COMPONENT OF THE ATMOSPHERIC STRESS OVER THE OCEAN 
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: VSTRA   ! V-COMPONENT OF THE ATMOSPHERIC STRESS OVER THE OCEAN 
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: UCUR    ! U COMPONENT OF CURRENT ON WAVE MODEL GRID
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: VCUR    ! V COMPONENT OF CURRENT ON WAVE MODEL GRID
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: WSWAVE  ! WIND SPEED (WAVE PARAMETER 245) ON WAVE MODEL GRID
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: WDWAVE  ! WIND DIRECTION (WAVE PARAMETER 249) ON WAVE MODEL GRID
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: UFRIC   ! FRICTION VELOCITY
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: TAUW    ! WAVE INDUCED KINEMATIC STRESS MAGNITUDE
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: TAUWDIR ! WAVE INDUCED KINEMATIC STRESS DIRECTION
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: Z0M     ! SURFACE ROUGHNESS LENGTH SCALE
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: Z0B     ! BACKGROUND SURFACE ROUGHNESS LENGTH SCALE
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: CHRNCK  ! CHARNOCK COEFFICIENT
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: XLON    ! LONGITUDE OF FORCING_FIELDS DATA THAT ARE NEEDED (i.e. all local points only !)
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: YLAT    ! LATITUDE  OF FORCING_FIELDS DATA THAT ARE NEEDED (i.e. all local points only !)
        CONTAINS
           PROCEDURE :: ALLOC=>FORCING_FIELDS_ALLOC
           PROCEDURE :: DEALLOC=>FORCING_FIELDS_DEALLOC
      END TYPE FORCING_FIELDS


!*    **  *VARIABLES USED FOR OUTPUT OF INTEGRATED/DERIVED 2D FIELDS
!          SEE *WVALLOC* AND *WVDEALLOC* FOR THEIR MANAGEMENT

      TYPE INTGT_PARAM_FIELDS 
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: WSEMEAN  ! WINDSEA VARIANCE.
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: WSFMEAN  ! WINDSEA MEAN FREQUENCY (1./MEAN PERIOD).

        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: USTOKES  ! U-COMP SURFACE STOKES DRIFT.
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: VSTOKES  ! V-COMP SURFACE STOKES DRIFT.

        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: PHIEPS   ! NORMALIZED ENERGY FLUX TO OCEAN.
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: PHIOCD   ! DIMENSIONAL TURBULENT ENERGY FLUX INTO OCEAN.

        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: PHIAW    ! NORMALIZED ENERGY FLUX FROM WIND TO WAVES.

        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: TAUOC    ! NORMALIZED MOMENTUM FLUX INTO OCEAN.
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: TAUXD    ! DIMENSIONAL U-COMPONENT OF MOMENTUM FLUX FROM ATMOSPHERE.
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: TAUYD    ! DIMENSIONAL V-COMPONENT OF MOMENTUM FLUX FROM ATMOSPHERE.
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: TAUOCXD  ! DIMENSIONAL U-COMPONENT OF MOMENTUM FLUX INTO OCEAN.
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: TAUOCYD  ! DIMENSIONAL V-COMPONENT OF MOMENTUM FLUX INTO OCEAN.

        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: STRNMS   ! MEAN SQUARE STRAIN INTO THE SEA ICE.

        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: ALTWH    ! ALTIMETER WAVE HEIGHT
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: CALTWH   ! CORRECTED ALTIMETER WAVE HEIGHT
        REAL(KIND=JWRB), DIMENSION(:,:), ALLOCATABLE :: RALTCOR  ! ALTIMETER RANGE CORRECTION
        CONTAINS
           PROCEDURE :: ALLOC=>INTGT_PARAM_FIELDS_ALLOC
           PROCEDURE :: DEALLOC=>INTGT_PARAM_FIELDS_DEALLOC
      END TYPE INTGT_PARAM_FIELDS 


!*    **  *VARIABLES USED FOR COLLECTING FIELDS PASSED FROM THE WAVES TO THE OCEAN

      TYPE WAVE2OCEAN
        REAL(KIND=JWRO), DIMENSION(:,:), ALLOCATABLE :: NSWH       ! SIGNIFICANT WAVE HEIGHT
        REAL(KIND=JWRO), DIMENSION(:,:), ALLOCATABLE :: NMWP       ! MEAN WAVE PERIOD
        REAL(KIND=JWRO), DIMENSION(:,:), ALLOCATABLE :: NPHIEPS    ! NORMALIZED TURBULENT KINETIC ENERGY FLUX INTO OCEAN
        REAL(KIND=JWRO), DIMENSION(:,:), ALLOCATABLE :: NEMOPHIF   ! TURBULENT KINETIC ENERGY FLUX INTO THE OCEAN
        REAL(KIND=JWRO), DIMENSION(:,:), ALLOCATABLE :: NTAUOC     ! NORMALIZED MOMENTUM FLUX INTO OCEAN
        REAL(KIND=JWRO), DIMENSION(:,:), ALLOCATABLE :: NEMOTAUX   ! U-COMPONENT OF OCEAN STRESS (MOMENTUM FLUX)
        REAL(KIND=JWRO), DIMENSION(:,:), ALLOCATABLE :: NEMOTAUY   ! V-COMPONENT OF OCEAN STRESS (MOMENTUM FLUX)
        REAL(KIND=JWRO), DIMENSION(:,:), ALLOCATABLE :: NEMOUSTOKES! U-COMPONENT OF SURFACE STOKES DRIFT
        REAL(KIND=JWRO), DIMENSION(:,:), ALLOCATABLE :: NEMOVSTOKES! V-COMPONENT OF SURFACE STOKES DRIFT
        REAL(KIND=JWRO), DIMENSION(:,:), ALLOCATABLE :: NEMOWSWAVE ! WIND SPEED USED TO GENERATE WAVES
        REAL(KIND=JWRO), DIMENSION(:,:), ALLOCATABLE :: NEMOSTRN   ! SEA ICE MEAN SQUARE WAVE STRAIN
        CONTAINS
           PROCEDURE :: ALLOC=>WAVE2OCEAN_ALLOC
           PROCEDURE :: DEALLOC=>WAVE2OCEAN_DEALLOC
      END TYPE WAVE2OCEAN


!*    **  *VARIABLES USED FOR COLLECTING FIELDS PASSED FROM THE OCEAN TO THE WAVES

      TYPE OCEAN2WAVE
        REAL(KIND=JWRO), DIMENSION(:,:), ALLOCATABLE :: NEMOSST     ! SEA SURFACE TEMPERATURE
        REAL(KIND=JWRO), DIMENSION(:,:), ALLOCATABLE :: NEMOCICOVER ! SEA ICE COVER
        REAL(KIND=JWRO), DIMENSION(:,:), ALLOCATABLE :: NEMOCITHICK ! SEA ICE THICKNESS
        REAL(KIND=JWRO), DIMENSION(:,:), ALLOCATABLE :: NEMOUCUR    ! ZONAL CURRENT
        REAL(KIND=JWRO), DIMENSION(:,:), ALLOCATABLE :: NEMOVCUR    ! MERIDIONAL CURRENT
        CONTAINS
           PROCEDURE :: ALLOC=>OCEAN2WAVE_ALLOC
           PROCEDURE :: DEALLOC=>OCEAN2WAVE_DEALLOC
      END TYPE OCEAN2WAVE

! ----------------------------------------------------------------------
      CONTAINS
           SUBROUTINE WVGRIDGLO_ALLOC(SELF, NIBLO)
              CLASS(WVGRIDGLO), INTENT(INOUT) :: SELF
              INTEGER(KIND=JWIM), INTENT(IN) :: NIBLO

              ALLOCATE(SELF%IXLG(NIBLO))
              ALLOCATE(SELF%KXLT(NIBLO))
           END SUBROUTINE

           SUBROUTINE WVGRIDGLO_DEALLOC(SELF)
              CLASS(WVGRIDGLO), INTENT(INOUT) :: SELF

              DEALLOCATE(SELF%IXLG)
              DEALLOCATE(SELF%KXLT)
           END SUBROUTINE

           SUBROUTINE WVGRIDLOC_ALLOC(SELF, NPROMA_WAM, NCHNK)
              CLASS(WVGRIDLOC), INTENT(INOUT) :: SELF
              INTEGER(KIND=JWIM), INTENT(IN) :: NPROMA_WAM, NCHNK

              ALLOCATE(SELF%IFROMIJ(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%JFROMIJ(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%KFROMIJ(NPROMA_WAM, NCHNK))
           END SUBROUTINE

           SUBROUTINE WVGRIDLOC_DEALLOC(SELF)
              CLASS(WVGRIDLOC), INTENT(INOUT) :: SELF

              DEALLOCATE(SELF%IFROMIJ)
              DEALLOCATE(SELF%JFROMIJ)
              DEALLOCATE(SELF%KFROMIJ)
           END SUBROUTINE

           SUBROUTINE ENVIRONMENT_ALLOC(SELF, NPROMA_WAM, NCHNK)
              CLASS(ENVIRONMENT), INTENT(INOUT) :: SELF
              INTEGER(KIND=JWIM), INTENT(IN) :: NPROMA_WAM, NCHNK

              ALLOCATE(SELF%INDEP(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%IODP(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%IOBND(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%DELLAM1(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%COSPHM1(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%DEPTH(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%EMAXDPT(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%UCUR(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%VCUR(NPROMA_WAM, NCHNK))
           END SUBROUTINE

           SUBROUTINE ENVIRONMENT_DEALLOC(SELF)
              CLASS(ENVIRONMENT), INTENT(INOUT) :: SELF

              DEALLOCATE(SELF%INDEP)
              DEALLOCATE(SELF%IODP)
              DEALLOCATE(SELF%IOBND)
              DEALLOCATE(SELF%DELLAM1)
              DEALLOCATE(SELF%COSPHM1)
              DEALLOCATE(SELF%DEPTH)
              DEALLOCATE(SELF%EMAXDPT)
              DEALLOCATE(SELF%UCUR)
              DEALLOCATE(SELF%VCUR)
           END SUBROUTINE

           SUBROUTINE FREQUENCY_LAND_ALLOC(SELF, NFRE)
              CLASS(FREQUENCY_LAND), INTENT(INOUT) :: SELF
              INTEGER(KIND=JWIM), INTENT(IN) :: NFRE

              ALLOCATE(SELF%WAVNUM(NFRE))
              ALLOCATE(SELF%CINV(NFRE))
              ALLOCATE(SELF%CGROUP(NFRE))
              ALLOCATE(SELF%XK2CG(NFRE))
              ALLOCATE(SELF%OMOSNH2KD(NFRE))
              ALLOCATE(SELF%STOKFAC(NFRE))
              ALLOCATE(SELF%CIWA(NFRE))
           END SUBROUTINE

           SUBROUTINE FREQUENCY_LAND_DEALLOC(SELF)
              CLASS(FREQUENCY_LAND), INTENT(INOUT) :: SELF

              DEALLOCATE(SELF%WAVNUM)
              DEALLOCATE(SELF%CINV)
              DEALLOCATE(SELF%CGROUP)
              DEALLOCATE(SELF%XK2CG)
              DEALLOCATE(SELF%OMOSNH2KD)
              DEALLOCATE(SELF%STOKFAC)
              DEALLOCATE(SELF%CIWA)
           END SUBROUTINE

           SUBROUTINE FREQUENCY_ALLOC(SELF, NPROMA_WAM, NFRE, NCHNK)
              CLASS(FREQUENCY), INTENT(INOUT) :: SELF
              INTEGER(KIND=JWIM), INTENT(IN) :: NPROMA_WAM, NFRE, NCHNK

              ALLOCATE(SELF%WAVNUM(NPROMA_WAM, NFRE, NCHNK))
              ALLOCATE(SELF%CINV(NPROMA_WAM, NFRE, NCHNK))
              ALLOCATE(SELF%CGROUP(NPROMA_WAM, NFRE, NCHNK))
              ALLOCATE(SELF%XK2CG(NPROMA_WAM, NFRE, NCHNK))
              ALLOCATE(SELF%OMOSNH2KD(NPROMA_WAM, NFRE, NCHNK))
              ALLOCATE(SELF%STOKFAC(NPROMA_WAM, NFRE, NCHNK))
              ALLOCATE(SELF%CIWA(NPROMA_WAM, NFRE, NCHNK))
           END SUBROUTINE

           SUBROUTINE FREQUENCY_DEALLOC(SELF)
              CLASS(FREQUENCY), INTENT(INOUT) :: SELF

              DEALLOCATE(SELF%WAVNUM)
              DEALLOCATE(SELF%CINV)
              DEALLOCATE(SELF%CGROUP)
              DEALLOCATE(SELF%XK2CG)
              DEALLOCATE(SELF%OMOSNH2KD)
              DEALLOCATE(SELF%STOKFAC)
              DEALLOCATE(SELF%CIWA)
           END SUBROUTINE

           SUBROUTINE FORCING_FIELDS_ALLOC(SELF, NPROMA_WAM, NCHNK, UBND0, UBND1)
              CLASS(FORCING_FIELDS), INTENT(INOUT) :: SELF
              INTEGER(KIND=JWIM), INTENT(IN) :: NPROMA_WAM, NCHNK
              INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL :: UBND0, UBND1

              IF (PRESENT(UBND0)) THEN
                IF (.NOT. PRESENT(UBND1))THEN
                   ERROR STOP
                ENDIF
                ALLOCATE(SELF%UWND(NPROMA_WAM:UBND0, NCHNK:UBND1))
                ALLOCATE(SELF%VWND(NPROMA_WAM:UBND0, NCHNK:UBND1))
                ALLOCATE(SELF%AIRD(NPROMA_WAM:UBND0, NCHNK:UBND1))
                ALLOCATE(SELF%WSTAR(NPROMA_WAM:UBND0, NCHNK:UBND1))
                ALLOCATE(SELF%CICOVER(NPROMA_WAM:UBND0, NCHNK:UBND1))
                ALLOCATE(SELF%CITHICK(NPROMA_WAM:UBND0, NCHNK:UBND1))
                ALLOCATE(SELF%LKFR(NPROMA_WAM:UBND0, NCHNK:UBND1))
                ALLOCATE(SELF%USTRA(NPROMA_WAM:UBND0, NCHNK:UBND1))
                ALLOCATE(SELF%VSTRA(NPROMA_WAM:UBND0, NCHNK:UBND1))
                ALLOCATE(SELF%UCUR(NPROMA_WAM:UBND0, NCHNK:UBND1))
                ALLOCATE(SELF%VCUR(NPROMA_WAM:UBND0, NCHNK:UBND1))
                ALLOCATE(SELF%WSWAVE(NPROMA_WAM:UBND0, NCHNK:UBND1))
                ALLOCATE(SELF%WDWAVE(NPROMA_WAM:UBND0, NCHNK:UBND1))
                ALLOCATE(SELF%UFRIC(NPROMA_WAM:UBND0, NCHNK:UBND1))
                ALLOCATE(SELF%TAUW(NPROMA_WAM:UBND0, NCHNK:UBND1))
                ALLOCATE(SELF%TAUWDIR(NPROMA_WAM:UBND0, NCHNK:UBND1))
                ALLOCATE(SELF%Z0M(NPROMA_WAM:UBND0, NCHNK:UBND1))
                ALLOCATE(SELF%Z0B(NPROMA_WAM:UBND0, NCHNK:UBND1))
                ALLOCATE(SELF%CHRNCK(NPROMA_WAM:UBND0, NCHNK:UBND1))
                ALLOCATE(SELF%XLON(NPROMA_WAM:UBND0, NCHNK:UBND1))
                ALLOCATE(SELF%YLAT(NPROMA_WAM:UBND0, NCHNK:UBND1))
              ELSE
                ALLOCATE(SELF%UWND(NPROMA_WAM, NCHNK))
                ALLOCATE(SELF%VWND(NPROMA_WAM, NCHNK))
                ALLOCATE(SELF%AIRD(NPROMA_WAM, NCHNK))
                ALLOCATE(SELF%WSTAR(NPROMA_WAM, NCHNK))
                ALLOCATE(SELF%CICOVER(NPROMA_WAM, NCHNK))
                ALLOCATE(SELF%CITHICK(NPROMA_WAM, NCHNK))
                ALLOCATE(SELF%LKFR(NPROMA_WAM, NCHNK))
                ALLOCATE(SELF%USTRA(NPROMA_WAM, NCHNK))
                ALLOCATE(SELF%VSTRA(NPROMA_WAM, NCHNK))
                ALLOCATE(SELF%UCUR(NPROMA_WAM, NCHNK))
                ALLOCATE(SELF%VCUR(NPROMA_WAM, NCHNK))
                ALLOCATE(SELF%WSWAVE(NPROMA_WAM, NCHNK))
                ALLOCATE(SELF%WDWAVE(NPROMA_WAM, NCHNK))
                ALLOCATE(SELF%UFRIC(NPROMA_WAM, NCHNK))
                ALLOCATE(SELF%TAUW(NPROMA_WAM, NCHNK))
                ALLOCATE(SELF%TAUWDIR(NPROMA_WAM, NCHNK))
                ALLOCATE(SELF%Z0M(NPROMA_WAM, NCHNK))
                ALLOCATE(SELF%Z0B(NPROMA_WAM, NCHNK))
                ALLOCATE(SELF%CHRNCK(NPROMA_WAM, NCHNK))
                ALLOCATE(SELF%XLON(NPROMA_WAM, NCHNK))
                ALLOCATE(SELF%YLAT(NPROMA_WAM, NCHNK))
              ENDIF
           END SUBROUTINE FORCING_FIELDS_ALLOC

           SUBROUTINE FORCING_FIELDS_DEALLOC(SELF)
              CLASS(FORCING_FIELDS), INTENT(INOUT) :: SELF

              DEALLOCATE(SELF%UWND)
              DEALLOCATE(SELF%VWND)
              DEALLOCATE(SELF%AIRD)
              DEALLOCATE(SELF%WSTAR)
              DEALLOCATE(SELF%CICOVER)
              DEALLOCATE(SELF%CITHICK)
              DEALLOCATE(SELF%LKFR)
              DEALLOCATE(SELF%USTRA)
              DEALLOCATE(SELF%VSTRA)
              DEALLOCATE(SELF%UCUR)
              DEALLOCATE(SELF%VCUR)
              DEALLOCATE(SELF%WSWAVE)
              DEALLOCATE(SELF%WDWAVE)
              DEALLOCATE(SELF%UFRIC)
              DEALLOCATE(SELF%TAUW)
              DEALLOCATE(SELF%TAUWDIR)
              DEALLOCATE(SELF%Z0M)
              DEALLOCATE(SELF%Z0B)
              DEALLOCATE(SELF%CHRNCK)
              DEALLOCATE(SELF%XLON)
              DEALLOCATE(SELF%YLAT)
           END SUBROUTINE FORCING_FIELDS_DEALLOC

           SUBROUTINE INTGT_PARAM_FIELDS_ALLOC(SELF, NPROMA_WAM, NCHNK)
              CLASS(INTGT_PARAM_FIELDS), INTENT(INOUT) :: SELF
              INTEGER(KIND=JWIM), INTENT(IN) :: NPROMA_WAM, NCHNK

              ALLOCATE(SELF%WSEMEAN(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%WSFMEAN(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%USTOKES(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%VSTOKES(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%PHIEPS(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%PHIOCD(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%PHIAW(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%TAUOC(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%TAUXD(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%TAUYD(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%TAUOCXD(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%TAUOCYD(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%STRNMS(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%ALTWH(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%CALTWH(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%RALTCOR(NPROMA_WAM, NCHNK))
           END SUBROUTINE INTGT_PARAM_FIELDS_ALLOC

           SUBROUTINE INTGT_PARAM_FIELDS_DEALLOC(SELF)
              CLASS(INTGT_PARAM_FIELDS), INTENT(INOUT) :: SELF

              DEALLOCATE(SELF%WSEMEAN)
              DEALLOCATE(SELF%WSFMEAN)
              DEALLOCATE(SELF%USTOKES)
              DEALLOCATE(SELF%VSTOKES)
              DEALLOCATE(SELF%PHIEPS)
              DEALLOCATE(SELF%PHIOCD)
              DEALLOCATE(SELF%PHIAW)
              DEALLOCATE(SELF%TAUOC)
              DEALLOCATE(SELF%TAUXD)
              DEALLOCATE(SELF%TAUYD)
              DEALLOCATE(SELF%TAUOCXD)
              DEALLOCATE(SELF%TAUOCYD)
              DEALLOCATE(SELF%STRNMS)
              DEALLOCATE(SELF%ALTWH)
              DEALLOCATE(SELF%CALTWH)
              DEALLOCATE(SELF%RALTCOR)
           END SUBROUTINE INTGT_PARAM_FIELDS_DEALLOC

           SUBROUTINE WAVE2OCEAN_ALLOC(SELF, NPROMA_WAM, NCHNK)
              CLASS(WAVE2OCEAN), INTENT(INOUT) :: SELF
              INTEGER(KIND=JWIM), INTENT(IN) :: NPROMA_WAM, NCHNK

              ALLOCATE(SELF%NSWH(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%NMWP(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%NPHIEPS(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%NEMOPHIF(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%NTAUOC(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%NEMOTAUX(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%NEMOTAUY(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%NEMOUSTOKES(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%NEMOVSTOKES(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%NEMOWSWAVE(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%NEMOSTRN(NPROMA_WAM, NCHNK))
           END SUBROUTINE WAVE2OCEAN_ALLOC

           SUBROUTINE WAVE2OCEAN_DEALLOC(SELF)
              CLASS(WAVE2OCEAN), INTENT(INOUT) :: SELF

              DEALLOCATE(SELF%NSWH)
              DEALLOCATE(SELF%NMWP)
              DEALLOCATE(SELF%NPHIEPS)
              DEALLOCATE(SELF%NEMOPHIF)
              DEALLOCATE(SELF%NTAUOC)
              DEALLOCATE(SELF%NEMOTAUX)
              DEALLOCATE(SELF%NEMOTAUY)
              DEALLOCATE(SELF%NEMOUSTOKES)
              DEALLOCATE(SELF%NEMOVSTOKES)
              DEALLOCATE(SELF%NEMOWSWAVE)
              DEALLOCATE(SELF%NEMOSTRN)
           END SUBROUTINE WAVE2OCEAN_DEALLOC

           SUBROUTINE OCEAN2WAVE_ALLOC(SELF, NPROMA_WAM, NCHNK)
              CLASS(OCEAN2WAVE), INTENT(INOUT) :: SELF
              INTEGER(KIND=JWIM), INTENT(IN) :: NPROMA_WAM, NCHNK

              ALLOCATE(SELF%NEMOSST(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%NEMOCICOVER(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%NEMOCITHICK(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%NEMOUCUR(NPROMA_WAM, NCHNK))
              ALLOCATE(SELF%NEMOVCUR(NPROMA_WAM, NCHNK))
           END SUBROUTINE OCEAN2WAVE_ALLOC

           SUBROUTINE OCEAN2WAVE_DEALLOC(SELF)
              CLASS(OCEAN2WAVE), INTENT(INOUT) :: SELF

              DEALLOCATE(SELF%NEMOSST)
              DEALLOCATE(SELF%NEMOCICOVER)
              DEALLOCATE(SELF%NEMOCITHICK)
              DEALLOCATE(SELF%NEMOUCUR)
              DEALLOCATE(SELF%NEMOVCUR)
           END SUBROUTINE OCEAN2WAVE_DEALLOC
      END MODULE YOWDRVTYPE
