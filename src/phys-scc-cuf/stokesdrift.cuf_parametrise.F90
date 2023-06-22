! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE STOKESDRIFT_CUF_PARAMETRISE_MOD
  CONTAINS
  ATTRIBUTES(DEVICE) SUBROUTINE STOKESDRIFT_CUF_PARAMETRISE (KIJS, KIJL, FL1, STOKFAC, WSWAVE, WDWAVE, CICOVER, USTOKES,  &
  & VSTOKES, IJ)
    
    !
    !***  *STOKESDRIFT*   DETERMINES THE STOKES DRIFT
    !
    !     PETER JANSSEN MARCH 2009
    !
    !     PURPOSE.
    !     --------
    !
    !              DETERMINATION OF STOKES DRIFT VECTOR
    !
    !     INTERFACE.
    !     ----------
    !              *CALL*  *STOKESDRIFT(KIJS, KIJL, FL1, STOKFAC, WSWAVE,WDWAVE,CICOVER,USTOKES,VSTOKES)*
    !
    !                       INPUT:
    !                            *KIJS*   - FIRST GRIDPOINT
    !                            *KIJL*   - LAST GRIDPOINT
    !                            *FL1*    - 2-D SPECTRUM
    !                            *STOKFAC*- FACTOR TO COMPUTE THE STOKES DRIFT
    !                            Auxilliary fields to specify Stokes when model sea ice cover the blocking threshold
    !                            as 0.016*WSWAVE, aligned in the wind direction
    !                            *WSWAVE* - WIND SPEED IN M/S.
    !                            *WDWAVE* - WIND DIRECTION IN RADIANS.
    !                            *CICOVER*- SEA ICE COVER.
    !
    !                       OUTPUT:
    !                            *USTOKES*   - U-COMPONENT STOKES DRIFT
    !                            *VSTOKES*   - V-COMPONENT STOKES DRIFT
    !
    !     METHOD.
    !     -------
    !              DETERMINE U- AND V-COMPONENT OF STOKES DRIFT FOLLOWING
    !              K.E. KENYON, J.G.R., 74, 6991-6994
    !
    !     EXTERNALS.
    !     ----------
    !              NONE
    !
    !
    !-----------------------------------------------------------------------
    
    USE cudafor
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWPCONS, ONLY: G_D, ZPI_D
    USE YOWFRED, ONLY: FR_D, DFIM_D, DELTH_D, TH_D, DFIM_SIM_D, FRATIO, COSTH_D, SINTH_D
    USE YOWICE, ONLY: LICERUN_D, LWAMRSETCI_D, CITHRSH_D
    USE YOWPARAM, ONLY: NANG_D, NFRE_D, NFRE_ODD_D
    
    
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NANG_LOKI_PARAM = 12
    INTEGER, PARAMETER :: NFRE_LOKI_PARAM = 36
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: KIJS, KIJL
    
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NANG_LOKI_PARAM, NFRE_LOKI_PARAM) :: FL1
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL, NFRE_LOKI_PARAM) :: STOKFAC
    REAL(KIND=JWRB), INTENT(IN), DEVICE, DIMENSION(KIJL) :: WSWAVE, WDWAVE, CICOVER
    REAL(KIND=JWRB), INTENT(OUT), DEVICE, DIMENSION(KIJL) :: USTOKES, VSTOKES
    
    
    INTEGER(KIND=JWIM) :: M, K
    INTEGER(KIND=JWIM), VALUE, INTENT(IN) :: IJ
    
    REAL(KIND=JWRB), PARAMETER :: STMAX = 1.5_JWRB    ! maximum magnitude (this is for safety when coupled)
    REAL(KIND=JWRB) :: CONST, FAC, FAC1, FAC2, FAC3
    REAL(KIND=JWRB) :: STFAC
    
    ! ----------------------------------------------------------------------
    
    
    
    !***  1. DETERMINE STOKE DRIFT VECTOR.
    !     --------------------------------
    
    CONST = 2.0_JWRB*DELTH_D*ZPI_D**3 / G_D*FR_D(NFRE_ODD_D)**4
    
    !***  1.1 PERFORM INTEGRATION.
    !     ------------------------
    
    USTOKES(IJ) = 0.0_JWRB
    VSTOKES(IJ) = 0.0_JWRB
    
    DO M=1,NFRE_ODD_D
      STFAC = STOKFAC(IJ, M)*DFIM_SIM_D(M)
      DO K=1,NANG_D
        FAC3 = STFAC*FL1(IJ, K, M)
        USTOKES(IJ) = USTOKES(IJ) + FAC3*SINTH_D(K)
        VSTOKES(IJ) = VSTOKES(IJ) + FAC3*COSTH_D(K)
      END DO
    END DO
    
    !***  1.2 ADD CONTRIBUTION OF UNRESOLVED WAVES.
    !     -----------------------------------------
    
    DO K=1,NANG_D
      FAC1 = CONST*SINTH_D(K)
      FAC2 = CONST*COSTH_D(K)
      USTOKES(IJ) = USTOKES(IJ) + FAC1*FL1(IJ, K, NFRE_ODD_D)
      VSTOKES(IJ) = VSTOKES(IJ) + FAC2*FL1(IJ, K, NFRE_ODD_D)
    END DO
    
    
    !***  1.3 Sea Ice exception
    !     ---------------------
    IF (LICERUN_D .and. LWAMRSETCI_D) THEN
      IF (CICOVER(IJ) > CITHRSH_D) THEN
        USTOKES(IJ) = 0.016_JWRB*WSWAVE(IJ)*SIN(WDWAVE(IJ))*(1.0_JWRB - CICOVER(IJ))
        VSTOKES(IJ) = 0.016_JWRB*WSWAVE(IJ)*COS(WDWAVE(IJ))*(1.0_JWRB - CICOVER(IJ))
      END IF
    END IF
    
    !***  1.4 Protection
    !     --------------
    
    USTOKES(IJ) = MIN(MAX(USTOKES(IJ), -STMAX), STMAX)
    VSTOKES(IJ) = MIN(MAX(VSTOKES(IJ), -STMAX), STMAX)
    
    
  END SUBROUTINE STOKESDRIFT_CUF_PARAMETRISE
END MODULE STOKESDRIFT_CUF_PARAMETRISE_MOD
