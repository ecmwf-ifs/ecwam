! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWICE

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*     ** *ICE* ICE POINTS


      INTEGER(KIND=JWIM) :: IPARAMCI
      INTEGER(KIND=JWIM) :: NICT, NICH

      REAL(KIND=JWRB), PARAMETER   :: FLMIN=0.00001_JWRB ! MINIMUM ENERGY IN SPECTRAL BINS
      REAL(KIND=JWRB), PARAMETER   :: HICMIN=0.2_JWRB
      REAL(KIND=JWRB)              :: CITHRSH
      REAL(KIND=JWRB)              :: CIBLOCK
      REAL(KIND=JWRB)              :: CITHRSH_SAT
      REAL(KIND=JWRB)              :: CITHRSH_TAIL
      REAL(KIND=JWRB)              :: CDICWA
      REAL(KIND=JWRB)              :: TICMIN
      REAL(KIND=JWRB)              :: DTIC, DHIC
      REAL(KIND=JWRB)              :: ZALPFACB
      REAL(KIND=JWRB)              :: ZALPFACX
      REAL(KIND=JWRB)              :: ZALPWRS
      REAL(KIND=JWRB)              :: ZIBRW_THRSH
      REAL(KIND=JWRB), ALLOCATABLE :: CIDEAC(:,:) 

      LOGICAL :: LICERUN 
      LOGICAL :: LICETH
      LOGICAL :: LMASKICE
      LOGICAL :: LWAMRSETCI
      LOGICAL :: LCIWA1
      LOGICAL :: LCIWA2
      LOGICAL :: LCIWA3
      LOGICAL :: LCISCAL

!--------------------------------------------------------------------

!*    VARIABLE     TYPE      PURPOSE
!     --------     ----      -------
!     IPARAMCI     INTEGER   GRIB PARAMETER VALUE OF SEA ICE FRACTION OR SST.
!     NICT         INTEGER   FIRST DIMENSION (WAVE PERIOD) OF TABLE CIDEAC.
!     NICH         INTEGER   SECOND DIMENSION (ICE THICKNESS) OF TABLE CIDEAC.
!     FLMIN        REAL      ABSOLUTE MINIMUM VALUE
!                            OF EACH SPECTRAL COMPONENTS.
!     CITHRSH      REAL      SEA ICE TRESHOLD, ALL SEA POINTS WITH CICOVER > CITHRSH
!                            WILL BE SET TO MISSING 
!                            WHEN LMASKICE IS TRUE IT IS DONE AT ALL TIME STEP
!                            WHEN IT IS FALSE IT IS !ONLY! DONE FOR THE OUTPUT OF
!                            THE WAVE INTEGRATED PARAMETERS.
!     CIBLOCK      REAL      IF LMASKICE : FULL BLOCKING SEA ICE TRESHOLD, ALL SEA POINTS WITH CICOVER > CIBLOCK
!                            AND THRESHOLD OVER WHICH FIELDS THAT ARE EXCHANGED WITH THE ATMOSPHERE AND THE OCEAN 
!                            ARE RESET.
!     CITHRSH_SAT  REAL      SEA ICE TRESHOLD FOR DATA ASSIMILATION, WHICH WILL ONLY BE DONE FOR
!                            ALL SEA POINTS WITH CICOVER < CITHRSH_SAT
!     CITHRSH_TAIL REAL      SEA ICE TRESHOLD FOR IMPOSITION OF SPECTRAL TAIL,
!                            FOR ALL SEA POINTS WITH CICOVER < CITHRSH_TAIL
!     CDICWA       REAL      DRAG COEFFICIENT ICE-WATER.
!     TICMIN       REAL      MINIMUM WAVE PERIOD IN TABLE CIDEAC.
!     HICMIN       REAL      MINIMUM ICE THICKNESS IN TABLE CIDEAC.
!     DTIC         REAL      WAVE PERIOD INCREMENT IN TABLE CIDEAC.
!     DHIC         REAL      ICE THICKNESS INCREMENT IN TABLE CIDEAC.
!     CIDEAC       REAL      SEA ICE DIMENSIONLESS ENERGY ATTENUATION COEFFICIENT TABLE. 
!     LICERUN      LOGICAL   SET TO TRUE IF SEA ICE IS TAKEN INTO ACCOUNT.
!     LICETH       LOGICAL   SET TO TRUE IF SEA ICE THICKNESS IS PART OF THE INPUT DATA.
!     LMASKICE     LOGICAL   SET TO TRUE IF ICE MASK IS APPLIED (SEE CITHRSH)
!                            NOTE THAT THE MASK IS ALWAYS APPLIED FOR SATELLITE
!                            DATA PROCESSING AS WELL AS FOR CHARNOCK THAT
!                            IS RETURNED TO THE IFS.
!     LWAMRSETCI   LOGICAL   SET TO TRUE IF FIELDS THAT ARE EXCHANGED WITH THE ATMOSPHERE AND THE OCEAN
!                            ARE RESET TO WHAT WOULD BE USED IF THERE WERE NO WAVE MODELS.
!                            IT IS ONLY ACTIVE FOR CICOVER>CIBLOCK OR CICOVER>CITHRSH
!     LCIWA1       LOGICAL   SET TO TRUE IF SEA ICE SCATTERING ATTENUATION IS USED.
!     LCIWA2       LOGICAL   SET TO TRUE IF SEA ICE BOTTOM FRICTION ATTENUATION IS USED.
!     LCIWA3       LOGICAL   SET TO TRUE IF SEA ICE VISCOUS FRICTION ATTENUATION IS USED.
!     LCISCAL      LOGICAL   SET TO TRUE TO SCALE INPUT AND DISSIPATION SOURCE TERMS BY SEA ICE CONCENTRATION
!     ZALPFACB     REAL      FACTOR TO SCALE ATTENUATION FOR ALL SEA ICE
!     ZALPFACX     REAL      FACTOR TO SCALE ATTENUATION UP/DOWN FOR SOLID/BROKEN ICE
!     ZALPWRS      REAL      PROPORTION OF ENERGY LOST FROM SLICE THAT GOES INTO WAVE RADIATIVE STRESS
!     ZIBRW_THRSH  REAL      THRESHOLD AT WHICH SEA ICE IS CONSIDERED BROKEN
!--------------------------------------------------------------------
      END MODULE YOWICE
