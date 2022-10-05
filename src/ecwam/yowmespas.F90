      MODULE YOWMESPAS

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!*    **  *YOWMESPAS* - SELECTS MACHINE SPECIFIC PART OF THE CODE 
!                    AND SOME LATER CONTROL SWITCHES
      LOGICAL :: LFDBIOOUT 
      LOGICAL :: LGRIBIN 
      LOGICAL :: LGRIBOUT 
      LOGICAL :: LNOCDIN 
      LOGICAL :: LWAVEWIND=.FALSE. 

!*     VARIABLE.   TYPE.     PURPOSE.
!      --------    ----      -------

!      *LFDBIOOUT  LOGICAL   IF TRUE THE PBIO SOFTWARE WILL BE USED FOR 
!                            OUTPUT.
!      *LGRIBIN    LOGICAL   SELECTS TYPE OF INPUT FILES
!                            IF TRUE INPUT GRIB SPECTRA , CD, WINDS
!                            ARE EXPECTED AS INPUTS. EXCEPT IN preset
!                            WHERE GRIB SPECTRA ARE ONLY EXPECTED.
!      *LGRIBOUT   LOGICAL   SELECTS TYPE OF OUTPUT FILES
!                            IF TRUE OUTPUT GRIB SPECTRA , CD, WINDS
!                            ARE EXPECTED AS OUTPUTS. EXCEPT IN preset
!                            WHERE GRIB SPECTRA ARE ONLY EXPECTED.
!      *LNOCDIN    LOGICAL   SELECTS WHETHER A DRAG COEFFICIENT FILE IS
!                            SUPPLIED AS INPUT. 
!      *LWAVEWIND  LOGICAL   TELLS WHETHER A WAVE 10M WIND FILE IS 
!                            SUPPLIED AS INPUT. 
! ----------------------------------------------------------------------
      END MODULE YOWMESPAS
