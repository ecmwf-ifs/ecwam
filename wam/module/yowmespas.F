      MODULE YOWMESPAS

      IMPLICIT NONE

!*    **  *YOWMESPAS* - SELECTS MACHINE SPECIFIC PART OF THE CODE 
!                    AND SOME LATER CONTROL SWITCHES
      LOGICAL LMESSPASS 
      LOGICAL LPBIOIN 
      LOGICAL LPBIOOUT 
      LOGICAL LFDBIOIN 
      LOGICAL LFDBIOOUT 
      LOGICAL LGRIBIN 
      LOGICAL LGRIBOUT 
      LOGICAL LNOCDIN 

!*     VARIABLE.   TYPE.     PURPOSE.
!      --------    ----      -------

!      *LMESSPASS  LOGICAL   SELECTS PART OF THE CODE THAT IS SPECIFIC
!                            TO A MESSAGE PASSING OR A MULTI-TASKED 
!                            ENVIRONMENT.
!      *LPBIOIN    LOGICAL   SELECTS TYPE OF INPUT OF BINARY FILES,
!                            IF TRUE THE PBIO SOFTWARE WILL BE USED FOR
!                            INPUT.
!      *LPBIOOUT   LOGICAL   SELECTS TYPE OF OUTPUT OF BINARY FILES,
!                            IF TRUE THE PBIO SOFTWARE WILL BE USED FOR
!                            OUTPUT.
!      *LFDBIOIN   LOGICAL   IF TRUE THE FDB SOFTWARE WILL BE USED FOR 
!                            INPUT.
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
! ----------------------------------------------------------------------
      END MODULE YOWMESPAS
