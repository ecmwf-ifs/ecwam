      SUBROUTINE WINDCOR(
!                DIMENSION
     &                NSPECW,MPART,
!                INPUT
     &                EW,XMEANFREW,XMEANANGW,
     &                EA,XFREA,XANGA,
     &                NINTW, LOTEST, IU06,
!                MODIFIED
     &                DUVSAR)

! ----------------------------------------------------------------------

!    S.HASSELMANN, MPM, HAMBURG 3/94.
!    R.BROKOPF, MPM, HAMBURG 3/94.
!    P.LIONELLO (UNIV. OF PADUA), HAMBURG 3/94
!    J. BIDLOT, ECMWF : USE USTAR INSTEAD

!     PURPOSE.
!     --------
!     CORRECT LOCAL WINDS FROM LOCAL WAVE CORRECTIONS.

!     EXTERNALS.
!     ----------
!      U10COR - COMPUTES RELATION BETWEEN U10 AND ENERGY , 
!               MEAN FREQUENCY.

!     METHOD:
!     -------
!     A DIMENSIONLESS RELATION IN THE FORM
!          E10= A * FR10 **B,                            (1)
!     IT IS ASSUMED WITH
!          E10=E G**2/USTAR**4  (DIMENSIONLESS ENERGY)     (2)
!     AND
!          FR10=FR * USTAR/G (DIMENSIONLESS FREQUENCY).    (3)
!     USTAR is THE FRICTION VELOCITY, I.E. THE INPUT OF WAM.

!     SUBSTITUITING E10 AND FR10 IN (1) ACCORDING TO (2) AND (3),
!     ONE HAS AN EQUATION FOR THE UNKNOWN USTAR. WITH GIVES

!      USTAR= XX**(-1/(4+B)) ; XX = A F**B / (E G**(B+2).

!     THE CONSTANT A AND B HAVE BEEN OBTAINED FROM THE WAM TIME
!     LIMITED GROWTH CURVE.

!     IN ORDER TO AVOID A TOO STRONG DEPENDANCE ON THE GROWTH CURVE
!     USED FOR THE ESTIMATE,THE DEVIATION DU10=U1NE-U10OE IS USED:
!     USAN=U10OLD+DU10
!     WHERE USAN IS THE ANALYSED WIND SPEED.

!     THE ANALYSIS WIND INCREMENT ARE PASSED TO MERGESARCOR WHERE
!     A SIMPLE INTERPOLATION IS USED TO SMOOTH OUT THE RESULTS
!    

! ----------------------------------------------------------------------

      IMPLICIT NONE

!     INTERFACE

      INTEGER 
     &   NSPECW, MPART

      REAL 
     &   EW(NSPECW,MPART),
     &   XMEANFREW(NSPECW,MPART),
     &   XMEANANGW(NSPECW,MPART),
     &   EA(NSPECW,MPART),
     &   XFREA(NSPECW,MPART),
     &   XANGA(NSPECW,MPART),
     &   DUVSAR(2,NSPECW)

      REAL 
     &   USFG(NSPECW),
     &   USA(NSPECW)

      INTEGER 
     &   NINTW(NSPECW)
      REAL
     &  UVXW,UVYW,UVXA,UVYA
 
      LOGICAL LOTEST

      INTEGER IU06

!     LOCAL VARIABLES:
      INTEGER ISPEC,JPART

! ----------------------------------------------------------------------

      IF(LOTEST) WRITE(IU06,*)' INSIDE WINDCOR '

!*    ESTIMATE OF THE FIRST GUESS USTAR.
!     ---------------------------------

      CALL U10COR (USFG, XMEANFREW, EW, NSPECW, MPART, NINTW)

!*    ESTIMATE OF THE ANALYSED USTAR.
!     -------------------------------

      CALL U10COR (USA, XFREA, EA, NSPECW, MPART, NINTW)

!*    COMPUTE WIND CORRECTIONS.
!     -------------------------

      DO ISPEC=1,NSPECW
        IF(NINTW(ISPEC).NE.0 .AND.
     &     USFG(ISPEC).GT.0. .AND. USA(ISPEC) .GT.0.) THEN
            JPART=NINTW(ISPEC)
            UVXW=USFG(ISPEC)*COS(XMEANANGW(ISPEC,JPART))
            UVYW=USFG(ISPEC)*SIN(XMEANANGW(ISPEC,JPART))
            UVXA=USA(ISPEC)*COS(XANGA(ISPEC,JPART))
            UVYA=USA(ISPEC)*SIN(XANGA(ISPEC,JPART))
            DUVSAR(1,ISPEC) = UVXA-UVXW
            DUVSAR(2,ISPEC) = UVYA-UVYW
        END IF
      END DO

      RETURN
      END SUBROUTINE WINDCOR
