      SUBROUTINE U10COR (USE, FPTMP, EPTMP, NSPECW, MPART, NINTW)
!     PURPOSE.
!     --------
!     COMPUTES RELATION BETWEEN USTAR AND ENERGY , MEAN FREQUENCY.

!     EXTERNALS.
!     ----------

!     METHOD:
!     -------
!     A DIMENSIONLESS RELATION IN THE FORM
!          E10 = A * FR10 **B,                            (1)
!     IT IS ASSUMED WITH
!          E10=E G**2/USTAR**4  (DIMENSIONLESS ENERGY)     (2)
!     AND
!          FR10=FR * USTAR/G (DIMENSIONLESS FREQUENCY).    (3)
!     USTAR IS THE FRICTION VELOCITY, 
!     DIMENSIONLESS FREQUENCY based on the 1/f F(f) formulation

!     SUBSTITUITING E10 AND FR10 IN (1) ACCORDING TO (2) AND (3),
!     ONE HAS AN EQUATION FOR THE UNKNOWN USTAR. WITH GIVES

!      USTAR= XX**(-1/(4+B)) ; XX = A F**B / (E G**(B+2).

!     THE CONSTANT A AND B HAVE BEEN OBTAINED FROM THE WAM TIME
!     LIMITED GROWTH CURVE. SEE WORK DONE FOR THE ALTIMETER 
!     ASSIMILATION.

! ----------------------------------------------------------------------

      USE YOWPCONS , ONLY : G        ,AFMRV    ,BFMRV
      USE YOWTABL  , ONLY : USTARM

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL A,B,B2,B4

      REAL 
     &   FPTMP(NSPECW,MPART),
     &   EPTMP(NSPECW,MPART),
     &   USE(NSPECW),
     &   XX

      INTEGER
     &   NINTW(NSPECW),
     &   NSPECW,
     &   MPART,
     &   ISPEC,
     &   JPART

! ----------------------------------------------------------------------

      A = AFMRV
      B = BFMRV
      B2=-(B+2.)
      B4=-1./(B+4.)

      DO ISPEC=1,NSPECW
        IF(NINTW(ISPEC).NE.0) THEN
          JPART=NINTW(ISPEC)
  
          IF(EPTMP(ISPEC,JPART).GT.0.) THEN

!*          USTAR IS DERIVED FROM THE WINDSEA E-FM CURVE.
!           -------------------------------------------

            XX=A*G**B2*FPTMP(ISPEC,JPART)**B/EPTMP(ISPEC,JPART)
            USE(ISPEC)=XX**B4

            USE(ISPEC)=MIN(USE(ISPEC),USTARM)
          ELSE
            USE(ISPEC)=-1.
          ENDIF
        ENDIF
      ENDDO

      RETURN

      END SUBROUTINE U10COR
