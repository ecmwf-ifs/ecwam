      SUBROUTINE U10COR (U10E, FPTMP, EPTMP, NSPECW, MPART, NINTW, UW10,
     &                   WSPMIN)
!     PURPOSE.
!     --------
!     COMPUTES RELATION BETWEEN U10 AND ENERGY , MEAN FREQUENCY.

!     EXTERNALS.
!     ----------

!     METHOD:
!     -------
!     A DIMENSIONLESS RELATION IN THE FORM
!          E10= A * FR10 **B,                            (1)
!     IT IS ASSUMED WITH
!          E10=E G**2/U10**4  (DIMENSIONLESS ENERGY)     (2)
!     AND
!          FR10=FR * U10/G (DIMENSIONLESS FREQUENCY).    (3)
!     U10 IS THE WIND SPEED AT 10M., I.E. THE INPUT OF WAM.

!     SUBSTITUITING E10 AND FR10 IN (1) ACCORDING TO (2) AND (3),
!     ONE HAS AN EQUATION FOR THE UNKNOWN U10. WITH GIVES

!      U10= XX**(-1/(4+B)) ; XX = A F**B / (E G**(B+2).

!     THE CONSTANT A AND B HAVE BEEN OBTAINED FROM THE WAM TIME
!     LIMITED GROWTH CURVE. A IS AN EVALUATION OBTAINED AT
!     20M/S, AM AT 10M/S. THE WIND SPEED IS ESTIMATED FROM A
!     WEIGHTED AVERAGE OF THE TWO EVALUATION.

! ----------------------------------------------------------------------

      IMPLICIT NONE

!     LOCAL PARAMETERS

      REAL A,B,B2,B4,AM,G

      PARAMETER( A=0.0000225, B=-2.694, B2=-(B+2.), B4=-1./(B+4.),
     &           AM= 0.0000176, G=9.806 )

      REAL 
     &   FPTMP(NSPECW,MPART),
     &   EPTMP(NSPECW,MPART),
     &   U10E(NSPECW),
     &   UW10(NSPECW),
     &   U10,U10M,XX,WD,WDM,WT,WSPMIN

      INTEGER
     &   NINTW(NSPECW),
     &   NSPECW,
     &   MPART,
     &   ISPEC,
     &   JPART

! ----------------------------------------------------------------------

      DO ISPEC=1,NSPECW
        IF(NINTW(ISPEC).NE.0) THEN
          JPART=NINTW(ISPEC)
  
          IF(EPTMP(ISPEC,JPART).GT.0.) THEN

!*          U10 IS DERIVED FROM THE WINDSEA E-FM CURVE.
!           -------------------------------------------

            XX=A*G**B2*FPTMP(ISPEC,JPART)**B/EPTMP(ISPEC,JPART)
            U10=XX**B4

!*          AUXILIARY EVALUATION.
!           ---------------------

            XX=AM*G**B2*FPTMP(ISPEC,JPART)**B/EPTMP(ISPEC,JPART)
            U10M=XX**B4

!*          WEIGHTED AVERAGE.
!           -----------------

            WD=MAX( ABS(U10-20.),0.001)**(-2)
            WDM=MAX( ABS(U10M-10.),0.001)**(-2)
            WT=WD+WDM
            U10E(ISPEC)=(WD*U10+WDM*U10M)/WT
          ELSE
!           in case the assimilation scheme has reduced the windsea
!           contribution to 0. set wind speed to WSPMIN 
            U10E(ISPEC)=WSPMIN
          ENDIF
        ELSE
          U10E(ISPEC)=UW10(ISPEC)
        ENDIF
      ENDDO

      RETURN 

      END SUBROUTINE U10COR
