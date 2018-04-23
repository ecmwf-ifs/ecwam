!--------------------------------------------------------------------

      SUBROUTINE SKEWNESS(IU06,F1,NCOLL,XKAPPA1,DELH_ALT)

!--------------------------------------------------------------------

!*****SKEWNESS** COMPUTES PARAMETERS OF THE NEARLY-GAUSSIAN
!             DISTRIBUTION OF OCEAN WAVES AT A FIXED GRID POINT.

!     P.JANSSEN JULY 1997

!     PURPOSE
!     -------
!             DETERMINES SKEWNESS PARAMETERS IN ORDER TO OBTAIN
!             CORRECTION ON ALTIMETER WAVE HEIGHT.

!     INTERFACE
!     ---------
!             *CALL* *SKEWNESS(IU06,F1,NCOLL,XKAPPA1,DELH_ALT)*


!     PARAMETER   TYPE      PURPOSE.
!     ---------   ----      -------
!
!       IU06      INTEGER   PRINTER OUTPUT UNIT.
!       F1        REAL      TWO DIMENSIONAL SPECTRUM
!       NCOLL     INTEGER   NUMBER OF COLLOCATED SPECTRA
!       XKAPPA1   REAL      CORRECTED KAPPA1 FROM ALTIMETER WAVE HEIGHT
!                           ALGORITHM
!       DELH_ALT  REAL      RELATIVE ALTIMETER RANGE CORRECTION,
!                           I.E. DELH_ALT*HS GIVES ACTUAL RANGE
!                           CORRECTION

!     METHOD
!     ------
!             EVALUATE DEVIATIONS FROM GAUSSIANITY FOLLOWING THE WORK
!             OF SROKOSZ AND LONGUET-HIGGINS. FOR SECOND ORDER
!             CORRECTIONS TO SURFACE ELEVATION THE APPROACH OF
!             ZAKHAROV HAS BEEN USED.

!     EXTERNALS
!     ---------
!             NONE

!     REFERENCES
!     ----------
!             M.A. SROKOSZ, J.G.R.,91,995-1006(1986)
!             V.E. ZAKHAROV, HAMILTONIAN APPROACH(1967)

!--------------------------------------------------------------------

      USE YOWFRED  , ONLY : TH       ,COSTH    ,SINTH 
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWTABL  , ONLY : NFREHF   ,FAC0     ,FAC1     ,FAC2     ,
     &            FAC3     ,FAK      ,FRHF     ,DFIMHF

!--------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IU06, NCOLL

      REAL, DIMENSION(NCOLL), INTENT(OUT) :: XKAPPA1, DELH_ALT
      REAL, DIMENSION(NCOLL,NANG,NFRE), INTENT(IN) :: F1


      INTEGER :: M, K, IJ, M1, K1, M2, K2, I, J
      INTEGER :: MSTART
   
      REAL :: FH, DELF, XK1
      REAL :: XPI, XPJ, XPK, XN, DELTA, XFAC
      REAL, DIMENSION(NCOLL,NANG,NFREHF) :: F2
      REAL, DIMENSION(NCOLL,0:3,0:2,0:2) :: XMU, XLAMBDA


!     1. COMPUTATION OF FREQUENCY-DIRECTION INCREMENT
!     -----------------------------------------------

      MSTART = 1

      XMU(:,:,:,:) = 0.

      DO K=1,NANG
        DO M=1,NFRE
          DO IJ=1,NCOLL
            F2(IJ,K,M)=F1(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO

      DO M=NFRE+1,NFREHF
        FH=(FRHF(NFRE)/FRHF(M))**5
        DO K=1,NANG
          DO IJ=1,NCOLL
            F2(IJ,K,M)=F1(IJ,K,NFRE)*FH
          ENDDO
        ENDDO
      ENDDO

!     2. COMPUTATION OF THE SKEWNESS COEFFICIENTS
!     --------------------------------------------

      DO M1=MSTART,NFREHF
        DO M2=MSTART,NFREHF
          DO K1=1,NANG
            DO K2=1,NANG
              DO IJ=1,NCOLL
                DELF = DFIMHF(M1)*DFIMHF(M2)*F2(IJ,K1,M1)*F2(IJ,K2,M2)
                XMU(IJ,3,0,0) = XMU(IJ,3,0,0)+3.*FAC0(K1,K2,M1,M2)*DELF
                XMU(IJ,1,2,0) = XMU(IJ,1,2,0)+FAC1(K1,K2,M1,M2)*DELF
                XMU(IJ,1,0,2) = XMU(IJ,1,0,2)+FAC2(K1,K2,M1,M2)*DELF
                XMU(IJ,1,1,1) = XMU(IJ,1,1,1)+FAC3(K1,K2,M1,M2)*DELF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO K1=1,NANG
        DO M1=MSTART,NFREHF
          XK1 = FAK(M1)**2
          DO IJ=1,NCOLL
            DELF = DFIMHF(M1)*F2(IJ,K1,M1)
            XMU(IJ,2,0,0) = XMU(IJ,2,0,0) + DELF
            XMU(IJ,0,2,0) = XMU(IJ,0,2,0) + XK1*COSTH(K1)**2*DELF
            XMU(IJ,0,0,2) = XMU(IJ,0,0,2) + XK1*SINTH(K1)**2*DELF
            XMU(IJ,0,1,1) = XMU(IJ,0,1,1) +
     &                      XK1*COSTH(K1)*SINTH(K1)*DELF
          ENDDO
        ENDDO
      ENDDO


!     3. COMPUTATION OF THE NORMALISED SKEWNESS COEFFICIENTS
!     ------------------------------------------------------

      DO I=0,3
        XPI = 0.5*FLOAT(I)
        DO J=0,2
          XPJ = 0.5*FLOAT(J)
          DO K=0,2
            XPK = 0.5*FLOAT(K)
            DO IJ=1,NCOLL
              XN = XMU(IJ,2,0,0)**XPI*XMU(IJ,0,2,0)**XPJ
     &            *XMU(IJ,0,0,2)**XPK
              XLAMBDA(IJ,I,J,K) = XMU(IJ,I,J,K)/XN
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!     4. CORRECTION TO KAPPA1
!     -----------------------

      DO IJ=1,NCOLL
         DELTA = ( XLAMBDA(IJ,1,2,0) + XLAMBDA(IJ,1,0,2)
     &             - 2.*XLAMBDA(IJ,0,1,1)*XLAMBDA(IJ,1,1,1) )/
     &             (1. - XLAMBDA(IJ,0,1,1)**2)
         XFAC  = 2.*(XLAMBDA(IJ,3,0,0)/3. + DELTA)*
     &           (5.*XLAMBDA(IJ,3,0,0)/24. + 0.125*DELTA)
         XKAPPA1(IJ) = 1. + XFAC
         DELH_ALT(IJ) = -0.125*(XLAMBDA(IJ,3,0,0)/3.+DELTA)
      ENDDO

!--------------------------------------------------------------------

      RETURN
      END SUBROUTINE SKEWNESS
