!--------------------------------------------------------------------

      SUBROUTINE SKEWNESS(IU06, F1,NCOLL,NANG,NFRE,FR1,TH,COSTH,SINTH,
     1                   XKAPPA1,DELH_ALT)

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
!             *CALL* *SKEWNESS(IU06,F1,NCOLL,NANG,NFRE,FR1,TH,COSTH,
!                             SINTH,XKAPPA1,DELH_ALT)*


!     PARAMETER   TYPE      PURPOSE.
!     ---------   ----      -------
!
!       IU06      INTEGER   PRINTER OUTPUT UNIT.
!       F1        REAL      TWO DIMENSIONAL SPECTRUM
!       NCOLL     INTEGER   NUMBER OF COLLOCATED SPECTRA
!       NANG      INTEGER   NUMBER OF DIRECTIONS
!       NFRE      INTEGER   NUMBER OF FREQUENCIES
!       FR1       REAL      FIRST FREQUENCY
!       TH        REAL      DIRECTION (RADIANS) ARRAY
!       COSTH     REAL      COS(TH)
!       SINTH     REAL      SIN(TH)
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

      USE YOWTABL  , ONLY : NFREHF   ,A        ,B        ,FAC1     ,
     &            FAC2     ,FAC3

!--------------------------------------------------------------------

      DIMENSION F1(NCOLL,NANG,NFRE),TH(NANG), COSTH(NANG),
     &          SINTH(NANG), XKAPPA1(NCOLL), DELH_ALT(NCOLL)

      DIMENSION DF(NFREHF),F2(NCOLL,NANG,NFREHF),
     &          XMU(NCOLL,0:3,0:2,0:2),XLAMBDA(NCOLL,0:3,0:2,0:2)
      DIMENSION FR(NFREHF),FK(NFREHF)
      LOGICAL FRSTIME
      DATA FRSTIME/.TRUE./


!     1. COMPUTATION OF FREQUENCY-DIRECTION INCREMENT
!     -----------------------------------------------

      PI  = 4.*ATAN(1.)
      ZPI = 2.*PI
      G   = 9.806

      MSTART = 1
      DELTH  = TH(2)-TH(1)
      FR(1)  = FR1

      DO M=2,NFREHF
        FR(M) = 1.1*FR(M-1)
      ENDDO

      CO1   = 0.1*DELTH/2.
      DF(1) = CO1*FR(1)

      DO M=2,NFREHF-1
         DF(M)=CO1*(FR(M)+FR(M-1))
      ENDDO
      DF(NFREHF)=CO1*FR(NFREHF-1)

      DO I=0,3
        DO J=0,2
          DO K=0,2
            DO IJ=1,NCOLL
              XMU(IJ,I,J,K) = 0.
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO K=1,NANG
        DO M=1,NFRE
          DO IJ=1,NCOLL
            F2(IJ,K,M)=F1(IJ,K,M)
          ENDDO
        ENDDO
      ENDDO

      DO M=NFRE,NFREHF
        FH=(FR(NFRE)/FR(M))**5
        DO K=1,NANG
          DO IJ=1,NCOLL
            FLH=F1(IJ,K,NFRE)
            F2(IJ,K,M)=FLH*FH
          ENDDO
        ENDDO
      ENDDO

      DO M=1,NFREHF
         FK(M) = (ZPI*FR(M))**2/G
      ENDDO

!     IF THEN DETERMINE SECOND ORDER COEFFICIENTS

      IF (FRSTIME) THEN

        ALLOCATE(A(NANG,NANG,NFREHF,NFREHF))
        ALLOCATE(B(NANG,NANG,NFREHF,NFREHF))
        ALLOCATE(FAC1(NANG,NANG,NFREHF,NFREHF))
        ALLOCATE(FAC2(NANG,NANG,NFREHF,NFREHF))
        ALLOCATE(FAC3(NANG,NANG,NFREHF,NFREHF))

        CALL SECONDHH(NFREHF,NANG,TH,COSTH,SINTH,FK,A,B,FAC1,
     &                 FAC2,FAC3)
        PRINT*,' SECOND ORDER COEFFICIENTS DETERMINED '
        FRSTIME = .FALSE.
      END IF


!     2. COMPUTATION OF THE SKEWNESS COEFFICIENTS
!     --------------------------------------------

      DO M1=MSTART,NFREHF
        DO M2=MSTART,NFREHF
          XK1 = FK(M1)
          XK2 = FK(M2)
          DO K1=1,NANG
            DO K2=1,NANG
              C22 = A(K1,K2,M1,M2)
              S22 = B(K1,K2,M1,M2)
              FC1 = FAC1(K1,K2,M1,M2)
              FC2 = FAC2(K1,K2,M1,M2)
              FC3 = FAC3(K1,K2,M1,M2)
              DO IJ=1,NCOLL
                DELF = DF(M1)*DF(M2)*F2(IJ,K1,M1)*F2(IJ,K2,M2)
                XMU(IJ,3,0,0) = XMU(IJ,3,0,0) + 3.*C22*DELF

                XMU(IJ,1,2,0) = XMU(IJ,1,2,0) + FC1*DELF

                XMU(IJ,1,0,2) = XMU(IJ,1,0,2) + FC2*DELF

                XMU(IJ,1,1,1) = XMU(IJ,1,1,1) + FC3*DELF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO K1=1,NANG
        COS1 = COSTH(K1)
        SIN1 = SINTH(K1)
        COS12 = COS1**2
        SIN12 = SIN1**2
        DO M1=MSTART,NFREHF
          XK1 = FK(M1)
          DO IJ=1,NCOLL
            DELF = DF(M1)*F2(IJ,K1,M1)
            XMU(IJ,2,0,0) = XMU(IJ,2,0,0) + DELF
            XMU(IJ,0,2,0) = XMU(IJ,0,2,0) + XK1**2*COS12*DELF
            XMU(IJ,0,0,2) = XMU(IJ,0,0,2) + XK1**2*SIN12*DELF
            XMU(IJ,0,1,1) = XMU(IJ,0,1,1) +
     &                      XK1**2*COS1*SIN1*DELF
          ENDDO
        ENDDO
      ENDDO

!     3. COMPUTATION OF THE NORMALISED SKEWNESS COEFFICIENTS
!     ------------------------------------------------------

      DO I=0,3
        DO J=0,2
          DO K=0,2
            XPI = 0.5*FLOAT(I)
            XPJ = 0.5*FLOAT(J)
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
     &           (5.*XLAMBDA(IJ,3,0,0)/24. + DELTA/8.)
         XKAPPA1(IJ) = 1. + XFAC
         DELH_ALT(IJ) = -1./8.*(XLAMBDA(IJ,3,0,0)/3.+DELTA)
      ENDDO

!--------------------------------------------------------------------

      RETURN
      END SUBROUTINE SKEWNESS
