      REAL FUNCTION AKI_ICE (G,XK,DEPTH,RHOW,CITH)

! ----------------------------------------------------------------------

!**** *AKI_ICE* - FUNCTION TO COMPUTE WAVE NUMBER UNDER THE ICE.
!                 FOR A GIVEN OPEN WATER WAVE NUMBER.

!*    PURPOSE.
!     -------

!       *AKI_ICE* COMPUTES THE REAL WAVE NUMBER UNDER THE ICE AS FUNCTION OF
!                 OPEN OCEAN WAVE NUMBER, THE WATER DEPTH AND THE ICE THICKNESS.

!**   INTERFACE.
!     ----------

!       *FUNCTION* *AKI_ICE (G,XK,DEPTH,RHOW,CITH))*
!          *G*        - ACCELERATION OF GRAVITY (m/s**2).
!          *XK*       - OPEN OCEAN WAVE NUMBER (1/m).
!          *DEPTH*    - WATER DEPTH (m).
!          *RHOW*     - WATER DENSITY (kg/m**3).
!          *CITH*     - SEA ICE THICKNESS (m).

!     METHOD.
!     -------

!       NEWTONS METHOD TO SOLVE THE LINEAR DISPERSION RELATION IN 
!       SHALLOW WATER UNDER AN INFINITELY LONG ELASTIC FLOATING PLATE
!       REPRESENTING THE SEA ICE. THE MECHANICAL PROPERTIES OF THE SEA
!       ICE IS GIVEn BY ITS YOUNG MODULUS, THE POISSON'S RATIO AND ITS
!       DENSITY (these are fixed for now, see below). 

!       IF F(x)=0, then solve iteratively
!       x(n+1) = x(n) - F(x(n))/F'(x(n))

!       WHERE F'(x) IS THE FIRST DERIVATIVE OF F WITH RESPECT TO x.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOX AND SQUIRE, 1991, JGR 96, C3, 4531-4547.

!       NONE.

! ----------------------------------------------------------------------
      IMPLICIT NONE

!     ICE PROPERTIES (assumed fixed for now) 
      REAL, PARAMETER :: YMICE=5.5E+9   ! typical value of Young modulus of sea ice
      REAL, PARAMETER :: RMUICE=0.3      ! Poisson's ratio of sea ice
      REAL, PARAMETER :: RHOI=922.5     ! typical value of the sea ice density

!     RELATIVE ERROR LIMIT OF NEWTONS METHOD.
      REAL, PARAMETER :: EBS = 0.000001
!     MAXIMUM WAVE NUMBER
      REAL, PARAMETER :: AKI_MAX = 20.

      REAL :: G, XK, DEPTH, RHOW, CITH
      REAL :: FICSTF, RDH
      REAL :: OM2, AKI, AKIOLD, F, FPRIME, AKID 


      IF(CITH.LE.0.0) THEN
        AKI=XK
      ELSE
!       BENDING STIFFNESS / WATER DENSITY
        FICSTF=(YMICE*CITH**3/(12*(1-RMUICE**2)))/RHOW

!       DENSITY RATIO * ICE THICKNESS
        RDH=(RHOI/RHOW)*CITH

!       SQUARE OF THE OPEN OCEAN ANGULAR FREQUENCY
        OM2=G*XK*TANH(XK*DEPTH)

!*      2. ITERATION LOOP.
!          ---------------

        AKIOLD=0.
        AKI=MIN(XK ,(OM2/MAX(FICSTF,1.))**0.2)

        DO WHILE (ABS(AKI-AKIOLD).GT.EBS*AKIOLD .AND.
     &            AKI.LT.AKI_MAX
     &           )
          AKIOLD=AKI 
          AKID=MIN(DEPTH*AKI,50.)
          F=FICSTF*AKI**5+G*AKI-OM2*(RDH*AKI+1./TANH(AKID))
          FPRIME=5.*FICSTF*AKI**4+G-OM2*(RDH-DEPTH/(SINH(AKID)**2))
          AKI = AKI-F/FPRIME
!         in case of overshoot because it is trying to find a very large wave number
          IF(AKI.LE.0.) AKI=AKI_MAX
        ENDDO

      ENDIF

      AKI_ICE = AKI

      RETURN
      END
