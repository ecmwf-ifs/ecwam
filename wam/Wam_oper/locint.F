      SUBROUTINE LOCINT (IG, IJS, IJL, NC, NR, KCOL, KROW, IGPER,
     &                  DLAM, DPHI, RLONL, RLATS, U, V, US, DS)

! ----------------------------------------------------------------------

!**** *LOCINT* - LOCATE WAM POINT IN INPUT GRID AND INTERPOLATE.

!     H. GUNTHER    ECMWF/GKSS  DECEMBER 1990  MODIFIED FOR CYCLE_4.

!*    PURPOSE.
!     --------

!        *LOCINT* - LOCATE AND INTERPOLATE IN INPUT GRID.

!**   INTERFACE.
!     ----------

!       *CALL* *LOCINT (IG, IJS, IJL, NC, NR, KCOL, KROW, IGPER,
!                       DLAM, DPHI, RLONL, RLATS, U, V, US, DS)*
!         *IG*      - BLOCK NUMBER.
!         *IJS*     - FIRST POINT IN BLOCK.
!         *IJL*     - LAST POINT IN BLOCK.
!         *NC*      - DIMENSION OF INPUT ARRAY, NUMBER OF COLUMNES.
!         *NR*      - DIMENSION OF INPUT ARRAY, NUMBER OF ROWS.
!         *KCOL*    - NUMBER OF COLUMNES IN INPUT (USED).
!         *KROW*    - NUMBER OF ROWS     IN INPUT (USED).
!         *IGPER*   - INDICATOR PERIODICAL INPUT GRID OR NOT
!                     1 = PERIODICAL
!                     OTHERWISE NON-PERIODICAL
!         *DLAM*    - STEPSIZE BETWEEN LONGITUDES IN INPUT (DEG).
!         *DPHI*    - STEPSIZE BETWEEN LATITUDES  IN INPUT (DEG).
!         *RLATS*   - MOST SOUTHERN LATITUDE OF INPUT (DEG).
!         *RLONL*   - MOST WESTERN LONGITUDE OF INPUT (DEG).
!         *U(I,J)*  - "HORIZONTAL" COMPONENT (INPUT).
!         *V(I,J)*  - "VERTICAL  " COMPONENT (INPUT).
!         *US*      - SPACE INTERPOLATED "HORIZONTAL" COMPONENT.
!         *DS*      - SPACE INTERPOLATED "VERTICAL  " COMPONENT.

!     METHOD.
!     -------

!       DOUBLE LINEAR INTERPOLATION WITHIN A MESH FOR U AND V.

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,AMOWEP   ,AMOSOP   ,
     &            XDELLA   ,ZDELLO
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      REAL,DIMENSION(IJS:IJL)  :: US,DS

      DIMENSION U(NC,NR), V(NC,NR)

! ----------------------------------------------------------------------

!*    1. LOOP OVER POINTS IN WAMODEL BLOCKS.
!        -----------------------------------

      IF(DLAM == 0.) THEN
        DLAMM1 = 1.
      ELSE
        DLAMM1 = 1./DLAM
      ENDIF

      IF(DPHI == 0.) THEN
        DPHIM1 = 1.
      ELSE
        DPHIM1 = 1./DPHI
      ENDIF

      DO IJ=IJS,IJL

!*    1.1 TRANSFORM WAM COORDINATE TO INPUT GRID.
!         ---------------------------------------

        KX = KXLT(IJ,IG)
        XI = AMOWEP + REAL(IXLG(IJ,IG)-1)*ZDELLO(KX)-RLONL
        XI = MOD(XI+720.,360.)
        XI = XI*DLAMM1+1.00001
        XJ = AMOSOP + REAL(KXLT(IJ,IG)-1)*XDELLA-RLATS
        XJ = XJ*DPHIM1+1.00001

!*    1.2 COMPUTE CORNER POINT INDICES OF INPUT GRID.
!         -------------------------------------------

        I1  = XI
        J1  = XJ
        J2  = MIN(KROW,J1+1)
        I2  = I1+1

!*    1.3 DISTANCES OF INTERPOLATION POINT FROM CORNER POINTS.
!         ----------------------------------------------------

        DI1 = XI-I1
        DI2 = 1.-DI1
        DJ1 = XJ-J1
        DJ2 = 1.-DJ1

!*    1.4. CORRECTIONS FOR FIRST AND LAST GRID LINES
!          PERIODIC OR UNPERIODIC GRIDS ARE CONSIDERED.
!          --------------------------------------------

        IF (IGPER.EQ.1) THEN
          IF (I1.EQ.KCOL) I2 = 1
          IF (I1.EQ.0   ) I1 = KCOL
        ELSE
          IF (I1.EQ.KCOL) I2 = KCOL
        ENDIF

!*    1.5 CHECK WHETHER POINTS ARE IN GRID.
!         ---------------------------------

        IF (I1.LT.1 .OR. I1.GT.KCOL .OR. J1.LT.1 .OR. J1.GT.KROW) THEN
          WRITE(IU06,*) ' *******************************************'
          WRITE(IU06,*) ' *                                         *'
          WRITE(IU06,*) ' *          FATAL ERROR IN LOCINT          *'
          WRITE(IU06,*) ' *          =====================          *'
          WRITE(IU06,*) ' * POINT IS OUTSIDE OF INPUT GRID          *'
          WRITE(IU06,*) ' * DIMENSION OF INPUT GRID IS KCOL = ', KCOL
          WRITE(IU06,*) ' *                            KROW = ', KROW
          WRITE(IU06,*) ' * POINT REQUESTED IS         XI   = ', XI
          WRITE(IU06,*) ' *                            XJ   = ', XJ
          WRITE(IU06,*) ' *                                         *'
          WRITE(IU06,*) ' *  PROGRAM ABORTS     PROGRAM ABORTS      *'
          WRITE(IU06,*) ' *                                         *'
          WRITE(IU06,*) ' *******************************************'
          CALL ABORT1
        ENDIF

!*    1.6 LINEAR INTERPOLATION.
!         ---------------------

        U1 = U(I1,J1)*DI2 + U(I2,J1)*DI1
        U2 = U(I1,J2)*DI2 + U(I2,J2)*DI1
        US(IJ) = U1*DJ2 + U2*DJ1
        V1 = V(I1,J1)*DI2 + V(I2,J1)*DI1
        V2 = V(I1,J2)*DI2 + V(I2,J2)*DI1
        DS(IJ) = V1*DJ2 + V2*DJ1
      ENDDO

      RETURN
      END SUBROUTINE LOCINT
