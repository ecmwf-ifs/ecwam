      SUBROUTINE WAMWND (US, DS, ADS, ZIDL, KINF, KSUP,IJS, IJL, IG)

! ----------------------------------------------------------------------

!     MODIFIED  S. ABDALLA OCTOBER 2001  INCLUSION OF AIR DENSITY & Zi/L

!**** *WAMWND* - TRANSFORMS INPUT WINDS TO BLOCKED WAM POINTS.


!*    PURPOSE.
!     --------

!       INTERPOLATE AND CONVERT INPUT WINDS TO WAM WINDS FOR ALL
!       POINTS IN THE GRID.

!**   INTERFACE.
!     ----------

!       *CALL WAMWND (US, DS, ADS, ZIDL, KINF, KSUP,IJS, IJL, IG)*
!          *US*   - INTERPOLATED WINDS AT ALL POINTS AND BLOCKS.
!          *DS*   - INTERPOLATED WIND DIRECTION AT ALL POINTS.
!          *ADS*  - INTERPOLATED AIR DENSITY AT ALL POINTS.
!          *ZIDL* - INTERPOLATED Zi/L AT ALL POINTS
!                   (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!          *KINF* - LOWER BOUND ON FIRST DIMENSION OF FIELD
!          *KSUP* - UPPER BOUND ON FIRST DIMENSION OF FIELD
!          *IJS*    - INDEX OF FIRST GRIDPOINT
!          *IJL*    - INDEX OF LAST GRIDPOINT
!          *IG*     - BLOCK NUMBER

!     METHOD.
!     -------

!       THE INPUT WINDS ARE INTERPOLATED TO THE WAVE MODEL GRID POINTS
!       IN BLOCKED FORMAT. THE INTERPOLATED VALUES ARE TRANSFORMED TO
!       MAGNITUDE AND DIRECTION. INPUT MAY BE WIND IN 10M HEIGHT ,
!       SURFACE WINDS OR FRICTION VELOCITIES. THE INPUT GRID HAS TO BE
!       A LATITUDE/LONGITUDE GRID EITHER PERIODIC OR NON PERIODIC.

!     EXTERNALS.
!     ----------

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NGY      ,NBLO
      USE YOWPCONS , ONLY : G        ,ZPI      ,ROAIR
      USE YOWCOUP  , ONLY : XKAPPA
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMAP   , ONLY : IXLG     ,KXLT
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWWNDG  , ONLY : DLAM     ,DPHI     ,RLATS    ,RLONL    ,
     &            KCOL     ,KROW     ,ICODE    ,IWPER
      USE YOWWIND  , ONLY : WSPMIN   ,NC       ,NR       ,UWND     ,
     &            VWND     ,AIRD     ,PZIDL

! ----------------------------------------------------------------------

      REAL, PARAMETER :: ALPHACH = 0.0185

! ----------------------------------------------------------------------

      REAL, DIMENSION (KINF:KSUP,NBLO) :: US,DS,ADS,ZIDL

! ----------------------------------------------------------------------

!*    2. TRANSFORM GRIDDED WIND INPUT INTO BLOCK
!*       AND TRANSFORM TO MAGNITUDE AND DIRECTION.
!*       AND SET LOWER WIND SPEED TO WSPMIN m/s
!        ----------------------------------------

      DO IJ = IJS,IJL
        IX = IXLG(IJ,IG)
        JY = NGY-KXLT(IJ,IG) + 1

        UU = UWND(IX,JY)
        VV = VWND(IX,JY)
        US(IJ,IG) = SQRT(UU**2 + VV**2)
        IF (US(IJ,IG).NE.0.) THEN 
          DS(IJ,IG) = ATAN2(UU,VV)
        ELSE
          DS(IJ,IG) = 0.
        ENDIF
        IF (DS(IJ,IG).LT.0.) DS(IJ,IG) = DS(IJ,IG) + ZPI
        US(IJ,IG) = MAX(US(IJ,IG),WSPMIN)

        ADS(IJ,IG) = AIRD(IX,JY)
        ZIDL(IJ,IG) = PZIDL(IX,JY)
      ENDDO

! ----------------------------------------------------------------------

!*    3. PROCESS WINDS ACCORDING TO TYPE
!        NOTHING TO DO FOR WIND SPEED U10 (ICODE = 3).
!        ---------------------------------------------

      IF (ICODE.EQ.1) THEN

!*    3.2  INPUT IS FRICTION VELOCITY.
!          ---------------------------

        DO IJ = IJS,IJL
          USTAR = MAX(0.01,US(IJ,IG))
          Z0  = ALPHACH/G*USTAR**2
          CD  = XKAPPA/LOG(10./Z0)
          US(IJ,IG) = USTAR/CD
          US(IJ,IG) = MAX(US(IJ,IG),WSPMIN)
        ENDDO

      ELSEIF (ICODE.EQ.2) THEN

!*    3.3 INPUT WINDS ARE SURFACE STRESSES.
!         ---------------------------------

        DO IJ = IJS,IJL
          USTAR = MAX (0.01, SQRT(US(IJ,IG)/ADS(IJ,IG)))
          Z0  = ALPHACH/G*USTAR**2
          CD  = XKAPPA/LOG(10./Z0)
          US(IJ,IG) = USTAR/CD
          US(IJ,IG) = MAX(US(IJ,IG),WSPMIN)
        ENDDO
      ENDIF

!*    4. TEST OUTPUT OF WAVE MODEL BLOCKS
!        ---------------------------------

      IF (ITEST.GE.3) THEN
        WRITE (IU06,*) ' '
        WRITE (IU06,*) '      SUB. WAMWND:',
     &   ' WINDFIELDS CONVERTED TO BLOCKS'
        WRITE (IU06,*) ' '
        IJA = IJS
        IJE = MIN(IJS+4,IJL)
        WRITE (IU06,*) ' BLOCK= ',IG
        WRITE (IU06,*) ' '
        WRITE (IU06,*) ' US (IJS - IJS+4): ',(US(I,IG),I=IJA,IJE)
        WRITE (IU06,*) ' DS (IJS - IJS+4): ',(DS(I,IG),I=IJA,IJE)
        WRITE (IU06,*) ' ADS(IJS - IJS+4): ',(ADS(I,IG),I=IJA,IJE)
        WRITE (IU06,*) ' ZIDL(IJS - IJS+4): ',(ZIDL(I,IG),I=IJA,IJE)
        CALL FLUSH(IU06)
      ENDIF

      RETURN
      END SUBROUTINE WAMWND
