      SUBROUTINE WAMWND (US, DS, KINF, KSUP , KEND)

! ----------------------------------------------------------------------

!**** *WAMWND* - TRANSFORMS INPUT WINDS TO BLOCKED WAM POINTS.


!*    PURPOSE.
!     --------

!       INTERPOLATE AND CONVERT INPUT WINDS TO WAM WINDS FOR ALL
!       POINTS IN THE GRID.

!**   INTERFACE.
!     ----------

!       *CALL WAMWND (US, DS, KINF, KSUP , KEND)*
!          *US*   - INTERPOLATED WINDS IN ALL POINTS AND BLOCKS.
!          *DS*   - INTERPOLATED WIND DIRECTION IN ALL POINTS.
!          *KINF* - LOWER BOUND ON FIRST DIMENSION OF FIELD
!          *KSUP* - UPPER BOUND ON FIRST DIMENSION OF FIELD
!          *KEND* - USED UPPER LIMIT ON FIRST DIMENSION OF FIELD

!     METHOD.
!     -------

!       THE INPUT WINDS ARE INTERPOLATED TO THE WAVE MODEL GRID POINTS
!       IN BLOCKED FORMAT. THE INTERPOLATED VALUES ARE TRANSFORMED TO
!       MAGNITUDE AND DIRECTION. INPUT MAY BE WIND IN 10M HEIGHT ,
!       SURFACE WINDS OR FRICTION VELOCITIES. THE INPUT GRID HAS TO BE
!       A LATITUDE/LONGITUDE GRID EITHER PERIODIC OR NON PERIODIC.

!     EXTERNALS.
!     ----------

!       *LOCINT*    - LOCATES WAM POINT IN INPUT GRID AND INTERPOLATES.

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NGY      ,NBLO     ,NIBLO
      USE YOWPCONS , ONLY : G        ,ZPI      ,ROAIR
      USE YOWCOUP  , ONLY : XKAPPA
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMAP   , ONLY : IXLG     ,KXLT
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWWNDG  , ONLY : DLAM     ,DPHI     ,RLATS    ,RLONL    ,
     &            KCOL     ,KROW     ,ICODE    ,IWPER
      USE YOWWIND  , ONLY : WSPMIN   ,NC       ,NR       ,UWND     ,
     &            VWND

! ----------------------------------------------------------------------

      REAL, PARAMETER :: ALPHACH = 0.0185

! ----------------------------------------------------------------------

      REAL, DIMENSION (KINF:KSUP,NBLO) :: US,DS
      LOGICAL LFRSTIME

      DATA LFRSTIME /.TRUE./

! ----------------------------------------------------------------------

!*    2. TRANSFORM GRIDDED WIND INPUT INTO BLOCK
!        ----------------------------------------

      DO IG=1,IGL
        IF(KINF.EQ.1.AND.KSUP.EQ.NIBLO.AND.LMESSPASS) THEN
          DO IJ = KINF,KSUP
            IX = IXLG(IJ,IG)
            JY = NGY-KXLT(IJ,IG) + 1
            US(IJ,IG) = UWND(IX,JY)
            DS(IJ,IG) = VWND(IX,JY)
          ENDDO
        ELSE
          DO IJ = IJS(IG),IJL(IG)
            IX = IXLG(IJ,IG)
            JY = NGY-KXLT(IJ,IG) + 1
            US(IJ,IG) = UWND(IX,JY)
            DS(IJ,IG) = VWND(IX,JY)
          ENDDO
        ENDIF
      ENDDO

! ----------------------------------------------------------------------

!*    3. PROCESS WINDS ACCORDING TO TYPE
!        NOTHING TO DO FOR WIND SPEED U10 (ICODE = 3).
!        ---------------------------------------------


!*    3.1 TRANSFORM TO MAGNITUDE AND DIRECTION.
!         -------------------------------------

      DO IG = 1,IGL
        IF(KINF.EQ.1.AND.KSUP.EQ.NIBLO.AND.LMESSPASS) THEN
          DO  IJ = 1,KEND
            UU = US(IJ,IG)
            VV = DS(IJ,IG)
            US(IJ,IG) = SQRT(UU**2 + VV**2)
            IF (US(IJ,IG).NE.0.) DS(IJ,IG) = ATAN2(UU,VV)
            IF (DS(IJ,IG).LT.0.) DS(IJ,IG) = DS(IJ,IG) + ZPI
          ENDDO
        ELSE
          DO  IJ = IJS(IG),IJL(IG)
            UU = US(IJ,IG)
            VV = DS(IJ,IG)
            US(IJ,IG) = SQRT(UU**2 + VV**2)
            IF (US(IJ,IG).NE.0.) DS(IJ,IG) = ATAN2(UU,VV)
            IF (DS(IJ,IG).LT.0.) DS(IJ,IG) = DS(IJ,IG) + ZPI
          ENDDO
        ENDIF
      ENDDO 

      IF (ICODE.EQ.1) THEN

!*    3.2  INPUT IS FRICTION VELOCITY.
!          ---------------------------

        DO IG = 1,IGL
          IF(KINF.EQ.1.AND.KSUP.EQ.NIBLO.AND.LMESSPASS) THEN
            DO IJ = 1,KEND
              USTAR = MAX(0.01,US(IJ,IG))
              Z0  = ALPHACH/G*USTAR**2
              CD  = XKAPPA/LOG(10./Z0)
              US(IJ,IG) = USTAR/CD
            ENDDO
          ELSE
            DO IJ = IJS(IG),IJL(IG)
              USTAR = MAX(0.01,US(IJ,IG))
              Z0  = ALPHACH/G*USTAR**2
              CD  = XKAPPA/LOG(10./Z0)
              US(IJ,IG) = USTAR/CD
            ENDDO
          ENDIF
        ENDDO

      ELSEIF (ICODE.EQ.2) THEN

!*    3.3 INPUT WINDS ARE SURFACE STRESSES.
!         ---------------------------------

        DO IG = 1,IGL
          IF(KINF.EQ.1.AND.KSUP.EQ.NIBLO.AND.LMESSPASS) THEN
            DO IJ = 1,KEND
              USTAR = MAX (0.01, SQRT(US(IJ,IG)/ROAIR))
              Z0  = ALPHACH/G*USTAR**2
              CD  = XKAPPA/LOG(10./Z0)
              US(IJ,IG) = USTAR/CD
            ENDDO
          ELSE
            DO IJ = IJS(IG),IJL(IG)
              USTAR = MAX (0.01, SQRT(US(IJ,IG)/ROAIR))
              Z0  = ALPHACH/G*USTAR**2
              CD  = XKAPPA/LOG(10./Z0)
              US(IJ,IG) = USTAR/CD
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!     3.4 SET LOWER WIND SPEED TO WSPMIN m/s
!         ----------------------------------
      DO IG = 1,IGL
        IF(KINF.EQ.1.AND.KSUP.EQ.NIBLO.AND.LMESSPASS) THEN
          DO IJ = 1,KEND
            US(IJ,IG) = MAX(US(IJ,IG),WSPMIN)
          ENDDO
        ELSE
          DO IJ = IJS(IG),IJL(IG)
            US(IJ,IG) = MAX(US(IJ,IG),WSPMIN)
          ENDDO
        ENDIF
      ENDDO
      IF(LFRSTIME) THEN
        WRITE (IU06,*) ' '
        WRITE (IU06,*) ' SUB. WAMWND : '
        WRITE (IU06,*) ' '
        WRITE (IU06,*) ' WIND SPEEDS LOWER THAN ',WSPMIN, ' M/S'
        WRITE (IU06,*) ' WERE RESET TO  ',WSPMIN, ' M/S'
        WRITE (IU06,*) ' '
        WRITE (IU06,*) ' '
        LFRSTIME=.FALSE.
      ENDIF


!*    4. TEST OUTPUT OF WAVE MODEL BLOCKS
!        ---------------------------------

      IF (ITEST.GE.3) THEN
        WRITE (IU06,*) ' '
        WRITE (IU06,*) '      SUB. WAMWND:',
     &   ' WINDFIELDS CONVERTED TO BLOCKS'
        WRITE (IU06,*) ' '
        DO IG = 1,MIN(ITESTB,NBLO)
          IJA = IJS(IG)
          IJE = IJS(IG)+4
          WRITE (IU06,*) ' BLOCK= ',IG
          WRITE (IU06,*) ' '
          WRITE (IU06,*) ' US(IJS - IJS+4): ',(US(I,IG),I=IJA,IJE)
          WRITE (IU06,*) ' DS(IJS - IJS+4): ',(DS(I,IG),I=IJA,IJE)
        ENDDO
      ENDIF

      IF(ALLOCATED(UWND)) DEALLOCATE(UWND)
      IF(ALLOCATED(VWND)) DEALLOCATE(VWND)

      RETURN
      END SUBROUTINE WAMWND
