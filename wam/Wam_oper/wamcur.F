      SUBROUTINE WAMCUR (IUCUR)

! ----------------------------------------------------------------------

!**** *WAMCUR* - ROUTINE PROCESS A CURRENT FIELD.

!     SUSANNE HASSELMANN        JUNE     1990.
!     H. GUNTHER    ECMWF/GKSS  DECEMBER 1990  MODIFIED FOR CYCLE_4.

!*    PURPOSE.
!     --------

!       INPUT AND INTERPOLATION TO BLOCKED CURRENT ARRAYS.

!**   INTERFACE.
!     ----------

!       *CALL* *WAMCUR (IUCUR)*
!          *IUCUR*  - INPUT UNIT FOR CURRENTS.

!     METHOD.
!     -------

!       LINEAR INTERPOLATION.

!     EXTERNALS.
!     ----------

!       *LOCINT*    - INTERPOLATE TO BLOCKS.
!       *READCUR*   - READ A CURRENT FIELD.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NIBLC    ,NBLC
      USE YOWCURR  , ONLY : U        ,V
      USE YOWCURG  , ONLY : DUCLO    ,DUCLA    ,UCLO     ,UCLW     ,
     &            UCUR     ,VCUR     ,KRCI     ,KCCI     ,IPCUR
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJLT
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB

! ----------------------------------------------------------------------

      CHARACTER*12 CDTCUR

! ----------------------------------------------------------------------

!*    1. READ CURRENT FILE.
!        ------------------

      CALL READCUR (IUCUR, CDTCR)

! ----------------------------------------------------------------------

!*    2. INTERPOLATION OF INPUT CURRENT ON MODEL GRID.
!        ---------------------------------------------

      IJA =1

      ALLOCATE(U(0:NIBLC,NBLC))
      ALLOCATE(V(0:NIBLC,NBLC))

      DO IG=1,IGL
        CALL LOCINT (IG, IJA, IJLT(IG), KCCI, KRCI, KCCI, KRCI, IPCUR,
     &   DUCLO, DUCLA, UCLW, UCLO, UCUR, VCUR, U(IJS(IG),IG),
     &   V(IJS(IG),IG))
      ENDDO

      DEALLOCATE(UCUR)
      DEALLOCATE(VCUR)

! ----------------------------------------------------------------------

!*    3. TEST OUTPUT OF WAVE MODEL BLOCKS
!        ---------------------------------

      IF (ITEST.GE.3) THEN
        WRITE (IU06,*) ' '
        WRITE (IU06,*) '      SUB. WAMCUR..',
     &   ' CURRENTS CONVERTED TO BLOCKS'
        WRITE (IU06,*) ' '
        DO IG = 1,ITESTB
          IJA = IJS(IG)
          IJE = IJS(IG)+4
          WRITE (IU06,*) ' BLOCK= ',IG
          WRITE (IU06,*) ' '
          WRITE (IU06,*) ' U(IJS - IJS+4): ',(U(I,IG),I=IJA,IJE)
          WRITE (IU06,*) ' V(IJS - IJS+4): ',(V(I,IG),I=IJA,IJE)
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE WAMCUR
