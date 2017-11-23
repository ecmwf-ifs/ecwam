      SUBROUTINE CHECK (ML, KL, IINPC)

! ----------------------------------------------------------------------

!**** *CHECK* - ROUTINE TO CHECK CONSISTENCY BETWEEN COMPUTED BLOCKS.

!     H.GUNTHER            ECMWF       04/04/1990

!*    PURPOSE.
!     -------

!       *CHECK* CHECKS CONSISTENCY BETWEEN BLOCK INDICES.

!**   INTERFACE.
!     ----------

!       *CALL* *CHECK (ML, KL, IINPC)*
!          *ML*      - NUMBER OF FREQUENCIES.
!          *KL*      - NUMBER OF DIRECTIONS.
!          *IINPC*   - NUMBER INPUT POINTS FROM A PREVIOUS COARSE GRID.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.
!       *OUTPP*     - WRITE OUT A GRID.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      ,
     &            NBLO     ,NIBLO    ,NOVER    ,NIBL1    ,NIBLD    ,
     &            NBLD     ,NIBLC    ,NBLC
      USE YOWPCONS , ONLY : DEG
      USE YOWCPBO  , ONLY : IBOUNC   ,NBOUNC   ,IJARC    ,IGARC
      USE YOWFPBO  , ONLY : IBOUNF   ,NBOUNF   ,IJARF    ,IGARF
      USE YOWCOUT  , ONLY : NGOUT    ,IGAR     ,IJAR
      USE YOWCURR  , ONLY : U        ,V
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL2     ,IJLS     ,
     &            IJL      ,IJLT
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,NX       ,NY       ,
     &            AMOWEP   ,AMOSOP   ,AMOEAP   ,AMONOP   ,XDELLO
      USE YOWSHAL  , ONLY : NDEPTH   ,DEPTH
      USE YOWSTAT  , ONLY : IREFRA
      USE YOWTEST  , ONLY : IU06, ITEST

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "outpp.intfb.h"

      INTEGER, INTENT(IN) :: ML, KL, IINPC

      INTEGER :: IG
      INTEGER :: IU1, IU2, IO, IO1, IO2
      INTEGER :: I, K, IJ, IERR
      INTEGER :: ILEN, IPAGE, LAST, L, IA, IE
      INTEGER :: IJFLAT, IJLLAT, IJMAX, ISEA, IPOI, IOV 

      REAL :: BMOWEP, BMOEAP
      REAL,  ALLOCATABLE :: GRID(:,:)

      CHARACTER(LEN=100) :: TITL
      CHARACTER(LEN=1), ALLOCATABLE :: LST(:,:)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *LST*       CHARACTER  LAND SEA TABLE  L = LAND
!                                             S = SEA
!                                             + = SEA AND OUTPUT POINT.
!      *GRID*      REAL       ARRAY FOR GRIDDED PRINT OUTPUT.

! ----------------------------------------------------------------------

      ALLOCATE(LST(NGX,NGY))
      ALLOCATE(GRID(NGX,NGY))

!*    1. COMPARE LENGTH OF OVERLAPPING LAT.
!        -----------------------------------

      DO IG=1,IGL-1
        IU1 = IJS(IG+1)-1
        IO2 = IJL(IG)-IJLS(IG)+1
        IF (IU1.NE.IO2) THEN
          WRITE (IU06,*) ' *****************************************'
          WRITE (IU06,*) ' *                                       *'
          WRITE (IU06,*) ' *      FATAL ERROR IN SUB. CHECK        *'
          WRITE (IU06,*) ' *      =========================        *'
          WRITE (IU06,*) ' *                                       *'
          WRITE (IU06,*) ' * LENGTH OF FIRST LAT. IN BLOCK IG+1    *'
          WRITE (IU06,*) ' * IS NOT EQUAL TO SECOND TO LAST OF     *'
          WRITE (IU06,*) ' * BLOCK IG                              *'
          WRITE (IU06,*) ' * BLOCK  NUMBER IS IG = ', IG
          WRITE (IU06,*) ' * LENGTH IN BLOCK IG   IS IU1 = ', IU1
          WRITE (IU06,*) ' * LENGTH IN BLOCK IG+1 IS IO2 = ', IO2
          WRITE (IU06,*) ' *                                       *'
          WRITE (IU06,*) ' *****************************************'
        ENDIF
        IU2 = IJL2(IG+1)-IJS(IG+1)+1
        IO1 = IJLT(IG)-IJL(IG)
        IF (IU2.NE.IO1) THEN
          WRITE (IU06,*) ' *****************************************'
          WRITE (IU06,*) ' *                                       *'
          WRITE (IU06,*) ' *      FATAL ERROR IN SUB. CHECK        *'
          WRITE (IU06,*) ' *      =========================        *'
          WRITE (IU06,*) ' *                                       *'
          WRITE (IU06,*) ' * LENGTH OF SECOND LAT. IN BLOCK IG+1   *'
          WRITE (IU06,*) ' * IS NOT EQUAL TO LAST OF BLOCK IG      *'
          WRITE (IU06,*) ' * BLOCK  NUMBER IS IG = ', IG
          WRITE (IU06,*) ' * LENGTH IN BLOCK IG   IS IU2 = ', IU2
          WRITE (IU06,*) ' * LENGTH IN BLOCK IG+1 IS IO1 = ', IO1
          WRITE (IU06,*) ' *                                       *'
          WRITE (IU06,*) ' *****************************************'
        ENDIF
      ENDDO

! ----------------------------------------------------------------------

!*    2. GENERATE LAND SEA TABLE FROM INDEX ARRAYS.
!        ------------------------------------------

      DO K=1,NY
        DO I=1,NX
          LST(I,K) = 'L'
        ENDDO
      ENDDO

      IERR = 0
      DO IG=1,IGL
        DO IJ=IJS(IG),IJL(IG)
          IF (IXLG(IJ,IG).NE.0.OR.KXLT(IJ,IG).NE.0)
     &     LST(IXLG(IJ,IG),KXLT(IJ,IG)) = 'S'
        ENDDO
      ENDDO

!*    2.1 INCLUDE OUTPUT POINTS.
!         ----------------------

      IF (NGOUT.GT.0) THEN
        DO IO=1,NGOUT
          IG = IGAR(IO)
          IJ = IJAR(IO)
          IF (IG.LT.1.OR.IG.GT.IGL) THEN
            IERR = IERR+1
            WRITE (IU06,*) ' ***************************************'
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' *      FATAL ERROR IN SUB. CHECK      *'
            WRITE (IU06,*) ' *      =========================      *'
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' * BLOCK NUMBER OF OUTPUT POINT IS     *'
            WRITE (IU06,*) ' * OUT OF RANGE.                       *'
            WRITE (IU06,*) ' * OUTPUT POINT NUMBER IS  IO = ', IO
            WRITE (IU06,*) ' * BLOCK NUMBER IS         IG = ', IG
            WRITE (IU06,*) ' * MAX. BLOCK NUMBER IS   IGL = ', IGL
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' ***************************************'
            IF (IERR.GT.20) CALL ABORT1
          ENDIF
          IF (IJ.LT.IJS(IG).OR.IJ.GT.IJL(IG)) THEN
            IERR = IERR+1
            WRITE (IU06,*) ' ***************************************'
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' *      FATAL ERROR IN SUB. CHECK      *'
            WRITE (IU06,*) ' *      =========================      *'
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' * GRID POINT NUMBER OF OUTPUT POINT IS*'
            WRITE (IU06,*) ' * OUT OF RANGE.                       *'
            WRITE (IU06,*) ' * OUTPUT POINT NUMBER IS IO = ', IO
            WRITE (IU06,*) ' * GRID POINT NUMBER IS   IJ = ', IJ
            WRITE (IU06,*) ' * MIN. NUMBER IS        IJS = ', IJS(IG)
            WRITE (IU06,*) ' * MAX. NUMBER IS        IJL = ', IJL(IG)
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' ***************************************'
            IF (IERR.GT.20) CALL ABORT1
          ENDIF
          IF (IXLG(IJ,IG).NE.0.OR.KXLT(IJ,IG).NE.0)
     &     LST(IXLG(IJ,IG),KXLT(IJ,IG)) = '+'
        ENDDO
      ENDIF

!*    2.2 INCLUDE COARSE GRID NEST OUTPUT POINTS.
!         ---------------------------------------

      IF (IBOUNC.EQ.1) THEN
        DO IO=1,NBOUNC
          IG = IGARC(IO)
          IJ = IJARC(IO)
          IF (IG.LT.1.OR.IG.GT.IGL) THEN
            IERR = IERR+1
            WRITE (IU06,*) ' ***************************************'
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' *      FATAL ERROR IN SUB. CHECK      *'
            WRITE (IU06,*) ' *      =========================      *'
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' * BLOCK NUMBER OF OUTPUT POINT IS     *'
            WRITE (IU06,*) ' * OUT OF RANGE.                       *'
            WRITE (IU06,*) ' * COARSE BOUNDARY POINT NUMBER IS IO= ',IO
            WRITE (IU06,*) ' * BLOCK NUMBER IS         IG = ', IG
            WRITE (IU06,*) ' * MAX. BLOCK NUMBER IS   IGL = ', IGL
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' ***************************************'
            IF (IERR.GT.20) CALL ABORT1
          ENDIF
          IF (IJ.LT.IJS(IG).OR.IJ.GT.IJL(IG)) THEN
            IERR = IERR+1
            WRITE (IU06,*) ' ***************************************'
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' *      FATAL ERROR IN SUB. CHECK      *'
            WRITE (IU06,*) ' *      =========================      *'
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' * GRID POINT NUMBER OF OUTPUT POINT IS*'
            WRITE (IU06,*) ' * OUT OF RANGE.                       *'
            WRITE (IU06,*) ' * COARSE BOUNDARY POINT NUMBER IS IO= ',IO
            WRITE (IU06,*) ' * GRID POINT NUMBER IS   IJ = ', IJ
            WRITE (IU06,*) ' * MIN. NUMBER IS        IJS = ', IJS(IG)
            WRITE (IU06,*) ' * MAX. NUMBER IS        IJL = ', IJL(IG)
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' ***************************************'
            IF (IERR.GT.20) CALL ABORT1
          ENDIF
          IF (IXLG(IJ,IG).NE.0.OR.KXLT(IJ,IG).NE.0)
     &     LST(IXLG(IJ,IG),KXLT(IJ,IG)) = '/'
        ENDDO
      ENDIF

!*    2.3 INCLUDE FINE GRID NEST INPUT POINTS.
!         ------------------------------------

      IF (IBOUNF.EQ.1) THEN
        DO IO=1,NBOUNF
          IG = IGARF(IO)
          IJ = IJARF(IO)
          IF (IG.LT.1.OR.IG.GT.IGL) THEN
            IERR = IERR+1
            WRITE (IU06,*) ' ***************************************'
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' *      FATAL ERROR IN SUB. CHECK      *'
            WRITE (IU06,*) ' *      =========================      *'
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' * BLOCK NUMBER OF OUTPUT POINT IS     *'
            WRITE (IU06,*) ' * OUT OF RANGE.                       *'
            WRITE (IU06,*) ' * FINE BOUNDARY POINT NUMBER IS  IO =',IO
            WRITE (IU06,*) ' * BLOCK NUMBER IS         IG = ', IG
            WRITE (IU06,*) ' * MAX. BLOCK NUMBER IS   IGL = ', IGL
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' ***************************************'
            IF (IERR.GT.20) CALL ABORT1
          ENDIF
          IF (IJ.LT.IJS(IG).OR.IJ.GT.IJL(IG)) THEN
            IERR = IERR+1
            WRITE (IU06,*) ' ***************************************'
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' *      FATAL ERROR IN SUB. CHECK      *'
            WRITE (IU06,*) ' *      =========================      *'
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' * GRID POINT NUMBER OF OUTPUT POINT IS*'
            WRITE (IU06,*) ' * OUT OF RANGE.                       *'
            WRITE (IU06,*) ' * FINE BOUNDARY POINT NUMBER IS IO= ', IO
            WRITE (IU06,*) ' * GRID POINT NUMBER IS   IJ = ', IJ
            WRITE (IU06,*) ' * MIN. NUMBER IS        IJS = ', IJS(IG)
            WRITE (IU06,*) ' * MAX. NUMBER IS        IJL = ', IJL(IG)
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' ***************************************'
            IF (IERR.GT.20) CALL ABORT1
          ENDIF
          IF (IXLG(IJ,IG).NE.0.OR.KXLT(IJ,IG).NE.0)
     &     LST(IXLG(IJ,IG),KXLT(IJ,IG)) = 'B'
        ENDDO
      ENDIF

!*    2.4 PRINT LAND SEA MAP.
!         -------------------


      IF (ITEST.GE.3) THEN
        ILEN = 120
        IPAGE = (NX+ILEN-1)/ILEN
        IF (IPAGE.GT.1) THEN
          LAST = (NX-ILEN*(IPAGE-1)+IPAGE-2)/(IPAGE-1)
          IF (LAST.LE.10) THEN
            ILEN = ILEN + 10
            IPAGE = (NX+ILEN-1)/ILEN
          ENDIF
        ENDIF
        DO L=1,IPAGE
          IA = (L-1)*ILEN
          IE = MIN(IA+ILEN,NX)
          IA = IA + 1
          BMOWEP = AMOWEP +REAL(IA-1)*XDELLO
          BMOEAP = AMOWEP +REAL(IE-1)*XDELLO
          WRITE (IU06,'(1H1,'' LAND SEA MAP OF FULL GRID '',
     &     ''   L = LAND  S = SEA  + = OUTPUT POINT'',
     &     ''                PAGE: '',I2)') L
          WRITE (IU06,'(2X,''LONGITUDE IS FROM '',F7.2,'' TO '',F7.2)')
     &     BMOWEP, BMOEAP
          WRITE (IU06,'(2X,''LATITUDE  IS FROM '',F7.2,'' TO '',F7.2)')
     &     AMONOP, AMOSOP
          WRITE (IU06,'(2X,130I1)') (MOD(I,10),I=IA,IE)
          DO K=NY,1,-1
            WRITE (IU06,'(1X,I1,130A1)') MOD(K,10),(LST(I,K),I=IA,IE)
          ENDDO
          WRITE (IU06,'(2X,130I1)') (MOD(I,10),I=IA,IE)
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    3. OUTPUT OF DEPTH FIELD.
!        ----------------------

      IF (ITEST.GE.3) THEN
        DO K=1,NY
          DO I=1,NX
            GRID(I,K) = 99999
          ENDDO
        ENDDO

        DO IG=1,IGL
          DO IJ=IJS(IG),IJL(IG)
            GRID(IXLG(IJ,IG),NY+1-KXLT(IJ,IG)) = MIN(DEPTH(IJ,IG),999.)
          ENDDO
        ENDDO

        TITL = 'WATER DEPTH IN METERS. '//
     &         '(DEPTH DEEPER THAN 999M ARE PRINTED AS 999)'
        CALL OUTPP('00000000000000',IU06,NGX,NGY,NX,NY,TITL,1.,GRID,
     &              AMOWEP,AMOSOP,AMOEAP,AMONOP)
      ENDIF

! ----------------------------------------------------------------------

!*    4. OUTPUT OF CURRENT FIELD.
!        ------------------------


      IF (ITEST.GE.3) THEN

        IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN
 
        DO K=1,NY
          DO I=1,NX
            GRID(I,K) = 99999
          ENDDO
        ENDDO

        DO IG=1,IGL
          DO IJ=IJS(IG),IJL(IG)
            GRID(IXLG(IJ,IG),NY+1-KXLT(IJ,IG)) =
     &       SQRT(U(IJ,IG)**2+V(IJ,IG)**2)
          ENDDO
        ENDDO

        TITL = 'CURRENT SPEED IN 0.01 METRES/SECOND '
        CALL OUTPP('00000000000000',IU06,NGX,NGY,NX,NY,TITL,100.,GRID,
     &              AMOWEP,AMOSOP,AMOEAP,AMONOP)

        DO K=1,NY
          DO I=1,NX
            GRID(I,K) = 99999
          ENDDO
        ENDDO

        DO IG=1,IGL
          DO IJ=IJS(IG),IJL(IG)
            GRID(IXLG(IJ,IG),NY+1-KXLT(IJ,IG)) =
     &       MOD(ATAN2(U(IJ,IG),V(IJ,IG)+0.1E-10)*DEG+360.,360.)
          ENDDO
        ENDDO

        TITL = 'CURRENT DIRECTIOM IN DEGREES (CLOCKWISE FROM NORTH)'
        CALL OUTPP('00000000000000',IU06,NGX,NGY,NX,NY,TITL,1.,GRID,
     &              AMOWEP,AMOSOP,AMOEAP,AMONOP)

        ENDIF

      ENDIF

! ----------------------------------------------------------------------

!*    5. OUTPUT OF OVERALL GRID INFORMATION.
!        -----------------------------------

 5000 CONTINUE
      WRITE (IU06,'(1H1,'' GRID SUMMERY:'')')
      WRITE (IU06,*) ' NUMBER OF BLOCKS GENERATED IS IGL ....: ', IGL
      IJFLAT = 0
      IJLLAT = 0
      IJMAX = 0
      ISEA = 0
      IPOI = 0
      DO IG=1,IGL
        IJFLAT= MAX(IJFLAT,IJS(IG)-1)
        IJLLAT= MAX(IJLLAT,IJLT(IG)-IJL(IG))
        IJMAX = MAX(IJMAX,IJLT(IG))
        IPOI  = IPOI + IJLT(IG)
        ISEA  = ISEA + IJL(IG)-IJS(IG) + 1
      ENDDO
      IOV = IPOI-ISEA
      WRITE (IU06,*) ' MAXIMUM NUMBER OF POINTS IN A BLOCK ..: ', IJMAX
      WRITE (IU06,*) ' TOTAL NUMBER OF POINT IN ALL BLOCKS ..: ', IPOI
      WRITE (IU06,*) ' TOTAL NUMBER OF SEA POINTS ...........: ', ISEA
      WRITE (IU06,*) ' TOTAL NUMBER OF POINTS IN OVERLAP.....: ', IOV
      WRITE (IU06,*) ' MAXIMUM LENGTH OF FIRST LAT OF A BLOCK: ', IJFLAT
      WRITE (IU06,*) ' MAXIMUM LENGTH OF LAST  LAT OF A BLOCK: ', IJLLAT

! ----------------------------------------------------------------------

!*    6. OUTPUT OF OPTIMAL DIMENSIONS.
!        -----------------------------

      WRITE (IU06,'(//,'' DIMENSIONS OF ARRAYS, WHICH ARE USED'',
     &             '' IN PRESET AND CHIEF '',/)')
      WRITE (IU06,'(''                                     DEFINED'',
     &           ''      USED'',''  REQUIRED'')')
      WRITE (IU06,'('' NUMBER OF DIRECTIONS        NANG '', 3I10)')
     &           NANG, KL, KL
      WRITE (IU06,'('' NUMBER OF FREQUENCIES       NFRE '', 3I10)')
     &           NFRE, ML, ML
      WRITE (IU06,'('' NUMBER LONGITUDE GRID POINTS NGX '', 3I10)')
     &           NGX, NX, NX
      WRITE (IU06,'('' NUMBER LATITUDE GRID POINTS  NGY '', 3I10)')
     &           NGY, NY, NY
      WRITE (IU06,'('' NUMBER OF BLOCKS            NBLO '', 3I10)')
     &           NBLO, IGL, NBLO
      WRITE (IU06,'('' MAXIMUM BLOCK LENGTH       NIBLO '', 3I10)')
     &           NIBLO, IJMAX, NIBLO
      WRITE (IU06,'('' NUMBER POINTS IN FIRST LAT NOVER '', 3I10)')
     &           NOVER, MAX(1,IJFLAT), MAX(1,IJFLAT)
      WRITE (IU06,'('' NUMBER OF OUTPUT POINTS    NGOUT '', 3I10)')
     &           NGOUT, MAX(1,NGOUT), NGOUT 

      IF (IGL.EQ.1) THEN
        WRITE (IU06,'('' 1 BLOCK VERSION            NIBL1 '', 3I10)')
     &   NIBL1, 1, 1
      ELSE
        WRITE (IU06,'('' MULTI BLOCK VERSION        NIBL1 '', 3I10)')
     &   NIBL1, IJMAX, NIBLO
      ENDIF

      WRITE (IU06,'('' SHALLOW WATER TABLE LEN.  NDEPTH '', 3I10)')
     &           NDEPTH, NDEPTH, NDEPTH

      WRITE (IU06,'(/,'' THE DIMENSIONS IN PRESET AND CHIEF HAVE TO '',
     &             '' BE THE VALUES IN COLUMN - REQUIRED - '')')
      WRITE (IU06,'(  '' IF YOU WANT TO USE THE OPTIMAL DIMENSION'',
     &             '' LENGTH IN THE WAMODEL, THEN  '',/,
     &             '' RERUN PREPROC WITH THE DIMENSION'',
     &             '' GIVEN AS -USED-'')')

      DEALLOCATE(LST)
      DEALLOCATE(GRID)

      RETURN
      END SUBROUTINE CHECK
