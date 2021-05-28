      SUBROUTINE CHECK (IINPC)

! ----------------------------------------------------------------------

!**** *CHECK* - ROUTINE TO CHECK CONSISTENCY BETWEEN COMPUTED BLOCKS.

!     H.GUNTHER            ECMWF       04/04/1990

!*    PURPOSE.
!     -------

!       *CHECK* CHECKS CONSISTENCY BETWEEN BLOCK INDICES.

!**   INTERFACE.
!     ----------

!       *CALL* *CHECK (IINPC)*
!          *IINPC*   - NUMBER INPUT POINTS FROM A PREVIOUS COARSE GRID.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_RED ,              &
     &            NGX      ,NGY      ,NIBLO
      USE YOWPCONS , ONLY : DEG
      USE YOWCPBO  , ONLY : IBOUNC   ,NBOUNC   ,IJARC
      USE YOWFPBO  , ONLY : IBOUNF   ,NBOUNF   ,IJARF
      USE YOWCOUT  , ONLY : NGOUT    ,IJAR
      USE YOWGRID  , ONLY : IJS      ,IJL
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,NX       ,NY       ,    &
     &            AMOWEP   ,AMOSOP   ,AMOEAP   ,AMONOP   ,XDELLO
      USE YOWSHAL  , ONLY : NDEPTH   ,DEPTH
      USE YOWTEST  , ONLY : IU06, ITEST

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IINPC

      INTEGER(KIND=JWIM) :: IO
      INTEGER(KIND=JWIM) :: I, K, IJ, IERR
      INTEGER(KIND=JWIM) :: ILEN, IPAGE, LAST, L, IA, IE

      REAL(KIND=JWRB) :: BMOWEP, BMOEAP
      REAL(KIND=JWRB) :: GRID(NGX,NGY)

      CHARACTER(LEN=100) :: TITL
      CHARACTER(LEN=1) :: LST(NGX,NGY)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *LST*       CHARACTER  LAND SEA TABLE  L = LAND
!                                             S = SEA
!                                             + = SEA AND OUTPUT POINT.
!      *GRID*      REAL       ARRAY FOR GRIDDED PRINT OUTPUT.

! ----------------------------------------------------------------------


!*    2. GENERATE LAND SEA TABLE FROM INDEX ARRAYS.
!        ------------------------------------------

      DO K=1,NY
        DO I=1,NX
          LST(I,K) = 'L'
        ENDDO
      ENDDO

      IERR = 0
      DO IJ=IJS,IJL
        IF (IXLG(IJ).NE.0.OR.KXLT(IJ).NE.0)                     &
     &   LST(IXLG(IJ),KXLT(IJ)) = 'S'
      ENDDO

!*    2.1 INCLUDE OUTPUT POINTS.
!         ----------------------

      IF (NGOUT.GT.0) THEN
        DO IO=1,NGOUT
          IJ = IJAR(IO)
          IF (IJ.LT.IJS .OR. IJ.GT.IJL) THEN
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
            WRITE (IU06,*) ' * MIN. NUMBER IS        IJS = ', IJS
            WRITE (IU06,*) ' * MAX. NUMBER IS        IJL = ', IJL
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' ***************************************'
            IF (IERR.GT.20) CALL ABORT1
          ENDIF
          IF (IXLG(IJ).NE.0.OR.KXLT(IJ).NE.0) 
     &     LST(IXLG(IJ),KXLT(IJ)) = '+'
        ENDDO
      ENDIF

!*    2.2 INCLUDE COARSE GRID NEST OUTPUT POINTS.
!         ---------------------------------------

      IF (IBOUNC.EQ.1) THEN
        DO IO=1,NBOUNC
          IJ = IJARC(IO)
          IF (IJ.LT.IJS.OR.IJ.GT.IJL) THEN
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
            WRITE (IU06,*) ' * MIN. NUMBER IS        IJS = ', IJS
            WRITE (IU06,*) ' * MAX. NUMBER IS        IJL = ', IJL
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' ***************************************'
            IF (IERR.GT.20) CALL ABORT1
          ENDIF
          IF (IXLG(IJ).NE.0.OR.KXLT(IJ).NE.0)                     &
     &     LST(IXLG(IJ),KXLT(IJ)) = '/'
        ENDDO
      ENDIF

!*    2.3 INCLUDE FINE GRID NEST INPUT POINTS.
!         ------------------------------------

      IF (IBOUNF.EQ.1) THEN
        DO IO=1,NBOUNF
          IJ = IJARF(IO)
          IF (IJ.LT.IJS.OR.IJ.GT.IJL) THEN
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
            WRITE (IU06,*) ' * MIN. NUMBER IS        IJS = ', IJS
            WRITE (IU06,*) ' * MAX. NUMBER IS        IJL = ', IJL
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' ***************************************'
            IF (IERR.GT.20) CALL ABORT1
          ENDIF
          IF (IXLG(IJ).NE.0.OR.KXLT(IJ).NE.0)                     &
     &     LST(IXLG(IJ),KXLT(IJ)) = 'B'
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
          WRITE (IU06,'(1H1,'' LAND SEA MAP OF FULL GRID '',            &
     &     ''   L = LAND  S = SEA  + = OUTPUT POINT'',                  &
     &     ''                PAGE: '',I2)') L
          WRITE (IU06,'(2X,''LONGITUDE IS FROM '',F7.2,'' TO '',F7.2)') &
     &     BMOWEP, BMOEAP
          WRITE (IU06,'(2X,''LATITUDE  IS FROM '',F7.2,'' TO '',F7.2)') &
     &     AMONOP, AMOSOP
          WRITE (IU06,'(2X,130I1)') (MOD(I,10),I=IA,IE)
          DO K=NY,1,-1
            WRITE (IU06,'(1X,I1,130A1)') MOD(K,10),(LST(I,K),I=IA,IE)
          ENDDO
          WRITE (IU06,'(2X,130I1)') (MOD(I,10),I=IA,IE)
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    5. OUTPUT OF OVERALL GRID INFORMATION.
!        -----------------------------------

 5000 CONTINUE
      WRITE (IU06,'(1H1,'' GRID SUMMERY:'')')

! ----------------------------------------------------------------------

!*    6. OUTPUT OF OPTIMAL DIMENSIONS.
!        -----------------------------

      WRITE (IU06,'(//,'' DIMENSIONS OF ARRAYS, WHICH ARE USED'',       &
     &             '' IN PRESET AND CHIEF '',/)')
      WRITE (IU06,'(''                                     DEFINED'',   &
     &           ''      USED'',''  REQUIRED'')')
      WRITE (IU06,'('' NUMBER OF DIRECTIONS        NANG '', 3I10)')     &
     &           NANG, NANG, NANG 
      WRITE (IU06,'('' NUMBER OF FREQUENCIES       NFRE '', 3I10)')     &
     &           NFRE, NFRE, NFRE 
      WRITE (IU06,'('' NUMBER OF FREQUENCIES   NFRE_RED '', 3I10)')     &
     &           NFRE_RED, NFRE_RED, NFRE_RED 
      WRITE (IU06,'('' NUMBER LONGITUDE GRID POINTS NGX '', 3I10)')     &
     &           NGX, NX, NX
      WRITE (IU06,'('' NUMBER LATITUDE GRID POINTS  NGY '', 3I10)')     &
     &           NGY, NY, NY
      WRITE (IU06,'('' MAXIMUM BLOCK LENGTH       NIBLO '', 3I10)')     &
     &           NIBLO
      WRITE (IU06,'('' NUMBER OF OUTPUT POINTS    NGOUT '', 3I10)')     &
     &           NGOUT, MAX(1,NGOUT), NGOUT 

      WRITE (IU06,'('' SHALLOW WATER TABLE LEN.  NDEPTH '', 3I10)')     &
     &           NDEPTH, NDEPTH, NDEPTH

      WRITE (IU06,'(/,'' THE DIMENSIONS IN PRESET AND CHIEF HAVE TO '', &
     &             '' BE THE VALUES IN COLUMN - REQUIRED - '')')
      WRITE (IU06,'(  '' IF YOU WANT TO USE THE OPTIMAL DIMENSION'',    &
     &             '' LENGTH IN THE WAMODEL, THEN  '',/,                &
     &             '' RERUN PREPROC WITH THE DIMENSION'',               &
     &             '' GIVEN AS -USED-'')')


      END SUBROUTINE CHECK
