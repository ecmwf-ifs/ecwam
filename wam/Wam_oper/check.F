      SUBROUTINE CHECK (IREFRA, ML, KL, IINPC)

C ----------------------------------------------------------------------
C
C**** *CHECK* - ROUTINE TO CHECK CONSISTENCY BETWEEN COMPUTED BLOCKS.
C
C     H.GUNTHER            ECMWF       04/04/1990
C
C*    PURPOSE.
C     -------
C
C       *CHECK* CHECKS CONSISTENCY BETWEEN BLOCK INDICES.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *CHECK (IREFRA, ML, KL, IINPC)*
C          *IREFRA*  - REFRACTION OPTION.
C          *ML*      - NUMBER OF FREQUENCIES.
C          *KL*      - NUMBER OF DIRECTIONS.
C          *IINPC*   - NUMBER INPUT POINTS FROM A PREVIOUS COARSE GRID.
C
C     METHOD.
C     -------
C
C       NONE.
C
C     EXTERNALS.
C     ----------
C
C       *ABORT1*     - TERMINATES PROCESSING.
C       *OUTPP*     - WRITE OUT A GRID.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "parameter_preproc.h"
C
#include "comcpbo.h"
C
#include "comfpbo.h"
C
#include "comcout.h"
C
#include "comcurr.h"
C
#include "comgrid.h"
C
#include "commap.h"
C
#include "comshal.h"
C
#include "comtest.h"
C
C ----------------------------------------------------------------------
C
#include "parcons.h"
C
      CHARACTER*1 LST(NGX,NGY)
      CHARACTER TITL*100
      DIMENSION GRID(NGX,NGY)
C
C*     VARIABLE.   TYPE.     PURPOSE.
C      ---------   -------   --------
C      *LSTAB*     CHARACTER  LAND SEA TABLE  L = LAND
C                                             S = SEA
C                                             + = SEA AND OUTPUT POINT.
C      *GRID*      REAL       ARRAY FOR GRIDDED PRINT OUTPUT.
C
C ----------------------------------------------------------------------
C
C*    1. COMPARE LENGTH OF OVERLAPPING LAT.
C        -----------------------------------
C
 1000 CONTINUE
      DO 1001 IG=1,IGL-1
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
 1001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    2. GENERATE LAND SEA TABLE FROM INDEX ARRAYS.
C        ------------------------------------------
C
 2000 CONTINUE
      DO K=1,NY
        DO I=1,NX
          LST(I,K) = 'L'
        ENDDO
      ENDDO

      IERR = 0
      DO 2002 IG=1,IGL
         DO 2003 IJ=IJS(IG),IJL(IG)
            IF (IXLG(IJ,IG).NE.0.OR.KXLT(IJ,IG).NE.0)
     1          LST(IXLG(IJ,IG),KXLT(IJ,IG)) = 'S'
 2003    CONTINUE
 2002 CONTINUE
C
C*    2.1 INCLUDE OUTPUT POINTS.
C         ----------------------
C
 2100 CONTINUE
      IF (NGOUT.GT.0) THEN
         DO 2101 IO=1,NGOUT
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
     1          LST(IXLG(IJ,IG),KXLT(IJ,IG)) = '+'
 2101    CONTINUE
      ENDIF
C
C*    2.2 INCLUDE COARSE GRID NEST OUTPUT POINTS.
C         ---------------------------------------
C
 2200 CONTINUE
      IF (IBOUNC.EQ.1) THEN
         DO 2201 IO=1,NBOUNC
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
               WRITE (IU06,*) ' * BOUNDARY POINT NUMBER IS  IO = ', IO
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
               WRITE (IU06,*) ' * BOUNDARY POINT NUMBER IS IO = ', IO
               WRITE (IU06,*) ' * GRID POINT NUMBER IS   IJ = ', IJ
               WRITE (IU06,*) ' * MIN. NUMBER IS        IJS = ', IJS(IG)
               WRITE (IU06,*) ' * MAX. NUMBER IS        IJL = ', IJL(IG)
               WRITE (IU06,*) ' *                                     *'
               WRITE (IU06,*) ' ***************************************'
               IF (IERR.GT.20) CALL ABORT1
            ENDIF
            IF (IXLG(IJ,IG).NE.0.OR.KXLT(IJ,IG).NE.0)
     1          LST(IXLG(IJ,IG),KXLT(IJ,IG)) = '/'
 2201    CONTINUE
      ENDIF
C
C*    2.3 INCLUDE FINE GRID NEST INPUT POINTS.
C         ------------------------------------
C
 2300 CONTINUE
      IF (IBOUNF.EQ.1) THEN
         DO 2301 IO=1,NBOUNF
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
               WRITE (IU06,*) ' * BOUNDARY POINT NUMBER IS  IO = ', IO
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
               WRITE (IU06,*) ' * BOUNDARY POINT NUMBER IS IO = ', IO
               WRITE (IU06,*) ' * GRID POINT NUMBER IS   IJ = ', IJ
               WRITE (IU06,*) ' * MIN. NUMBER IS        IJS = ', IJS(IG)
               WRITE (IU06,*) ' * MAX. NUMBER IS        IJL = ', IJL(IG)
               WRITE (IU06,*) ' *                                     *'
               WRITE (IU06,*) ' ***************************************'
               IF (IERR.GT.20) CALL ABORT1
            ENDIF
            IF (IXLG(IJ,IG).NE.0.OR.KXLT(IJ,IG).NE.0)
     1          LST(IXLG(IJ,IG),KXLT(IJ,IG)) = 'B'
 2301    CONTINUE
      ENDIF
C
C*    2.4 PRINT LAND SEA MAP.
C         -------------------
C
 2400 CONTINUE
      ILEN = 120
      IPAGE = (NX+ILEN-1)/ILEN
      IF (IPAGE.GT.1) THEN
         LAST = (NX-ILEN*(IPAGE-1)+IPAGE-2)/(IPAGE-1)
         IF (LAST.LE.10) THEN
            ILEN = ILEN + 10
            IPAGE = (NX+ILEN-1)/ILEN
         ENDIF
      ENDIF
      DO 2401 L=1,IPAGE
         IA = (L-1)*ILEN
         IE = MIN(IA+ILEN,NX)
         IA = IA + 1
         BMOWEP = AMOWEP +REAL(IA-1)*XDELLO
         BMOEAP = AMOWEP +REAL(IE-1)*XDELLO
         WRITE (IU06,'(1H1,'' LAND SEA MAP OF FULL GRID '',
     1               ''   L = LAND  S = SEA  + = OUTPUT POINT'',
     2               ''                PAGE: '',I2)') L
         WRITE (IU06,'(2X,''LONGITUDE IS FROM '',F7.2,'' TO '',F7.2)')
     1              BMOWEP, BMOEAP
         WRITE (IU06,'(2X,''LATITUDE  IS FROM '',F7.2,'' TO '',F7.2)')
     1              AMONOP, AMOSOP
         WRITE (IU06,'(2X,130I1)') (MOD(I,10),I=IA,IE)
         DO 2402 K=NY,1,-1
            WRITE (IU06,'(1X,I1,130A1)') MOD(K,10),(LST(I,K),I=IA,IE)
 2402    CONTINUE
         WRITE (IU06,'(2X,130I1)') (MOD(I,10),I=IA,IE)
 2401 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    3. OUTPUT OF DEPTH FIELD.
C        ----------------------
C
 3000 CONTINUE
      DO K=1,NY
        DO I=1,NX
          GRID(I,K) = 99999
        ENDDO
      ENDDO

      DO 3002 IG=1,IGL
         DO 3003 IJ=IJS(IG),IJL(IG)
             GRID(IXLG(IJ,IG),NY+1-KXLT(IJ,IG)) = MIN(DEPTH(IJ,IG),999.)
 3003    CONTINUE
 3002 CONTINUE

      TITL = 'WATER DEPTH IN METERS. '//
     1       '(DEPTH DEEPER THAN 999M ARE PRINTED AS 999)'
      CALL OUTPP ('000000000000',IU06, NGX, NGY, NX, NY, TITL, 1., GRID,
     1            AMOWEP, AMOSOP, AMOEAP, AMONOP)
C
C ----------------------------------------------------------------------
C
C*    4. OUTPUT OF CURRENT FIELD.
C        ------------------------
C
 4000 CONTINUE

      IF (IREFRA.LT.2) GOTO 5000

      DO K=1,NY
        DO I=1,NX
          GRID(I,K) = 99999
        ENDDO
      ENDDO

      DO 4002 IG=1,IGL
         DO 4003 IJ=IJS(IG),IJL(IG)
             GRID(IXLG(IJ,IG),NY+1-KXLT(IJ,IG)) =
     1          SQRT(U(IJ,IG)**2+V(IJ,IG)**2)
 4003    CONTINUE
 4002 CONTINUE

      TITL = 'CURRENT SPEED IN 0.01 METRES/SECOND '
      CALL OUTPP ('000000000000',IU06, NGX, NGY,NX,NY,TITL, 100., GRID,
     1            AMOWEP, AMOSOP, AMOEAP, AMONOP)

      DO K=1,NY
        DO I=1,NX
          GRID(I,K) = 99999
        ENDDO
      ENDDO

      DO 4102 IG=1,IGL
         DO 4103 IJ=IJS(IG),IJL(IG)
             GRID(IXLG(IJ,IG),NY+1-KXLT(IJ,IG)) =
     1           MOD(ATAN2(U(IJ,IG),V(IJ,IG)+0.1E-10)*DEG+360.,360.)
 4103    CONTINUE
 4102 CONTINUE

      TITL = 'CURRENT DIRECTIOM IN DEGREES (CLOCKWISE FROM NORTH)'
      CALL OUTPP ('000000000000',IU06, NGX, NGY, NX, NY, TITL, 1., GRID,
     1            AMOWEP, AMOSOP, AMOEAP, AMONOP)
C
C ----------------------------------------------------------------------
C
C*    5. OUTPUT OF OVERALL GRID INFORMATION.
C        -----------------------------------
C
 5000 CONTINUE
      WRITE (IU06,'(1H1,'' GRID SUMMERY:'')')
      WRITE (IU06,*) ' NUMBER OF BLOCKS GENERATED IS IGL ....: ', IGL
      IJFLAT = 0
      IJLLAT = 0
      IJMAX = 0
      ISEA = 0
      IPOI = 0
      DO 5001 IG=1,IGL
         IJFLAT= MAX(IJFLAT,IJS(IG)-1)
         IJLLAT= MAX(IJLLAT,IJLT(IG)-IJL(IG))
         IJMAX = MAX(IJMAX,IJLT(IG))
         IPOI  = IPOI + IJLT(IG)
         ISEA  = ISEA + IJL(IG)-IJS(IG) + 1
 5001 CONTINUE
      IOV = IPOI-ISEA
      WRITE (IU06,*) ' MAXIMUM NUMBER OF POINTS IN A BLOCK ..: ', IJMAX
      WRITE (IU06,*) ' TOTAL NUMBER OF POINT IN ALL BLOCKS ..: ', IPOI
      WRITE (IU06,*) ' TOTAL NUMBER OF SEA POINTS ...........: ', ISEA
      WRITE (IU06,*) ' TOTAL NUMBER OF POINTS IN OVERLAP.....: ', IOV
      WRITE (IU06,*) ' MAXIMUM LENGTH OF FIRST LAT OF A BLOCK: ', IJFLAT
      WRITE (IU06,*) ' MAXIMUM LENGTH OF LAST  LAT OF A BLOCK: ', IJLLAT
C
C ----------------------------------------------------------------------
C
C*    6. OUTPUT OF OPTIMAL DIMENSIONS.
C        -----------------------------
C
 6000 CONTINUE
      WRITE (IU06,'(//,'' DIMENSIONS OF ARRAYS, WHICH ARE USED'',
     1             '' IN PRESET AND CHIEF '',/)')
      WRITE (IU06,'(''                                     DEFINED'',
     1           ''      USED'',''  REQUIRED'')')
      WRITE (IU06,'('' NUMBER OF DIRECTIONS        NANG '', 3I10)')
     1           NANG, KL, KL
      WRITE (IU06,'('' NUMBER OF FREQUENCIES       NFRE '', 3I10)')
     1           NFRE, ML, ML
      WRITE (IU06,'('' NUMBER LONGITUDE GRID POINTS NGX '', 3I10)')
     1           NGX, NX, NX
      WRITE (IU06,'('' NUMBER LATITUDE GRID POINTS  NGY '', 3I10)')
     1           NGY, NY, NY
      WRITE (IU06,'('' NUMBER OF BLOCKS            NBLO '', 3I10)')
     1           NBLO, IGL, NBLO
      WRITE (IU06,'('' MAXIMUM BLOCK LENGTH       NIBLO '', 3I10)')
     1           NIBLO, IJMAX, NIBLO
      WRITE (IU06,'('' NUMBER POINTS IN FIRST LAT NOVER '', 3I10)')
     1           NOVER, MAX(1,IJFLAT), MAX(1,IJFLAT)
      WRITE (IU06,'('' NUMBER OF OUTPUT POINTS    MOUTP '', 3I10)')
     1           MOUTP, MAX(1,NGOUT), MOUTP
      WRITE (IU06,'('' BOUNDARY POINTS (COURSE).  NMAXC '', 3I10)')
     1           NMAXC, MAX(1,NBOUNC), MAX(1,NBOUNC)
      WRITE (IU06,'('' BOUNDARY POINTS (FINE).    NMAXF '', 3I10)')
     1           NMAXF, MAX(1,NBOUNF), MAX(1,NBOUNF)
      WRITE (IU06,'('' BOUNDARY POINTS (INPUT)    NBINP '', 3I10)')
     1           NBINP, MAX(1,IINPC), MAX(1,IINPC)

      IF (IGL.EQ.1) THEN
         WRITE (IU06,'('' 1 BLOCK VERSION           NIBL1 '', 3I10)')
     1           NIBL1, 1, 1
      ELSE
         WRITE (IU06,'('' MULTI BLOCK VERSION        NIBL1 '', 3I10)')
     1           NIBL1, IJMAX, NIBLO
      ENDIF
      IF (IREFRA.EQ.0) THEN
         WRITE (IU06,'('' DEPTH REFRAC. BLOCK LENGTH NIBLD '', 3I10)')
     1           NIBLD, 1, 1
         WRITE (IU06,'('' DEPTH REFRAC. BLOCK NUMBER  NBLD '', 3I10)')
     1           NBLD, 1, 1
         WRITE (IU06,'('' CURR. REFRAC. BLOCK LENGTH NIBLC '', 3I10)')
     1           NIBLC, 1, 1
         WRITE (IU06,'('' CURR. REFRAC. BLOCK NUMBER  NBLC '', 3I10)')
     1           NBLC, 1, 1
      ELSE IF (IREFRA.EQ.1) THEN
         WRITE (IU06,'('' DEPTH REFRAC. BLOCK LENGTH NIBLD '', 3I10)')
     1           NIBLD, IJMAX, NIBLO
         WRITE (IU06,'('' DEPTH REFRAC. BLOCK NUMBER  NBLD '', 3I10)')
     1           NBLD, IGL, NBLO
         WRITE (IU06,'('' CURR. REFRAC. BLOCK LENGTH NIBLC '', 3I10)')
     1           NIBLC, 1, 1
         WRITE (IU06,'('' CURR. REFRAC. BLOCK NUMBER  NBLC '', 3I10)')
     1           NBLC, 1, 1
      ELSE IF (IREFRA.EQ.2) THEN
         WRITE (IU06,'('' DEPTH REFRAC. BLOCK LENGTH NIBLD '', 3I10)')
     1           NIBLD, IJMAX, NIBLO
         WRITE (IU06,'('' DEPTH REFRAC. BLOCK NUMBER  NBLD '', 3I10)')
     1           NBLD, IGL, NBLO
         WRITE (IU06,'('' CURR. REFRAC. BLOCK LENGTH NIBLC '', 3I10)')
     1           NIBLC, IJMAX, NIBLO
         WRITE (IU06,'('' CURR. REFRAC. BLOCK NUMBER  NBLC '', 3I10)')
     1           NBLC, IGL, NBLO
      ENDIF

      WRITE (IU06,'('' SHALLOW WATER TABLE LEN.  NDEPTH '', 3I10)')
     1           NDEPTH, NDEPTH, NDEPTH

      WRITE (IU06,'(/,'' THE DIMENSIONS IN PRESET AND CHIEF HAVE TO '',
     1             '' BE THE VALUES IN COLUMN - REQUIRED - '')')
      WRITE (IU06,'(  '' IF YOU WANT TO USE THE OPTIMAL DIMENSION'',
     1             '' LENGTH IN THE WAMODEL, THEN  '',/,
     2             '' RERUN PREPROC WITH THE DIMENSION'',
     3             '' GIVEN AS -USED-'')')

      RETURN
      END
