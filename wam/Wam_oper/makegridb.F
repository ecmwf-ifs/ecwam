      SUBROUTINE MAKEGRIDB (BLOCK, GRID, IJS, IJL, IG,
     &                     NGX, NGY, NIBLO, NBLO,
     &                     IXLG, KXLT)

C ----------------------------------------------------------------------
C
C**** *MAKEGRID* - MAKE GRIDDED WAM MODEL FIELDS FROM BLOCKED FIELDS.
C
C     H. GUNTHER       ECMWF    NOVEMBER 1989
C
C*    PURPOSE.
C     --------
C
C       GRIDDED WAVE FIELDS ARE CREATED FROM BLOCKED WAVE FIELDS.
C
C**   INTERFACE.
C     ----------
C
C        *CALL MAKEGRID*
C           *BLOCK*   REAL   DATA IN BLOCKED FORMAT
C           *GRID*    REAL   DATA IN GRID FORMAT
C           *IJS*     INTEGER INDEX OF FIRST GRIDPOINT IN BLOCK
C           *IJL*     INTEGER INDEX OF LAST  GRIDPOINT IN BLOCK
C           *IG*      INTEGER BLOCK NUMBER
C
C     METHOD.
C     -------
C
C       THE PARAMETER, WHICH IS GIVEN IN BLOCKED FORMAT, IS
C       DISTRIBUTED IN GRID FORMAT. ONE BLOCK IS DONE IN ONE
C       CALL. BEFORE THE FIRST BLOCK IS TRANSFORMED THE GRID ARRAY IS
C       INITIALISED WITH ZERO.
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
      PARAMETER(ZMISS=-999.)
C
C ----------------------------------------------------------------------
C
      DIMENSION BLOCK(NIBLO), GRID(NGX,NGY)
      DIMENSION IXLG(NIBLO,NBLO), KXLT(NIBLO,NBLO)
      LOGICAL LICERUN ; DATA LICERUN / .FALSE. /
C
C ----------------------------------------------------------------------
C
C*    1. INITIALIZATION
C     -----------------
C
      IF (IG.EQ.1) THEN
         DO J = 1,NGY
           DO I = 1,NGX
             GRID(I,J) = ZMISS
           ENDDO
         ENDDO
      ENDIF
C
C*    2. MAKE GRIDDED FIELD
C     ---------------------
C
      DO 2000 IJ = IJS,IJL
         IX = IXLG(IJ,IG)
         IY = NGY- KXLT(IJ,IG) +1
         GRID(IX,IY) = BLOCK(IJ)
 2000 CONTINUE
C
C*    3. SET ICE POINTS TO MISSING DATA.
C     ----------------------------------
C
C     IF (LICERUN) THEN
C        DO 3100 KICE = 1,IGNO
C           IF (IGICE(KICE).EQ.IG) THEN
C              DO 3000 I = KSTART(KICE),KEND(KICE)
C                 IX = IXLG(ICEIJ(I),IG)
C                 IY = NGY- KXLT(ICEIJ(I),IG) +1
C                 GRID(IX,IY) = ZMISS
 3000          CONTINUE
C           ENDIF
 3100    CONTINUE
C     ENDIF

      RETURN
      END
