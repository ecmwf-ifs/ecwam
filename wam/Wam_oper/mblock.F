      SUBROUTINE MBLOCK (IA1, KA, KE, IPP)

C ----------------------------------------------------------------------
C
C**** *MBLOCK* - ROUTINE TO ARRANGE WAMODEL GRID FOR ONE BLOCK.
C
C     H.GUNTHER            ECMWF       04/04/1990
C
C*    PURPOSE.
C     -------
C
C       *MBLOCK* ARRANGES WAMODEL GRID FOR A BLOCK AND
C                COMPUTES VARIOUS MODEL CONSTANTS
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *MBLOCK (IA1, KA, KE, IPP)*
C          *IA1*     - TOPOGRAPHIC DATA.
C          *KA*      - NUMBER OF FIRST LAT IN BLOCK.
C          *KE*      - NUMBER OF LAST LAT IN BLOCK.
C          *IPP*     - NUMBER OF SEA POINTS PER LAT.
C
C     METHOD.
C     -------
C
C       THE LAND POINTS ARE REMOVED. ALL MODEL CONSTANTS WHICH ARE
C       GRID DEPENDENT ARE COMPUTED AND STORED IN THE COMMON BLOCKS.
C
C     EXTERNALS.
C     ----------
C
C       *ABORT1*     - TERMINATES PROCESSING.
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
#include "parameter_topo.h"
C
#include "comcinp.h"
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
      DIMENSION IA1(NGX, NGY), IPP(NGY)
C
C ----------------------------------------------------------------------
C
#include "parcons.h"
C
C ----------------------------------------------------------------------
C
C*    1. UPDATE BLOCK NUMBER AND INITIALIZES ARRAYS.
C        -------------------------------------------
C
 1000 CONTINUE
      IGL = IGL + 1
      IF (IGL.GT.NBLO) THEN
         WRITE (IU06,*) '**********************************************'
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '*        FATAL ERROR IN SUB. MBLOCK          *'
         WRITE (IU06,*) '*        ==========================          *'
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '* MORE BLOCKS THAN DIMENSION ALLOWS.         *'
         WRITE (IU06,*) '* BLOCK NUMBER IS                 IGL = ', IGL
         WRITE (IU06,*) '* DIMENSION IS                   NBLO = ', NBLO
         WRITE (IU06,*) '* NUMBER OF FIRST LATITUDE IS      KA = ', KA
         WRITE (IU06,*) '* NUMBER OF LAST  LATITUDE IS      KE = ', KE
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '*  PROGRAM WILL BE ABORTED                   *'
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '**********************************************'
         CALL ABORT1
      ENDIF
      DO 1001 IJ=1,NIBLO
         DEPTH(IJ,IGL) = 0.
         IXLG(IJ,IGL) = 0
         KXLT(IJ,IGL) = 0
 1001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    2. THE FIRST AND LAST BLOCK MUST CONTAIN MORE THAN 2
C*       ALL OTHER BLOCKS MORE  THAN 3 LATITUDES.
C        -------------------------------------------------
C
 2000 CONTINUE
      IF ((KE.EQ.1) .OR. (KA.EQ.NY) .OR.
     1    ((KA.NE.1) .AND. (KE.EQ.NY) .AND. (KE-KA.LT.2))) THEN
         WRITE (IU06,*) '**********************************************'
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '*        FATAL ERROR IN SUB. MBLOCK          *'
         WRITE (IU06,*) '*        ==========================          *'
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '* BLOCK LENGTH IS TO SHORT.                  *'
         WRITE (IU06,*) '* LESS THAN 2 LATITUDES IN FIRST OR LAST, OR *'
         WRITE (IU06,*) '* LESS THAN 3 LATITUDES IN OTHER BLOCKS.     *'
         WRITE (IU06,*) '* BLOCK NUMBER IS               IGL = ', IGL
         WRITE (IU06,*) '* BLOCK LENGTH IS             NIBLO = ', NIBLO
         WRITE (IU06,*) '* NUMBER OF FIRST LATITUDE IS    KA = ', KA
         WRITE (IU06,*) '* NUMBER OF LAST  LATITUDE IS    KE = ', KE
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '*  PROGRAM WILL BE ABORTED                   *'
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '**********************************************'
         CALL ABORT1
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    3. COMPUTE INDICES OF FIRST, SECOND, BEFORE LAST, AND LAST LAT.
C        -----------------------------------------------------------
C
 3000 CONTINUE
      IF (KA.EQ.1) THEN
         IJS (IGL) = 1
         IJL2(IGL) = IPP(1)
      ELSE
        IJS (IGL) = IPP(KA)+1
        IJL2(IGL) = IPP(KA)+IPP(KA+1)
      ENDIF
      IJLT(IGL) = 0
      DO 3001 K=KA,KE
         IJLT(IGL) = IJLT(IGL)+IPP(K)
 3001 CONTINUE
      IF (KE.EQ.NY) THEN
         IJL (IGL) = IJLT(IGL)
      ELSE
         IJL (IGL) = IJLT(IGL)-IPP(KE)
      ENDIF
      IJLS(IGL) = IJL(IGL)-IPP(KE-1)+1
C
C ----------------------------------------------------------------------
C
C*    4. REMOVE LAND POINTS AND STORE COS AND SIN OF LAT.
C        ------------------------------------------------
C
 4000 CONTINUE
      IP = 0
      DO 4001 K=KA,KE
         DO 4002 I=1,NLONRGG(K)
            IF (IA1(I,K).GT.-990) THEN
              IP = IP+1
              DEPTH(IP,IGL) = REAL(IA1(I,K))
              IXLG(IP,IGL) = I
              KXLT(IP,IGL) = K
            ENDIF
 4002    CONTINUE
 4001 CONTINUE
      IF (IP.NE.IJLT(IGL)) THEN
         WRITE (IU06,*) '**********************************************'
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '*        FATAL ERROR IN SUB. MBLOCK          *'
         WRITE (IU06,*) '*        ==========================          *'
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '* TOTAL NUMBER OF SEAPOINTS DO NOT MATCH.    *'
         WRITE (IU06,*) '* BLOCK NUMBER                    IGL = ', IGL
         WRITE (IU06,*) '* NO. OF SEAPOINTS COUNTED         IP = ', IP
         WRITE (IU06,*) '* NO. OF SEAPOINTS EXPECTED IJLT(IGL) = ',
     1                                            IJLT(IGL)
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '*  PROGRAM WILL BE ABORTED                   *'
         WRITE (IU06,*) '*                                            *'
         WRITE (IU06,*) '**********************************************'
         CALL ABORT1
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    5. PRINTER PROTOCOL OF BLOCK.
C        --------------------------
C
 5000 CONTINUE
      IF (IGL.EQ.1) THEN
         WRITE (IU06,'(1H0,'' BLOCKING INFORMATION:'')')
         WRITE (IU06,'(1H ,''            LATITUDES   '',
     1                     ''   SECOND LAT. INDEX '',
     2                     '' SECOND TO LAST LAT  '',
     3                     ''   TOTAL'')')
         WRITE (IU06,'(1H ,''  NO     SOUTH     NORTH'',
     1                     ''     START       END'',
     2                     ''     START       END'',
     3                     ''    POINTS'')')
      ENDIF
      WRITE (IU06,'(1X,I4,2F10.2,5I10)')
     1        IGL, AMOSOP+(KA-1)*XDELLA, AMOSOP+(KE-1)*XDELLA,
     2        IJS(IGL), IJL2(IGL), IJLS(IGL), IJL(IGL), IJLT(IGL)

      RETURN
      END
