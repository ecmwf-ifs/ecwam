      SUBROUTINE MPDECOMP(NPR,NSTART,NEND,KLENTOP,KLENBOT,MAXLEN)
C
C****  *MPDECOMP* - EVEN DECOMPOSITION OF THE GRID DOMAIN
C****               AMONG PROCESSES
C
!     J. BIDLOT    ECMWF   MARCH 1996  MESSAGE PASSING
!     J. BIDLOT    ECMWF   JANUARY 1998 introduce NPR 
C
C     PURPOSE.
C     --------
C     DETERMINE AN EVEN DECOMPOSITION OF THE GRID ARRAYS FOR USE ON
C     A DISTRIBUTED MEMORY COMPUTER USING A MESSAGE PASSING PROTOCOL
C     FOR THE EXCHANGE OF INFORMATION ACROSS THE DIFFERENT PROCESSES
C     THE LENGTH OF THE MESSAGES BETWEEN DIFFERENT PROCESS IS ALSO
C     COMPUTED.
C     IT IS ASSUMED THAT EACH GRID SUB DOMAIN WILL ONLY SEND OR
C     RECEIVE MESSAGES FROM THE PREVIOUS AND NEXT PROCESS 
C
C*    INTERFACE.
C     ----------
C     CALL *MPDECOMP*(NPR,NSTART,NEND,KLENTOP,KLENBOT,MAXLEN)
C
C      *NPR*       NUMBER OF SUBDOMAINS (USUALLY SET TO BE THE NUMBER OF PE'S 
C      *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
C      *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
C      *KLENTOP*   LENGTH OF THE MESSAGE SENT BY THE NEXT
C                  SUB GRID DOMAIN (IN GRID POINT UNIT)
C      *KLENBOT*   LENGTH OF THE MESSAGE SENT BY THE PREVIOUS
C                  SUB GRID DOMAIN (IN GRID POINT UNIT)
C      *MAXLEN*    MAXIMUM NUMBER OF POINTS IN ANY SUB DOMAIN

C     METHOD.
C     -------
C     SINCE THE GLOBAL GRID DOMAIN IS MAPPED ONTO A ONE DIMENSIONAL ARRAY
C     FOLLOWING INCREASING LATITUDE LINES, SUB GRID DOMAINS ARE SEGMENTS
C     OF THAT LONG ARRAY. THE LENGTH OF THOSE SEGMENTS IS DETERMINED BY
C     DISTRIBUTING AS EVENLY AS POSSIBLE THE TOTAL NUMBER OF POINTS.
C
C     THE ADVECTION SCHEME IS SUCH THAT IT ONLY USES
C     NEIGHBORING GRID POINTS IN THE 2D GRID. THEREFORE,(see figure) 
C     FOR A GIVEN SEGMENT (PROCESS) INFORMATION FROM
C     PREVIOUS GRID POINTS ARE ONLY NEEDED FROM FEW GRID POINTS
C     STARTING FROM THE GRID POINT LOCATED, IN THE 2-D GRID,
C     JUST BELOW (ON THE PREVIOUS LATITUDE LINE) THE STARTING POINT
C     OF THAT PARTICULAR SEGMENT PROVIDED THIS PREVIOUS GRID POINT
C     IS A VALID SEA POINT. IF NOT A VALID GRID POINT, FIND THE NEXT ONE
C     WHICH IS A VALID SEA POINT BUT DOES NOT BELONG TO THE SEGEMENT IN QUESTION
C     NATURALLY, THE FIRST SEGMENT DOES NOT REQUIRE ANY INFORMATION
C     FROM PREVIOUS SEGMENT. 
C
C
C
C     FIGURE: 2-D GRID
C     -------
C
C       +    +    +    +    +    +    +    +
C       +    +    +    +    +    +    +    +
C       +    +    +    I    I    I    I    I
C       I    I    I    *    *    *    *    *
C       *    *    *    *    *    *    *    *
C       *    *    *    *    *    *    *    *
C       *    *    *    *    *    *    *    *
C       *    *    *    *    *    *    *    *
C       *    *    *    *    *    *    I    I
C       I    I    I    I    I    I    +    +
C       +    +    +    +    +    +    +    +
C       +    +    +    +    +    +    +    +

C     FIGURE LEGEND :
C     --------------
C     *  GRID POINTS OF THE SELECTED SEGMENT (PROCESS)
C     I  GRID POINTS OF THE PREVIOUS/FOLLOWING SEGEMENT FOR WHICH
C        INFORMATION IS NECESSARY IF THE ADVECTION SCHEME IS ONLY COMPUTED
C        AT THE * POINTS
C
C     +  OTHER GRID POINTS THAT DO NOT AFFECT THE * POINTS 

C
C     SIMILARLY, INFORMATION FROM FOLLOWING GRID POINTS ARE ONLY NECESSARY
C     UP TO THE GRID POINT JUST ABOVE (IN THE 2D GRID) THE LAST POINT
C     OF THE SELECTED SEGMENT, PROVIDED THIS FOLLOWING GRID POINT IS A
C     VALID SEA POINT, OTHERWISE, FIND THE FIRST PREVIOUS ONE WHICH IS
C     A VALID ONE BUT DOES NOT BELONG TO THE SEGMENT IN QUESTION.
C
C     EXTERNALS.
C     ----------
C          NONE
C
C     REFERENCES.
C     -----------
C          NONE
C
C -------------------------------------------------------------------
C
#include "param.h"
C
#include "comcoup.h"
C
#include "comtest.h"
C
#include "comfred.h"
C
#include "comgrid.h"
C
#include "comubuf.h"
C
#include "txtubuf.h"
C
#include "commesspass.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
!----------------------------------------------------------------------
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      INTEGER,DIMENSION(NPR) :: NSTART,NEND,KLENTOP,KLENBOT

C ----------------------------------------------------------------------
C
      INTEGER,ALLOCATABLE :: ICOMBUF(:)

C
C     0. READ THE ACTUAL LENGTH OF THE READ ARRAYS
C        -----------------------------------------
        IF(LWCOU) THEN
          IU07=29
        ELSE
          IU07=7
        ENDIF
        REWIND (UNIT=IU07)
        READ(IU07) FR, DFIM, GOM, C,
     1             DELTH, DELTR, TH, COSTH, SINTH
        READ(IU07) DELPHI, DELLAM, NLONRGG, SINPH, COSPH,
     1             IGL, IJS, IJL2, IJLS, IJL, IJLT
C
        REWIND (UNIT=IU07)

C
C*    1. INPUT NEIGHBOURING GRID POINT INDICES 
C        -------------------------------------
C
        IF(LWCOU) THEN
          IU08=31
        ELSE
          IU08=8
        ENDIF

        IF(LMESSPASS) THEN

          IREAD=1
          KTAG=1

C         READ IOU8 ON PE IREAD
          IF(IRANK.EQ.IREAD) THEN
            READ (IU08) KLAT, KLON
            REWIND (UNIT=IU08)
          ENDIF

C         SEND KLAT AND KLON TO OTHER PE'S
          IF(NPR.GT.1) THEN
            ITAG=KTAG
            MPLENGTH=4*NIBLO
            ALLOCATE(ICOMBUF(MPLENGTH))

            IF(IRANK.EQ.IREAD) THEN
              KCOUNT=0
              DO IC=1,2
                DO IJ=1,NIBLO
                 KCOUNT=KCOUNT+1
                 ICOMBUF(KCOUNT)=KLAT(IJ,IC)
                END DO
              END DO
              DO IC=1,2
                DO IJ=1,NIBLO
                   KCOUNT=KCOUNT+1
                   ICOMBUF(KCOUNT)=KLON(IJ,IC)
                END DO
              END DO

              CALL MPE_BROADCAST(ICOMBUF,MPLENGTH,1,IRANK,
     &                           ITAG,0,0,0,IERR)
              IF(IERR.LT.0) CALL MPEI_ABORT
     1                    ('MPE_BROADCAST ERROR in MPDECOMP' )
            ELSE

              CALL MPE_RECV(ICOMBUF,MPLENGTH,1,IREAD,ITAG,0,0,0,
     $                      KRCOUNT,KRFROM,KRTAG,IERR)
              IF(IERR.LT.0) CALL MPEI_ABORT
     1                         ('MPE_RECV ERROR in MPDECOMP' )
              IF(KRCOUNT.NE.MPLENGTH) CALL MPEI_ABORT
     1        ('MPE_RECV ERROR in MPDECOMP:MISMATCHED MSG LENGTH')
              IF(KRTAG.NE.ITAG) CALL MPEI_ABORT
     1        ('MPE_RECV ERROR in MPRECOMP MISMATCHED TAGS' )


              KCOUNT=0
              DO IC=1,2
                DO IJ=1,NIBLO
                  KCOUNT=KCOUNT+1
                  KLAT(IJ,IC)=ICOMBUF(KCOUNT)
                END DO
              END DO
              DO IC=1,2
                DO IJ=1,NIBLO
                  KCOUNT=KCOUNT+1
                  KLON(IJ,IC)=ICOMBUF(KCOUNT)
                END DO
              END DO

            ENDIF

            DEALLOCATE(ICOMBUF)

          ENDIF

        ELSE IF(.NOT.LMESSPASS) THEN
          READ (IU08) KLAT, KLON
          REWIND (UNIT=IU08)
        ENDIF



C*    2. FIND THE NUMBER OF POINTS PER PROCESS, THE START AND END INDEX
C        --------------------------------------------------------------
      NMEAN=IJL(NBLO)/NPR
      NREST=IJL(NBLO)-NMEAN*NPR
      
      NSTART(1)=1
      IF(NREST.GT.0) THEN
        NPTS=NMEAN+1
        NREST=NREST-1 
      ELSE
        NPTS=NMEAN
      ENDIF
      NEND(1)=NSTART(1)+NPTS-1
c
      DO IP=2,NPR
          NSTART(IP)=NSTART(IP-1)+NPTS
          IF(NREST.GT.0) THEN
            NPTS=NMEAN+1
            NREST=NREST-1 
          ELSE
            NPTS=NMEAN
          ENDIF
         NEND(IP)=NSTART(IP)+NPTS-1
       END DO
c
c
C     3. DETERMINE THE LENGTH OF THE MESSAGE THAT WILL BE EXCHANGED 
C        BETWEEN NEIGHBORING SUB GRID DOMAIN
C        -----------------------------------------------------------
c
         KLENBOT(1)=0
      DO IP=2,NPR
         I=-1
20       I=I+1
         ISB=KLAT(NSTART(IP)+I,1)
         IF(ISB.LE.0) GOTO 20
c
         KLENBOT(IP)=MAX(0,NSTART(IP)-ISB)
      END DO

         KLENTOP(NPR)=0
      DO IP=1,NPR-1
         I=-1
40       I=I+1
         IET=KLAT(NEND(IP)-I,2)
         IF(IET.LE.0) GOTO 40
c
         KLENTOP(IP)=MAX(0,IET-NEND(IP))
      END DO
C
      MAXLEN=0
      DO IP=1,NPR
        MAXLEN=MAX(MAXLEN,NEND(IP)-NSTART(IP)+1)
      ENDDO
C
      RETURN
      END
