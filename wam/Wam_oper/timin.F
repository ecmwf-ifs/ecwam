      SUBROUTINE TIMIN (CDTWIS, CDTWIE, NC, NR,
     1                   U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     2                   IUGB, IVGB, ILEN)

C ----------------------------------------------------------------------
!     MODIFIED  J. BIDLOT  FEBRARY 1996  MESSAGE PASSING

C
C**** *TIMIN* - STEERING MODULE IF TIME INTERPOLATION IS WANTED.
C
C     H. GUNTHER    ECMWF    MAY 1990     MODIFIED FOR SUB VERSION.
C     H. GUNTHER    ECMWF    DECEMBER 90  MODIFIED FOR CYCLE_4.
C
C*    PURPOSE.
C     --------
C
C       TIME INTERPOLATION: PROCESS WINDFIELDS.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *TIMIN (CDTWIS, CDTWIE, NC, NR)*
C          *CDTWIS* - DATE OF FIRST WIND FIELD
C          *CDTWIE* - DATE OF LAST FIRST WIND FIELD
C          *NC*     - NUMBER OF COLUMNES IN INPUT WIND
C          *NR*     - NUMBER OF ROWS     IN INPUT WIND
C
C     METHOD.
C     -------
C
C       WINDFIELDS ARE READ IN EVERY IDELWI SECONDS (U,V),
C       INTERPOLATED IN SPACE, BLOCKED AND SAVED ON SCRATCH UNITS.
C       MAGNITUDE AND DIRECTION ARE INTERPOLATED LINEARLY IN TIME AND
C       WRITTEN TO THE OUTPUT FILES.
C
C     EXTERNALS.
C     ----------
C
C       *ABORT1*     - TERMINATES PROCESSING.
C       *AIRSEA*    - TOTAL STRESS IN SURFACE LAYER.
C       *CREWFN*    - ASSIGN A FILE NAME.
C       *GETWND*    - READ A WINDFIELD (UWND,VWND) AND COMPUTE WIND
C                     FOR ALL BLOCKS (US,DS).
C       *INCDATE*   - INCREMENT DATE.
C       *LOCINT*    - SPACE INTERPOLATION.
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
#include "comcoup.h"
C
#include "comgrid.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
#include "comstat.h"
C
#include "comtest.h"
C
#include "comunit.h"
C
#include "comwind.h"
C
      CHARACTER*10 CDTWIE, CDTH, CDT1, CDT2, CDTWIS, ZERO
      LOGICAL LWNDFILE, LCLOSEWND
C
C ----------------------------------------------------------------------
      INTEGER IUGB(ILEN), IVGB(ILEN)
C ----------------------------------------------------------------------
C
C     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10OLD,THWOLD,USOLD,Z0OLD,TAUW
C ----------------------------------------------------------------------
C
#include "parcons.h"
C
      REAL, ALLOCATABLE :: UWND(:,:),VWND(:,:)

      REAL,DIMENSION(NINF:NSUP,NBLO) :: US,DS
      REAL,DIMENSION(NINF:NSUP,NBLO) :: US2,DS2
      REAL,DIMENSION(NINF:NSUP) :: US3,DS3
C
C        *UWND*   REAL   INPUT WIND FIELD ARRAY (U COMPONENT)
C        *VWND    REAL   INPUT WIND FIELD ARRAY (V COMPONENT)
C        *US*     REAL   OUTPUT WIND FIELD (FRICTION VELOCITIES).
C        *DS*     REAL   OUTPUT WIND FIELD (DIRECTIONS).
C
C ----------------------------------------------------------------------
C
C*    1. INITIALIZE TIMECOUNTER.
C        -----------------------
C
 1000 CONTINUE

      IF(.NOT.ALLOCATED(UWND)) ALLOCATE(UWND(NC,NR))
      IF(.NOT.ALLOCATED(VWND)) ALLOCATE(VWND(NC,NR))
 
      ZERO= ' '
C
      LCLOSEWND=.FALSE.
      IF (LWCOU) THEN
        LWNDFILE=.FALSE.
        ILEV=2
      ELSE
        LWNDFILE=.TRUE.
        ILEV=1
      ENDIF
C
      IF (CDA.EQ.ZERO) THEN
         CDT1 = CDTWIS
         CALL GETWND (US, DS, NINF, NSUP, NSUP, CDT1, UWND, VWND, NC,NR,
     1                IUGB, IVGB, ILEN, LWNDFILE, LCLOSEWND)
         DO IG=1,IGL
           DO IJ=IJS(IG),IJL(IG)
             U10OLD(IJ,IG) = US(IJ,IG)
             THWOLD(IJ,IG) = DS(IJ,IG)
           ENDDO
         ENDDO
         DO 1002 IG = 1,IGL
            CALL AIRSEA (U10OLD(IJS(IG),IG), TAUW(IJS(IG),IG),
     1                   USOLD(IJS(IG),IG), Z0OLD(IJS(IG),IG),
     2                   IJS(IG), IJL(IG), ILEV)

 1002    CONTINUE
         IF (ITEST.GE.3) THEN
            WRITE(IU06,'(''       SUB. TIMIN: FIRST WIND FIELD '',
     1                ''SAVED IN COMMON WIND'')')
         ENDIF
      ELSE
         CDT1 = CDA
         DO IG=1,IGL
           DO IJ=IJS(IG),IJL(IG)
             US(IJ,IG) = U10OLD(IJ,IG)
             DS(IJ,IG) = THWOLD(IJ,IG)
           ENDDO
         ENDDO
         CDTH = CDT1
         CALL INCDATE (CDTH,IDELWO)
         IF (CDTWIS.NE.CDTH) THEN
            WRITE(IU06,*) ' *******************************************'
            WRITE(IU06,*) ' *                                         *'
            WRITE(IU06,*) ' *        FATAL ERROR IN --TIMIN--         *'
            WRITE(IU06,*) ' *        =========================        *'
            WRITE(IU06,*) ' * DATES DO NOT MATCH.                     *'
            WRITE(IU06,*) ' * START DATE FOR WIND IS       CDTWIS = ',
     1                        CDTWIS
            WRITE(IU06,*) ' * LAST DATE SAVED IN COM WIND IS CDT1 = ',
     1                       CDT1
            WRITE(IU06,*) ' * PROCESSING WILL BE ABORTED              *'
            WRITE(IU06,*) ' *                                         *'
            WRITE(IU06,*) ' *******************************************'
            CALL ABORT1
         ENDIF
      ENDIF
      CDT2 = CDT1
      CALL INCDATE(CDT2,IDELWI)

      NTS = IDELWI/IDELWO
      DEL = REAL(IDELWO)/REAL(IDELWI)
      CALL CREWFN (IUVELO, CDTWIE)
      IF (ITEST.GE.3) THEN
         WRITE(IU06,*) '       SUB. TIMIN: NEW WIND FILE '
         WRITE(IU06,*) '       UNIT IS IUVELO = ', IUVELO,
     1                 ' WIND FILE DATE IS CDTWIE = ', CDTWIE
      ENDIF
      REWIND (UNIT=IUVELO)
C
C*    2. LOOP OVER INPUT WINDFIELDS.
C        ---------------------------
C
 2000 CONTINUE
C
C*    2.1 READ ONE WINDFIELD AND TRANSFORM ALL BLOCKS.
C         -------------------------------------------
C
      CDT2 = CDT1
      CALL INCDATE(CDT2,IDELWI)
      CALL GETWND (US2, DS2, NINF, NSUP, NSUP, CDT2, UWND, VWND, NC, NR,
     1             IUGB, IVGB, ILEN, LWNDFILE, LCLOSEWND)
C
C*    2.2 SAVE BLOCKED WIND FIELD ON SCRATCH UNITS.
C         -----------------------------------------
C
      DO 2201 IG=1,IGL
         IUNIT =IUSCR(IG)
         CDTH = CDT1
         DO 2202 N=1,NTS
            CALL INCDATE(CDTH,IDELWO)
            DO 2203 IJ=IJS(IG),IJL(IG)
               US3(IJ) = US(IJ,IG)+REAL(N)*DEL*(US2(IJ,IG)-US(IJ,IG))
               D = DS2(IJ,IG)-DS(IJ,IG)
               IF (ABS(D).GT.PI) D = D-ZPI*SIGN(1.,D)
               D = DS(IJ,IG)+REAL(N)*DEL*D
               DS3(IJ) = MOD(D+ZPI,ZPI)
 2203       CONTINUE
            WRITE (IUNIT,ERR=6100,IOSTAT=IOS)
     1             CDTH, IG, (US3(I),I=IJS(IG),IJL(IG)),
     2             (DS3(I),I=IJS(IG),IJL(IG))
 2202    CONTINUE
         DO 2204 IJ=IJS(IG),IJL(IG)
            US(IJ,IG) = US2(IJ,IG)
            DS(IJ,IG) = DS2(IJ,IG)
 2204    CONTINUE
 2201 CONTINUE
      CDT1 = CDT2
      IF (ITEST.GE.3) THEN
         WRITE(IU06,*) '       SUB. TIMIN: LAST WIND FIELD AT ',
     1                 'CDTH = ', CDTH,' WRITTEN TO SCRATCH UNITS'
      ENDIF
C
C*    2.3 UPDATE WIND FIELD REQUEST TIME AND READ NEXT IF REQUESTED.
C         ----------------------------------------------------------
C
         CALL INCDATE (CDTH,IDELWI)
         IF (CDTH.LE.CDTWIE) GOTO 2000
C
C ----------------------------------------------------------------------
C
C*    3. REWIND SCRATCH UNITS.
C        ---------------------
C
 3000 CONTINUE
      DO 3001 IG=1,IGL
         IUNIT = IUSCR(IG)
         REWIND (UNIT=IUNIT)
 3001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    4. RE-ARRANGE THE BLOCKS.
C        ----------------------
C
 4000 CONTINUE
      MP = MAX(IDELWI,IDELPRO)/IDELWO
      MSTEP = MAX(1,IDELPRO/IDELWO)
C
C*    4.1 LOOP OVER WIND FILES.
C         ---------------------
C
      DO 4101 M=1,MP,MSTEP
C
C*    4.1.1 LOOP OVER BLOCKS.
C           -----------------
C
         DO 4102 IG=1,IGL
            IUNIT = IUSCR(IG)
C
C*    4.1.1.1 LOOP OVER WIND FIELDS.
C             ----------------------
C
            DO 4103 MM=1,MSTEP
               READ (IUNIT,END=6200,IOSTAT=IOS)
     1            CDTH, IGG, (US(I,IG),I=IJS(IG),IJL(IG)),
     2             (DS(I,IG),I=IJS(IG),IJL(IG))
               IF (CDTH.GE.CDTWIS)
     1            WRITE (IUVELO,ERR=6300,IOSTAT=IOS)
     2            CDTH, IGG, (US(I,IG),I=IJS(IG),IJL(IG)),
     3             (DS(I,IG),I=IJS(IG),IJL(IG))
 4103       CONTINUE
 4102    CONTINUE
 4101 CONTINUE
      DO 4104 IG=1,IGL
         CLOSE (UNIT=IUSCR(IG) , STATUS='DELETE')
 4104 CONTINUE
      CLOSE (UNIT=IUVELO, STATUS='KEEP')
      IF (ITEST.GE.3) THEN
         WRITE(IU06,*) '       SUB. TIMIN: LAST WIND FIELD AT ',
     1                 'CDTH = ', CDTH,' WRITTEN TO OUTPUT UNIT'
      ENDIF

      IF(ALLOCATED(UWND)) DEALLOCATE(UWND)
      IF(ALLOCATED(VWND)) DEALLOCATE(VWND)

      RETURN
C
C ----------------------------------------------------------------------
C
C*    6. ERROR MESSAGES
C
 6100 CONTINUE
         WRITE (IU06,*) ' '
         WRITE (IU06,*) ' ********************************************'
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' *       FATAL ERROR IN SUB. TIMIN:         *'
         WRITE (IU06,*) ' *       ==========================         *'
         WRITE (IU06,*) ' * ERROR WHEN WRITTING ON SCRATCH UNIT      *'
         WRITE (IU06,*) ' * FOR BLOCK NUMBER        IG = ', IG
         WRITE (IU06,*) ' * NUMBER OF WIND FIELD IS  N = ', N
         WRITE (IU06,*) ' * DATE OF WIND FIELD IS CDTH = ', CDTH
         WRITE (IU06,*) ' * MESSAGE IS             IOS = ', IOS
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS         *'
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' ********************************************'
         CALL ABORT1
 6200    CONTINUE
         WRITE (IU06,*) ' '
         WRITE (IU06,*) ' ********************************************'
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' *       FATAL ERROR IN SUB. TIMIN:         *'
         WRITE (IU06,*) ' *       ==========================         *'
         WRITE (IU06,*) ' * ERROR WHEN READING FROM SCRATCH UNIT     *'
         WRITE (IU06,*) ' * FOR BLOCK NUMBER        IG = ', IG
         WRITE (IU06,*) ' * DATE OF WIND FIELD IS CDTH = ', CDTH
         WRITE (IU06,*) ' * MESSAGE IS             IOS = ', IOS
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS         *'
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' ********************************************'
         CALL ABORT1
 6300 CONTINUE
         WRITE (IU06,*) ' '
         WRITE (IU06,*) ' ********************************************'
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' *       FATAL ERROR IN SUB. TIMIN:         *'
         WRITE (IU06,*) ' *       ==========================         *'
         WRITE (IU06,*) ' * ERROR WHEN WRITTING ON WIND OUTPUT FILE  *'
         WRITE (IU06,*) ' * FOR BLOCK NUMBER       IG  = ', IG
         WRITE (IU06,*) ' * DATE OF WIND FIELD IS CDTH = ', CDTH
         WRITE (IU06,*) ' * MESSAGE IS            IOS  = ', IOS
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS         *'
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' ********************************************'
         CALL ABORT1

      END
