      SUBROUTINE NOTIM  (CDTWIS, CDTWIE, NC, NR,
     1                   U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     2                   IUGB, IVGB, ILEN)

C ----------------------------------------------------------------------
C
!     MODIFIED  J. BIDLOT  FEBRARY 1996  MESSAGE PASSING

C**** *NOTIM* - STEERING MODULE IF NO TIME INTERPOLATION WANTED.
C
C*    PURPOSE.
C     --------
C
C       NOTIM NO TIME INTERPOLATION: PROCESS WINDFIELDS.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *NOTIM (CDTWIS, CDTWIE, NC, NR)*
C          *CDTWIS* - DATE OF FIRST WIND FIELD.
C          *CDTWIE* - DATE OF LAST FIRST WIND FIELD.
C          *NC*     - FIRST DIMENSION OF INPUT WIND FIELD.
C          *NR*     - SECOND DIMENSION OF INPUT WIND FIELD.
C
C     METHOD.
C     -------
C
C       NO TIME INTERPOLATION:
C       WINDFIELDS ARE PROCESSED EVERY IDELWI SECONDS (U,V),
C       THE WINDS INTERPOLATED IN SPACE ONLY (US,DS)
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
      CHARACTER*10 CDTWIE, CDTWIH, CDTWIS, CDT, ZERO
C
C ----------------------------------------------------------------------
      INTEGER IUGB(ILEN), IVGB(ILEN)
C ----------------------------------------------------------------------
C
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10OLD,THWOLD,USOLD,Z0OLD,TAUW
C ----------------------------------------------------------------------
C
      REAL, ALLOCATABLE :: UWND(:,:),VWND(:,:)
      REAL,DIMENSION(NINF:NSUP,NBLO) :: US,DS
C
C        *UWND*    REAL   INPUT WIND FIELD ARRAY (U COMPONENT).
C        *VWND*    REAL   INPUT WIND FIELD ARRAY (V COMPONENT).
C        *US*      REAL   OUTPUT WIND FIELD ARRAY (FRICTION VELOCITY).
C        *DS*      REAL   OUTPUT WIND FIELD ARRAY (DIRECTION).
C
C ----------------------------------------------------------------------
C
C*    1. INITIALISE WIND REQUEST DATE AND PROCESS FIRST WINDFIELD
C*       IF COLD START.
C        --------------------------------------------------------
C
 1000 CONTINUE

      IF(.NOT.ALLOCATED(UWND)) ALLOCATE(UWND(NC,NR))
      IF(.NOT.ALLOCATED(VWND)) ALLOCATE(VWND(NC,NR))

      ZERO = ' '
      CDTWIH = CDTWIS
      IF (CDA.EQ.ZERO) THEN
         CDA = CDTWIS
         CALL GETWND (U10OLD, THWOLD, CDA, UWND, VWND, NC, NR,
     1                IUGB, IVGB, ILEN)
         DO 1001 IG = 1,IGL
            CALL AIRSEA (U10OLD(IJS(IG),IG), TAUW(IJS(IG),IG),
     1                   USOLD(IJS(IG),IG), Z0OLD(IJS(IG),IG),
     2                   IJS(IG), IJL(IG))
 1001    CONTINUE
         IF (ITEST.GE.3) THEN
            WRITE(IU06,'(''       SUB. NOTIM: FIRST WIND FIELD '',
     1                ''SAVED IN COMMON WIND'')')
         ENDIF
         IF (CDA.EQ.CDTWIE) THEN
            IF(ALLOCATED(UWND)) DEALLOCATE(UWND)
            IF(ALLOCATED(VWND)) DEALLOCATE(VWND)
            RETURN
         ENDIF
         CALL INCDATE (CDTWIH,IDELWO)
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    2. CREATE WIND FILE NAME FROM END DATE.
C        ------------------------------------
C
 2000 CONTINUE
ccccc      CALL CREWFN (IUVELO, CDTWIE)
ccccc      IF (ITEST.GE.3) THEN
ccccc         WRITE(IU06,*) '       SUB. NOTIM: NEW WIND FILE '
ccccc         WRITE(IU06,*) '       UNIT IS IUVELO = ', IUVELO,
ccccc     1                 ' WIND FILE DATE IS CDTWIE = ', CDTWIE
ccccc      ENDIF
ccccc      REWIND (UNIT=IUVELO)
      MP = 0
C
C ----------------------------------------------------------------------
C
C*    3. LOOP OVER OUTPUT WIND TIMES.
C        ----------------------------
C
 3000 CONTINUE
C
C*    3.1 READ ONE WIND FIELD AND TRANSFORM TO BLOCKS.
C         --------------------------------------------
C
 3100 CONTINUE
      CALL GETWND (US, DS, CDTWIH, UWND, VWND, NC, NR,
     1             IUGB, IVGB, ILEN)
      MP = MP + 1
C
C*    3.2 SAVE BLOCKED WIND FIELD.
C         ------------------------
C
      IF (IDELPRO.GT.IDELWO) THEN
C
C*        MORE THAN ONE FIELD PER PROPAGATION TIME STEP.
C*        SAVE ON SCRATCH UNITS:
C
         DO 3201 IG=1,IGL
            IUNIT =IUSCR(IG)
            WRITE (IUNIT,ERR=6100,IOSTAT=IOS)
     1            CDTWIH, IG, (US(I,IG),I=IJS(IG),IJL(IG)),
     2            (DS(I,IG),I=IJS(IG),IJL(IG))
 3201    CONTINUE
         IF (ITEST.GE.3) THEN
            WRITE(IU06,*) '       SUB. NOTIM: NEW WIND FILES AT ',
     1                    'CDTWIH = ', CDTWIH,' WRITTEN TO SCRATCH UNIT'
            WRITE(IU06,'(''+'', A8)') IUNIT
         ENDIF
C
C*        UPDATE WIND FIELD REQUEST TIME.
C
         CALL INCDATE (CDTWIH,IDELWO)
C
C*        IF TIME LEFT BRANCH BACK TO 3.0 FOR NEXT FIELD.
C
         IF (CDTWIH.LE.CDTWIE) GOTO 3000
      ELSE
C
C*        WIND FIELD IS CONSTANT FOR ONE PROPAGATION TIME STEP:
C*        WRITE TO OUTPUT.
C         -----------------------------------------------------
C
         DO 3202 IG=1,IGL

cccc            WRITE (IUVELO,ERR=6400,IOSTAT=IOS)
cccc     1            CDTWIH, IG, (US(I,IG),I=IJS(IG),IJL(IG)),
cccc     2            (DS(I,IG),I=IJS(IG),IJL(IG))
ccc 
cc save to array instead
         CDTNEXT=CDTWIH
         DO I=IJS(IG),IJL(IG)
            USNEXT(I,IG)=US(I,IG)
            DSNEXT(I,IG)=DS(I,IG)
         ENDDO

cccccc

 3202    CONTINUE
C
C*        UPDATE WIND FIELD REQUEST TIME.
C         -------------------------------
C
         CALL INCDATE (CDTWIH,IDELWO)
C
C*        IF TIME LEFT BRANCH BACK TO 3.0 FOR NEXT FIELD.
C
cccccccc
cccccc         IF (CDTWIH.LE.CDTWIE) GOTO 3000
ccc we don't allow to loop anymore
         IF (CDTWIH.LE.CDTWIE) GOTO 6500

         IF (ITEST.GE.3) THEN
            WRITE(IU06,*) '       SUB. NOTIM: NEW WIND FILE GENERATED '
         ENDIF

ccccccccccccc        CLOSE (UNIT=IUVELO, STATUS='KEEP')

         IF(ALLOCATED(UWND)) DEALLOCATE(UWND)
         IF(ALLOCATED(VWND)) DEALLOCATE(VWND)

         RETURN

      ENDIF
C
C ----------------------------------------------------------------------
C
C*    4. RE-ARRANGE THE BLOCKS.
C        ----------------------
C
 4000 CONTINUE
      DO 4001 IG=1,IGL
         REWIND (UNIT=IUSCR(IG))
 4001 CONTINUE

      MSTEP = IDELPRO/IDELWO
      DO 4002 M=1,MP,MSTEP
         DO 4003 IG=1,IGL
            IUNIT = IUSCR(IG)
            DO 4004 MM=1,MSTEP
               READ (IUNIT,ERR=6200,IOSTAT=IOS)
     1            CDT,IGG, (US(I,IG),I=IJS(IG),IJL(IG)),
     2            (DS(I,IG),I=IJS(IG),IJL(IG))
C
               WRITE (IUVELO,ERR=6300,IOSTAT=IOS)
     1            CDT, IGG, (US(I,IG),I=IJS(IG),IJL(IG)),
     2            (DS(I,IG),I=IJS(IG),IJL(IG))
 4004       CONTINUE
 4003    CONTINUE
 4002 CONTINUE
      REWIND (UNIT=IUVELO)
      DO 4005 IG=1,IGL
         CLOSE (UNIT=IUSCR(IG), STATUS='DELETE')
 4005 CONTINUE
      CLOSE (UNIT=IUVELO, STATUS='KEEP')
      IF (ITEST.GE.3) THEN
         WRITE(IU06,*) '       SUB. NOTIM: NEW WIND FILE ',
     1                 'GENERATED AS COPY FROM  SCRATCH UNITS'
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
         WRITE (IU06,*) ' *       FATAL ERROR IN SUB. NOTIM:         *'
         WRITE (IU06,*) ' *       ==========================         *'
         WRITE (IU06,*) ' * ERROR WHEN WRITTING ON SCRATCH UNIT      *'
         WRITE (IU06,*) ' * FOR BLOCK NUMBER          IG = ', IG
         WRITE (IU06,*) ' * NUMBER OF WIND FIELD IS   MP = ', MP
         WRITE (IU06,*) ' * DATE OF WIND FIELD IS CDTWIH = ', CDTWIH
         WRITE (IU06,*) ' * MESSAGE IS               IOS = ', IOS
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS         *'
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' ********************************************'
         CALL ABORT1
 6200    CONTINUE
         WRITE (IU06,*) ' '
         WRITE (IU06,*) ' ********************************************'
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' *       FATAL ERROR IN SUB. NOTIM:         *'
         WRITE (IU06,*) ' *       ==========================         *'
         WRITE (IU06,*) ' * ERROR WHEN READING FROM SCRATCH UNIT     *'
         WRITE (IU06,*) ' * FOR BLOCK NUMBER       IG = ', IG
         WRITE (IU06,*) ' * DATE OF WIND FIELD IS CDT = ', CDT
         WRITE (IU06,*) ' * MESSAGE IS            IOS = ', IOS
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS         *'
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' ********************************************'
         CALL ABORT1
 6300 CONTINUE
         WRITE (IU06,*) ' '
         WRITE (IU06,*) ' ********************************************'
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' *       FATAL ERROR IN SUB. NOTIM:         *'
         WRITE (IU06,*) ' *       ==========================         *'
         WRITE (IU06,*) ' * ERROR WHEN WRITTING ON WIND OUTPUT FILE  *'
         WRITE (IU06,*) ' * AFTER READING FROM SCRATCH UNIT          *'
         WRITE (IU06,*) ' * FOR BLOCK NUMBER       IG = ', IG
         WRITE (IU06,*) ' * DATE OF WIND FIELD IS CDT = ', CDT
         WRITE (IU06,*) ' * MESSAGE IS            IOS = ', IOS
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS         *'
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' ********************************************'
         CALL ABORT1
 6400 CONTINUE
         WRITE (IU06,*) ' '
         WRITE (IU06,*) ' ********************************************'
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' *       FATAL ERROR IN SUB. NOTIM:         *'
         WRITE (IU06,*) ' *       ==========================         *'
         WRITE (IU06,*) ' * ERROR WHEN WRITTING ON WIND OUTPUT FILE  *'
         WRITE (IU06,*) ' * SCRATCH UNITS ARE NOT USED               *'
         WRITE (IU06,*) ' * FOR BLOCK NUMBER          IG = ', IG
         WRITE (IU06,*) ' * DATE OF WIND FIELD IS CDTWIH = ',CDTWIH
         WRITE (IU06,*) ' * MESSAGE IS               IOS = ', IOS
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS         *'
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' ********************************************'
         CALL ABORT1
 6500 CONTINUE
         WRITE (IU06,*) ' '
         WRITE (IU06,*) ' ********************************************'
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' *       FATAL ERROR IN SUB. NOTIM:         *'
         WRITE (IU06,*) ' *       ==========================         *'
         WRITE (IU06,*) ' * ERROR WHEN WRITTING ON WIND OUTPUT FILE  *'
         WRITE (IU06,*) ' * MORE THAN ONE LOOP NOT ALLOWED ANY LONGER*'
         WRITE (IU06,*) ' * WHEN USING ARRAY INSTEAD OF WRITE TO FILE*'
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS         *'
         WRITE (IU06,*) ' *                                          *'
         WRITE (IU06,*) ' ********************************************'
         CALL ABORT1
      END
