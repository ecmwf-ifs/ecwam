      SUBROUTINE CLOSEND (IJS,IJL,IG,IGL,CDATE,CDATEWH,NEWREAD,NEWFILE,
     1                    U10OLD,THWOLD,U10NEW,THWNEW) 

C ----------------------------------------------------------------------
C
C**** *CLOSEND* - HANDLING OF WINDS OUTSIDE MULTITASKED AREA    
C
C     P.A.E.M. JANSSEN  KNMI/ECMWF  SEPTEMBER 1994
C     J. BIDLOT         ECMWF       FEBRUARY  1996  MESSAGE PASSING	
C
C*    PURPOSE.
C     --------
C
C       READ WINDS WHEN NEEDED.                                   
C
C**   INTERFACE.
C     ----------
C
C     *CALL* *CLOSEND (IJS, IJL, IG, IGL, CDATE,CDATEWH,NEWREAD,NEWFILE)*
C    1                  U10OLD,THWOLD,U10NEW,THWNEW) 
C      *IJS*    - INDEX OF FIRST GRIDPOINT
C      *IJL*    - INDEX OF LAST GRIDPOINT
C      *IG*     - BLOCK NUMBER
C      *IGL*    - NUMBER OF BLOCKS
C      *NEWREAD* - TRUE IF NEW WINDS HAVE BEEN READ
C      *NEWFILE* - TRUE IF NEW WIND FILE HAS BEEN OPENED
C      *U10NEW*    NEW WIND SPEED IN M/S.
C      *U10OLD*    INTERMEDIATE STORAGE OF MODULUS OF WIND
C                  VELOCITY.
C      *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
C                  NOTATION (POINTING ANGLE OF WIND VECTOR,
C                  CLOCKWISE FROM NORTH).
C      *THWOLD*    INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
C                  WIND VELOCITY.
C
C     METHOD.
C     -------
C
C       COPIES NEW TO OLD WINDS WHEN NEEDED. CLOSES FILES
C
C     EXTERNALS.
C     ---------
C
C       *CREWFN*    - CREATE A WIND FILE NAME.
C       *INCDATE*   - UPDATE DATE TIME GROUP.
C
C     REFERENCE.
C     ----------
C
C       NONE
C
C ----------------------------------------------------------------------
C
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.
C
#include "param.h"
C
#include "comstat.h"
C
#include "comtest.h"
C
#include "commpp.h"
C
#include "comunit.h"
C
#include "comwind.h"
C
C ----------------------------------------------------------------------
C
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10OLD,THWOLD
      REAL,DIMENSION(NINF:NSUP) :: U10NEW,THWNEW

C ----------------------------------------------------------------------
C
      LOGICAL NEWREAD, NEWFILE
C
      CHARACTER*12 CDATEWH, CDATE
C
C
C ----------------------------------------------------------------------
C
C*    1. SAVE WIND INTO INTERMEDIATE STORAGE.
C     ----------------------------------------
C
      IF (NEWREAD) THEN
         DO 1000 IJ=IJS,IJL
            U10OLD(IJ,IG) = U10NEW(IJ)
            THWOLD(IJ,IG) = THWNEW(IJ)
 1000    CONTINUE
         NEWREAD = .FALSE.
         IF (ITEST.GE.2) THEN
            IF (ITESTB.GE.IG) THEN
               WRITE(IU06,*) '       SUB. CLOSEND: STORE WIND ',
     1                       ' AT CDATE = ',CDATE
            ENDIF
         ENDIF
      ENDIF                                     
C
C ----------------------------------------------------------------------
C
C*    2. UPDATE WIND COUNTERS IF LAST BLOCK HAS BEEN DONE.
C        -------------------------------------------------
C
      IF (IG.EQ.IGL) THEN
         IF (NEWFILE) THEN
C
C*    2.1 A NEW WIND FILE HAS BEEN USED.
C         ------------------------------
C
C*    2.1.1 CLOSE WIND FILES SO THEY DO NOT ACCUMULATE.
C           -------------------------------------------
C
            CLOSE (UNIT=IU17, STATUS='DELETE')
            IF (ITEST.GE.2)
     1         WRITE(IU06,'(7X,''SUB. CLOSEND: CLOSED IU17 = '',
     2                      I8,'' AT CDATE = '',A12)') IU17, CDATE
C
C*    2.1.2 UPDATE WIND FILE TIME COUNTER AND UNITS.
C           ----------------------------------------
C
            IUH  = IU17
            IU17 = IU18
            IU18 = IUH
            NEWFILE = .FALSE.
            IDELWH = MAX(IDELWI,IDELPRO)
            CALL INCDATE(CDAWIFL,IDELWH)
            CALL INCDATE(CDATEFL,IDELWH)
            IF (ITEST.GE.2) THEN
               WRITE(IU06,*) '      SUB. CLOSEND: NEW WINDFILE DATES'
               WRITE(IU06,*) '        DATE OF NEXT WIND FILE IS ',
     1                       'CDAWIFL = ', CDAWIFL
               WRITE(IU06,*) '        FILE WILL BE ACCESSED  AT ',
     1                       'CDATEFL = ', CDATEFL
            ENDIF
         ENDIF
C
C*    2.2 UPDATE WIND FIELD COUNTER.
C         --------------------------
C
         CDATEWO=CDATEWH
      ENDIF
C
      RETURN
      END     
