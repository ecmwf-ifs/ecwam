      SUBROUTINE NEWWIND (IJS,IJL,IG,IGL,CDATE,CDATEWH,NEWREAD,NEWFILE,
     1                    U10OLD,THWOLD,U10NEW,THWNEW)
C ----------------------------------------------------------------------
C
C**** *NEWWIND* - HANDLING OF WINDS OUTSIDE MULTITASKED AREA    
C
C     P.A.E.M. JANSSEN  KNMI/ECMWF  SEPTEMBER 1994
!     J. BIDLOT ECMWF FEBRUARY 1996 MESSAGE PASSING 
C
C*    PURPOSE.
C     --------
C
C       READ WINDS WHEN NEEDED.                                   
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *NEWWIND (IJS, IJL, IG, IGL, CDATE,NEWREAD, NEWFILE)
C     1                    U10OLD,THWOLD,U10NEW,THWNEW)*
C      *IJS*     - INDEX OF FIRST GRIDPOINT
C      *IJL*     - INDEX OF LAST GRIDPOINT
C      *IG*      - BLOCK NUMBER
C      *IGL*     - NUMBER OF BLOCKS
C      *CDATE*   - START DATE OF SOURCE FUNCTION INTEGRATION
C      *NEWREAD* - TRUE IF NEW WINDS HAVE BEEN READ
C      *NEWFILE* - TRUE IF NEW WIND FILE HAS BEEN OPENED
C      *U10NEW*    REAL      NEW WIND SPEED IN M/S.
C      *U10OLD*    REAL      INTERMEDIATE STORAGE OF MODULUS OF WIND
C                            VELOCITY.
C      *THWNEW*    REAL      WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
C                            NOTATION (POINTING ANGLE OF WIND VECTOR,
C                            CLOCKWISE FROM NORTH).
C      *THWOLD*    REAL      INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
C                            WIND VELOCITY.
C
C     METHOD.
C     -------
C
C       READ NEW WINDS WHEN NEEDED.
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
C ----------------------------------------------------------------------
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10OLD,THWOLD
      REAL,DIMENSION(NINF:NSUP) :: U10NEW,THWNEW

C ----------------------------------------------------------------------
C
      LOGICAL NEWREAD, NEWFILE
C
      CHARACTER*10 CDATEWH, CDATE
C
C ----------------------------------------------------------------------
C
C*    1. INITIALISATION.
C        ---------------
C
 1000 CONTINUE
C
C*    1.1 INTERNAL TIME COUNTER AND WIND FILE UNIT.
C         -----------------------------------------
C
 1100 CONTINUE
      IDNIN = IU17    
C
C*    1.2 WINDS ARE TAKEN FROM INTERMEDIATE STORAGE.
C         ------------------------------------------
C
 1200 CONTINUE
      IF (CDATE.LT.CDATEWH) THEN
         DO 1201 IJ=IJS,IJL
            U10NEW(IJ) = U10OLD(IJ,IG)
            THWNEW(IJ) = THWOLD(IJ,IG)
 1201    CONTINUE 
      ENDIF
C
C
C ---------------------------------------------------------------------- 
C
C
C*    2. NEW WIND INPUT.
C        ---------------
C
 2000 CONTINUE
C
      IF (CDATE.GE.CDATEWH) THEN
C
C*    2.1 WINDS HAVE TO BE READ FROM FILE.
C         --------------------------------
C
         IF (CDATE.GE.CDATEFL) THEN
C
C*    2.1.1 A NEW WIND FILE HAS TO BE USED.
C           -------------------------------
C
            IF (IG.EQ.1) THEN
C
C*    2.1.1.1 A NEW WIND FILE HAS TO BE ACCESSED.
C             -----------------------------------
C
ccccc               CALL CREWFN (IU18, CDAWIFL)
               NEWFILE = .TRUE.
               IF (ITEST.GE.2) THEN
                  WRITE(IU06,*) '      SUB. NEWWIND: NEW WIND FILE',
     1                          ' AT CDATE = ', CDATE
                  WRITE(IU06,*) '     WINDFILE DATE IS CDAWIFL ',
     1                          CDAWIFL, ' UNIT IS IU18 = ',IU18
               ENDIF
            ENDIF
            IDNIN = IU18
         ENDIF
C
C*    2.2 NEW WINDS ARE READ IN.
C         ----------------------
C
         IF (ITEST.GE.2) THEN
            IF (ITESTB.GE.IG) THEN
               WRITE(IU06,*) '       SUB. NEWWIND: READ WIND ',
     1        'IDNIN = ', IDNIN, ' AT CDATE = ', CDATE
            ENDIF
         ENDIF
ccccc         READ(IDNIN) CDA, IIG, (U10NEW(L),L=IJS,IJL),
ccccc     1                         (THWNEW(L),L=IJS,IJL)
cccc
cccc read from file instead
c

         CDA=CDTNEXT
         IIG=IG
         DO L=IJS,IJL
            U10NEW(L)=USNEXT(L,IIG)
            THWNEW(L)=DSNEXT(L,IIG)
         ENDDO


         CALL INCDATE(CDATEWH,IDELWO)
         NEWREAD = .TRUE.   
      ENDIF   
C
      RETURN
      END  
