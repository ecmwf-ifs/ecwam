      SUBROUTINE PHYS (FL3,FL1,FL,IJS,IJL,IG,NSIZE,NWORK,CDATE,CDTPRA,
     1                 U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     2                 U10NEW,THWNEW,USNEW,Z0NEW,SL,FCONST)
C
C ----------------------------------------------------------------------
C
C**** *PHYS* - SET UP PROCESS TO RUN IMPLSCH AND PROPAGS
C
C     P.A.E.M. JANSSEN K.N.M.I./E.C.M.W.F. SEPTEMBER 1994
C     J. BIDLOT ECMWF FEBRUARY 1996-97 MESSAGE PASSING
C
C*    PURPOSE.
C     --------
C
C       BY SPLITTING UP INNER DO LOOP IN A NUMBER OF PARTS WE DISTRIBUTE
C       THE WORK OVER A NUMBER OF PROCESSORS.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *PHYS (FL3,FL,IJS,IJL,IG,NSIZE,NWORK,CDATE,CDTPRA,
C    1                 U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
C    2                 U10NEW,THWNEW,USNEW,Z0NEW,SL,FCONST)*
C      *FL3*    - FREQUENCY SPECTRUM(INPUT AND OUTPUT).
C      *FL1*    - FREQUENCY SPECTRUM(INPUT).
C      *IJS*    - INDEX OF FIRST GRIDPOINT
C      *IJL*    - INDEX OF LAST GRIDPOINT
C      *IG*     - BLOCK NUMBER
C      *NSIZE*  - SIZE OF SUB TASK
C      *NWORK*  - TOTAL AMOUNT OF SUB TASKS
C      *CDATE*  - DATE OF RUN
C      *CDTPRA* - PROPAGATION DATE
C      *U10NEW*    NEW WIND SPEED IN M/S.
C      *U10OLD*    INTERMEDIATE STORAGE OF MODULUS OF WIND
C                  VELOCITY.
C      *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
C                  NOTATION (POINTING ANGLE OF WIND VECTOR,
C                  CLOCKWISE FROM NORTH).
C      *THWOLD*    INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
C                  WIND VELOCITY.
C      *USNEW*     NEW FRICTION VELOCITY IN M/S.
C      *USOLD*     INTERMEDIATE STORAGE OF MODULUS OF FRICTION
C                  VELOCITY.
C      *Z0NEW*     ROUGHNESS LENGTH IN M.
C      *Z0OLD*     INTERMEDIATE STORAGE OF ROUGHNESS LENGTH IN
C                  M.
C      *TAUW*      WAVE STRESS IN (M/S)**2
C      *SL*        REAL      TOTAL SOURCE FUNCTION ARRAY.
C      *FCONST*    REAL      = 1 FOR PROGNOSTIC FREQUENCY BANDS.
C                            = 0 FOR DIAGNOSTIC FREQUENCY BANDS.
C
C     METHOD.
C     -------
C
C       NONE.
C
C     EXTERNALS.
C     ---------
C
C       *IMPLSCH*    - INTEGRATES SOURCE TERMS.     
C
C     REFERENCE.
C     ----------
C
C       NONE.                                                     
C

C ----------------------------------------------------------------------
C
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.
C
#include "param.h"
C
#include "commltsk.h"
C
#include "commesspass.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
C ----------------------------------------------------------------------
C
      CHARACTER*10 CDATE, CDTPRA
C     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL,FL3,FL1,SL
      REAL,DIMENSION(NINF:NSUP,NFRE) :: FCONST
      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10OLD,THWOLD,USOLD,Z0OLD,TAUW
      REAL,DIMENSION(NINF:NSUP) :: U10NEW,THWNEW,USNEW,Z0NEW

C ----------------------------------------------------------------------
C
C
C*    1. INITIALISATION.
C        ---------------
C              
C  
! DETERMINE WHAT FORM OF PARALELLISM IS USED BASED ON THE TYPE
! OF MACHINE ARCHTECTURE

       IF(LMESSPASS) THEN
C
         IF (CDATE.EQ.CDTPRA) THEN
            CALL PROPAGS(FL1, FL3, IJS, IJL, IG, SL, FCONST)
         ENDIF
         CALL IMPLSCH (FL3, FL, IJS, IJL, IG,
     1                 THWOLD,USOLD,TAUW,Z0OLD,
     2                 U10NEW,THWNEW,USNEW,Z0NEW,SL,FCONST)
C
       ELSE

 1000 CONTINUE  
C                
         NN = INCINT(NCHUNK,1)
         IF (NN.LE.NWORK) THEN
            MIJS = IJS + (NN-1)*NSIZE
            MIJL = MIJS + NSIZE-1
            MIJS = MIN(MIJS,IJL)
            MIJL = MIN(MIJL,IJL)
C

            IF (CDATE.EQ.CDTPRA) THEN
               CALL PROPAGS(FL1, FL3, MIJS, MIJL, IG, SL, FCONST)
            ENDIF
            CALL IMPLSCH (FL3, FL, MIJS, MIJL, IG,
     1                    THWOLD,USOLD,TAUW,Z0OLD,
     2                    U10NEW,THWNEW,USNEW,Z0NEW,SL,FCONST)
          GO TO 1000
          END IF

      ENDIF
C
C ----------------------------------------------------------------------
C
C
      RETURN
      END     
