      SUBROUTINE PREWIND (NC,NR,U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     1                    IUGB, IVGB, ILEN)

C ----------------------------------------------------------------------
C
C**** *PREWIND* - PREPARES WIND DATA FOR WAVE MODEL.
C
C     P.GROENWOUD     DELFT HYDRAULICS LABORATORY  OKTOBER 1986
C
C     E. BAUER        MPI       FEB 1987   VERSION FOR CDC 205 HAMBURG.
C
C     S. HASSELMANN   MPI       MAY 1987   COMBINED CDC 205 AND CRAY
C                                          VERSIONS.
C     W. BRUEGGEMANN  MPI    AUGUST 1988   SIMPLIFIED PROGRAM.
C
C     L. ZAMBRESKY    ECMWF    JUNE 1988   MODIFIED EXTENSIVELY FOR
C                                          COUPLING TO SPECTRAL MODEL.
C
C     H. GUNTHER      ECMWF    JUNE 1990   MODIFIED FOR CYCLE_4.

!     J. BIDLOT       ECMWF    FEBRUARY 1996 MESSAGE PASSING
C
C*    PURPOSE.
C     --------
C
C       EVALUATE WIND SPEED AND DIRECTION AT WAVE MODEL GRID POINTS.
C
C**   INTERFACE.
C     ----------
C
C     *CALL* *PREWIND (NC,NR,U10OLD,THWOLD,USOLD,TAUW,Z0OLD)
C
C      *NC*    FIRST  DIMENSION OF INPUT WIND FIELD.
C      *NR*    SECOND DIMENSION OF INPUT WIND FIELD.
C      *U10OLD*    REAL      INTERMEDIATE STORAGE OF MODULUS OF WIND
C                            VELOCITY.
C                            CLOCKWISE FROM NORTH).
C      *THWOLD*    REAL      INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
C                            WIND VELOCITY.
C      *USOLD*     REAL      INTERMEDIATE STORAGE OF MODULUS OF FRICTION
C                            VELOCITY.
C      *Z0OLD*     REAL      INTERMEDIATE STORAGE OF ROUGHNESS LENGTH IN
C                            M.
C      *TAUW*      REAL      WAVE STRESS IN (M/S)**2
C
C
C       *UNIT* *DESCRIPTION*
C
C          IU01    INPUT WIND DATA (SUB READWND).
C          IU06    PRINTER OUTPUT (SUB INITMDL).
C          IUVELO  OUTPUT OF BLOCKED WIND FIELDS. (SUB CREWFN).
C          IUSCR   SCRATCH UNITS FOR ALL BLOCKS (INTERMEDIATE STORAGE,
C                  INPUT/OUTPUT) (SUB INITMDL).
C
C     METHOD.
C     -------
C
C       INPUT WIND FIELDS WHICH CAN BE COMPONENTS OF
C                USTAR, U10, USTRESS
C       ARE TRANSFORMED TO FRICTION VELOCITIES.
C       THE INPUT FIELDS HAVE TO BE ON A LAT /LONG GRID.
C       SEE SUB READWND FOR FORMATS AND HEADER INFORMATION,
C       WHICH HAVE TO BE GIVEN TO THE PROGRAM.
C
C       A DOUBLE LINEAR INTERPOLATION IN SPACE IS PERFORMED
C       ONTO THE MODEL BLOCKS.
C       IF THE WIND OUTPUT TIMSTEP IS LESS THAN THE INPUT TIMESTEP
C       A LINEAR INTERPOLATION IN TIME IS PERFORMED.
C
C       THERE ARE TWO POSSIBILITIES WITH RESPECT TO THE WIND
C       OUTPUT FILES:
C
C           1. PROPAGATION TIMESTEP >= WIND INPUT STEP
C              ONE OUTPUT FILE CONTAINS IDELPRO/IDELWO WINDFIELDS
C              I.E. INFORMATION FOR ONE PROPAGATION TIMESTEP.
C              TIME FILE(I+1)= TIME FILE(I)+ IDELPRO
C
C           2. PROPAGATION TIMESTEP < INPUT WIND TIMESTEP
C              ONE OUTPUT FILE CONTAINS IDELWI/IDELWO WINDFIELDS
C              I.E. INFORMATION FOR ONE WIND INPUT TIMESTEP.
C              TIME FILE(I+1)= TIME FILE(I) + IDELWI
C
C     EXTERNALS.
C     ----------
C
C       *ABORT*     - TERMINATES PROCESSING.
C       *AIRSEA*    - SURFACE LAYER STRESS.
C       *CREWFN*    - CREATES A WIND FILE NAME.
C       *GETWND*    - PROCESSES ONE WIND FIELD.
C       *INCDAT*    - INCREMENTS DATE TIME GROUP.
C       *LOCINT*    - INTERPOLATES IN SPACE.
C       *NOTIM*     - STEERING SUB FOR INTERPOLATION IN SPACE ONLY.
C       *READWND*   - READS A WIND FIELD.
C       *TIMIN*     - STEERING SUB FOR INTERPOLATION IN SPACE AND TIME.
C       *WAMWND*    - BLOCKS A WIND FIELD AND CONVERTS TO USTAR.
C
C     REFERENCE.
C     -----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "comgrid.h"
C
#include "txtgrid.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
#include "comstat.h"
C
#include "txtstat.h"
C
#include "comtest.h"
C
#include "txttest.h"
C
#include "comunit.h"
C
#include "txtunit.h"
C
#include "comwind.h"
C
#include "txtwind.h"
C
#include "comwndg.h"
C
#include "txtwndg.h"
C
#include "comtext.h"
C
!
      CHARACTER*10 CDTWIE, CDTWIS, ZERO

      LOGICAL LLFRST
      DATA LLFRST / .TRUE. /
C
C ----------------------------------------------------------------------
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10OLD,THWOLD,USOLD,Z0OLD,TAUW
C
C ----------------------------------------------------------------------
      INTEGER IUGB(ILEN), IVGB(ILEN)
C ----------------------------------------------------------------------
C
C
C*    1. BEGIN AND END DATES OF WIND FIELDS TO BE PROCESSED.
C        ---------------------------------------------------
C
      ZERO = ' '
C      
      IF (CDA.EQ.ZERO) THEN
C
C        IF START FROM PRESET FIELDS DO FIRST FIELD IN ADDITION.
C
         WRITE(6,*) "  PREWIND AT CDATEA = ", CDATEA, 
     .              " AND CDA = ", CDA
         CALL FLUSH (6)
         CDTWIS = CDATEA
      ELSE
         CDTWIS = CDAWIFL
         IDELWH = -MAX(IDELPRO,IDELWI)+IDELWO
         CALL INCDATE (CDTWIS,IDELWH)
      ENDIF
C      
      IF(CDAWIFL.LT.CDATEE) THEN
        CDTWIE = CDAWIFL
         LLFRST = .FALSE.
      ELSE
        CDTWIE = CDATEE
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    2. PROCESS WIND FIELDS.
C        --------------------
C
 2000 CONTINUE
      IF (IDELWO.GE.IDELWI) THEN
C
C*    2.1 NO TIME INTERPOLATION.
C         ----------------------
C
         IF (ITEST.GE.2) THEN
            WRITE (IU06,*) '   SUB. PREWIND: WIND REQUEST'
            WRITE (IU06,*) '     NO TIME INTERPOLATION'
            WRITE (IU06,*) '     START OF PERIOD IS    CDTWIS = ',CDTWIS
            WRITE (IU06,*) '     END   OF PERIOD IS    CDTWIE = ',CDTWIE
            WRITE (IU06,*) '     WIND INPUT TIME STEP  IDELWI = ',IDELWI
            WRITE (IU06,*) '     WIND OUTPUT TIME STEP IDELWO = ',IDELWO
         ENDIF
         CALL NOTIM (CDTWIS, CDTWIE, NC, NR,
     1               U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     2               IUGB, IVGB, ILEN)
      ELSE
C
C*    2.2 TIME INTERPOLATION.
C         -------------------
C
         IF (ITEST.GE.2) THEN
            WRITE (IU06,*) '   SUB. PREWIND: WIND REQUEST'
            WRITE (IU06,*) '     TIME INTERPOLATION'
            WRITE (IU06,*) '     START OF PERIOD IS    CDTWIS = ',CDTWIS
            WRITE (IU06,*) '     END   OF PERIOD IS    CDTWIE = ',CDTWIE
            WRITE (IU06,*) '     WIND INPUT TIME STEP  IDELWI = ',IDELWI
            WRITE (IU06,*) '     WIND OUTPUT TIME STEP IDELWO = ',IDELWO
         ENDIF
         CALL TIMIN (CDTWIS, CDTWIE, NC, NR,
     1               U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     2               IUGB, IVGB, ILEN)
      ENDIF

      RETURN
      END
