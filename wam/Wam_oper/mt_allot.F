      SUBROUTINE MT_ALLOT (F,ETOI, USOI, USMO, THMO,IJS,IJL,NSIZE,NWORK)
C
C ----------------------------------------------------------------------
C
C**** *MT_ALLOT* - SET UP TASKS TO RUN UPDATE.
C
C     B.       HANSEN  MPIfM/E.C.M.W.F.    OCTOBER   1994
C     (  FROM ROUTINE PHYS BY
C        P.A.E.M. JANSSEN K.N.M.I./E.C.M.W.F. SEPTEMBER 1994 )
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
C       *CALL* *MT_ALLOT (F,ETOI, USOI, USMO, THMO, IJS,IJL,IG,NSIZE,
C                         NWORK)*
C          *F*      -  REAL      SPECTRUM ( FIRST GUESS IN INPUT ;
C                                ANALYSIS IN OUTPUT ).
C          *ETOI*   -  REAL      WAVE ENERGY (FROM O.I.).
C          *USOI*   -  REAL      USTAR ( FROM O.I.).
C          *USMO*   -  REAL      USTAR ( FIRST GUESS IN INPUT ;
C                                ANALYSIS IN OUTPUT ).
C          *THMO*   -  REAL      WIND DIRECTION.
C          *IJS*    -  INDEX OF FIRST GRIDPOINT
C          *IJL*    -  INDEX OF LAST GRIDPOINT
C          *IG*     -  BLOCK NUMBER
C          *NSIZE*  -  SIZE OF SUB TASK
C          *NWORK*  -  TOTAL AMOUNT OF SUB TASKS
C
C     METHOD.
C     -------
C
C       NONE.
C
C     EXTERNALS.
C     ---------
C
C       *UPDATE*    - ANALYSES THE WAM MODEL SPECTRA.
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
#include "commpp.h"
C
C ----------------------------------------------------------------------
C
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: F
      REAL,DIMENSION(NINF:NSUP) :: ETOI,USOI,USMO,THMO
C
C ----------------------------------------------------------------------
C
C
C*    1. INITIALISATION.
C        ---------------
C              
C  
 1000 CONTINUE  
      NN = INCINT(NCHUNK,1)
      IF (NN.LE.NWORK) THEN
         MIJS = IJS + (NN-1)*NSIZE
         MIJL = MIJS + NSIZE-1
         MIJS = MIN(MIJS,IJL)
         MIJL = MIN(MIJL,IJL)
C                
         CALL UPDATE (F, ETOI, USOI, USMO, THMO, MIJS, MIJL)
         GO TO 1000
      END IF
C
C ----------------------------------------------------------------------
C
C
      RETURN
      END     
