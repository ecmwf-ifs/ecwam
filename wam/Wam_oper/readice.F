      SUBROUTINE READICE (IU06,IU90,ITEST,NSTART,NEND,IREAD,ITAG)
C
C--------------------------------------------------------------------
C
C**** *READICE* *DEFINES SEA POINTS AS ICE IF THE CORRESPONDING
C                SST IS LT 271.5 KELVIN.
C
C     B. HANSEN    MPIFM/ECMWF     NOVEMBER 92
C     J. BIDLOT    ECMWF           JUNE 1996  MESSAGE PASSING
C     P. JANSSEN   ECMWF           OCTOBER 1997 SIMPLIFY BY INTRODUCING
C                                               ICEMASK
C     J. BIDLOT    ECMWF           JANURAY 1998 : CORRECT THRESHOLD FOR ICE
C                                                 TO BE LT 271.50 
C
C     PURPOSE.
C     --------
C            DYNAMICALLY ADJUST NORTH AND SOUTH ICE BOUNDARIES
C            OF THE WAMODEL.
C
C**   INTERFACE
C     ---------
C     *CALL* *READICE (IU06,IU90,ITEST,NSTART,NEND,IREAD,ITAG)*
C
C     *IU06*      OUTPUT UNIT FOR STANDARD OUTPUT.
C     *IU90*      UNIT 90   SEA SUFACE TEMPERATURES IN KELVIN (INPUT)
C     *ITEST*     OUTPUT FLAG IF LE 2 NO EXTRA OUTPUT IS GENERATED
C     *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
C     *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
C     *IREAD*     RANK OF THE PROCESS WHICH INPUTS THE SST FIELD
C
C
C     METHOD.
C     -------
C            READ GLOBAL TEMPRERATURES (UNIT 99).
C            DEFINE BLOCK(IG) AND INDEX(IJ) FOR EACH POINT
C            WITH TEMPERATURE < 271.50 K  -->  ICE.
C            EXTRACT WATER POINTS.
C            SORT ICE POINTS ACCORDING TO BLOCK AND INDEX IN A BLOCK.
C            KEEP INFORMATION IN COMMON COMICE FOR LATER USE.
C
C
C    EXTERNAL.
C    ---------
C             *MICEP*     DEFINE BLOCK AND INDEX IN A BLOCK FOR
C                         EACH ICE POINT (T < 271.5 KELVIN).
C             *READT*     READ GLOBAL TEMPERATURES.
C             *MPBCASTGRDFLD* BROADCAST SST FIELD TO ALL OTHER PE's
C
C---------------------------------------------------------------------
C
#include "comtemp.h"
C
#include "commesspass.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
C------------------------------------------------------------------------
C     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS
C
      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
C
C------------------------------------------------------------------------
C
      REAL,ALLOCATABLE :: TEMPG(:,:), WORK(:)
C
C*    1. READ SST FROM UNIT 90.
C        ----------------------
C
 1000 CONTINUE
C
      ALLOCATE(TEMPG(NXT,NYT),WORK(6))
      TEMPG=273.15

      IF (LMESSPASS) THEN
C
C        READ SST FIELD ON PE IREAD
C
         IF (IRANK.EQ.IREAD) CALL READT(IU06,IU90,ITEST,TEMPG)
C
C        BROADCAST  SST FIELD ARRAY INFORMATION TO OTHER PE'S
C
         IF (IRANK.EQ.IREAD) THEN
            WORK(1)=TMOWEP
            WORK(2)=TMOSOP
            WORK(3)=TMOEAP
            WORK(4)=TMONOP
            WORK(5)=DELLA
            WORK(6)=DELLO
         ENDIF

         CALL MPBCASTGRDFLD(IREAD,ITAG,6,1,WORK)
         ITAG=ITAG+1

         IF (IRANK.NE.IREAD) THEN
            TMOWEP=WORK(1)
            TMOSOP=WORK(2)
            TMOEAP=WORK(3)
            TMONOP=WORK(4)
            DELLA=WORK(5)
            DELLO=WORK(6)
         ENDIF
C
C        BROADCAST  SST FIELD TO OTHER PE'S
C
         CALL MPBCASTGRDFLD(IREAD,ITAG,NXT,NYT,TEMPG)
         ITAG=ITAG+1

      ELSE
         CALL READT(IU06,IU90,ITEST,TEMPG)
      ENDIF
C
C-------------------------------------------------------------------
C
C*    2. DEFINE ICE POINTS AND COMPUTE ICEFLAG AND WRITE COMMON BLOCK
C        COMICE.
C        ----------------------------------------------
C
 2000 CONTINUE
C
      CALL MICEP (IU06, NSTART, NEND, TEMPG)
C

      DEALLOCATE(TEMPG,WORK)
C
      RETURN
      END
