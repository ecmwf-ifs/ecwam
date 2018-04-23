C
C       USER INPUT OF BOUINT.            CREATED 19/4/91
C       =====================
C
C        USER INFORMATION IS BEING READ WITH THE PRESUMPTIONS THAT:
C         1. EVERY LINE STARTING WITH 'C' IS A COMMENT LINE
C         2. LINE 1 IS A COMMENT LINE
C         3. VALUES ARE PUT IN BELOW POSITIONS INDICATED WITH '-'
C            (RIGHT-JUSTIFIED, BUT CHARACTER LEFT JUSTIFIED)
C
C ----------------------------------------------------------------------
C
C           **** DATES AND TIME STEPS FOR OUTPUT  ****
C
C     IDATEA   INTERGER*10   START DATE OF OUTPUT    (YYMMDDHHMM)
C     IDATEE   INTERGER*10   END DATE OF OUTPUT      (YYMMDDHHMM)
C     IDELPRF  INTERGER*7    DATA OUTPUT TIME STEP   (SECOND)
C     U        CHARACTER*1   TIME UNIT  S SECONDS  OR  H HOURS
C
C    IDATEA      IDATEE  IDELPRF U
C----------  ----------  ------- -
 7809060610  7809080600      400 S
C
C           **** FILE NAMES AND FILE DATES ****
C
C     IDATES   INTERGER*10   DATE OF FIRST FILE NAME (YYMMDDHHMM)
C     IDELFI   INTERGER*7    FILE TIME STEP (SECONDS)
C     U        CHARACTER*1   TIME UNIT  S SECONDS  OR  H HOURS
C     USD      CHARACTER*3   USER ID
C     RDI      CHARACTER*3   RUN  ID (INPUT DATA)
C     FDI      CHARACTER*3   FILE ID (INPUT DATA)
C     PATHI    CHARACTER*60  PATH (INPUT DATA)
C     RDO      CHARACTER*3   RUN  ID (OUTPUT DATA)
C     FDO      CHARACTER*3   FILE ID (OUTPUT DATA)
C     PATHO    CHARACTER*60  PATH (OUTPUT DATA)
C
C    IDATES   IDELFI U  USD
C----------  ------- -  ---
 7809070600       24 H  NAB
C
CRDI  FDI  PATHI
C---  ---  ---------+---------+---------+---------+---------+---------+
 CAT  CBO  /tmp/nab/coarse
C
CRDO  FDO  PATHO
C---  ---  ---------+---------+---------+---------+---------+---------+
 FAT  FBI  /tmp/nab/fine
C
C ----------------------------------------------------------------------
C                END OF USER INPUT PROG. BOUINT.
C ----------------------------------------------------------------------
