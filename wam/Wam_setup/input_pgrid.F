C
C        INPUT OF PROGRAM PGRID (PRINT GRID)
C        ===================================
C
C        USER INFORMATION IS BEING READ WITH THE PRESUMPTIONS THAT:
C         1. EVERY LINE STARTING WITH 'C' IS A COMMENT LINE
C         2. VALUES ARE PUT IN BELOW POSITIONS INDICATED WITH '-'
C            (RIGHT-JUSTIFIED, BUT CHARACTER LEFT JUSTIFIED)
C
C           **** OUTPUT TIME IN FIXED INTERVALS ****
C
C           IDATEA  = START DATE OF OUTPUT    (YYMMDDHHMM)
C           IDATEE  = END DATE OF OUTPUT      (YYMMDDHHMM)
C           IDELDO  = DATA OUTPUT TIME STEP   (SECOND)
C           U       = TIME UNIT: H = HOURS, S = SECONDS
C
C    IDATEA      IDATEE   IDELDO  U
C----------  ----------  -------  -
 7809060900  7809060900    10800  S
C
C           **** OUTPUT TIME AT SPECIFIED TIMES ****
C
C   SPECIFY A LIST WITH UP TO 100 DATES (YYMMDDHHMM) AND CLOSE THE
C   LIST BY 'END       '. IF THE LIST INCLUDES ANYTHING ELSE THAN 'END'
C   ALL OUTPUT TIME STEPS SPECFIED ABOVE ARE IGNORED.
C
CYYMMDDHHMM  YYMMDDHHMM  YYMMDDHHMM  YYMMDDHHMM  YYMMDDHHMM  YYMMDDHHMM
C----------  ----------  ----------  ----------  ----------  ----------
 END
C
C           **** FILE REQUESTS ****
C
C           IDATEF  = DATE OF FIRST FILE NAME (YYMMDDHHMM)
C           IDELFI  = FILE TIME STEP (SECONDS)
C           U       = TIME UNIT: H = HOURS, S = SECONDS
C           USD     = 3 CHARACTER USER ID
C           RID     = 3 CHARACTER RUN  ID (FOR CYBER_205 ONLY)
C           FID     = 3 CHARACTER FILE ID (NORMALLY MAP)
C           PATH    = 60 CHARACTER PATH (WITHOUT FILE NAME)
C
C    IDATEF   IDELFI  U  USD
C----------  -------  -  ---
 7809070600       24  H  NAB
C
CRID  FID  PATH
C---  ---  ---------+---------+---------+---------+---------+---------+
 CAT  MAP  /tmp/nab/coarse
C
C           **** FIELDS TO BE PRINTED ****
C
C           NAME    = 50 CHARACTER PRINT HEADER
C                        (WILL BE EXTENTED BY DATE)
C           P       = T (TRUE) IF PRINT IS REQUESTED  ELSE F (FALSE)
C
CNAME                                                P
C---------+---------+---------+---------+---------+  -
 WAM WAVE HEIGHTS                                    T
 WAM WAVE DIRECTIONS                                 T
 10M WIND SPEED                                      F
 FRICTION VELOCITY                                   F
 WIND DIRECTION                                      F
 MEAN PERIOD                                         F
 PEAK PERIOD                                         F
 DRAG COEFFICIENT                                    F
 NORMALISED WAVE STRESS                              F
C
 END
C
C ----------------------------------------------------------------------
C         END OF USER INPUT PROGRAM PGRID.
C ----------------------------------------------------------------------
