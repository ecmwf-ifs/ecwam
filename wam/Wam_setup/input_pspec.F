C
C        INPUT OF PROGRAM PSPEC (PRINT SPECTRA)
C        ======================================
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
C           FID     = 3 CHARACTER FILE ID (NORMALLY SPE, CBO, OR FBI)
C           PATH    = 60 CHARACTER PATH (WITHOUT FILE NAME)
C
C    IDATEF   IDELFI  U  USD
C----------  -------  -  ---
 7809070600       24  H  NAB
C
CRID  FID  PATH
C---  ---  ---------+---------+---------+---------+---------+---------+
 CAT  SPE  /tmp/nab/coarse
C
C           **** OUTPUT SIDES ****
C
C           XLONG   = LONGITUDE IN DEGREE
C           XLAT    = LATITUDE IN DEGREE
C           NAME    = SIDE NAME (UP TO 20 CHARACTER)
C
C           ----------------------------------------------
C           THE PROGRAM DIMENSIONS ALLOW UP TO 26 SIDES
C           SPECTRA ARE NOT SORTED AS TIME SERIES.
C           ----------------------------------------------
C
C   XLONG      XLAT  NAME
C--------  --------  ---------+---------+
   10.000    -9.500  FETCH  0.5 DEG
   10.000    -9.000  FETCH  1.0 DEG
   10.000    -8.500  FETCH  1.5 DEG
   10.000    -8.000  FETCH  2.0 DEG
   10.000    -7.500  FETCH  2.5 DEG
   10.000    -7.000  FETCH  3.0 DEG
   10.000    -5.500  FETCH  4.5 DEG
   10.000    -5.000  FETCH  5.0 DEG
   10.000    -4.500  FETCH  5.5 DEG
   10.000     0.000  FETCH  10. DEG
   10.000     1.000  FETCH  11. DEG
   10.000     2.000  FETCH  12. DEG
   10.000     3.000  FETCH  13. DEG
   10.000     4.000  FETCH  14. DEG
   10.000     5.000  FETCH  15. DEG
   10.000     6.000  FETCH  16. DEG
   10.000    10.000  FETCH  20. DEG
C
C           THE INPUT LIST HAS TO BE CLOSED BY AN END LINE
C
    0.000     0.000  END
C
C ----------------------------------------------------------------------
C         END OF USER INPUT PROGRAM PSPEC.
C ----------------------------------------------------------------------
