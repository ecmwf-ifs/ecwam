      REAL FUNCTION WAM_USER_CLOCK()
C
C   Returns system clock converted to micro-seconds
C
      CALL SYSTEM_CLOCK(ICOUNT,ICOUNT_RATE,ICOUNT_MAX)
      IF(ICOUNT.NE.HUGE(0)) THEN
        WAM_USER_CLOCK = FLOAT(ICOUNT) * 1.E6 / FLOAT(ICOUNT_RATE)
      ELSE
        WAM_USER_CLOCK = 0.0
      ENDIF    
C
      RETURN
      END
