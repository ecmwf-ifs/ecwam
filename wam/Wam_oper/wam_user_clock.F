      REAL FUNCTION WAM_USER_CLOCK()

!   Returns system clock converted to micro-seconds

      CALL SYSTEM_CLOCK(ICOUNT,ICOUNT_RATE,ICOUNT_MAX)
      IF(ICOUNT.NE.HUGE(0)) THEN
        WAM_USER_CLOCK = REAL(ICOUNT) * 1.E6 / REAL(ICOUNT_RATE)
      ELSE
        WAM_USER_CLOCK = 0.0
      ENDIF    

      RETURN
      END




