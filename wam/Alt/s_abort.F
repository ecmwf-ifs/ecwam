C=======================================================================
C
      SUBROUTINE S_ABORT(pvalus,klen,k)
C
C=======================================================================
      CHARACTER*128 cl_what
      DATA cl_what / '@(#) s_abort.F   2.0 96/12/04 B.Hansen'/

      REAL pvalus(klen)
      IF ( k .EQ.0) THEN
        DO 1000 j=1,klen
           WRITE(*,'(i10,F20.7)')j,pvalus(j)
 1000   CONTINUE
        RETURN
      ENDIF
      DO 2000 j=1,klen
         WRITE(0,'(i10,F20.7)')j,pvalus(j)
 2000 CONTINUE
      CALL ABORT
      END
