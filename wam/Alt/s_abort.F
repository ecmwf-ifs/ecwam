C=======================================================================
C
      SUBROUTINE s_abort(pvalus,klen,k)
C
C=======================================================================
      CHARACTER*128 cl_what
      DATA cl_what / '@(#) s_abort.F   2.0 96/12/04 B.Hansen'/
C
      REAL pvalus(klen*2)
C
      WRITE(*,*) ' KLEN= ',klen
      IF ( k .EQ.0) THEN
        DO 1000 j=1,klen
           WRITE(*,'(i10,F20.7)')j,pvalus(j)
 1000   CONTINUE
        RETURN
      ENDIF
C
      DO 2000 j=1,klen
         WRITE(*,'(i10,2(2x,E20.12))')j,pvalus(j), pvalus(j+klen)
 2000 CONTINUE
C
      CALL ABORT
      END
