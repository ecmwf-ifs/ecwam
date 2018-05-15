      SUBROUTINE WAM_NPROMA(NTOT, MTHREADS, NPROMA)
! ----------------------------------------------------------------------

!**** *WAM_NPROMA* - ADJUSTS DOWN NPROMA BASED ON NTOT AND MTHREADS
!                    TO INSURE THAT NTOT IS DISTRIBUTED AS EQUALLY
!                    AS POSSIBLE ON MTHREADS THREADS

!**   INTERFACE.
!     ----------

!       *WAM_NPROMA(NTOT, MTHREADS, NPROMA)
!       *NTOT*    - TOTAL NUMBER OF GRID POINTS TO BE DISTRIBUTED
!       *MTHREADS*- NUMBER OF OPENMP THREADS
!       *NPROMA*  - MAXIMUM NUMBER OF POINTS PER THREAD PER CALL

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: NTOT, MTHREADS
      INTEGER(KIND=JWIM), INTENT(INOUT) :: NPROMA
      INTEGER(KIND=JWIM) :: NLOOP, NREMAIN, NPROMA_MAX

! ----------------------------------------------------------------------

      NPROMA=MAX(NPROMA,1)
      NPROMA_MAX=NPROMA
      NLOOP=NTOT/(MTHREADS*NPROMA)
      NREMAIN=NTOT-NLOOP*MTHREADS*NPROMA
      IF(NREMAIN.GT.0 .AND. MTHREADS.GT.1 .AND. NTOT.GT.0) THEN
        DO WHILE (NREMAIN.LT.(NPROMA*MTHREADS) .AND. NPROMA.GT.1)
          NPROMA=NPROMA-1
          NREMAIN=NTOT-NLOOP*MTHREADS*NPROMA
        ENDDO
        NPROMA=MIN(NPROMA+1,NPROMA_MAX)
      ENDIF

      END SUBROUTINE WAM_NPROMA
