      SUBROUTINE WVDEALLOC
C
C ----------------------------------------------------------------------
C
C**** *WVDEALLOC* - WAVE MODEL DEALLOCATION 
C
C     J. BIDLOT     ECMWF   JANUARY 1997 ATMOSPHERIC COUPLING
C
C
C ----------------------------------------------------------------------
C
      USE WAVE_MODULE
C
C ----------------------------------------------------------------------
C     1.  DEALLOCATE NECESSARY ARRAYS
C         -------------------------
      IF(ALLOCATED(U10NEW)) DEALLOCATE(U10NEW)
      IF(ALLOCATED(THWNEW)) DEALLOCATE(THWNEW)
      IF(ALLOCATED(USNEW)) DEALLOCATE(USNEW)
      IF(ALLOCATED(Z0NEW)) DEALLOCATE(Z0NEW)

      IF(ALLOCATED(FL2)) DEALLOCATE(FL2)
      IF(ALLOCATED(FL3)) DEALLOCATE(FL3)

      IF(ALLOCATED(SL)) DEALLOCATE(SL)
      IF(ALLOCATED(FCONST)) DEALLOCATE(FCONST)

C
      RETURN
      END
