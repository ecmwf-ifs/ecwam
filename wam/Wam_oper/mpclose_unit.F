      SUBROUTINE MPCLOSE_UNIT

C ----------------------------------------------------------------------
C
C**** *MPCLOSE_UNIT* -
C
C     J. BIDLOT     ECMWF    MAY 1996    MESSAGE PASSING
C
C
C*    PURPOSE.
C     --------
C     CLOSE ALL OPENED UNIT
C
C**   INTERFACE.
C     ----------
C     *CALL MPCLOSE_UNIT*
C
C     METHOD.
C     -------
C     LOOP ON ALL UNIT NUMBER AND INQUIRE IF UNIT IS OPENED AND THEN CLOSE IT
C
C     EXTERNALS.
C     ----------
C       NONE.
C     REFERENCE.
C     ----------
C       NONE.
C
C ----------------------------------------------------------------------
C
      LOGICAL LOPENED
C
      DO IU=1,99
        INQUIRE(UNIT=IU,OPENED=LOPENED)
        IF(LOPENED) THEN 
          CALL FLUSH(IU)
          CLOSE(IU)
        ENDIF
      ENDDO 
      END
