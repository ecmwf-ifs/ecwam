      SUBROUTINE WSTREAM_STRG(ISTREAM,CSTREAM) 

! ----------------------------------------------------------------------

!**** *WSTREAM_STRG* -

!     J. BIDLOT     ECMWF    AUGUST 2000.


!*    PURPOSE.
!     --------
!     TO RETURN A 4 CHARACTER STRING WHICH CORRESPONDS TO THE VALUE
!     ISTREAM WHEN DEALING WITH GRIB DATA AND THE FDB.
!     IF NO MATCH IS FOUND, CSTREAM IS SET TO '****'

!**   INTERFACE.
!     ----------
!     *CALL WSTREAM_STRG*(ISTREAM,CSTREAM)

! ----------------------------------------------------------------------
      INTEGER :: ISTREAM
      CHARACTER*4 CSTREAM

      IF(ISTREAM.EQ.1045) THEN
        CSTREAM = 'wave'
      ELSE IF(ISTREAM.EQ.1081) THEN
        CSTREAM = 'waef'
      ELSE IF(ISTREAM.EQ.1082) THEN
        CSTREAM = 'wasf'
      ELSE IF(ISTREAM.EQ.1083) THEN
        CSTREAM = 'mawv'
      ELSE IF(ISTREAM.EQ.1027) THEN
        CSTREAM = 'scwv'
      ELSE IF(ISTREAM.EQ.1029) THEN
        CSTREAM = 'dcwv'
      ELSE IF(ISTREAM.EQ.1095) THEN
        CSTREAM = 'wamf'
      ELSE
        CSTREAM='****'
      ENDIF

      RETURN
      END SUBROUTINE WSTREAM_STRG
