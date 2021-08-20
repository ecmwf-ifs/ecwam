      SUBROUTINE CLOSEND (CDATE, CDATEWH, LLNEWREAD, LLNEWFILE)

! ----------------------------------------------------------------------

!      obsolete ???
!**** *CLOSEND* - HANDLING OF WINDS OUTSIDE MULTITASKED AREA    

!     P.A.E.M. JANSSEN  KNMI/ECMWF  SEPTEMBER 1994
!     J. BIDLOT         ECMWF       FEBRUARY  1996  MESSAGE PASSING
!     S. ABDALLA        ECMWF       OCTOBER   1996  AIR DENSITY AND Zi/L

!*    PURPOSE.
!     --------

!       READ WINDS WHEN NEEDED.                                   

!**   INTERFACE.
!     ----------

!     *CALL* *CLOSEND*(CDATE,CDATEWH,LLNEWREAD,LLNEWFILE)
!      *LLNEWREAD* - TRUE IF NEW WINDS HAVE BEEN READ
!      *LLNEWFILE* - TRUE IF NEW WIND FILE HAS BEEN OPENED

!       *INCDATE*   - UPDATE DATE TIME GROUP.

!     REFERENCE.
!     ----------

!       NONE

! ----------------------------------------------------------------------

!*    *PARAMETER*  FOR ARRAY DIMENSIONS.

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWSTAT  , ONLY : IDELPRO  ,IDELWI   ,NPROMA_WAM
      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      CHARACTER(LEN=14), INTENT(IN) :: CDATEWH, CDATE
      LOGICAL, INTENT(INOUT) :: LLNEWREAD, LLNEWFILE

      INTEGER(KIND=JWIM) :: IDELWH

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('CLOSEND',0,ZHOOK_HANDLE)


        LLNEWREAD = .FALSE.

!*    2. UPDATE WIND COUNTERS IF LAST BLOCK HAS BEEN DONE.
!        -------------------------------------------------

        IF (LLNEWFILE) THEN
!*        UPDATE WIND FILE TIME COUNTER AND UNITS.
          LLNEWFILE = .FALSE.
          IDELWH = MAX(IDELWI,IDELPRO)
          CALL INCDATE(CDAWIFL,IDELWH)
          CALL INCDATE(CDATEFL,IDELWH)
        ENDIF

!*      UPDATE WIND FIELD COUNTER.
        CDATEWO=CDATEWH

      IF (LHOOK) CALL DR_HOOK('CLOSEND',1,ZHOOK_HANDLE)

      END SUBROUTINE CLOSEND
