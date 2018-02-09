      SUBROUTINE WAMSENDNORM(WNORM,ISEND,IRECV)

! ----------------------------------------------------------------------

!**** *WAMSENDNORM* -


!*    PURPOSE.
!     --------

!     SENDS WNORM FROM ISEND TO IRECV

!**   INTERFACE.
!     ----------
!     *CALL*  WAMSENDNORM(WNORM,ISEND,IRECV)
!          *WNORM*     REAL    WAM NORMS 
!          *ISEND*     INTEGER PE WHERE IT IS DEFINED
!          *IRECV*     INTEGER PE WHERE IT SHOULD BE WRITTEN OUT 

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------


!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMPP   , ONLY : IRANK
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE MPL_MODULE,ONLY : MPL_SEND, MPL_RECV

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: ISEND, IRECV
      REAL(KIND=JWRB), INTENT(INOUT) :: WNORM(4)

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('WAMSENDNORM',0,ZHOOK_HANDLE)

      IF(ISEND.NE.IRECV) THEN
        CALL GSTATS(675,0)
        IF(IRANK.EQ.ISEND) THEN
        ELSEIF(IRANK.EQ.IRECV) THEN
        ENDIF
        CALL GSTATS(675,1)
      ENDIF

      IF (LHOOK) CALL DR_HOOK('WAMSENDNORM',1,ZHOOK_HANDLE)

      END SUBROUTINE WAMSENDNORM
