      SUBROUTINE OUTWSPEC_IO_SERV_REC (YDIOS,YDFLDSC,KGRIB_HANDLE,KSTEP, PFLD)

!----------------------------------------------------------------------

!**** *OUTWSPEC*  SECOND PART OF OUTWSPEC_IO_SERV, PERFORMED ON THE IO SERVER

!     J. HAWKES   ECMWF  OCTOBER 2017 

!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!     SUBROUTINE OUTWSPEC_IO_SERV_REC ()

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!       NONE.

!-------------------------------------------------------------------

      USE PARKIND1,       ONLY : JPRB, JPIM
      USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
      USE YOMIO_SERV, ONLY : IO_SERV
      USE IOFLDDESC_MOD,  ONLY : IOFLDDESC

!-----------------------------------------------------------------------
      IMPLICIT NONE

      REAL :: ZHOOK_HANDLE

      TYPE (IO_SERV),      INTENT (INOUT) :: YDIOS
      TYPE (IOFLDDESC),    INTENT (IN)    :: YDFLDSC
      INTEGER (KIND=JPIM), INTENT (INOUT) :: KGRIB_HANDLE
      INTEGER (KIND=JPIM), INTENT (IN)    :: KSTEP
      REAL (KIND=JPRB),    INTENT (IN)    :: PFLD (:)

!-----------------------------------------------------------------------

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTWSPEC_IO_SERV_REC',0,ZHOOK_HANDLE)
#endif

      WRITE(*,*) "OUTWSPEC_IO_SERV_REC"      

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTWSPEC_IO_SERV_REC',1,ZHOOK_HANDLE)
#endif

      RETURN

      END SUBROUTINE OUTWSPEC_IO_SERV_REC
