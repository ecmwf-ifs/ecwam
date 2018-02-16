      SUBROUTINE OUTWSPEC_IO_SERV_REC (YDIOS,FIELD,FLDDESC,KGRIB_HANDLE)

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
      REAL,                INTENT (IN)    :: FIELD (:,:)
      TYPE (IOFLDDESC),    INTENT (IN)    :: FLDDESC
      INTEGER (KIND=JPIM), INTENT (INOUT) :: KGRIB_HANDLE

!-----------------------------------------------------------------------

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTWSPEC_IO_SERV_REC',0,ZHOOK_HANDLE)
#endif

      WRITE(*,*) "OUTWSPEC_IO_SERV_REC"    
      WRITE(*,*) "IANG, IFREQ RECEIVED", FLDDESC%IFREQ, FLDDESC%IANGLE
      

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTWSPEC_IO_SERV_REC',1,ZHOOK_HANDLE)
#endif

      RETURN

      END SUBROUTINE OUTWSPEC_IO_SERV_REC
