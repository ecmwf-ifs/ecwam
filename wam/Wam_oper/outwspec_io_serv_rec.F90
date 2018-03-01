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

      REAL                    :: ZMISS

!-----------------------------------------------------------------------

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTWSPEC_IO_SERV_REC',0,ZHOOK_HANDLE)
#endif

      WRITE(*,*) "OUTWSPEC_IO_SERV_REC"    
      WRITE(*,*) "IANG, IFREQ RECEIVED", FLDDESC%IFREQ, FLDDESC%IANGLE


      ! What do I need to do encoding?
      ! ZMISS (which should be applied earlier)
      ! CDATE, IFCST, MARSTYPE

      ! (got) NGX, NGY, FIELD, ITABLE, IPARAM, IK, IM 

      ! Fixed
      ! IPARAM=251
      ! ITABLE=140

      ! UPDATE GRIB HANDLES

      !CALL WGRIBENCODE_IO_SERV( NGX, NGY, FIELD, ITABLE, IPARAM, 0, IK, IM, CDATE, IFCST, MARSTYPE, KGRIB_HANDLE  )


      !            CALL WGRIBENCODE(IU06, ITEST, NGX, NGY, FIELD,
      ! &                     ITABLE, IPARAM, 0, IK , IM, 
      ! &                     CDATE, IFCST, MARSTYPE,
      ! &                     IGRIB_HANDLE)

      ! Copy to ISENDMSG

      !               CALL WGRIBOUT(IU06, ITEST,
      !        &                    LFDBIOOUT, CFDB2DSP, NWFDBREF, LFDBOPEN,
      !        &                    IUOUT,
     !         &                    IGRIB_HANDLE,MSGSIZE(IRANK),ISENDMSG(1))

      ! Deallocate ISENDMSG

      ! GRIBD RELEASE HANDLE

      ! Or this?
      !CALL IGRIB_NEW_FROM_MESSAGE(IGRBHNDL,IRECVMSG(1:MSIZE,IRCV))
      !CALL WGRIBOUT(IU06, ITEST,
      !&                    LFDBIOOUT, CFDB2DSP, NWFDBREF, LFDBOPEN,
      !&                    IUOUT,
      !&                    IGRBHNDL,MSIZE,IRECVMSG(1,IRCV))

      ! CALL IGRIB_CLOSE_FILE(IUOUT)

      

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTWSPEC_IO_SERV_REC',1,ZHOOK_HANDLE)
#endif

      RETURN

      END SUBROUTINE OUTWSPEC_IO_SERV_REC
