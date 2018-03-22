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
      USE GRIB_API_INTERFACE
      USE GRIB_UTILS_MOD

!-----------------------------------------------------------------------
      IMPLICIT NONE

      REAL :: ZHOOK_HANDLE

      TYPE (IO_SERV),      INTENT (INOUT) :: YDIOS
      REAL,                INTENT (IN)    :: FIELD (YDIOS%MODELPAR%YWAM%NGX,YDIOS%MODELPAR%YWAM%NGY)
      TYPE (IOFLDDESC),    INTENT (IN)    :: FLDDESC
      INTEGER (KIND=JPIM), INTENT (INOUT) :: KGRIB_HANDLE
      
      INTEGER(KIND=JPKSIZE_T)             :: KBYTES
      INTEGER(KIND=JPIM)                  :: ISIZE, IERR
      INTEGER, ALLOCATABLE                :: KGRIB_BUFR(:)
      
      INTEGER(KIND=JPIM) :: ILEV, IGRIBCD, IDUMMY(1)
      CHARACTER (LEN=2) :: CLREPR
      CHARACTER (LEN=3) :: CLTYPE
      REAL(KIND=JPRB) :: ZDUMMY(1)

!-----------------------------------------------------------------------

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTWSPEC_IO_SERV_REC',0,ZHOOK_HANDLE)
#endif

      ASSOCIATE( &
            WAMPAR => YDIOS%MODELPAR%YWAM, &
            ECPAR  => YDIOS%MODELPAR%YECGRIB, &
            WAMHDR => FLDDESC%YWAM &
      )

      WRITE(*,*) 'Writing wave model spectral field from IO server.'

      ! Update GRIB handles
      ! CALL WVCOUPLE_UPDATE_GRIB_HANDLES( .TRUE., WAMPAR%NSTPW, ECPAR%TSTEP, FLDDESC%YWAM%NSTEP, ECPAR%LPPSTEPS, KGRIB_HANDLE)

      IGRIBCD=165
      ILEV=0
      CLREPR='GG'
      CLTYPE='SFC'
      ! CALL GRIB_CODE_MESSAGE(TSTEP,KGRIB_HANDLE,IGRIBCD,ILEV,CLREPR,CLTYPE,&
      !       & LLGRAD_DUM,ITOP_DUM,IBOT_DUM)
      CALL GRIB_SET_TIME( &
            KGRIB_HANDLE, &
            ECPAR%LPPSTEPS, &
            FLDDESC%YWAM%NSTEP,&
            ECPAR%TSTEP, &
            ECPAR%NSTEPINI, &
            ECPAR%LVAREPS, ECPAR%NLEG, &
            ECPAR%NFCHO_TRUNC_INI, ECPAR%NFCLENGTH_INI, &
            ECPAR%NREFERENCE, ECPAR%NSTREAM,&
            & 'fc',FLDDESC%IPREVPP, 165)

      CALL GRIB_SET_PARAMETER(KGRIB_HANDLE,IGRIBCD,ILEV,&
            & IDUMMY,IDUMMY,IDUMMY,IDUMMY,ZDUMMY)

      ZDUMMY=0.0_JPRB
      CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'values',ZDUMMY)

      ! Do field-specific GRIB encoding
      CALL WGRIBENCODE_IO_SERV( WAMPAR%NGX, WAMPAR%NGY, FIELD, YDIOS, FLDDESC, KGRIB_HANDLE  )

      CALL IGRIB_GET_MESSAGE_SIZE( KGRIB_HANDLE, KBYTES )

      ISIZE = ( KBYTES + WAMPAR%NPRECI - 1 ) / WAMPAR%NPRECI

      ALLOCATE( KGRIB_BUFR( ISIZE ) )

      CALL IGRIB_GET_MESSAGE( KGRIB_HANDLE, KGRIB_BUFR )

      ! Use the the FDB opened by the IO server
      CALL WGRIB2FDB ( &
            WAMPAR%IU06, WAMPAR%ITEST, &
            KGRIB_HANDLE, ISIZE, KGRIB_BUFR, &
            WAMPAR%CFDB2DSP, ECPAR%NFDBREF, &
            .TRUE., &
            WAMPAR%IMDLGRBID_G, WAMPAR%IMDLGRBID_M, &
            IERR &
      )

      IF(IERR.NE.0)THEN
          WRITE(WAMPAR%IU06,*) ' ------------------------'
          WRITE(WAMPAR%IU06,*) ' ERROR ACCESSING FDB '
          WRITE(WAMPAR%IU06,*) ' FDB ERROR CODE IS ',IERR
          WRITE(WAMPAR%IU06,*) ' ------------------------'
          WRITE(*,*) ' ------------------------'
          WRITE(*,*) ' ERROR ACCESSING FDB '
          WRITE(*,*) ' FDB ERROR CODE IS ',IERR
          WRITE(*,*) ' ------------------------'
          CALL ABORT1
      ENDIF

      DEALLOCATE( KGRIB_BUFR )

      END ASSOCIATE

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTWSPEC_IO_SERV_REC',1,ZHOOK_HANDLE)
#endif

      RETURN

      END SUBROUTINE OUTWSPEC_IO_SERV_REC
