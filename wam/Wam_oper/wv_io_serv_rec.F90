SUBROUTINE WV_IO_SERV_REC (YDIOS,FIELD,FLDDESC,KGRIB_HANDLE)

!----------------------------------------------------------------------

!**** *WV_IO_SERV_REC*
!      SECOND PART OF OUTWSPEC_IO_SERV AND OUTINT_IO_SERV.
!      PERFORMED BY THE IO SERVER

!      J. HAWKES   ECMWF  OCTOBER 2017 

!-------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
      USE YOMIO_SERV, ONLY : IO_SERV
      USE IOFLDDESC_MOD,  ONLY : IOFLDDESC
      USE GRIB_API_INTERFACE
      USE GRIB_UTILS_MOD

!-----------------------------------------------------------------------
      IMPLICIT NONE
#include "wgribencode_io_serv.intfb.h"
#include "abort1.intfb.h"
#include "wgrib2fdb.intfb.h"


      TYPE (IO_SERV),      INTENT (INOUT) :: YDIOS
      REAL(KIND=JWRB),     INTENT (INOUT) :: FIELD (YDIOS%MODELPAR%YWAM%NGX,YDIOS%MODELPAR%YWAM%NGY)
      TYPE (IOFLDDESC),    INTENT (IN)    :: FLDDESC
      INTEGER(KIND=JWIM),  INTENT (INOUT) :: KGRIB_HANDLE
      
      INTEGER(KIND=JPKSIZE_T)             :: KBYTES
      INTEGER(KIND=JWIM)                  :: ISIZE, IERR
      INTEGER(KIND=JWIM)                  :: ILEV, IGRIBCD, IDUMMY(1), LPPSTEPS
      INTEGER(KIND=JWIM), ALLOCATABLE     :: KGRIB_BUFR(:)
      
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: ZDUMMY(1)

      CHARACTER(LEN=2) :: CLREPR
      CHARACTER(LEN=3) :: CLTYPE

      LOGICAL :: LLDUM
      LOGICAL :: LFDBOPEN_IO

!-----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WV_IO_SERV_REC',0,ZHOOK_HANDLE)

      ASSOCIATE( &
 &          WAMPAR => YDIOS%MODELPAR%YWAM, &
 &          ECPAR  => YDIOS%MODELPAR%YECGRIB, &
 &          WAMHDR => FLDDESC%YWAM &
 &    )

      ! This is a limited version of WVCOUPLE_UPDATE_GRIB_HANDLES/GRIB_CODE_MESSAGE
      IF(ECPAR%NLOCGRB == 1 .OR. ECPAR%NLOCGRB == 36) THEN
            CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'type',ECPAR%CTYPE)
      ENDIF

      IGRIBCD  = 165
      LPPSTEPS = -1
      CALL GRIB_SET_TIME( &
 &          KGRIB_HANDLE, &
 &          ECPAR%LPPSTEPS, &
 &          WAMHDR%NSTEP, &
 &          ECPAR%TSTEP, &
 &          ECPAR%NSTEPINI, &
 &          ECPAR%LVAREPS, ECPAR%NLEG, &
 &          ECPAR%NFCHO_TRUNC_INI, ECPAR%NFCLENGTH_INI, &
 &          ECPAR%NREFERENCE, ECPAR%NSTREAM, &
 &          'fc', -1, IGRIBCD, LLDUM )

      CALL GRIB_SET_PARAMETER(KGRIB_HANDLE,IGRIBCD,ILEV,&
            & IDUMMY,IDUMMY,IDUMMY,IDUMMY,ZDUMMY)

      ZDUMMY=0.0_JWRB
      CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'values',ZDUMMY)

      CALL WGRIBENCODE_IO_SERV( WAMPAR%NGX, WAMPAR%NGY, FIELD, YDIOS, FLDDESC, KGRIB_HANDLE )

      CALL IGRIB_GET_MESSAGE_SIZE( KGRIB_HANDLE, KBYTES )

      ISIZE = ( KBYTES + WAMPAR%NPRECI - 1 ) / WAMPAR%NPRECI

      ALLOCATE( KGRIB_BUFR( ISIZE ) )

      CALL IGRIB_GET_MESSAGE( KGRIB_HANDLE, KGRIB_BUFR )

      ! Use the the FDB opened by the IO server
      LFDBOPEN_IO=.TRUE.
      CALL WGRIB2FDB ( &
 &           WAMPAR%IU06, WAMPAR%ITEST, &
 &           KGRIB_HANDLE, ISIZE, KGRIB_BUFR, &
 &           WAMPAR%CFDB2DSP, ECPAR%NFDBREF, &
 &          LFDBOPEN_IO, &
 &          WAMPAR%IMDLGRBID_G, WAMPAR%IMDLGRBID_M, &
 &          IERR &
 &    )

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

      IF (LHOOK) CALL DR_HOOK('WV_IO_SERV_REC',1,ZHOOK_HANDLE)

END SUBROUTINE WV_IO_SERV_REC
