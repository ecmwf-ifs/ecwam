SUBROUTINE WGRIBENCODE_IO_SERV (I1, I2, FIELD, &
                             YDIOS, FLDDESC, &
                             IGRIB_HANDLE)

! ----------------------------------------------------------------------

!****  *WGRIBENCODE_IO_SERV*  ENCODES WAM MODEL FIELD INTO GRIB CODE AND OUTPUT

!       J. BIDLOT    ECMWF JULY 2009: USE GRIB API 

!      PURPOSE.
!      --------
!      
!      A wrapper around wgribencode which runs on the IO server
!      Uses variables baked into the IO server at initializations
!      (see suiosctmpl.F90)
!      and variables passed on a per-field basis in the field descriptor.

!      INTERFACE.
!      ----------

!      METHOD.
!      -------

!      EXTERNALS.
!      ----------

!      REFERENCES.
!      -----------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOMIO_SERV, ONLY : IO_SERV
      USE IOFLDDESC_MOD,  ONLY : IOFLDDESC

      USE GRIB_API_INTERFACE
      USE YOMHOOK  , ONLY : LHOOK, DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "wgribencode.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN)     :: I1, I2
      INTEGER(KIND=JWIM), INTENT(INOUT)  :: IGRIB_HANDLE
      REAL(KIND=JWRB), INTENT(INOUT)     :: FIELD(I1,I2)
      REAL(KIND=JWRB)                    :: ZHOOK_HANDLE
      TYPE (IO_SERV), INTENT (IN)   :: YDIOS
      TYPE (IOFLDDESC), INTENT (IN) :: FLDDESC

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('WGRIBENCODE_IO_SERV',0,ZHOOK_HANDLE)

      CALL GSTATS(1709,0)

      ! Static variables, stored in IO server ( set in IO_SERV_SUIOSCTMPL_WAM )
      ASSOCIATE( WAMPAR => YDIOS%MODELPAR%YWAM )
      ! Field-specific data ( set in OUTWSPEC_IO_SERV )
      ASSOCIATE ( WAMFLD => FLDDESC%YWAM )

      CALL WGRIBENCODE( WAMPAR%IU06, WAMPAR%ITEST, &
 &                      I1, I2, &
 &                      FIELD, &
 &                      WAMFLD%ITABLE, WAMFLD%IPARAM, &
 &                      WAMFLD%KLEV, &
 &                      WAMFLD%IANGLE, WAMFLD%IFREQ, &
 &                      WAMFLD%CDATE, WAMFLD%IFCST, WAMFLD%MARSTYPE, &
 &                      WAMPAR%PPMISS, WAMPAR%PPEPS, WAMPAR%PPREC, WAMPAR%PPRESOL, WAMPAR%PPMIN_RESET, WAMPAR%NTENCODE, &
 &                      .TRUE., &
 &                      WAMFLD%DATE_TIME_WINDOW_END, &
 &                      WAMPAR%NGRBRESS, WAMPAR%LNEWLVTP, WAMPAR%LPADPOLES, &
 &                      SIZE(WAMPAR%NLONRGG), WAMPAR%NLONRGG, WAMPAR%IRGG, &
 &                      WAMPAR%AMONOP, WAMPAR%AMOSOP, WAMPAR%XDELLA, WAMPAR%CLDOMAIN, &
 &                      WAMFLD%KCOUSTEP, WAMFLD%LRSTST0, &
 &                      WAMPAR%ZMISS, &
 &                      IGRIB_HANDLE)
      
      END ASSOCIATE
      END ASSOCIATE
!     

      CALL GSTATS(1709,1)

      IF (LHOOK) CALL DR_HOOK('WGRIBENCODE_IO_SERV',1,ZHOOK_HANDLE)

END SUBROUTINE WGRIBENCODE_IO_SERV
