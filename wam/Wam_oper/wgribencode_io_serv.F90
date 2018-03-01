      SUBROUTINE WGRIBENCODE_IO_SERV (I1, I2, FIELD, &
                             IK, IM, &
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

      USE YOMIO_SERV, ONLY : IO_SERV
      USE IOFLDDESC_MOD,  ONLY : IOFLDDESC

      USE GRIB_API_INTERFACE
      USE YOMHOOK  , ONLY : LHOOK, DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN)           :: I1, I2
      INTEGER, INTENT(IN)           :: ITABLE, IPARAM, KLEV, IK, IM, IFCST
      INTEGER, INTENT(INOUT)        :: IGRIB_HANDLE
      CHARACTER(LEN=2), INTENT(IN)  :: MARSTYPE
      CHARACTER(LEN=14), INTENT(IN) :: CDATE
      INTEGER, INTENT(IN)           :: DATE_TIME_WINDOW_END
      INTEGER, INTENT(IN)           :: KCOUSTEP
      LOGICAL, INTENT(IN)           :: LRSTST0
      REAL, INTENT(INOUT)           :: FIELD(I1,I2)
      REAL                          :: ZHOOK_HANDLE
      TYPE (IO_SERV), INTENT (IN)   :: YDIOS
      TYPE (IOFLDDESC), INTENT (IN) :: FLDDESC

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('WGRIBENCODE_IO_SERV',0,ZHOOK_HANDLE)
#endif

      CALL GSTATS(1709,0)

      ASSOCIATE( IOS => YDIOS%MODELPAR )

      CALL WGRIBENCODE( IOS%WV_IU06, IOS%WV_ITEST, &
                    IOS%WV_NGX, IOS%WV_AMONOP_NGY, &
                    FIELD, &
                    FLDDESC%ITABLE, FLDDESC%IPARAM, &
                    FLDDESC%KLEV, &
                    FLDDESC%WV_ANG, FLDDESC%WV_FREQ, &
                    FLDDESC%WV_CDATE, FLDDESC%IFCST, FLDDESC%WV_MARSTYPE, &
                    IOS%WV_PPMISS, IOS%WV_PPEPS, IOS%WV_PPREC, IOS%WV_PPRESOL, IOS%WV_PPMIN_RESET, IOS%WV_NTENCODE, &
                    FLDDESC%DATE_TIME_WINDOW_END, &
                    IOS%WV_NGRBRESS, IOS%WV_LNEWLVTP, IOS%WV_LPADPOLES, &
                    IOS%WV_NLONRGG, IOS%WV_IRGG, &
                    IOS%WV_AMONOP, IOS%WV_AMOSOP, IOS%WV_XDELLA, IOS%WV_CLDOMAIN, &
                    FLDDESC%KCOUSTEP, FLDDESC%LRSTST0, &
                    IOS%WV_ZMISS, &
                    IGRIB_HANDLE)
      
      END ASSOCIATE
!     

      CALL GSTATS(1709,1)

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('WGRIBENCODE_IO_SERV',1,ZHOOK_HANDLE)
#endif

      RETURN
      END SUBROUTINE WGRIBENCODE_IO_SERV
