SUBROUTINE GET_PRESET_WGRIB_TEMPLATE(CT,IGRIB_HANDLE)

!----------------------------------------------------------------------

!**** *GET_PRESET_WGRIB_TEMPLATE* GETS THE WAVE MODEL GRIB TEMPLATE
!     WHICH HAS BEEN CREATED BY THE WAVE MODEL.

!     J. HAWKES    ECMWF NVOEMBER 20017

!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!     SUBROUTINE GET_PRESET_WGRIB_TEMPLATE(CT, IGRIB_HANDLE)
!                INPUT:
!                CT           : "I" for INTEGRATED PARAMETERS AND
!                               "S" for SPECTRA
!                OUTPUT:
!                IGRIB_HANDLE : GRIB HANDLE THAT WILL BE CREATED. 

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!       NONE.

!-------------------------------------------------------------------

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
      USE YOWGRIB_HANDLES, ONLY : NGRIB_HANDLE_WAM_S, NGRIB_HANDLE_WAM_I

      IMPLICIT NONE

      CHARACTER(LEN=1) :: CT 
      INTEGER :: IGRIB_HANDLE
      REAL :: ZHOOK_HANDLE

!-------------------------------------------------------------------

#ifdef ECMWF 
      IF (LHOOK) CALL DR_HOOK('GET_PRESET_WGRIB_TEMPLATE',0,ZHOOK_HANDLE)
#endif

    IF(CT.EQ."S") THEN
        IGRIB_HANDLE = NGRIB_HANDLE_WAM_S
    ELSE IF (CT.EQ."I") THEN
        IGRIB_HANDLE = NGRIB_HANDLE_WAM_I
    ELSE
      ![majh] Need to throw error here
      WRITE(*,*) "[majh] NO SUCH DOMAIN"
    END IF

#ifdef ECMWF 
      IF (LHOOK) CALL DR_HOOK('GET_PRESET_WGRIB_TEMPLATE',1,ZHOOK_HANDLE)
#endif

      RETURN
END SUBROUTINE GET_PRESET_WGRIB_TEMPLATE 
