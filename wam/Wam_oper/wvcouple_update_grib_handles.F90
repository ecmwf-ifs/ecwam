SUBROUTINE WVCOUPLE_UPDATE_GRIB_HANDLES(LWINIT, NSTPW, TSTEP, NSTEP, LPPSTEPS, NGRIB_HANDLE_FOR_WAM)

!----------------------------------------------------------------------

!**** *WVCOUPLE_UPDATE_GRIB_HANDLES* 

!     J. HAWKES    ECMWF MARCH 2018

!*    PURPOSE.
!     --------
!     	Initialize and update the GRIB handle (NGRIB_HANDLE_FOR_WAM) for WAM when coupled to IFS.
!		If LWINIT == false, a new grib handle will be created.
!		If LWINIT == true, the existing grib handle will be updated.


!**   INTERFACE.
!     ----------

!     SUBROUTINE WVCOUPLE_UPDATE_GRIB_HANDLES
!                INPUT: YDGEOMETRY
!                OUTPUT:

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------
!     WVCOUPLE_INIT_EARLY
!     WVXF2GB

!     REFERENCE.
!     ----------
!     EXTRACTED FROM WVXF2FB
!

!-------------------------------------------------------------------

	USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
	USE PARKIND1 , ONLY : JPIM, JPRB
    USE YOMDIM   , ONLY : TDIM
	USE IOSTREAM_MIX, ONLY : NGRIB_HANDLE_GG

	USE GRIB_API_INTERFACE, ONLY : IGRIB_CLONE
	USE GRIB_API_INTERFACE, ONLY : IGRIB_SET_VALUE
	USE GRIB_UTILS_MOD,ONLY : GRIB_SET_PARAMETER

	IMPLICIT NONE

	
	LOGICAL,	 		INTENT(IN) :: LWINIT	! Is wave model already initialized?
	INTEGER(KIND=JPIM), INTENT(IN) :: NSTPW 	! Frequency of call to the wave model (from YREWCOU)
	REAL(KIND=JPRB), 	INTENT(IN) :: TSTEP 	! IFS timestep (from YDRIP)
	INTEGER(KIND=JPIM), INTENT(IN) :: NSTEP 	! Current IFS step (from YOMCT3)
	LOGICAL, 			INTENT(IN) :: LPPSTEPS	! (from YOMDYNCORE)
	INTEGER(KIND=JPIM), INTENT(INOUT) :: NGRIB_HANDLE_FOR_WAM ! The grib handle to be created/updated

	REAL :: ZHOOK_HANDLE

	INTEGER(KIND=JPIM) :: ILEV, IGRIBCD, IDUMMY(1)
	CHARACTER (LEN=2) :: CLREPR
	CHARACTER (LEN=3) :: CLTYPE
	INTEGER(KIND=JPIM) :: ITOP_DUM,IBOT_DUM
	LOGICAL :: LLGRAD_DUM
	REAL(KIND=JPRB) :: ZDUMMY(1)
	INTEGER(KIND=JPIM) :: ISTEP_OUT

!-------------------------------------------------------------------

#ifdef ECMWF 
      IF (LHOOK) CALL DR_HOOK('WVCOUPLE_UPDATE_GRIB_HANDLES',0,ZHOOK_HANDLE)
#endif

	CALL GSTATS(1714,0)
	IF (.NOT.LWINIT) THEN
		! INITIALISE GRIB HANDLE FOR WAM
		IF(NGRIB_HANDLE_GG < 0 ) THEN
			CALL ABOR1("WVXF2GB: PROBLEM NGRIB_HANDLE_GG < 0 ")
		ELSE
			NGRIB_HANDLE_FOR_WAM=-99
			CALL IGRIB_CLONE(NGRIB_HANDLE_GG,NGRIB_HANDLE_FOR_WAM)
		ENDIF
	ENDIF

	! UPDATE GRIB HANDLE FOR WAM
	! ONLY THE FIRST TIME AND EVERY HOUR (3 or 6) for long runs AFTER THAT 
	! The wave model does not output more often than every hour
	! (it will be longer it TSEP is not a multiple of one hour)
	! We could do it all the time but then there is a potential
	! problem in *GRIB_CODE_MESSAGE* whereby grib_api is unable
	! to find an appropriate time unit when coding the step (it is a GRIB-1 limitation)  
	! There is another GRIB-1 limitation that limit the forecast to be less than 2**16=65536
	! whatever the time unit is. This limits the hourly output to (2**16-1) hours
	! Instead, one can hope to have 3 hourly output up to (2**16-1)*3 hours
	! Followed by 6 hourly output up to 2**16-1)*6 hours and so on.
	IF(MOD(3600,NINT(TSTEP)) == 0 .AND. NINT(NSTEP*(TSTEP/3600._JPRB)) <= 65535) THEN 
		ISTEP_OUT=MAX(NINT(3600._JPRB/TSTEP),1)
	ELSEIF(MOD(10800,NINT(TSTEP)) == 0 .AND. NINT(NSTEP*(TSTEP/10800._JPRB)) <= 65535 ) THEN 
		ISTEP_OUT=MAX(NINT(10800._JPRB/TSTEP),1)
	ELSEIF(MOD(21600,NINT(TSTEP)) == 0 .AND. NINT(NSTEP*(TSTEP/21600._JPRB)) <= 65535 ) THEN 
		ISTEP_OUT=MAX(NINT(21600._JPRB/TSTEP),1)
	ELSEIF(MOD(43200,NINT(TSTEP)) == 0 .AND. NINT(NSTEP*(TSTEP/43200._JPRB)) <= 65535 ) THEN 
		ISTEP_OUT=MAX(NINT(43200._JPRB/TSTEP),1)
	ELSEIF(MOD(86400,NINT(TSTEP)) == 0 .AND. NINT(NSTEP*(TSTEP/86400._JPRB)) <= 65535 ) THEN 
		ISTEP_OUT=MAX(NINT(86400._JPRB/TSTEP),1)
	ELSE
		ISTEP_OUT=1
	ENDIF

	IF ((NSTEP == NSTPW) .OR. MOD(NSTEP,ISTEP_OUT) == 0 .OR. LPPSTEPS ) THEN
		IGRIBCD=165
		ILEV=0
		CLREPR='GG'
		CLTYPE='SFC'
		CALL GRIB_CODE_MESSAGE(TSTEP,NGRIB_HANDLE_FOR_WAM,IGRIBCD,ILEV,CLREPR,CLTYPE,&
			& LLGRAD_DUM,ITOP_DUM,IBOT_DUM)
		CALL GRIB_SET_PARAMETER(NGRIB_HANDLE_FOR_WAM,IGRIBCD,ILEV,&
			& IDUMMY,IDUMMY,IDUMMY,IDUMMY,ZDUMMY)
		ZDUMMY=0.0_JPRB
		CALL IGRIB_SET_VALUE(NGRIB_HANDLE_FOR_WAM,'values',ZDUMMY)
	ENDIF
	CALL GSTATS(1714,1)

#ifdef ECMWF 
    IF (LHOOK) CALL DR_HOOK('WVCOUPLE_UPDATE_GRIB_HANDLES',1,ZHOOK_HANDLE)
#endif

      RETURN
END SUBROUTINE WVCOUPLE_UPDATE_GRIB_HANDLES 
