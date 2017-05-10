MODULE UNSTRUCT_WIND
  !
  ! Rule of conversion of angle in WAM.
  ! theta_{WAM} = 90 - theta_{trigonometric}
  !
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWPARAM , ONLY : NANG, NFRE
      USE YOWSTAT,  ONLY : IREFRA
      USE YOWUNPOOL, ONLY : LCALC
      USE yowpd, only: MNE=>ne, INE, MNP=>npa, NP_RES => np
      USE yowpd, only: XP=>x, YP=>y, DEP=>z
      USE yowpd, only: exchange, np_global
      USE yownodepool, only : iplg
      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal
      USE WAV_NETCDF_FCT
      IMPLICIT NONE
      REAL(KIND=JWRU), ALLOCATABLE :: WIND (:,:)
      REAL(KIND=JWRU), ALLOCATABLE :: WIND1(:,:)
      REAL(KIND=JWRU), ALLOCATABLE :: WIND2(:,:)
      integer recTime1, recTime2
      integer, allocatable :: Indexes_boundary(:)
      TYPE(TIMEPERIOD) RecTimeWind
      character(len=*), parameter :: eFileWind = 'wwm_wind_format.nc'
      real*8 :: WAV_WindTime = 0
      LOGICAL :: WIND_UseSingleFile = .FALSE.
      LOGICAL :: WIND_IsInitializedSingleFile = .FALSE.
      CONTAINS
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_UNSTRUCT_SINGLEFILE_WIND
      USE NETCDF
      USE WAV_NETCDF_FCT, ONLY : WAV_GENERIC_NETCDF_ERROR
      USE WAV_NETCDF_FCT, ONLY : CF_EXTRACT_TIME
      IMPLICIT NONE
      character(len=*), parameter :: CallFct = "INIT_FILE_BOUNDARY"
      real*8, allocatable :: ListTime_mjd(:)
      character (len=100) :: eStrUnitTime
      real(KIND=JWRU) ConvertToDay, eTimeStart, eTimeBnd
      INTEGER varid, ncid
      integer, dimension(nf90_max_var_dims) :: dimids
      integer istat
      integer nbtime_mjd
      integer iTime
      real(KIND=JWRU) eWD1, eWD2, eDiff, eDiff1, eDiff2
      real(KIND=JWRU) DeltaDiff, eDir
      logical IsAssigned
      integer ID, ID1, ID2
      !
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'eFileWind=', TRIM(eFileWind)
#endif
      !
      CALL TEST_FILE_EXIST_DIE(eFileWind, "Needed for currents")
      ISTAT = NF90_OPEN(TRIM(eFileWind), NF90_NOWRITE, ncid)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)
      !
      ! reading the time
      !
      CALL READ_LIST_TIME(ncid, RecTimeWind)
      !
      ! closing the file
      !
      ISTAT = NF90_CLOSE(ncid)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 15, ISTAT)
!
! allocate wave boundary arrays ... 
!
      recTime1=-1
      recTime2=-1
      ALLOCATE(WIND(2,MNP), WIND1(2,MNP), WIND2(2,MNP), stat=istat)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SINGLE_READ_WIND(eFile, WIND, IT)
      USE WAV_NETCDF_FCT, ONLY : WAV_GENERIC_NETCDF_ERROR
      USE NETCDF
      IMPLICIT NONE
      character(len=*), intent(in) :: eFile
      REAL(KIND=JWRU), intent(out) :: WIND(2,MNP)
      integer, intent(in) :: IT
      REAL(KIND=JWRU) WIND_GL(2,np_global)
      REAL(KIND=JWRU) singleRead(np_global)
      integer istat, ncid, var_id
      integer IP, IS, ID, idx
      integer IPglob
      character(len=*), parameter :: CallFct = "SINGLE_READ_WIND"
      !
      ! We have this inversion of the order because that is so in WWM at
      ! the present time
      !
      ISTAT = NF90_OPEN(TRIM(eFile), NF90_NOWRITE, ncid)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)
      !
      ISTAT = nf90_inq_varid(ncid, 'Uwind', var_id)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)
      !
      ISTAT = NF90_GET_VAR(ncid, var_id, singleRead, start=(/1,IT/), count = (/np_global,1/))
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)
      WIND_GL(1,:) = singleRead
      !
      ISTAT = nf90_inq_varid(ncid, 'Vwind', var_id)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)
      !
      ISTAT = NF90_GET_VAR(ncid, var_id, singleRead, start=(/1,IT/), count = (/np_global,1/))
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)
      WIND_GL(2,:) = singleRead
      !
      ISTAT = NF90_CLOSE(ncid)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)
      !
      ! Now reassigning
      !
      DO IP=1,MNP
        IPglob=iplg(IP)
        WIND(:,IP) = WIND_GL(:,IPglob)
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_WIND_SINGLEFILE
      USE YOWSTAT  , ONLY : IDELT
      USE WAV_NETCDF_FCT, ONLY : WAV_GET_ETIMEDAY
      IMPLICIT NONE
      REAL*8 eTimeDay
      CHARACTER(LEN=15) STIME
      REAL(KIND=JWRU) w1, w2
      integer iTime1, iTime2
      WIND_UseSingleFile = .TRUE.
      IF (WIND_IsInitializedSingleFile .eqv. .FALSE.) THEN
        WIND_IsInitializedSingleFile=.FALSE.
        CALL INIT_UNSTRUCT_SINGLEFILE_WIND
      END IF
      CALL WAV_GET_ETIMEDAY(eTimeDay, WAV_WindTime)
      CALL FIND_MATCH_TIME(RecTimeWind, eTimeDay, iTime1, w1, iTime2, w2)
      IF (iTime1 .ne. recTime1) THEN
        CALL SINGLE_READ_WIND(eFileWind, WIND1, iTime1)
        recTime1 = iTime1
      END IF
      IF (iTime2 .ne. recTime2) THEN
        CALL SINGLE_READ_WIND(eFileWind, WIND2, iTime2)
        recTime2 = iTime2
      END IF
      WIND = w1*WIND1 + w2*WIND2

      WAV_WindTime = WAV_WindTime + DBLE(IDELT)
      CALL MJD2CT(WAV_WindTime, STIME)
      WRITE(740+MyRankGlobal,*) 'WIND(min/max)=', minval(WIND), maxval(WIND)
      WRITE(740+MyRankGlobal,*) 'SET_WIND_SINGLEFILE, STIME=', TRIM(STIME)
      FLUSH(740+MyRankGlobal)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_VARIABLES_UNSTRUCTURED(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,  &
     &                    ROAIRO, ZIDLOLD, CICOVER, CITHICK)
      USE YOWWIND  , ONLY : WSPMIN
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NBLO
      USE YOWGRID  , ONLY : IJS      ,IJL, IGL
      USE YOWWIND  , ONLY : NSTORE   , CDA, CDTNEXT, CDATEFL,           &
     &            FF_NEXT  ,RWFAC
      IMPLICIT NONE
      REAL,intent(inout), DIMENSION(NINF:NSUP,NBLO) :: USOLD, Z0OLD, TAUW
      REAL,intent(out), DIMENSION(NINF:NSUP,NBLO) :: U10OLD, THWOLD
      REAL,intent(out), DIMENSION(NINF:NSUP,NBLO) :: ROAIRO, ZIDLOLD
      REAL,intent(out), DIMENSION(NINF:NSUP,NBLO) :: CICOVER, CITHICK
      !
      integer IP, IJ
      integer ISTORE, IG
      real eAIRD, eDS, UU, VV, eWindSpeedRaw
      real eWindDir, eWindSpeed
      real UU_wind, VV_wind
      real eCICOVER, eCITHICK, eZIDL, eUS
      logical UseCurrent
      logical HaveUV
      !
      NSTORE=1
      IF(.NOT.ALLOCATED(CDTNEXT)) ALLOCATE(CDTNEXT(NSTORE))
      IF(.NOT.ALLOCATED(FF_NEXT)) ALLOCATE(FF_NEXT(IJS(1):IJL(1),NSTORE))
      !
      CALL SET_WIND_SINGLEFILE
      ISTORE = 1
      IG = 1
      DO IP=1,MNP
        IJ = IP + NINF -1
        UU_wind = WIND(1,IP)
        VV_wind = WIND(2,IP)
        eAIRD=1
        eZIDL=0
        eCICOVER = 0
        eCITHICK = 0
        !
        UseCurrent = .FALSE.
        HaveUV = .FALSE.
        IF (UseCurrent .and. HaveUV) THEN
!          UU = UU_wind - RWFAC*U(IJ,IG)
!          VV = VV_wind - RWFAC*V(IJ,IG)
        ELSE
          UU = UU_wind
          VV = VV_wind
        END IF
        !
        eWindSpeedRaw = SQRT(UU**2 + VV**2)
        IF (eWindSpeedRaw .NE. 0.) THEN
          eWindDir = ATAN2(UU,VV)
        ELSE
          eWindDir = 0.
        ENDIF
        eWindSpeed=MAX(eWindSpeedRaw, WSPMIN)
        !
        U10OLD(IJ,IG)=eWindSpeed
        THWOLD(IJ,IG)=eWindDir
        ROAIRO(IJ,IG)=eAIRD
        ZIDLOLD(IJ,IG)=eZIDL
        CICOVER(IJ,IG)=eCICOVER
        CITHICK(IJ,IG)=eCITHICK
        !
        FF_NEXT(IJ,ISTORE)%WSWAVE=eWindSpeed
        FF_NEXT(IJ,ISTORE)%WDWAVE=eWindDir
        FF_NEXT(IJ,ISTORE)%AIRD=eAIRD
        FF_NEXT(IJ,ISTORE)%ZIDL=eZIDL
        FF_NEXT(IJ,ISTORE)%CIFR=eCICOVER
        FF_NEXT(IJ,ISTORE)%CITH=eCITHICK
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_WIND_UNSTRUCTURED
      USE YOWWIND  , ONLY : NXFF     ,NYFF     ,FIELDG
      IMPLICIT NONE
      INTEGER :: I, J, IP
      CALL SET_WIND_SINGLEFILE
      DO J=1,NYFF
        DO I=1,NXFF
          IP = I
          FIELDG(I,J)%UWND = WIND(1,IP)
          FIELDG(I,J)%VWND = WIND(2,IP)
        ENDDO
      ENDDO
      END SUBROUTINE
END MODULE
