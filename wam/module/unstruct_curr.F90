MODULE UNSTRUCT_CURR
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
      REAL(KIND=JWRU), ALLOCATABLE :: CURTXY (:,:)
      REAL(KIND=JWRU), ALLOCATABLE :: CURTXY1(:,:)
      REAL(KIND=JWRU), ALLOCATABLE :: CURTXY2(:,:)
      INTEGER(KIND=JWIM) :: recTime1, recTime2
      TYPE(TIMEPERIOD) RecTimeCurr
      character(len=*), parameter :: eFileBnd = 'wwm_curr_format.nc'
      REAL(KIND=JWRU) :: WAV_CurrTime = 0.
      LOGICAL :: UseSingleFile = .FALSE.
      LOGICAL :: IsInitializedSingleFile = .FALSE.
#ifdef NETCDF_OUTPUT_WAM

      CONTAINS
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_UNSTRUCT_SINGLEFILE_CURRENT
      USE NETCDF
      USE WAV_NETCDF_FCT, ONLY : WAV_GENERIC_NETCDF_ERROR
      USE WAV_NETCDF_FCT, ONLY : CF_EXTRACT_TIME
      IMPLICIT NONE
      character(len=*), parameter :: CallFct = "INIT_FILE_BOUNDARY"
      REAL(KIND=JWRU), allocatable :: ListTime_mjd(:)
      character (len=100) :: eStrUnitTime
      REAL(KIND=JWRU) ConvertToDay, eTimeStart, eTimeBnd
      INTEGER(KIND=JWIM) :: varid, ncid
      INTEGER(KIND=JWIM), dimension(nf90_max_var_dims) :: dimids
      INTEGER(KIND=JWIM) :: istat
      INTEGER(KIND=JWIM) :: nbtime_mjd
      INTEGER(KIND=JWIM) :: iTime
      REAL(KIND=JWRU) eWD1, eWD2, eDiff, eDiff1, eDiff2
      REAL(KIND=JWRU) DeltaDiff, eDir
      logical IsAssigned
      INTEGER(KIND=JWIM) :: ID, ID1, ID2
      !
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'eFileBnd=', TRIM(eFileBnd)
#endif
      !
      CALL TEST_FILE_EXIST_DIE(eFileBnd, "Needed for currents")
      ISTAT = NF90_OPEN(TRIM(eFileBnd), NF90_NOWRITE, ncid)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)
      !
      ! reading the time
      !
      CALL READ_LIST_TIME(ncid, RecTimeCurr)
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
      ALLOCATE(CURTXY1(2,MNP), CURTXY2(2,MNP), stat=istat)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SINGLE_READ_CURRENT(eFile, CURT, IT)
      USE WAV_NETCDF_FCT, ONLY : WAV_GENERIC_NETCDF_ERROR
      USE NETCDF
      IMPLICIT NONE
      character(len=*), intent(in) :: eFile
      REAL(KIND=JWRU), intent(out) :: CURT(2,MNP)
      INTEGER(KIND=JWIM), intent(in) :: IT
      REAL(KIND=JWRU) CURT_GL(2,np_global)
      REAL(KIND=JWRU) singleRead(np_global)
      INTEGER(KIND=JWIM) :: istat, ncid, var_id
      INTEGER(KIND=JWIM) :: IP, IS, ID, idx
      INTEGER(KIND=JWIM) :: IPglob
      character(len=*), parameter :: CallFct = "SINGLE_READ_CURRENT"
      !
      ! We have this inversion of the order because that is so in WWM at
      ! the present time
      !
      ISTAT = NF90_OPEN(TRIM(eFile), NF90_NOWRITE, ncid)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)
      !
      ISTAT = nf90_inq_varid(ncid, 'CURTX', var_id)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)
      !
      ISTAT = NF90_GET_VAR(ncid, var_id, singleRead, start=(/1,IT/), count = (/np_global,1/))
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)
      CURT_GL(1,:) = singleRead
      !
      ISTAT = nf90_inq_varid(ncid, 'CURTY', var_id)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)
      !
      ISTAT = NF90_GET_VAR(ncid, var_id, singleRead, start=(/1,IT/), count = (/np_global,1/))
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)
      CURT_GL(2,:) = singleRead
      !
      ISTAT = NF90_CLOSE(ncid)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)
      !
      ! Now reassigning
      !
      DO IP=1,MNP
        IPglob=iplg(IP)
        CURT(:,IP) = CURT_GL(:,IPglob)
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_CURTXY_SINGLEFILE
      USE YOWSTAT  , ONLY : IDELT
      USE WAV_NETCDF_FCT, ONLY : WAV_GET_ETIMEDAY
      IMPLICIT NONE
      REAL(KIND=JWRU) eTimeDay
      REAL(KIND=JWRU) w1, w2
      INTEGER(KIND=JWIM) :: iTime1, iTime2
      UseSingleFile = .TRUE.
      IF (IsInitializedSingleFile .eqv. .FALSE.) THEN
        IsInitializedSingleFile=.FALSE.
        CALL INIT_UNSTRUCT_SINGLEFILE_CURRENT
      END IF
      CALL WAV_GET_ETIMEDAY(eTimeDay, WAV_CurrTime)
      CALL FIND_MATCH_TIME(RecTimeCurr, eTimeDay, iTime1, w1, iTime2, w2)
      IF (iTime1 .ne. recTime1) THEN
        CALL SINGLE_READ_CURRENT(eFileBnd, CURTXY1, iTime1)
        recTime1 = iTime1
      END IF
      IF (iTime2 .ne. recTime2) THEN
        CALL SINGLE_READ_CURRENT(eFileBnd, CURTXY2, iTime2)
        recTime2 = iTime2
      END IF
      CURTXY = w1*CURTXY1 + w2*CURTXY2
      WAV_CurrTime = WAV_CurrTime + DBLE(IDELT)
      LCALC=.TRUE.
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_CURTXY
      USE YOWCURR, ONLY : U, V
      IMPLICIT NONE
      INTEGER(KIND=JWIM) :: IP, IG
      IF (UseSingleFile .eqv. .FALSE.) THEN
        IG=1
        DO IP=1,MNP
          CURTXY(1,IP)=REAL(U(IP,IG),JWRU)
          CURTXY(2,IP)=REAL(V(IP,IG),JWRU)
        END DO
        LCALC=.TRUE.
      END IF
      END SUBROUTINE SET_CURTXY
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
END MODULE
