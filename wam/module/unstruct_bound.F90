MODULE UNSTRUCT_BOUND
  !
  ! Rule of conversion of angle in WAM.
  ! theta_{WAM} = 90 - theta_{trigonometric}
  !
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWPARAM , ONLY : NANG, NFRE
      USE YOWSTAT,  ONLY : IREFRA
      USE yowpd, only: MNE=>ne, INE, MNP=>npa, NP_RES => np
      USE yowpd, only: XP=>x, YP=>y, DEP=>z
      USE yowpd, only: exchange
      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal
      USE WAV_NETCDF_FCT
      IMPLICIT NONE
      REAL(KIND=JWRU), ALLOCATABLE :: SPSIG(:)
      REAL(KIND=JWRU), ALLOCATABLE :: WBAC (:,:,:)
      REAL(KIND=JWRU), ALLOCATABLE :: WBAC1(:,:,:)
      REAL(KIND=JWRU), ALLOCATABLE :: WBAC2(:,:,:)
      integer nbDirWWM, nbFreqWWM
      REAL(KIND=JWRU), ALLOCATABLE :: SPDIR_WWM(:)
      REAL(KIND=JWRU), ALLOCATABLE :: SPSIG_WWM(:)
      REAL(KIND=JWRU), ALLOCATABLE :: SPDIR_WAM(:)
      integer, allocatable :: WWM_ID1(:), WWM_ID2(:)
      REAL(KIND=JWRU), allocatable :: WWM_WD1(:), WWM_WD2(:)
      integer recTime1, recTime2
      integer, allocatable :: Indexes_boundary(:)
      TYPE(TIMEPERIOD) RecTimeBnd
      real*8 :: WAV_BoucTime = 0
      character(len=*), parameter :: eFileBnd = 'wwm_bouc_format.nc'
      INTEGER                :: IWBMNP ! number of wave boundary points
      INTEGER                :: IWBMNPGL
      INTEGER, ALLOCATABLE   :: IWBNDLC(:) ! local wave boundary index
      INTEGER, ALLOCATABLE   :: IWBNDLC_REV(:) ! local wave boundary index
      INTEGER, ALLOCATABLE   :: IOBPD(:,:) ! boundary direction pointer
      INTEGER, ALLOCATABLE   :: IOBWB(:)   ! gl. wave boundary index ... will vanish in the decomp.
      INTEGER, ALLOCATABLE   :: IOBP(:)    ! boundary points index
      LOGICAL                :: LBCWA = .FALSE.
      CONTAINS
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_BOUNDARY
      IMPLICIT NONE
      ALLOCATE(IOBP(MNP)  ); IOBP = 0
      ALLOCATE(IOBPD(NANG,MNP)); IOBPD = 0
      ALLOCATE(IOBWB(MNP) ); IOBWB = 1 ! for boundary nodes we set this to 0 since in this case we omit advection at these nodes
      CALL SET_IOBP          ! boundary point marker
      CALL SET_IOBPD         ! boundary directional marker
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE RHEADER_NODE(IFILE,NKR,NKG)

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     XFN header nodes ...

!     Externals. 
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF* based on Aron Roland & Mathieu Dutour 

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IFILE
      INTEGER(KIND=JWIM), INTENT(OUT):: NKR, NKG

      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*) NKR
      READ(IFILE,*)
      READ(IFILE,*) NKG
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)

      END SUBROUTINE RHEADER_NODE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE RHEADER_ELEMENT(IFILE,NELEM)

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     XFN header elements

!     Externals.  
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF* based on Aron Roland & Mathieu Dutour 

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

      IMPLICIT NONE
      INTEGER(KIND=JWIM), INTENT(IN) :: IFILE
      INTEGER(KIND=JWIM), INTENT(OUT):: NELEM
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*) NELEM
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)
      END SUBROUTINE RHEADER_ELEMENT
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_STATUS(STATUS)
      IMPLICIT NONE
      integer, intent(out) :: STATUS(MNP)
      integer COLLECTED(MNP)
      integer PREVVERT(MNP)
      integer NEXTVERT(MNP)
      integer IE, I, IPREV, INEXT, IP
      integer IPPREV, IPNEXT
      integer ZNEXT
      LOGICAL IsFinished
      STATUS(:) = 0
      DO IE=1,MNE
        DO I=1,3
          IF (I.EQ.1) THEN
            IPREV=3
          ELSE
            IPREV=I-1
          END IF
          IF (I.EQ.3) THEN
            INEXT=1
          ELSE
            INEXT=I+1
          END IF
          IP=INE(I,IE)
          IPNEXT=INE(INEXT,IE)
          IPPREV=INE(IPREV,IE)
          IF (STATUS(IP).EQ.0) THEN
            STATUS(IP)=1
            PREVVERT(IP)=IPPREV
            NEXTVERT(IP)=IPNEXT
          END IF
        END DO
      END DO
      STATUS(:)=0
      DO
        COLLECTED(:)=0
        DO IE=1,MNE
          DO I=1,3
            IF (I.EQ.1) THEN
              IPREV=3
            ELSE
              IPREV=I-1
            END IF
            IF (I.EQ.3) THEN
              INEXT=1
            ELSE
              INEXT=I+1
            END IF
            IP=INE(I,IE)
            IPNEXT=INE(INEXT,IE)
            IPPREV=INE(IPREV,IE)
            IF (STATUS(IP).eq.0) THEN
              ZNEXT=NEXTVERT(IP)
              IF (ZNEXT.eq.IPPREV) THEN
                COLLECTED(IP)=1
                NEXTVERT(IP)=IPNEXT
                IF (NEXTVERT(IP).eq.PREVVERT(IP)) THEN
                  STATUS(IP)=1
                END IF
              END IF
            END IF
          END DO
        END DO
        ISFINISHED=.TRUE.
        DO IP=1,MNP
          IF ((COLLECTED(IP).eq.0).and.(STATUS(IP).eq.0)) THEN
            STATUS(IP)=-1
          END IF
          IF (STATUS(IP).eq.0) THEN
            ISFINISHED=.FALSE.
          END IF
        END DO
        IF (ISFINISHED) THEN
          EXIT
        END IF
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_IOBP
      USE YOWUNPOOL
      USE yowpd, only: ipgl, iplg, np_global
      IMPLICIT NONE
      INTEGER(KIND=JWIM) :: I, IWILD(MNP)
      INTEGER(KIND=JWIM) :: I1, I2, I3, IE, IP, ID, IFSTAT, IPglob
      INTEGER(KIND=JWIM) :: ZNEXT, ITMP
      INTEGER(KIND=JWIM) :: IOBPcopy(MNP), eDiff, nbDiff
      integer STATUS(MNP)
      integer istat
      integer, allocatable :: Indexes(:)
      REAL(KIND=JWRB) :: BNDTMP
      REAL(KIND=JWRU) :: x1, y1, x2, y2
      REAL(KIND=JWRU) :: EVX, EVY
      REAL(KIND=JWRU) :: eDet1, eDet2
      REAL(KIND=JWRU) :: ATMP, BTMP
      REAL(KIND=JWRU) :: rtemp(MNP)
      LOGICAL :: LFLIVE
      integer IOBPglobal(np_global)
      integer idx
      integer eIOBP
!
! open and read boundary nodes file ...
!
      IOBP    = 0
      IOBPD   = 0


      OPEN(BND%FHNDL, FILE = BND%FNAME, STATUS = 'OLD')
      CALL RHEADER_NODE(BND%FHNDL,ITMP,ITMP)
      DO IPglob = 1, np_global
        READ(BND%FHNDL,*) ITMP, BNDTMP, BNDTMP, BNDTMP
        eIOBP=INT(BNDTMP)
        IF (.NOT.  LBCWA) THEN
          IF ((eIOBP .EQ. 2).or.(eIOBP .EQ. 3)) eIOBP = 1
        END IF
        IOBPglobal(IPglob) = eIOBP
        IP = ipgl(IPglob)
        IF(IP /= 0) THEN
          IOBP(IP) = eIOBP
        END IF
      END DO
      CLOSE(BND%FHNDL)
      IOBPcopy=IOBP
!
! find islands and domain boundary ....
!
      CALL COMPUTE_STATUS(STATUS)
      DO IP=1,MNP
        IF (STATUS(IP).eq.-1 .AND. IOBP(IP) .EQ. 0) THEN
          IOBP(IP)=1
        END IF
      END DO
!
! reporting differences found
!
#ifdef DEBUG
      nbDiff=0
      DO IP=1,MNP
        IF (IOBP(IP) .ne. IOBPcopy(IP)) THEN
          nbDiff=nbDiff+1
        END IF
      END DO
      WRITE(740+MyRankGlobal,*) 'nbDiff=', nbDiff
#endif
!
! Determining number of boundary nodes
!
      IWBMNP = 0
      DO IP = 1, MNP
        IF (IOBP(IP) == 2) IWBMNP = IWBMNP + 1 ! Local number of boundary nodes ...
      END DO
!
! map boundary nodes ... needed later for the decomposition ...
!
      ALLOCATE(IWBNDLC(IWBMNP), IWBNDLC_REV(MNP), stat=istat)
      IWBNDLC_REV = 0
      idx = 0
      DO IP = 1, MNP
        IF (IOBP(IP) == 2) THEN
          idx = idx + 1
          IWBNDLC(idx)    = IP
          IWBNDLC_REV(IP) = idx
        END IF
      END DO
!
! allocate wave boundary arrays ... 
!
      ALLOCATE(Indexes(np_global), stat=istat)
      Indexes = 0
      idx = 0
      DO IPglob = 1, np_global
        IF (IOBPglobal(IPglob) == 2) THEN
          idx = idx + 1
          Indexes(IPglob)=idx
        END IF
      END DO
      IWBMNPGL = idx
      allocate(Indexes_boundary(IWBMNP), stat=istat)
      DO idx=1,IWBMNP
        IP=IWBNDLC(idx)
        IPglob=iplg(IP)
        Indexes_boundary(idx) = Indexes(IPglob)
      END DO
      IF (LBCWA) THEN
        CALL INIT_FILE_BOUNDARY
      END IF      
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_IOBPD

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Estimate the max. integration time step and amount of iterations ...

!     Externals.  Estimate spectral direction that are pointing into the domain from the boundary ...
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF* based on Aron Roland & Mathieu Dutour 

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

        USE YOWUNPOOL
        USE YOWFRED, ONLY : DFIM, COSTH, SINTH

        IMPLICIT NONE

        INTEGER(KIND=JWIM) :: I1, I2, I3, IE, IP, ID
        INTEGER(KIND=JWIM) :: I, IWILD(MNP)
        REAL(KIND=JWRU) :: DXP1, DXP2, DXP3, DYP1, DYP2, DYP3
        REAL(KIND=JWRU) :: x1, y1, x2, y2
        REAL(KIND=JWRU) :: EVX, EVY
        REAL(KIND=JWRU) :: eDet1, eDet2
        REAL(KIND=JWRU) :: rtemp(MNP)
!
! SET IOBPD ...
!
        DO IE=1,MNE
          I1   =   INE(1,IE)
          I2   =   INE(2,IE)
          I3   =   INE(3,IE)
          DXP1 =   IEN(6,IE)
          DYP1 = - IEN(5,IE)
          DXP2 =   IEN(2,IE)
          DYP2 = - IEN(1,IE)
          DXP3 =   IEN(4,IE)
          DYP3 = - IEN(3,IE)
!2do ... modifly wave direction by currents ...
          DO ID=1,MDC
            EVX=SINTH(ID)
            EVY=COSTH(ID)
            DO I=1,3
              IF (I.eq.1) THEN
                x1=   DXP1
                y1=   DYP1
                x2= - DXP3
                y2= - DYP3
                IP=   I1
              END IF
              IF (I.eq.2) THEN
                x1 =   DXP2
                y1 =   DYP2
                x2 = - DXP1
                y2 = - DYP1
                IP =   I2
              END IF
              IF (I.eq.3) THEN
                x1 =   DXP3
                y1 =   DYP3
                x2 = - DXP2
                y2 = - DYP2
                IP =   I3
              END IF
              eDet1 = SMALL-x1*EVY+y1*EVX
              eDet2 = SMALL+x2*EVY-y2*EVX
              IF ((eDet1.gt.0.0_JWRU).and.(eDet2.gt.0.0_JWRU)) THEN
                IOBPD(ID,IP)=1
              END IF
            END DO
          END DO
        END DO

        DO IP = 1, MNP
          IF ( LBCWA ) THEN
            IF ( IOBP(IP) == 2 ) THEN
              IOBWB(IP) = 0
              IOBPD(:,IP) = 1
            ENDIF
          END IF
          IF ( IOBP(IP) == 3 ) THEN ! If Neumann boundary condition is given set IOBP to 4
            IOBPD(:,IP) = 1 ! Update Neumann nodes ...
          END IF
        END DO

        DO ID=1, MDC
          rtemp = IOBPD(ID,:)
          call exchange(rtemp)
          IOBPD(ID,:) = INT(rtemp)
        END DO
      END SUBROUTINE SET_IOBPD
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SINGLE_READ_BOUNDARY(eFile, WBAC, IT)
      USE YOWPCONS , ONLY : ZPI
      USE NETCDF
      IMPLICIT NONE
      character(len=*), intent(in) :: eFile
      REAL(KIND=JWRU), intent(out) :: WBAC(NANG,NFRE,IWBMNP)
      integer, intent(in) :: IT
      REAL(KIND=JWRU) WBAC_GL(NFRE,NANG,IWBMNPGL)
      REAL(KIND=JWRU) eAC, eFL
      integer istat, ncid, var_id
      integer IP, IS, ID, idx
      character(len=*), parameter :: CallFct = "SINGLE_READ_BOUNDARY"
      integer ID1, ID2
      REAL(KIND=JWRU) eWD1, eWD2
      !
      ! We have this inversion of the order because that is so in WWM at 
      ! the present time
      !
      ISTAT = NF90_OPEN(TRIM(eFile), NF90_NOWRITE, ncid)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)
      ISTAT = nf90_inq_varid(ncid, 'WBAC', var_id)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)
      ISTAT = NF90_GET_VAR(ncid, var_id, WBAC_GL, start=(/1,1,1,IT/), count = (/NFRE,NANG,IWBMNPGL,1/))
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)
      ISTAT = NF90_CLOSE(ncid)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)
      !
      ! Now reassigning 
      !
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'IWBMNPGL = ', IWBMNPGL
      WRITE(740+MyRankGlobal,*) 'IWBMNP   = ', IWBMNP
#endif
      DO IP=1,IWBMNP
        idx=Indexes_boundary(IP)
!        WRITE(740+MyRankGlobal,*) 'IP=', IP, ' idx=', idx
        DO ID=1,NANG
           DO IS=1,NFRE
            ID1=WWM_ID1(ID)
            ID2=WWM_ID2(ID)
            eWD1=WWM_WD1(ID)
            eWD2=WWM_WD2(ID)
            eAC=WBAC_GL(IS,ID1,idx) * eWD1 + WBAC_GL(IS,ID2,idx) * eWD2
            eFL=eAC * SPSIG(IS) * ZPI
            WBAC(ID,IS,IP) = eFL
!#ifdef DEBUG
!            WRITE(740+MyRankGlobal,*) 'ID=', ID, ' IS=', IS, ' IP=', IP
!            WRITE(740+MyRankGlobal,*) '   wbac=', WBAC(ID,IS,IP)
!#endif
          END DO
        END DO
      END DO
!#ifdef DEBUG
!      FLUSH(740+MyRankGlobal)
!#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_UP_WBAC
      USE YOWSTAT  , ONLY : IDELT
      IMPLICIT NONE
      REAL*8 eTimeDay
      REAL(KIND=JWRU) w1, w2
      integer iTime1, iTime2
      CALL WAV_GET_ETIMEDAY(eTimeDay, WAV_BoucTime)
      CALL FIND_MATCH_TIME(RecTimeBnd, eTimeDay, iTime1, w1, iTime2, w2)
      IF (iTime1 .ne. recTime1) THEN
        CALL SINGLE_READ_BOUNDARY(eFileBnd, WBAC1, iTime1)
        recTime1 = iTime1
      END IF
      IF (iTime2 .ne. recTime2) THEN
        CALL SINGLE_READ_BOUNDARY(eFileBnd, WBAC2, iTime2)
        recTime2 = iTime2
      END IF
      WBAC = w1*WBAC1 + w2*WBAC2
      WAV_BoucTime = WAV_BoucTime + DBLE(IDELT)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE APPLY_BOUNDARY_CONDITION(FL)
      USE YOWMPP   , ONLY : NINF, NSUP
      IMPLICIT NONE
      REAL(KIND=JWRB), INTENT(INOUT) :: FL(NINF-1:NSUP,NANG,NFRE)
      integer ID, IS, IP, idx
      DO IS=1,NFRE
        DO ID=1,NANG
          DO idx=1,IWBMNP
            IP=IWBNDLC(idx)
            FL(IP,ID,IS) = REAL(WBAC(ID,IS,idx),JWRB)
          END DO
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_FILE_BOUNDARY
      USE YOWPCONS, ONLY : PI, ZPI
      USE YOWFRED, ONLY : FR, TH
      USE NETCDF
      IMPLICIT NONE
      character(len=*), parameter :: CallFct = "INIT_FILE_BOUNDARY"
      INTEGER varid, ncid
      integer, dimension(nf90_max_var_dims) :: dimids
      integer istat
      integer iTime
      real(KIND=JWRU) eWD1, eWD2, eDiff, eDiff1, eDiff2
      real(KIND=JWRU) DeltaDiff, eDir
      logical IsAssigned
      integer ID, ID1, ID2
      !
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'eFileBnd=', TRIM(eFileBnd)
#endif
      !
      CALL TEST_FILE_EXIST_DIE(eFileBnd, "Need boundary file for boundary forcing")
      ISTAT = NF90_OPEN(TRIM(eFileBnd), NF90_NOWRITE, ncid)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)
      !
      ! reading direction
      !
      ISTAT = nf90_inq_varid(ncid, "SPDIR", varid)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)
      !
      ISTAT = nf90_inquire_variable(ncid, varid, dimids=dimids)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)
      !
      ISTAT = nf90_inquire_dimension(ncid, dimids(1), len = nbDirWWM)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)
      allocate(SPDIR_WWM(nbDirWWM), stat=istat)
      !
      ISTAT = nf90_get_var(ncid, varid, SPDIR_WWM)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)
      !
      ! reading frequencies
      !
      ISTAT = nf90_inq_varid(ncid, "SPSIG", varid)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)
      !
      ISTAT = nf90_inquire_variable(ncid, varid, dimids=dimids)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 7, ISTAT)
      !
      ISTAT = nf90_inquire_dimension(ncid, dimids(1), len = nbFreqWWM)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 8, ISTAT)
      allocate(SPSIG_WWM(nbFreqWWM), stat=istat)
      !
      ISTAT = nf90_get_var(ncid, varid, SPSIG_WWM)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 9, ISTAT)
      !
      ! reading the time
      !
      CALL READ_LIST_TIME(ncid, RecTimeBnd)
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
      ALLOCATE(WBAC(NANG,NFRE,IWBMNP), WBAC1(NANG,NFRE,IWBMNP), WBAC2(NANG,NFRE,IWBMNP), stat=istat)
      !
      ! allocating the interpolation arrays
      !
      allocate(WWM_ID1(NANG), WWM_ID2(NANG), WWM_WD1(NANG), WWM_WD2(NANG), stat=istat)
      WWM_ID1=0
      WWM_ID2=0
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'NANG=', NANG, 'nbDirWWM=', nbDirWWM
#endif
      DO ID=1,NANG
        eDir = PI *0.5_JWRU - TH(ID)
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) '--------------------------------------------------'
        WRITE(740+MyRankGlobal,*) 'ID=', ID
        WRITE(740+MyRankGlobal,*) 'eDir=', eDir, ' TH=', TH(ID)
#endif
        IsAssigned=.FALSE.
        DO ID1=1,nbDirWWM
          IF (ID1 .lt. nbDirWWM) THEN
            ID2=ID1 + 1
          ELSE
            ID2=1
          END IF
          IF (IsAssigned .eqv. .FALSE.) THEN
            eDiff = SPDIR_WWM(ID2) - SPDIR_WWM(ID1)
            CALL RenormalizeAngle(eDiff)
            !
            eDiff1 = eDir - SPDIR_WWM(ID1)
            CALL RenormalizeAngle(eDiff1)
            !
            eDiff2 = SPDIR_WWM(ID2) - eDir
            CALL RenormalizeAngle(eDiff2)
            !
            DeltaDiff=abs(eDiff) - abs(eDiff1) - abs(eDiff2)
#ifdef DEBUG
            WRITE(740+MyRankGlobal,*) 'ID1=', ID1, ' ID2=', ID2
            WRITE(740+MyRankGlobal,*) 'SPDIR12=', SPDIR_WWM(ID1), SPDIR_WWM(ID2)
            WRITE(740+MyRankGlobal,*) 'eDiff=', eDiff
            WRITE(740+MyRankGlobal,*) 'eDiff1=', eDiff1, 'eDiff2=', eDiff2
            WRITE(740+MyRankGlobal,*) 'DeltaDiff=', DeltaDiff
#endif
            IF (abs(DeltaDiff) .lt. 0.0001_JWRU) THEN
#ifdef DEBUG
              WRITE(740+MyRankGlobal,*) 'M A T C H I N G'
#endif
              eWD1 = eDiff2 / eDiff
              eWD2 = eDiff1 / eDiff
              IsAssigned=.TRUE.
              WWM_ID1(ID) = ID1
              WWM_ID2(ID) = ID2
              WWM_WD1(ID) = eWD1
              WWM_WD2(ID) = eWD2
#ifdef DEBUG
              WRITE(740+MyRankGlobal,*) 'eWD1=', eWD1, ' eWD2=', eWD2
#endif
            END IF
          END IF
        END DO
        IF (IsAssigned .eqv. .FALSE.) THEN
          FLUSH(740+MyRankGlobal)
          Print *, 'Failed in the directional interpolation'
          STOP
        END IF
#ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'ID12=', WWM_ID1(ID), WWM_ID2(ID)
        WRITE(740+MyRankGlobal,*) 'WD12=', WWM_WD1(ID), WWM_WD2(ID)
        WRITE(740+MyRankGlobal,*) '--------------------------------------------------'
#endif
      END DO
      CONTAINS
      SUBROUTINE RenormalizeAngle(eDiff)
      REAL(KIND=JWRU), intent(inout) :: eDiff
      IF (eDiff .gt. PI) THEN
        eDiff = eDiff - ZPI
      END IF
      IF (eDiff .lt. -PI) THEN
        eDiff = eDiff + ZPI
      END IF
      END SUBROUTINE
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
END MODULE
