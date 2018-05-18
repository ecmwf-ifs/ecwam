MODULE WAV_netcdf
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
#ifdef NETCDF_OUTPUT_WAM
      USE WAV_netcdf_var
      USE WAV_NETCDF_FCT, ONLY : WAV_GENERIC_NETCDF_ERROR
      CONTAINS
!**********************************************************************
!*                                                                    *
!**********************************************************************
# ifdef DEBUG
      SUBROUTINE PRINT_MINMAX_U10_NEWOLD
      USE YOWGRID  , ONLY : IJS      ,IJL, IJSLOC, IJLLOC
      USE YOWSPEC, ONLY   : U10NEW   ,U10OLD
      implicit none
      REAL(KIND=JWRB) :: siz, avgHS, avgWindSpeed
      logical IsFirst
      REAL(KIND=JWRB) :: minU10new, maxU10new
      REAL(KIND=JWRB) :: minU10old, maxU10old
      INTEGER(KIND=JWIM) :: IG, IJ
      IG=1
      IsFirst=.TRUE.
      DO IJ=IJSLOC,IJLLOC
        IF (IsFirst) THEN
          minU10new=U10NEW(IJ)
          maxU10new=U10NEW(IJ)
          minU10old=U10OLD(IJ,IG)
          maxU10old=U10OLD(IJ,IG)
        ELSE
          IF (minU10new .gt. U10NEW(IJ) ) THEN
            minU10new=U10NEW(IJ)
          END IF
          IF (maxU10new .lt. U10NEW(IJ) ) THEN
            maxU10new=U10NEW(IJ)
          END IF
          IF (minU10old .gt. U10OLD(IJ,IG) ) THEN
            minU10old=U10OLD(IJ,IG)
          END IF
          IF (maxU10old .lt. U10OLD(IJ,IG) ) THEN
            maxU10old=U10OLD(IJ,IG)
          END IF
        END IF
        IsFirst=.FALSE.
      END DO
      WRITE(740+MyRankGlobal,*) 'IJSLOC/IJLLOC=', IJSLOC, IJLLOC
      WRITE(740+MyRankGlobal,*) 'U10NEW(min/max)sel=', minU10new, maxU10new
      WRITE(740+MyRankGlobal,*) 'U10OLD(min/max)sel=', minU10old, maxU10old
      FLUSH(740+MyRankGlobal)
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAV_netcdf_init
      USE YOWGRID  , ONLY : IJS      ,IJL, IJSLOC, IJLLOC, IJGLOBAL_OFFSET
      USE YOWMPP   , ONLY : NPROC, NPRECR, NINF, NSUP, IRANK
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWMAP   , ONLY : IXLG     ,KXLT
      USE yowpd    , only : np_global, MNP=>npa, NP_RES => np
      USE YOWUNPOOL, ONLY : LLUNSTR, LCFL
      USE yownodepool, ONLY : iplg
!# if !defined MODEL_COUPLING_ATM_WAV && !defined MODEL_COUPLING_OCN_WAV
      USE MPL_MPIF 
!# endif
      implicit none
      INTEGER(KIND=JWIM) :: rbuf_int(3), ierr, istat
      INTEGER(KIND=JWIM) :: IG, idx_loc
      INTEGER(KIND=JWIM) :: IX, IY, IXY, IPROC, NB_loc, IJ
      REAL(KIND=JWRB), allocatable :: dspl_recv(:)
      INTEGER(KIND=JWIM) :: status(MPI_STATUS_SIZE)
      INTEGER(KIND=JWIM), allocatable :: AttainedMatrix(:,:), eInt(:)
      INTEGER(KIND=JWIM) :: MNPloc, NP_RESloc, IP, idx
# if defined MODEL_COUPLING_ATM_WAV || defined MODEL_COUPLING_OCN_WAV
      MPI_COMM_NETCDF=WAV_COMM_WORLD
# else
      MPI_COMM_NETCDF=MPI_COMM_WORLD
# endif
      IG=1
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'MyRankLocal=', MyRankLocal
      WRITE(740+MyRankGlobal,*) 'IRANK=', IRANK
      FLUSH(740+MyRankGlobal)
# endif
      idx=11
# if defined MODEL_COUPLING_OCN_WAV
      NETCDF_nbVar=NETCDF_nbVar + 3
      idx=idx+1
      idxUcurr=idx
      !
      idx=idx+1
      idxVcurr=idx
      !
      idx=idx+1
      idxZeta=idx
# endif
      IF (LLUNSTR .and. LCFL) THEN
        NETCDF_nbVar=NETCDF_nbVar + 3
        idx=idx+1
        idxcfl1=idx
        idx=idx+1
        idxcfl2=idx
        idx=idx+1
        idxcfl3=idx
      END IF
      IF (LLUNSTR) THEN
        NETCDF_X=np_global
        NETCDF_Y=1
        IF (IRANK .eq. 1) THEN
          allocate(ListMNPloc(nproc), ListNP_RESloc(nproc))
          ListMNPloc(1)=MNP
          ListNP_RESloc(1)=NP_RES
          DO iProc=2,nproc
            allocate(eInt(2), stat=istat)
            CALL MPI_RECV(eInt,2,MPI_INTEGER(KIND=JWIM),iProc-1, 53, MPI_COMM_NETCDF, status,ierr)
# ifdef DEBUG
            WRITE(740+MyRankGlobal,*) 'After MPI_RECV 53, ierr=', ierr
            FLUSH(740+MyRankGlobal)
# endif
            MNPloc=eInt(1)
            NP_RESloc=eInt(2)
            deallocate(eInt)
            ListMNPloc(iProc)=MNPloc
            ListNP_RESloc(iProc)=NP_RESloc
          END DO
          MaxMNPloc=maxval(ListMNPloc)
          allocate(ListIPLGloc(MaxMNPloc, nproc))
          DO IP=1,MNP
            ListIPLGloc(IP,1)=iplg(IP)
          END DO
          DO iProc=2,nproc
            MNPloc=ListMNPloc(iProc)
            allocate(eInt(MNPloc), stat=istat)
            CALL MPI_RECV(eInt,MNPloc,MPI_INTEGER(KIND=JWIM),iProc-1, 54, MPI_COMM_NETCDF, status,ierr)
# ifdef DEBUG
            WRITE(740+MyRankGlobal,*) 'After MPI_RECV 54, ierr=', ierr
            FLUSH(740+MyRankGlobal)
# endif
            DO IP=1,MNPloc
              ListIPLGloc(IP,iProc)=eInt(IP)
            END DO
            deallocate(eInt)
          END DO
        ELSE
          allocate(eInt(2))
          eInt(1)=MNP
          eInt(2)=NP_RES
          CALL MPI_SEND(eInt,2,MPI_INTEGER(KIND=JWIM), 0, 53, MPI_COMM_NETCDF, ierr)
# ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'After MPI_SEND 53, ierr=', ierr
          FLUSH(740+MyRankGlobal)
# endif
          deallocate(eInt)
          !
          allocate(eInt(MNP))
          DO IP=1,MNP
            eInt(IP)=iplg(IP)
          END DO
          CALL MPI_SEND(eInt,MNP,MPI_INTEGER(KIND=JWIM), 0, 54, MPI_COMM_NETCDF, ierr)
# ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'After MPI_SEND 54, ierr=', ierr
          FLUSH(740+MyRankGlobal)
# endif
          deallocate(eInt)
        END IF
        NB_loc=NP_RES
      ELSE
        NETCDF_X=NGX
        NETCDF_Y=NGY
        allocate(AttainedMatrix(NETCDF_X, NETCDF_Y))
        allocate(ListIJS(NPROC), ListIJL(NPROC), ListIJ_OFFSET(NPROC))
        IF (IRANK .eq. 1) THEN
          ListIJS(1)=IJSLOC
          ListIJL(1)=IJLLOC
          ListIJ_OFFSET(1)=IJGLOBAL_OFFSET
          DO iProc=2,NPROC
            CALL MPI_RECV(rbuf_int,3,MPI_INT, iProc-1, 1958, MPI_COMM_NETCDF, status, ierr)
# ifdef DEBUG
            WRITE(740+MyRankGlobal,*) 'After MPI_RECV 1958, ierr=', ierr
            FLUSH(740+MyRankGlobal)
# endif
            ListIJS(iProc)=rbuf_int(1)
            ListIJL(iProc)=rbuf_int(2)
            ListIJ_OFFSET(iProc)=rbuf_int(3)
          END DO
          DO iProc=2,NPROC
            CALL MPI_SEND(ListIJS,NPROC,MPI_INT, iProc-1, 1957, MPI_COMM_NETCDF, ierr)
# ifdef DEBUG
            WRITE(740+MyRankGlobal,*) 'After MPI_SEND 1957, ierr=', ierr
            FLUSH(740+MyRankGlobal)
# endif
            CALL MPI_SEND(ListIJL,NPROC,MPI_INT, iProc-1, 1959, MPI_COMM_NETCDF, ierr)
# ifdef DEBUG
            WRITE(740+MyRankGlobal,*) 'After MPI_SEND 1959, ierr=', ierr
            FLUSH(740+MyRankGlobal)
# endif
            CALL MPI_SEND(ListIJ_OFFSET,NPROC,MPI_INT, iProc-1, 1961, MPI_COMM_NETCDF, ierr)
# ifdef DEBUG
            WRITE(740+MyRankGlobal,*) 'After MPI_SEND 1961, ierr=', ierr
            FLUSH(740+MyRankGlobal)
# endif
          END DO
        ELSE
          rbuf_int(1)=IJSLOC
          rbuf_int(2)=IJLLOC
          rbuf_int(3)=IJGLOBAL_OFFSET
          CALL MPI_SEND(rbuf_int,3,MPI_INT, 0, 1958, MPI_COMM_NETCDF, ierr)
# ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'After MPI_SEND 1958, ierr=', ierr
          FLUSH(740+MyRankGlobal)
# endif
          CALL MPI_RECV(ListIJS,NPROC,MPI_INT, 0, 1957, MPI_COMM_NETCDF, status, ierr)
# ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'After MPI_RECV 1957, ierr=', ierr
          FLUSH(740+MyRankGlobal)
# endif
          CALL MPI_RECV(ListIJL,NPROC,MPI_INT, 0, 1959, MPI_COMM_NETCDF, status, ierr)
# ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'After MPI_RECV 1959, ierr=', ierr
          FLUSH(740+MyRankGlobal)
# endif
          CALL MPI_RECV(ListIJ_OFFSET,NPROC,MPI_INT, 0, 1961, MPI_COMM_NETCDF, status, ierr)
# ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'After MPI_RECV 1961, ierr=', ierr
          FLUSH(740+MyRankGlobal)
# endif
        END IF
        IF (MyRankLocal .eq. 0) THEN
          allocate(NETCDF_var_rqst(NPROC-1), NETCDF_var_stat(MPI_STATUS_SIZE,NPROC-1), NETCDF_var_type(NPROC-1))
          AttainedMatrix=0
          DO IPROC=2,NPROC
            NB_loc=1+ListIJL(IPROC) - ListIJS(iPROC)
            allocate(dspl_recv(NB_loc))
            idx_loc=0
            DO IJ=ListIJS(IPROC),ListIJL(IPROC)
              idx_loc=idx_loc+1
              IX = IXLG(IJ,IG)
              IY = NGY- KXLT(IJ,IG) +1
              IXY= IX + NGX*(IY-1)
              AttainedMatrix(IX,IY)=AttainedMatrix(IX,IY) + 1
              dspl_recv(idx_loc)=NETCDF_nbVar*(IXY - 1)
            END DO
            call mpi_type_create_indexed_block(NB_loc,NETCDF_nbVar,dspl_recv,MPI_REAL4,NETCDF_var_type(IPROC-1), ierr)
            call mpi_type_commit(NETCDF_var_type(IPROC-1), ierr)
            deallocate(dspl_recv)
          END DO
        END IF
        NB_loc=1 + IJLLOC - IJSLOC
      END IF
      allocate(NETCDF_var_gridded(NETCDF_nbVar, NETCDF_X,NETCDF_Y))
      NETCDF_var_gridded=0
      allocate(NETCDF_var(NETCDF_nbVar, NB_loc))
      NETCDF_var=0
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'End WAV_netcdf_init'
      WRITE(740+MyRankGlobal,*) 'NETCDF_X=', NETCDF_X
      WRITE(740+MyRankGlobal,*) 'NETCDF_Y=', NETCDF_Y
      FLUSH(740+MyRankGlobal)
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAV_netcdf_setup_array_variables_V1
      USE YOWGRID  , ONLY : IJS      ,IJL, IJSLOC, IJLLOC, IJGLOBAL_OFFSET
      USE YOWMPP   , ONLY : NPROC, NPRECR, NINF, NSUP, IRANK
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWMAP   , ONLY : IXLG     ,KXLT
      USE YOWUNPOOL, ONLY : LLUNSTR
!#if !defined MODEL_COUPLING_ATM_WAV && !defined MODEL_COUPLING_OCN_WAV
      USE MPL_MPIF
!#endif
      implicit none
      INTEGER(KIND=JWIM) :: IG, idx_loc, IX, IY, IXY, IJ
      INTEGER(KIND=JWIM) :: IPROC, NB_loc, IJglob
      INTEGER(KIND=JWIM) :: ierr, istatus, nbWait
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'Begin WAV_netcdf_setup_array_variables'
      FLUSH(740+MyRankGlobal)
# endif
      IF (LLUNSTR) THEN
        Print *, 'Need to write the unstructured code'
        Print *, 'anyway the structured code is bugged as well'
        STOP
      END IF
      IG=1
      NB_loc=1 + IJLLOC - IJSLOC
      IF (IRANK .eq. 1) THEN
        nbWait=NPROC-1
        NETCDF_var_rqst=0
        DO IPROC=1,nbWait
# ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'before IPROC=', IPROC, 'rqst=', NETCDF_var_rqst(IPROC)
          FLUSH(740+MyRankGlobal)
# endif
          call mpi_irecv(NETCDF_var_gridded,1,NETCDF_var_type(IPROC),iPROC,1968,MPI_COMM_NETCDF,NETCDF_var_rqst(IPROC),ierr)
# ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'mpi_irecv, ierr=', ierr
          WRITE(740+MyRankGlobal,*) 'after IPROC=', IPROC, 'rqst=', NETCDF_var_rqst(IPROC)
          FLUSH(740+MyRankGlobal)
# endif
        END DO
        idx_loc=0
        DO IJ=IJSLOC,IJLLOC
          idx_loc=idx_loc+1
          IJglob = IJ + IJGLOBAL_OFFSET
          IX = IXLG(IJ,IG)
          IY = NGY- KXLT(IJ,IG) +1
          IXY= IX + NGX*(IY-1)
          NETCDF_var_gridded(:,IX,IY)=NETCDF_var(:,idx_loc)
        END DO
# ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'allocated(NETCDF_var_rqst)=', allocated(NETCDF_var_rqst)
        WRITE(740+MyRankGlobal,*) 'allocated(NETCDF_var_stat)=', allocated(NETCDF_var_stat)
        WRITE(740+MyRankGlobal,*) 'size(NETCDF_var_rqst)=', size(NETCDF_var_rqst)
        WRITE(740+MyRankGlobal,*) 'size(NETCDF_var_stat)=', size(NETCDF_var_stat)
        WRITE(740+MyRankGlobal,*) 'MPI_STATUS_SIZE=', MPI_STATUS_SIZE
        FLUSH(740+MyRankGlobal)
# endif
        IF (nbWait > 0) THEN
# ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'Before MPI_waitall'
          WRITE(740+MyRankGlobal,*) 'NPROC=', NPROC
          WRITE(740+MyRankGlobal,*) 'nbWait=', nbWait
          FLUSH(740+MyRankGlobal)
# endif
          call mpi_waitall(nbWait, NETCDF_var_rqst, MPI_STATUS_IGNORE, ierr)
!          call mpi_waitall(nbWait, NETCDF_var_rqst, NETCDF_var_stat, ierr)
# ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'mpi_waitall, ierr=', ierr
          FLUSH(740+MyRankGlobal)
# endif
        END IF
      ELSE
# ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'Before MPI_SEND data 1958'
        WRITE(740+MyRankGlobal,*) 'NB_loc=', NB_loc
        FLUSH(740+MyRankGlobal)
# endif
        CALL MPI_SEND(NETCDF_var,NETCDF_nbVar*NB_loc,MPI_REAL4, 0, 1968, MPI_COMM_NETCDF, ierr)
# ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'After MPI_SEND data 1958, ierr=', ierr
        FLUSH(740+MyRankGlobal)
# endif
      END IF
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'After WAV_netcdf_setup_array_variables'
      FLUSH(740+MyRankGlobal)
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAV_netcdf_setup_array_variables
      USE YOWGRID  , ONLY : IJS      ,IJL, IJSLOC, IJLLOC, IJGLOBAL_OFFSET
      USE YOWMPP   , ONLY : NPROC, NPRECR, NINF, NSUP, IRANK
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWMAP   , ONLY : IXLG     ,KXLT
      USE YOWUNPOOL, ONLY : LLUNSTR
      USE yownodepool, ONLY : iplg
      USE yowpd    , only : np_global, MNP=>npa, NP_RES => np
!#if !defined MODEL_COUPLING_ATM_WAV && !defined MODEL_COUPLING_OCN_WAV
      USE MPL_MPIF 
!#endif
      implicit none
      INTEGER(KIND=JWIM) :: IG, idx_loc, IX, IY, IXY, IJ, IP
      INTEGER(KIND=JWIM) :: IPROC, NB_loc, NB_loc_side
      INTEGER(KIND=JWIM) :: ierr, istatus, nbWait
      REAL(KIND=JWRB), allocatable :: TheRecv(:,:)
      INTEGER(KIND=JWIM) :: status(MPI_STATUS_SIZE)
      INTEGER(KIND=JWIM) :: iVar, IPglob, NP_RESloc, IJglob
      INTEGER(KIND=JWIM) :: VertStatus(np_global)
      REAL(KIND=JWRB) :: eMin, eMax, eAvg
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'Before WAV_netcdf_setup_array_variables'
      WRITE(740+MyRankGlobal,*) 'minWind min(NETCDF_var(3,:))=', minval(NETCDF_var(3,:))
      WRITE(740+MyRankGlobal,*) 'NPROC=', NPROC
      WRITE(740+MyRankGlobal,*) 'NETCDF_nbVar=', NETCDF_nbVar
      FLUSH(740+MyRankGlobal)
# endif
      IF (LLUNSTR) THEN
        IF (IRANK .eq. 1) THEN
          NETCDF_var_gridded=0
          VertStatus=0
          DO IP=1,NP_RES
            IPglob=iplg(IP)
            VertStatus(IPglob)=1
            NETCDF_var_gridded(:,IPglob,1)=NETCDF_var(:,IP)
          END DO
          DO IPROC=2,NPROC
            NP_RESloc=ListNP_RESloc(IPROC)
            allocate(TheRecv(NETCDF_nbVar,NP_RESloc))
            CALL MPI_RECV(TheRecv,NETCDF_nbVar*NP_RESloc,MPI_REAL4,IPROC-1,1968,MPI_COMM_NETCDF,status,ierr)
# ifdef DEBUG
            WRITE(740+MyRankGlobal,*) 'IPROC=', IPROC
            WRITE(740+MyRankGlobal,*) 'NP_RESloc=', NP_RESloc
            WRITE(740+MyRankGlobal,*) 'MPI_RECV, 1968, ierr=', ierr
            WRITE(740+MyRankGlobal,*) 'minval(TheRecv(3,:))=', minval(TheRecv(3,:))
            FLUSH(740+MyRankGlobal)
# endif
            DO IP=1,NP_RESloc
              IPglob=ListIPLGloc(IP,IPROC)
              VertStatus(IPglob)=1
              NETCDF_var_gridded(:,IPglob,1)=TheRecv(:,IP)
            END DO
            deallocate(TheRecv)
          END DO
# ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'min(VertStatus)=', minval(VertStatus)
          FLUSH(740+MyRankGlobal)
# endif
# ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'minWind min(NETCDF_var_gridded(3,:))=', minval(NETCDF_var_gridded(3,:,:))
          FLUSH(740+MyRankGlobal)
# endif
        ELSE
# ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'Before MPI_SEND of NETCDF_var'
          WRITE(740+MyRankGlobal,*) 'minval(NETCDF_var(3,:))=', minval(NETCDF_var(3,:))
          WRITE(740+MyRankGlobal,*) 'NP_RES=', NP_RES
          FLUSH(740+MyRankGlobal)
# endif
          CALL MPI_SEND(NETCDF_var,NETCDF_nbVar*NP_RES,MPI_REAL4, 0, 1968, MPI_COMM_NETCDF, ierr)
# ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'MPI_SEND, 1968, ierr=', ierr
          FLUSH(740+MyRankGlobal)
# endif
        END IF
      ELSE
        IG=1
        NB_loc=1+IJLLOC - IJSLOC
        IF (IRANK .eq. 1) THEN
          NETCDF_var_gridded=0
          idx_loc=0
          DO IJ=IJSLOC,IJLLOC
            idx_loc=idx_loc+1
            IJglob=IJ + IJGLOBAL_OFFSET
            IX = IXLG(IJglob,IG)
            IY = NGY- KXLT(IJglob,IG) +1
            NETCDF_var_gridded(:,IX,IY)=NETCDF_var(:,idx_loc)
          END DO
          DO IPROC=2,NPROC
            NB_loc_side=1 + ListIJL(IPROC) - ListIJS(IPROC)
            allocate(TheRecv(NETCDF_nbVar,NB_loc_side))
            CALL MPI_RECV(TheRecv,NETCDF_nbVar*NB_loc_side,MPI_REAL4,IPROC-1,1968,MPI_COMM_NETCDF,status,ierr)
            idx_loc=0
            DO IJ=ListIJS(IPROC),ListIJL(IPROC)
              idx_loc=idx_loc+1
              IJglob=IJ + ListIJ_OFFSET(IPROC)
              IX = IXLG(IJglob,IG)
              IY = NGY- KXLT(IJglob,IG) +1
              NETCDF_var_gridded(:,IX,IY)=TheRecv(:,idx_loc)
            END DO
            deallocate(TheRecv)
          END DO
        ELSE
          CALL MPI_SEND(NETCDF_var,NETCDF_nbVar*NB_loc,MPI_REAL4, 0, 1968, MPI_COMM_NETCDF, ierr)
        END IF
      END IF
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'After WAV_netcdf_setup_array_variables'
      FLUSH(740+MyRankGlobal)
# endif
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'MyRankLocal=', MyRankLocal
      IF (MyRankLocal .eq. 0) THEN
        WRITE(740+MyRankGlobal,*) 'setup_array_variable : NETCDF_nbVar=', NETCDF_nbVar
        DO iVar=1,NETCDF_nbVar
          eMin=minval(NETCDF_var_gridded(iVar,:,:))
          eMax=maxval(NETCDF_var_gridded(iVar,:,:))
          eAvg=sum(NETCDF_var_gridded(iVar,:,:))/(REAL(NETCDF_X)*REAL(NETCDF_Y))
          WRITE(740+MyRankGlobal,*) 'iV=', iVar, ' v(min/max/sum)=', eMin, eMax, eAvg
        END DO
        WRITE(740+MyRankGlobal,*) 'Begin WAV_netcdf_setup_array_variables'
        WRITE(740+MyRankGlobal,*) 'From NETCDF_var_gridded'
        WRITE(740+MyRankGlobal,*) 'MaxWind=', maxval(NETCDF_var_gridded(3,:,:))
        WRITE(740+MyRankGlobal,*) 'MinWind=', minval(NETCDF_var_gridded(3,:,:))
        FLUSH(740+MyRankGlobal)
      END IF
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAV_netcdf_export
      implicit none
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'WAV_netcdf_export, begin'
      WRITE(740+MyRankGlobal,*) 'DoNETCDF_sync=', DoNETCDF_sync
      FLUSH(740+MyRankGlobal)
# endif
      IF (DoNETCDF_sync) THEN
# ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'Before WAV_netcdf_init'
        WRITE(740+MyRankGlobal,*) 'NETCDF_initialized=', NETCDF_initialized
        FLUSH(740+MyRankGlobal)
# endif
        IF (NETCDF_initialized .eqv. .FALSE.) THEN
          CALL WAV_netcdf_init
        END IF
# ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'After WAV_netcdf_init'
        FLUSH(740+MyRankGlobal)
# endif
        NETCDF_initialized=.TRUE.
        CALL WAV_assign_netcdf_array
# ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'After WAV_assign_netcdf_array'
        FLUSH(740+MyRankGlobal)
# endif
        CALL WAV_netcdf_setup_array_variables
      END IF
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'WAV_netcdf_export, end'
      FLUSH(740+MyRankGlobal)
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAV_assign_netcdf_array
      USE YOWGRID  , ONLY : IJS      ,IJL, IJSLOC, IJLLOC
      USE YOWSPEC, ONLY   : NSTART   ,NEND     ,                     &
     &            U10NEW   ,U10OLD   ,THWNEW   ,THWOLD   ,USNEW    , &
     &            USOLD    ,Z0NEW    ,Z0OLD    ,TAUW     ,           &
     &            ROAIRN   ,ROAIRO   ,ZIDLNEW  ,ZIDLOLD  ,           &
     &            FL3
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,THQ,FPMEAN
      USE YOWINTP  , ONLY : WHGTTG   ,WDIRTG   ,WPKFTG   ,WMNFTG,     &
     &            USTARG, CDG
      USE YOWPCONS , ONLY : EPSUS    ,EPSU10
      USE YOWUNPOOL, ONLY : LLUNSTR, LCFL, CFLCXY
      USE yowpd    , only : NP_RES => np
#if defined MODEL_COUPLING_ATM_WAV
      USE pgmcl_lib_WAM, only : Uwind_atm, Vwind_atm
#endif
#if defined MODEL_COUPLING_OCN_WAV
      USE pgmcl_lib_WAM, only : Uwind_ocn, Vwind_ocn, ZETA_ocn, Ucurr_ocn, Vcurr_ocn
      USE coupling_var, only : NlevelVert
#endif
      implicit none
      INTEGER(KIND=JWIM) :: IG, idx_loc, IJ
      INTEGER(KIND=JWIM) :: idx
      REAL(KIND=JWRB) :: eSingZ0, eUstarE, gzValue, eAlpha, eHS, eDir
      REAL(KIND=JWRB) :: ePeakFreq, eMeanFreq, eTAU, eCD
      REAL(KIND=JWRB) :: eU_10, eV_10, eWndMag, gValue
      REAL(KIND=JWRB) :: eWindSpeed, eWindDir
      REAL(KIND=JWRB) :: diffU, diffV, diffUV
      REAL(KIND=JWRB) :: eU_10_wam, eV_10_wam
      INTEGER(KIND=JWIM) :: k, IJfirst, IJlast
      REAL(KIND=JWRB) :: SumError, MaxError, TotalSum, MinValue
# ifdef DEBUG
      REAL(KIND=JWRB) :: siz, avgHS, avgWindSpeed
      logical IsFirst
      REAL(KIND=JWRB) :: minU10new, maxU10new
# endif
      gValue=9.806
      IG=1
      idx_loc=0
      IF (LLUNSTR) THEN
        IJfirst=1
        IJlast=NP_RES
      ELSE
        IJfirst=IJSLOC
        IJlast =IJLLOC
      END IF
# ifdef DEBUG
      IsFirst=.true.
# endif
      DO IJ=IJfirst,IJlast
        idx_loc=idx_loc+1
        eSingZ0=Z0NEW(IJ)
        eUstarE=USNEW(IJ)
        gzValue=gValue*eSingZ0
        eAlpha=gzValue/(eUstarE*eUstarE)
        eHS=4*SQRT(EMEAN(IJ))
        eDir=THQ(IJ)
        ePeakFreq=FPMEAN(IJ)
        eMeanFreq=FMEAN(IJ)
        eTAU = MAX(USNEW(IJ)**2,EPSUS)
        eCD = eTAU/MAX(U10NEW(IJ)**2,EPSU10)
        eWindSpeed=U10NEW(IJ)
        eWindDir=THWNEW(IJ)
# ifdef DEBUG
        IF (IsFirst) THEN
          minU10new=eWindSpeed
          maxU10new=eWindSpeed
        ELSE
          IF (minU10new .gt. eWindSpeed) THEN
            minU10new=eWindSpeed
          END IF
          IF (maxU10new .lt. eWindSpeed) THEN
            maxU10new=eWindSpeed
          END IF
        END IF
        IsFirst=.FALSE.
# endif
        eU_10_wam=eWindSpeed*SIN(eWindDir)
        eV_10_wam=eWindSpeed*COS(eWindDir)
# if defined MODEL_COUPLING_ATM_WAV && defined DEBUG
        eU_10=Uwind_atm(IJ)
        eV_10=Vwind_atm(IJ)
        diffU=abs(eU_10 - eU_10_wam)
        diffV=abs(eV_10 - eV_10_wam)
        diffUV=diffU + diffV
        IF (diffU .gt. 1) THEN
          WRITE(740+MyRankGlobal,*) 'IJ / idx_loc=', IJ, idx_loc
          WRITE(740+MyRankGlobal,*) 'U_10(wam/cosmo)=', eU_10_wam, eU_10
          WRITE(740+MyRankGlobal,*) 'V_10(wam/cosmo)=', eV_10_wam, eV_10
          FLUSH(740+MyRankGlobal)
        END IF
# endif
        eU_10=eU_10_wam
        eV_10=eV_10_wam
        eWindSpeed=SQRT(eU_10**2 + eV_10**2)
        !
        NETCDF_var(1, idx_loc)=eU_10
        NETCDF_var(2, idx_loc)=eV_10
        NETCDF_var(3, idx_loc)=eWindSpeed
        NETCDF_var(4, idx_loc)=eAlpha
        NETCDF_var(5, idx_loc)=eUstarE
        NETCDF_var(6, idx_loc)=eSingZ0
        NETCDF_var(7, idx_loc)=eHS
        NETCDF_var(8, idx_loc)=eDir
        NETCDF_var(9, idx_loc)=ePeakFreq
        NETCDF_var(10,idx_loc)=eMeanFreq
        NETCDF_var(11,idx_loc)=eCD
        idx=11
# if defined MODEL_COUPLING_OCN_WAV
        idx=idx+1
        NETCDF_var(idx,idx_loc)=Ucurr_ocn(NlevelVert, IJ)
        idx=idx+1
        NETCDF_var(idx,idx_loc)=Vcurr_ocn(NlevelVert, IJ)
        idx=idx+1
        NETCDF_var(idx,idx_loc)=ZETA_ocn(IJ)
# endif
        IF (LLUNSTR .and. LCFL) THEN
          idx=idx+1
          NETCDF_var(idx,idx_loc)=CFLCXY(1,IJ)
          idx=idx+1
          NETCDF_var(idx,idx_loc)=CFLCXY(2,IJ)
          idx=idx+1
          NETCDF_var(idx,idx_loc)=CFLCXY(3,IJ)
        END IF
      END DO
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'U10NEW(min/max)=', minval(U10NEW), maxval(U10NEW)
      WRITE(740+MyRankGlobal,*) 'U10NEW(min/max)sel=', minU10new, maxU10new
      WRITE(740+MyRankGlobal,*) 'From NETCDF_var'
      WRITE(740+MyRankGlobal,*) 'U_10(min/max)=', minval(NETCDF_var(1,:)), maxval(NETCDF_var(1,:))
      WRITE(740+MyRankGlobal,*) 'V_10(min/max)=', minval(NETCDF_var(2,:)), maxval(NETCDF_var(2,:))
#  if defined MODEL_COUPLING_ATM_WAV
      WRITE(740+MyRankGlobal,*) 'Uwind_atm(min/max)=', minval(Uwind_atm), maxval(Uwind_atm)
      WRITE(740+MyRankGlobal,*) 'Vwind_atm(min/max)=', minval(Vwind_atm), maxval(Vwind_atm)
#  endif
      WRITE(740+MyRankGlobal,*) 'MaxWind=', maxval(NETCDF_var(3,:))
      WRITE(740+MyRankGlobal,*) 'MinWind=', minval(NETCDF_var(3,:))
#  if defined MODEL_COUPLING_OCN_WAV
      WRITE(740+MyRankGlobal,*) 'Usurf(min/max)=', minval(NETCDF_var(12,:)), maxval(NETCDF_var(12,:))
      WRITE(740+MyRankGlobal,*) 'Vsurf(min/max)=', minval(NETCDF_var(13,:)), maxval(NETCDF_var(13,:))
      DO k=1,NlevelVert
        WRITE(740+MyRankGlobal,*) 'k=', k
        WRITE(740+MyRankGlobal,*) 'Ucurr_ocn(min/max)=', minval(Ucurr_ocn(k,IJSLOC:IJLLOC)), maxval(Ucurr_ocn(k,IJSLOC:IJLLOC))
        WRITE(740+MyRankGlobal,*) 'Vcurr_ocn(min/max)=', minval(Vcurr_ocn(k,IJSLOC:IJLLOC)), maxval(Vcurr_ocn(k,IJSLOC:IJLLOC))
      END DO
#  endif
      siz=REAL(1 + IJLLOC - IJSLOC)
      avgHS=sum(NETCDF_var(7,:))/siz
      avgWindSpeed=sum(NETCDF_var(3,:))/siz
      WRITE(740+MyRankGlobal,*) 'avg(HS/U10)=', avgHS, avgWindSpeed
      FLUSH(740+MyRankGlobal)
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GET_DEPTHG_tot(DEPTHG_tot)
      USE YOWPARAM , ONLY : NGX      ,NGY, NIBLO
      USE YOWSPEC, ONLY   : NSTART   ,NEND
      USE YOWMPP   , ONLY : NPROC, IRANK
      USE YOWGRID, ONLY : IJS      ,IJL
      USE YOWSHAL  , ONLY : DEPTH
      USE yownodepool,     ONLY : nodes_global
      USE YOWUNPOOL, ONLY : LLUNSTR
      USE yowpd    , only : np_global, MNP=>npa
      IMPLICIT NONE
      INTEGER(KIND=JWIM) :: i, j, IX
      INTEGER(KIND=JWIM) :: IRANKdepth, IG, IJ
      REAL(KIND=JWRB), intent(inout) :: DEPTHG_tot(NETCDF_X, NETCDF_Y)
      REAL(KIND=JWRB) :: TEMP(NIBLO)
      logical LSQRT
      REAL(KIND=JWRU) :: eREAL_8
      REAL(KIND=JWRB) :: eREAL_4
      REAL(KIND=JWRB) :: ZMISS
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'Begin of function GET_DEPTHG_tot'
      WRITE(740+MyRankGlobal,*) 'LLUNSTR=', LLUNSTR
      WRITE(740+MyRankGlobal,*) 'NETCDF_X=', NETCDF_X
      WRITE(740+MyRankGlobal,*) 'NETCDF_Y=', NETCDF_Y
      WRITE(740+MyRankGlobal,*) 'size(DEPTHG_tot,1)=', size(DEPTHG_tot,1)
      WRITE(740+MyRankGlobal,*) 'size(DEPTHG_tot,2)=', size(DEPTHG_tot,2)
      WRITE(740+MyRankGlobal,*) 'allocated(nodes_global)=', allocated(nodes_global)
      WRITE(740+MyRankGlobal,*) 'size(nodes_global)=', size(nodes_global)
      FLUSH(740+MyRankGlobal)
# endif
      IF (LLUNSTR) THEN
# ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'Before IX loop'
        WRITE(740+MyRankGlobal,*) 'np_global=', np_global
        FLUSH(740+MyRankGlobal)
# endif
        DO IX=1,np_global
          eREAL_8=nodes_global(IX)%z
          eREAL_4=SNGL(eREAL_8)
          DEPTHG_tot(IX,1)=eREAL_4
        END DO
# ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'After IX loop'
        FLUSH(740+MyRankGlobal)
# endif
      ELSE
        ZMISS=-999
        DEPTHG_tot=0
        IRANKdepth=1
        IG=1
          DO IJ = IJS(IG),IJL(IG)
            TEMP(IJ) = DEPTH(IJ,IG)
          END DO
          LSQRT=.FALSE.
          CALL MPGATHERSCFLD(IRANKdepth,NSTART, NEND,TEMP,NIBLO)
          IF(IRANKdepth.EQ.IRANK) THEN
            CALL MAKEGRID (TEMP,DEPTHG_tot,IG,LSQRT,ZMISS)
          END IF
      END IF
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'End of function GET_DEPTHG_tot'
      FLUSH(740+MyRankGlobal)
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAV_DEFINE_VAR(ncid, nx_wav_dims, ny_wav_dims, ntime_dims, eVarName, eVarUnit, LongName)
      USE netcdf
      USE YOWUNPOOL ,ONLY : LLUNSTR
      IMPLICIT NONE
      REAL(KIND=JWRB) :: TheFillVal
      INTEGER(KIND=JWIM) :: var_id, iret
      INTEGER(KIND=JWIM), intent(in) :: ncid, nx_wav_dims, ny_wav_dims, ntime_dims
      character (len = *), intent(in) :: eVarName
      character (len = *), intent(in) :: eVarUnit
      character (len = *), intent(in) :: LongName
      !
      character (len = *), parameter :: UNITS = "units"
      character (len = *), parameter :: LONG_NAME = "longname"
      character (len = *), parameter :: FILLVAL = "_FillValue"
      character (len = *), parameter :: CallFct="WAV_create_netcdf"
      TheFillVal=-999
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*)  'WAV_DEFINE_VAR, begin'
      WRITE(740+MyRankGlobal,*)  'LLUNSTR=', LLUNSTR
      WRITE(740+MyRankGlobal,*)  'eVarUnit=', TRIM(eVarUnit)
      WRITE(740+MyRankGlobal,*)  'eVarName=', TRIM(eVarName)
      WRITE(740+MyRankGlobal,*)  'LongName=', TRIM(LongName)
      FLUSH(740+MyRankGlobal)
# endif
      IF (LLUNSTR) THEN
        iret=nf90_def_var(ncid,TRIM(eVarName),NF90_REAL,(/ nx_wav_dims, ntime_dims /),var_id)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 1, iret)
      ELSE
        iret=nf90_def_var(ncid,TRIM(eVarName),NF90_REAL,(/ nx_wav_dims, ny_wav_dims, ntime_dims /),var_id)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 1, iret)
      END IF
      iret=nf90_put_att(ncid,var_id,UNITS,TRIM(eVarUnit))
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 2, iret)
      iret=nf90_put_att(ncid,var_id,LONG_NAME,TRIM(LongName))
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 3, iret)
      iret=nf90_put_att(ncid,var_id,FILLVAL,TheFillVal)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 4, iret)
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*)  'WAV_DEFINE_VAR, end'
      FLUSH(740+MyRankGlobal)
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAV_create_netcdf(idxFile, DEPTHG_tot)
      USE netcdf
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NBLO     ,NIBLO    , &
     &            CLDOMAIN, LWDINTS, NFRE, NANG
      USE YOWFRED, ONLY : FR, TH
      USE YOWSPEC, ONLY   : NSTART   ,NEND
      USE YOWMAP, ONLY : AMOWEP, AMONOP, XDELLA, XDELLO
      USE YOWGRID, ONLY : IJS      ,IJL
      USE YOWMPP   , ONLY : NPROC, IRANK
      USE YOWPD,     ONLY : nodes_global
      USE yowelementpool, only : ne_global, INE_global
      USE YOWUNPOOL, ONLY : LLUNSTR, LCFL
# if defined MODEL_COUPLING_OCN_WAV
      USE pgmcl_lib_WAM, only : MSK_att_ocn_uv, MSK_att_ocn_rho
# endif
      implicit none
      INTEGER(KIND=JWIM), intent(in) :: idxFile
      REAL(KIND=JWRB), intent(in) :: DEPTHG_tot(NETCDF_X,NETCDF_Y)
      character (len = 400) :: FILE_NAME
      character(len=4) eStr
      character (len = *), parameter :: CallFct="WAV_create_netcdf"
      character (len = *), parameter :: UNITS = "units"
      character (len = *), parameter :: LONG_NAME = "longname"
      character (len = *), parameter :: FILLVAL = "_FillValue"
      INTEGER(KIND=JWIM) :: iret, ncid, ntime_dims, fifteen_dims
      INTEGER(KIND=JWIM) :: mne_dims, three_dims
      INTEGER(KIND=JWIM) :: nx_wav_dims, ny_wav_dims
      INTEGER(KIND=JWIM) :: irec_dim, var_id, nfre_dims, nang_dims, one_dims
      REAL(KIND=JWRB) :: eFieldOne(1)
      INTEGER(KIND=JWIM) :: eInt(1)
      REAL(KIND=JWRB) :: TheFillVal
      REAL(KIND=JWRB) :: LON_wav_loc(NETCDF_X,NETCDF_Y), LAT_wav_loc(NETCDF_X,NETCDF_Y), MSK_wav_loc(NETCDF_X,NETCDF_Y)
      INTEGER(KIND=JWIM) :: i, j, IX
      REAL(KIND=JWRB) :: ZMISS
      ZMISS=-999
      TheFillVal=-999
      !
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'Beginning of WAV_create_netcdf'
      FLUSH(740+MyRankGlobal)
# endif
      CALL WAV_GetString(idxFile, eStr)
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) '1: eStr=', eStr
      FLUSH(740+MyRankGlobal)
# endif
      IF (LLUNSTR) THEN
        DO IX=1,NETCDF_X
          LON_wav_loc(IX,1)=nodes_global(IX) % x
          LAT_wav_loc(IX,1)=nodes_global(IX) % y
          MSK_wav_loc(IX,1)=1
        END DO
      ELSE
        DO i=1,NGX
          DO j=1,NGY
            LON_wav_loc(i,j)=AMOWEP+(i-1)*XDELLO
            LAT_wav_loc(i,j)=AMONOP-(j-1)*XDELLA
            IF (DEPTHG_tot(i,j).eq.ZMISS) THEN
              MSK_wav_loc(i,j)=0
            ELSE
              MSK_wav_loc(i,j)=1
            END IF
          END DO
        END DO
      END IF
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'xSize_WAV=', NETCDF_X
      WRITE(740+MyRankGlobal,*) 'ySize_WAV=', NETCDF_Y
      WRITE(740+MyRankGlobal,*) 'NFRE=', NFRE
      WRITE(740+MyRankGlobal,*) 'NANG=', NANG
      FLUSH(740+MyRankGlobal)
# endif
      FILE_NAME= "WAM_output_" // eStr // ".nc"
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'FILE_NAME=', TRIM(FILE_NAME)
      FLUSH(740+MyRankGlobal)
# endif
      iret = nf90_create(TRIM(FILE_NAME), NF90_CLOBBER, ncid)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 1, iret)
      !
      ! Time related definitions
      !
      CALL CPL_WRITE_NETCDF_TIME_HEADER(ncid, ntime_dims)
      !
      ! Definition of dimensions
      !
      iret = nf90_def_dim(ncid, 'nx_wav', NETCDF_X, nx_wav_dims)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 5, iret)
      IF (LLUNSTR) THEN
        iret = nf90_def_dim(ncid, 'ne_global', ne_global, mne_dims)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 6, iret)
        iret = nf90_def_dim(ncid, 'three', 3, three_dims)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 6, iret)
      END IF
      IF (.NOT. LLUNSTR) THEN
        iret = nf90_def_dim(ncid, 'ny_wav', NETCDF_Y, ny_wav_dims)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 6, iret)
      END IF
      iret = nf90_def_dim(ncid, 'nfre', NFRE, nfre_dims)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 7, iret)
      iret = nf90_def_dim(ncid, 'one', 1, one_dims)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 8, iret)
      iret = nf90_def_dim(ncid, 'nang', NANG, nang_dims)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 9, iret)
      !
      ! Definition of var idxoutput
      !
      iret=nf90_def_var(ncid,'idxOutput',NF90_INT,(/ ntime_dims/),var_id)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 11, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'nondimensional')
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 12, iret)
      iret=nf90_put_att(ncid,var_id,LONG_NAME,'model time step')
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 13, iret)
      !
      ! Definition of LLUNSTR
      !
      iret=nf90_def_var(ncid,'LLUNSTR',NF90_INT,(/ one_dims/),var_id)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 14, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'nondimensional')
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 15, iret)
      iret=nf90_put_att(ncid,var_id,LONG_NAME,'1 for unstructured, 0 otherwise')
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 16, iret)
      !
      ! Definition of var idxfile
      !
      iret=nf90_def_var(ncid,'idxFile',NF90_INT,(/ ntime_dims/),var_id)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 14, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'nondimensional')
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 15, iret)
      iret=nf90_put_att(ncid,var_id,LONG_NAME,'index of file')
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 16, iret)
      !
      ! Definition of var idxInfile
      !
      iret=nf90_def_var(ncid,'idxInfile',NF90_INT,(/ ntime_dims/),var_id)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 17, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'nondimensional')
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 18, iret)
      iret=nf90_put_att(ncid,var_id,LONG_NAME,'index in the file')
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 19, iret)
      IF (idxFile .eq. 1) THEN
        IF (LLUNSTR) THEN
          !
          ! Definition of var LON_wav
          !
          iret=nf90_def_var(ncid,'LON_wav',NF90_REAL,(/ nx_wav_dims/),var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 38, iret)
          iret=nf90_put_att(ncid,var_id,UNITS,'deg')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 39, iret)
          iret=nf90_put_att(ncid,var_id,LONG_NAME,'longitude of wave model')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 40, iret)
          !
          ! Definition of var LAT_wav
          !
          iret=nf90_def_var(ncid,'LAT_wav',NF90_REAL,(/ nx_wav_dims/),var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 41, iret)
          iret=nf90_put_att(ncid,var_id,UNITS,'deg')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 42, iret)
          iret=nf90_put_att(ncid,var_id,LONG_NAME,'latitude of wave model')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 43, iret)
          !
          ! Definition of var MSK_wav
          !
          iret=nf90_def_var(ncid,'MSK_wav',NF90_INT,(/ nx_wav_dims/),var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 44, iret)
          iret=nf90_put_att(ncid,var_id,UNITS,'nondimensional')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 45, iret)
          iret=nf90_put_att(ncid,var_id,LONG_NAME,'mask of wave model')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 46, iret)
#if defined MODEL_COUPLING_OCN_WAV
          !
          ! Definition of var MSK_att_ocn_uv
          !
          iret=nf90_def_var(ncid,'MSK_att_ocn_uv',NF90_INT,(/ nx_wav_dims/),var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 44, iret)
          iret=nf90_put_att(ncid,var_id,UNITS,'nondimensional')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 45, iret)
          iret=nf90_put_att(ncid,var_id,LONG_NAME,'mask of oceanic model')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 46, iret)
          !
          ! Definition of mask of oceanic model rho
          !
          iret=nf90_def_var(ncid,'MSK_att_ocn_rho',NF90_INT,(/ nx_wav_dims/),var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 44, iret)
          iret=nf90_put_att(ncid,var_id,UNITS,'nondimensional')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 45, iret)
          iret=nf90_put_att(ncid,var_id,LONG_NAME,'mask of oceanic model')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 46, iret)
#endif
          !
          ! Definition of var DEP_wav
          !
          iret=nf90_def_var(ncid,'DEP_wav',NF90_REAL,(/ nx_wav_dims/),var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 47, iret)
          iret=nf90_put_att(ncid,var_id,UNITS,'deg')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 48, iret)
          iret=nf90_put_att(ncid,var_id,LONG_NAME,'depth of wave model')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 49, iret)
          !
          ! Definition of elements
          !
          iret=nf90_def_var(ncid,'ele',NF90_INT,(/three_dims, mne_dims/),var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 27, iret)
          iret=nf90_put_att(ncid,var_id,UNITS,'non-dimensional')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 28, iret)
          iret=nf90_put_att(ncid,var_id,LONG_NAME,'connecting table elements')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 49, iret)
        ELSE
          !
          ! Definition of var LON_wav
          !
          iret=nf90_def_var(ncid,'LON_wav',NF90_REAL,(/ nx_wav_dims, ny_wav_dims/),var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 38, iret)
          iret=nf90_put_att(ncid,var_id,UNITS,'deg')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 39, iret)
          iret=nf90_put_att(ncid,var_id,LONG_NAME,'longitude of wave model')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 40, iret)
          !
          ! Definition of var LAT_wav
          !
          iret=nf90_def_var(ncid,'LAT_wav',NF90_REAL,(/ nx_wav_dims, ny_wav_dims/),var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 41, iret)
          iret=nf90_put_att(ncid,var_id,UNITS,'deg')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 42, iret)
          iret=nf90_put_att(ncid,var_id,LONG_NAME,'latitude of wave model')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 43, iret)
          !
          ! Definition of var MSK_wav
          !
          iret=nf90_def_var(ncid,'MSK_wav',NF90_INT,(/ nx_wav_dims, ny_wav_dims/),var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 44, iret)
          iret=nf90_put_att(ncid,var_id,UNITS,'nondimensional')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 45, iret)
          iret=nf90_put_att(ncid,var_id,LONG_NAME,'mask of wave model')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 46, iret)
#if defined MODEL_COUPLING_OCN_WAV
          !
          ! Definition of var MSK_att_ocn_uv
          !
          iret=nf90_def_var(ncid,'MSK_att_ocn_uv',NF90_INT,(/ nx_wav_dims, ny_wav_dims/),var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 44, iret)
          iret=nf90_put_att(ncid,var_id,UNITS,'nondimensional')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 45, iret)
          iret=nf90_put_att(ncid,var_id,LONG_NAME,'mask of oceanic model')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 46, iret)
          !
          ! Definition of mask of oceanic model rho
          !
          iret=nf90_def_var(ncid,'MSK_att_ocn_rho',NF90_INT,(/ nx_wav_dims, ny_wav_dims/),var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 44, iret)
          iret=nf90_put_att(ncid,var_id,UNITS,'nondimensional')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 45, iret)
          iret=nf90_put_att(ncid,var_id,LONG_NAME,'mask of oceanic model')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 46, iret)
#endif
          !
          ! Definition of var DEP_wav
          !
          iret=nf90_def_var(ncid,'DEP_wav',NF90_REAL,(/ nx_wav_dims, ny_wav_dims/),var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 47, iret)
          iret=nf90_put_att(ncid,var_id,UNITS,'deg')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 48, iret)
          iret=nf90_put_att(ncid,var_id,LONG_NAME,'depth of wave model')
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 49, iret)
        END IF
        !
        !
        ! Definition of var FR
        !
        iret=nf90_def_var(ncid,'FR',NF90_REAL,(/ nfre_dims/),var_id)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 50, iret)
        iret=nf90_put_att(ncid,var_id,UNITS,'Hz')
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 51, iret)
        !
        ! Definition of var TH
        !
        iret=nf90_def_var(ncid,'TH',NF90_REAL,(/ nang_dims/),var_id)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 52, iret)
        iret=nf90_put_att(ncid,var_id,UNITS,'rad')
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 53, iret)
        !
        iret=nf90_def_var(ncid,'AMOWEP',NF90_REAL,(/ one_dims/),var_id)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 54, iret)
        iret=nf90_put_att(ncid,var_id,UNITS,'degree')
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 55, iret)
        iret=nf90_def_var(ncid,'XDELLO',NF90_REAL,(/ one_dims/),var_id)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 56, iret)
        iret=nf90_put_att(ncid,var_id,UNITS,'degree')
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 57, iret)
        iret=nf90_def_var(ncid,'AMONOP',NF90_REAL,(/ one_dims/),var_id)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 58, iret)
        iret=nf90_put_att(ncid,var_id,UNITS,'degree')
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 59, iret)
        iret=nf90_def_var(ncid,'XDELLA',NF90_REAL,(/ one_dims/),var_id)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 60, iret)
        iret=nf90_put_att(ncid,var_id,UNITS,'degree')
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 61, iret)
      END IF
      !
      CALL WAV_DEFINE_VAR(ncid, nx_wav_dims, ny_wav_dims, ntime_dims, 'Z0wave', 'm', 'roughness length')
      CALL WAV_DEFINE_VAR(ncid, nx_wav_dims, ny_wav_dims, ntime_dims, 'Ustar', 'm/s', 'friction velocity')
      CALL WAV_DEFINE_VAR(ncid, nx_wav_dims, ny_wav_dims, ntime_dims, 'Cd', 'nondim.', 'drag coefficient')
      CALL WAV_DEFINE_VAR(ncid, nx_wav_dims, ny_wav_dims, ntime_dims, 'Hwave', 'm', 'significant wave height')
      CALL WAV_DEFINE_VAR(ncid, nx_wav_dims, ny_wav_dims, ntime_dims, 'Dwave', 'deg', 'mean wave direction')
      CALL WAV_DEFINE_VAR(ncid, nx_wav_dims, ny_wav_dims, ntime_dims, 'PwaveFreq', 'Hz', 'peak wave frequency')
      CALL WAV_DEFINE_VAR(ncid, nx_wav_dims, ny_wav_dims, ntime_dims, 'MwaveFreq', 'Hz', 'mean wave frequency')
      CALL WAV_DEFINE_VAR(ncid, nx_wav_dims, ny_wav_dims, ntime_dims, 'AlphaWave', 'nondim.', 'Charnock coefficient')
      CALL WAV_DEFINE_VAR(ncid, nx_wav_dims, ny_wav_dims, ntime_dims, 'WndMag', 'm/s', 'wind velocity')
      CALL WAV_DEFINE_VAR(ncid, nx_wav_dims, ny_wav_dims, ntime_dims, 'U_10', 'm/s', 'U wind')
      CALL WAV_DEFINE_VAR(ncid, nx_wav_dims, ny_wav_dims, ntime_dims, 'V_10', 'm/s', 'V wind')
# if defined MODEL_COUPLING_OCN_WAV
      CALL WAV_DEFINE_VAR(ncid, nx_wav_dims, ny_wav_dims, ntime_dims, 'ucurr', 'm/s', 'U current')
      CALL WAV_DEFINE_VAR(ncid, nx_wav_dims, ny_wav_dims, ntime_dims, 'vcurr', 'm/s', 'V current')
      CALL WAV_DEFINE_VAR(ncid, nx_wav_dims, ny_wav_dims, ntime_dims, 'ZetaOcean', 'm', 'free surface')
# endif
      IF (LLUNSTR .and. LCFL) THEN
        CALL WAV_DEFINE_VAR(ncid, nx_wav_dims, ny_wav_dims, ntime_dims, 'cfl1', 'nondim.', 'CFL number 1')
        CALL WAV_DEFINE_VAR(ncid, nx_wav_dims, ny_wav_dims, ntime_dims, 'cfl2', 'nondim.', 'CFL number 2')
        CALL WAV_DEFINE_VAR(ncid, nx_wav_dims, ny_wav_dims, ntime_dims, 'cfl3', 'nondim.', 'CFL number 3')
      END IF
      iret=nf90_close(ncid)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 90, iret)
      !
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'FILE_NAME=', TRIM(FILE_NAME)
      FLUSH(740+MyRankGlobal)
# endif
      iret=nf90_open(TRIM(FILE_NAME), nf90_write, ncid)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 91, iret)
      IF (idxFile .eq. 1) THEN
        IF (LLUNSTR) THEN
          eInt(1)=1
        ELSE
          eInt(1)=0
        END IF
        iret=nf90_inq_varid(ncid, "LLUNSTR", var_id)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 98, iret)
        iret=nf90_put_var(ncid,var_id,eInt,start = (/1/), count=(/ 1/))
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 99, iret)
        IF (LLUNSTR) THEN
          iret=nf90_inq_varid(ncid, "LON_wav", var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 98, iret)
          iret=nf90_put_var(ncid,var_id,LON_wav_loc(:,1),start = (/1/), count=(/ NETCDF_X/))
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 99, iret)
          !
          iret=nf90_inq_varid(ncid, "LAT_wav", var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 100, iret)
          iret=nf90_put_var(ncid,var_id,LAT_wav_loc(:,1),start = (/1/), count=(/ NETCDF_X/))
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 101, iret)
          !
          iret=nf90_inq_varid(ncid, "MSK_wav", var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 102, iret)
          iret=nf90_put_var(ncid,var_id,MSK_wav_loc(:,1),start = (/1/), count=(/ NETCDF_X/))
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 103, iret)
# if defined MODEL_COUPLING_OCN_WAV
          iret=nf90_inq_varid(ncid, "MSK_att_ocn_uv", var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 102, iret)
          iret=nf90_put_var(ncid,var_id,MSK_att_ocn_uv(:,1),start = (/1/), count=(/ NETCDF_X/))
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 103, iret)
          iret=nf90_inq_varid(ncid, "MSK_att_ocn_rho", var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 102, iret)
          iret=nf90_put_var(ncid,var_id,MSK_att_ocn_rho(:,1),start = (/1/), count=(/ NETCDF_X/))
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 103, iret)
# endif

          iret=nf90_inq_varid(ncid, "ele", var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 7, iret)
          iret=nf90_put_var(ncid,var_id,INE_global, start = (/1,1/), count = (/ 3, ne_global/))
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 8, iret)
        ELSE
          iret=nf90_inq_varid(ncid, "LON_wav", var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 98, iret)
          iret=nf90_put_var(ncid,var_id,LON_wav_loc,start = (/1, 1/), count=(/ NETCDF_X, NETCDF_Y/))
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 99, iret)
          iret=nf90_inq_varid(ncid, "LAT_wav", var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 100, iret)
          iret=nf90_put_var(ncid,var_id,LAT_wav_loc,start = (/1, 1/), count=(/ NETCDF_X, NETCDF_Y/))
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 101, iret)
          iret=nf90_inq_varid(ncid, "MSK_wav", var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 102, iret)
          iret=nf90_put_var(ncid,var_id,MSK_wav_loc,start = (/1, 1/), count=(/ NETCDF_X, NETCDF_Y/))
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 103, iret)
# if defined MODEL_COUPLING_OCN_WAV
          iret=nf90_inq_varid(ncid, "MSK_att_ocn_uv", var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 102, iret)
          iret=nf90_put_var(ncid,var_id,MSK_att_ocn_uv,start = (/1, 1/), count=(/ NETCDF_X, NETCDF_Y/))
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 103, iret)
          iret=nf90_inq_varid(ncid, "MSK_att_ocn_rho", var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 102, iret)
          iret=nf90_put_var(ncid,var_id,MSK_att_ocn_rho,start = (/1, 1/), count=(/ NETCDF_X, NETCDF_Y/))
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 103, iret)
# endif
        END IF
        iret=nf90_inq_varid(ncid, "FR", var_id)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 104, iret)
        iret=nf90_put_var(ncid,var_id,FR,start = (/1/), count=(/ NFRE/))
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 105, iret)
        iret=nf90_inq_varid(ncid, "TH", var_id)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 106, iret)
        iret=nf90_put_var(ncid,var_id,TH,start = (/1/), count=(/ NANG/))
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 107, iret)
        iret=nf90_inq_varid(ncid, "DEP_wav", var_id)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 108, iret)
        iret=nf90_put_var(ncid,var_id,DEPTHG_tot,start = (/1, 1/), count=(/ NETCDF_X, NETCDF_Y/))
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 109, iret)
        !
        iret=nf90_inq_varid(ncid, "AMOWEP", var_id)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 110, iret)
        eFieldOne(1)=AMOWEP
        iret=nf90_put_var(ncid,var_id,eFieldOne,start = (/1/), count=(/ 1/))
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 111, iret)
        iret=nf90_inq_varid(ncid, "XDELLO", var_id)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 112, iret)
        eFieldOne(1)=XDELLO
        iret=nf90_put_var(ncid,var_id,eFieldOne,start = (/1/), count=(/ 1/))
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 113, iret)
        iret=nf90_inq_varid(ncid, "AMONOP", var_id)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 114, iret)
        eFieldOne(1)=AMONOP
        iret=nf90_put_var(ncid,var_id,eFieldOne,start = (/1/), count=(/ 1/))
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 115, iret)
        iret=nf90_inq_varid(ncid, "XDELLA", var_id)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 116, iret)
        eFieldOne(1)=XDELLA
        iret=nf90_put_var(ncid,var_id,eFieldOne,start = (/1/), count=(/ 1/))
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 117, iret)
      END IF
      !
      iret=nf90_close(ncid)
      CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 118, iret)
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'WAV: preamble of netcdf file created'
      FLUSH(740+MyRankGlobal)
# endif
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'Ending of WAV_create_netcdf'
      FLUSH(740+MyRankGlobal)
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAV_WRITE_VAR(ncid, idx, idVar, VarName)
      USE netcdf
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWUNPOOL ,ONLY : LLUNSTR
      IMPLICIT NONE
      INTEGER(KIND=JWIM), intent(in) :: ncid, idx, idVar
      character(len = *), intent(in) :: VarName
      INTEGER(KIND=JWIM) :: var_id, iret
      REAL(KIND=JWRB) :: TheFieldNC(NETCDF_X,NETCDF_Y)
      character (len = *), parameter :: CallFct="WAV_create_netcdf"
      REAL(KIND=JWRB) :: TheFieldUnstruct(NETCDF_X)
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*)  'WAV_WRITE_VAR, begin'
      WRITE(740+MyRankGlobal,*)  'LLUNSTR=', LLUNSTR
      WRITE(740+MyRankGlobal,*)  'idVar=', idVar
      WRITE(740+MyRankGlobal,*)  'VarName=', TRIM(VarName)
      FLUSH(740+MyRankGlobal)
# endif
      IF (MyRankLocal.eq.0) THEN
        IF (LLUNSTR) THEN
          TheFieldUnstruct=NETCDF_var_gridded(idVar,:,1)
          iret=nf90_inq_varid(ncid, TRIM(VarName), var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 1, iret)
          iret=nf90_put_var(ncid,var_id,TheFieldUnstruct,start = (/1, idx/), count=(/NETCDF_X,1/))
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 2, iret)
        ELSE
          TheFieldNC=NETCDF_var_gridded(idVar,:,:)
          iret=nf90_inq_varid(ncid, TRIM(VarName), var_id)
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 1, iret)
          iret=nf90_put_var(ncid,var_id,TheFieldNC,start = (/1, 1, idx/), count=(/NETCDF_X,NETCDF_Y,1/))
          CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 2, iret)
        END IF
      END IF
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*)  'WAV_WRITE_VAR, end'
      FLUSH(740+MyRankGlobal)
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE wav_single_write_netcdf(idxOutput, idxFile, idxInFile, Increment)
      USE netcdf
      USE YOWINTP  , ONLY : WHGTTG   ,WDIRTG   ,WPKFTG   ,WMNFTG,     &
     &            USTARG, CDG
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWCOUT  , ONLY : IPFGTBL
      USE YOWWIND  , ONLY : FIELDG_coupl, FIELDG
      USE YOWUNPOOL, ONLY : LLUNSTR, LCFL
      USE WAV_NETCDF_FCT, ONLY : WAV_GET_ETIMEDAY
      implicit none
      character (len = 400) :: FILE_NAME
      character (len=4) :: eStr
      INTEGER(KIND=JWIM), intent(in) :: idxOutput, idxFile, idxInfile
      REAL(KIND=JWRU), intent(in) :: Increment
      INTEGER(KIND=JWIM) :: iret, ncid, var_id, idx
      INTEGER(KIND=JWIM) :: eVarInt(1)
      character (len = *), parameter :: CallFct="wav_single_write_netcdf"
      INTEGER(KIND=JWIM) :: IX, IY
      REAL(KIND=JWRB) :: eU, eV, eSpeed
      INTEGER(KIND=JWIM) :: I, J
      REAL(KIND=JWRU):: eTimeDay
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*)  'allocated(FIELDG)=', allocated(FIELDG)
      WRITE(740+MyRankGlobal,*)  'allocated(FIELDG_coupl)=', allocated(FIELDG_coupl)
      FLUSH(740+MyRankGlobal)
# endif
      idx=idxInfile
      CALL WAV_GetString(idxFile, eStr)
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) '2: eStr=', eStr
      FLUSH(740+MyRankGlobal)
# endif
      FILE_NAME= "WAM_output_" // eStr // ".nc"
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'FILE_NAME=', TRIM(FILE_NAME)
      FLUSH(740+MyRankGlobal)
# endif
      IF (MyRankLocal.eq.0) THEN
        iret=nf90_open(TRIM(FILE_NAME), nf90_write, ncid)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 1, iret)
      END IF
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*)  'wav_single_write_netcdf, step 1'
      FLUSH(740+MyRankGlobal)
# endif
      IF (MyRankLocal.eq.0) THEN
        CALL WAV_GET_ETIMEDAY(eTimeDay, Increment)
        CALL CPL_WRITE_NETCDF_TIME(ncid, idx, eTimeDay)
        !
        iret=nf90_inq_varid(ncid, "idxOutput", var_id)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 2, iret)
        eVarInt(1)=idxOutput
        iret=nf90_put_var(ncid,var_id,eVarInt,start = (/idx/), count=(/ 1 /))
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 3, iret)
        iret=nf90_inq_varid(ncid, "idxFile", var_id)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 4, iret)
        eVarInt(1)=idxFile
        iret=nf90_put_var(ncid,var_id,eVarInt,start = (/idx/), count=(/ 1 /))
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 5, iret)
        iret=nf90_inq_varid(ncid, "idxInfile", var_id)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 6, iret)
        eVarInt(1)=idxInfile
        iret=nf90_put_var(ncid,var_id,eVarInt,start = (/idx/), count=(/ 1 /))
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 7, iret)
      END IF
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*)  'wav_single_write_netcdf, step 2'
      FLUSH(740+MyRankGlobal)
# endif
# ifdef MODEL_COUPLING_ATM_WAV
      CALL WAV_WRITE_VAR(ncid, idx, 4, "AlphaWave")
      CALL WAV_WRITE_VAR(ncid, idx, 5, "Ustar")
      CALL WAV_WRITE_VAR(ncid, idx, 6, "Z0wave")
# endif
      CALL WAV_WRITE_VAR(ncid, idx, 1, "U_10")
      CALL WAV_WRITE_VAR(ncid, idx, 2, "V_10")
      CALL WAV_WRITE_VAR(ncid, idx, 3, "WndMag")
      CALL WAV_WRITE_VAR(ncid, idx, 7, "Hwave")
      IF (maxval(NETCDF_VAR_GRIDDED(7, :, :)) > 100) THEN
        Print *, 'HS is larger than 100 m'
        Print *, 'We stop the program now'
        STOP
      END IF
      CALL WAV_WRITE_VAR(ncid, idx, 8, "Dwave")
      CALL WAV_WRITE_VAR(ncid, idx, 9, "PwaveFreq")
      CALL WAV_WRITE_VAR(ncid, idx, 10, "MwaveFreq")
      CALL WAV_WRITE_VAR(ncid, idx, 11, "Cd")
      !
# if defined MODEL_COUPLING_OCN_WAV
      CALL WAV_WRITE_VAR(ncid, idx, idxUcurr, "ucurr")
      CALL WAV_WRITE_VAR(ncid, idx, idxVcurr, "vcurr")
      CALL WAV_WRITE_VAR(ncid, idx, idxZeta, "ZetaOcean")
# endif
      IF (LLUNSTR .and. LCFL) THEN
        CALL WAV_WRITE_VAR(ncid, idx, idxcfl1, "cfl1")
        CALL WAV_WRITE_VAR(ncid, idx, idxcfl2, "cfl2")
        CALL WAV_WRITE_VAR(ncid, idx, idxcfl3, "cfl3")
      END IF
      IF (MyRankLocal.eq.0) THEN
        iret=nf90_close(ncid)
        CALL WAV_GENERIC_NETCDF_ERROR(CallFct, 8, iret)
      END IF
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*)  'wav_single_write_netcdf, step 3'
      FLUSH(740+MyRankGlobal)
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAV_netcdf_output
      USE YOWSTAT  , ONLY : IDELT
      USE YOWGRID  , ONLY : IJS      ,IJL, IJSLOC, IJLLOC, IJGLOBAL_OFFSET
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWUNIT  , ONLY : IU25     ,IU26
      USE YOWSPEC, ONLY   : FL3
      IMPLICIT NONE
      INTEGER(KIND=JWIM) :: IG
      LOGICAL DO_OUTPUT, TEST_CON, TEST_NETCDF
      REAL(KIND=JWRB), allocatable :: DEPTHG_tot(:,:)
      LOGICAL LLOUTBS
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'wav_netcdf_output, NETCDF_X=', NETCDF_X
      WRITE(740+MyRankGlobal,*) 'wav_netcdf_output, NETCDF_Y=', NETCDF_Y
      FLUSH(740+MyRankGlobal)
# endif
      IF (DeltaTimeNetcdfWAV .eq. -1) THEN
        DeltaTimeNetcdfWAV=DBLE(IDELT)
      END IF
      IF (FileSizeTimeWAV .eq. -1) THEN
        FileSizeTimeWAV=86400 ! daily output is the default
      END IF
      CALL TEST_DIV(FileSizeTimeWAV, DeltaTimeNetcdfWAV, TEST_NETCDF)
      IF (.NOT. TEST_NETCDF) THEN
        Print *, 'DeltaTimeNetcdfWAV=', DeltaTimeNetcdfWAV
        Print *, 'FileSizeTimeWAV=', FileSizeTimeWAV
        Print *, 'FileSizeTimeWAV should be a multiple of DeltaTimeNetcdfWAV'
        STOP
      END IF
      NDEFHIS_WAV=INT(FileSizeTimeWAV / DeltaTimeNetcdfWAV)
      CALL TEST_DIV(DeltaTimeNetcdfWAV, DBLE(IDELT), TEST_CON)
      IF (.NOT. TEST_CON) THEN
        Print *, 'DeltaTimeNetcdfWAV=', DeltaTimeNetcdfWAV
        Print *, 'IDELT=', IDELT
        Print *, 'DeltaTimeNetcdfWAV should be a multiple of IDELT'
        STOP
      END IF
      CALL TEST_DIV(WAV_NetcdfPresTime, DeltaTimeNetcdfWAV, DO_OUTPUT)
      IF (DO_OUTPUT) THEN
        IG=1
        DoNETCDF_sync=.TRUE.
# ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'minval(FL3)=', minval(FL3(IJSLOC:IJLLOC,:,:))
        FLUSH(740+MyRankGlobal)
# endif
        CALL OUTBS (FL3(IJSLOC:IJLLOC,:,:), IJSLOC, IJLLOC, IJGLOBAL_OFFSET, IG, IU25, IU26, LLOUTBS)
        DoNETCDF_sync=.FALSE.
        CALL WAV_netcdf_export
# ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'WAVw: beginning'
        WRITE(740+MyRankGlobal,*) 'WAVw: 1 idxOutput=', idxOutput, idxFile, idxInfile
        FLUSH(740+MyRankGlobal)
# endif
        IF (idxInFile .eq. 1) THEN
# ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'WAV_create_netcdf, MRL=', MyRankLocal, 'idxFile=', idxFile
          FLUSH(740+MyRankGlobal)
# endif
          allocate(DEPTHG_tot(NETCDF_X, NETCDF_Y))
# ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'Before GET_DEPTHG_TOT function'
          WRITE(740+MyRankGlobal,*) 'NETCDF_X=', NETCDF_X
          WRITE(740+MyRankGlobal,*) 'NETCDF_Y=', NETCDF_Y
          WRITE(740+MyRankGlobal,*) 'size(DEPTHG_tot,1)=', size(DEPTHG_tot,1)
          WRITE(740+MyRankGlobal,*) 'size(DEPTHG_tot,2)=', size(DEPTHG_tot,2)
          FLUSH(740+MyRankGlobal)
# endif
          CALL GET_DEPTHG_tot(DEPTHG_tot)
# ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'After GET_DEPTHG_tot'
          FLUSH(740+MyRankGlobal)
# endif
          IF (MyRankLocal .eq. 0) THEN
            CALL WAV_create_netcdf(idxFile, DEPTHG_tot)
          END IF
          deallocate(DEPTHG_tot)
# ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'After NETCDF creation'
          FLUSH(740+MyRankGlobal)
# endif
        ENDIF
# ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'WAVw: WRITE TO THE NETCDF FILE, idxInfile=', idxInfile
        FLUSH(740+MyRankGlobal)
# endif
        CALL wav_single_write_netcdf(idxOutput, idxFile, idxInfile, WAV_NetcdfPresTime)
# ifdef DEBUG
        WRITE(740+MyRankGlobal,*) 'WAVw: ending'
        FLUSH(740+MyRankGlobal)
# endif
        idxInfile=idxInfile+1
        IF (idxInFile .eq. NDEFHIS_WAV + 1) THEN
          idxFile=idxFile+1
          idxInfile=1
        END IF
# ifdef DEBUG
        WRITE(740+MyRankGlobal,*)  'WAVw: 2 idxOutput=', idxOutput, idxFile, idxInfile
        FLUSH(740+MyRankGlobal)
# endif
      ENDIF
      WAV_NetcdfPresTime = WAV_NetcdfPresTime + DBLE(IDELT)
      idxOutput=idxOutput+1
# ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'End of WAV_netcdf_output'
      FLUSH(740+MyRankGlobal)
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAV_GetString (TheNb, eStr)
      implicit none
      character*4, intent(out) :: eStr
      INTEGER(KIND=JWIM), intent(in) :: TheNb
      INTEGER(KIND=JWIM) :: STAT_VALUE
      character*1 eStr1
      character*2 eStr2
      character*3 eStr3
      IF (TheNb.le.9) THEN
         WRITE (FMT=10, UNIT=eStr1, IOSTAT=STAT_VALUE) TheNb
         eStr='000' // eStr1
      ELSE IF (TheNb.le.99) THEN
         WRITE (FMT=20, UNIT=eStr2, IOSTAT=STAT_VALUE) TheNb
         eStr='00' // eStr2
      ELSE IF (TheNb.le.999) THEN
         WRITE (FMT=30, UNIT=eStr3, IOSTAT=STAT_VALUE) TheNb
         eStr='0' // eStr3
      ELSE
         WRITE (FMT=40, UNIT=eStr, IOSTAT=STAT_VALUE) TheNb
      END IF
 10   FORMAT (i1)
 20   FORMAT (i2)
 30   FORMAT (i3)
 40   FORMAT (i4)
      END SUBROUTINE
#else
  contains
     SUBROUTINE WAV_NETCDF_OUTPUT()
        write(*,*) 'calling WAV_NETCDF_OUTPUT '
        write(*,*) '!!!!!!!!!!!! this is a dummy version of it !!!'
        call abort
     END SUBROUTINE
     SUBROUTINE WAV_NETCDF_EXPORT()
        write(*,*) 'calling WAV_NETCDF_EXPORT '
        write(*,*) '!!!!!!!!!!!! this is a dummy version of it !!!'
        call abort
     END SUBROUTINE
#endif
END MODULE WAV_netcdf
