      SUBROUTINE WAM2ODB(IJS, IJL, HSOIB, HSAN, U10FG, U10AN)

! ----------------------------------------------------------------------

!****  *WAM2ODB* - WRITE WAVE DATA BACK INTO ODB


!     PURPOSE.                                                          
!     --------                                                          

!*    INTERFACE.                                                        
!     ----------                                                        

!     *CALL* *WAM2ODB

!     METHOD.                                                           
!     -------                                                           

!     EXTERNALS.                                                        
!     ----------                                                        

!     REFERENCES.                                                       
!     -----------                                                       

!          NONE                                                         

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWALTAS , ONLY : NIJALT   ,NOBSPE, IJALT ,DIFFALTFG, ALTDATA
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : LL1D
      USE YOWSTAT  , ONLY : CDTPRO   ,IDELALT
      USE YOWSPEC  , ONLY : NSTART   ,NEND          ,IJ2NEWIJ
      USE YOWTEST  , ONLY : IU06

      USE MPL_MODULE

#ifdef WITH_ODB
      USE YOWODB   , ONLY : WAM_YDODB
      USE DBASE_VIEW_MOD, ONLY: DBASE_VIEW
#endif

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "incdate.intfb.h"
#include "mpgatherscfld.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: HSOIB, HSAN, U10FG, U10AN

#ifdef WITH_ODB
      TYPE(DBASE_VIEW) :: ROBSU

      INTEGER(KIND=JWIM) :: IJ, IR, IRN, IRK, IC
      INTEGER(KIND=JWIM) :: ICOUNT
      INTEGER(KIND=JWIM) :: IRECV
      INTEGER(KIND=JWIM) :: I_PARAM, NUM_PARAM
      INTEGER(KIND=JWIM) :: IZCOMLEN, NENTRY, NOBSPETOT
      INTEGER(KIND=JWIM) :: JOBS, IOBS
      INTEGER(KIND=JWIM) :: IDATE, ITIME
      INTEGER(KIND=JWIM) :: IJALT_MDL, IJALT_ODB, ISAT_MDL, ISAT_ODB
      INTEGER(KIND=JWIM),DIMENSION(NPROC,2) :: NPEOBS
      INTEGER(KIND=JWIM) :: IRET
      INTEGER(KIND=JWIM) :: IJALT_LOC(NIJALT)
      INTEGER(KIND=JWIM) :: IPREV_SEQNO
      INTEGER(KIND=JWIM) :: I_START_OF_THIS_REPORT

      REAL(KIND=JWRB), PARAMETER :: ODB_RMDI = -2147483647.0_JWRB  ! Missing value in ODB
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: ZFG_DEPART, ZAN_DEPART
      REAL(KIND=JWRB) :: ZU10FG, ZU10AN
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:) :: ZCOMBUF
      REAL(KIND=JWRU), POINTER :: ZODB_SEQNO(:) => NULL()
      REAL(KIND=JWRU), POINTER :: ZODB_DATE(:) => NULL()
      REAL(KIND=JWRU), POINTER :: ZODB_TIME(:) => NULL()
      REAL(KIND=JWRU), POINTER :: ZODB_VARNO(:) => NULL()
      REAL(KIND=JWRU), POINTER :: ZODB_OBSVALUE(:) => NULL()
      REAL(KIND=JWRU), POINTER :: ZODB_DATUM_STATUS(:) => NULL()
      REAL(KIND=JWRU), POINTER :: ZODB_REPORT_STATUS(:) => NULL()
      REAL(KIND=JWRU), POINTER :: ZODB_DATUM_EVENT2(:) => NULL()
      REAL(KIND=JWRU), POINTER :: ZODB_GP_NUMBER(:) => NULL()
      REAL(KIND=JWRU), POINTER :: ZODB_SATELLITE_IDENTIFIER(:) => NULL()
      REAL(KIND=JWRU), POINTER :: ZODB_FG_DEPAR(:) => NULL()
      REAL(KIND=JWRU), POINTER :: ZODB_AN_DEPAR(:) => NULL()

      REAL(KIND=JWRU) :: ZINVAR, ZSWITCH_BIT
      INTEGER(KIND=JWIM) :: INVAL, INBIT

      ZSWITCH_BIT(ZINVAR,INVAL,INBIT) = IOR(INT(ZINVAR),ISHFT(INVAL,INBIT))
      CHARACTER(LEN=14) :: CBEGINDT, CENDDT, COBSDT

      LOGICAL :: LLEPSMIN
      LOGICAL(KIND=1) :: LLIN 
      LOGICAL(KIND=1), ALLOCATABLE, DIMENSION(:) :: LUNSET_FLAGS
      LOGICAL(KIND=1), ALLOCATABLE, DIMENSION(:) :: LIN_TIMWIN
      LOGICAL(KIND=1) :: L_ANY_ACTIVE_THIS_REPORT

!     -------------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('WAM2ODB',0,ZHOOK_HANDLE)

      NUM_PARAM=2  ! Significant wave height and Surface wind speed
      LLEPSMIN=.TRUE.

      IRECV=1

!     OUTPUT TO ODB:
!     -------------
      ! Define the start and end of time window
      CBEGINDT = CDTPRO
      IF (IDELALT <= 10800) THEN
        CALL INCDATE (CBEGINDT,-IDELALT)
      ELSE 
        CALL INCDATE (CBEGINDT,-IDELALT/2)
      ENDIF
      CENDDT = CBEGINDT
      CALL INCDATE (CENDDT,IDELALT)

      NOBSPETOT=SUM(NOBSPE)

      NENTRY=NIJALT+4
      IZCOMLEN = NPROC+NENTRY*NOBSPETOT
      ALLOCATE(ZCOMBUF(IZCOMLEN))
      ! ZCOMBUF needs to be initialised to 0
      ZCOMBUF(:)=0.0_JWRB

      IR=1
      NPEOBS(IR,1)=1
      NPEOBS(IR,2)=NPEOBS(IR,1)+1+NENTRY*NOBSPE(IR)-1
      DO IR = 2,NPROC
        NPEOBS(IR,1)=NPEOBS(IR-1,2)+1
        NPEOBS(IR,2)=NPEOBS(IR,1)+1+NENTRY*NOBSPE(IR)-1
      ENDDO

      ICOUNT = NPEOBS(IRANK,1)
      ZCOMBUF(ICOUNT) = IRANK 
      DO IOBS=1,NOBSPE(IRANK)
        IJALT_MDL=IJALT(IOBS,1)

        IF (IJALT_MDL >= NSTART(IRANK) .AND. IJALT_MDL <= NEND(IRANK)) THEN
          ! Observations overlap, so take only those belonging to IRANK
          LLIN=.TRUE.
        ELSE
          LLIN=.FALSE.
        ENDIF

        DO IC=1,NIJALT
          ICOUNT = ICOUNT + 1
          IF (LLIN) ZCOMBUF(ICOUNT) = FLOAT(IJALT(IOBS,IC))
        ENDDO

          ! SWH FG DEPAR, OVER-WRITING DIFFALTFG(IOBS)
          ICOUNT = ICOUNT + 1
          ZCOMBUF(ICOUNT) = DIFFALTFG(IOBS) 

          ! SWH AN DEPAR
          ICOUNT = ICOUNT + 1
          IF (LLIN) THEN
            IF (HSOIB(IJALT_MDL) <= 0.0_JWRB) THEN
              ZCOMBUF(ICOUNT)=ODB_RMDI
            ELSE
              ZCOMBUF(ICOUNT)=ALTDATA(IOBS,1)-HSAN(IJALT_MDL)
            ENDIF
          ENDIF

          ! U10FG
          ICOUNT = ICOUNT + 1
          IF (LLIN) ZCOMBUF(ICOUNT)=U10FG(IJALT_MDL)

          ! U10AN
          ICOUNT = ICOUNT + 1
          IF (LLIN) ZCOMBUF(ICOUNT)=U10AN(IJALT_MDL)

          ! increase NENTRY if new entry is added
      ENDDO

!       SEND ALL CONTRIBUTIONS OF ZCOMBUF TO PE IRECV

        IF (IZCOMLEN > 0) CALL MPGATHERSCFLD(IRECV,NPEOBS(1,1),NPEOBS(1,2),ZCOMBUF,IZCOMLEN)

!       SETTING UP ODB ON PE IRECV
        IF (IRANK == IRECV) THEN
          !*        GET RALT ODB TO UPDATE IT
          IRET = WAM_YDODB%SELECT('ralt_wam',  ROBSU)

          ZODB_DATE => ROBSU%GET_COLUMN_PTR('date@hdr')
          ZODB_TIME => ROBSU%GET_COLUMN_PTR('time@hdr')
          ZODB_FG_DEPAR => ROBSU%GET_COLUMN_PTR('fg_depar@body')
          ZODB_AN_DEPAR => ROBSU%GET_COLUMN_PTR('an_depar@body')
          ZODB_GP_NUMBER => ROBSU%GET_COLUMN_PTR('gp_number@hdr')
          ZODB_SATELLITE_IDENTIFIER => &
            & ROBSU%GET_COLUMN_PTR('satellite_identifier@sat')
          ZODB_REPORT_STATUS => ROBSU%GET_COLUMN_PTR('report_status@hdr')
          ZODB_DATUM_STATUS => ROBSU%GET_COLUMN_PTR('datum_status@body')
          ZODB_DATUM_EVENT2 => ROBSU%GET_COLUMN_PTR('datum_event2@body')
          ZODB_VARNO => ROBSU%GET_COLUMN_PTR('varno@body')
          ZODB_OBSVALUE => ROBSU%GET_COLUMN_PTR('obsvalue@body')
          ZODB_SEQNO => ROBSU%GET_COLUMN_PTR('seqno@hdr')

          IF (ROBSU%NROWS /= 0) THEN
            IF (ALLOCATED(LUNSET_FLAGS)) DEALLOCATE (LUNSET_FLAGS)
            ALLOCATE (LUNSET_FLAGS(ROBSU%NROWS))
            LUNSET_FLAGS=.TRUE.
            IF (ALLOCATED(LIN_TIMWIN)) DEALLOCATE (LIN_TIMWIN)
            ALLOCATE (LIN_TIMWIN(ROBSU%NROWS))

            DO JOBS = 1, ROBSU%NROWS
              IDATE=NINT(ZODB_DATE(JOBS))
              ITIME=NINT(ZODB_TIME(JOBS))
              IF (IDATE >= 0 .AND. IDATE <= 99991231 .AND.              &
     &            ITIME >= 0 .AND. ITIME <= 999999         ) THEN
                WRITE(COBSDT,'(I8.8,I6.6)') IDATE, ITIME
                LIN_TIMWIN(JOBS)= LGE(COBSDT,CBEGINDT) .AND. LLE(COBSDT,CENDDT)
              ELSE
                LIN_TIMWIN(JOBS)=.FALSE.
              ENDIF

              IF (LIN_TIMWIN(JOBS)) THEN
                ZODB_FG_DEPAR(JOBS)=ODB_RMDI  
                ZODB_AN_DEPAR(JOBS)=ODB_RMDI

              ENDIF
            ENDDO
          ENDIF


!          GET CONTRIBUTION FROM EVERY PEs
           DO IRN = 1, NPROC
             ICOUNT = NPEOBS(IRN,1)
             IRK=NINT(ZCOMBUF(ICOUNT))

             DO IOBS=1,NOBSPE(IRN)

               DO IC=1,NIJALT
                ICOUNT = ICOUNT + 1
                IJALT_LOC(IC)=NINT(ZCOMBUF(ICOUNT))
               ENDDO
               IJALT_MDL=IJALT_LOC(1)
               ISAT_MDL =IJALT_LOC(2)

               IF (IJALT_MDL >= NSTART(IRN) .AND. IJALT_MDL <= NEND(IRN)) THEN
                 LLIN=.TRUE.
               ELSE
                 LLIN=.FALSE.
               ENDIF

               ICOUNT = ICOUNT + 1
               IF (IJALT_LOC(3) == -1) THEN
                 ZFG_DEPART=ODB_RMDI
               ELSE
                 ZFG_DEPART=ZCOMBUF(ICOUNT)
               ENDIF

               ICOUNT = ICOUNT + 1
               IF (IJALT_LOC(3) == -1) THEN
                 ZAN_DEPART=ODB_RMDI
               ELSE
                 ZAN_DEPART=ZCOMBUF(ICOUNT)
               ENDIF

               ICOUNT = ICOUNT + 1
               ZU10FG =ZCOMBUF(ICOUNT)

               ICOUNT = ICOUNT + 1
               ZU10AN=ZCOMBUF(ICOUNT)

               IF (LLIN) THEN
                 I_PARAM=0
                 ! FIND THE RECORD THAT CORRESPONDS TO IOBS
                 DO JOBS = 1, ROBSU%NROWS 
                 IF (LIN_TIMWIN(JOBS)) THEN
                   IJALT_ODB=NINT(ZODB_GP_NUMBER(JOBS))
                   IF (.NOT.LL1D .AND. NPROC > 1) IJALT_ODB=IJ2NEWIJ(IJALT_ODB)
                   IF (IJALT_ODB == IJALT_MDL) THEN
                     ISAT_ODB =NINT(ZODB_SATELLITE_IDENTIFIER(JOBS))
                     IF (ISAT_ODB == ISAT_MDL) THEN
                       IF (NINT(ZODB_VARNO(JOBS)) == 220) THEN
                         ! Select wave height Var.Num.=220
                         ! If $ralt_swh = 220 in file odb/ddl/varno.h changes, you need to change it here
                         ZODB_FG_DEPAR(JOBS)=ZFG_DEPART
                         ZODB_AN_DEPAR(JOBS)=ZAN_DEPART

                         !The report is active (used in assimilation system):
                         ZODB_REPORT_STATUS(JOBS)  = 0 
                         ZODB_REPORT_STATUS(JOBS)  = ZSWITCH_BIT(ZODB_REPORT_STATUS(JOBS),1,0)  ! swith active bit

                         ! datum_status
                         ZODB_DATUM_STATUS(JOBS) = 0 
                         IF (ABS(ZFG_DEPART-ODB_RMDI) > 1.D-5 .AND. ABS(ZAN_DEPART-ODB_RMDI) > 1.D-5) THEN
                           IF (IJALT_LOC(3) < 0) THEN
                             !The observation was rejected
                             ZODB_DATUM_STATUS(JOBS) = ZSWITCH_BIT(ZODB_DATUM_STATUS(JOBS), 1, 2) ! switch rejected bit
                           ELSE IF (IJALT_LOC(3) == 0) THEN
                             !The observation was blacklisted 
                             ZODB_DATUM_STATUS(JOBS) = ZSWITCH_BIT(ZODB_DATUM_STATUS(JOBS), 1, 3) ! switch blacklisted bit
                           ELSE
                             !The observation is active (used in assimilation system):
                             ZODB_DATUM_STATUS(JOBS) = ZSWITCH_BIT(ZODB_DATUM_STATUS(JOBS),1,0)  ! swith active bit
                           ENDIF
                         ELSE
                           !The observation is rejected as there is no FG or AN value
                           ZODB_DATUM_STATUS(JOBS) = ZSWITCH_BIT(ZODB_DATUM_STATUS(JOBS), 1, 2) ! switch rejected bit
                         ENDIF

                         ! datum_event2
                         ZODB_DATUM_EVENT2(JOBS)=IJALT_LOC(3)

                         LUNSET_FLAGS(JOBS)=.FALSE.

                         I_PARAM=I_PARAM+1
                         IF (I_PARAM >= NUM_PARAM) EXIT
                       ELSEIF (NINT(ZODB_VARNO(JOBS)) == 221) THEN
                         ! Select 10-m wind speed Var.Num.=221
                         ! If $ralt_sws = 221 in file odb/ddl/varno.h changes, you need to change it here

                         !Report status is active
                         ZODB_REPORT_STATUS(JOBS)  = 0 
                         ZODB_REPORT_STATUS(JOBS)  = ZSWITCH_BIT(ZODB_REPORT_STATUS(JOBS),1,0)  ! swith active bit

                         ZODB_DATUM_STATUS(JOBS) = 0 
                         IF (ABS(ZODB_OBSVALUE(JOBS)-ODB_RMDI) > 1.D-5) THEN
                           ZODB_FG_DEPAR(JOBS)=ZODB_OBSVALUE(JOBS)-ZU10FG
                           ZODB_AN_DEPAR(JOBS)=ZODB_OBSVALUE(JOBS)-ZU10AN
                           !The observation is passive (for monitoring):
                           ZODB_DATUM_STATUS(JOBS) = ZSWITCH_BIT(ZODB_DATUM_STATUS(JOBS), 1, 1) ! switch passive bit
                         ELSE
                           !The observation is rejected as there is no observation value
                           ZODB_FG_DEPAR(JOBS)=ODB_RMDI
                           ZODB_AN_DEPAR(JOBS)=ODB_RMDI
                           ZODB_DATUM_STATUS(JOBS) = ZSWITCH_BIT(ZODB_DATUM_STATUS(JOBS), 1, 2) ! switch rejected bit
                         ENDIF

                         LUNSET_FLAGS(JOBS)=.FALSE.

                         I_PARAM=I_PARAM+1
                         IF (I_PARAM >= NUM_PARAM) EXIT
                       ENDIF
                     ENDIF
                   ENDIF
                 ENDIF
                 ENDDO
                 IF (I_PARAM < NUM_PARAM) THEN
                    WRITE(IU06,*)'WARNING [WAM2ODB]: NO MATCH-UP MODEL '// &
     &                    'VALUES FOR IOBS,ISAT,PE,INDX,IJS,IJL: ',        &
     &                    IOBS,ISAT_MDL,IRN,IJALT_MDL,IJS,IJL
                 ENDIF
               ENDIF  ! LLIN
             ENDDO
           ENDDO

           ! Setting all observations with unset flags to rejected (but the report is always active!):
           DO JOBS = 1, ROBSU%NROWS
             IF (LIN_TIMWIN(JOBS) .AND. LUNSET_FLAGS(JOBS)) THEN
               ZODB_REPORT_STATUS(JOBS)  = 0 
               ZODB_REPORT_STATUS(JOBS)  = ZSWITCH_BIT(ZODB_REPORT_STATUS(JOBS),1,0)  ! swith active bit
               ZODB_DATUM_STATUS(JOBS) = 0 
               ZODB_DATUM_STATUS(JOBS) = ZSWITCH_BIT(ZODB_DATUM_STATUS(JOBS), 1, 2) ! switch rejected bit
             ENDIF
           ENDDO



           ! Make report_status consistent with datum_status
           ! i.e. if all datum_status for a report are rejected then set report_status as rejected as well
           IPREV_SEQNO = 0
           I_START_OF_THIS_REPORT = 0
           L_ANY_ACTIVE_THIS_REPORT = .FALSE.
           DO JOBS = 1, ROBSU%NROWS
             ! If new report detected
             IF (ZODB_SEQNO(JOBS) /= IPREV_SEQNO) THEN
               IF (IPREV_SEQNO /= 0) THEN
                 IF (.NOT.L_ANY_ACTIVE_THIS_REPORT) THEN
                    ZODB_REPORT_STATUS(I_START_OF_THIS_REPORT:JOBS-1)  = 0
                 ENDIF
               ENDIF

               I_START_OF_THIS_REPORT = JOBS
               IPREV_SEQNO = ZODB_SEQNO(JOBS)
               L_ANY_ACTIVE_THIS_REPORT = .FALSE.
             ENDIF

             ! if datum_status.active == 1 
             IF (IBITS(INT(ZODB_DATUM_STATUS(JOBS)), 0, 1) == 1) L_ANY_ACTIVE_THIS_REPORT=.TRUE.

           ENDDO

           ! Check final report
           IF ((.NOT.L_ANY_ACTIVE_THIS_REPORT).AND.(ROBSU%NROWS>0)) &
             & ZODB_REPORT_STATUS(I_START_OF_THIS_REPORT:JOBS-1)  = 0


           IF(ALLOCATED(ROBSU%DATA)) IRET = WAM_YDODB%PUT(ROBSU)

           IF (ALLOCATED(LUNSET_FLAGS)) DEALLOCATE (LUNSET_FLAGS)
           IF (ASSOCIATED(ZODB_DATE)) NULLIFY(ZODB_DATE)
           IF (ASSOCIATED(ZODB_TIME)) NULLIFY(ZODB_TIME)
           IF (ASSOCIATED(ZODB_FG_DEPAR)) NULLIFY(ZODB_FG_DEPAR)
           IF (ASSOCIATED(ZODB_AN_DEPAR)) NULLIFY(ZODB_AN_DEPAR)
           IF (ASSOCIATED(ZODB_GP_NUMBER)) NULLIFY(ZODB_GP_NUMBER)
           IF (ASSOCIATED(ZODB_SATELLITE_IDENTIFIER)) NULLIFY(ZODB_SATELLITE_IDENTIFIER)
           IF (ASSOCIATED(ZODB_REPORT_STATUS)) NULLIFY(ZODB_REPORT_STATUS)
           IF (ASSOCIATED(ZODB_DATUM_STATUS)) NULLIFY(ZODB_DATUM_STATUS)
           IF (ASSOCIATED(ZODB_DATUM_EVENT2)) NULLIFY(ZODB_DATUM_EVENT2)
           IF (ASSOCIATED(ZODB_VARNO)) NULLIFY(ZODB_VARNO)
           IF (ASSOCIATED(ZODB_OBSVALUE)) NULLIFY(ZODB_OBSVALUE)
           IF (ASSOCIATED(ZODB_SEQNO)) NULLIFY(ZODB_SEQNO)

         ENDIF ! end if on PE IRECV

         IF (ALLOCATED(ZCOMBUF)) DEALLOCATE(ZCOMBUF)

      IF (LHOOK) CALL DR_HOOK('WAM2ODB',1,ZHOOK_HANDLE)
!     -------------------------------------------------------------------------
#endif
      END SUBROUTINE WAM2ODB
