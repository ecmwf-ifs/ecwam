      SUBROUTINE OUTGRID (IJS, IJL, BOUT)
! ----------------------------------------------------------------------

!**** *OUTGRID* - OUTPUT OF WAVE ENERGY,MEAN ANGLE, MEAN FREQUENCY
!****             FRICTION VELOCITY AND WIND DIRECTION IN GRIDDED
!****             FORMAT.

!     P.JANSSEN      KNMI


!     P.LIONELLO  FEB. 87
!                 OUTPUT OF SWELL ENERGY ,MEAN SWELL DIRECTION ,
!                 MEAN SWELL FREQUENCY AND MEAN WIND-SEA WAVE
!                 DIRECTION AT ALL ACTIVE GRID POINTS

!     H.GUNTHER   GKSS/ECMWF     NOVEMBER 1989
!     J. BIDLOT   ECMWF          APRIL 1996
!                 MESSAGE PASSING

!                 OUTPUT FORMAT IS NOW GRIDDED FIELDS.

!*    PURPOSE.
!     --------

!       OUTPUT OF WAVE AND WIND FIELDS INTO COMMON BLOCK.


!**   INTERFACE.
!     ----------

!        *CALL* *OUTGRID (IJS, IJL, BOUT)
!         *IJS*    - INDEX OF FIRST LOCAL GRIDPOINT.
!         *IJL*    - INDEX OF LAST LOCAL GRIDPOINT.
!         *BOUT*   - BLOCk OF OUTPUT PARAMETERS


!     EXTERNALS.
!     ----------

!       *MAKEGRID*  - GENERATES GRIDDED FIELDS.
!       *MPGATHERSCFLD* - GATHER SCALAR FIELD DISTRIBUTED ON DIFFERENT
!                         PROCESSES

!     METHOD.
!     -------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE MPL_MPIF
      USE yowDatapool, only: comm
      USE OUTPUT_STRUCT, only : INITIAL_OUTPUT_INITS_NEXTGEN
      USE OUTPUT_STRUCT, only : SET_UP_ARR_OUT_RECV
      USE YOWPD, only : INE, npa, NP_RES =>np
      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal
      USE YOWCOUT  , ONLY : NTRAIN   ,JPPFLAG  ,IPFGTBL,                &
     &            NIPRMOUT, ITOBOUT
      USE YOWCURR  , ONLY : U, V 
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE ,CICOVER  ,CITHRSH  ,
     &            CIBLOCK ,CITHICK 
      USE YOWINTP  , ONLY : GOUT
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,MPMAXLENGTH        ,
     &            KTAG
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO     ,NIBLO    ,    &
     &                      NGX      ,NGY
      USE YOWPCONS , ONLY : ZMISS    ,EPSUS    ,EPSU10
      USE YOWSHAL  , ONLY : DEPTH
      USE YOWSPEC  , ONLY : NBLKS   ,NBLKE     , 
     &            U10NEW   ,THWNEW   ,USNEW    ,
     &            ROAIRN   ,ZIDLNEW
      USE YOWSWEL  , ONLY : ESWELL   ,FSWELL   ,THSWELL  ,ESEA     ,
     &            FSEA     ,THWISEA  ,P1SWELL  ,P2SWELL  ,SPRDSWELL,
     &            P1SEA    ,P2SEA    ,SPRDSEA
      USE YOWTRAINS, ONLY : EMTRAIN  ,THTRAIN  ,PMTRAIN
      USE YOWTEST  , ONLY : IU06     
      USE YOWNEMOFLDS,ONLY: NEMOSST, NEMOCICOVER, NEMOCITHICK, 
     &                      NEMOUCUR, NEMOVCUR, LNEMOCITHICK
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE OUTPUT_STRUCT, ONLY : ARR_OUT_RECV, LocalPosICT
      USE YOWUNPOOL, ONLY : NIBLO_FD, NIBLO_OUT
      USE YOWUNPOOL, ONLY : OUT_METHOD
      USE YOWUNPOOL, ONLY : LLUNSTR
      USE MPL_MODULE

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL
      REAL, INTENT(IN) :: BOUT(IJS:IJL,JPPFLAG)

      INTEGER, PARAMETER :: IG=1
      INTEGER :: IJ
      INTEGER :: ITT,ICT
      INTEGER :: NFLDTOT, NFLDPPEMAX
      INTEGER :: ITG
      INTEGER :: IFLD
      INTEGER :: IR, IPR
      INTEGER :: ICOUNT
      INTEGER, DIMENSION(NPROC) :: ICNT
      INTEGER, DIMENSION(NPROC+JPPFLAG) :: IREQ
      INTEGER, DIMENSION(NPROC) :: NFLDPPE 
      INTEGER, DIMENSION(NPROC) :: ISENDCOUNTS,IRECVCOUNTS

      REAL :: ZHOOK_HANDLE
      REAL :: ZDUM(2)
      REAL :: ZCOMBUFR(MPMAXLENGTH)
      REAL, ALLOCATABLE, DIMENSION(:) :: GTEMP 
      REAL, ALLOCATABLE, DIMENSION(:,:) :: ZSENDBUF, ZRECVBUF 

      LOGICAL :: LSQRT
      LOGICAL, SAVE :: HaveMPI_arrays = .FALSE.
!----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTGRID',0,ZHOOK_HANDLE)
#endif

!     FIND THE NUMBER OF FIELDS EACH PE HAS TO DEAL WITH
      NFLDPPE(:)=0
      NFLDTOT=0
      NFLDPPEMAX=1
      DO ICT=1,JPPFLAG
        IF(IPFGTBL(ICT).NE.0) THEN 
          NFLDPPE(IPFGTBL(ICT))=NFLDPPE(IPFGTBL(ICT))+1
          NFLDTOT=NFLDTOT+1
          NFLDPPEMAX=MAX(NFLDPPEMAX,NFLDPPE(IPFGTBL(ICT)))
        ENDIF
      ENDDO
      IF (HaveMPI_arrays .eqv. .FALSE.) THEN
        IF (LLUNSTR .and. (OUT_METHOD .eq. 2)) THEN
           CALL INITIAL_OUTPUT_INITS_NEXTGEN
        END IF
      END IF
      HaveMPI_arrays=.TRUE.
      IF(NFLDTOT.EQ.0) THEN
#ifdef ECMWF
        IF (LHOOK) CALL DR_HOOK('OUTGRID',1,ZHOOK_HANDLE)
#endif
        RETURN
      ENDIF
 


!     SENDING TO RELEVANT PE'S
!     ------------------------

      DO IPR=1,NPROC
        ISENDCOUNTS(IPR)=NFLDPPE(IPR)*(NBLKE(IRANK)-NBLKS(IRANK)+1)
      ENDDO

      DO IPR=1,NPROC
        IRECVCOUNTS(IPR)=NFLDPPE(IRANK)*(NBLKE(IPR)-NBLKS(IPR)+1)
      ENDDO

      IF ((.NOT.LLUNSTR).OR.(OUT_METHOD.eq.1)) THEN
      
!     LOADING THE COMMUNICATION BUFFER
        ALLOCATE(ZSENDBUF(NFLDPPEMAX * MPMAXLENGTH,NPROC))

        ICNT(:)=0
        DO ICT=1,JPPFLAG
          IF(IPFGTBL(ICT).NE.0) THEN 
            IPR=IPFGTBL(ICT)
            DO IJ=IJS,IJL
              ICNT(IPR) = ICNT(IPR) + 1
              ZSENDBUF(ICNT(IPR),IPR) = BOUT(IJ,ITOBOUT(ICT))
            ENDDO
          ENDIF
        ENDDO
      ELSE
!!!!! this will need to be adapted to use BOUT
        WRITE(0,*) '!!! ********************************* !!'
        WRITE(0,*) '!!! in outgrid. Not yet ready !!!' 
        WRITE(0,*) '!!! ********************************* !!'
        CALL ABORT1

        CALL SET_UP_ARR_OUT_RECV(IJS, IJL, BOUT, NFLDPPE)
      END IF


!     GLOBAL EXCHANGE

      IF ((.NOT.LLUNSTR).OR.(OUT_METHOD.eq.1)) THEN
        CALL GSTATS(693,0)
        IR=0
        IF(NFLDPPE(IRANK).GT.0) THEN
          ALLOCATE(ZRECVBUF(NFLDPPE(IRANK)*MPMAXLENGTH,NPROC))
          DO IPR=1,NPROC
            IR=IR+1
            CALL MPL_RECV(ZRECVBUF(1:IRECVCOUNTS(IPR),IPR),KSOURCE=IPR,
     &        KTAG=KTAG,
     &        KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ(IR),
     &        CDSTRING='OUTGRID:')
          ENDDO
        ENDIF
        DO IPR=1,NPROC
          IF(NFLDPPE(IPR).GT.0) THEN
            IR=IR+1
            CALL MPL_SEND(ZSENDBUF(1:ICNT(IPR),IPR),KDEST=IPR,
     &        KTAG=KTAG,
     &        KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ(IR),
     &        CDSTRING='OUTGRID:')
          ENDIF
        ENDDO
!       NOW WAIT FOR ALL TO COMPLETE
        CALL MPL_WAIT(KREQUEST=IREQ(1:IR),CDSTRING='OUTGRID:')
        CALL GSTATS(693,1)
      END IF


!     RETRIEVE THE INFORMATION
!     ------------------------

      IF (LLUNSTR .and. (OUT_METHOD .eq. 2)) THEN
        NIBLO_OUT = NIBLO_FD
      ELSE
        NIBLO_OUT = NIBLO
      END IF
      ALLOCATE(GTEMP(NIBLO_OUT))


!     RECEIVING AND TRANSFERING FROM BLOCK TO GRID
!     --------------------------------------------

      IF(NFLDPPE(IRANK).GT.0) THEN
        IF(ALLOCATED(GOUT)) DEALLOCATE(GOUT)
        ALLOCATE(GOUT(NFLDPPE(IRANK),NGX,NGY))
      ENDIF

      ICNT(:)=0
      ICT=1
      IFLD=0
      DO WHILE ( IFLD.LT.NFLDPPE(IRANK) .AND. ICT .LE. JPPFLAG )
        IF(IPFGTBL(ICT).EQ.IRANK ) THEN
          IF (.NOT.(LLUNSTR) .OR. (OUT_METHOD .eq. 1)) THEN
            DO IPR=1,NPROC
              ICOUNT=ICNT(IPR)
              DO IJ=NBLKS(IPR),NBLKE(IPR)
                ICOUNT=ICOUNT+1
                GTEMP(IJ)=ZRECVBUF(ICOUNT,IPR)
              ENDDO
              ICNT(IPR)=ICOUNT
            ENDDO
          ELSE
            GTEMP=ARR_OUT_RECV(LocalPosICT(ICT),:)
          END IF

          LSQRT=.FALSE.
          CALL MAKEGRID (GTEMP,GOUT(ITOBOUT(ICT),:,:),IG,ZMISS)
          IFLD=IFLD+1

        ENDIF ! (IPFGTBL) 

        ICT=ICT+1
      ENDDO  ! WHILE

      DEALLOCATE(GTEMP)
      IF(ALLOCATED(ZRECVBUF)) DEALLOCATE(ZRECVBUF)
      IF(ALLOCATED(ZSENDBUF)) DEALLOCATE(ZSENDBUF)

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTGRID',1,ZHOOK_HANDLE)
#endif
      END SUBROUTINE OUTGRID
