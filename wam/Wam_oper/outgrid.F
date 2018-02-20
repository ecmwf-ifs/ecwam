      SUBROUTINE OUTGRID (IJS, IJL, IJ_OFFSET, BOUT)
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

!        *CALL* *OUTGRID (IJS, IJL, IJ_OFFSET, BOUT)
!         *IJS*    - INDEX OF FIRST LOCAL GRIDPOINT.
!         *IJL*    - INDEX OF LAST LOCAL GRIDPOINT.
!         *IJ_OFFSET* OFFSET to point IJSLOC and IJLLOC to the global block of data
!                     only meaningful if unstructured grid
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
      USE YOWCOUT  , ONLY : NTRAIN   ,JPPFLAG  ,IPFGTBL
      USE YOWCURR  , ONLY : U, V 
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE ,CICOVER  ,CITHRSH  ,
     &            CIBLOCK ,CITHICK 
      USE YOWINTP  , ONLY : WHGTTG   ,WDIRTG   ,WPKFTG   ,WMNFTG   ,
     &            USTARG   ,UDIRG    ,TAUWG    ,CDG      ,SMEANG   ,
     &            U10G     ,MWP1G    ,MWP2G    ,WSPRDG   ,C4G      ,
     &            C3G      ,WHGTAG   ,CWHTAG   ,RANGAG   ,
     &            BFG      ,QPG      ,DEPTHG   ,HMAXG    ,TMAXG    ,
     &            USTOKESG ,VSTOKESG ,UCURG    ,VCURG    ,PHIEPSG  ,
     &            PHIAWG   ,TAUOCG   ,STRNMSG  ,RHOAG    ,WSTARG   ,
     &            CICG     ,CIHG     ,
     &            E10G     ,WEFLXMG  ,WEFLXDG  ,E1012G   , 
     &            E1214G   ,E1417G   ,E1721G   ,E2125G   ,E2530G   ,
     &            WX1G     ,WX2G     ,WX3G     ,WX4G     ,WX5G     ,
     &            NSSTG    ,NICCG    ,NICTG    ,NUCUG    ,NVCUG
      USE YOWINTS  , ONLY : WHGTSG   ,WDIRSG   ,WMNFSG   ,WHGTWG   ,
     &            WDIRWG   ,WMNFWG   ,MWP1SG   ,MWP2SG   ,WSPRDSG  ,
     &            MWP1WG   ,MWP2WG   ,WSPRDWG
      USE YOWINTT  , ONLY : WHGTTRG  ,WDIRTRG  ,MWPMTRG
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,THQ      ,CDMEAN   ,
     &            WVSTRMEAN,SMEAN    ,ALTWH    ,CALTWH   ,RALTCOR  ,
     &            P1MEAN   ,P2MEAN   ,SPRDMEAN ,C3MEAN   ,C4MEAN   ,
     &            BFMEAN   ,QPMEAN   ,HMAXMEAN ,FPMEAN   ,TMAXMEAN ,
     &            USTMEAN  ,VSTMEAN  ,PHIEPS   ,PHIAW    ,TAUOC    ,
     &            STRNMS   ,E10MEAN  ,WEFLXM   ,WEFLXD   ,E1012    ,
     &            E1214    ,E1417    ,E1721    ,E2125    ,E2530    ,
     &            WX1      ,WX2      ,WX3      ,WX4      ,WX5
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,MPMAXLENGTH        ,
     &            KTAG
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO     ,NIBLO
      USE YOWPCONS , ONLY : ZMISS    ,EPSUS    ,EPSU10
      USE YOWSHAL  , ONLY : DEPTH
      USE YOWSPEC  , ONLY : NBLKS   ,NBLKE     , 
     &            U10NEW   ,THWNEW   ,USNEW    ,
     &            ROAIRN   ,ZIDLNEW
      USE YOWSTAT  , ONLY : CDTPRO   ,CDTINTS
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

      INTEGER, INTENT(IN) :: IJS, IJL, IJ_OFFSET
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
      REAL, ALLOCATABLE, DIMENSION(:,:) :: TEMP

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
        IF(ICT.LE.16 .OR. ICT.GE.22) THEN
          IF(IPFGTBL(ICT).NE.0) THEN 
            NFLDPPE(IPFGTBL(ICT))=NFLDPPE(IPFGTBL(ICT))+1
            NFLDTOT=NFLDTOT+1
            NFLDPPEMAX=MAX(NFLDPPEMAX,NFLDPPE(IPFGTBL(ICT)))
          ENDIF
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
          IF(ICT.LE.16 .OR. ICT.GE.22) THEN
            IF(IPFGTBL(ICT).NE.0) THEN 
              IPR=IPFGTBL(ICT)
              DO IJ=IJS,IJL
                ICNT(IPR) = ICNT(IPR) + 1
                ZSENDBUF(ICNT(IPR),IPR) = TEMP(IJ,ICT)
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ELSE
        CALL SET_UP_ARR_OUT_RECV(IJS, IJL, TEMP, NFLDPPE)
      END IF

      DEALLOCATE(TEMP)

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

      ICNT(:)=0
      ICT=1
      IFLD=0
      DO WHILE ( IFLD.LT.NFLDPPE(IRANK) .AND. ICT .LE. JPPFLAG )
        IF((ICT.LE.16 .OR. ICT.GE.22).AND. IPFGTBL(ICT).EQ.IRANK ) THEN
          IF (.NOT.(LLUNSTR) .or. (OUT_METHOD .eq. 1)) THEN
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

!!! until I move all grid arrays to one work array, need to loop on all
!!! special case
          IF(ICT.EQ.1) THEN
            LSQRT=.TRUE.
            CALL MAKEGRID (GTEMP,WHGTTG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.2) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WDIRTG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.3) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WMNFTG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.4) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,USTARG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.5) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,UDIRG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.6) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WPKFTG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.7) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,CDG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.8) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,TAUWG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.9) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,SMEANG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.10) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,U10G,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.11) THEN
            LSQRT=.TRUE.
            CALL MAKEGRID (GTEMP,WHGTWG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.12) THEN
            LSQRT=.TRUE.
            CALL MAKEGRID (GTEMP,WHGTSG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.13) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WDIRWG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.14) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WDIRSG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.15) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WMNFWG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.16) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WMNFSG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.22) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WHGTAG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.23) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,CWHTAG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.24) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,RANGAG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.25) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,MWP1G,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.26) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,MWP2G,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.27) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WSPRDG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.28) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,MWP1WG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.29) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,MWP1SG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.30) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,MWP2WG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.31) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,MWP2SG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.32) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WSPRDWG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.33) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WSPRDSG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.34) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,C4G,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.35) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,BFG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.36) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,QPG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.37) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,DEPTHG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.38) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,HMAXG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.39) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,TMAXG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.40) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,USTOKESG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.41) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,VSTOKESG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.42) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,UCURG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.43) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,VCURG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.44) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,PHIEPSG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.45) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,PHIAWG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.46) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,TAUOCG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ENDIF

          DO ITT = 1, NTRAIN
            IF(ICT.EQ.(46+(ITT-1)*3+1) ) THEN
              LSQRT=.TRUE.
              CALL MAKEGRID (GTEMP,WHGTTRG(:,:,ITT),
     &                   IG,LSQRT,ZMISS)
              IFLD=IFLD+1

            ELSEIF(ICT.EQ.(46+(ITT-1)*3+2) ) THEN
              LSQRT=.FALSE.
              CALL MAKEGRID (GTEMP,WDIRTRG(:,:,ITT),
     &                   IG,LSQRT,ZMISS)
              IFLD=IFLD+1

            ELSEIF(ICT.EQ.(46+(ITT-1)*3+3) ) THEN
              LSQRT=.FALSE.
              CALL MAKEGRID (GTEMP,MWPMTRG(:,:,ITT),
     &                   IG,LSQRT,ZMISS)
              IFLD=IFLD+1
            ENDIF
          ENDDO

          IF(ICT.EQ.(46+3*NTRAIN+1)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,STRNMSG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+2)) THEN
            LSQRT=.TRUE.
            CALL MAKEGRID (GTEMP,E10G,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+3)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,RHOAG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+4)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WSTARG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+5)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,CICG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+6)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,CIHG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+7)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,C3G,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+8)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,NSSTG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+9)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,NICCG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+10)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,NICTG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+11)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,NUCUG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+12)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,NVCUG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+13)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WEFLXMG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+14)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WEFLXDG,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+15)) THEN
            LSQRT=.TRUE.
            CALL MAKEGRID (GTEMP,E1012G,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+16)) THEN
            LSQRT=.TRUE.
            CALL MAKEGRID (GTEMP,E1214G,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+17)) THEN
            LSQRT=.TRUE.
            CALL MAKEGRID (GTEMP,E1417G,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+18)) THEN
            LSQRT=.TRUE.
            CALL MAKEGRID (GTEMP,E1721G,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+19)) THEN
            LSQRT=.TRUE.
            CALL MAKEGRID (GTEMP,E2125G,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+20)) THEN
            LSQRT=.TRUE.
            CALL MAKEGRID (GTEMP,E2530G,IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ENDIF

!       EXTRA OUTPUT FIELDS:
!       See *OUTBS* where they should have been calculated or pick them up from
!       the forcing fields (see above if the sea ice mask has to be removed)
          IF(ICT.EQ.(JPPFLAG-4)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WX1G,IG,LSQRT,0.)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(JPPFLAG-3)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WX2G,IG,LSQRT,0.)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(JPPFLAG-2)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WX3G,IG,LSQRT,0.)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(JPPFLAG-1)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WX4G,IG,LSQRT,0.)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(JPPFLAG)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WX5G,IG,LSQRT,0.)
            IFLD=IFLD+1
          ENDIF

        ENDIF ! ( <=16 >=22  and IPFGTBL) 

        ICT=ICT+1
      ENDDO  ! WHILE

      DEALLOCATE(GTEMP)
      IF(ALLOCATED(ZRECVBUF)) DEALLOCATE(ZRECVBUF)
      IF(ALLOCATED(ZSENDBUF)) DEALLOCATE(ZSENDBUF)

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTGRID',1,ZHOOK_HANDLE)
#endif
      END SUBROUTINE OUTGRID
