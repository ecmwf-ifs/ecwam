      SUBROUTINE OUTGRID (IJS, IJL, IJ_OFFSET, IODP)
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

!        *CALL* *OUTGRID (IJS, IJL, IJ_OFFSET, IODP)
!         *IJS*    - INDEX OF FIRST LOCAL GRIDPOINT.
!         *IJL*    - INDEX OF LAST LOCAL GRIDPOINT.
!         *IJ_OFFSET* OFFSET to point IJSLOC and IJLLOC to the global block of data
!                     only meaningful if unstructured grid
!         *IODP*     LAND MASK IF IODP(IJ)=0


!     EXTERNALS.
!     ----------

!       *SETWMASK*   - TRANSFERS A LOCAL BLOCK ARRAY TO A GLOBAL BLOCK 
!                     ARRAY AND IMPOSES ICE MASK.
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
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,THQ      ,
     &            SMEAN    ,ALTWH    ,CALTWH   ,RALTCOR  ,
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
     &            U10NEW   ,THWNEW   ,USNEW    ,TAUW     ,
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

      INTEGER, PARAMETER :: IG=1
      INTEGER, INTENT(IN) :: IJS, IJL, IJ_OFFSET
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
      INTEGER, DIMENSION(IJS:IJL) :: IODP
      INTEGER, DIMENSION(IJS:IJL) :: NIODP

      REAL :: TAU 
      REAL :: ZHOOK_HANDLE
      REAL :: ZDUM(2)
      REAL :: CD(IJS:IJL)
      REAL :: ZCOMBUFR(MPMAXLENGTH)
      REAL, ALLOCATABLE, DIMENSION(:) :: GTEMP 
      REAL, ALLOCATABLE, DIMENSION(:) :: WVSTR
      REAL, ALLOCATABLE, DIMENSION(:,:) :: ZSENDBUF, ZRECVBUF 
      REAL, ALLOCATABLE, DIMENSION(:,:) :: TEMP

      LOGICAL :: LSQRT
      LOGICAL, SAVE :: HaveMPI_arrays = .FALSE.
!----------------------------------------------------------------------
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'Begin OUTGRID'
      FLUSH(740+MyRankGlobal)
#endif
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
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'OUT_METHOD=', OUT_METHOD
      FLUSH(740+MyRankGlobal)
#endif
      IF (HaveMPI_arrays .eqv. .FALSE.) THEN
        IF (LLUNSTR .and. (OUT_METHOD .eq. 2)) THEN
#ifdef DEBUG
           WRITE(740+MyRankGlobal,*) 'Before call to NEXTGEN'
#endif
           CALL INITIAL_OUTPUT_INITS_NEXTGEN
#ifdef DEBUG
           WRITE(740+MyRankGlobal,*) 'After call to NEXTGEN'
#endif
        END IF
      END IF
      HaveMPI_arrays=.TRUE.
      IF(NFLDTOT.EQ.0) THEN
#ifdef ECMWF
        IF (LHOOK) CALL DR_HOOK('OUTGRID',1,ZHOOK_HANDLE)
#endif
        RETURN
      ENDIF
 
      ALLOCATE(TEMP(IJS:IJL,JPPFLAG))
      NIODP(:)=1

!*    1. INTEGRATED PARAMETERS OF TOTAL SEA.
!        -----------------------------------

        IF(IPFGTBL(8).NE.0) ALLOCATE(WVSTR(IJS:IJL))
        DO IJ = IJS,IJL
          TAU = MAX(USNEW(IJ)**2,EPSUS)
          IF(IPFGTBL(8).NE.0) WVSTR(IJ) = TAUW(IJ,IG)/TAU
!!       if the numerical computation of TAU and CD change, a similar
!!       modification has to be put in buildstress where the friction
!!       velocity is determined from U10 and CD.
          CD(IJ) = TAU/MAX(U10NEW(IJ)**2,EPSU10)
        ENDDO

!*    1.3 APPLY SEA ICE MASK (IF NEEDED) AND SEND TO DEDICATED OUTPUT PE
!         -------------------------------------------------

        IF(IPFGTBL(1).NE.0) THEN 
         CALL SETWMASK(EMEAN(IJS),TEMP(IJS,1),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(2).NE.0) THEN 
          CALL SETWMASK(THQ(IJS),TEMP(IJS,2),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(3).NE.0) THEN 
         CALL SETWMASK(FMEAN(IJS),TEMP(IJS,3),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(4).NE.0) THEN
!         The friction velocity is NOT set to missing over ice 
          CALL SETWMASK(USNEW(IJS),TEMP(IJS,4),IJS,IJL,1.1,IODP(IJS))
        ENDIF

        IF(IPFGTBL(5).NE.0) THEN 
!         The wind direction is NOT set to missing over ice or land 
          CALL SETWMASK(THWNEW(IJS),TEMP(IJS,5),IJS,IJL,1.1,NIODP(IJS))
        ENDIF

        IF(IPFGTBL(6).NE.0) THEN
          CALL SETWMASK(FPMEAN,TEMP(IJS,6),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(7).NE.0) THEN 
!         The drag coefficient is NOT set to missing over ice or land 
          CALL SETWMASK(CD,TEMP(IJS,7),IJS,IJL,1.1,NIODP(IJS))
        ENDIF

        IF(IPFGTBL(8).NE.0) THEN
!         The wave stress is NOT set to missing over ice 
          CALL SETWMASK(WVSTR,TEMP(IJS,8),IJS,IJL,1.1,IODP(IJS))
          DEALLOCATE(WVSTR)
        ENDIF

        IF(IPFGTBL(9).NE.0) THEN
         CALL SETWMASK(SMEAN(IJS),TEMP(IJS,9),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(10).NE.0) THEN
!         The wind speed is NOT set to missing over ice or land 
          CALL SETWMASK(U10NEW(IJS),TEMP(IJS,10),IJS,IJL,1.1,NIODP(IJS))
#ifdef DEBUG
          WRITE(740+MyRankGlobal,*) 'minval(U10NEW)=',                  &
     &       minval(U10NEW(IJS:IJL))
          WRITE(740+MyRankGlobal,*) 'maxval(U10NEW)=',                  &
     &       maxval(U10NEW(IJS:IJL))
          FLUSH(740+MyRankGlobal)
#endif          
        ENDIF

        IF(IPFGTBL(11).NE.0) THEN
          CALL SETWMASK(ESEA,TEMP(IJS,11),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(12).NE.0) THEN
          CALL SETWMASK(ESWELL,TEMP(IJS,12),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(13).NE.0) THEN
          CALL SETWMASK(THWISEA,TEMP(IJS,13),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(14).NE.0) THEN
          CALL SETWMASK(THSWELL,TEMP(IJS,14),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(15).NE.0) THEN
          CALL SETWMASK(FSEA,TEMP(IJS,15),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(16).NE.0) THEN
          CALL SETWMASK(FSWELL,TEMP(IJS,16),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(22).NE.0) THEN
          CALL SETWMASK(ALTWH,TEMP(IJS,22),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(23).NE.0) THEN
          CALL SETWMASK(CALTWH,TEMP(IJS,23),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(24).NE.0) THEN
          CALL SETWMASK(RALTCOR,TEMP(IJS,24),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(25).NE.0) THEN
          CALL SETWMASK(P1MEAN,TEMP(IJS,25),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(26).NE.0) THEN
          CALL SETWMASK(P2MEAN,TEMP(IJS,26),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(27).NE.0) THEN
          CALL SETWMASK(SPRDMEAN,TEMP(IJS,27),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(28).NE.0) THEN
          CALL SETWMASK(P1SEA,TEMP(IJS,28),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(29).NE.0) THEN
          CALL SETWMASK(P1SWELL,TEMP(IJS,29),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(30).NE.0) THEN
          CALL SETWMASK(P2SEA,TEMP(IJS,30),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(31).NE.0) THEN
          CALL SETWMASK(P2SWELL,TEMP(IJS,31),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(32).NE.0) THEN
          CALL SETWMASK(SPRDSEA,TEMP(IJS,32),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(33).NE.0) THEN
         CALL SETWMASK(SPRDSWELL,TEMP(IJS,33),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(34).NE.0) THEN
          CALL SETWMASK(C4MEAN,TEMP(IJS,34),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(35).NE.0) THEN
          CALL SETWMASK(BFMEAN,TEMP(IJS,35),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(36).NE.0) THEN
          CALL SETWMASK(QPMEAN,TEMP(IJS,36),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(37).NE.0) THEN
!         The bathymetry is NOT set to missing over ice. 
         CALL SETWMASK(DEPTH(IJS,IG),TEMP(IJS,37),IJS,IJL,1.1,IODP(IJS))
        ENDIF

        IF(IPFGTBL(38).NE.0) THEN
          CALL SETWMASK(HMAXMEAN,TEMP(IJS,38),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(39).NE.0) THEN
          CALL SETWMASK(TMAXMEAN,TEMP(IJS,39),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(40).NE.0) THEN
          CALL SETWMASK(USTMEAN,TEMP(IJS,40),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(41).NE.0) THEN
          CALL SETWMASK(VSTMEAN,TEMP(IJS,41),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(42).NE.0) THEN
!         the ice mask does not apply to the imposed ocean current 
          IF(ALLOCATED(U)) THEN
            CALL SETWMASK(U(IJS,IG),TEMP(IJS,42),IJS,IJL,1.1,IODP(IJS))
          ELSE
            TEMP(:,42)=0.0
          ENDIF
        ENDIF

        IF(IPFGTBL(43).NE.0) THEN
!         the ice mask does not apply to the imposed ocean current 
          IF(ALLOCATED(V)) THEN
            CALL SETWMASK(V(IJS,IG),TEMP(IJS,43),IJS,IJL,1.1,IODP(IJS))
          ELSE
            TEMP(:,43)=0.0
          ENDIF
        ENDIF

        IF(IPFGTBL(44).NE.0) THEN
!         The energy flux into ocean is NOT set to missing over ice 
          CALL SETWMASK(PHIEPS(IJS),TEMP(IJS,44),IJS,IJL,1.1,IODP(IJS))
        ENDIF

        IF(IPFGTBL(45).NE.0) THEN
!         The energy flux into waves is NOT set to missing over ice 
          CALL SETWMASK(PHIAW(IJS),TEMP(IJS,45),IJS,IJL,1.1,IODP(IJS))
        ENDIF

        IF(IPFGTBL(46).NE.0) THEN
!         The momentum flux into ocean is NOT set to missing over ice 
          CALL SETWMASK(TAUOC(IJS),TEMP(IJS,46),IJS,IJL,1.1,IODP(IJS))
        ENDIF

        ICT=46
        DO ITT = 1, NTRAIN
          ICT=ICT+1
          IF(IPFGTBL(ICT).NE.0) THEN
          CALL SETWMASK(EMTRAIN(:,ITT),TEMP(IJS,ICT),
     &                    IJS,IJL,CITHRSH,IODP(IJS))
          ENDIF

          ICT=ICT+1
          IF(IPFGTBL(ICT).NE.0) THEN
          CALL SETWMASK(THTRAIN(:,ITT),TEMP(IJS,ICT),
     &                    IJS,IJL,CITHRSH,IODP(IJS))
          ENDIF

          ICT=ICT+1
          IF(IPFGTBL(ICT).NE.0) THEN
          CALL SETWMASK(PMTRAIN(:,ITT),TEMP(IJS,ICT),
     &                    IJS,IJL,CITHRSH,IODP(IJS))
          ENDIF

        ENDDO


        ITG=46+3*NTRAIN+1
        IF(IPFGTBL(ITG).NE.0) THEN
!         the ice mask does not apply to the mean square strain in the ice
          CALL SETWMASK(STRNMS(IJS),TEMP(IJS,ITG),IJS,IJL,1.1,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
         CALL SETWMASK(E10MEAN(IJS),TEMP(IJS,ITG),
     &                 IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
!         the ice mask does not apply to air density or land 
         CALL SETWMASK(ROAIRN(IJS),TEMP(IJS,ITG),IJS,IJL,1.1,NIODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
!         the ice mask does not apply to the imposed convective velocity scale
        CALL SETWMASK(ZIDLNEW(IJS),TEMP(IJS,ITG),IJS,IJL,1.1,NIODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
!         the ice mask does not apply
          CALL SETWMASK(CICOVER(IJS,IG),TEMP(IJS,ITG),
     &                  IJS,IJL,1.1,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETWMASK(CITHICK(IJS,IG),TEMP(IJS,ITG),
     &                  IJS,IJL,1.1,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETWMASK(C3MEAN,TEMP(IJS,ITG),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          DO IJ = IJS,IJL
            TEMP(IJ,ITG) = NEMOSST(IJ) 
          ENDDO
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          DO IJ = IJS,IJL
            TEMP(IJ,ITG) = NEMOCICOVER(IJ) 
          ENDDO
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          IF (LNEMOCITHICK) THEN
            DO IJ = IJS,IJL
              TEMP(IJ,ITG) = NEMOCITHICK(IJ) 
            ENDDO
          ELSE
            DO IJ = IJS,IJL
              TEMP(IJ,ITG) = ZMISS
            ENDDO
          ENDIF
        ENDIF
          
        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          DO IJ = IJS,IJL
            TEMP(IJ,ITG) = NEMOUCUR(IJ) 
          ENDDO
        ENDIF
        
        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          DO IJ = IJS,IJL
            TEMP(IJ,ITG) = NEMOVCUR(IJ) 
          ENDDO
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
         CALL SETWMASK(WEFLXM(IJS),TEMP(IJS,ITG),
     &                 IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
         CALL SETWMASK(WEFLXD(IJS),TEMP(IJS,ITG),
     &                 IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
         CALL SETWMASK(E1012(IJS),TEMP(IJS,ITG),
     &                 IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
         CALL SETWMASK(E1214(IJS),TEMP(IJS,ITG),
     &                 IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
         CALL SETWMASK(E1417(IJS),TEMP(IJS,ITG),
     &                 IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
         CALL SETWMASK(E1721(IJS),TEMP(IJS,ITG),
     &                 IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
         CALL SETWMASK(E2125(IJS),TEMP(IJS,ITG),
     &                 IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
         CALL SETWMASK(E2530(IJS),TEMP(IJS,ITG),
     &                 IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF


!       EXTRA OUTPUT FIELDS:
!       See *OUTBS* where they should have been calculated or pick them up from
!       the forcing fields (see above if the sea ice mask has to be removed)
        ITG=JPPFLAG-4
        IF(IPFGTBL(ITG).NE.0) THEN
         CALL SETWMASK(WX1(IJS),TEMP(IJS,ITG),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF
        ITG=JPPFLAG-3
        IF(IPFGTBL(ITG).NE.0) THEN
         CALL SETWMASK(WX2(IJS),TEMP(IJS,ITG),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF
        ITG=JPPFLAG-2
        IF(IPFGTBL(ITG).NE.0) THEN
         CALL SETWMASK(WX3(IJS),TEMP(IJS,ITG),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF
        ITG=JPPFLAG-1
        IF(IPFGTBL(ITG).NE.0) THEN
         CALL SETWMASK(WX4(IJS),TEMP(IJS,ITG),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF
        ITG=JPPFLAG
        IF(IPFGTBL(ITG).NE.0) THEN
         CALL SETWMASK(WX5(IJS),TEMP(IJS,ITG),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'After construction of TEMP'
      FLUSH(740+MyRankGlobal)
#endif          


!     SENDING TO RELEVANT PE'S
!     ------------------------

      DO IPR=1,NPROC
        ISENDCOUNTS(IPR)=NFLDPPE(IPR)*(NBLKE(IRANK)-NBLKS(IRANK)+1)
      ENDDO

      DO IPR=1,NPROC
        IRECVCOUNTS(IPR)=NFLDPPE(IRANK)*(NBLKE(IPR)-NBLKS(IPR)+1)
      ENDDO

#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'Before exchanges'
      FLUSH(740+MyRankGlobal)
#endif          
      
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
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'After construction of ARR_OUT_RECV'
      WRITE(740+MyRankGlobal,*) 'allocated(TEMP)=', allocated(TEMP)
      FLUSH(740+MyRankGlobal)
#endif          
      DEALLOCATE(TEMP)
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'After deallocation'
      FLUSH(740+MyRankGlobal)
#endif          

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

#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'Now leaving OUTGRID'
      FLUSH(740+MyRankGlobal)
#endif          

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTGRID',1,ZHOOK_HANDLE)
#endif
      RETURN
      END SUBROUTINE OUTGRID
