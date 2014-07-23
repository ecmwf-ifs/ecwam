      SUBROUTINE OUTGRID (IJS, IJL, IG)
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

!        *CALL* *OUTGRID (IJS, IJL, IG)
!         *IJS*    - INDEX OF FIRST GRIDPOINT.
!         *IJL*    - INDEX OF LAST GRIDPOINT.
!         *IG*     - BLOCK NUMBER.

!     EXTERNALS.
!     ----------

!       *SETICEMASK*   - TRANSFERS A LOCAL BLOCK ARRAY TO A GLOBAL BLOCK 
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
     &            STRNMS   ,
     &            WX1      ,WX2      ,WX3      ,WX4      ,WX5
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,MPMAXLENGTH
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO     ,NIBLO
      USE YOWPCONS , ONLY : ZMISS    ,EPSUS    ,EPSU10
      USE YOWSHAL  , ONLY : DEPTH
      USE YOWSPEC  , ONLY : NSTART   ,NEND     , 
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
      USE MPL_MODULE

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: IJ,IJS,IJL,IG
      INTEGER :: ITT,ICT
      INTEGER :: NFLDTOT, NFLDPPEMAX
      INTEGER :: ITG
      INTEGER :: ICNT, IFLD
      INTEGER :: IR, IPR
      INTEGER, DIMENSION(NPROC) :: NFLDPPE 
      INTEGER, DIMENSION(NPROC) :: ISENDCOUNTS,IRECVCOUNTS

      REAL :: TAU 
      REAL :: ZHOOK_HANDLE
      REAL :: CD(IJS:IJL)
      REAL :: ZCOMBUFR(MPMAXLENGTH)
      REAL, ALLOCATABLE, DIMENSION(:) :: GTEMP 
      REAL, ALLOCATABLE, DIMENSION(:) :: WVSTR
      REAL, ALLOCATABLE, DIMENSION(:) :: ZSENDBUF,ZRECVBUF
      REAL, ALLOCATABLE, DIMENSION(:,:) :: TEMP

      LOGICAL :: LSQRT

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('OUTGRID',0,ZHOOK_HANDLE)

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
      IF(NFLDTOT.EQ.0) THEN
        IF (LHOOK) CALL DR_HOOK('OUTGRID',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
 
      ALLOCATE(TEMP(IJS:IJL,JPPFLAG))

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
          CALL SETICEMASK(EMEAN(IJS),TEMP(IJS,1),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(2).NE.0) THEN 
          CALL SETICEMASK(THQ(IJS),TEMP(IJS,2),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(3).NE.0) THEN 
          CALL SETICEMASK(FMEAN(IJS),TEMP(IJS,3),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(4).NE.0) THEN
!         The friction velocity is NOT set to missing over ice 
          CALL SETICEMASK(USNEW(IJS),TEMP(IJS,4),IJS,IJL,IG,1.1)
        ENDIF

        IF(IPFGTBL(5).NE.0) THEN 
!         The wind direction is NOT set to missing over ice 
          CALL SETICEMASK(THWNEW(IJS),TEMP(IJS,5),IJS,IJL,IG,1.1)
        ENDIF

        IF(IPFGTBL(6).NE.0) THEN
          CALL SETICEMASK(FPMEAN,TEMP(IJS,6),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(7).NE.0) THEN 
!         The drag coefficient is NOT set to missing over ice 
          CALL SETICEMASK(CD,TEMP(IJS,7),IJS,IJL,IG,1.1)
        ENDIF

        IF(IPFGTBL(8).NE.0) THEN
!         The wave stress is NOT set to missing over ice 
          CALL SETICEMASK(WVSTR,TEMP(IJS,8),IJS,IJL,IG,1.1)
          DEALLOCATE(WVSTR)
        ENDIF

        IF(IPFGTBL(9).NE.0) THEN
          CALL SETICEMASK(SMEAN(IJS),TEMP(IJS,9),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(10).NE.0) THEN
!         The wind speed is NOT set to missing over ice. 
          CALL SETICEMASK(U10NEW(IJS),TEMP(IJS,10),IJS,IJL,IG,1.1)
        ENDIF

        IF(IPFGTBL(11).NE.0) THEN
          CALL SETICEMASK(ESEA,TEMP(IJS,11),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(12).NE.0) THEN
          CALL SETICEMASK(ESWELL,TEMP(IJS,12),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(13).NE.0) THEN
          CALL SETICEMASK(THWISEA,TEMP(IJS,13),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(14).NE.0) THEN
          CALL SETICEMASK(THSWELL,TEMP(IJS,14),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(15).NE.0) THEN
          CALL SETICEMASK(FSEA,TEMP(IJS,15),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(16).NE.0) THEN
          CALL SETICEMASK(FSWELL,TEMP(IJS,16),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(22).NE.0) THEN
          CALL SETICEMASK(ALTWH,TEMP(IJS,22),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(23).NE.0) THEN
          CALL SETICEMASK(CALTWH,TEMP(IJS,23),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(24).NE.0) THEN
          CALL SETICEMASK(RALTCOR,TEMP(IJS,24),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(25).NE.0) THEN
          CALL SETICEMASK(P1MEAN,TEMP(IJS,25),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(26).NE.0) THEN
          CALL SETICEMASK(P2MEAN,TEMP(IJS,26),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(27).NE.0) THEN
          CALL SETICEMASK(SPRDMEAN,TEMP(IJS,27),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(28).NE.0) THEN
          CALL SETICEMASK(P1SEA,TEMP(IJS,28),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(29).NE.0) THEN
          CALL SETICEMASK(P1SWELL,TEMP(IJS,29),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(30).NE.0) THEN
          CALL SETICEMASK(P2SEA,TEMP(IJS,30),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(31).NE.0) THEN
          CALL SETICEMASK(P2SWELL,TEMP(IJS,31),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(32).NE.0) THEN
          CALL SETICEMASK(SPRDSEA,TEMP(IJS,32),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(33).NE.0) THEN
          CALL SETICEMASK(SPRDSWELL,TEMP(IJS,33),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(34).NE.0) THEN
          CALL SETICEMASK(C4MEAN,TEMP(IJS,34),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(35).NE.0) THEN
          CALL SETICEMASK(BFMEAN,TEMP(IJS,35),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(36).NE.0) THEN
          CALL SETICEMASK(QPMEAN,TEMP(IJS,36),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(37).NE.0) THEN
!         The bathymetry is NOT set to missing over ice. 
          CALL SETICEMASK(DEPTH(IJS,IG),TEMP(IJS,37),IJS,IJL,IG,1.1)
        ENDIF

        IF(IPFGTBL(38).NE.0) THEN
          CALL SETICEMASK(HMAXMEAN,TEMP(IJS,38),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(39).NE.0) THEN
          CALL SETICEMASK(TMAXMEAN,TEMP(IJS,39),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(40).NE.0) THEN
          CALL SETICEMASK(USTMEAN,TEMP(IJS,40),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(41).NE.0) THEN
          CALL SETICEMASK(VSTMEAN,TEMP(IJS,41),IJS,IJL,IG,CITHRSH)
        ENDIF

        IF(IPFGTBL(42).NE.0) THEN
!         the ice mask does not apply to the imposed ocean current 
          CALL SETICEMASK(U(IJS,IG),TEMP(IJS,42),IJS,IJL,IG,1.1)
        ENDIF

        IF(IPFGTBL(43).NE.0) THEN
!         the ice mask does not apply to the imposed ocean current 
          CALL SETICEMASK(V(IJS,IG),TEMP(IJS,43),IJS,IJL,IG,1.1)
        ENDIF

        IF(IPFGTBL(44).NE.0) THEN
!         The energy flux into ocean is NOT set to missing over ice 
          CALL SETICEMASK(PHIEPS(IJS),TEMP(IJS,44),IJS,IJL,IG,1.1)
        ENDIF

        IF(IPFGTBL(45).NE.0) THEN
!         The energy flux into waves is NOT set to missing over ice 
          CALL SETICEMASK(PHIAW(IJS),TEMP(IJS,45),IJS,IJL,IG,1.1)
        ENDIF

        IF(IPFGTBL(46).NE.0) THEN
!         The momentum flux into ocean is NOT set to missing over ice 
          CALL SETICEMASK(TAUOC(IJS),TEMP(IJS,46),IJS,IJL,IG,1.1)
        ENDIF

        ICT=46
        DO ITT = 1, NTRAIN
          ICT=ICT+1
          IF(IPFGTBL(ICT).NE.0) THEN
          CALL SETICEMASK(EMTRAIN(:,ITT),TEMP(IJS,ICT),
     &                    IJS,IJL,IG,CITHRSH)
          ENDIF

          ICT=ICT+1
          IF(IPFGTBL(ICT).NE.0) THEN
          CALL SETICEMASK(THTRAIN(:,ITT),TEMP(IJS,ICT),
     &                    IJS,IJL,IG,CITHRSH)
          ENDIF

          ICT=ICT+1
          IF(IPFGTBL(ICT).NE.0) THEN
          CALL SETICEMASK(PMTRAIN(:,ITT),TEMP(IJS,ICT),
     &                    IJS,IJL,IG,CITHRSH)
          ENDIF

        ENDDO


        ITG=46+3*NTRAIN+1
        IF(IPFGTBL(ITG).NE.0) THEN
!         the ice mask does not apply to the mean square strain in the ice
          CALL SETICEMASK(STRNMS(IJS),TEMP(IJS,ITG),IJS,IJL,IG,1.1)
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
!         the ice mask does not apply to air density 
          CALL SETICEMASK(ROAIRN(IJS),TEMP(IJS,ITG),IJS,IJL,IG,1.1)
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
!         the ice mask does not apply to the imposed convective velocity scale
          CALL SETICEMASK(ZIDLNEW(IJS),TEMP(IJS,ITG),IJS,IJL,IG,1.1)
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
!         the ice mask does not apply
          CALL SETICEMASK(CICOVER(IJS,IG),TEMP(IJS,ITG),IJS,IJL,IG,1.1)
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETICEMASK(CITHICK(IJS,IG),TEMP(IJS,ITG),IJS,IJL,IG,1.1)
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETICEMASK(C3MEAN,TEMP(IJS,ITG),IJS,IJL,IG,CITHRSH)
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

!       EXTRA OUTPUT FIELDS:
!       See *OUTBS* where they should have been calculated or pick them up from
!       the forcing fields (see above if the sea ice mask has to be removed)
        ITG=JPPFLAG-4
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETICEMASK(WX1(IJS),TEMP(IJS,ITG),IJS,IJL,IG,CITHRSH)
        ENDIF
        ITG=JPPFLAG-3
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETICEMASK(WX2(IJS),TEMP(IJS,ITG),IJS,IJL,IG,CITHRSH)
        ENDIF
        ITG=JPPFLAG-2
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETICEMASK(WX3(IJS),TEMP(IJS,ITG),IJS,IJL,IG,CITHRSH)
        ENDIF
        ITG=JPPFLAG-2
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETICEMASK(WX4(IJS),TEMP(IJS,ITG),IJS,IJL,IG,CITHRSH)
        ENDIF
        ITG=JPPFLAG
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETICEMASK(WX5(IJS),TEMP(IJS,ITG),IJS,IJL,IG,CITHRSH)
        ENDIF


!     SENDING TO RELEVANT PE'S
!     ------------------------


      DO IPR=1,NPROC
        ISENDCOUNTS(IPR)=NFLDPPE(IPR)*(NEND(IRANK)-NSTART(IRANK)+1)
      ENDDO

      DO IPR=1,NPROC
        IRECVCOUNTS(IPR)=NFLDPPE(IRANK)*(NEND(IPR)-NSTART(IPR)+1)
      ENDDO

!     LOADING THE COMMUNICATION BUFFER
      ALLOCATE(ZSENDBUF(NFLDTOT*MPMAXLENGTH))
      ICNT=0
      DO IPR=1,NPROC
        DO ICT=1,JPPFLAG
          IF(IPFGTBL(ICT).EQ.IPR) THEN 
            DO IJ=IJS,IJL
              ICNT = ICNT + 1
              ZSENDBUF(ICNT) = TEMP(IJ,ICT)
            ENDDO
          ENDIF
        ENDDO
      ENDDO

      DEALLOCATE(TEMP)
      ALLOCATE(ZRECVBUF(NFLDPPEMAX*NIBLO))

!     GLOBAL EXCHANGE
      CALL GSTATS(674,0)
      CALL MPL_ALLTOALLV(ZSENDBUF,ISENDCOUNTS,
     &                   ZRECVBUF,IRECVCOUNTS,
     &                   CDSTRING='OUTGRID:')
      CALL GSTATS(674,1)
      DEALLOCATE(ZSENDBUF)


!     RETRIEVE THE INFORMATION
!     ------------------------

      IF(NFLDPPE(IRANK).GT.0) THEN
        ALLOCATE(GTEMP(NIBLO))


!     RECEIVING AND TRANSFERING FROM BLOCK TO GRID
!     --------------------------------------------

      IFLD=0
      DO ICT=1,JPPFLAG
        IF(IFLD.GE.NFLDPPE(IRANK)) EXIT
        IF(IPFGTBL(ICT).EQ.IRANK) THEN
          ICNT=0
          DO IPR=1,NPROC
            ICNT=ICNT+IFLD*(NEND(IPR)-NSTART(IPR)+1)
            DO IJ=NSTART(IPR),NEND(IPR)
              ICNT=ICNT+1
              GTEMP(IJ)=ZRECVBUF(ICNT)
            ENDDO
            ICNT=ICNT+(NFLDPPE(IRANK)-1-IFLD)*(NEND(IPR)-NSTART(IPR)+1)
          ENDDO

!!! until I move all grid arrays to one work array, need to loop on all
!!! special case
          IF(ICT.EQ.1) THEN
            LSQRT=.TRUE.
            CALL MAKEGRID (GTEMP,WHGTTG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.2) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WDIRTG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.3) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WMNFTG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.4) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,USTARG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.5) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,UDIRG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.6) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WPKFTG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.7) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,CDG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.8) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,TAUWG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.9) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,SMEANG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.10) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,U10G,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.11) THEN
            LSQRT=.TRUE.
            CALL MAKEGRID (GTEMP,WHGTWG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.12) THEN
            LSQRT=.TRUE.
            CALL MAKEGRID (GTEMP,WHGTSG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.13) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WDIRWG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.14) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WDIRSG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.15) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WMNFWG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.16) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WMNFSG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.22) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WHGTAG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.23) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,CWHTAG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.24) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,RANGAG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.25) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,MWP1G,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.26) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,MWP2G,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.27) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WSPRDG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.28) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,MWP1WG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.29) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,MWP1SG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.30) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,MWP2WG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.31) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,MWP2SG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.32) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WSPRDWG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.33) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WSPRDSG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.34) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,C4G,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.35) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,BFG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.36) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,QPG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.37) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,DEPTHG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.38) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,HMAXG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.39) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,TMAXG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.40) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,USTOKESG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.41) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,VSTOKESG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.42) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,UCURG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.43) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,VCURG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.44) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,PHIEPSG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.45) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,PHIAWG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.46) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,TAUOCG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ENDIF

          DO ITT = 1, NTRAIN
            IF(ICT.EQ.(46+(ITT-1)*3+1) ) THEN
              LSQRT=.TRUE.
              CALL MAKEGRID (GTEMP,WHGTTRG(:,:,ITT),
     &                   1,NEND(NPROC),IG,LSQRT,ZMISS)
              IFLD=IFLD+1

            ELSEIF(ICT.EQ.(46+(ITT-1)*3+2) ) THEN
              LSQRT=.FALSE.
              CALL MAKEGRID (GTEMP,WDIRTRG(:,:,ITT),
     &                   1,NEND(NPROC),IG,LSQRT,ZMISS)
              IFLD=IFLD+1

            ELSEIF(ICT.EQ.(46+(ITT-1)*3+3) ) THEN
              LSQRT=.FALSE.
              CALL MAKEGRID (GTEMP,MWPMTRG(:,:,ITT),
     &                   1,NEND(NPROC),IG,LSQRT,ZMISS)
              IFLD=IFLD+1
            ENDIF
          ENDDO

          IF(ICT.EQ.(46+3*NTRAIN+1)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,STRNMSG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+2)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,RHOAG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+3)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,WSTARG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+4)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,CICG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+5)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,CIHG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+6)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,C3G,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+7)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,NSSTG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+8)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,NICCG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+9)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,NICTG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+10)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,NUCUG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ELSEIF(ICT.EQ.(46+3*NTRAIN+11)) THEN
            LSQRT=.FALSE.
            CALL MAKEGRID (GTEMP,NVCUG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IFLD=IFLD+1
          ENDIF

!     EXTRA OUTPUT FIELDS:
!     See *OUTBS* where they should have been calculated or pick them up from
!     the forcing fields (see above if the sea ice mask has to be removed)
        IF(ICT.EQ.(JPPFLAG-4)) THEN
          LSQRT=.FALSE.
          CALL MAKEGRID (GTEMP,WX1G,1,NEND(NPROC),IG,LSQRT,0.)
          IFLD=IFLD+1
        ELSEIF(ICT.EQ.(JPPFLAG-3)) THEN
          LSQRT=.FALSE.
          CALL MAKEGRID (GTEMP,WX2G,1,NEND(NPROC),IG,LSQRT,0.)
          IFLD=IFLD+1
        ELSEIF(ICT.EQ.(JPPFLAG-2)) THEN
          LSQRT=.FALSE.
          CALL MAKEGRID (GTEMP,WX3G,1,NEND(NPROC),IG,LSQRT,0.)
          IFLD=IFLD+1
        ELSEIF(ICT.EQ.(JPPFLAG-1)) THEN
          LSQRT=.FALSE.
          CALL MAKEGRID (GTEMP,WX4G,1,NEND(NPROC),IG,LSQRT,0.)
          IFLD=IFLD+1
        ELSEIF(ICT.EQ.(JPPFLAG)) THEN
          LSQRT=.FALSE.
          CALL MAKEGRID (GTEMP,WX5G,1,NEND(NPROC),IG,LSQRT,0.)
          IFLD=IFLD+1
        ENDIF


        ENDIF ! IPFGTBL
      ENDDO  ! ICT

        DEALLOCATE(GTEMP)
      ENDIF

      DEALLOCATE(ZRECVBUF)


      IF (LHOOK) CALL DR_HOOK('OUTGRID',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE OUTGRID
