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

!       *LOC2GLO*   - TRANSFERS A LOCAL BLOCK ARRAY TO A GLOBAL BLOCK 
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

      USE YOWCOUT  , ONLY : JPPFLAG  ,IPFGTBL
      USE YOWCURR  , ONLY : U, V 
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE ,CICOVER  ,CITHRSH  ,
     &            CICOVER  ,CITHICK 
      USE YOWINTP  , ONLY : WHGTTG   ,WDIRTG   ,WPKFTG   ,WMNFTG   ,
     &            USTARG   ,UDIRG    ,TAUWG    ,CDG      ,SMEANG   ,
     &            U10G     ,MWP1G    ,MWP2G    ,WSPRDG   ,C4G      ,
     &            BFG      ,QPG      ,DEPTHG   ,HMAXG    ,TMAXG    ,
     &            USTOKESG ,VSTOKESG ,UCURG    ,VCURG    ,PHIEPSG  ,
     &            PHIAWG   ,TAUOCG   ,STRNMSG  ,RHOAG    ,WSTARG   ,
     &            CICG     ,CIHG     ,
     &            WX1G     ,WX2G     ,WX3G     ,WX4G     ,WX5G
      USE YOWINTS  , ONLY : WHGTSG   ,WDIRSG   ,WMNFSG   ,WHGTWG   ,
     &            WDIRWG   ,WMNFWG   ,MWP1SG   ,MWP2SG   ,WSPRDSG  ,
     &            MWP1WG   ,MWP2WG   ,WSPRDWG
      USE YOWINTT  , ONLY : WHGTT1G  ,WDIRT1G  ,MWP1T1G  ,
     &                      WHGTT2G  ,WDIRT2G  ,MWP1T2G  ,
     &                      WHGTT3G  ,WDIRT3G  ,MWP1T3G
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,THQ      ,
     &            SMEAN    ,P1MEAN   ,P2MEAN   ,SPRDMEAN ,C4MEAN   ,
     &            BFMEAN   ,QPMEAN   ,HMAXMEAN ,FPMEAN   ,TMAXMEAN ,
     &            USTMEAN  ,VSTMEAN  ,PHIEPS   ,PHIAW    ,TAUOC    ,
     &            STRNMS   ,
     &            WX1      ,WX2      ,WX3      ,WX4      ,WX5
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
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
      USE YOWTRAINS , ONLY: ETRAIN1  ,THTRAIN1 ,P1TRAIN1 ,
     &                      ETRAIN2  ,THTRAIN2 ,P1TRAIN2 ,
     &                      ETRAIN3  ,THTRAIN3 ,P1TRAIN3
      USE YOWTEST  , ONLY : IU06     
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE MPL_MODULE

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER ::  IJ,IJS,IJL,IG
      INTEGER :: ITG

      REAL :: TAU 
      REAL :: ZHOOK_HANDLE
      REAL,ALLOCATABLE,DIMENSION(:) :: CD,WVSTR
      REAL :: TEMP(NIBLO)

      LOGICAL :: LSQRT

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('OUTGRID',0,ZHOOK_HANDLE)

!*    1. INTEGRATED PARAMETERS OF TOTAL SEA.
!        -----------------------------------

        ALLOCATE (CD(IJS:IJL))
        IF(IPFGTBL(8).NE.0) ALLOCATE(WVSTR(IJS:IJL))
        DO IJ = IJS,IJL
          TAU = MAX(USNEW(IJ)**2,EPSUS)
          IF(IPFGTBL(8).NE.0) WVSTR(IJ) = TAUW(IJ,IG)/TAU
!!       if the numerical computation of TAU and CD change, a similar
!!       modification has to be put in buildstress where the friction
!!       velocity is determined from U10 and CD.
          CD(IJ) = TAU/MAX(U10NEW(IJ)**2,EPSU10)
        ENDDO

!*    1.3 STORE PARAMETERS FOR ACTUAL BLOCK IN GRID ARRAYS.
!         -------------------------------------------------

          IF(IPFGTBL(1).NE.0) THEN 
            CALL LOC2GLO(EMEAN(IJS),TEMP,IJS,IJL,IG)
            LSQRT=.TRUE.
            CALL MPGATHERSCFLD(IPFGTBL(1),1,NSTART,NEND,TEMP)
            IF(IPFGTBL(1).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,WHGTTG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(2).NE.0) THEN 
            CALL LOC2GLO(THQ(IJS),TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(2),2,NSTART,NEND,TEMP)
            IF(IPFGTBL(2).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,WDIRTG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(3).NE.0) THEN 
            CALL LOC2GLO(FMEAN(IJS),TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(3),3,NSTART,NEND,TEMP)
            IF(IPFGTBL(3).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,WMNFTG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(4).NE.0) THEN
            DO IJ = IJS,IJL
              TEMP(IJ) = USNEW(IJ)
!             cannot set USNEW to ZMISS since it can still be used
!             in the next time step
              IF (LICERUN.AND.LMASKICE.AND.
     &            CICOVER(IJ,IG).GT.CITHRSH) THEN
                TEMP(IJ)=ZMISS
              ELSE
                TEMP(IJ) = USNEW(IJ)
              ENDIF
            ENDDO
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(4),4,NSTART,NEND,TEMP)
            IF(IPFGTBL(4).EQ.IRANK) 
     &        CALL MAKEGRID (TEMP,USTARG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(5).NE.0) THEN 
            DO IJ = IJS,IJL
!             The wind direction is NOT set to missing over ice 
              TEMP(IJ) = THWNEW(IJ) 
            ENDDO
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(5),5,NSTART,NEND,TEMP)
            IF(IPFGTBL(5).EQ.IRANK) 
     &        CALL MAKEGRID (TEMP,UDIRG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(6).NE.0) THEN
            CALL LOC2GLO(FPMEAN,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(6),6,NSTART,NEND,TEMP)
            IF(IPFGTBL(6).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,WPKFTG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(7).NE.0) THEN 
            CALL LOC2GLO(CD,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(7),7,NSTART,NEND,TEMP)
            IF(IPFGTBL(7).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,CDG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            IF(ALLOCATED(CD)) DEALLOCATE(CD)
          ENDIF

          IF(IPFGTBL(8).NE.0) THEN
            CALL LOC2GLO(WVSTR,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(8),8,NSTART,NEND,TEMP)
            IF(IPFGTBL(8).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,TAUWG,1,NEND(NPROC),IG,LSQRT,ZMISS)
            DEALLOCATE(WVSTR)
          ENDIF

          IF(IPFGTBL(9).NE.0) THEN
            CALL LOC2GLO(SMEAN(IJS),TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(9),9,NSTART,NEND,TEMP)
            IF(IPFGTBL(9).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,SMEANG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(10).NE.0) THEN
            DO IJ = IJS,IJL
!             The wind speed is NOT set to missing over ice. 
              TEMP(IJ) = U10NEW(IJ)
            ENDDO
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(10),10,NSTART,NEND,TEMP)
            IF(IPFGTBL(10).EQ.IRANK) 
     &        CALL MAKEGRID (TEMP,U10G,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(25).NE.0) THEN
            CALL LOC2GLO(P1MEAN,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(25),25,NSTART,NEND,TEMP)
            IF(IPFGTBL(25).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,MWP1G,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(26).NE.0) THEN
            CALL LOC2GLO(P2MEAN,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(26),26,NSTART,NEND,TEMP)
            IF(IPFGTBL(26).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,MWP2G,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(27).NE.0) THEN
            CALL LOC2GLO(SPRDMEAN,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(27),27,NSTART,NEND,TEMP)
            IF(IPFGTBL(27).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,WSPRDG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(34).NE.0) THEN
            CALL LOC2GLO(C4MEAN,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(34),34,NSTART,NEND,TEMP)
            IF(IPFGTBL(34).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,C4G,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(35).NE.0) THEN
            CALL LOC2GLO(BFMEAN,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(35),35,NSTART,NEND,TEMP)
            IF(IPFGTBL(35).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,BFG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(36).NE.0) THEN
            CALL LOC2GLO(QPMEAN,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(36),36,NSTART,NEND,TEMP)
            IF(IPFGTBL(36).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,QPG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(37).NE.0) THEN
            DO IJ = IJS,IJL
!             The bathymetry is NOT set to missing over ice. 
              TEMP(IJ) = DEPTH(IJ,IG)
            ENDDO
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(37),37,NSTART,NEND,TEMP)
            IF(IPFGTBL(37).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,DEPTHG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(38).NE.0) THEN
            CALL LOC2GLO(HMAXMEAN,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(38),38,NSTART,NEND,TEMP)
            IF(IPFGTBL(38).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,HMAXG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(39).NE.0) THEN
            CALL LOC2GLO(TMAXMEAN,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(39),39,NSTART,NEND,TEMP)
            IF(IPFGTBL(39).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,TMAXG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(40).NE.0) THEN
            CALL LOC2GLO(USTMEAN,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(40),40,NSTART,NEND,TEMP)
            IF(IPFGTBL(40).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,USTOKESG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(41).NE.0) THEN
            CALL LOC2GLO(VSTMEAN,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(41),41,NSTART,NEND,TEMP)
            IF(IPFGTBL(41).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,VSTOKESG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(42).NE.0) THEN
            DO IJ = IJS,IJL
!             the ice mask does not apply to the imposed ocean current 
              TEMP(IJ) = U(IJ,IG)
            ENDDO
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(42),42,NSTART,NEND,TEMP)
            IF(IPFGTBL(42).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,UCURG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(43).NE.0) THEN
            DO IJ = IJS,IJL
!             the ice mask does not apply to the imposed ocean current 
              TEMP(IJ) = V(IJ,IG)
            ENDDO
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(43),43,NSTART,NEND,TEMP)
            IF(IPFGTBL(43).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,VCURG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(44).NE.0) THEN
            CALL LOC2GLO(PHIEPS(IJS),TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(44),44,NSTART,NEND,TEMP)
            IF(IPFGTBL(44).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,PHIEPSG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(45).NE.0) THEN
            CALL LOC2GLO(PHIAW(IJS),TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(45),45,NSTART,NEND,TEMP)
            IF(IPFGTBL(45).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,PHIAWG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(46).NE.0) THEN
            CALL LOC2GLO(TAUOC(IJS),TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(46),46,NSTART,NEND,TEMP)
            IF(IPFGTBL(46).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,TAUOCG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          ITG=56
          IF(IPFGTBL(ITG).NE.0) THEN
            DO IJ = IJS,IJL
!             the ice mask does not apply to the mean square strain in the ice
              TEMP(IJ) = STRNMS(IJ)
            ENDDO
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(ITG),ITG,NSTART,NEND,TEMP)
            IF(IPFGTBL(ITG).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,STRNMSG,1,NEND(NPROC),IG,LSQRT,0.)
          ENDIF

          ITG=57
          IF(IPFGTBL(ITG).NE.0) THEN
            DO IJ = IJS,IJL
!             the ice mask does not apply to the imposed air density
              TEMP(IJ) = ROAIRN(IJ)
            ENDDO
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(ITG),ITG,NSTART,NEND,TEMP)
            IF(IPFGTBL(ITG).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,RHOAG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          ITG=58
          IF(IPFGTBL(ITG).NE.0) THEN
            DO IJ = IJS,IJL
!             the ice mask does not apply to the imposed convective velocity scale
              TEMP(IJ) = ZIDLNEW(IJ)
            ENDDO
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(ITG),ITG,NSTART,NEND,TEMP)
            IF(IPFGTBL(ITG).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,WSTARG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          ITG=59
          IF(IPFGTBL(ITG).NE.0) THEN
            DO IJ = IJS,IJL
!             the ice mask does not apply
              TEMP(IJ) = CICOVER(IJ,IG)
            ENDDO
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(ITG),ITG,NSTART,NEND,TEMP)
            IF(IPFGTBL(ITG).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,CICG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          ITG=60
          IF(IPFGTBL(ITG).NE.0) THEN
            DO IJ = IJS,IJL
!             the ice mask does not apply
              TEMP(IJ) = CITHICK(IJ,IG)
            ENDDO
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(ITG),ITG,NSTART,NEND,TEMP)
            IF(IPFGTBL(ITG).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,CIHG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

!         EXTRA OUTPUT FIELDS:
!         See *OUTBS* where they should have been calculated or pick them up from
!         the forcing fields (see above if the sea ice mask has to be removed)
          ITG=JPPFLAG-4
          IF(IPFGTBL(ITG).NE.0) THEN
            CALL LOC2GLO(WX1(IJS),TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(ITG),ITG,NSTART,NEND,TEMP)
            IF(IPFGTBL(ITG).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,WX1G,1,NEND(NPROC),IG,LSQRT,0.)
          ENDIF
          ITG=JPPFLAG-3
          IF(IPFGTBL(ITG).NE.0) THEN
            CALL LOC2GLO(WX2(IJS),TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(ITG),ITG,NSTART,NEND,TEMP)
            IF(IPFGTBL(ITG).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,WX2G,1,NEND(NPROC),IG,LSQRT,0.)
          ENDIF
          ITG=JPPFLAG-2
          IF(IPFGTBL(ITG).NE.0) THEN
            CALL LOC2GLO(WX3(IJS),TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(ITG),ITG,NSTART,NEND,TEMP)
            IF(IPFGTBL(ITG).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,WX3G,1,NEND(NPROC),IG,LSQRT,0.)
          ENDIF
          ITG=JPPFLAG-2
          IF(IPFGTBL(ITG).NE.0) THEN
            CALL LOC2GLO(WX4(IJS),TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(ITG),ITG,NSTART,NEND,TEMP)
            IF(IPFGTBL(ITG).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,WX4G,1,NEND(NPROC),IG,LSQRT,0.)
          ENDIF
          ITG=JPPFLAG
          IF(IPFGTBL(ITG).NE.0) THEN
            CALL LOC2GLO(WX5(IJS),TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(ITG),ITG,NSTART,NEND,TEMP)
            IF(IPFGTBL(ITG).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,WX5G,1,NEND(NPROC),IG,LSQRT,0.)
          ENDIF


! ----------------------------------------------------------------------

!*    2. INTEGRATED PARAMETERS OF SWELL AND WIND SEA.
!        --------------------------------------------

      IF (CDTINTS.EQ.CDTPRO) THEN

!*    2.1 STORE PARAMETERS FOR ACTUAL BLOCK IN GRID ARRAYS.
!         -------------------------------------------------

          IF(IPFGTBL(11).NE.0) THEN
            CALL LOC2GLO(ESEA,TEMP,IJS,IJL,IG)
            LSQRT=.TRUE.
            CALL MPGATHERSCFLD(IPFGTBL(11),11,NSTART,NEND,TEMP)
            IF(IPFGTBL(11).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,WHGTWG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(12).NE.0) THEN
            CALL LOC2GLO(ESWELL,TEMP,IJS,IJL,IG)
            LSQRT=.TRUE.
            CALL MPGATHERSCFLD(IPFGTBL(12),12,NSTART,NEND,TEMP)
            IF(IPFGTBL(12).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,WHGTSG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(13).NE.0) THEN
            CALL LOC2GLO(THWISEA,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(13),13,NSTART,NEND,TEMP)
            IF(IPFGTBL(13).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,WDIRWG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(14).NE.0) THEN
            CALL LOC2GLO(THSWELL,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(14),14,NSTART,NEND,TEMP)
            IF(IPFGTBL(14).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,WDIRSG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(15).NE.0) THEN
            CALL LOC2GLO(FSEA,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(15),15,NSTART,NEND,TEMP)
            IF(IPFGTBL(15).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,WMNFWG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(16).NE.0) THEN
            CALL LOC2GLO(FSWELL,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(16),16,NSTART,NEND,TEMP)
            IF(IPFGTBL(16).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,WMNFSG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(28).NE.0) THEN
            CALL LOC2GLO(P1SEA,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(28),28,NSTART,NEND,TEMP)
            IF(IPFGTBL(28).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,MWP1WG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(29).NE.0) THEN
            CALL LOC2GLO(P1SWELL,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(29),29,NSTART,NEND,TEMP)
            IF(IPFGTBL(29).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,MWP1SG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(30).NE.0) THEN
            CALL LOC2GLO(P2SEA,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(30),30,NSTART,NEND,TEMP)
            IF(IPFGTBL(30).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,MWP2WG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(31).NE.0) THEN
            CALL LOC2GLO(P2SWELL,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(31),31,NSTART,NEND,TEMP)
            IF(IPFGTBL(31).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,MWP2SG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(32).NE.0) THEN
            CALL LOC2GLO(SPRDSEA,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(32),32,NSTART,NEND,TEMP)
            IF(IPFGTBL(32).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,WSPRDWG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(33).NE.0) THEN
            CALL LOC2GLO(SPRDSWELL,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(33),33,NSTART,NEND,TEMP)
            IF(IPFGTBL(33).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,WSPRDSG,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

      ENDIF

! ----------------------------------------------------------------------

!*    3. INTEGRATED PARAMETERS OF 3 WAVE TRAINS
!        --------------------------------------------

      IF (CDTINTS.EQ.CDTPRO) THEN

!*    3.1 STORE PARAMETERS FOR ACTUAL BLOCK IN GRID ARRAYS.
!         -------------------------------------------------

          IF(IPFGTBL(47).NE.0) THEN
            CALL LOC2GLO(ETRAIN1,TEMP,IJS,IJL,IG)
            LSQRT=.TRUE.
            CALL MPGATHERSCFLD(IPFGTBL(47),47,NSTART,NEND,TEMP)
            IF(IPFGTBL(47).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,WHGTT1G,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(48).NE.0) THEN
            CALL LOC2GLO(THTRAIN1,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(48),48,NSTART,NEND,TEMP)
            IF(IPFGTBL(48).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,WDIRT1G,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(49).NE.0) THEN
            CALL LOC2GLO(P1TRAIN1,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(49),49,NSTART,NEND,TEMP)
            IF(IPFGTBL(49).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,MWP1T1G,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(50).NE.0) THEN
            CALL LOC2GLO(ETRAIN2,TEMP,IJS,IJL,IG)
            LSQRT=.TRUE.
            CALL MPGATHERSCFLD(IPFGTBL(50),50,NSTART,NEND,TEMP)
            IF(IPFGTBL(50).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,WHGTT2G,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(51).NE.0) THEN
            CALL LOC2GLO(THTRAIN2,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(51),51,NSTART,NEND,TEMP)
            IF(IPFGTBL(51).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,WDIRT2G,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(52).NE.0) THEN
            CALL LOC2GLO(P1TRAIN2,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(52),52,NSTART,NEND,TEMP)
            IF(IPFGTBL(52).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,MWP1T2G,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(53).NE.0) THEN
            CALL LOC2GLO(ETRAIN3,TEMP,IJS,IJL,IG)
            LSQRT=.TRUE.
            CALL MPGATHERSCFLD(IPFGTBL(53),53,NSTART,NEND,TEMP)
            IF(IPFGTBL(53).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,WHGTT3G,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(54).NE.0) THEN
            CALL LOC2GLO(THTRAIN3,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(54),54,NSTART,NEND,TEMP)
            IF(IPFGTBL(54).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,WDIRT3G,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

          IF(IPFGTBL(55).NE.0) THEN
            CALL LOC2GLO(P1TRAIN3,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(55),55,NSTART,NEND,TEMP)
            IF(IPFGTBL(55).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,MWP1T3G,1,NEND(NPROC),IG,LSQRT,ZMISS)
          ENDIF

      ENDIF

      CALL GSTATS(750,0)
      CALL MPL_BARRIER(CDSTRING='OUTGRID:')
      CALL GSTATS(750,1)

      IF (LHOOK) CALL DR_HOOK('OUTGRID',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE OUTGRID
