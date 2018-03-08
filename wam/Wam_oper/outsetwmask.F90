      SUBROUTINE OUTSETWMASK (IJS, IJL, IODP, BOUT)
! ----------------------------------------------------------------------

!**** *OUTSETWMASK* -


!*    PURPOSE.
!     --------

!       OUTPUT OF WAVE AND WIND FIELDS INTO COMMON BLOCK.


!**   INTERFACE.
!     ----------

!        *CALL* *OUTSETWMASK (IJS, IJL, IODP, BOUT)
!         *IJS*    - INDEX OF FIRST LOCAL GRIDPOINT.
!         *IJL*    - INDEX OF LAST LOCAL GRIDPOINT.
!         *IODP*     LAND MASK IF IODP(IJ)=0
!         *BOUT*


!     EXTERNALS.
!     ----------

!       *SETWMASK*   - TRANSFERS A LOCAL BLOCK ARRAY TO A GLOBAL BLOCK 
!                      ARRAY AND IMPOSES ICE MASK.

!     METHOD.
!     -------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : NTRAIN   ,JPPFLAG  ,IPFGTBL,IFRSTPARTI
      USE YOWCURR  , ONLY : U, V 
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE ,CICOVER  ,CITHRSH  ,
     &            CIBLOCK ,CITHICK 
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,THQ      ,CDMEAN   ,
     &            WVSTRMEAN,SMEAN    ,ALTWH    ,CALTWH   ,RALTCOR  ,
     &            P1MEAN   ,P2MEAN   ,SPRDMEAN ,C3MEAN   ,C4MEAN   ,
     &            BFMEAN   ,QPMEAN   ,HMAXMEAN ,FPMEAN   ,TMAXMEAN ,
     &            USTMEAN  ,VSTMEAN  ,PHIEPS   ,PHIAW    ,TAUOC    ,
     &            STRNMS   ,E10MEAN  ,WEFLXM   ,WEFLXD   ,E1012    ,
     &            E1214    ,E1417    ,E1721    ,E2125    ,E2530    ,
     &            WX1      ,WX2      ,WX3      ,WX4      ,WX5
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO     ,NIBLO
      USE YOWPCONS , ONLY : ZMISS    ,EPSUS    ,EPSU10
      USE YOWSHAL  , ONLY : DEPTH
      USE YOWSPEC  , ONLY :  U10NEW   ,THWNEW   ,USNEW    ,
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

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL), INTENT(IN) :: IODP
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,JPPFLAG), INTENT(INOUT) :: BOUT

      INTEGER(KIND=JWIM), PARAMETER :: IG=1
      INTEGER(KIND=JWIM) :: IJ, ITT, ICT, ITG, IFLD, IR
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL) :: NIODP

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

!----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTSETWMASK',0,ZHOOK_HANDLE)
#endif


      DO IR=1,JPPFLAG
         ITG=ITOBOUT(IR)
         IF(ITG.GT.0) THEN
         CALL SETWMASK(BOUT(IJS,ITG),IJS,IJL,CITHRSH,IODP(IJS))
         ENDIF
      ENDDO


      NIODP(:)=1

!*    1. INTEGRATED PARAMETERS OF TOTAL SEA.
!        -----------------------------------

!*    1.3 APPLY SEA ICE MASK (IF NEEDED)
!         ------------------------------

        IF(IPFGTBL(1).NE.0) THEN 
         CALL SETWMASK(EMEAN(IJS),BOUT(IJS,1),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(2).NE.0) THEN 
          CALL SETWMASK(THQ(IJS),BOUT(IJS,2),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(3).NE.0) THEN 
         CALL SETWMASK(FMEAN(IJS),BOUT(IJS,3),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(4).NE.0) THEN
!         The friction velocity is NOT set to missing over ice 
          CALL SETWMASK(USNEW(IJS),BOUT(IJS,4),IJS,IJL,1.1,IODP(IJS))
        ENDIF

        IF(IPFGTBL(5).NE.0) THEN 
!         The wind direction is NOT set to missing over ice or land 
          CALL SETWMASK(THWNEW(IJS),BOUT(IJS,5),IJS,IJL,1.1,NIODP(IJS))
        ENDIF

        IF(IPFGTBL(6).NE.0) THEN
          CALL SETWMASK(FPMEAN,BOUT(IJS,6),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(7).NE.0) THEN 
!         The drag coefficient is NOT set to missing over ice or land 
          CALL SETWMASK(CDMEAN,BOUT(IJS,7),IJS,IJL,1.1,NIODP(IJS))
        ENDIF

        IF(IPFGTBL(8).NE.0) THEN
!         The wave stress is NOT set to missing over ice 
          CALL SETWMASK(WVSTRMEAN,BOUT(IJS,8),IJS,IJL,1.1,IODP(IJS))
        ENDIF

        IF(IPFGTBL(9).NE.0) THEN
         CALL SETWMASK(SMEAN(IJS),BOUT(IJS,9),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(10).NE.0) THEN
!         The wind speed is NOT set to missing over ice or land 
          CALL SETWMASK(U10NEW(IJS),BOUT(IJS,10),IJS,IJL,1.1,NIODP(IJS))
        ENDIF

        IF(IPFGTBL(11).NE.0) THEN
          CALL SETWMASK(ESEA,BOUT(IJS,11),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(12).NE.0) THEN
          CALL SETWMASK(ESWELL,BOUT(IJS,12),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(13).NE.0) THEN
          CALL SETWMASK(THWISEA,BOUT(IJS,13),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(14).NE.0) THEN
          CALL SETWMASK(THSWELL,BOUT(IJS,14),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(15).NE.0) THEN
          CALL SETWMASK(FSEA,BOUT(IJS,15),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(16).NE.0) THEN
          CALL SETWMASK(FSWELL,BOUT(IJS,16),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(17).NE.0) THEN
          CALL SETWMASK(ALTWH,BOUT(IJS,22),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(18).NE.0) THEN
          CALL SETWMASK(CALTWH,BOUT(IJS,23),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(19).NE.0) THEN
          CALL SETWMASK(RALTCOR,BOUT(IJS,24),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(20).NE.0) THEN
          CALL SETWMASK(P1MEAN,BOUT(IJS,25),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(21).NE.0) THEN
          CALL SETWMASK(P2MEAN,BOUT(IJS,26),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(22).NE.0) THEN
          CALL SETWMASK(SPRDMEAN,BOUT(IJS,27),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(23).NE.0) THEN
          CALL SETWMASK(P1SEA,BOUT(IJS,28),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(24).NE.0) THEN
          CALL SETWMASK(P1SWELL,BOUT(IJS,29),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(25).NE.0) THEN
          CALL SETWMASK(P2SEA,BOUT(IJS,30),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(26).NE.0) THEN
          CALL SETWMASK(P2SWELL,BOUT(IJS,31),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(27).NE.0) THEN
          CALL SETWMASK(SPRDSEA,BOUT(IJS,32),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(28).NE.0) THEN
         CALL SETWMASK(SPRDSWELL,BOUT(IJS,33),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(29).NE.0) THEN
          CALL SETWMASK(C4MEAN,BOUT(IJS,34),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(30).NE.0) THEN
          CALL SETWMASK(BFMEAN,BOUT(IJS,35),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(31).NE.0) THEN
          CALL SETWMASK(QPMEAN,BOUT(IJS,36),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(32).NE.0) THEN
!         The bathymetry is NOT set to missing over ice. 
         CALL SETWMASK(DEPTH(IJS,IG),BOUT(IJS,37),IJS,IJL,1.1,IODP(IJS))
        ENDIF

        IF(IPFGTBL(33).NE.0) THEN
          CALL SETWMASK(HMAXMEAN,BOUT(IJS,38),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(34).NE.0) THEN
          CALL SETWMASK(TMAXMEAN,BOUT(IJS,39),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(35).NE.0) THEN
          CALL SETWMASK(USTMEAN,BOUT(IJS,40),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(36).NE.0) THEN
          CALL SETWMASK(VSTMEAN,BOUT(IJS,41),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        IF(IPFGTBL(37).NE.0) THEN
!         the ice mask does not apply to the imposed ocean current 
          IF(ALLOCATED(U)) THEN
            CALL SETWMASK(U(IJS,IG),BOUT(IJS,42),IJS,IJL,1.1,IODP(IJS))
          ELSE
            BOUT(:,37)=0.0_JWRB
          ENDIF
        ENDIF

        IF(IPFGTBL(38).NE.0) THEN
!         the ice mask does not apply to the imposed ocean current 
          IF(ALLOCATED(V)) THEN
            CALL SETWMASK(V(IJS,IG),BOUT(IJS,43),IJS,IJL,1.1,IODP(IJS))
          ELSE
            BOUT(:,38)=0.0_JWRB
          ENDIF
        ENDIF

        IF(IPFGTBL(39).NE.0) THEN
!         The energy flux into ocean is NOT set to missing over ice 
          CALL SETWMASK(PHIEPS(IJS),BOUT(IJS,44),IJS,IJL,1.1,IODP(IJS))
        ENDIF

        IF(IPFGTBL(40).NE.0) THEN
!         The energy flux into waves is NOT set to missing over ice 
          CALL SETWMASK(PHIAW(IJS),BOUT(IJS,45),IJS,IJL,1.1,IODP(IJS))
        ENDIF

        IF(IPFGTBL(41).NE.0) THEN
!         The momentum flux into ocean is NOT set to missing over ice 
          CALL SETWMASK(TAUOC(IJS),BOUT(IJS,46),IJS,IJL,1.1,IODP(IJS))
        ENDIF

        ICT=IFRSTPARTI-1
        DO ITT = 1, NTRAIN
          ICT=ICT+1
          IF(IPFGTBL(ICT).NE.0) THEN
            CALL SETWMASK(EMTRAIN(:,ITT),BOUT(IJS,ICT),IJS,IJL,CITHRSH,IODP(IJS))
          ENDIF

          ICT=ICT+1
          IF(IPFGTBL(ICT).NE.0) THEN
            CALL SETWMASK(THTRAIN(:,ITT),BOUT(IJS,ICT),IJS,IJL,CITHRSH,IODP(IJS))
          ENDIF

          ICT=ICT+1
          IF(IPFGTBL(ICT).NE.0) THEN
            CALL SETWMASK(PMTRAIN(:,ITT),BOUT(IJS,ICT),IJS,IJL,CITHRSH,IODP(IJS))
          ENDIF

        ENDDO


        ITG=IFRSTPARTI-1+3*NTRAIN+1
        IF(IPFGTBL(ITG).NE.0) THEN
!         the ice mask does not apply to the mean square strain in the ice
          CALL SETWMASK(STRNMS(IJS),BOUT(IJS,ITG),IJS,IJL,1.1,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETWMASK(E10MEAN(IJS),BOUT(IJS,ITG),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
!         the ice mask does not apply to air density or land 
          CALL SETWMASK(ROAIRN(IJS),BOUT(IJS,ITG),IJS,IJL,1.1,NIODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
!         the ice mask does not apply to the imposed convective velocity scale
          CALL SETWMASK(ZIDLNEW(IJS),BOUT(IJS,ITG),IJS,IJL,1.1,NIODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
!         the ice mask does not apply
          CALL SETWMASK(CICOVER(IJS,IG),BOUT(IJS,ITG),IJS,IJL,1.1,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETWMASK(CITHICK(IJS,IG),BOUT(IJS,ITG),IJS,IJL,1.1,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETWMASK(C3MEAN,BOUT(IJS,ITG),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          DO IJ = IJS,IJL
            BOUT(IJ,ITG) = NEMOSST(IJ) 
          ENDDO
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          DO IJ = IJS,IJL
            BOUT(IJ,ITG) = NEMOCICOVER(IJ) 
          ENDDO
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          IF (LNEMOCITHICK) THEN
            DO IJ = IJS,IJL
              BOUT(IJ,ITG) = NEMOCITHICK(IJ) 
            ENDDO
          ELSE
            DO IJ = IJS,IJL
              BOUT(IJ,ITG) = ZMISS
            ENDDO
          ENDIF
        ENDIF
          
        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          DO IJ = IJS,IJL
            BOUT(IJ,ITG) = NEMOUCUR(IJ) 
          ENDDO
        ENDIF
        
        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          DO IJ = IJS,IJL
            BOUT(IJ,ITG) = NEMOVCUR(IJ) 
          ENDDO
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETWMASK(WEFLXM(IJS),BOUT(IJS,ITG),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETWMASK(WEFLXD(IJS),BOUT(IJS,ITG),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETWMASK(E1012(IJS),BOUT(IJS,ITG),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETWMASK(E1214(IJS),BOUT(IJS,ITG),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETWMASK(E1417(IJS),BOUT(IJS,ITG),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETWMASK(E1721(IJS),BOUT(IJS,ITG),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETWMASK(E2125(IJS),BOUT(IJS,ITG),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

        ITG=ITG+1
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETWMASK(E2530(IJS),BOUT(IJS,ITG),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF


!       EXTRA OUTPUT FIELDS:
!       See *OUTBS* where they should have been calculated or pick them up from
!       the forcing fields (see above if the sea ice mask has to be removed)
        ITG=JPPFLAG-4
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETWMASK(WX1(IJS),BOUT(IJS,ITG),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF
        ITG=JPPFLAG-3
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETWMASK(WX2(IJS),BOUT(IJS,ITG),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF
        ITG=JPPFLAG-2
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETWMASK(WX3(IJS),BOUT(IJS,ITG),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF
        ITG=JPPFLAG-2
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETWMASK(WX4(IJS),BOUT(IJS,ITG),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF
        ITG=JPPFLAG
        IF(IPFGTBL(ITG).NE.0) THEN
          CALL SETWMASK(WX5(IJS),BOUT(IJS,ITG),IJS,IJL,CITHRSH,IODP(IJS))
        ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTSETWMASK',1,ZHOOK_HANDLE)
#endif
      END SUBROUTINE OUTSETWMASK
