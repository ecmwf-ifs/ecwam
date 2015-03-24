      SUBROUTINE CLEAN_OUTBS(LLCLEAN_OUTBS)

! ----------------------------------------------------------------------

!****  *CLEAN_OUTBS* - CLEAN ALL ARRAYS THAT WERE LEFT OVER FROM 
!                      CALLING OUTBS.


!     PURPOSE.                                                          
!     --------                                                          

!     DEALLOCATES ARRAYS THAT WERE CREATED BY CALLING OUTBS


!*    INTERFACE.                                                        
!     ----------                                                        

!     *CALL* *CLEAN_OUTBS(LLCLEAN_OUTBS)*

!          *LLCLEAN_OUTBS*  LOGICAl IF TRUE ALL ARRAYS ARE DELETE
!                           SET TO FALSE ON LEAVING

!     METHOD.                                                           
!     -------                                                           

!     EXTERNALS.                                                        
!     ----------                                                        

!     REFERENCES.                                                       
!     -----------                                                       

! ----------------------------------------------------------------------

      USE YOWINTP  , ONLY : WHGTTG   ,WDIRTG   ,WPKFTG   ,WMNFTG   ,
     &            USTARG   ,UDIRG    ,TAUWG    ,CDG      ,SMEANG   ,
     &            U10G     ,WHGTAG   ,CWHTAG   ,RANGAG   ,MWP1G    ,
     &            MWP2G    ,WSPRDG   ,C4G      ,BFG      ,QPG      ,
     &            DEPTHG   ,HMAXG    ,TMAXG    ,USTOKESG ,VSTOKESG ,
     &            UCURG    ,VCURG    ,PHIEPSG  ,PHIAWG   ,TAUOCG   ,
     &            STRNMSG  ,RHOAG    ,WSTARG   ,CICG     ,CIHG     ,
     &            C3G      ,
     &            WX1G     ,WX2G     ,WX3G     ,WX4G     ,WX5G     ,
     &            NSSTG    ,NICCG    ,NICTG    ,NUCUG    ,NVCUG
      USE YOWINTS  , ONLY : WHGTSG   ,WDIRSG   ,WMNFSG   ,WHGTWG   ,
     &            WDIRWG   ,WMNFWG
      USE YOWINTT  , ONLY : WHGTTRG  ,WDIRTRG  ,MWPMTRG 
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      REAL :: ZHOOK_HANDLE

      LOGICAL, INTENT(INOUT) :: LLCLEAN_OUTBS

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('CLEAN_OUTBS',0,ZHOOK_HANDLE)
#endif

      IF(LLCLEAN_OUTBS) THEN
!     DEALLOCATE OUTPUT ARRAYS

        IF(ALLOCATED(WHGTTG)) DEALLOCATE(WHGTTG)
        IF(ALLOCATED(WDIRTG)) DEALLOCATE(WDIRTG)
        IF(ALLOCATED(WPKFTG)) DEALLOCATE(WPKFTG)
        IF(ALLOCATED(WMNFTG)) DEALLOCATE(WMNFTG)
        IF(ALLOCATED(USTARG)) DEALLOCATE(USTARG)
        IF(ALLOCATED(UDIRG)) DEALLOCATE(UDIRG)
        IF(ALLOCATED(TAUWG)) DEALLOCATE(TAUWG)
        IF(ALLOCATED(CDG)) DEALLOCATE(CDG)
        IF(ALLOCATED(SMEANG)) DEALLOCATE(SMEANG)
        IF(ALLOCATED(U10G)) DEALLOCATE(U10G)
        IF(ALLOCATED(MWP1G)) DEALLOCATE(MWP1G)
        IF(ALLOCATED(MWP2G)) DEALLOCATE(MWP2G)
        IF(ALLOCATED(WSPRDG)) DEALLOCATE(WSPRDG)
        IF(ALLOCATED(C4G)) DEALLOCATE(C4G)
        IF(ALLOCATED(BFG)) DEALLOCATE(BFG)
        IF(ALLOCATED(QPG)) DEALLOCATE(QPG)
        IF(ALLOCATED(DEPTHG)) DEALLOCATE(DEPTHG)
        IF(ALLOCATED(HMAXG)) DEALLOCATE(HMAXG)
        IF(ALLOCATED(TMAXG)) DEALLOCATE(TMAXG)
        IF(ALLOCATED(USTOKESG)) DEALLOCATE(USTOKESG)
        IF(ALLOCATED(VSTOKESG)) DEALLOCATE(VSTOKESG)
        IF(ALLOCATED(UCURG)) DEALLOCATE(UCURG)
        IF(ALLOCATED(VCURG)) DEALLOCATE(VCURG)
        IF(ALLOCATED(PHIEPSG)) DEALLOCATE(PHIEPSG)
        IF(ALLOCATED(PHIAWG)) DEALLOCATE(PHIAWG)
        IF(ALLOCATED(TAUOCG)) DEALLOCATE(TAUOCG)
        IF(ALLOCATED(STRNMSG)) DEALLOCATE(STRNMSG)
        IF(ALLOCATED(RHOAG)) DEALLOCATE(RHOAG)
        IF(ALLOCATED(WSTARG)) DEALLOCATE(WSTARG)
        IF(ALLOCATED(CICG)) DEALLOCATE(CICG)
        IF(ALLOCATED(CIHG)) DEALLOCATE(CIHG)
        IF(ALLOCATED(C3G)) DEALLOCATE(C3G)
        IF(ALLOCATED(NSSTG)) DEALLOCATE(NSSTG)
        IF(ALLOCATED(NICCG)) DEALLOCATE(NICCG)
        IF(ALLOCATED(NICTG)) DEALLOCATE(NICTG)
        IF(ALLOCATED(NUCUG)) DEALLOCATE(NUCUG)
        IF(ALLOCATED(NVCUG)) DEALLOCATE(NVCUG)
        IF(ALLOCATED(WX1G)) DEALLOCATE(WX1G)
        IF(ALLOCATED(WX2G)) DEALLOCATE(WX2G)
        IF(ALLOCATED(WX3G)) DEALLOCATE(WX3G)
        IF(ALLOCATED(WX4G)) DEALLOCATE(WX4G)
        IF(ALLOCATED(WX5G)) DEALLOCATE(WX5G)

        IF(ALLOCATED(WHGTSG)) DEALLOCATE(WHGTSG)
        IF(ALLOCATED(WDIRSG)) DEALLOCATE(WDIRSG)
        IF(ALLOCATED(WMNFSG)) DEALLOCATE(WMNFSG)
        IF(ALLOCATED(WHGTWG)) DEALLOCATE(WHGTWG)
        IF(ALLOCATED(WDIRWG)) DEALLOCATE(WDIRWG)
        IF(ALLOCATED(WMNFWG)) DEALLOCATE(WMNFWG)

        IF(ALLOCATED(WHGTTRG)) DEALLOCATE(WHGTTRG)
        IF(ALLOCATED(WDIRTRG)) DEALLOCATE(WDIRTRG)
        IF(ALLOCATED(MWPMTRG)) DEALLOCATE(MWPMTRG)

        IF(ALLOCATED(WHGTAG)) DEALLOCATE(WHGTAG)
        IF(ALLOCATED(CWHTAG)) DEALLOCATE(CWHTAG)
        IF(ALLOCATED(RANGAG)) DEALLOCATE(RANGAG)

        LLCLEAN_OUTBS=.FALSE.
      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('CLEAN_OUTBS',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE CLEAN_OUTBS
