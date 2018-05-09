      SUBROUTINE OUTWNORM(IJS, IJL, BOUT)

! ----------------------------------------------------------------------

!**** *WAMNORM* - PRINTS AVERAGE, MINIMUM AND MAXIMUM VALUES OF
!                 INTEGRATED PARAMETER FIELDS 

!     J. BIDLOT       ECMWF   JULY  1998 
!     D. PETTENUZZO   ADDED 9 NEW PARAMETERS FOR 3 WAVE TRAINS

!*    PURPOSE.
!     --------


!**   INTERFACE.
!     ----------

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWCOUT  , ONLY : FFLAG    ,PFLAG    ,GFLAG    ,NFLAG    ,
     &            NFLAGALL ,NTRAIN   ,IPFGTBL  ,JPPFLAG  ,COUTNAME ,    &
     &            NIPRMOUT, ITOBOUT
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWINTP  , ONLY : WHGTTG   ,WDIRTG   ,WPKFTG   ,WMNFTG   ,
     &            USTARG   ,UDIRG    ,TAUWG    ,CDG      ,SMEANG   ,
     &            U10G     ,MWP1G    ,MWP2G    ,WSPRDG   ,C4G      ,
     &            BFG      ,QPG      ,DEPTHG   ,HMAXG    ,TMAXG    ,
     &            USTOKESG ,VSTOKESG ,UCURG    ,VCURG    ,PHIEPSG  ,
     &            PHIAWG   ,TAUOCG   ,STRNMSG  ,RHOAG    ,WSTARG   ,
     &            CICG     ,CIHG     ,C3G      ,E10G     ,
     &            NSSTG    ,NICCG    ,NICTG    ,NUCUG    ,NVCUG    ,
     &            WEFLXMG  ,WEFLXDG  ,E1012G   , 
     &            E1214G   ,E1417G   ,E1721G   ,E2125G   ,E2530G   ,
     &            WX1G     ,WX2G     ,WX3G     ,WX4G     ,WX5G
      USE YOWINTS  , ONLY : WHGTSG   ,WDIRSG   ,WMNFSG   ,WHGTWG   ,
     &            WDIRWG   ,WMNFWG   ,MWP1SG   ,MWP2SG   ,WSPRDSG  ,
     &            MWP1WG   ,MWP2WG   ,WSPRDWG
      USE YOWINTT  , ONLY : WHGTTRG  ,WDIRTRG  ,MWPMTRG
      USE YOWMPP   , ONLY : IRANK    ,KTAG
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWTEST  , ONLY : IU06
      USE MPL_MODULE
! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL, NIPRMOUT
      REAL, INTENT(IN) :: BOUT(IJS:IJL,JPPFLAG)

      INTEGER :: ITG, ITT, I
      INTEGER :: ISEND, IRECV, ITAG
      INTEGER, DIMENSION(JPPFLAG) :: NAVG
      REAL :: ZHOOK_HANDLE
      REAL, DIMENSION(4,JPPFLAG) :: WNORM
      CHARACTER(LEN=47), DIMENSION(JPPFLAG) :: CNORMNAME
      LOGICAL, DIMENSION(JPPFLAG) :: LLPRNORM

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTWNORM',0,ZHOOK_HANDLE)
#endif
      IF(NFLAGALL) THEN

        IRECV=1
        LLPRNORM(:)=.FALSE.
        DO ITG = 1,JPPFLAG
        IF(NFLAG(ITG).AND.(GFLAG(ITG).OR.FFLAG(ITG).OR.PFLAG(ITG))) THEN
          LLPRNORM(ITG)=.TRUE.
!!!       all other cases are omitted until the output arrays are fused
          WNORM(:,ITG)=0.
          CNORMNAME(ITG)='not available '//COUTNAME(ITG)

          SELECT CASE (ITG)
          CASE(1)
            CNORMNAME(ITG)='   WAMNORM FOR'//COUTNAME(ITG)
       IF(IPFGTBL(ITG).EQ.IRANK) CALL WAMNORM(WHGTTG,NGX,NGY,WNORM(:,1))
          CASE(6)
            CNORMNAME(ITG)='   WAMNORM FOR'//COUTNAME(ITG)
       IF(IPFGTBL(ITG).EQ.IRANK) CALL WAMNORM(WPKFTG,NGX,NGY,WNORM(:,6))
          CASE(7)
            CNORMNAME(ITG)='   WAMNORM FOR'//COUTNAME(ITG)
       IF(IPFGTBL(ITG).EQ.IRANK) CALL WAMNORM(CDG,NGX,NGY,WNORM(:,7))
          CASE(10)
            CNORMNAME(ITG)='   WAMNORM FOR'//COUTNAME(ITG)
       IF(IPFGTBL(ITG).EQ.IRANK) CALL WAMNORM(U10G,NGX,NGY,WNORM(:,10))
          END SELECT 

          ISEND=IPFGTBL(ITG)

          IF(ISEND.NE.IRECV) THEN
            ITAG=KTAG+ITG
            CALL GSTATS(675,0)
            IF(IRANK.EQ.ISEND) THEN
              CALL MPL_SEND(WNORM(:,ITG),KDEST=IRECV,KTAG=ITAG,         &
     &                      CDSTRING='OUTWNORM:')
            ELSEIF(IRANK.EQ.IRECV) THEN
              CALL MPL_RECV(WNORM(:,ITG),KSOURCE=ISEND,KTAG=ITAG,       &
     &                      CDSTRING='OUTWNORM:')
            ENDIF
            CALL GSTATS(675,1)
          ENDIF

        ENDIF
        ENDDO

!       WRITE NORM TO LOGFILE
        IF(IRANK.EQ.IRECV) THEN
          WRITE(IU06,*) ' ' 
          WRITE(IU06,*) '   WAMNORM ON ',CDTPRO
          WRITE(IU06,*) '          AVERAGE            MINIMUM  ',
     &     '           MAXIMUM        SEA POINTS' 
          DO ITG = 1, JPPFLAG
            IF(LLPRNORM(ITG)) THEN
              WRITE(IU06,*) CNORMNAME(ITG)
              WRITE(IU06,*) '  ',(WNORM(I,ITG),I=1,4)
              WRITE(IU06,111) (WNORM(I,ITG),I=1,4)
            ENDIF
          ENDDO
        ELSE
          WRITE(IU06,*) ' ' 
          WRITE(IU06,*) '   WAMNORM ON ',CDTPRO,' SEE OUTPUT PE',IRECV
        ENDIF
111     FORMAT(4x,'HEX: ',4(Z16.16,2x))
        CALL FLUSH(IU06)

      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTWNORM',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE OUTWNORM
