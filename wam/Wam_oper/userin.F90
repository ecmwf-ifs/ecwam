      SUBROUTINE USERIN (IFORCA, LWCUR)

! ----------------------------------------------------------------------

!**** *USERIN* - ROUTINE TO READ AND WRITE NAMELIST INPUT.

!     H. GUNTHER     GKSS/ECMWF   NOVEMBER 1989
!     J. BIDLOT      ECMWF        JUNE  1996
!     J. BIDLOT      ECMWF        JUNE  1997: coupled/uncoupled
!                                             restarted/normal
!     B. HANSEN      ECMWF        JULY  1997  COUPLED AND ANALYSIS
!                                             (IE. PREPAN)
!     B. HANSEN      ECMWF        JANUARY 1998 NAMELIST INPUT.
!     B. HANSEN      ECMWF        NOVEMBER 1998 WPOSNAM USED TO POSITION
!                                                 NAMLIST FOR READING.
!     S. ABDALLA     ECMWF        OCTOBER 2001  LGUST & LADEN ADDED TO
!                                               CONTROL GUSTINESS & AIR DENSITY 
!     J BIDLOT       ECMWF        MARCH 2008    LLWSWAVE AND LLWDWAVE
!                                               ADDED TO CONTROL WHETHER OR
!                                               NOT PARAMETER 245 AND/OR
!                                               PARAMETER 249 ARE USED TO
!                                               GENERATE THE WIND FORCING
!                                               (ONLY FOR UNCOUPLED RUNS).
!     D. PETTENUZZO               MAY 2012      ADDED 9 NEW PARAMETERS

!*    PURPOSE.
!     --------

!       READ USER INPUT CONCERNING PERIOD OF INTEREST,TIMESTEPS AND
!       MODEL OPTIONS TO INITIALIZE COMMON USERD. A CONSISTENCY CHECK
!       IS DONE TOO.

!**   INTERFACE.
!     ----------

!       *CALL* *USERIN (IFORCA, LWCUR)*
!          *IFORCA* -  FORCAST START OPTION.
!          *LWCUR* -  INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS


!     METHOD.
!     -------

!        USER INFORMATION IS BEING READ WITH THE PRESUMPTIONS THAT:
!         1. EVERY LINE STARTING WITH 'C' IS A COMMENT LINE
!         2. VALUES ARE PUT IN BELOW POSITIONS INDICATED WITH '-'
!            (RIGHT-JUSTIFIED)

!     EXTERNALS.
!     ----------

!       *ABORT1*    - TERMINATES PROCESSING.
!       *DIFDATE*   - COMPUTES A TIME DIFFERENCE.
!       *IINITFDB*  - INITIALISE DATA BASE.
!       *WPOSNAM*   - POSITION NAMELIST FOR READING AND CONTROLLED
!                     TRAPPING OF EOF.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWALTAS , ONLY : NUMALT   ,IBUFRSAT  ,ALTSDTHRSH,ALTBGTHRSH, &
     &            HSALTCUT, LALTGRDOUT, LALTPAS, LALTPASSIV,            &
     &            XKAPPA2  ,HSCOEFCOR,HSCONSCOR ,LALTCOR   ,LALTLRGR,   &
     &            LODBRALT ,CSATNAME
      USE YOWCOUP  , ONLY : LWCOU    ,KCOUSTEP  ,LWFLUX, LWVFLX_SNL,    &
     &            LWNEMOCOU, LWNEMOCOUSEND, LWNEMOCOURECV 
      USE YOWCOUT  , ONLY : COUTT    ,COUTS    ,CASS     ,FFLAG    ,    &
     &            FFLAG20  ,                                            &
     &            GFLAG    ,                                            &
     &            GFLAG20  ,NFLAG    ,                                  &
     &            IPRMINFO ,                                            &
     &            IRWDIR, IRCD ,IRU10  , IRALTHS ,IRALTHSC ,IRALTRC ,   &
     &            IRHS     ,IRTP     ,IRT1     ,IRPHIOC  ,IRTAUOC   ,   &
     &            IRBATHY  ,                                            &
     &            NFLAGALL ,UFLAG    ,LFDB     ,NOUTT    ,NOUTS    ,    &
     &            IRCD     ,IRU10    ,                                  &
     &            NASS     ,JPPFLAG  ,                                  &
     &            LRSTPARALR, LRSTPARALW,                               &
     &            COUTNAME,                                             &
     &            NTRAIN   ,                                            &
     &            IPFGTBL  ,                                            &
     &            LWAMANOUT,                                            &
     &            NWRTOUTWAM,                                           &
     &            LWFLUXOUT,                                            &
     &            LSECONDORDER,                                         &
     &            LWAM_USE_IO_SERV
      USE YOWCPBO  , ONLY : IBOUNC
      USE YOWCURR  , ONLY : IDELCUR  ,CDATECURA
      USE YOWFPBO  , ONLY : IBOUNF
      USE YOWGRIBHD, ONLY : LGRHDIFS ,LNEWLVTP ,IMDLGRBID_G,IMDLGRBID_M 
      USE YOWGRIB_HANDLES , ONLY : NGRIB_HANDLE_IFS
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE ,LCIWABR  ,              &
     &            CITHRSH  ,CIBLOCK  ,LICETH   ,                        &
     &            CITHRSH_SAT, CITHRSH_TAIL    ,CDICWA
      USE YOWMESPAS, ONLY : LMESSPASS,                                  &
     &            LFDBIOOUT,LGRIBIN  ,LGRIBOUT ,LNOCDIN
      USE YOWMPP   , ONLY : NPROC
      USE YOWPARAM , ONLY : SWAMPWIND,SWAMPWIND2,DTNEWWIND,LTURN90 ,    &
     &            SWAMPCIFR,SWAMPCITH,LWDINTS  ,LL1D     ,CLDOMAIN
      USE YOWPHYS  , ONLY : ALPHAPMAX
      USE YOWSTAT  , ONLY : CDATEE   ,CDATEF   ,CDATER   ,CDATES   ,    &
     &            CFDBSF   ,CFDB2DSP ,IDELPRO  ,IDELT    ,IDELWI   ,    &
     &            IDELWO   ,IDELALT  ,IREST    ,IDELRES  ,IDELINT  ,    &
     &            IDELBC   ,                                            &
     &            ICASE    ,ISHALLO  ,                                  &
     &            ISNONLIN ,                                            &
     &            IPHYS    ,                                            &
     &            IDAMPING ,                                            &
     &            LBIWBK   ,                                            &
     &            IREFRA   ,IPROPAGS ,IASSI    ,NTASKS   ,NSIZE    ,    &
     &            NENSFNB  ,NTOTENS  ,NSYSNB   ,NMETNB   ,CDATEA   ,    &
     &            YCLASS   ,YEXPVER  ,L4VTYPE  ,LFRSTFLD ,LALTAS   ,    &
     &            LSARAS   ,LSARINV  ,ISTREAM  ,NLOCGRB  ,NCONSENSUS,   &
     &            NDWD     ,NMFR     ,NNCEP    ,NUKM     ,IREFDATE ,    &
     &            LGUST    ,LADEN    ,LSUBGRID ,LLSOURCE ,              &
     &            LNSESTART,                                            &
     &            LSMSSIG_WAM,CMETER ,CEVENT   ,                        &
     &            LRELWIND ,                                            &
     &            IDELWI_LST, IDELWO_LST, CDTW_LST, NDELW_LST
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWTEXT  , ONLY : LRESTARTED,ICPLEN   ,USERID   ,RUNID    ,   &
     &            PATH     ,CPATH    ,CWI
      USE YOWUNIT  , ONLY : IU04     ,IU20     ,IU30     ,IU32
      USE YOWUNPOOL, ONLY : LLUNSTR  ,LPREPROC, LVECTOR, IVECTOR   ,    &
     &            LLUNBINOUT
      USE YOWWAMI  , ONLY : CBEGDT   ,CENDDT   ,CBPLTDT  ,CEPLTDT  ,    &
     &            CLSPDT   ,CRSTDT   ,IANALPD  ,IFOREPD  ,IDELWIN  ,    &
     &            IASSIM   ,NFCST    ,ISTAT
      USE YOWWIND  , ONLY : CWDFILE  ,LLWSWAVE ,LLWDWAVE ,RWFAC
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE FDBSUBS_MOD, ONLY : IINITFDBSUBS
      USE GRIB_API_INTERFACE

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "difdate.intfb.h"
#include "mpcrtbl.intfb.h"
#include "outxt.intfb.h"
#include "readsta.intfb.h"
#include "set_wflags.intfb.h"
#include "wam_u2l1cr.intfb.h"
#include "wstream_strg.intfb.h"

      INTEGER(KIND=JWIM) :: IFORCA 
      INTEGER(KIND=JWIM) :: ITG, IC, I, J, ISAT
      INTEGER(KIND=JWIM) :: LEN, ILCFDB2DSP, ILCFDBSF
      INTEGER(KIND=JWIM) :: IFS_STREAM, KSTREAM
      INTEGER(KIND=JWIM) :: IDELT_NEW
      INTEGER(KIND=JWIM) :: ISHIFT
      INTEGER(KIND=JWIM) :: IDELPRO_NEW
      INTEGER(KIND=JWIM) :: MINBUFRSAT, MAXBUFRSAT
      INTEGER(KIND=JWIM) :: INITFDB_STAT
      INTEGER(KIND=JWIM) :: IUWDFILE
      INTEGER(KIND=JWIM) :: IWTIME, IWTIME_old
      INTEGER(KIND=JWIM) :: IWAM_GET_UNIT
      INTEGER(KIND=JWIM) :: NPROC_RST

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: WSPEED, WTHETA

      CHARACTER(LEN=2) :: MARSFCTYPE
      CHARACTER(LEN=3) :: CITG
      CHARACTER(LEN=4) :: CSTREAM
      CHARACTER(LEN=256) :: CLFORM

      LOGICAL :: LERROR, LASTREAM
      LOGICAL :: LLNALTGO
      LOGICAL :: LWCUR
      LOGICAL :: LRSTPARAL

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('USERIN',0,ZHOOK_HANDLE)

! !!!!!! THE READING OF NALINE HAS BEEN MOVED TO MPUSERIN !!!!!!!

      IFORCA=1
      NPROC_RST=0


!*    1. ANALYSE NAMELIST NALINE.
!        ------------------------

      IF (IDAMPING.NE.0 .AND. IDAMPING.NE.1) THEN
        WRITE(IU06,*) ' WRONG VALUE FOR IDAMPING !!!'
        WRITE(IU06,*) ' IDAMPING =',IDAMPING
        CALL ABORT1
      ENDIF

      IF (.NOT.LWCOU .AND. LODBRALT) THEN
        WRITE(IU06,*)'WARNING IN SUBROUTINE USERIN:'
        WRITE(IU06,*)'        ODB IS NOT READY FOR STAND-ALONE MODEL'
        WRITE(IU06,*)'        ODB WILL NOT BE USED!'
        LODBRALT = .FALSE.  ! ODB is not ready for stand-alone wave model
                            ! this needs to be changed later.
      ENDIF

      IF(.NOT.LWCOU) LGRHDIFS = .FALSE.  ! by definition
      IF(LGRHDIFS) THEN
!       GET ISTREAM THAT CORRESPONDS TO IFS_STREAM
        CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'stream',IFS_STREAM)
        IF(.NOT.LNEWLVTP) THEN
          CALL WSTREAM_STRG(IFS_STREAM, CSTREAM, NENSFNB, NTOTENS,      &
     &                      MARSFCTYPE, ISTREAM, LASTREAM)
          IF(CSTREAM.EQ.'****') THEN
            WRITE(IU06,*) '*****************************************'
            WRITE(IU06,*) ''
            WRITE(IU06,*) ' ERROR IN USERIN !!!!'
            WRITE(IU06,*) ' IFS STREAM UNKNOWN '
            WRITE(IU06,*) ' IFS STREAM = ', IFS_STREAM
            WRITE(IU06,*) ' BUT NOT DEFINED IN WSTREAM_STRG !!!!'
            WRITE(IU06,*) ''
            WRITE(IU06,*) '*****************************************'
            CALL ABORT1
          ENDIF
        ELSE
          ISTREAM=IFS_STREAM
        ENDIF
      ELSEIF(ISTREAM.LE.0) THEN
        WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*)'+                                          +'
        WRITE(IU06,*)'+ SUBROUTINE USERIN :                      +'
        WRITE(IU06,*)'+ READ NAMELIST FAILED                     +' 
        WRITE(IU06,*)'+ ISTREAM MUST BE SPECIFIED > 0 !!!!       +' 
        WRITE(IU06,*)'+ PROGRAM WILL ABORT                       +'
        WRITE(IU06,*)'+                                          +'
        WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
        CALL ABORT1
      ENDIF

      IF(ISTREAM.LE.0) THEN
        WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*)'+                                          +'
        WRITE(IU06,*)'+ SUBROUTINE USERIN :                      +'
        WRITE(IU06,*)'+ READ NAMELIST FAILED                     +' 
        WRITE(IU06,*)'+ ISTREAM MUST BE SPECIFIED > 0 !!!!       +' 
        WRITE(IU06,*)'+ PROGRAM WILL ABORT                       +'
        WRITE(IU06,*)'+                                          +'
        WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
        CALL ABORT1
      ENDIF

      IF (LWNEMOCOU.AND..NOT.LWNEMOCOUSEND.AND..NOT.LWNEMOCOURECV) THEN
        WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*) '+                                             +'
        WRITE(IU06,*) '+ SUBROUTINE USERIN :                         +'
        WRITE(IU06,*) '+ WAM<->NEMO COUPLING ACTIVE                  +'
        WRITE(IU06,*) '+ BUT BOTH RECV AND SEND IS FALSE             +'
        WRITE(IU06,*) '+ THIS IS NOT SUPPOSED TO WORK                +'
        WRITE(IU06,*) '+                                             +'
        WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
        CALL ABORT1
      ENDIF

      IF(.NOT.LWCOU) THEN
        CDATEA = CBPLTDT
        CDATEE = CEPLTDT
      ELSE
        CBPLTDT = CDATEA
        CEPLTDT = CDATEF
      ENDIF
      CALL WAM_U2L1CR( YCLASS )


!           **** CHECK LENGTH OF YEXPVER AND PUT IT RIGHT JUSTIFIED ****
      LEN=LEN_TRIM(YEXPVER)
      YEXPVER(5-LEN:4)=YEXPVER(1:LEN)
      YEXPVER(1:4-LEN)='0000'

      IF(LWCOU) THEN
!       TIME STEP SELECTION: 
        IF(KCOUSTEP.LE.IDELPRO .AND. KCOUSTEP.LE.IDELT) THEN
!         TIGHT COUPLING
          IDELPRO=KCOUSTEP
          IDELT=KCOUSTEP
        ELSE IF(KCOUSTEP.GE.IDELPRO .AND. KCOUSTEP.LE.IDELT) THEN
!         IDELT IS SET BY COUPLING
          IDELT=KCOUSTEP
          DO IC=1,IDELT
            IDELPRO_NEW=IDELT/IC
            IF(IDELPRO_NEW*IC.EQ.IDELT.AND.IDELPRO_NEW.LE.IDELPRO )EXIT
          ENDDO
          IDELPRO=IDELPRO_NEW
        ELSE IF(KCOUSTEP.LE.IDELPRO.AND.KCOUSTEP.GE.IDELT) THEN
!         IDELPRO IS SET BY COUPLING
          IDELPRO=KCOUSTEP
          DO IC=1,IDELPRO
            IDELT_NEW=IDELPRO/IC
            IF(IDELT_NEW*IC.EQ.IDELPRO.AND.IDELT_NEW.LE.IDELT) EXIT
          ENDDO
          IDELT=IDELT_NEW
        ELSE
          IF(IDELPRO.LE.IDELT) THEN
            DO IC=1,KCOUSTEP
              IDELT_NEW=KCOUSTEP/IC
              IF(IDELT_NEW*IC.EQ.KCOUSTEP.AND.IDELT_NEW.LE.IDELT) EXIT
            ENDDO
            IDELT=IDELT_NEW
  
            DO IC=1,IDELT
            IDELPRO_NEW=IDELT/IC
              IF(IDELPRO_NEW*IC.EQ.IDELT .AND. IDELPRO_NEW.LE.IDELPRO) EXIT
            ENDDO
            IDELPRO=IDELPRO_NEW

          ELSE
            DO IC=1,KCOUSTEP
            IDELPRO_NEW=KCOUSTEP/IC
              IF(IDELPRO_NEW*IC.EQ.KCOUSTEP .AND. IDELPRO_NEW.LE.IDELPRO) EXIT
            ENDDO
            IDELPRO=IDELPRO_NEW

            DO IC=1,IDELPRO
              IDELT_NEW=IDELPRO/IC
              IF(IDELT_NEW*IC.EQ.IDELPRO.AND.IDELT_NEW.LE.IDELT) EXIT
            ENDDO
            IDELT=IDELT_NEW
          ENDIF
        ENDIF

!       IF COUPLED MODEL
!       SET TEMPORARILY IDELWO and IDELWI TO IDELPRO
!       TO AVOID ALL THE CHECKS. SEE IN INITMDL WHERE THEY ARE RESET
        IDELWO=IDELPRO
        IDELWI=IDELPRO

      ELSE
!           **** FORCING TIME STEPPING AT SPECIFIED TIMES ****
!        SEE MPUSERIN

      ENDIF

!* CHECK FLAG FOR GRIBING AS SOME OPTIONS ARE NOT IMPLEMENTED
!* IN OUTINT
      IF (IREST.GT.0 .AND. LGRIBOUT .AND. .NOT.GFLAG(IRWDIR)) THEN
        GFLAG(IRWDIR) = .TRUE.
        WRITE(IU06,*) ' ******** NOTE *** NOTE **** NOTE *************' 
        WRITE(IU06,*) ' YOU HAVE REQUESTED RESTART FILE IN GRIB FORMAT' 
        WRITE(IU06,*) ' BUT DID NOT REQUEST OUTPUT OF WIND DIRECTION. '
        WRITE(IU06,*) ' THE MODEL WILL RESET THE OUTPUT OPTION TO GET ' 
        WRITE(IU06,*) ' IT OUT AT TIME PRESCRIBED BY THE OUTPUT OF ALL'
        WRITE(IU06,*) ' OTHER INTEGRATED PARAMETERS.'
        WRITE(IU06,*) ' ******** NOTE *** NOTE **** NOTE *************' 
        WRITE(IU06,*) ' '
      ENDIF
      IF (IREST.GT.0 .AND. LGRIBOUT .AND. .NOT.GFLAG(IRCD)) THEN
        GFLAG(IRCD) = .TRUE.
        WRITE(IU06,*) ' ******** NOTE *** NOTE **** NOTE *************' 
        WRITE(IU06,*) ' YOU HAVE REQUESTED RESTART FILE IN GRIB FORMAT' 
        WRITE(IU06,*) ' BUT FAILED TO ASK FOR OUTPUT OF THE DRAG COEF.'
        WRITE(IU06,*) ' THE MODEL WILL RESET THE OUTPUT OPTION TO GET ' 
        WRITE(IU06,*) ' CD OUT AT TIME PRESCRIBED BY THE OUTPUT OF ALL'
        WRITE(IU06,*) ' OTHER INTEGRATED PARAMETERS.'
        WRITE(IU06,*) ' ******** NOTE *** NOTE **** NOTE *************' 
        WRITE(IU06,*) ' '
      ENDIF
      IF (IREST.GT.0 .AND. LGRIBOUT .AND. .NOT.GFLAG(IRU10)) THEN
        GFLAG(IRU10) = .TRUE.
        WRITE(IU06,*) ' ******** NOTE *** NOTE **** NOTE *************' 
        WRITE(IU06,*) ' YOU HAVE REQUESTED RESTART FILE IN GRIB FORMAT' 
        WRITE(IU06,*) ' BUT FAILED TO ASK FOR OUTPUT OF U10wave'
        WRITE(IU06,*) ' THE MODEL WILL RESET THE OUTPUT OPTION TO GET ' 
        WRITE(IU06,*) ' u10 OUT AT TIME PRESCRIBED BY THE OUTPUT OF ALL'
        WRITE(IU06,*) ' OTHER INTEGRATED PARAMETERS.'
        WRITE(IU06,*) ' ******** NOTE *** NOTE **** NOTE *************' 
        WRITE(IU06,*) ' '
      ENDIF
      IF (GFLAG(IRALTHS) .AND. (IASSI.NE.1.OR..NOT.LALTAS)) THEN
        GFLAG(IRALTHS) = .FALSE.
      ENDIF
      IF (GFLAG(IRALTHSC) .AND. (IASSI.NE.1.OR..NOT.LALTAS)) THEN
        GFLAG(IRALTHSC) = .FALSE.
      ENDIF
      IF (GFLAG(IRALTRC) .AND. (IASSI.NE.1.OR..NOT.LALTAS)) THEN
        GFLAG(IRALTRC) = .FALSE.
      ENDIF


!     SET INTEGRATED OUTPUT PARAMETER TABLE FOR MESSAGE PASSING
!     IN CASE OF NO MESSAGE PASSING THE TABLE IS STILL USED TO
!     POINT TO THE PARAMETRE WHICH ARE OUTPUT.
      CALL MPCRTBL
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '    SUB. INITMDL: MPCRTBL DONE'
        CALL FLUSH (IU06)
      ENDIF

      FFLAG20 = SET_WFLAGS(FFLAG,JPPFLAG)
      GFLAG20 = SET_WFLAGS(GFLAG,JPPFLAG)
      NFLAGALL= SET_WFLAGS(NFLAG,JPPFLAG) 

!     ARE THE OCEAN FLUXES OUTPUT PARAMETERS:
      LWFLUXOUT = .TRUE.

!     DEFINE GRIDDED OUTPUT ARRAYS WHICH ARE USED IN THE ASSIMILATION
      UFLAG=.FALSE.
      IF(IASSI.EQ.1) THEN
        UFLAG(1)=.TRUE.
      ENDIF

      ICPLEN=LEN_TRIM(CPATH)
      IF(ICPLEN.GT.0.AND.CPATH(ICPLEN:ICPLEN).EQ.'/') THEN
        CPATH=CPATH(1:ICPLEN-1)
        ICPLEN=ICPLEN-1
      ENDIF
      CWI=CPATH(1:ICPLEN)//'/waminfo'

      IF(LFDBIOOUT.AND..NOT.LGRIBOUT) THEN
        WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*)'+                                          +'
        WRITE(IU06,*)'+ LFDBIOOUT = TRUE AND LGRIBOUT = FALSE    +'
        WRITE(IU06,*)'+ IS AN OBSOLETE OPTION:                   +'
        WRITE(IU06,*)'+ LFDBIOOUT IS RESET TO FALSE              +'
        WRITE(IU06,*)'+                                          +'
        WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
        LFDBIOOUT=.FALSE.
      ENDIF

!*    1.1  READ THE WAMINFO FILE AND OVERWRITE INPUT.
!          ------------------------------------------

      INQUIRE(FILE=CWI,EXIST=LRESTARTED)
      IF (LRESTARTED) THEN
        IU04 =  IWAM_GET_UNIT (IU06, CWI(1:ICPLEN+8) , 'r', 'f', 0)
        CALL READSTA(IU04, CBEGDT, CENDDT, IANALPD, IFOREPD, IDELWIN,   &
     &               CRSTDT, CLSPDT, CBPLTDT, CEPLTDT, IASSIM, NFCST,   &
     &               ISTAT, CDATECURA, LRSTPARAL, NPROC_RST)
        CLOSE (IU04)
        WRITE(IU06,*)'+++++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*)'+                                           +'
        WRITE(IU06,*)'+    WARNING MODEL RESTART FROM RC FILE     +'
        WRITE(IU06,*)'+    ==================================     +'
        WRITE(IU06,*)'+                                           +'
        WRITE(IU06,*)'+ RC: ', CWI(1:ICPLEN+8)
        WRITE(IU06,*)'+                                           +'
        WRITE(IU06,*)'+ START DATE OF RUN         CDATEA  = ', CDATEA
        WRITE(IU06,*)'+ START DATE OF RESTART     CBEGDT  = ', CBEGDT
        WRITE(IU06,*)'+ ANALYSIS PERIOD (SECONDS) IANALPD = ', IANALPD
        WRITE(IU06,*)'+ FORECAST PERIOD (SECONDS) IFOREPD = ', IFOREPD
        WRITE(IU06,*)'+ END DATE FROM WAMINFO     CENDDT  = ', CENDDT
        WRITE(IU06,*)'+ END DATE FROM INPUT       CDATEE  = ', CDATEE
        WRITE(IU06,*)'+ START DATE FOR CURRENTS CDATECURA = ',CDATECURA
        WRITE(IU06,*)'+                                           +'
!       OVER-WRITE LRSTPARALR WITH WHAT THE RESTART FILE WANTS
        LRSTPARALR=LRSTPARAL
        IF(LRSTPARALR) THEN
          IF(NPROC_RST .NE. NPROC) THEN
            WRITE(IU06,*)'+                                           +'
            WRITE(IU06,*)'+++++++++++++++++++++++++++++++++++++++++++++'
            WRITE(IU06,*)'+                                           +'
            WRITE(IU06,*)'+   F A T A L   E R R O R  IN SUB. USERIN   +'
            WRITE(IU06,*)'+   =====================================   +'
            WRITE(IU06,*)'+                                           +'
            WRITE(0,*)'+   RESTART FAILED.                         +'
            WRITE(0,*)'+                                           +'
            WRITE(0,*)'+   PARALLEL READ OPTION WAS SELECTED, BUT  +'
            WRITE(0,*)'+   THE NUMBER OF MPI TASKS THAT WERE USED  +'
            WRITE(0,*)'+   NPROC_RST = ',NPROC_RST
            WRITE(0,*)'+   IS NOT EQUAL THE CURRENT NUMBER !       +'
            WRITE(0,*)'+   NPROC = ',NPROC
            WRITE(IU06,*)'+   RESTART FAILED.                         +'
            WRITE(IU06,*)'+                                           +'
            WRITE(IU06,*)'+   PARALLEL READ OPTION WAS SELECTED, BUT  +'
            WRITE(IU06,*)'+   THE NUMBER OF MPI TASKS THAT WERE USED  +'
            WRITE(IU06,*)'+   NPROC_RST = ',NPROC_RST
            WRITE(IU06,*)'+   IS NOT EQUAL THE CURRENT NUMBER !       +'
            WRITE(IU06,*)'+   NPROC = ',NPROC
            WRITE(IU06,*)'+                                           +'
            WRITE(IU06,*)'+   ABORT SERVICE ROUTINE CALLED BY USERIN  +'
            WRITE(IU06,*)'+   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  +'
            CALL ABORT1
          ELSE
            WRITE(IU06,*)'+ PARALLEL READ OPTION WAS SELECTED         +'
          ENDIF
        ELSE
          WRITE(IU06,*)'+ NON PARALLEL READ OPTION WAS SELECTED     +'
        ENDIF  
        WRITE(IU06,*)'+ MODEL CONTINUES USING RESTART DATE.       +'
        WRITE(IU06,*)'+                                           +'
        WRITE(IU06,*)'+++++++++++++++++++++++++++++++++++++++++++++'

        CDATEF = CEPLTDT
        IF (.NOT. LWCOU) THEN
          CDATEA = CBEGDT
          CDATEE = CENDDT
        ENDIF

        CDATER = CRSTDT
        CDATES = CLSPDT
        IDELWI = IDELWIN
        IASSI = IASSIM
        IFORCA = NFCST

        IF (.NOT. LWCOU) THEN
          IF (CDATEE.NE.CENDDT) THEN
            WRITE(IU06,*)'+++++++++++++++++++++++++++++++++++++++++++++'
            WRITE(IU06,*)'+                                           +'
            WRITE(IU06,*)'+    WARNING ERROR IN SUB. USERIN           +'
            WRITE(IU06,*)'+    ============================           +'
            WRITE(IU06,*)'+ END DATE OF RUN AND RUN LENGTH DO NOT     +'
            WRITE(IU06,*)'+ MATCH IN THE WAM INFO FILE.               +'
            WRITE(IU06,*)'+ START DATE OF RUN         CBEGDT =', CBEGDT
            WRITE(IU06,*)'+ ANALYSIS PERIOD (SECONDS) IANALPD=', IANALPD
            WRITE(IU06,*)'+ FORECAST PERIOD (SECONDS) IFOREPD=', IFOREPD
            WRITE(IU06,*)'+ END DATE FROM WAMINFO     CENDDT =', CENDDT
            WRITE(IU06,*)'+ END DATE COMPUTED         CDATEE =', CDATEE
            WRITE(IU06,*)'+                                           +'
            WRITE(IU06,*)'+ MODEL CONTINUES USING COMPUTED END DATE.  +'
            WRITE(IU06,*)'+                                           +'
            WRITE(IU06,*)'+++++++++++++++++++++++++++++++++++++++++++++'
            CENDDT = CDATEE
          ENDIF
        ELSE
          IF ( CDATEA .NE. CBEGDT ) THEN
            WRITE(IU06,*)'+                                           +'
            WRITE(IU06,*)'+++++++++++++++++++++++++++++++++++++++++++++'
            WRITE(IU06,*)'+                                           +'
            WRITE(IU06,*)'+   F A T A L   E R R O R  IN SUB. USERIN   +'
            WRITE(IU06,*)'+   =====================================   +'
            WRITE(IU06,*)'+                                           +'
            WRITE(IU06,*)'+   COUPLED RESTART FAILED.                 +'
            WRITE(IU06,*)'+                                           +'
            WRITE(IU06,*)'+ START DATE IN WAMINFO DOES NOT CORRESPOND +'
            WRITE(IU06,*)'+ WITH CURRENT DATE OF THE ATMOSPHERE       +'
            WRITE(IU06,*)'+                                           +'
            WRITE(IU06,*)'+   ABORT SERVICE ROUTINE CALLED BY USERIN  +'
            WRITE(IU06,*)'+   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  +'
            CALL ABORT1
          ENDIF
        ENDIF

        IF(LFDBIOOUT.AND..NOT.LGRIBOUT) THEN
          WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
          WRITE(IU06,*)'+                                          +'
          WRITE(IU06,*)'+ LFDBIOOUT = TRUE AND LGRIBOUT = FALSE    +'
          WRITE(IU06,*)'+ IS AN OBSOLETE OPTION                    +'
          WRITE(IU06,*)'+ PROGRAM WILL ABORT                       +'
          WRITE(IU06,*)'+                                          +'
          WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
          CALL ABORT1
        ENDIF

        CALL FLUSH(IU06)
      ENDIF

      IF(CBPLTDT.GT.CDATEA.AND..NOT.LRESTARTED) THEN
        WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*) '+                                           +'
        WRITE(IU06,*) '+    WARNING IN SUB. USERIN                 +'
        WRITE(IU06,*) '+    ======================                 +'
        WRITE(IU06,*) '+    CBPLTDT GREATER THAN CDATEA            +'
        WRITE(IU06,*) '+    CBPLTDT = ',CBPLTDT
        WRITE(IU06,*) '+    CDATEA  = ',CDATEA
        WRITE(IU06,*) '+    CBPLTDT WAS RESET TO CDATEA            +'
        WRITE(IU06,*) '+                                           +'
        WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++'
        CBPLTDT=CDATEA
      ENDIF

      IF(CDATER .NE. '00000000000000') THEN
        IF(CDATER .LT. CDATEA) CDATER = CDATEE
        IF(CDATER .GT. CDATEE) CDATER = CDATEE
      ENDIF
      IF(CDATES .LT. CDATEA) CDATES = CDATEE
      IF(CDATES .GT. CDATEE) CDATES = CDATEE

!     IF IDELBC IS SET TO ZERO THEN IT WILL BE RESET TO THE LENGTH
!     OF THE RUN (I.E. THE FILE(S) WILL DISPOSED AT THE END OF THE
!     RUN.
      IF(IDELBC.LE.0) CALL DIFDATE(CDATEA,CDATEE,IDELBC)

!     1.3 CHECK IF IDELALT WAS SET ELSE SET IT TO IDELWO
!         ----------------------------------------------

      IF (IDELALT.EQ.0 .AND. IASSI.EQ.1) THEN
        IDELALT=IDELWO
        WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*) '+                                             +'
        WRITE(IU06,*) '+         WARNING IN SUB. USERIN              +'
        WRITE(IU06,*) '+         ======================              +'
        WRITE(IU06,*) '+THE ALTIMETER DATA ASSIMILATION WINDOW LENGTH+'
        WRITE(IU06,*) '+WAS NOT PROVIDED IN WAMINPUT, IT WAS SET TO  +'
        WRITE(IU06,*) '+IDELWO = ', IDELWO
        WRITE(IU06,*) '+                                             +'
        WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
      ENDIF
      

!*    1.4  PRINT INITIAL CONDITIONS.
!          -------------------------
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' WAVE MODEL'
      WRITE(IU06,*) ' '
      WRITE(IU06,*) ' COUPLING WITH ATMOS. MODEL (LWCOU) : ',LWCOU
      WRITE(IU06,*) ' COUPLING WITH NEMO MOD.(LWNEMOCOU) : ',LWNEMOCOU
      IF(.NOT.LMESSPASS) THEN
        WRITE(IU06,*) ' USING ',NTASKS,' PROCESSORS'
        WRITE(IU06,*) ' WITH A WORKLOAD OF ',NSIZE,' PER PROCESSOR'
      ENDIF
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' STARTING DATE ........... : ',CDATEA
      WRITE(IU06,*) ' END DATE ................ : ',CDATEE
      IF(CDATEF.LT.CDATEE) THEN
        WRITE(IU06,*) ' FORECAST STARTING DATE    : ',CDATEF
      ENDIF
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' MODEL TIME STEPS:'
      WRITE(IU06,*) ' SOURCE TERM INTEGRATION TIME STEP : ', &
     & IDELT,' SECS'
      WRITE(IU06,*) ' PROPAGATION TIME STEP ............: ', &
     & IDELPRO,' SECS'
      IF(.NOT.LWCOU) THEN
        IF(NDELW_LST.LE.0) THEN
          WRITE(IU06,*) ' MODEL WIND INPUT TIME STEP .......: ', &
     &     IDELWI,' SECS'
          WRITE(IU06,*) ' MODEL WIND OUTPUT TIME STEP.......: ', &
     &     IDELWO,' SECS'
        ELSE
          DO IC=1,NDELW_LST
            WRITE(IU06,*) ' MODEL WIND INPUT  TIME STEP UNTIL ', &
     &      CDTW_LST(IC),' IS ',IDELWI_LST(IC),' SECS'
            WRITE(IU06,*) ' MODEL WIND OUTPUT TIME STEP UNTIL ', &
     &      CDTW_LST(IC),' IS ',IDELWO_LST(IC),' SECS'
          ENDDO
        ENDIF
      ELSE
        WRITE(IU06,*) ' IFS COUPLING TIME STEP . .........: ',   &
     & KCOUSTEP,' SECS'
      ENDIF
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' MODEL OPTIONS:'
      WRITE(IU06,*) '  '
      IF (LLUNSTR) THEN
        WRITE(IU06,*) ' THE UNSTRUCTED GRID OPTION HAS BEEN SELECTED' 
        WRITE(IU06,*) '   LPREPROC= ', LPREPROC
        WRITE(IU06,*) '   LVECTOR = ', LVECTOR
        WRITE(IU06,*) '   IVECTOR = ', IVECTOR
        WRITE(IU06,*) '   LLUNBINOUT= ',LLUNBINOUT
        WRITE(IU06,*) ' '
      ENDIF
      IF (ICASE.EQ.1) THEN
        WRITE(IU06,*) ' PROPAGATION GRID SPHERICAL LAT/LON COORDINATES'
      ELSE
        WRITE(IU06,*) ' PROPAGATION GRID CARTESIAN COORDINATES'
      ENDIF
      IF(LL1D) THEN
        WRITE(IU06,*) ' 1D DECOMPOSITION OF THE DOMAIN ' 
      ELSE
        WRITE(IU06,*) ' 2D DECOMPOSITION OF THE DOMAIN ' 
      ENDIF
      WRITE(IU06,*) ' MODEL PHYSICS: IPHYS = ', IPHYS
      IF (ISHALLO.EQ.1) THEN
        WRITE(IU06,*) ' THIS IS A DEEP WATER RUN '
      ELSE
        WRITE(IU06,*) ' THIS IS A SHALLOW WATER RUN '
        IF(ISNONLIN.EQ.0) THEN
          WRITE(IU06,*) ' THE OLD SHALLOW WATER SNONLIN IS USED '
        ELSE IF (ISNONLIN.EQ.1) THEN
          WRITE(IU06,*) ' THE NEW SHALLOW WATER SNONLIN IS USED '
        ELSE
          WRITE(IU06,*) ' WARNING: DO NOT KNOW WHICH SNONLIN TO USE !'
          WRITE(IU06,*) ' !!!! THIS ISNONLIN IS NOT YET AVAILABLE !!!'
          WRITE(IU06,*) ' !!!! WARNING: ISNONLIN = ',ISNONLIN
          WRITE(IU06,*) ' '
          CALL ABORT1
        ENDIF
      ENDIF
      WRITE(IU06,*) ' MAXIMUM PHILLIPS PARAMETER ALLOWED: ',ALPHAPMAX
      IF (IDAMPING.EQ.1 .AND. IPHYS.EQ.0) THEN
        WRITE(IU06,*) ' SWELL DAMPING FORMULATION IS USED'
      ENDIF
      IF (ISHALLO.EQ.0) THEN
        IF(LBIWBK) THEN
          WRITE(IU06,*) ' BOTTOM INDUCED WAVE BREAKING IS USED'
        ENDIF
      ENDIF
      IF (IPROPAGS.EQ.1) THEN
        WRITE(IU06,*) ' PROPAGATION: DUAL ROTATED QUADRANTS SCHEME'
!!!! PROPAGS1 has not yet been adapted for anything else then
!!! propgation on a spherical coordinates without currents.
        IF(ICASE.NE.1 .OR.                                              &
     &     (ICASE.EQ.1 .AND. (IREFRA.EQ.2 .OR. IREFRA.EQ.3))            &
     &    ) THEN
          WRITE(IU06,*) ' !!!! WARNING: '
          WRITE(IU06,*) ' !!!! THIS OPTION IS NOT YET AVAILABLE !!!'
          WRITE(IU06,*) ' !!!! IT WILL BE AS IF IPROPAGS=0 ' 
          WRITE(IU06,*) ' '
        ENDIF
      ELSE IF (IPROPAGS.EQ.2) THEN
        WRITE(IU06,*) ' PROPAGATION: CORNER TRANSPORT UPSTREAM SCHEME'
!!!! PROPAGS2 has not yet been adapted for anything else then
        IF(ICASE.NE.1) THEN
          WRITE(IU06,*) ' !!!! WARNING: '
          WRITE(IU06,*) ' !!!! THIS OPTION IS NOT YET AVAILABLE !!!'
          WRITE(IU06,*) ' '
          CALL ABORT1
        ENDIF
      ELSE
        WRITE(IU06,*) ' PROPAGATION: THE SINGLE QUADRANT SCHEME'
      ENDIF
      IF (LSUBGRID) THEN
        WRITE(IU06,*) ' WITH SUB-GRID PARAMETRISATION '
      ELSE
        WRITE(IU06,*) ' WITHOUT SUB-GRID PARAMETRISATION !!! '
      ENDIF
      IF (IREFRA.EQ.0) THEN
        WRITE(IU06,*) ' MODEL RUNS WITHOUT REFRACTION'
      ELSEIF (IREFRA.EQ.1) THEN
        WRITE(IU06,*) ' MODEL RUNS WITH DEPTH REFRACTION ONLY'
      ELSEIF (IREFRA.EQ.2) THEN
        WRITE(IU06,*) ' MODEL RUNS WITH CURRENT REFRACTION ONLY'
        IF(.NOT.LWCOU) THEN
        WRITE(IU06,*) ' WITH A CURRENT INPUT TIME STEP OF ',IDELCUR,    &
     &                ' SECONDS.'
        WRITE(IU06,*) ' STARTING FROM ',CDATECURA
        ENDIF
      ELSEIF (IREFRA.EQ.3) THEN
        WRITE(IU06,*) ' MODEL RUNS WITH DEPTH AND CURRENT REFRACTION'
        IF(.NOT.LWCOU) THEN
        WRITE(IU06,*) ' WITH A CURRENT INPUT TIME STEP OF ',IDELCUR,    &
     &                ' SECONDS.'
        WRITE(IU06,*) ' STARTING FROM ',CDATECURA
        ENDIF
      ENDIF

      IF (LLSOURCE) THEN
        IF(NDELW_LST.LE.0) THEN
          IF (IDELWO.GE.IDELWI) THEN
            WRITE(IU06,*) ' WIND FIELDS ARE NOT INTERPOLATED IN TIME.'
          ELSE
            WRITE(IU06,*) ' WIND FIELDS ARE INTERPOLATED IN TIME.'
          ENDIF
        ELSE
          DO IC=1,NDELW_LST
            IF (IDELWO_LST(IC).GE.IDELWI_LST(IC)) THEN
              WRITE(IU06,*)' WIND FIELDS ARE NOT INTERPOLATED IN TIME', &
     &                     ' UNTIL ', CDTW_LST(IC)
            ELSE
              WRITE(IU06,*)' WIND FIELDS ARE INTERPOLATED IN TIME',     &
     &                     ' UNTIL ', CDTW_LST(IC)
            ENDIF
          ENDDO
        ENDIF
        IF (LLWSWAVE .AND. .NOT. LWCOU) THEN
          WRITE(IU06,*) ' WIND SPEED FROM WAVE MODEL USED AS FORCING.'
        ENDIF
        IF (LLWDWAVE .AND. .NOT. LWCOU) THEN
          WRITE(IU06,*) ' WIND DIRECTION FROM WAVE MODEL USED AS WELL.'
        ENDIF
        IF (LGUST) THEN
          WRITE(IU06,*) ' GUSTINESS EFFECT IS INCLUDED.'
        ENDIF
        IF (LADEN) THEN
          WRITE(IU06,*) ' VARIABLE AIR DENSITY EFFECT IS INCLUDED.'
        ENDIF
        IF ( (LWCOU .AND. LWCUR) .OR.                                   &
     &        IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN
          WRITE(IU06,*) ' SURFACE CURRENTS ARE PROVIDED,'
          IF ( LWCOU ) THEN 
            IF (LRELWIND) THEN
              WRITE(IU06,*) ' THE WINDS ARE RELATIVE TO THE CURRENT.'
            ELSE
              WRITE(IU06,*) ' BUT ABSOLUTE WINDS WILL BE MODIFIED'
              WRITE(IU06,*) ' WITH A FACTOR OF ',-RWFAC 
            ENDIF
          ELSE
            IF (LRELWIND) THEN
              WRITE(IU06,*) ' HENCE, RELATIVE WINDS WILL BE USED, BUT'
              WRITE(IU06,*) ' WITH A REDUCTION FACTOR OF ',RWFAC 
            ELSE
              WRITE(IU06,*) ' BUT ABSOLUTE WINDS WILL BE USED.'
            ENDIF
          ENDIF
        ENDIF
      ELSE
        WRITE(IU06,*) ' ADVECTION ONLY RUN '
        WRITE(IU06,*) ' NO CONTRIBUTION FROM SOURCE TERMS'
      ENDIF

!     WHEN IMPOSING THE ICE MASK SET THRESHOLD TO 0.3
      IF(LMASKICE) THEN
        CITHRSH=0.3_JWRB
        CITHRSH_SAT=CITHRSH
        CIBLOCK=1.0_JWRB
        CITHRSH_TAIL=CITHRSH
        CDICWA=0.0_JWRb
      ELSE
!     RELAX IT A BIT WHEN THE WAVES ARE ALLOWED TO PROPAGATE INTO THE ICE.
        CITHRSH=0.70_JWRB
!     EXCEPT FOR DATA ASSIMILATION
        CITHRSH_SAT=0.01_JWRB
!     BUT ENFORCE FULL BLOCKING CI > CIBLOCK
        CIBLOCK=0.70_JWRB
!     HIGH FREQUENCY SPECTRAL WILL ONLY BE IMPOSED IF SEA ICE COVER <=CITHRSH_TAIL
        CITHRSH_TAIL=0.01_JWRB
!     ICE WATER DRAG COEFFICIENT
        CDICWA=0.01_JWRB
      ENDIF

      IF (LICERUN) THEN
        WRITE(IU06,*) ' VARYING SEA ICE BOUNDARY WILL BE DETERMINED.'
        IF(LMASKICE) THEN
          WRITE(IU06,*) ' AT A FIXED SEA ICE THRESHOLD OF ',CITHRSH
          WRITE(IU06,*) ' WITH FULL BLOCKING THRESHOLD OF ',CIBLOCK
          WRITE(IU06,*) ' AT A FIXED SEA ICE ASSIMLATION THRESHOLD OF ',&
     &                     CITHRSH_SAT
          WRITE(IU06,*) ' WITH A SPECTRAL TAIL SEA ICE THRESHOLD OF ',  &
     &                    CITHRSH_TAIL
        ELSE
          WRITE(IU06,*) ' WITH AN OUPUT SEA ICE THRESHOLD OF ',CITHRSH
          WRITE(IU06,*) ' WITH FULL BLOCKING THRESHOLD OF ',CIBLOCK
          WRITE(IU06,*) ' WITH AN ASSIMILATION SEA ICE THRESHOLD OF ',  &
     &                    CITHRSH_SAT
          WRITE(IU06,*) ' WITH A SPECTRAL TAIL SEA ICE THRESHOLD OF ',  &
     &                    CITHRSH_TAIL
          IF(LCIWABR) THEN
          WRITE(IU06,*) ' WITH AN ICE-WATER DRAG COEFFICENT OF ',CDICWA
          ENDIF
        ENDIF
      ENDIF

      IF (LICETH) THEN
        WRITE(IU06,*) ' SEA ICE THICKNESS IS PART OF THE INPUT.'
      ENDIF

      IF(.NOT.LWCOU) THEN
        IF (IBOUNC.EQ.1) THEN
          WRITE(IU06,*) ' MODEL PRODUCES BOUNDARY DATA (COARSE GRID)'
        ELSE
          IF(ITEST.GT.0) THEN
          WRITE(IU06,*) ' MODEL RUNS WITHOUT BOUNDARY POINTS',          &
     &     ' (COARSE GRID)'
          ENDIF
        ENDIF
        IF (IBOUNF.EQ.1) THEN
          WRITE(IU06,*) ' MODEL RUNS WITH BOUNDARY POINTS (FINE GRID)'
        ELSE
          IF(ITEST.GT.0) THEN
          WRITE(IU06,*) ' MODEL RUNS WITHOUT BOUNDARY POINTS',          &
     &     ' (FINE GRID)'
          ENDIF
        ENDIF
        WRITE(IU06,*) ' '
        IF (IFORCA.EQ.1) THEN
          WRITE(IU06,*) ' MODEL STARTS FROM ANALYSIS FIELDS'
        ELSE
          WRITE(IU06,*) ' MODEL STARTS FROM FORECAST FIELDS'
        ENDIF
      ENDIF
      IF(LNSESTART) THEN
        WRITE(IU06,*) ' INITIAL SPECTRA ARE RESET TO NOISE.' 
      ENDIF
      IF (IASSI.EQ.1) THEN
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' WAVE DATA ASSIMILATION IS CARRIED OUT'
        IF(NASS.GT.0) THEN
          WRITE(IU06,*) ' AT DATE(S) '
          WRITE(IU06,'(2X,A14)') (CASS(I),I=1,NASS)
        ELSE
          IF(LWCOU) THEN
            WRITE(IU06,*) ' AT DATE ', CDATEF
          ELSE
            WRITE(IU06,*) ' UNTIL THE END OF THE ANALYSIS PERIOD'
            WRITE(IU06,*) ' AT DATE ', CDATEF
          ENDIF
          CALL FLUSH(IU06)
        ENDIF
        IF(LALTAS) THEN
          WRITE(IU06,*) ' WITH ALTIMETER DATA IN TIME WINDOW(S) OF '
          WRITE(IU06,*) ' IDELALT = ', IDELALT,' SECONDS'
          IF(IDELALT.LE.10800) THEN
            WRITE(IU06,*) ' ENDING AT ASSIMILATION TIME(S) '
          ELSE
            WRITE(IU06,*) ' CENTERED AROUND THE ASSIMILATION TIME(S) '
          ENDIF
          LLNALTGO = .TRUE.

          MINBUFRSAT=IBUFRSAT(1)
          MAXBUFRSAT=IBUFRSAT(1)
          DO ISAT=2,NUMALT
            MINBUFRSAT=MIN(IBUFRSAT(ISAT),MINBUFRSAT)
            MAXBUFRSAT=MAX(IBUFRSAT(ISAT),MAXBUFRSAT)
          ENDDO
          IF(ALLOCATED(LALTPASSIV)) DEALLOCATE(LALTPASSIV)
          ALLOCATE(LALTPASSIV(MINBUFRSAT:MAXBUFRSAT))
          LALTPASSIV=.FALSE.
          DO ISAT=1,NUMALT
            LALTPASSIV(IBUFRSAT(ISAT))=LALTPAS(ISAT)
          ENDDO

          WRITE(IU06,*) '  '
          DO ISAT=1,NUMALT
             WRITE(IU06,*) ' THE ALTIMETER DATA FROM ',                 &
     &                     TRIM(CSATNAME(ISAT)),' (',IBUFRSAT(ISAT),')'
            IF(LALTPAS(ISAT)) THEN
              WRITE(IU06,*) ' THE DATA WILL ONLY BE USED PASSIVELY !!! '
              WRITE(IU06,*) '  '
            ENDIF
            IF(LALTLRGR(ISAT)) THEN
              WRITE(IU06,*) ' THE DATA WILL BE CORRECTED ' 
              WRITE(IU06,*) ' ACCORDING TO FOLLOWING LINEAR REGRESSION'
              WRITE(IU06,*) ' Hsnew= ',HSCOEFCOR(ISAT),' Hs + ',        &
     &                                 HSCONSCOR(ISAT) 
            ENDIF
            IF(LALTCOR(ISAT)) THEN
              WRITE(IU06,*) ' THE DATA WILL BE CORRECTED ' 
              WRITE(IU06,*) ' ACCORDING TO THE MODEL SEA STATE.' 
            ENDIF
            WRITE(IU06,*) ' THE THRESHOLD FOR BACKGROUND CHECK IS ',    &
     &                      ALTBGTHRSH(ISAT)
            IF(HSALTCUT(ISAT).LT.999999.) THEN
              WRITE(IU06,*) ' THE INPUT MINIMUM WAVE HEIGHT IS ',       &
     &                        HSALTCUT(ISAT)
            ELSE
              WRITE(IU06,*) ' THE MINIMUM WAVE HEIGHT WILL BE',         &
     &                      ' THE OBSERVATION ERROR.'
            ENDIF
            IF(LALTGRDOUT(ISAT)) THEN
              WRITE(IU06,*) ' GRIDDED ALTIMETER FIELDS WILL BE' 
              WRITE(IU06,*) ' PRODUCED FOR THIS ALTIMETER.' 
              LLNALTGO = .FALSE.
            ENDIF
            WRITE(IU06,*) '  '
          ENDDO
          IF (LLNALTGO) THEN
            WRITE(IU06,*) '  '
            WRITE(IU06,*) ' WARNING   WARNING   WARNING   WARNING'
            WRITE(IU06,*) ' GRIDDED ALTIMETER FIELDS WILL NOT BE' 
            WRITE(IU06,*) ' PRODUCED FOR ANY INSTRUMENT.'
            WRITE(IU06,*) '  '
          ENDIF
        ENDIF
        IF(LSARAS) THEN
          WRITE(IU06,*) ' '
          WRITE(IU06,*) ' WITH SAR DATA IN TIME WINDOW(S) OF ' 
          WRITE(IU06,*) ' IDELALT = ', IDELALT,' SECONDS' 
          WRITE(IU06,*) ' CENTERED AROUND THE ASSIMILATION TIME(S) '
        ENDIF
      ELSE
        WRITE(IU06,*) ' DATA ASSIMILATION IS NOT CARRIED OUT'
      ENDIF
      IF(LSARINV) THEN
        WRITE(IU06,*) '  '
        WRITE(IU06,*) ' SAR INVERSION IS CARRIED OUT'
        WRITE(IU06,*) '  '
      ENDIF
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' MODEL OUTPUT SELECTION:'
      WRITE(IU06,*) '  '
      IF (NOUTT.GT.0) THEN
        WRITE(IU06,*) ' NUMBER OF OUTPUT TIMES IS NOUTT = ', NOUTT
        WRITE(IU06,*) ' OUTPUT WILL BE PROCESSED AT:'
        WRITE(IU06,'(6(2X,A14))') (COUTT(I),I=1,NOUTT)
        WRITE(IU06,*) '  '
      ENDIF
      IF(FFLAG20) CALL OUTXT(IU06,IU20,IU30,IU20,IU32,IDELINT,IASSI,    &
     & NOUTT,'IU20','IU30','IU20','IU32')

      IF(LSECONDORDER) THEN
        WRITE(IU06,*) ' SECOND ORDER CORRECTION WILL BE APPLIED TO ',   &
     &                ' OUTPUT INTEGRATED PARAMETERS BASED ON MOMENTS'
        WRITE(IU06,*) '  '
      ENDIF


      WRITE(IU06,*) '###  NAME       OUTPUT OPTION :       FILE  GRIB ',&
     &               ' OUT PE  PARAMID'
      WRITE(IU06,*) '                            F = FALSE   T = TRUE '
      IF(LWAM_USE_IO_SERV) THEN
        WRITE(IU06,*) ''
        WRITE (IU06,*) ' OUTPUT TASK WILL USE THE IFS IOSERVER' 
        WRITE (IU06,*) ' INFORMATION ON OUT PE (below) HAS NO MEANING' 
        WRITE(IU06,*) ''
      ENDIF
      DO ITG = 1,JPPFLAG
        IF(FFLAG(ITG) .OR. GFLAG(ITG)) THEN
        WRITE(CITG,'(I3.3)') ITG
        WRITE(IU06,*) CITG,COUTNAME(ITG),                               &
     &   '...', FFLAG(ITG),'...',GFLAG(ITG),'...',IPFGTBL(ITG),'...',   &
     &   IPRMINFO(ITG,1)*1000+IPRMINFO(ITG,2)
        ENDIF
      ENDDO
      WRITE(IU06,*) ''

      IF(LWCOU .AND. LWFLUX) THEN
        WRITE(IU06,*) ''
        WRITE(IU06,*) ' OCEAN FLUXES WILL ALSO BE RETURNED TO IFS.'
        WRITE(IU06,*) ''
      ENDIF

      IF(LWVFLX_SNL) THEN
        WRITE(IU06,*) ' OCEAN FLUXES WILL INCLUDE SNL CONTRIBUTION'
      ELSE
        WRITE(IU06,*) ' OCEAN FLUXES WILL NOT INCLUDE SNL CONTRIBUTION'
      ENDIF
      WRITE(IU06,*) ''

      WRITE(IU06,*) ' ACCESS TO THE FIELD DATA BASE: '                  &
     & ,'    F = DISABLED   T = ENABLED    ', LFDB
      IF (LFDB) THEN
        ILCFDBSF = LEN_TRIM(CFDBSF)
        IF(ILCFDBSF.GT.0) THEN
          CLFORM='("  ITS ROOT DIRECTORY IS ", A256)'
          WRITE (CLFORM(31:33), '(I3.3)') ILCFDBSF
          WRITE(IU06,CLFORM) CFDBSF(1:ILCFDBSF)
        ELSE
          WRITE(IU06,*)' **********************************************'
          WRITE(IU06,*)'   WARNING : FDB ROOT DIRECTORY NOT SPECIFIED ' 
          WRITE(IU06,*)' **********************************************'
        ENDIF
      ENDIF
      WRITE(IU06,*) '  '


      IF ( LFDB .AND. .NOT. GFLAG20) THEN
        WRITE(IU06,*)' ************************************************'
        WRITE(IU06,*)' *                                              *'
        WRITE(IU06,*)' * ACCESS TO THE FIELD DATA BASE REQUIRES GRIB  *'
        WRITE(IU06,*)' * CODED DATA BUT NO GFLAG WAS SET TO TRUE      *'
        WRITE(IU06,*)' * THIS IS CHANGED AUTOMATICALLY AND THE        *'
        WRITE(IU06,*)' * STANDARD SET OF PARAMETERS IS ENABLED        *'
        WRITE(IU06,*)' * FOR PACKING                                  *'
        WRITE(IU06,*)' *                                              *'
        WRITE(IU06,*)' ************************************************'
        GFLAG(IRHS) = .TRUE.
        GFLAG(IRTP) = .TRUE.
        GFLAG(IRT1) = .TRUE.
        GFLAG(IRCD) = .TRUE.
        GFLAG(IRU10)= .TRUE.
        GFLAG(IRBATHY) = .TRUE.
        GFLAG20 = .TRUE.
      ENDIF
      IF (LFDB) THEN
        WRITE(IU06,*) ' OUTPUT OF GRIB INTEGRATED PARAMETERS REDIRECTED' &
     &   ,' TO THE FIELD DATA BASE'
        WRITE(IU06,*) '                    '
      ENDIF

      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' BINARY RESTART READ  PARALLEL = ', LRSTPARALR 
      WRITE(IU06,*) ' BINARY RESTART WRITE PARALLEL = ', LRSTPARALW 
      WRITE(IU06,*) '  '
      IF (NOUTS.GT.0.AND.IREST.EQ.1) THEN
        WRITE(IU06,*) ' NUMBER OF SPECTRA OUTPUT TIMES IS NOUTS = ', NOUTS
        WRITE(IU06,*) ' SPECTRA OUTPUT WILL BE PROCESSED AT:'
        WRITE(IU06,'(6(2X,A14))') (COUTS(I),I=1,NOUTS)
        WRITE(IU06,*) '  '

      ELSEIF (IREST.EQ.1) THEN
        IF(.NOT.LGRIBOUT) THEN
          WRITE(IU06,*) ' SPECTRA FILES WILL BE WRITTEN OUT TO DISK '   &
     &    ,'EVERY ...', IDELRES, ' SECONDS AND AT THE END OF THE RUN'
        ELSE
          WRITE(IU06,*) ' WAVE SPECTRA WILL BE DISPOSED '               &
     &    ,'EVERY ...', IDELRES
          WRITE(IU06,*) '  AND AT THE END OF THE RUN.'
        ENDIF
        IF(CDATER.LT.CDATEE.AND..NOT.LGRIBOUT) WRITE(IU06,*)            &
     &   ' !! HOWEVER BOTH RESTART FILES WILL ONLY BE SAVED',           &
     &   ' AT ...', CDATER 
        IF(CDATES.LT.CDATEE) WRITE(IU06,*)                              &
     &   ' BUT SPECTRA FILES ALONE  WILL BE SAVED UNTIL '               &
     &   ,'...', CDATES 
      ELSE
        WRITE(IU06,*) ' SPECTRA FILES WILL NOT BE WRITTEN OUT TO DISK'
      ENDIF

      WRITE(IU06,*) '  '
      IF(LNOCDIN.AND.LGRIBIN) THEN
        WRITE(IU06,*) '  '
        WRITE (IU06,*) ' NO DRAG COEFFICIENT FIELD IS PROVIDED AS INPUT'
        WRITE (IU06,*) ' THE FIELD WILL BE INITIALISED BY TAKING' 
        WRITE (IU06,*) ' ZERO WAVE STRESS (TAUW)' 
      ENDIF
      IF(LGRIBIN.AND.ITEST.GT.0) THEN
        WRITE(IU06,*) '  '
        WRITE (IU06,*) ' GRIB SPECTRA FIELD ARE USED AS INPUT'
      ENDIF
      WRITE(IU06,*) '  '
      IF(LFDBIOOUT) THEN
        WRITE (IU06,*) ' FDB SOFTWARE IS USED TO WRITE OUTPUT SPECTRA FILES' 
        IF(LWAM_USE_IO_SERV) THEN
          WRITE (IU06,*) ' OUTPUT TASK WILL USE THE IFS IOSERVER' 
        ELSE
        WRITE (IU06,*) ' OUTPUT TASK SELECTED WITH STRIDE = ',NWRTOUTWAM
        ENDIF
      ELSE
        WRITE (IU06,*) ' OUTPUT SPECTRA FILES ARE WRITTEN OUT TO DISK' 
      ENDIF
!     CFDB2DSP is reset to null since the FDB will be initialised
!     with CFDBSF, provided it is used (ie grib output of integrated
!     parameters.
      IF(LFDB.AND.GFLAG20) CFDB2DSP  = "" 
      IF (LFDBIOOUT) THEN
        CLFORM='("  ITS ROOT DIRECTORY IS ", A256)'
        ILCFDB2DSP = LEN_TRIM(CFDB2DSP)
        IF ( ILCFDB2DSP .EQ. 0 ) THEN
          ILCFDBSF = LEN_TRIM(CFDBSF)
          IF(ILCFDBSF.GT.0) THEN
            WRITE (CLFORM(31:33), '(I3.3)') ILCFDBSF
            WRITE(IU06,CLFORM) CFDBSF(1:ILCFDBSF)
          ENDIF
        ELSE
          WRITE (CLFORM(31:33), '(I3.3)') ILCFDB2DSP
          WRITE(IU06,CLFORM) CFDB2DSP(1:ILCFDB2DSP)
        ENDIF
      ENDIF
      WRITE(IU06,*) '  '
      CALL FLUSH(IU06)

      CALL WSTREAM_STRG(ISTREAM,CSTREAM,NENSFNB,NTOTENS,MARSFCTYPE,     &
     &                  KSTREAM, LASTREAM)

      WRITE(IU06,'("  HARD DRIVE PATH NAME : ",/,5X,A70)') CPATH
      WRITE(IU06,*) '  '
      WRITE(IU06,'("  CURRENT RUN:")')
      WRITE(IU06,'("  GRIB TABLE..: ", I4)') NLOCGRB
      WRITE(IU06,'("  STREAM .....: ", A4)') CSTREAM 
      WRITE(IU06,'("  CLASS.......: ", A4)') YCLASS
      WRITE(IU06,'("  EXPERIMENT..: ", A4)') YEXPVER
      IF( CLDOMAIN == 'g' ) THEN
        WRITE(IU06,'("  MODEL NUMBER: ", I4)') IMDLGRBID_G
      ELSE IF( CLDOMAIN == 'm' ) THEN
        WRITE(IU06,'("  MODEL NUMBER: ", I4)') IMDLGRBID_M
      ELSE IF( CLDOMAIN == 's' ) THEN
        WRITE(IU06,*) '  ONE GRIDPOINT OR SWAMP CASE !!!'
        WRITE(IU06,*) '  SWAMPWIND = ',SWAMPWIND
        WRITE(IU06,*) '  SWAMPWIND2= ',SWAMPWIND2
        WRITE(IU06,*) '  DTNEWWIND = ',DTNEWWIND
        WRITE(IU06,*) '  LTURN90   = ',LTURN90 
        WRITE(IU06,*) '  SWAMPCIFR = ',SWAMPCIFR
        WRITE(IU06,*) '  SWAMPCITH = ',SWAMPCITH
      ENDIF

      IF ( ISTREAM .EQ. 1082 .OR. ISTREAM .EQ. 1095 .OR.                &
     &     ISTREAM .EQ. 1203 .OR. ISTREAM .EQ. 1204 ) THEN
        WRITE(IU06,*) '  '
        IF(ISTREAM.EQ.1095 .OR. ISTREAM.EQ.1203 ) THEN
          WRITE(IU06,'("  MONTHLY FORECAST RUN : ")')
        ELSEIF(ISTREAM.EQ.1204 ) THEN
          WRITE(IU06,'("  MONTHLY FORECAST HINDCAST RUN : ")')
          WRITE(IU06,'("  WITH REFERENCE DATE:      ", I8  )') IREFDATE
        ELSE
          WRITE(IU06,'("  SEASONAL FORECAST RUN : ")')
        ENDIF
        WRITE(IU06,'("  ************************************ ")')
        WRITE(IU06,'("  ENSEMBLE NUMBER:         ", I4  )') NENSFNB
        WRITE(IU06,'("  TOTAL NUMBER OF ENSEMBLE:", I4,/)') NTOTENS
        WRITE(IU06,'("  SYSTEM NUMBER:           ", I4  )') NSYSNB
        WRITE(IU06,'("  METHOD NUMBER:           ", I4  )') NMETNB 
        WRITE(IU06,*) '  '
      ELSE IF ( ISTREAM .EQ. 1083 ) THEN
        WRITE(IU06,'("  MULTI ANALYSIS FORECAST HINDCAST RUN : ")')
        WRITE(IU06,'("  ************************************ ")')
        WRITE(IU06,'("  ENSEMBLE NUMBER:         ", I4  )') NENSFNB
        WRITE(IU06,'("  TOTAL NUMBER OF ENSEMBLE:", I4,/)') NTOTENS
      ELSE IF ( ISTREAM .EQ. 1084 ) THEN
        WRITE(IU06,'("  ENSEMBLE FORECAST HINDCAST RUN : ")')
        WRITE(IU06,'("  WITH REFERENCE DATE:      ", I8  )') IREFDATE
        WRITE(IU06,'("  ENSEMBLE FORECAST RUN : ")')
        WRITE(IU06,'("  ************************************ ")')
        WRITE(IU06,'("  ENSEMBLE NUMBER:         ", I4  )') NENSFNB
        WRITE(IU06,'("  TOTAL NUMBER OF ENSEMBLE:", I4,/)') NTOTENS
      ELSE IF ( ISTREAM .EQ. 1085 ) THEN
        WRITE(IU06,'("  FORECAST HINDCAST RUN : ")')
        WRITE(IU06,'("  WITH REFERENCE DATE:      ", I8  )') IREFDATE
      ELSE IF ( ISTREAM .EQ. 1078 ) THEN
        WRITE(IU06,'("  NEW ENSEMBLE FORECAST HINDCAST RUN OVERLAP: ")')
        WRITE(IU06,'("  WITH REFERENCE DATE:      ", I8  )') IREFDATE
        WRITE(IU06,'("  ENSEMBLE FORECAST RUN : ")')
        WRITE(IU06,'("  ************************************ ")')
        WRITE(IU06,'("  ENSEMBLE NUMBER:         ", I4  )') NENSFNB
        WRITE(IU06,'("  TOTAL NUMBER OF ENSEMBLE:", I4,/)') NTOTENS
      ELSE IF ( ISTREAM .EQ. 1079 ) THEN
        WRITE(IU06,'("  NEW ENSEMBLE FORECAST HINDCAST RUN : ")')
        WRITE(IU06,'("  WITH REFERENCE DATE:      ", I8  )') IREFDATE
        WRITE(IU06,'("  ENSEMBLE FORECAST RUN : ")')
        WRITE(IU06,'("  ************************************ ")')
        WRITE(IU06,'("  ENSEMBLE NUMBER:         ", I4  )') NENSFNB
        WRITE(IU06,'("  TOTAL NUMBER OF ENSEMBLE:", I4,/)') NTOTENS
      ELSE IF ( ISTREAM .EQ. 1088 ) THEN
          WRITE(IU06,'("  ENSEMBLE DATA ASSIMILATION RUN : ")')
          WRITE(IU06,'("  ************************************ ")')
          WRITE(IU06,'("  ENSEMBLE NUMBER:         ", I4  )') NENSFNB
          WRITE(IU06,'("  TOTAL NUMBER OF ENSEMBLE:", I4,/)') NTOTENS
      ELSE
        IF (NENSFNB.NE.0.OR.NTOTENS.NE.0) THEN
          WRITE(IU06,'("  ENSEMBLE FORECAST RUN : ")')
          WRITE(IU06,'("  ************************************ ")')
          WRITE(IU06,'("  ENSEMBLE NUMBER:         ", I4  )') NENSFNB
          WRITE(IU06,'("  TOTAL NUMBER OF ENSEMBLE:", I4,/)') NTOTENS
        ENDIF
      ENDIF

      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' TEST OUTPUT LEVEL IS .............. ITEST = ',    &
     & ITEST
      WRITE(IU06,*) ' TEST OUTPUT OF IN BLOCK LOOPS UPTO ITESTB = ',    &
     & ITESTB
      WRITE(IU06,*) '  '

      IF (.NOT. LWCOU .AND. (LFDB .OR. LFDBIOOUT ) ) THEN
        WRITE(IU06,*) ' ACCESS TO FIELD DATA BASE REQUESTED '
        WRITE(IU06,*) ' DATA BASE WILL BE INITIALIZED '
        WRITE(IU06,*) ' '
        CALL IINITFDBSUBS(INITFDB_STAT)
        IF (INITFDB_STAT .NE. 0) THEN
          WRITE(IU06,'("****************************************")')
          WRITE(IU06,'("*                                      *")')
          WRITE(IU06,'("*    FATAL ERROR IN SUB. USERIN        *")')
          WRITE(IU06,'("*    ==========================        *")')
          WRITE(IU06,'("*    INITIALIZATIO OF THE              *")')
          WRITE(IU06,'("*    FIELD DATA BASE                   *")')
          WRITE(IU06,'("*     F A I L E D .                    *")')
          WRITE(IU06,'("*                                      *")')
          WRITE(IU06,'("*  ERROR CODE (ISTAT) =",i7, "         *")')    &
     &     INITFDB_STAT
          WRITE(IU06,'("*                                      *")')
          WRITE(IU06,'("****************************************")')
          CALL ABORT1
        ENDIF
      ENDIF

      IF(LSMSSIG_WAM) THEN
        WRITE(IU06,*) '  '
        WRITE(IU06,*) ' SIGNALLING IS ACTIVE with '
        WRITE(IU06,*) ' ',TRIM(CMETER),' and ',TRIM(CEVENT)
        WRITE(IU06,*) '  '
      ENDIF

! ----------------------------------------------------------------------

!*    2. CHECK INTEGER RATIOS BETWEEN TIMESTEPS.
!        ---------------------------------------

      LERROR = .FALSE.

!*    2.1 WIND OUTPUT AND PROPAGATION TIME STEP.
!         --------------------------------------

      IF ((IDELWO.LT.IDELPRO.AND.MOD(IDELPRO,IDELWO).NE.0) .OR.         &
     &    (IDELWO.GE.IDELPRO.AND.MOD(IDELWO,IDELPRO).NE.0)) THEN
        WRITE(IU06,*) '*******************************************'
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
        WRITE(IU06,*) '*    ==========================           *'
        WRITE(IU06,*) '* WIND OUTPUT TIMSTEP AND PROPAGATION     *'
        WRITE(IU06,*) '* TIME STEP DO NOT HAVE INTEGER RATIO.    *'
        WRITE(IU06,*) '* WIND OUTPUT TIMSTEP    IDELWO = ', IDELWO
        WRITE(IU06,*) '* PROPAGATION TIME STEP IDELPRO = ', IDELPRO
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '*******************************************'
        LERROR = .TRUE.
      ENDIF

!*    2.2 SOURCE FUNCTION AND PROPAGATION TIMESTEP.
!         -----------------------------------------

      IF (MOD(IDELPRO,IDELT).NE.0 .AND. MOD(IDELT,IDELPRO).NE.0 ) THEN
        WRITE(IU06,*) '*******************************************'
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
        WRITE(IU06,*) '*    ==========================           *'
        WRITE(IU06,*) '* SOURCE FUNCTION  AND PROPAGATION        *'
        WRITE(IU06,*) '* TIME STEP DO NOT HAVE INTEGER RATIO     *'
        WRITE(IU06,*) '* SOURCE FUNCTION TIMESTEP IDELT = ', IDELT
        WRITE(IU06,*) '* PROPAGATION TIMESTEP   IDELPRO = ', IDELPRO
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '*******************************************'
        LERROR = .TRUE.
      ENDIF

!*    2.3 SOURCE FUNCTION AND WIND OUTPUT TIMESTEP.
!         -----------------------------------------

      IF (MOD(IDELWO,IDELT).NE.0 .AND. MOD(IDELT,IDELWO).NE.0 ) THEN
        WRITE(IU06,*) '*******************************************'
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
        WRITE(IU06,*) '*    ==========================           *'
        WRITE(IU06,*) '* THE SOURCE FUNCTION TIME STEP IS NOT    *'
        WRITE(IU06,*) '* AN INTEGER MULTIPLE OF THE WIND OUTPUT  *'
        WRITE(IU06,*) '* TIME STEP !                             *'
        WRITE(IU06,*) '* SOURCE FUNCTION TIMESTEP IDELT = ', IDELT
        WRITE(IU06,*) '* WIND OUTPUT TIMESTEP    IDELWO = ', IDELWO
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '*******************************************'
        LERROR = .TRUE.
      ENDIF

!*    2.4 WIND INPUT AND WIND OUTPUT TIMESTEP.
!         ------------------------------------

      IF(NDELW_LST.LE.0) THEN
        IF (IDELWO.GT.IDELWI) THEN
          WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++'
          WRITE(IU06,*) '+                                         +'
          WRITE(IU06,*) '+   WARNING ERROR IN SUB. USERIN          +'
          WRITE(IU06,*) '+   ============================          +'
          WRITE(IU06,*) '+ WIND INPUT TIME STEP IS LESS THAN       +'
          WRITE(IU06,*) '+ WIND OUTPUT STEP                        +'
          WRITE(IU06,*) '+ WIND INPUT TIMESTEP   IDELWI = ', IDELWI
          WRITE(IU06,*) '+ WIND OUTPUT TIMESTEP  IDELWO = ', IDELWO
          WRITE(IU06,*) '+                                         +'
          WRITE(IU06,*) '+ WIND INPUT CHANGED TO WIND OUTPUT       +'
          WRITE(IU06,*) '+ MODEL WILL USE A NEW WIND FIELD EVERY   +'
          WRITE(IU06,*) '+ WIND OUTPUT TIME STEP AND IGNORE FIELDS +'
          WRITE(IU06,*) '+ IN BETWEEN.                             +'
          WRITE(IU06,*) '+                                         +'
          WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++'
          IDELWI = IDELWO
        ENDIF
        IF ((IDELWO.LE.IDELWI.AND.MOD(IDELWI,IDELWO).NE.0)) THEN
          WRITE(IU06,*) '*******************************************'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
          WRITE(IU06,*) '*    ==========================           *'
          WRITE(IU06,*) '* WIND INPUT  AND WIND OUTPUT             *'
          WRITE(IU06,*) '* TIME STEP DO NOT HAVE INTEGER RATIO OR  *'
          WRITE(IU06,*) '* WIND INPUT TIMESTEP   IDELWI = ', IDELWI
          WRITE(IU06,*) '* WIND OUTPUT TIMESTEP  IDELWO = ', IDELWO
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*******************************************'
          LERROR = .TRUE.
        ENDIF
      ELSE
        DO IC=1,NDELW_LST
          IF (IDELWO_LST(IC).GT.IDELWI_LST(IC)) THEN
            WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++'
            WRITE(IU06,*) '+                                         +'
            WRITE(IU06,*) '+   WARNING ERROR IN SUB. USERIN          +'
            WRITE(IU06,*) '+   ============================          +'
            WRITE(IU06,*) '+ WIND INPUT TIME STEP IS LESS THAN       +'
            WRITE(IU06,*) '+ WIND OUTPUT STEP UNTIL ',CDTW_LST(IC)
            WRITE(IU06,*) '+ WIND INPUT TIMESTEP  = ', IDELWI_LST(IC)
            WRITE(IU06,*) '+ WIND OUTPUT TIMESTEP = ', IDELWO_LST(IC)
            WRITE(IU06,*) '+                                         +'
            WRITE(IU06,*) '+ WIND INPUT CHANGED TO WIND OUTPUT       +'
            WRITE(IU06,*) '+ MODEL WILL USE A NEW WIND FIELD EVERY   +'
            WRITE(IU06,*) '+ WIND OUTPUT TIME STEP AND IGNORE FIELDS +'
            WRITE(IU06,*) '+ IN BETWEEN.                             +'
            WRITE(IU06,*) '+                                         +'
            WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++'
            IDELWI_LST(IC) = IDELWO_LST(IC)
          ENDIF
          IF ((IDELWO_LST(IC).LE.IDELWI_LST(IC).AND.                    &
     &         MOD(IDELWI_LST(IC),IDELWO_LST(IC)).NE.0)) THEN
            WRITE(IU06,*) '*******************************************'
            WRITE(IU06,*) '*                                         *'
            WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
            WRITE(IU06,*) '*    ==========================           *'
            WRITE(IU06,*) '* WIND INPUT AND WIND OUTPUT UNTIL ',        &
     &                       CDTW_LST(IC)
            WRITE(IU06,*) '* TIME STEP DO NOT HAVE INTEGER RATIO OR  *'
            WRITE(IU06,*) '* WIND INPUT TIMESTEP   = ', IDELWI_LST(IC)
            WRITE(IU06,*) '* WIND OUTPUT TIMESTEP  = ', IDELWO_LST(IC)
            WRITE(IU06,*) '*                                         *'
            WRITE(IU06,*) '*******************************************'
            LERROR = .TRUE.
           ENDIF
          IF(IC.GT.1) THEN
            IF(CDTW_LST(IC).LE.CDTW_LST(IC-1)) THEN
              WRITE(IU06,*) '***************************************'
              WRITE(IU06,*) '*                                     *'
              WRITE(IU06,*) '*   FATAL ERROR IN SUB. USERIN        *'
              WRITE(IU06,*) '*   ==========================        *'
              WRITE(IU06,*) '* TIMES FOR WIND INPUT STEPS MUST BE  *'
              WRITE(IU06,*) '* CHRONOLOGICALLY ORDERED IN NAMELIST *' 
              WRITE(IU06,*) '* NAWI !!!!!                          *'
              WRITE(IU06,*) '* IC-1, CDTW_LST(IC-1): ',IC-1,' ',        &
     &                         CDTW_LST(IC-1) 
              WRITE(IU06,*) '* IC, CDTW_LST(IC): ',IC,' ',CDTW_LST(IC) 
              WRITE(IU06,*) '*                                     *'
              WRITE(IU06,*) '***************************************'
              LERROR = .TRUE.
            ENDIF
          ENDIF
        ENDDO
      ENDIF

      CWDFILE='wind_forcing_time_series'
      INQUIRE(FILE=CWDFILE,EXIST=LWDINTS)
      IF(LWDINTS) THEN
        LEN=LEN_TRIM(CWDFILE)
        WRITE(IU06,*) ' '
        WRITE(IU06,*) '  FILE ',CWDFILE(1:LEN),' WAS FOUND.'
        WRITE(IU06,*) '  WIND INPUT WILL BE SET BY IT.'
        WRITE(IU06,*) ' '
!       CHECK THAT THE FILE CONTAINS WINDS GIVEN EVERY IDELWI 
!       (also see readwind)
        IUWDFILE=IWAM_GET_UNIT(IU06,CWDFILE, 'r', 'f', 0)
        OPEN(IUWDFILE,FILE=CWDFILE,FORM='FORMATTED')
        READ(IUWDFILE,*,END=110) IWTIME,WSPEED,WTHETA
        IWTIME_old=IWTIME
        DO WHILE (.TRUE.)
          READ(IUWDFILE,*,END=110) IWTIME,WSPEED,WTHETA
          IF((IWTIME-IWTIME_old).NE.IDELWI) THEN
            WRITE(IU06,*) '****************************************'
            WRITE(IU06,*) '*                                      *'
            WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN        *'
            WRITE(IU06,*) '*    ==========================        *'
            WRITE(IU06,*) '*   THE WIND INPUT FILE DOES NOT HAVE  *'
            WRITE(IU06,*) '*   THE SAME INPUT TIME STEP AS MODEL  *' 
            WRITE(IU06,*) '*   IWTIME-IWTIME_old =',IWTIME-IWTIME_old
            WRITE(IU06,*) '*   IDELWI= ',IDELWI 
            WRITE(IU06,*) '*                                      *'
            WRITE(IU06,*) '****************************************'
            LERROR = .TRUE.
            EXIT 
          ENDIF
          IWTIME_old=IWTIME
        ENDDO
!100     FORMAT(i8,1x,f6.2,1x,f6.1)
110     CLOSE(IUWDFILE)

!       FORCE NO TIME INTERPOLATION
        IDELWO=IDELWI
      ENDIF

!*    2.5 FILE DISPOSE TIMESTEP.
!         ----------------------

      IF (MOD(IDELRES,IDELPRO).NE.0 ) THEN
        WRITE(IU06,*) '*******************************************'
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
        WRITE(IU06,*) '*    ==========================           *'
        WRITE(IU06,*) '* NEW OUTPUT FILES ARE REQUESTED EVERY    *'
        WRITE(IU06,*) '*    IDELRES = ',IDELRES,' SECONDS'
        WRITE(IU06,*) '* IDELRES MUST BE MULTIPLES OF            *'
        WRITE(IU06,*) '* THE PROPAGATION TIMESTEP IDELPRO = ',IDELPRO
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '*******************************************'
        LERROR = .TRUE.
      ENDIF

      IF (MOD(IDELBC,IDELPRO).NE.0 .AND. MOD(IDELBC,IDELWI).NE.0) THEN
        WRITE(IU06,*) '*******************************************'
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
        WRITE(IU06,*) '*    ==========================           *'
        WRITE(IU06,*) '* NEW OUTPUT FILES ARE REQUESTED EVERY    *'
        WRITE(IU06,*) '*    IDELRES = ',IDELBC,' SECONDS'
        WRITE(IU06,*) '* IDELRES MUST BE MULTIPLES OF            *'
        WRITE(IU06,*) '* THE WIND INPUT TIMESTEP   IDELWI = ', IDELWI
        WRITE(IU06,*) '* THE PROPAGATION TIMESTEP IDELPRO = ',IDELPRO
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '*******************************************'
        LERROR = .TRUE.
      ENDIF

!*    2.5 OUTPUT OPTION.
!         --------------

      IF (.NOT. LWAMANOUT) THEN
        WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*) '+ INFO FROM SUB. USERIN                     +'
        WRITE(IU06,*) '+   OUTPUT AT ANALYSIS TIME IS DISABLED.    +'
        WRITE(IU06,*) '+   REASON:  LWAMANOUT IS SET TO .FALSE.    +'
        WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++++'
      ENDIF

      IF (NOUTT.GT.0) THEN
        DO J=1,NOUTT
          CALL DIFDATE (CDATEA, COUTT(J), ISHIFT)
          IF (ISHIFT.LT.0) THEN
            WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++'
            WRITE(IU06,*) '+                                         +'
            WRITE(IU06,*) '+    WARNING INFO IN SUB. USERIN          +'
            WRITE(IU06,*) '+    ============================         +'
            WRITE(IU06,*) '+ OUTPUT DATE IS BEFORE THE STARTING DATE +'
            WRITE(IU06,*) '+ DATE IS : ', COUTT(J)
            WRITE(IU06,*) '+ PROGRAM WILL IGNORE THIS OUTPUT TIME    +'
            WRITE(IU06,*) '+                                         +'
            WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++++'
          ELSE IF (MOD(ISHIFT,IDELPRO).NE.0) THEN
            WRITE(IU06,*) '++++++++++++++++++++++++++++++++++++++++'
            WRITE(IU06,*) '+                                      +'
            WRITE(IU06,*) '+    WARNING ERROR IN SUB. USERIN      +'
            WRITE(IU06,*) '+    ============================      +'
            WRITE(IU06,*) '+ OUTPUT DATE IS NOT AT THE END OF A   +'
            WRITE(IU06,*) '+ PROPAGATION TIMESTEP.                +'
            WRITE(IU06,*) '+ DATE IS : ', COUTT(J)
            WRITE(IU06,*) '+ PROGRAM WILL IGNORE THIS OUTPUT TIME +'
            WRITE(IU06,*) '+                                      +'
            WRITE(IU06,*) '++++++++++++++++++++++++++++++++++++++++'
          ENDIF
        ENDDO
      ELSE
        IF ((FFLAG20.OR.GFLAG20) .AND. IDELINT.EQ.0) THEN
          WRITE(IU06,*) '*******************************************'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
          WRITE(IU06,*) '*    ==========================           *'
          WRITE(IU06,*) '* OUTPUT OF INTEGRATED DATA (TOTAL SEA)   *'
          WRITE(IU06,*) '* IS REQUESTED.                           *'
          WRITE(IU06,*) '* OUTPUT TIME STEP IDELINT                *'
          WRITE(IU06,*) '* HAS TO BE SPECIFIED IN NAMELIST         *'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*******************************************'
          LERROR = .TRUE.
        ENDIF
        IF ((FFLAG20.OR.GFLAG20) .AND. MOD(IDELINT,IDELPRO).NE.0) THEN
          WRITE(IU06,*) '*******************************************'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*    FATAL ERROR IN SUB. USERIN           *'
          WRITE(IU06,*) '*    ==========================           *'
          WRITE(IU06,*) '* OUTPUT OF INTEGRATED DATA (TOTAL SEA)   *'
          WRITE(IU06,*) '* IS REQUESTED.                           *'
          WRITE(IU06,*) '* OUTPUT TIME STEP HAS TO BE A MULTIPLE   *'
          WRITE(IU06,*) '* OF THE PROPAGATION TIME STEP.           *'
          WRITE(IU06,*) '* OUTPUT TIME STEP IS      IDELINT = ',        &
     &     IDELINT
          WRITE(IU06,*) '* PROPAGATION TIME STEP IS IDELPRO = ',        &
     &     IDELPRO
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*******************************************'
          LERROR = .TRUE.
        ENDIF
      ENDIF

      IF (NOUTS.GT.0) THEN
        DO J=1,NOUTS
          CALL DIFDATE (CDATEA, COUTS(J), ISHIFT)
          IF (ISHIFT.LE.0 .OR. MOD(ISHIFT,IDELPRO).NE.0) THEN
            WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++'
            WRITE(IU06,*) '+                                       +'
            WRITE(IU06,*) '+    WARNING ERROR IN SUB. USERIN       +'
            WRITE(IU06,*) '+    ============================       +'
            WRITE(IU06,*) '+ SPECTRA OUTPUT DATE IS NOT AT THE END +'
            WRITE(IU06,*) '+ OF A PROPAGATION TIMESTEP.            +'
            WRITE(IU06,*) '+ DATE IS : ', COUTS(J)
            WRITE(IU06,*) '+ PROGRAM WILL IGNORE THIS OUTPUT TIME  +'
            WRITE(IU06,*) '+                                       +'
            WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++'
          ENDIF
        ENDDO
      ENDIF

      IF (NASS.GT.0.AND.IASSI.EQ.1) THEN
        DO J=1,NASS
          CALL DIFDATE (CDATEA, CASS(J), ISHIFT)
          IF (ISHIFT.LE.0 .OR. MOD(ISHIFT,IDELPRO).NE.0) THEN
            WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++'
            WRITE(IU06,*) '+                                       +'
            WRITE(IU06,*) '+    WARNING ERROR IN SUB. USERIN       +'
            WRITE(IU06,*) '+    ============================       +'
            WRITE(IU06,*) '+ ASSIMILATION DATE IS NOT AT THE END   +'
            WRITE(IU06,*) '+ OF A PROPAGATION TIMESTEP.            +'
            WRITE(IU06,*) '+ DATE IS : ', CASS(J)
            WRITE(IU06,*) '+ PROGRAM WILL ABORT ' 
            WRITE(IU06,*) '+                                       +'
            WRITE(IU06,*) '+++++++++++++++++++++++++++++++++++++++++'
            CALL ABORT1
          ENDIF
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    3. ERROR CHECK.
!        ------------

      IF (LERROR) THEN
        WRITE(IU06,*) '*******************************************'
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '*    FATAL ERROR(S) IN SUB. USERIN        *'
        WRITE(IU06,*) '*    =============================        *'
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '* CORRECT USER INPUT AS INDICATED ABOVE   *'
        WRITE(IU06,*) '* AND TRY AGAIN!!!!!!!!!!!!!!!!!!!!!!!!   *'
        WRITE(IU06,*) '*                                         *'
        WRITE(IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.       *'
        WRITE(IU06,*) '* ---------------   --------------        *'
        WRITE(IU06,*) '*******************************************'
        CALL ABORT1
      ELSE

        IF (LHOOK) CALL DR_HOOK('USERIN',1,ZHOOK_HANDLE)

      ENDIF

! ----------------------------------------------------------------------

      END SUBROUTINE USERIN
