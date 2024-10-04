! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE WAMODEL (NADV, LINIONLY, LFRSTRST, LDSTOP, LDWRRE, BLK2GLO,&
 &                  WVENVI, WVPRPT, FF_NOW, FF_NEXT, INTFLDS,  &
 &                  WAM2NEMO, NEMO2WAM, VARS_4D)

! ----------------------------------------------------------------------

!**** *WAMODEL* - 3-G WAM MODEL - WRAPPER FOR TIME INTEGRATION OF WAVE FIELDS
!                                 AND OUTPUTS

!*    PURPOSE.
!     --------

!       COMPUTATION OF THE 2-D FREQUENCY-DIRECTION WAVE SPECTRUM AT ALL
!       GRID POINTS FOR A GIVEN INITIAL SPECTRUM AND FORCING SURFACE
!       STRESS FIELD.

!**   INTERFACE.
!     ----------

!     *CALL* *WAMODEL (NADV, LINIONLY, LFRSTRST, LDSTOP, LDWRRE, BLK2GLO,
!    &                 WVENVI, WVPRPT, FF_NOW, FF_NEXT, INTFLDS,
!    &                 WAM2NEMO, NEMO2WAM, FL1)
!        *NADV*      NUMBER OF ADVECTION ITERATIONS
!                    PER CALL OF WAMODEL, OUTPUT PARAMETER.
!        *LINIONLY*  INITIALISATION ONLY CALL (i.e. NO FOWARD TIME INTEGRATION) 
!        *LFRSTRST*  FIRST TIME INTEGRATION AFTER RESTART
!        *LDSTOP*    SET .TRUE. IF STOP SIGNAL RECEIVED.
!        *LDWRRE*    SET .TRUE. IF RESTART SIGNAL RECEIVED.
!        *BLK2GLO*   BLOCK TO GRID TRANSFORMATION
!        *WVENVI*    WAVE ENVIRONMENT FIELDS
!        *WVPRPT*    WAVE PROPERTIES FIELDS
!        *FF_NOW*    FORCING FIELDS AT CURRENT TIME.
!        *FF_NEXT*   DATA STRUCTURE WITH THE NEXT FORCING FIELDS
!        *INTFLDS*   INTEGRATED/DERIVED PARAMETERS
!        *WAM2NEMO*  WAVE FIELDS PASSED TO NEMO
!        *NEMO2WAM*  FIELDS FRON OCEAN MODEL to WAM

! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO, ENVIRONMENT, FREQUENCY, FORCING_FIELDS,  &
     &                         INTGT_PARAM_FIELDS, WAVE2OCEAN, OCEAN2WAVE, TYPE_4D, &
                               MIJ_TYPE

      USE YOWCPBO  , ONLY : IBOUNC   ,GBOUNC  , IPOGBO  , CBCPREF
      USE YOWCOUP  , ONLY : LWCOU    ,                                  &
     &                      LWNEMOCOU,                                  &
     &                      NEMOWSTEP, NEMOFRCO     ,                   &
     &                      NEMOCSTEP, NEMONSTEP    , KCOUSTEP
      USE YOWCOUT  , ONLY : COUTT    ,COUTS    ,FFLAG20  ,GFLAG20  ,    &
     &                      NGOUT    ,                                  &
     &                      NIPRMOUT ,                                  &
     &                      LFDB     ,NOUTT    ,NOUTS    ,              &
     &                      CASS     ,NASS     ,LOUTINT  ,              &
     &                      LRSTPARALW, LRSTINFDAT,                     &
     &                      LRSTST0  ,LWAMANOUT
      USE YOWCURR  , ONLY : CDTCUR
      USE YOWFPBO  , ONLY : IBOUNF
      USE YOWFRED  , ONLY : FR       ,TH, WVPRPT_LAND
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE
      USE YOWMESPAS, ONLY : LFDBIOOUT, LGRIBOUT , LNOCDIN, LWAVEWIND 
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,KTAG 
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSTAT  , ONLY : CDATEA   ,CDATEE   ,CDATEF   ,CDTPRO   ,CDTRES   ,    &
     &                      CDATER   ,CDATES   ,CDTINTT  ,IDELPRO  ,IDELT    ,    &
     &                      IDELWI   ,IREST    ,IDELRES  ,IDELINT  ,              &
     &                      CDTBC    ,IDELBC   ,                                  &
     &                      IASSI    ,MARSTYPE ,                                  &
     &                      LLSOURCE ,LANAONLY ,LFRSTFLD ,IREFDATE, LUPDATE_GPU_GLOBALS
      USE YOWSPEC, ONLY   : NBLKS    ,NBLKE, MIJ
      USE YOWTEST  , ONLY : IU06
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH    ,CWI      ,LRESTARTED
      USE YOWUNIT  , ONLY : IU02     ,IU19     ,IU20
      USE YOWWAMI  , ONLY : CBPLTDT  ,CEPLTDT  ,IANALPD  ,IFOREPD  ,    &
     &                      IDELWIN  ,NFCST    ,ISTAT
      USE YOWWIND  , ONLY : CDATEWO

      USE MPL_MODULE, ONLY : MPL_BARRIER
      USE WAM_MULTIO_MOD, ONLY : WAM_MULTIO_FLUSH
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE YOWABORT , ONLY : WAM_ABORT
      USE FIELD_ASYNC_MODULE, ONLY : WAIT_FOR_ASYNC_QUEUE


! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "outwint.intfb.h"
#include "outwpsp.intfb.h"
#include "abort1.intfb.h"
#include "bouinpt.intfb.h"
#include "difdate.intfb.h"
#include "gsfile_new.intfb.h"
#include "headbc.intfb.h"
#include "iwam_get_unit.intfb.h"
#include "incdate.intfb.h"
#include "outbc.intfb.h"
#include "outspec.intfb.h"
#include "outstep0.intfb.h"
#include "savspec.intfb.h"
#include "savstress.intfb.h"
#include "updnemofields.intfb.h"
#include "updnemostress.intfb.h"
#include "writsta.intfb.h"

#ifdef WAM_GPU
#include "outbs_loki_gpu.intfb.h"
#include "wamintgr_loki_gpu.intfb.h"
#else
#include "outbs.intfb.h"
#include "wamintgr.intfb.h"
#endif

      INTEGER(KIND=JWIM), INTENT(IN)                                           :: NADV
      LOGICAL, INTENT(IN)                                                      :: LINIONLY      
      LOGICAL, INTENT(INOUT)                                                   :: LFRSTRST
      LOGICAL, INTENT(INOUT)                                                   :: LDSTOP, LDWRRE
      TYPE(WVGRIDGLO), INTENT(IN)                                              :: BLK2GLO
      TYPE(ENVIRONMENT), INTENT(INOUT)                                         :: WVENVI
      TYPE(FREQUENCY), INTENT(INOUT)                                           :: WVPRPT
      TYPE(FORCING_FIELDS), INTENT(INOUT)                                      :: FF_NOW
      TYPE(FORCING_FIELDS), INTENT(IN)                                         :: FF_NEXT
      TYPE(INTGT_PARAM_FIELDS), INTENT(INOUT)                                  :: INTFLDS
      TYPE(WAVE2OCEAN), INTENT(INOUT)                                          :: WAM2NEMO
      TYPE(OCEAN2WAVE), INTENT(IN)                                             :: NEMO2WAM
      TYPE(TYPE_4D), INTENT(INOUT)                                             :: VARS_4D


      INTEGER(KIND=JWIM) :: IJ, K, M, J, IRA, KADV, ICH
      INTEGER(KIND=JWIM) :: IFIL, IC, ICL, ICR, II, ILOOP
      INTEGER(KIND=JWIM) :: ICHNK
      INTEGER(KIND=JWIM) :: JSTPNEMO, IDATE, ITIME
      INTEGER(KIND=JWIM) :: IU04

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE, ZHOOK_HANDLE_DATA_OFFLOAD, &
      &                    ZHOOK_HANDLE_ADVECTION_LOOP, ZHOOK_HANDLE_IO
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NIPRMOUT, NCHNK) :: BOUT

      CHARACTER(LEN= 2) :: MARSTYPEBAK
      CHARACTER(LEN=14) :: CDATEWH, CZERO
      CHARACTER(LEN=14) :: CDATE, CDTPRA, CDTIMP, CDTIMPNEXT, CDTRCF

      LOGICAL :: LLFLUSH
      LOGICAL :: LSV, LRST, LOUT
      LOGICAL :: LLNONASSI

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WAMODEL',0,ZHOOK_HANDLE)

!     0.0 INITIALISATION
!         --------------

      CZERO = ' '
      LLFLUSH = .FALSE.
      KTAG = 200
      LRSTST0 = .FALSE.

!     TIME FOR THE NEXT SOURCE TERM INTEGRATION
      CDTIMPNEXT = CDTPRO
      CALL INCDATE(CDTIMPNEXT, IDELT)
!     TIME FOR WIND INPUT UPDATE (SEE NEWWIND)
      CDTIMP = CDTPRO

!     0.2 FORECAST STEP 0 IF ANALYSIS IS FOLLOWED BY FORECAST (uncoupled only)
!         --------------------------------------------------------------------
      IF (.NOT.LWCOU .AND. CDTPRO /= CDATEA .AND. CDTPRO == CDATEF) THEN
         CALL OUTSTEP0 (WVENVI, WVPRPT, FF_NOW, INTFLDS,  &
 &                      WAM2NEMO, NEMO2WAM, VARS_4D%FL1, LINIONLY)
      ENDIF

!*    1. ADVECTION/PHYSICS TIME LOOP.
!        ----------------------------

      IF (LHOOK) CALL DR_HOOK('ADVECTION_LOOP',0,ZHOOK_HANDLE_ADVECTION_LOOP)
#ifdef WAM_GPU
      IF (LHOOK) CALL DR_HOOK('DATA_OFFLOAD',0,ZHOOK_HANDLE_DATA_OFFLOAD)
      CALL WVPRPT_LAND%SYNC_DEVICE_RDONLY(QUEUE=0)
      CALL VARS_4D%F_FL1%SYNC_DEVICE_RDWR(QUEUE=0)
      CALL BLK2GLO%SYNC_DEVICE_RDONLY(QUEUE=0)
      CALL WVPRPT%SYNC_DEVICE_RDWR(QUEUE=0)
      CALL WVENVI%SYNC_DEVICE_RDWR(DEPTH=.TRUE., DELLAM1=.TRUE., COSPHM1=.TRUE., UCUR=.TRUE., VCUR=.TRUE., &
      &                            EMAXDPT=.TRUE., IOBND=.TRUE., IODP=.TRUE., QUEUE=0)
      CALL FF_NOW%SYNC_DEVICE_RDWR(AIRD=.TRUE., WDWAVE=.TRUE., CICOVER=.TRUE., WSWAVE=.TRUE.,  &
      & WSTAR=.TRUE., UFRIC=.TRUE., TAUW=.TRUE., TAUWDIR=.TRUE., Z0M=.TRUE., Z0B=.TRUE.,  &
      & CHRNCK=.TRUE., CITHICK=.TRUE., USTRA=.TRUE., VSTRA=.TRUE., QUEUE=1)
      CALL FF_NEXT%SYNC_DEVICE_RDONLY(AIRD=.TRUE., WDWAVE=.TRUE., CICOVER=.TRUE., WSWAVE=.TRUE.,  &
      & WSTAR=.TRUE., UFRIC=.TRUE., TAUW=.TRUE., TAUWDIR=.TRUE., Z0M=.TRUE., Z0B=.TRUE.,  &
      & CHRNCK=.TRUE., CITHICK=.TRUE., USTRA=.TRUE., VSTRA=.TRUE., QUEUE=1)
      CALL WAM2NEMO%SYNC_DEVICE_RDWR(NEMOUSTOKES=.TRUE., NEMOVSTOKES=.TRUE., NEMOSTRN=.TRUE.,  &
      & NPHIEPS=.TRUE., NTAUOC=.TRUE., NSWH=.TRUE., NMWP=.TRUE., NEMOTAUX=.TRUE.,  &
      & NEMOTAUY=.TRUE., NEMOWSWAVE=.TRUE., NEMOPHIF=.TRUE., QUEUE=2)
      CALL INTFLDS%SYNC_DEVICE_RDWR(WSEMEAN=.TRUE., WSFMEAN=.TRUE., USTOKES=.TRUE.,  &
      & VSTOKES=.TRUE., STRNMS=.TRUE., TAUXD=.TRUE., TAUYD=.TRUE., TAUOCXD=.TRUE.,  &
      & TAUOCYD=.TRUE., TAUOC=.TRUE., PHIOCD=.TRUE., PHIEPS=.TRUE., PHIAW=.TRUE., QUEUE=2)
      CALL VARS_4D%F_XLLWS%SYNC_DEVICE_RDWR(QUEUE=2)
      CALL MIJ%SYNC_DEVICE_RDWR(QUEUE=2)
      IF (LHOOK) CALL DR_HOOK('DATA_OFFLOAD',1,ZHOOK_HANDLE_DATA_OFFLOAD)
#endif

      ADVECTION : DO KADV = 1,NADV

!*      1.1 FIX END DATE OF THIS PROPAGATION STEP AND OUTPUT TIMES.
!           -------------------------------------------------------
        CDTPRA = CDTPRO
        CALL INCDATE(CDTPRO, IDELPRO)

!       UPDATE OUTPUT TIMES.
        IF (NOUTT > 0) THEN
          CDTINTT = CZERO
          DO J=1,NOUTT
            IF (CDTPRO == COUTT(J)) THEN
              IF (FFLAG20 .OR. GFLAG20) CDTINTT = COUTT(J)
            ENDIF
          ENDDO
        ELSE
          IF ((FFLAG20.OR.GFLAG20))  THEN
            IF (LWCOU .AND. LRESTARTED .AND. LFRSTRST ) THEN
              LFRSTRST = .FALSE.
              CALL INCDATE (CDTINTT, KCOUSTEP)
            ENDIF
            IF (CDTINTT.LT.CDTPRO) THEN
              CALL INCDATE (CDTINTT, IDELINT)
            ENDIF
          ENDIF
        ENDIF

!       UPDATE SPECTRA OUTPUT DATE
        IF (NOUTS > 0) THEN
!         reset CDATES to insure that spectra output is only controlled
!         by list COUTS
          CDATES='000000000000'
          DO J=1,NOUTS
            IF (CDTPRO == COUTS(J)) THEN
              CDTRES=CDTPRO
              CDATER=CDTRES
              CDATES=CDTRES
              EXIT
            ENDIF
          ENDDO
        ELSE
          IF (CDTRES < CDTPRO) CALL INCDATE(CDTRES, IDELRES)
        ENDIF

!NEST (not used at ECMWF)
        IF ((IBOUNC == 1 .OR. IBOUNF == 1) .AND. CDTBC < CDTPRO) CALL INCDATE(CDTBC, IDELBC)
!NEST



!            -------------------------------------------------------
!*      1.2  THIS IS THE CORE OF THE WAVE MODEL:
!*           COMPUTATION OF PROPAGATION
!*           INTEGRATION OF SOURCE TERMS OVER SUB TIME STEPS BETWEEN
!*           PROPAGATION TIME STEPS.
!            -------------------------------------------------------
!       SET TIME COUNTER.
        CDATE   = CDTPRA
        CDATEWH = CDATEWO
        ILOOP = 1
        DO WHILE ( ILOOP == 1 .OR. CDTIMPNEXT <= CDTPRO)
#ifdef WAM_GPU
          CALL WAMINTGR_LOKI_GPU(CDTPRA, CDATE, CDATEWH, CDTIMP, CDTIMPNEXT, &
 &                       BLK2GLO,                                    &
 &                       WVENVI, WVPRPT, FF_NOW, FF_NEXT, INTFLDS,   &
 &                       WAM2NEMO, MIJ, VARS_4D)
#else
          CALL WAMINTGR (CDTPRA, CDATE, CDATEWH, CDTIMP, CDTIMPNEXT, &
 &                       BLK2GLO,                                    &
 &                       WVENVI, WVPRPT, FF_NOW, FF_NEXT, INTFLDS,   &
 &                       WAM2NEMO, MIJ, VARS_4D)
#endif
          ILOOP = ILOOP +1
        ENDDO


!       1.3 CHECK WHETHER OUTPUT(s) NEEDED
!           ------------------------------
        LRST = (LDWRRE .AND. KADV == NADV )
        IF (LRST) THEN
          WRITE(IU06,*) ' '
          WRITE(IU06,*) '  ******************************************'
          IF (LDSTOP) THEN
          WRITE(IU06,*) '  AN INTERRUPT SIGNAL HAS BEEN RECEIVED '
          ENDIF
          WRITE(IU06,*) '  THE NECESSARY BINARY RESTART FILES WILL BE'
          WRITE(IU06,*) '  GENERATED.'
          WRITE(IU06,*) '  ******************************************'
          WRITE(IU06,*) ' '
          CALL FLUSH (IU06)
        ENDIF


!NEST (not used at ECMWF)
!*      1.4.1 INPUT OF BOUNDARY VALUES.
!           -------------------------
#ifdef _OPENACC
        IF(IBOUNF == 1)THEN
            CALL WAM_ABORT("WAMODEL: IBOUNF==1 NOT SUPPORTED FOR GPU OFFLOAD")
        ENDIF
#endif
        IF (IBOUNF == 1) CALL BOUINPT (IU02, VARS_4D%FL1, NBLKS, NBLKE)
!*      1.4.2 OUTPUT OF BOUNDARY POINTS.
!           --------------------------
        IF (IBOUNC == 1) CALL OUTBC (VARS_4D%FL1, BLK2GLO, IU19)
!NEST

!       1.6 COMPUTE OUTPUT PARAMETERS FIELDS AND PRINT OUT NORMS
!           ----------------------------------------------------
        IF ( (CDTINTT == CDTPRO .OR. LRST) .AND. NIPRMOUT > 0 ) THEN

#ifdef WAM_GPU
          CALL OUTBS_LOKI_GPU (MIJ%PTR, VARS_4D%FL1, VARS_4D%XLLWS, &
     &                WVPRPT, WVENVI, FF_NOW, INTFLDS, NEMO2WAM,   &
     &                BOUT)
#else
          CALL OUTBS (MIJ%PTR, VARS_4D%FL1, VARS_4D%XLLWS, &
     &                WVPRPT, WVENVI, FF_NOW, INTFLDS, NEMO2WAM,   &
     &                BOUT)
#endif

        ENDIF

!*      1.5 POINT OUTPUT (not usually used at ECMWF)
!           ----------------------------------------
        IF ( NGOUT > 0 .AND. (CDTINTT == CDTPRO .OR. LRST) ) THEN
!           OUTPUT POINT SPECTRA (not usually used at ECMWF)
#ifdef WAM_GPU
            IF (LHOOK) CALL DR_HOOK('DATA_OFFLOAD',0,ZHOOK_HANDLE_DATA_OFFLOAD)
            CALL WAIT_FOR_ASYNC_QUEUE(QUEUE=3)
            CALL VARS_4D%F_FL1%GET_HOST_DATA_RDONLY(VARS_4D%FL1)
            !$acc exit data detach(VARS_4D%FL1)
            CALL FF_NOW%GET_HOST_DATA_RDONLY(AIRD=.TRUE., WDWAVE=.TRUE., CICOVER=.TRUE., WSWAVE=.TRUE.,  &
            & WSTAR=.TRUE., UFRIC=.TRUE., TAUW=.TRUE., TAUWDIR=.TRUE., Z0M=.TRUE., Z0B=.TRUE.,  &
            & CHRNCK=.TRUE., CITHICK=.TRUE., USTRA=.TRUE., VSTRA=.TRUE.)
            IF (LHOOK) CALL DR_HOOK('DATA_OFFLOAD',1,ZHOOK_HANDLE_DATA_OFFLOAD)
#endif

            CALL OUTWPSP (VARS_4D%FL1, FF_NOW)
        ENDIF


!*      1.7 ONE PROPAGATION TIMESTEP DONE
!           -----------------------------
        WRITE(IU06,*) ' !!!!!!!!!!!!!! WAVE FIELDS INTEGRATED FOR DATE : ', CDTPRO

!       IS THIS AN ANALYSIS TIME FROM THE INPUT LIST ?
        LLNONASSI=.TRUE.
        IF (IASSI == 1) THEN
          DO J=1,NASS
            IF (CDTPRO == CASS(J)) THEN
              LLNONASSI=.FALSE.
              EXIT
            ENDIF
          ENDDO
        ENDIF


!*      1.8 SAVE FIELD OF 2D SPECTRA AND/OR BINARY RESTART FILES.
!           -----------------------------------------------------
!         THIS WILL HAPPEN WHEN IT IS NOT IN DATA ASSIMILATION MODE AND
!         IT IS EITHER A DETERMINED OUTPUT TIME
!         OR THE INTERUPT SIGNAL HAS BEEN TRIGGERED and it will wait
!         until the end of the advection loop.
!         OTHERWISE THE OUTPUT WILL OCCUR IN WAMASSI.

        LSV=(CDTRES == CDTPRO .OR. CDATEE == CDTPRO .OR. CDTPRO == CDATER)

        IF (LSV .OR. LRST) THEN

          LOUT = ((IREST == 1) .AND. (CDTPRO == CDATER .OR. CDTPRO <= CDATES)) .AND. LSV .AND. LWAMANOUT

          IF ( LOUT .OR. LRST ) THEN

!           1.8.1 SAVE SPECTRUM IN GRIB
!                 ---------------------
            IF (LOUT .AND. LGRIBOUT) THEN
!             we have insured that the spectra will be written to FDB
!             even when the restart option is triggered and it is an
!             output step for the spectra.

!             IF THE OUTPUT TIME IS NOT AN ANALYSIS TIME THEN TYPE FG or 4V
!             BECOMES TYPE AN (i.e. pseudo analysis)
              MARSTYPEBAK=MARSTYPE
              IF ((MARSTYPE == 'fg' .AND. KADV < NADV) .OR.             &
     &            (MARSTYPE == '4v' .AND. LLNONASSI) ) THEN
                MARSTYPE='an'
              ENDIF

#ifdef WAM_GPU
              IF (LHOOK) CALL DR_HOOK('DATA_OFFLOAD',0,ZHOOK_HANDLE_DATA_OFFLOAD)
              CALL WAIT_FOR_ASYNC_QUEUE(QUEUE=3)
              CALL VARS_4D%F_FL1%GET_HOST_DATA_RDONLY(VARS_4D%FL1)
              !$acc exit data detach(VARS_4D%FL1)
              CALL FF_NOW%GET_HOST_DATA_RDONLY(AIRD=.TRUE., WDWAVE=.TRUE., CICOVER=.TRUE., WSWAVE=.TRUE.,  &
              & WSTAR=.TRUE., UFRIC=.TRUE., TAUW=.TRUE., TAUWDIR=.TRUE., Z0M=.TRUE., Z0B=.TRUE.,  &
              & CHRNCK=.TRUE., CITHICK=.TRUE., USTRA=.TRUE., VSTRA=.TRUE.)
              IF (LHOOK) CALL DR_HOOK('DATA_OFFLOAD',1,ZHOOK_HANDLE_DATA_OFFLOAD)
#endif

              CALL OUTSPEC(VARS_4D%FL1, FF_NOW)
              LLFLUSH = .TRUE.

              MARSTYPE=MARSTYPEBAK

              WRITE(IU06,*) ' '
              WRITE(IU06,*) '  GRIB WAVE SPECTRA DISPOSED AT........ CDTPRO  = ', CDTPRO
              WRITE(IU06,*) ' '
            ENDIF

!           1.8.2 SAVE RESTART FILES IN PURE BINARY FORM (in needed)
!                 --------------------------------------
            IF ( .NOT.LGRIBOUT .OR. LDWRRE ) THEN
#ifdef WAM_GPU
              IF (LHOOK) CALL DR_HOOK('DATA_OFFLOAD',0,ZHOOK_HANDLE_DATA_OFFLOAD)
              CALL WAIT_FOR_ASYNC_QUEUE(QUEUE=3)
              CALL WAIT_FOR_ASYNC_QUEUE(QUEUE=4)

              CALL VARS_4D%F_FL1%GET_HOST_DATA_RDONLY(VARS_4D%FL1)
              !$acc exit data detach(VARS_4D%FL1)
              CALL FF_NOW%GET_HOST_DATA_RDONLY(AIRD=.TRUE., WDWAVE=.TRUE., CICOVER=.TRUE., WSWAVE=.TRUE.,  &
              & WSTAR=.TRUE., UFRIC=.TRUE., TAUW=.TRUE., TAUWDIR=.TRUE., Z0M=.TRUE., Z0B=.TRUE.,  &
              & CHRNCK=.TRUE., CITHICK=.TRUE., USTRA=.TRUE., VSTRA=.TRUE.)
              CALL WVENVI%GET_HOST_DATA_RDONLY(DEPTH=.TRUE., DELLAM1=.TRUE., COSPHM1=.TRUE., UCUR=.TRUE., VCUR=.TRUE., &
              &                                EMAXDPT=.TRUE., IOBND=.TRUE., IODP=.TRUE.)
              IF (LHOOK) CALL DR_HOOK('DATA_OFFLOAD',1,ZHOOK_HANDLE_DATA_OFFLOAD)
#endif

              IF (LHOOK) CALL DR_HOOK('IO_TIME',0,ZHOOK_HANDLE_IO)
              CALL SAVSTRESS(WVENVI, FF_NOW, NBLKS, NBLKE, CDTPRO, CDATEF)
              WRITE(IU06,*) ' '
              WRITE(IU06,*) ' BINARY STRESS FILE DISPOSED AT........ CDTPRO  = ', CDTPRO
              WRITE(IU06,*) ' '

              CALL SAVSPEC(VARS_4D%FL1, NBLKS, NBLKE, CDTPRO, CDATEF, CDATER)
              WRITE(IU06,*) '  BINARY WAVE SPECTRA DISPOSED AT........ CDTPRO  = ', CDTPRO
              WRITE(IU06,*) ' '
              CALL FLUSH(IU06)
              IF (LHOOK) CALL DR_HOOK('IO_TIME',1,ZHOOK_HANDLE_IO)
            ENDIF


!*          1.8.3 UPDATE, WRITE AND SAVE WAMINFO FILE.
!                 -----------------------------------
            IF (LRST .AND. IRANK == 1) THEN
              ICH = 7
              CALL DIFDATE (CDATEF, CDATEE, IFOREPD)
              IF (CDTPRO <= CDATEF) THEN
                CALL DIFDATE (CDTPRO, CDATEF, IANALPD)
                CBPLTDT = CDTPRO
                NFCST = 1
              ELSE
                NFCST = 0
                IANALPD = 0
                CBPLTDT = CDATEF
                CALL DIFDATE (CDTPRO, CDATEE, IFOREPD)
              ENDIF
              ISTAT(:) = 0
              IF (CDATE == CDATEE) ISTAT(1) = 1
              IDELWIN = IDELWI

              CEPLTDT = CDATEF

              IU04 = IWAM_GET_UNIT (IU06,CWI(1:ICPLEN+8) , 'w', 'f', 0, 'READWRITE')

              CALL WRITSTA (IU04, CDTPRO, CDATEE, IANALPD, IFOREPD,     &
     &                      IDELWIN, CDATER, CDATES, CBPLTDT, CEPLTDT,  &
     &                      IASSI, NFCST, ISTAT, CDTCUR,                &
     &                      LRSTPARALW, NPROC)

              CLOSE (IU04)
              WRITE(IU06,*) ' WAMINFO FILE WRITTEN FOR RESTART... CDTPRO  = ', CDTPRO
              WRITE(IU06,*) '                                     CDATEF  = ', CDATEF
              WRITE(IU06,*) ' TO ', CWI(1:ICPLEN+8)
              CALL FLUSH(IU06)

              IF (LRSTINFDAT) THEN
!               WRITE AN ADDITIONAL wamfile WITH DATE/TIME INFO added to filename
                CDTRCF=CDTPRO
                CALL INCDATE(CDTRCF, IDELPRO)
                IU04 =  IWAM_GET_UNIT (IU06,CWI(1:ICPLEN+8)//'.'//      &
     &                              CDTRCF(1:8)//'_'//CDTRCF(9:14),     &
     &                              'w', 'f', 0, 'READWRITE')
                CALL WRITSTA (IU04, CDTPRO, CDATEE, IANALPD, IFOREPD,   &
     &                        IDELWIN, CDATER, CDATES, CBPLTDT, CEPLTDT,&
     &                        IASSI, NFCST, ISTAT, CDTCUR,              &
     &                        LRSTPARALW, NPROC)
                CLOSE (IU04)
              ENDIF

            ENDIF  ! waminfo

          ENDIF  ! END SAVE RESTART FIELDS

        ENDIF  ! END SAVE 1.8


!       1.9 WRITE INTEGRATED PARAMETER DATA TO FDB OR TO FILE
!           -------------------------------------------------
        IF (LWAMANOUT .AND. CDTINTT == CDTPRO .AND. NIPRMOUT > 0 ) THEN
!         IF THE OUTPUT TIME IS NOT AN ANALYSIS TIME THEN TYPE FG or 4V
!         BECOMES TYPE AN (i.e. speudo analysis)
          MARSTYPEBAK=MARSTYPE
          IF ((MARSTYPE == 'fg' .AND. KADV < NADV) .OR.                 &
     &        (MARSTYPE == '4v' .AND. LLNONASSI) ) THEN
            MARSTYPE='an'
          ENDIF

          IF (LHOOK) CALL DR_HOOK('IO_TIME',0,ZHOOK_HANDLE_IO)
          CALL OUTWINT(BOUT)
          IF (LHOOK) CALL DR_HOOK('IO_TIME',1,ZHOOK_HANDLE_IO)
          LLFLUSH = .TRUE.

          MARSTYPE=MARSTYPEBAK

          CALL GSTATS(753,0)
          CALL MPL_BARRIER(CDSTRING='WAMODEL:')
          CALL GSTATS(753,1)
        ENDIF


#ifdef WAM_GPU
        IF (LHOOK) CALL DR_HOOK('DATA_OFFLOAD',0,ZHOOK_HANDLE_DATA_OFFLOAD)
        CALL VARS_4D%F_FL1%SYNC_DEVICE_RDWR(QUEUE=0)
        CALL WVENVI%SYNC_DEVICE_RDWR(DEPTH=.TRUE., DELLAM1=.TRUE., COSPHM1=.TRUE., UCUR=.TRUE., VCUR=.TRUE., &
        &                            EMAXDPT=.TRUE., IOBND=.TRUE., IODP=.TRUE., QUEUE=0)
        CALL FF_NOW%SYNC_DEVICE_RDWR(AIRD=.TRUE., WDWAVE=.TRUE., CICOVER=.TRUE., WSWAVE=.TRUE.,  &
        & WSTAR=.TRUE., UFRIC=.TRUE., TAUW=.TRUE., TAUWDIR=.TRUE., Z0M=.TRUE., Z0B=.TRUE.,  &
        & CHRNCK=.TRUE., CITHICK=.TRUE., USTRA=.TRUE., VSTRA=.TRUE., QUEUE=1)
        IF (LHOOK) CALL DR_HOOK('DATA_OFFLOAD',1,ZHOOK_HANDLE_DATA_OFFLOAD)
#endif

!*      1.10 FLUSH FDB IF IT HAS BEEN USED AND IT IS NOT AN ANALYSIS (it will be done in *wamassi*)
!            -------------------------------------------------------

        IF ( LFDB .AND. LLFLUSH .AND. (IASSI /= 1 .OR. CDTPRO > CDATEF) ) THEN
          CALL GSTATS(1976,0)
          CALL WAM_MULTIO_FLUSH()
          CALL GSTATS(1976,1)
          WRITE(IU06,*) ' '
          WRITE(IU06,*) '  FDB FLUSHED AT ',  CDTPRO, ' FROM WAMODEL. '
          CALL FLUSH (IU06)
          LLFLUSH=.FALSE.
        ENDIF


!       1.11 FOR PURE BINARY DATA (obsolete option at ECMWF !!!!!!!):
!*           OUTPUT FILES AND RECOVERY FILES ARE DISPOSED WHEN
!            TIME REACHES THE DISPOSE DATE OR WHEN THE MODEL
!            HAS BEEN SIGNALLED TO DO SO.
!            -------------------------------------------------
        IF (FFLAG20) THEN
          IF (CDATEE == CDTPRO .AND. LOUTINT .AND. IRANK == 1  ) THEN
            CALL GSFILE (IU06, IU20, 0, CDTPRO, CDATEF, 'MAP', 'S')
          ENDIF
        ENDIF

!NEST (not used at ECMWF)
!       SAVE BOUNDARY VALUE FILE.
        IF (CDTBC == CDTPRO) THEN
          IF (IBOUNC == 1 .AND. IRANK == 1 ) THEN
            DO II=1,GBOUNC
            CALL GSFILE(IU06, IU19(II), 0, CDTBC, CDTBC,                &
     &        CBCPREF(II), 'S')
            IF (CDTBC < CDATEE)                                         &
     &        CALL HEADBC (IPOGBO(II)-IPOGBO(II-1), IDELPRO,            &
     &                     TH(1), FR(1), IU19(II), IU06)
            ENDDO
          ENDIF
        ENDIF
!NEST


!*      1.12 WAM-NEMO COUPLING (!!!!! WHEN NO atmospheric model !!!!!!) (currently not used at ECMWF)
!       (when coupled see cnt4 in ifs)
!       ----------------------------------------------------
        IF (LWNEMOCOU .AND. (.NOT.LWCOU)) THEN
          NEMOWSTEP=NEMOWSTEP+1

          IF (MOD(NEMOWSTEP,NEMOFRCO) == 0) THEN

#ifdef WAM_GPU
            IF (LHOOK) CALL DR_HOOK('DATA_OFFLOAD',0,ZHOOK_HANDLE_DATA_OFFLOAD)
            CALL WAIT_FOR_ASYNC_QUEUE(QUEUE=5)
            CALL WAM2NEMO%GET_HOST_DATA_RDONLY(NEMOUSTOKES=.TRUE., NEMOVSTOKES=.TRUE., NEMOSTRN=.TRUE.,  &
            & NPHIEPS=.TRUE., NTAUOC=.TRUE., NSWH=.TRUE., NMWP=.TRUE., NEMOTAUX=.TRUE.,  &
            & NEMOTAUY=.TRUE., NEMOWSWAVE=.TRUE., NEMOPHIF=.TRUE.)
            IF (LHOOK) CALL DR_HOOK('DATA_OFFLOAD',1,ZHOOK_HANDLE_DATA_OFFLOAD)
#endif

            CALL UPDNEMOFIELDS
            CALL UPDNEMOSTRESS

#ifdef WITH_NEMO
            DO JSTPNEMO = NEMOCSTEP, NEMOCSTEP+NEMONSTEP-1
               ! Advance the NEMO model 1 time step
               CALL NEMOGCMCOUP_STEP( JSTPNEMO, IDATE, ITIME )
               WRITE(IU06,*)'NEMO TIME IS : ',JSTPNEMO, IDATE, ITIME
            ENDDO
#endif
            NEMOCSTEP = NEMOCSTEP + NEMONSTEP

#ifdef WAM_GPU
            IF (LHOOK) CALL DR_HOOK('DATA_OFFLOAD',0,ZHOOK_HANDLE_DATA_OFFLOAD)
            CALL WAM2NEMO%SYNC_DEVICE_RDWR(NEMOUSTOKES=.TRUE., NEMOVSTOKES=.TRUE., NEMOSTRN=.TRUE.,  &
            & NPHIEPS=.TRUE., NTAUOC=.TRUE., NSWH=.TRUE., NMWP=.TRUE., NEMOTAUX=.TRUE.,  &
            & NEMOTAUY=.TRUE., NEMOWSWAVE=.TRUE., NEMOPHIF=.TRUE., QUEUE=2)
            IF (LHOOK) CALL DR_HOOK('DATA_OFFLOAD',1,ZHOOK_HANDLE_DATA_OFFLOAD)
#endif
          ENDIF
        ENDIF


        LUPDATE_GPU_GLOBALS = .FALSE.
!*    BRANCHING BACK TO 1.0 FOR NEXT PROPAGATION STEP.
      ENDDO ADVECTION

#ifdef WAM_GPU
      IF (LHOOK) CALL DR_HOOK('DATA_OFFLOAD',0,ZHOOK_HANDLE_DATA_OFFLOAD)
      CALL WVPRPT_LAND%GET_HOST_DATA_RDWR()
      CALL WVPRPT%GET_HOST_DATA_RDWR()
      CALL WVENVI%GET_HOST_DATA_RDWR()
      CALL FF_NOW%GET_HOST_DATA_RDWR()
      CALL FF_NEXT%GET_HOST_DATA_RDWR()
      CALL WAM2NEMO%GET_HOST_DATA_RDWR()
      CALL INTFLDS%GET_HOST_DATA_RDWR()
      CALL VARS_4D%GET_HOST_DATA_RDWR()
      CALL MIJ%GET_HOST_DATA_RDWR()
      CALL BLK2GLO%GET_HOST_DATA_RDWR()

      CALL WVPRPT_LAND%DELETE_DEVICE_DATA()
      CALL WVPRPT%DELETE_DEVICE_DATA()
      CALL WVENVI%DELETE_DEVICE_DATA()
      CALL FF_NOW%DELETE_DEVICE_DATA()
      CALL FF_NEXT%DELETE_DEVICE_DATA()
      CALL WAM2NEMO%DELETE_DEVICE_DATA()
      CALL INTFLDS%DELETE_DEVICE_DATA()
      CALL VARS_4D%DELETE_DEVICE_DATA()
      CALL MIJ%DELETE_DEVICE_DATA()
      CALL BLK2GLO%DELETE_DEVICE_DATA()
      IF (LHOOK) CALL DR_HOOK('DATA_OFFLOAD',1,ZHOOK_HANDLE_DATA_OFFLOAD)
#endif
      IF (LHOOK) CALL DR_HOOK('ADVECTION_LOOP',1,ZHOOK_HANDLE_ADVECTION_LOOP)

IF (LHOOK) CALL DR_HOOK('WAMODEL',1,ZHOOK_HANDLE)

END SUBROUTINE WAMODEL
