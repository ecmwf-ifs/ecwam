SUBROUTINE WAMASSI(LDSTOP, LDWRRE, BLK2GLO, WVENVI,   &
 &                 WVPRPT, FF_NOW, INTFLDS,  &
 &                 WAM2NEMO, NEMO2WAM, FL1)

! ----------------------------------------------------------------------

!****  *WAMASSI* - SUPERVISES EXECUTION OF MAIN MODULES                 
!****              OF THE WAVE DATA ASSIMILATION.                       

!      P. LIONELLO     ECMWF     APRIL   1990

!     PURPOSE.                                                          
!     --------                                                          

!         TO ANALYSE THE WAVE FIELD SUBSTITUTING THE FIRST GUESS        
!         SPECTRA WITH ANALYSED SPECTRA AND FIRST GUESS VALUES          
!         WITH ANALYSED VALUES IN THE GLOBAL GRID                       

!*    INTERFACE.                                                        
!     ----------                                                        

!     *CALL* *WAMASSI(LDSTOP, LDWRRE, BLK2GLO,
!                     WVENVI, WVPRPT, FF_NOW, INTFLDS, 
!                     WAM2NEMO, NEMO2WAM, FL1)
!        *LDSTOP*    SET .TRUE. IF STOP SIGNAL RECEIVED.
!        *LDWRRE*    SET .TRUE. IF RESTART SIGNAL RECEIVED.
!        *BLK2GLO*   BLOCK TO GRID TRANSFORMATION
!        *WVENVI*    WAVE ENVIRONMENT FIELDS
!        *WVPRPT*    WAVE PROPERTIES FIELDS
!        *FF_NOW*    FORCING FIELDS AT CURRENT TIME.
!        *INTFLDS*   INTEGRATED/DERIVED PARAMETERS
!        *WAM2NEMO*  WAVE FIELDS PASSED TO NEMO
!        *NEMO2WAM*  FIELDS FRON OCEAN MODEL to WAM
!        *FL1*       SPECTRUM


!     REFERENCES.                                                       
!     -----------                                                       

!          NONE                                                         

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO, ENVIRONMENT, FREQUENCY, FORCING_FIELDS,  &
     &                         INTGT_PARAM_FIELDS, WAVE2OCEAN, OCEAN2WAVE

      USE YOWCOUP  , ONLY : LWCOU    ,LWNEMOCOUSTRN
      USE YOWCOUT  , ONLY : FFLAG20  ,GFLAG20  ,LFDB     ,LOUTINT  ,             &
     &                      LWAMANOUT, NGOUT   ,NIPRMOUT ,LRSTPARALW
      USE YOWCURR  , ONLY : CDTCUR
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, IJFROMCHNK, KIJL4CHNK 
      USE YOWMESPAS, ONLY : LFDBIOOUT, LGRIBOUT
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : ZMISS    ,DEG
      USE YOWSPEC,   ONLY : NBLKS    ,NBLKE
      USE YOWSTAT  , ONLY : CDATEE   ,CDATEF   ,CDTPRO   ,CDTRES   ,              &
     &                      CDATER   ,CDATES   ,CDTINTT  ,                        &
     &                      IDELWI   ,IREST    ,IREFRA   ,IASSI    ,              &
     &                      NENSFNB  ,NTOTENS  ,NTOTENS  ,NSYSNB   ,CDATEA   ,    &
     &                      MARSTYPE ,YCLASS   ,YEXPVER  ,LALTAS   ,LSARAS   ,    &
     &                      LSARINV  ,IREFDATE
      USE YOWTEST  , ONLY : IU06
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH    ,CWI
      USE YOWWAMI  , ONLY : CBPLTDT  ,CEPLTDT  ,IANALPD  ,IFOREPD  ,    &
     &                      IDELWIN  ,NFCST    ,ISTAT

      USE MPL_MODULE, ONLY : MPL_BARRIER
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE WAM_MULTIO_MOD, ONLY : WAM_MULTIO_FLUSH

! ----------------------------------------------------------------------
      IMPLICIT NONE

#include "outwint.intfb.h"
#include "outwpsp.intfb.h"
#include "altas.intfb.h"
#include "difdate.intfb.h"
#include "gsfile_new.intfb.h"
#include "iwam_get_unit.intfb.h"
#include "outbs.intfb.h"
#include "outint.intfb.h"
#include "outspec.intfb.h"
#include "saras.intfb.h"
#include "savspec.intfb.h"
#include "savstress.intfb.h"
#include "writsta.intfb.h"

      LOGICAL, INTENT(IN) :: LDSTOP, LDWRRE
      TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO 
      TYPE(ENVIRONMENT), INTENT(INOUT) :: WVENVI
      TYPE(FREQUENCY), INTENT(IN) :: WVPRPT
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FF_NOW
      TYPE(INTGT_PARAM_FIELDS), INTENT(INOUT) :: INTFLDS
      TYPE(WAVE2OCEAN), INTENT(INOUT) :: WAM2NEMO
      TYPE(OCEAN2WAVE), INTENT(IN) :: NEMO2WAM
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(INOUT) :: FL1


      INTEGER(KIND=JWIM) :: IFIL, IJ, K, M
      INTEGER(KIND=JWIM) :: IU04
      INTEGER(KIND=JWIM) :: IJSG, IJLG
      INTEGER(KIND=JWIM), DIMENSION(NPROMA_WAM, NCHNK) :: MIJ

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, MAX(NIPRMOUT,1), NCHNK) :: BOUT
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK) :: XLLWS

      LOGICAL :: LLFLUSH
      LOGICAL :: LSV, LOUT

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WAMASSI',0,ZHOOK_HANDLE)


      WRITE (IU06,*) '  '
      WRITE (IU06,*) '  '
      WRITE (IU06,*) ' START OF DATA ASSIMILATION: DATE IS CDTPRO: ', CDTPRO 
      WRITE (IU06,*) '  '
      LLFLUSH = .TRUE.

!     1. CALLING THE MAIN DATA ASSIMILATION ROUTINE
!        ------------------------------------------

      IF (LSARINV) THEN
!!!!!!!!!!!!! inversion of SAR spectra has not happened since ENVISAT (2002-2012).
!!!!!!!!!!!!! This part of the code has not been tested recently
!!!!!!!        CALL SARINVERT(BLK2GLO%IXLG, BLK2GLO%KXLT, FF_NOW%WSWAVE, FF_NOW%WDWAVE, FF_NOW%CICOVER, FL1)
        WRITE(IU06,*) '   SUB. WAMASSI: SARINVERT DONE'
      ENDIF

      IF (LSARAS) THEN
!!!!!!!!!!!!! assimilation SAR spectra has not happened since ENVISAT (2002-2012).
!!!!!!!!!!!!! This part of the code has not been tested recently
        CALL SARAS(BLK2GLO, FF_NOW, FL1)
        WRITE(IU06,*) '   SUB. WAMASSI: SARAS DONE'
      ENDIF

      IF (LALTAS) THEN
        IJSG = IJFROMCHNK(1,1)
        IJLG = IJSG + SUM(KIJL4CHNK) - 1

        CALL ALTAS(IJSG, IJLG, MIJ, BLK2GLO,           &
 &                 WVENVI, FF_NOW, INTFLDS, WAM2NEMO,  &
 &                 WVPRPT, FL1, XLLWS)
        WRITE(IU06,*) '   SUB. WAMASSI: ALTAS DONE'
      ELSE
         MIJ(:,:)=NFRE
         XLLWS(:,:,:,:) = 0.0_JWRB
      ENDIF


!*    2. PREPARING OUTPUT AND UPDATE TO CHARNOCK PARAMETER
!        -------------------------------------------------


!       MODEL OUTPUT INTEGRATED DATA ARE SAVED IN MODULES.
!       --------------------------------------------------
        IF (CDTINTT == CDTPRO) THEN

!         OUTPUT POINT SPECTRA
          IF (NGOUT > 0 ) CALL OUTWPSP (FL1, FF_NOW)

!         COMPUTE OUTPUT PARAMETERS AND PRINT OUT NORMS
          IF (NIPRMOUT > 0) THEN
            CALL OUTBS (MIJ, FL1, XLLWS,                             &
     &                  WVPRPT, WVENVI, FF_NOW, INTFLDS, NEMO2WAM,   &
     &                  BOUT)
          ENDIF
        ENDIF

!*    3.1 SAVE FILES.
!         -----------

!     3.1.1 WRITE INTEGRATED DATA TO FILE AND/OR PRINTER
!           DATA WERE COLLECTED INSIDE THE BLOCK LOOP.
!           --------------------------------------------

      IF (LWAMANOUT .AND. CDTINTT == CDTPRO .AND. NIPRMOUT > 0 ) THEN

        CALL OUTWINT(BOUT)

        CALL GSTATS(752,0)
        CALL MPL_BARRIER(CDSTRING='WAMASSI 1:')
        CALL GSTATS(752,1)
      ENDIF


!*    3.1.4 SAVE RESTART FIELDS.
!           --------------------
      LSV=(CDTRES == CDTPRO .OR. CDATEE == CDTPRO .OR. CDTPRO == CDATER)

      IF (LSV .OR. LDWRRE) THEN
        LOUT=( IREST == 1  .AND. (CDTPRO == CDATER .OR. CDTPRO <= CDATES) .AND. LWAMANOUT )
        IF ( LOUT .OR. LDWRRE ) THEN 

!           SAVE SPECTRUM
!           -------------

          IF ( LOUT .AND. LGRIBOUT ) THEN
!           SAVE SPECTRUM IN GRIB
            CALL OUTSPEC(FL1, FF_NOW)

            CALL GSTATS(1976,0)
            IF (LFDB) CALL WAM_MULTIO_FLUSH()
            CALL GSTATS(1976,1)

            CALL GSTATS(752,0)
            CALL MPL_BARRIER(CDSTRING='WAMASSI 2:')
            CALL GSTATS(752,1)

            LLFLUSH = .FALSE.
            WRITE(IU06,*) ' '
            WRITE(IU06,*) '  GRIB WAVE SPECTRA DISPOSED AT...CDTPRO = ', CDTPRO
            WRITE(IU06,*) ' '
            CALL FLUSH(IU06)
          ENDIF

!         SAVE RESTART FILES IN PURE BINARY FORM
          IF ( .NOT.LGRIBOUT .OR. LDWRRE ) THEN
            CALL SAVSTRESS(WVENVI, FF_NOW, NBLKS, NBLKE, CDTPRO, CDATEF)
            WRITE(IU06,*) '  BINARY STRESS FILE DISPOSED AT... CDTPRO = ', CDTPRO
            WRITE(IU06,*) ' '

            CALL SAVSPEC(FL1, NBLKS, NBLKE, CDTPRO, CDATEF, CDATER)
            WRITE(IU06,*) '  BINARY WAVE SPECTRA DISPOSED AT... CDTPRO  = ', CDTPRO
            WRITE(IU06,*) ' '
            CALL FLUSH(IU06)
          ENDIF


!*    3.1.5 UPDATE, WRITE AND SAVE WAMINFO FILE.
!           ------------------------------------

          IF ((LDSTOP .OR. LDWRRE) .AND.  IRANK == 1 ) THEN
            CALL DIFDATE (CDATEF,CDATEE,IFOREPD)
            CALL DIFDATE (CDTPRO,CDATEF,IANALPD)
            IF (CDTPRO <= CDATEF) THEN
              NFCST = 1
            ELSE
              NFCST = 0
              CALL DIFDATE (CDTPRO,CDATEE,IFOREPD)
            ENDIF
            ISTAT(:) = 0
            IF (CDTPRO == CDATEE) ISTAT(1) = 1
            IDELWIN = IDELWI
            CBPLTDT = CDATEF
            CEPLTDT = CDATEE

            IU04 = IWAM_GET_UNIT (IU06,CWI(1:ICPLEN+8) , 'w', 'f', 0, 'READWRITE')

            CALL WRITSTA (IU04, CDTPRO, CDATEE, IANALPD, IFOREPD,       &
     &                    IDELWIN, CDATER, CDATES, CBPLTDT, CEPLTDT,    &
     &                    IASSI, NFCST, ISTAT, CDTCUR,                  &
     &                    LRSTPARALW, NPROC)

            CLOSE (IU04)
            WRITE(IU06,*) ' WAMINFO FILE WRITTEN FOR RESTART... CDTPRO  = ', CDTPRO
            WRITE(IU06,*) '                                     CDATEF  = ', CDATEF
            WRITE(IU06,*) ' TO ', CWI(1:ICPLEN+8)
          ENDIF
        ENDIF

        CALL GSTATS(752,0)
        CALL MPL_BARRIER(CDSTRING='WAMASSI 3:')
        CALL GSTATS(752,1)
      ENDIF

      IF (LFDB .AND. LLFLUSH .AND. CDTINTT == CDTPRO) THEN
        CALL GSTATS(1976,0)
        CALL WAM_MULTIO_FLUSH()
        CALL GSTATS(1976,1)
        WRITE(IU06,*) '   FDB FLUSHED AT ', CDTPRO, ' FROM WAMASSI. '
        CALL FLUSH (IU06)
      ENDIF

      CALL GSTATS(752,0)
      CALL MPL_BARRIER(CDSTRING='WAMASSI 4:')
      CALL GSTATS(752,1)

IF (LHOOK) CALL DR_HOOK('WAMASSI',1,ZHOOK_HANDLE)

END SUBROUTINE WAMASSI
