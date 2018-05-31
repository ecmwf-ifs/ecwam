      SUBROUTINE WAVEMDL (CBEGDAT, PSTEP, KSTOP, KSTPW,
     &                    NFIELDS, NGPTOTG, NC, NR,
     &                    IGRIB_HANDLE, RMISS, PRCHAR, FIELDS,
     &                    NATMFLX, 
     &                    LWCUR, LWSTOKES,
     &                    NWVFIELDS, WVFLDG,
     &                    NLONW, NLATW, LDSTOP, LDWRRE,
     &                    LDRESTARTED, ZDELATM, KQGAUSS,
     &                    LDWCOUNORMS, MASK_IN, MASK_OUT,
     &                    NFDBREF,
     &                    FRSTIME, NADV, PRPLRADI, PRPLRG,              &
     &                    RNU_ATM, RNUM_ATM,                            &
     &                    IDATE_TIME_WINDOW_END)

!****  *WAVEMDL* - SUPERVISES EXECUTION OF MAIN MODULES                 
!****              OF THE WAVE MODEL                                    

!      LIANA ZAMBRESKY    GKSS/ECMWF    OCTOBER 1988                    

!      MODIFICATION.
!      -------------
!         H. GUNTHER   ECMWF  MARCH    1990                      
!         P. LIONELLO  ECMWF  APRIL    1990  DATA ASSIMILATION   
!                                            MODULE WAMASSI ADDED.
!         J. BIDLOT    ECMWF  FEBRUARY 1996  MESSAGE PASSING.
!         J. DOYLE     ECMWF  OCTOBER  1996  ATMOSPHERIC COUPLING  .
!         J. BIDLOT    ECMWF  FEBRUARY 1997  MESSAGE PASSING.
!         B. HANSEN    ECMWF  MARCH    1997  SIGNAL HANDLING.
!         J. BIDLOT    ECMWF  April    1997  ADD ZDELATM IN PARAM. LIST 
!         S. ABDALLA   ECMWF  OCTOBER  2001  INCLUSION OF AIR DENSITY & Zi/L
!                                            GENERALIZE THE INTERFACE WITH
!                                            THE ATMOSPHERIC MODEL
!         G.Mozdzynski ECMWF  January  2005  OPTIMISE COUPLING COMMS
!         J BIDLOT     ECMWF  August 2006  PASS RMISS TO SPECIFY MISSING DATA.
!         M. Drusch    ECMWF  Sep 2007  Re-initialize through FRSTIME and NADV
!         J BIDLOT     ECMWF  August 2008 ADD LWCUR, MOVE PREWIND.
!         J BIDLOT     ECMWF  June 2009 ADD LWSTOKES. 
!         J BIDLOT     ECMWF  August 2010 ADD IGRIB_HANDLE.
!         P Bechtold   ECMWF  March 2012 ADD SMALL PLANET RADIUS&GRAVITY FACTOR
!                                            

!     PURPOSE.                                                          
!     --------                                                          

!          THIS SUBROUTINE SUPERVISES THE EXECUTION OF                  
!          MAIN MODULES FOR WAM MODEL INITIALIZATION,                   
!          WIND FIELD PREPROCESSING, WAM MODEL EXECUTION,               
!          AND WAVE DATA ASSIMILATION.                                  

!*    INTERFACE.                                                        
!     ----------                                                        

!          SEE MAIN MODULES SUB INITMDL, PREWIND, WAMODEL, WAMASSI.     

!     METHOD.                                                           
!     -------                                                           

!          THE FIRST TIME WAVEMDL IS CALLED, THE WAM MODEL IS           
!          INITIALIZED. THIS INITIALIZATION INCLUDES GETTING            
!          FROM ECFILE THE INITIAL SEA STATE FILES, FILLING             
!          COMMON BLOCKS DEFINING THE GRID AND SETTING GENERAL          
!          PARAMETERS. IN THE FIRST AND ALL SUBSEQUENT CALLS TO         
!          WAVEMDL  PREWIND REFORMATS THE WINDS INTO THE WAM            
!          MODEL BLOCKED STRUCTURE AND THE WAM MODEL IS EXECUTED.       
!          EACH CALL TO WAMODEL INTEGRATES THE WAVE SPECTRA FORWARD     
!          IN TIME BY ONE INPUT WIND TIME STEP OR PROPAGATION TIME      
!          STEP, WHAT EVER IS GREATER.                                  

!     EXTERNALS.                                                        
!     ----------                                                        

!          INITMDL  -  INITIALIZES THE WAM MODEL.                       
!                      GETS RECOVERY FILES OUT OF ECFILE,               
!                      SETS COMMON BLOCKS NECESSARY TO DEFINE           
!                      THE GRID AND BLOCKING STRUCTURE.                 
!                      DEFINES GENERAL PARAMETERS.                      

!          PREWIND  -  REFORMATS WINDS ON THE GAUSSIAN GRID             
!                      INTO THE WAM MODEL BLOCKED STRUCTURE.            

!          WAMODEL  -  INTEGRATES THE WAVE SPECTRA FORWARD IN TIME BY   
!                      ONE WIND INPUT TIME STEP OR ONE PROPAGATION      
!                      TIME STEP, WHATEVER IS GREATER.                  

!          WAMASSI  -  SUPERVISES DATA ASSIMILATION:                    
!                      PREPROCESSES DATA; PRODUCES ANALYSED             
!                      INTEGRATED QUANTITIES BY OPTIMAL INTERPOLATION;  
!                      ANALYSES WAVE SPECTRA ; SAVE ANALYSIS FOR        
!                      OUTPUT AND NEXT TIME STEP MODEL COMPUTATION.     

!     REFERENCES.                                                       
!     -----------                                                       

!          NONE                                                         

! -------------------------------------------------------------------   
      USE PARKIND_WAVE, ONLY : JWIM,JWRB,JWRU
#ifdef MODEL_COUPLING_ATM_WAV
      USE pgmcl_lib_WAM, ONLY : US_coupl, Z0_coupl
#endif
      USE YOWCOUT  , ONLY : CASS     ,NASS
      USE YOWCOUP  , ONLY : LWCOU    ,LWCOU2W  ,LWFLUX   ,LWCOUNORMS,
     &         RNU      ,RNUM        ,
     &         KCOUSTEP, LMASK_OUT_NOT_SET, LMASK_TASK_STR,
     &         I_MASK_OUT, J_MASK_OUT, N_MASK_OUT, LWNEMOCOU, 
     &         LWNEMOCOUSEND, LWNEMOCOURECV, LWNEMOCOUSTK, LWNEMOCOUSTRN
      USE YOWNEMOFLDS, ONLY : NEMOSTRN, NEMOUSTOKES, NEMOVSTOKES
      USE YOWGRIBHD, ONLY : LNEWLVTP,
     &         DATE_TIME_WINDOW_END
      USE YOWGRIB_HANDLES , ONLY : NGRIB_HANDLE_IFS
      USE YOWCURR  , ONLY : IDELCUR  ,LLCHKCFL
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWICE   , ONLY : CICOVER  ,CITHICK  ,CIWA
      USE YOWMAP   , ONLY : ZDELLO   ,IQGAUSS
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,PHIEPS   ,PHIAW   ,
     &            TAUOC
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NANG     ,NFRE    ,
     &                      NBLO     ,LL1D
      USE YOWPCONS , ONLY : ZMISS    ,G
      USE YOWSTAT  , ONLY : MARSTYPE ,CDATEA   ,CDATEE   ,CDATEF   ,
     &            CDTPRO   ,IDELPRO  ,IDELWI   ,IDELWO   ,IASSI    ,
     &            LSMSSIG_WAM,CMETER ,CEVENT   ,
     &            NWFDBREF ,LFDBOPEN ,LSARINV  ,NPROMA_WAM,
     &            IDELWI_LST,IDELWO_LST,CDTW_LST,NDELW_LST
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWNDG  , ONLY : ICODE_CPL
      USE YOWTEXT  , ONLY : LRESTARTED
      USE YOWSPEC, ONLY : U10OLD   ,THWOLD   ,USOLD    ,Z0OLD    ,
     &            TAUW   ,BETAOLD   ,ROAIRO   ,ZIDLOLD  ,
     &            NSTART ,NEND     ,FL1      ,FL3
      USE YOWWIND  , ONLY : CDAWIFL  ,IUNITW ,CDATEWO  ,CDATEFL
      USE FDBSUBS_MOD
      USE GRIB_API_INTERFACE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWUNPOOL ,ONLY : LLUNSTR
      USE MPL_MODULE
! --------------------------------------------------------------------- 

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "cimsstrn.intfb.h"
#include "difdate.intfb.h"
#include "incdate.intfb.h"
#include "initmdl.intfb.h"
#include "mpfldtoifs.intfb.h"
#include "outbeta.intfb.h"
#include "prewind.intfb.h"
#include "setmarstype.intfb.h"
#include "stokesdrift.intfb.h"
#include "wamassi.intfb.h"
#include "wamodel.intfb.h"

!     INITIAL DATE OF THE RUN
      CHARACTER(LEN=14), INTENT(IN) :: CBEGDAT
!     ATMOSPHERIC TIME STEP
      REAL(KIND=JWRB), INTENT(IN) :: PSTEP
!     NUMBER OF ATMOSPHERIC TIME STEPS UNTIL THE END OF THE RUN
      INTEGER(KIND=JWIM), INTENT(IN) :: KSTOP
!     NUMBER OF ATMOSPHERIC TIME STEP(S) FOR EACH WAVE MODEL CALL
      INTEGER(KIND=JWIM), INTENT(IN) :: KSTPW
!     NUMBER OF FIELDS HOLDING ATMOSPHERIC DATA
      INTEGER(KIND=JWIM), INTENT(IN) :: NFIELDS
!     NUMBER OF ATMOSPHERIC GRID POINTS IN THOSE FIELDS
      INTEGER(KIND=JWIM), INTENT(IN) :: NGPTOTG
!     NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR
      INTEGER(KIND=JWIM), INTENT(IN) :: NC
!     NUMBER OF ATM. ROWS OF LATITUDES
      INTEGER(KIND=JWIM), INTENT(IN) :: NR
!     IFS GRIB HANDLE
      INTEGER(KIND=JWIM), INTENT(IN) :: IGRIB_HANDLE
!     GRIB MISSING DATA INDICATOR
      REAL(KIND=JWRB), INTENT(IN) :: RMISS
!     DEFAULT VALUE FOR CHARNOCK
      REAL, INTENT(IN) :: PRCHAR 
!     FIELDS CONTAINING ATMOSPHERIC DATA
      REAL(KIND=JWRB), INTENT(INOUT) :: FIELDS(NGPTOTG,NFIELDS)
!     INDICATES WHICH PHYSICAL QUANTITY IS USED FOR WIND INPUT
!      0: NOT DECIDED WILL BE DETERMINED FROM INPUT
!      1: FRICTION VELOCITY, 2: SURFACE STRESS, 3: 10m WIND
      INTEGER(KIND=JWIM), INTENT(IN) :: NATMFLX
!     INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS FROM IFS
      LOGICAL, INTENT(INOUT) :: LWCUR
!     INDICATES WHETHER THE PRODUCTION OF THE STOKES DRIFT IS REQUIRED FOR THE IFS
      LOGICAL, INTENT(INOUT) :: LWSTOKES
!     NUMBER OF FIELDS RETURNED TO ATMOSPHERIC MODEL
      INTEGER(KIND=JWIM), INTENT(IN) :: NWVFIELDS
!     FIELDS RETURNED TO ATMOSPHERIC MODEL
      REAL(KIND=JWRB), INTENT(OUT) :: WVFLDG(NLONW,NLATW,NWVFIELDS)
!     FIRST DIMENSION OF WVFLDG
      INTEGER(KIND=JWIM), INTENT(IN) :: NLONW
!     SECOND DIMENSION OF WVFLDG
      INTEGER(KIND=JWIM), INTENT(IN) :: NLATW
!     SET .TRUE. IF STOP SIGNAL RECEIVED.
      LOGICAL, INTENT(INOUT) :: LDSTOP
!     SET .TRUE. IF RESTART FILE SIGNAL RECEIVED.
      LOGICAL, INTENT(INOUT) :: LDWRRE
!     TELLS ATMOSPHERIC MODEL THAT IT WAS A RESTART
      LOGICAL, INTENT(OUT) :: LDRESTARTED
!     WAVE MODEL GRID SPACING FOR EACH LATITUDE
      REAL(KIND=JWRB), INTENT(OUT) :: ZDELATM(NLATW)
!     WAVE MODEL GRID IS GAUSSIAN
      INTEGER(KIND=JWIM), INTENT(OUT) :: KQGAUSS
!     TELL ATMOS MODEL THE WAM REQUIREMENT FOR GLOBAL NORMS
      LOGICAL, INTENT(INOUT) :: LDWCOUNORMS
!     MASK TO INDICATE WHICH PART OF ARRAY FIELDS IS RELEVANT 
      INTEGER(KIND=JWIM), INTENT(INOUT) :: MASK_IN(NGPTOTG)
!     MASK TO INDICATE WHICH PART OF ARRAY WVFLDG IS RELEVANT 
      INTEGER(KIND=JWIM), INTENT(INOUT) :: MASK_OUT(NLONW,NLATW)
!     ATMOSPHERIC MODEL FDB ADDRESS 
      INTEGER(KIND=JWIM), INTENT(IN) :: NFDBREF
!     CONTROLS FIRST CALL
      LOGICAL, INTENT(INOUT) :: FRSTIME
!     NUMBER OF ADVECTION STEPS
      INTEGER(KIND=JWIM), INTENT(INOUT) :: NADV
!     MODIFICATION FACTOR FOR EARTH RADIUS FOR SMALL PLANET RUNS
      REAL(KIND=JWRB), INTENT(IN) :: PRPLRADI
!     MODIFICATION FACTOR FOR GRAVITY FOR SMALL PLANET RUNS
      REAL(KIND=JWRB), INTENT(IN) :: PRPLRG
!     KINEMATIC AIR DENSITY
      REAL(KIND=JWRB), INTENT(IN) :: RNU_ATM 
!     REDUCED KINEMATIC AIR DENSITY FOR MOEMNTUM TRANSFER
      REAL(KIND=JWRB), INTENT(IN) :: RNUM_ATM 
!     USED TO SPECIFY THE END OF THE 4DVAr ANALYSIS WINDOW WHEN COUPLED 
      INTEGER(KIND=JWIM), INTENT(IN) :: IDATE_TIME_WINDOW_END

      INTEGER(KIND=JWIM) :: IG
      INTEGER(KIND=JWIM) :: I, J, K, ICPLEN,ICPLEN_ECF
      INTEGER(KIND=JWIM) :: KDELWI, IDURAT
      INTEGER(KIND=JWIM) :: NDUR, KSTOP_BY, ISTOP
      INTEGER(KIND=JWIM) :: N_MASK_OUT_GLOBAL
      INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA
      INTEGER(KIND=JWIM) :: IREAD, ISTAT, IRECV
      INTEGER(KIND=JWIM) :: IYYYYMMDD, IHHMM, ISTEP, IRET, KRET
      INTEGER(KIND=JWIM) :: IC, IL, IST, IED, ICOUNT, JF, IP, IFLD
      INTEGER(KIND=JWIM) :: NCOMBUF, NCOMLOC, NTOT, NMASK
      INTEGER(KIND=JWIM) :: IFCST, IFCSTEP_HOUR
      INTEGER(KIND=JWIM), ALLOCATABLE :: ZCOMCNT(:)

      
      REAL(KIND=JWRB) :: VAL
      REAL(KIND=JWRB) :: STEP
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      REAL(KIND=JWRB) :: DURATION, DURATION_MAX, PSTEP8

      REAL(KIND=JWRB), DIMENSION(NWVFIELDS) :: FAVG,FMIN,FMAX
      REAL(KIND=JWRB), DIMENSION(NWVFIELDS) :: DEFVAL
      REAL(KIND=JWRB), ALLOCATABLE :: ZCOMBUFS(:), ZCOMBUFR(:)
      REAL(KIND=JWRB), ALLOCATABLE :: WVBLOCK(:,:)

#ifdef PARKIND1_SINGLE
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM) ::
     &     ZNEMOUSTOKES,ZNEMOVSTOKES,ZNEMOSTRN
#endif

      CHARACTER(LEN=7), PARAMETER :: CL_CPENV="SMSNAME"
      CHARACTER(LEN=8), PARAMETER :: CL_CPENV_ECF="ECF_NAME"

      CHARACTER(LEN=2) :: CCLASS, CTYPE
      CHARACTER(LEN=4) :: CSTREAM
      CHARACTER(LEN=4) :: CEXPVER
      CHARACTER(LEN=12) :: C12
      CHARACTER(LEN=12) :: FLABEL(NWVFIELDS)
      CHARACTER(LEN=14) :: CDUM
      CHARACTER(LEN=14) :: CDTASS
      CHARACTER(LEN=40) :: CLSETEV
      CHARACTER(LEN=256) :: CLSMSNAME,CLECFNAME

      LOGICAL, SAVE :: LFRST
      LOGICAL, SAVE :: LLGRAPI
      LOGICAL :: L1STCALL
      LOGICAL :: LLGLOBAL_WVFLDG
      LOGICAL :: LLINIT
      LOGICAL :: LLALLOC_FIELDG_ONLY

      DATA LFRST /.TRUE./
      DATA LLGRAPI /.TRUE./

! --------------------------------------------------------------------- 

!*    1.  THE FIRST CALL TO WAVEMDL PERFORMS INITIALIZATION.            
!         --------------------------------------------------  


#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('WAVEMDL',0,ZHOOK_HANDLE)
#endif



      DATE_TIME_WINDOW_END=IDATE_TIME_WINDOW_END
      ZMISS=RMISS
      G=G*PRPLRG ! modified for small planet

      RNU=RNU_ATM
      RNUM=RNUM_ATM

      KDELWI=NINT(PSTEP)*KSTPW

      ICODE_CPL=NATMFLX

      IG=1

      IF(LWCOU) THEN

        NGRIB_HANDLE_IFS=IGRIB_HANDLE

        IF(NGRIB_HANDLE_IFS < 0 ) THEN
          WRITE(IU06,*)' SUB: WAVEMDL:  NGRIB_HANDLE_IFS < 0 !'
          WRITE(IU06,*)' CALL ABORT1 '
          WRITE(IU06,*)'  '
          CALL ABORT1
        ENDIF

        IF (LLGRAPI .OR. ITEST.GE.1) THEN
           WRITE(IU06,*)'  '
           WRITE (IU06,*) ' WAVEMDL: GRIB HANDLE FROM IFS'
           CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'dataDate',IYYYYMMDD)
           CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'time',IHHMM)
           CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'step',STEP)
           CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'endStep',ISTEP)
           CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'expver',C12,KRET=IRET)
           IF(IRET.NE.0) THEN
             CEXPVER='****'
           ELSE
             CEXPVER=C12(1:4)
           ENDIF
           CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'class',C12)
           CCLASS=C12(1:2)
           CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'stream',C12)
           CSTREAM=C12(1:4)
           CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'type',C12)
           CTYPE=C12(1:2)
           WRITE(IU06,*)' EXPVER=', CEXPVER,
     &                  ' CLASS=', CCLASS,
     &                  ' STREAM=', CSTREAM,
     &                  ' TYPE=', CTYPE
           LLGRAPI=.FALSE.
        ENDIF

      ENDIF

      KCOUSTEP=KDELWI

      IREAD=1
      IF(NPROC.EQ.1) IREAD=1

      IF(ALLOCATED(FL3)) DEALLOCATE(FL3)


      IF (FRSTIME) THEN

!     !!!! FIRST TIME AROUND !!!!

        IF (ITEST.GE.1) THEN
          WRITE (IU06,*) ' WAVEMDL:'
          WRITE (IU06,*) ' FRSTIME', FRSTIME
        ENDIF
        CALL FLUSH(IU06)

        LMASK_OUT_NOT_SET=.TRUE.
        LMASK_TASK_STR=.TRUE.

        IF(LWCOU) THEN
          IU06=20
          IL = LEN_TRIM(CBEGDAT)
          IF(IL.NE.14) THEN
            WRITE (IU06,*) ' NON-Y2K COMPLIANT DATE IN CALLING WAVEMDL'
            WRITE (IU06,*) ' CBEGDAT =  ',CBEGDAT
            WRITE (IU06,*) ' PROGRAM WILL ABORT '
            CALL FLUSH(IU06)
            CALL ABORT1
          ENDIF
          CDATEA = CBEGDAT
          CDATEE = CDATEA

          DURATION_MAX=HUGE(KSTOP)
          PSTEP8=PSTEP
          DURATION=REAL(KSTOP,KIND(DURATION))*PSTEP8
          NDUR=INT(DURATION/DURATION_MAX)+1
          KSTOP_BY=KSTOP/NDUR
          DO ISTOP=1,KSTOP,KSTOP_BY
            IDURAT=MIN(KSTOP_BY,KSTOP-ISTOP+1)*NINT(PSTEP)
           CALL INCDATE (CDATEE,IDURAT)
          ENDDO

          WRITE (IU06,1011)
          WRITE (IU06,1012)
          WRITE (IU06,1013)
          IF (LWCOU2W) THEN
            WRITE (IU06,1014)
          ELSE
            WRITE (IU06,1015)
          ENDIF
          WRITE (IU06,1011)
          WRITE (IU06,1016) CDATEA
          WRITE (IU06,1018) CDATEE
          WRITE (IU06,1019) NINT(PSTEP) 
          WRITE (IU06,1020) KSTOP
          WRITE (IU06,1021) KDELWI
          WRITE (IU06,1011)
          WRITE (IU06,1010)

 1010 FORMAT ('  **************************************************')
 1011 FORMAT ('  *                                                *')
 1012 FORMAT ('  *    WAVEMDL: INITMDL STARTS FOR (RE-)START.     *')
 1013 FORMAT ('  *    ========================================    *')
 1014 FORMAT ('  *    TWO-WAY INTERACTION WIND AND WAVES          *')
 1015 FORMAT ('  *    ONE-WAY INTERACTION WIND AND WAVES          *')
 1016 FORMAT ('  * START DATE OF RUN      ', A14,  ' (CDATEA) *')
 1018 FORMAT ('  * END DATE OF RUN        ', A14,  ' (CDATEE) *')
 1019 FORMAT ('  * IFS MODEL TIME STEP    ', I10,  ' (PSTEP)     *')
 1020 FORMAT ('  * TOTAL NUMBER OF PSTEP  ', I10,  ' (KSTOP)     *')
 1021 FORMAT ('  * INTERACTION TIME STEP  ', I10,  ' (KDELWI)    *')

        ENDIF

!       INQUIRE IF IUNITW IS ALREADY OPEN THEN CLOSE IT
        IF(IUNITW.NE.0) CLOSE(IUNITW)

        CALL INITMDL (NADV,
     &                IREAD,
     &                NFIELDS, NGPTOTG, NC, NR,
     &                FIELDS, LWCUR, MASK_IN, PRPLRADI)



        NPROMA=NPROMA_WAM
        CALL GSTATS(1443,0)
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO=IJS(IG),IJL(IG),NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
          CALL OUTBETA (KIJS, KIJL, U10OLD(KIJS,IG),
     &                  USOLD(KIJS,IG), Z0OLD(KIJS,IG),
     &                  BETAOLD(KIJS))
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1443,1)
!

#ifdef ECMWF
        IF(IRANK==1.AND.LSMSSIG_WAM) CALL SIGMASTER()
#endif

        LLCHKCFL=.FALSE.

        KQGAUSS=IQGAUSS
!       if new levtype is used in coupled run than it assumed
!       that WAM and IFS will be sharing the same stream, therefore
!       the FDB should be open on the IFS side.
        IF(LWCOU .AND. LNEWLVTP) THEN
          IF(NFDBREF .GT. 0) THEN
            NWFDBREF=NFDBREF 
            LFDBOPEN=.TRUE.
          ELSE
            NWFDBREF = -5
            LFDBOPEN=.FALSE.
          ENDIF          
        ELSE
!       if you change the -5 value, do the same in wamodel
          NWFDBREF = -5
          LFDBOPEN=.FALSE.
        ENDIF

        FRSTIME = .FALSE.                                              
        L1STCALL = .TRUE.

        IF (ITEST.GE.1) THEN
          WRITE(IU06,*) ' SUB. WAVEMDL: INITMDL DONE'                 
        ENDIF

        IF(LWCOU) THEN
          IF(NLONW.NE.NGX.OR.NLATW.NE.NGY) THEN
            WRITE (IU06,*) ' *********************************'
            WRITE (IU06,*) ' *                               *'
            WRITE (IU06,*) ' * PROBLEM IN WAVEMDL..........  *'
            WRITE (IU06,*) ' * PROBLEM WITH NLONW AND NLATW  *'
            WRITE (IU06,*) ' * NOT EQUAL TO NGX   AND NGY  : *'
            WRITE (IU06,*) ' * ============================= *'
            WRITE (IU06,*) ' *                               *'
            WRITE (IU06,*) ' * NLONW=',NLONW
            WRITE (IU06,*) ' * NLATW=',NLATW
            WRITE (IU06,*) ' * NGX=',NGX
            WRITE (IU06,*) ' * NGY=',NGY
            WRITE (IU06,*) ' *                               *'
            WRITE (IU06,*) ' *                               *'
            WRITE (IU06,*) ' *********************************'
            CALL ABORT1
          ENDIF
          IF (LRESTARTED)  THEN
            LDRESTARTED = .TRUE.
!            RETURN
          ELSE
            LDRESTARTED = .FALSE.
          ENDIF
        ENDIF

      ELSE

!     !!!! ANY OTHER TIMES !!!!

!       CHECK IF THE INTERACTION/WIND-INPUT TIME STEP HAS NOT CHANGED
!       ALSO FOR SURFACE CURRENT INPUT,
!       OTHERWISE RESET NADV
        IF(LWCOU) THEN
          IF(IDELWI.NE.KDELWI) THEN
            IDELWO = KDELWI
            IDELWI = KDELWI
            IDELCUR= KDELWI
            NADV = IDELWI/IDELPRO
            CDAWIFL=CDTPRO
            CALL INCDATE(CDAWIFL,IDELWI)
            WRITE (IU06,*) ' ***********************************'
            WRITE (IU06,*) ' *                                 *'
            WRITE (IU06,*) ' * IN WAVEMDL :                    *'
            WRITE (IU06,*) ' * INTERACTION TIME STEP WAS RESET *' 
            WRITE (IU06,*) ' * KDELWI = ', KDELWI 
            WRITE (IU06,*) ' *                                 *'
            WRITE (IU06,*) ' ***********************************'
            CALL FLUSH(IU06)
            IF(IDELPRO*NADV.NE.IDELWI) THEN
              WRITE (IU06,*) ' ***************************************'
              WRITE (IU06,*) ' *                                     *'
              WRITE (IU06,*) ' * PROBLEM IN WAVEMDL :                *'
              WRITE (IU06,*) ' * THE NEW INTERACTION TIME STEP IS NOT*'
              WRITE (IU06,*) ' * A MULTIPLE OF IDELPRO  !!!          *'
              WRITE (IU06,*) ' * KDELWI = ', KDELWI 
              WRITE (IU06,*) ' * IDELPRO = ', IDELPRO 
              WRITE (IU06,*) ' *                                     *'
              WRITE (IU06,*) ' ***************************************'
              CALL ABORT1
            ENDIF
          ENDIF
        ELSE
          IF(NDELW_LST.GT.0) THEN
            DO IC=1,NDELW_LST
              IF(CDTPRO.LT.CDTW_LST(IC)) THEN
                IF(IDELWI.NE.IDELWI_LST(IC) .OR.
     &             IDELWO.NE.IDELWO_LST(IC)) THEN

                  CALL INCDATE(CDATEWO,-IDELWO/2)
                  IDELWI=IDELWI_LST(IC)
                  IDELWO=IDELWO_LST(IC)
                  CALL INCDATE(CDATEWO,IDELWO/2)
                  CDATEFL=CDATEWO
                  NADV = IDELWI/IDELPRO
                  CDAWIFL=CDTPRO
                  CALL INCDATE(CDAWIFL,IDELWI)
                  WRITE (IU06,*) ' **********************************'
                  WRITE (IU06,*) ' *                                *'
                  WRITE (IU06,*) ' * IN WAVEMDL :                   *'
                  WRITE (IU06,*) ' * WIND INPUT TIME STEP WAS RESET *' 
                  WRITE (IU06,*) ' * IDELWI = ', IDELWI 
                  WRITE (IU06,*) ' * IDELWO = ', IDELWO 
                  WRITE (IU06,*) ' *                                *'
                  WRITE (IU06,*) ' **********************************'
                  CALL FLUSH(IU06)
                ENDIF

                EXIT

              ENDIF
            ENDDO
          ENDIF

        ENDIF

!       KEEP ATMOSPHERIC MODEL INFORMED ABOUT OUR GRID
        KQGAUSS=IQGAUSS

        IF(LWCOU .AND. LNEWLVTP .AND. NFDBREF .GT. 0) THEN
          NWFDBREF=NFDBREF 
          LFDBOPEN=.TRUE.
        ENDIF

        L1STCALL = .FALSE.

!*      REFORMAT FORCING FIELDS FROM INPUT GRID TO BLOCKED.                     
!       ---------------------------------------------------
        LLINIT=.FALSE.
        LLALLOC_FIELDG_ONLY=LWCOU
!       !!!! PREWIND IS CALLED THE FIRST TIME IN INITMDL !!!!
        CALL PREWIND (U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     &              ROAIRO, ZIDLOLD,
     &              CICOVER, CITHICK, CIWA,
     &              LLINIT, LLALLOC_FIELDG_ONLY,
     &              IREAD,
     &              NFIELDS, NGPTOTG, NC, NR,
     &              FIELDS, LWCUR, MASK_IN)


        IF (ITEST.GE.1) THEN
          WRITE(IU06,*) ' SUB. WAVEMDL: PREWIND DONE'                   
        ENDIF

      ENDIF

!!test: reallocating the memory that was freed

      IF(.NOT.ALLOCATED(FL3)) ALLOCATE(FL3(NINF-1:NSUP,NANG,NFRE))


! --------------------------------------------------------------------  

!*    2.1  INTEGRATE THE WAVE SPECTRA FORWARD IN TIME.                  
!          -------------------------------------------                  

      CALL SETMARSTYPE

      CALL WAMODEL (NADV, LDSTOP, LDWRRE, L1STCALL)

      IF (ITEST.GE.1) THEN
        WRITE(IU06,*) ' SUB. WAVEMDL: WAMODEL DONE'
      ENDIF


!*    2.2  DATA ASSIMILATION AND/OR SAR INVERSION.
!*         IF REQUESTED AND MODEL IS IN ANALYSIS PERIOD.
!          THE DATA ASSIMILATION SOFTWARE IS NOT AVAILABLE
!          FOR GENERAL DISSIMINATION !
!          ---------------------------------------------

#ifdef ECMWF
      IF (IASSI.EQ.1) THEN

        MARSTYPE = 'an'

!       UPDATE ANALYSIS TIME
        DO J=1,NASS
          IF(CDTPRO.EQ.CASS(J)) THEN
            CDTASS=CDTPRO
            EXIT
          ENDIF
        ENDDO

        IF (NASS.GT.0 ) THEN
          IF ( CDTPRO.EQ.CDTASS ) THEN
            CALL WAMASSI (LDSTOP, LDWRRE)
            IF (ITEST.GE.1) THEN
              WRITE(IU06,*) ' SUB. WAVEMDL : WAMASSI DONE'
            ENDIF
          ENDIF
        ELSEIF ((.NOT.LWCOU .AND. CDTPRO.LE.CDATEF ) 
     &      .OR. 
     &      (LWCOU .AND. CDTPRO.EQ.CDATEF)) THEN
          CALL WAMASSI (LDSTOP, LDWRRE)
          IF (ITEST.GE.1) THEN
            WRITE(IU06,*) ' SUB. WAVEMDL: WAMASSI DONE'
          ENDIF
        ENDIF
      ELSE IF (IASSI.NE.1 .AND. LSARINV) THEN

!!!!
!!!!! we may have to introduce a list for sar inversion time !!!!
!       UPDATE ANALYSIS TIME (which is now the inversion time)
        DO J=1,NASS
          IF(CDTPRO.EQ.CASS(J)) THEN
            CDTASS=CDTPRO
            EXIT
          ENDIF
        ENDDO

        IF (NASS.GT.0 ) THEN
          IF ( CDTPRO.EQ.CDTASS ) THEN
            CALL SARINVERT
            IF(ITEST.GE.1)
     &         WRITE(IU06,*) '   SUB. WAVEMDL: SARINVERT DONE'
          ENDIF
        ELSEIF ((.NOT.LWCOU .AND. CDTPRO.LE.CDATEF ) 
     &      .OR. 
     &      (LWCOU .AND. CDTPRO.EQ.CDATEF)) THEN
          CALL SARINVERT
          IF(ITEST.GE.1) WRITE(IU06,*) '   SUB. WAVEMDL: SARINVERT done'
        ENDIF
      ENDIF
#endif


#ifdef ECMWF
!     SMS METER:
      IF(LSMSSIG_WAM) THEN

        CALL MPL_BARRIER(CDSTRING='WAVEMDL:')

        CALL DIFDATE (CDATEF, CDTPRO, IFCST)
        IFCSTEP_HOUR=IFCST/3600
        IF(IRANK==1) THEN
          WRITE(CLSETEV,' (A20,'' step '',I5,''&'') ') 
     &        CMETER,IFCSTEP_HOUR
          CLSMSNAME="                                             "
          CLECFNAME="                                             "
          CALL UTIL_CGETENV(CL_CPENV,'NOSMS',CLSMSNAME, ICPLEN)
          CALL UTIL_CGETENV(CL_CPENV_ECF,'NOECF',CLECFNAME,ICPLEN_ECF)
          IF ((ICPLEN > 0.AND.CLSMSNAME(1:5) /= 'NOSMS') .OR. 
     &        (ICPLEN_ECF > 0.AND.CLECFNAME(1:5) /= 'NOECF') ) THEN
            CALL SYSTEM(CLSETEV)
            WRITE(IU06,'(2X,A20,I5,'' posted '')')
     &       CMETER,IFCSTEP_HOUR
          ELSE
            WRITE(IU06,'(A20,I5)') CMETER,IFCSTEP_HOUR
            WRITE(IU06,*) 'not posted  because neither SMSNAME'
            WRITE(IU06,*) ICPLEN, CLSMSNAME
            WRITE(IU06,*) 'nor ECF_NAME  is defined. '
            WRITE(IU06,*) ICPLEN_ECF, CLECFNAME  
          ENDIF
        ENDIF
      ENDIF
#endif



!----------------------------------------------------------------------

!     3. PREPARE FIELDS THAT ARE RETURNED TO IFS (if coupled).
!        -----------------------------------------------------

      IF(LWCOU) THEN


        ALLOCATE(WVBLOCK(IJS(IG):IJL(IG),NWVFIELDS))

!       GRID LAYOUT
        DO K=1,NLATW
          ZDELATM(K) = ZDELLO(K)
        ENDDO

!       FIELDS TO BE PASSED TO THE ATMOSPHERIC MODEL ARE:

!       1. CHARNOCK FIELD
        FLABEL(1)=' Charnock'
        DEFVAL(1)=PRCHAR ! DEFAULT VALUE FOR GRID POINTS NOT COVERED BY
                         ! THE WAVE MODEL ICE FREE SEA POINTS.
                         !  !!! NO LONGER USED OVER SEA ICE SEE *BETADEF*
                         ! BUT WELL OVER LAND POINTS.

!       2. U-STOKESDRIFT
        IF(NWVFIELDS.GT.1) THEN
          FLABEL(2)=' U-Stokes'
          DEFVAL(2)=0.0  ! DEFAULT VALUE FOR GRID POINTS NOT COVERED BY
                         ! THE WAVE MODEL ICE FREE SEA POINTS.
        ENDIF

!       3. V-STOKESDRIFT
        IF(NWVFIELDS.GT.2) THEN
          FLABEL(3)=' V-Stokes'
          DEFVAL(3)=0.0  ! DEFAULT VALUE FOR GRID POINTS NOT COVERED BY
                         ! THE WAVE MODEL ICE FREE SEA POINTS.
        ENDIF

!       4. ENERGY FLUX TO OCEAN
        IF(NWVFIELDS.GT.3) THEN
          FLABEL(4)=' Phi_oc'
          DEFVAL(4)=0.0  ! DEFAULT VALUE FOR GRID POINTS NOT COVERED BY
                         ! THE WAVE MODEL ICE FREE SEA POINTS.
        ENDIF

!       5. ENERGY FLUX TO WAVES
        IF(NWVFIELDS.GT.4) THEN
          FLABEL(5)=' Phi_aw'
          DEFVAL(5)=0.0  ! DEFAULT VALUE FOR GRID POINTS NOT COVERED BY
                         ! THE WAVE MODEL ICE FREE SEA POINTS.
        ENDIF

!       6. MOMENTUM FLUX TO OCEAN
        IF(NWVFIELDS.GT.5) THEN
          FLABEL(6)=' Tau_oc'
          DEFVAL(6)=0.0  ! DEFAULT VALUE FOR GRID POINTS NOT COVERED BY
                         ! THE WAVE MODEL ICE FREE SEA POINTS.
        ENDIF

!       7. WAVE VARIANCE
        IF(NWVFIELDS.GT.6) THEN
          FLABEL(7)=' Emean'
          DEFVAL(7)=0.0  ! DEFAULT VALUE FOR GRID POINTS NOT COVERED BY
                         ! THE WAVE MODEL ICE FREE SEA POINTS.
        ENDIF
  
!       8. MEAN FREQUENCY
        IF(NWVFIELDS.GT.7) THEN
          FLABEL(8)=' Fmean'
          DEFVAL(8)=0.0  ! DEFAULT VALUE FOR GRID POINTS NOT COVERED BY
                         ! THE WAVE MODEL ICE FREE SEA POINTS.
        ENDIF

        NPROMA=NPROMA_WAM
        CALL GSTATS(1443,0)
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO=IJS(IG),IJL(IG),NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
          CALL OUTBETA (KIJS, KIJL, U10OLD(KIJS,IG),
     &                  USOLD(KIJS,IG), Z0OLD(KIJS,IG),
     &                  WVBLOCK(KIJS,1))

          BETAOLD(KIJS:KIJL)=WVBLOCK(KIJS:KIJL,1)

!         SURFACE STOKES DRIFT NEEDED FOR THE IFS
!         IT MIGHT ALSO BE USED FOR NEMO !!!!!!!!!! 

          IF(LWSTOKES) THEN
             CALL STOKESDRIFT(FL1(KIJS:KIJL,1:NANG,1:NFRE),KIJS,KIJL,
     &                        WVBLOCK(KIJS,2),WVBLOCK(KIJS,3))
          ELSE
             IF(NWVFIELDS.GT.1) THEN
               WVBLOCK(KIJS:KIJL,2)=DEFVAL(2)
             ENDIF
             IF(NWVFIELDS.GT.2) THEN
               WVBLOCK(KIJS:KIJL,3)=DEFVAL(3)
             ENDIF
          ENDIF

          IF(LWFLUX) THEN
             WVBLOCK(KIJS:KIJL,4)=PHIEPS(KIJS:KIJL)
             WVBLOCK(KIJS:KIJL,5)=PHIAW(KIJS:KIJL)
             WVBLOCK(KIJS:KIJL,6)=TAUOC(KIJS:KIJL)
             WVBLOCK(KIJS:KIJL,7)=EMEAN(KIJS:KIJL)
             WVBLOCK(KIJS:KIJL,8)=FMEAN(KIJS:KIJL)
          ELSE
             IF(NWVFIELDS.GT.3) THEN
               WVBLOCK(KIJS:KIJL,4)=DEFVAL(4)
             ENDIF
             IF(NWVFIELDS.GT.4) THEN
               WVBLOCK(KIJS:KIJL,5)=DEFVAL(5)
             ENDIF
             IF(NWVFIELDS.GT.5) THEN
               WVBLOCK(KIJS:KIJL,6)=DEFVAL(6)
             ENDIF             
             IF(NWVFIELDS.GT.6) THEN
               WVBLOCK(KIJS:KIJL,7)=DEFVAL(7)
             ENDIF    
             IF(NWVFIELDS.GT.7) THEN
               WVBLOCK(KIJS:KIJL,8)=DEFVAL(8)
             ENDIF    
          ENDIF
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1443,1)


!        GRIDDED FIELDS THAT WILL RETURNED TO NEMO (if coupled).
!        -----------------------------------------------------

        IF(LWNEMOCOU.AND.LWNEMOCOUSEND.AND.LWCOU) THEN
          CALL GSTATS(1443,0)
#ifdef PARKIND1_SINGLE
!$OMP     PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
!$OMP+    PRIVATE(ZNEMOUSTOKES,ZNEMOVSTOKES,ZNEMOSTRN)
#else
!$OMP     PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
#endif
          DO JKGLO=IJS(IG),IJL(IG),NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
            IF(LWNEMOCOUSTK) THEN
              IF(LWSTOKES) THEN
!               IT WAS ALREADY COMPUTED FOR THE IFS
                NEMOUSTOKES(KIJS:KIJL)=WVBLOCK(KIJS:KIJL,2)
                NEMOVSTOKES(KIJS:KIJL)=WVBLOCK(KIJS:KIJL,3)
              ELSE
!               IT NEEDS TO BE COMPUTED
#ifdef PARKIND1_SINGLE
!     Single precision -- needs tmp copies
                CALL STOKESDRIFT(FL1(KIJS:KIJL,1:NANG,1:NFRE),KIJS,KIJL,
     &                           ZNEMOUSTOKES(1),ZNEMOVSTOKES(1))
                NEMOUSTOKES(KIJS:KIJL) = ZNEMOUSTOKES(1:KIJL-KIJS+1)
                NEMOVSTOKES(KIJS:KIJL) = ZNEMOVSTOKES(1:KIJL-KIJS+1)
#else
!     Double precision
                CALL STOKESDRIFT(FL1(KIJS:KIJL,1:NANG,1:NFRE),KIJS,KIJL,
     &                           NEMOUSTOKES(KIJS),NEMOVSTOKES(KIJS))
#endif
              ENDIF
            ELSE
              NEMOUSTOKES(KIJS:KIJL)=0.
              NEMOVSTOKES(KIJS:KIJL)=0.
            ENDIF

            IF(LWNEMOCOUSTRN) THEN
#ifdef PARKIND1_SINGLE
!     Single precision -- needs a tmp copy
              CALL CIMSSTRN(FL1(KIJS:KIJL,1:NANG,1:NFRE),KIJS,KIJL,
     &                      ZNEMOSTRN(1))
              NEMOSTRN(KIJS:KIJL) = ZNEMOSTRN(1:KIJL-KIJS+1)
#else
!     Double precision
              CALL CIMSSTRN(FL1(KIJS:KIJL,1:NANG,1:NFRE),KIJS,KIJL,
     &                      NEMOSTRN(KIJS))
#endif
            ENDIF

          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1443,1)
        ENDIF

!       GRIDDED FIELDS ARE NEEDED FOR IFS, GATHERING OF THE NECESSARY 
!       INFORMATION

        CALL MPFLDTOIFS(IG, IJS(IG), IJL(IG), NWVFIELDS, WVBLOCK,
     &                  WVFLDG, DEFVAL, MASK_OUT, LLGLOBAL_WVFLDG)


        DEALLOCATE(WVBLOCK)


!       COMPUTATION OF THE NORMS OF OUTPUT FIELDS
        WRITE(IU06,*) ' '
        IF(LLGLOBAL_WVFLDG)THEN
          WRITE(IU06,*) ' GLOBAL NORM OF FIELDS RETURNED TO IFS :'

          N_MASK_OUT_GLOBAL=NLATW*NLONW
          CALL GSTATS(1443,0)
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(JF,J,I)
          DO IFLD=1,NWVFIELDS
            JF=IFLD
            FAVG(JF)=0.
            FMIN(JF)=WVFLDG(1,1,JF)
            FMAX(JF)=WVFLDG(1,1,JF)
            DO J=1,NLATW
              DO I=1,NLONW
                FAVG(JF)=FAVG(JF)+WVFLDG(I,J,JF)
                FMIN(JF)=MIN(FMIN(JF),WVFLDG(I,J,JF))
                FMAX(JF)=MAX(FMAX(JF),WVFLDG(I,J,JF))
              ENDDO
            ENDDO
            FAVG(JF)=FAVG(JF)/N_MASK_OUT_GLOBAL
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1443,1)

          DO IFLD=1,NWVFIELDS
            WRITE(IU06,*) FLABEL(IFLD),
     &               FAVG(IFLD),FMIN(IFLD),FMAX(IFLD),N_MASK_OUT_GLOBAL
            WRITE(IU06,111) 
     &               FAVG(IFLD),FMIN(IFLD),FMAX(IFLD)
          ENDDO
111       FORMAT(14x,'HEX: ',3(Z16.16,2x))
          WRITE(IU06,*) ' '
          CALL FLUSH(IU06)

        ELSE


!         Local NORM OF FIELDS RETURNED TO IFS :
!         INITIALISE I,J POINTER TO MASK_OUT
          IF(LFRST) THEN
            LFRST=.FALSE.
            JF=1
            N_MASK_OUT=0
            DO J=1,NLATW
              DO I=1,NLONW
                IF(MASK_OUT(I,J).EQ.1) THEN
                  N_MASK_OUT=N_MASK_OUT+1
                ENDIF
              ENDDO
            ENDDO
            ALLOCATE(I_MASK_OUT(MAX(1,N_MASK_OUT)))
            I_MASK_OUT(1)=1
            ALLOCATE(J_MASK_OUT(MAX(1,N_MASK_OUT)))
            J_MASK_OUT(1)=1
            JF=1
            IC=0
            DO J=1,NLATW
              DO I=1,NLONW
                IF(MASK_OUT(I,J).EQ.1) THEN
                  IC=IC+1
                  I_MASK_OUT(IC)=I
                  J_MASK_OUT(IC)=J
                ENDIF
              ENDDO
            ENDDO
          ENDIF

!         COMPUTE local NORMS
          DO JF=1,NWVFIELDS
            VAL=WVFLDG(I_MASK_OUT(1),J_MASK_OUT(1),JF)
            FAVG(JF)=VAL
            FMIN(JF)=VAL
            FMAX(JF)=VAL
            DO IC=2,N_MASK_OUT
              VAL=WVFLDG(I_MASK_OUT(IC),J_MASK_OUT(IC),JF)
              FAVG(JF)=FAVG(JF)+VAL
              FMIN(JF)=MIN(FMIN(JF),VAL)
              FMAX(JF)=MAX(FMAX(JF),VAL)
            ENDDO
            FAVG(JF)=FAVG(JF)/MAX(N_MASK_OUT,1)
          ENDDO

!         FOR PRIMARY PE (IRECV), COLLECT ALL THE NORMS TO PRODUCE A
!         PSEUDO GLOBAL NORM.
 
          IRECV=1
          NCOMLOC=1+3*NWVFIELDS
          NCOMBUF=NCOMLOC*NPROC
          ALLOCATE(ZCOMBUFS(NCOMBUF))

          IST=1+(IRANK-1)*NCOMLOC
          IED=IST+NCOMLOC-1
          ICOUNT=IST
          ZCOMBUFS(ICOUNT)=N_MASK_OUT
          DO IFLD=1,NWVFIELDS
            ICOUNT=ICOUNT+1
            ZCOMBUFS(ICOUNT)=FAVG(IFLD)
            ICOUNT=ICOUNT+1
            ZCOMBUFS(ICOUNT)=FMIN(IFLD)
            ICOUNT=ICOUNT+1
            ZCOMBUFS(ICOUNT)=FMAX(IFLD)
          ENDDO

          IF(IRANK.EQ.IRECV) THEN
            ALLOCATE(ZCOMBUFR(NCOMBUF))
          ENDIF
          ALLOCATE(ZCOMCNT(NPROC))
          ZCOMCNT=NCOMLOC
          CALL MPL_GATHERV(PSENDBUF=ZCOMBUFS(IST:IED),KROOT=IRECV,
     &                    PRECVBUF=ZCOMBUFR(:),KRECVCOUNTS=ZCOMCNT,
     &                    CDSTRING='WAVEMDL:')
          DEALLOCATE(ZCOMCNT)

!         COMPUTE PSEUDO GLOBAL NORM
          IF(IRANK.EQ.IRECV) THEN
          WRITE(IU06,*) ' Local NORM OF global FIELDS RETURNED TO IFS :'
            IST=1+(IRANK-1)*NCOMLOC
            ICOUNT=IST
            NMASK=ZCOMBUFR(ICOUNT)
            NTOT=NMASK
            DO JF=1,NWVFIELDS
              ICOUNT=ICOUNT+1
              FAVG(JF)=NMASK*ZCOMBUFR(ICOUNT)
              ICOUNT=ICOUNT+2
            ENDDO
            DO IP=1,NPROC
              IF(IP.NE.IRECV) THEN
                IST=1+(IP-1)*NCOMLOC
                ICOUNT=IST
                NMASK=ZCOMBUFR(ICOUNT)
                NTOT=NTOT+NMASK
                DO JF=1,NWVFIELDS
                  ICOUNT=ICOUNT+1
                  FAVG(JF)=FAVG(JF)+NMASK*ZCOMBUFR(ICOUNT)
                  ICOUNT=ICOUNT+1
                  FMIN(JF)=MIN(FMIN(JF),ZCOMBUFR(ICOUNT))
                  ICOUNT=ICOUNT+1
                  FMAX(JF)=MAX(FMAX(JF),ZCOMBUFR(ICOUNT))
                ENDDO
              ENDIF
            ENDDO
            DEALLOCATE(ZCOMBUFR)
            DO JF=1,NWVFIELDS
              FAVG(JF)=FAVG(JF)/MAX(NTOT,1)
            ENDDO
          ELSE
            NTOT=N_MASK_OUT
            WRITE(IU06,*) ' Local NORM OF FIELDS RETURNED TO IFS :'
          ENDIF

          DO IFLD=1,NWVFIELDS
            WRITE(IU06,*) FLABEL(IFLD),
     &                    FAVG(IFLD),FMIN(IFLD),FMAX(IFLD),NTOT,
     &                    IRANK, NPROC, LL1D
            WRITE(IU06,111) 
     &                    FAVG(IFLD),FMIN(IFLD),FMAX(IFLD)
          ENDDO
          IF (ITEST.GE.1)  CALL FLUSH(IU06)

          DEALLOCATE(ZCOMBUFS)
        ENDIF


!       TELL ATMOS MODEL THE WAM REQUIREMENT FOR GLOBAL NORMS
        LDWCOUNORMS=LWCOUNORMS
      ENDIF

#ifdef MODEL_COUPLING_ATM_WAV
      Z0_coupl=Z0OLD
      US_coupl=USOLD
#endif


!     4. END OF RUN ?
!        -----------
      IF (CDATEE .EQ. CDTPRO) THEN
        CALL MPL_BARRIER(CDSTRING='WAVEMDL: END')
        CALL FLUSH(IU06)
      ENDIF
      IF ( NWFDBREF .GT. 0 .AND. CDATEE .EQ. CDTPRO) THEN
        IF(LWCOU .AND. LNEWLVTP) THEN
!         the fdb should be closed by IFS
          IF (ITEST.GE.1) THEN
            WRITE(IU06,*) ' SUB. WAVEMDL: END OF WAVE MODEL RUN'
            CALL FLUSH(IU06)
          ENDIF
        ELSE
        WRITE(IU06,*) ' SUB. WAVEMDL: END OF RUN: CLOSING FDB'
        CALL GSTATS(1787,0)
        ISTAT = ICLOSEFDBSUBS (NWFDBREF)
        CALL GSTATS(1787,1)
          IF (ITEST.GE.1) THEN
            WRITE(IU06,*) ' SUB. WAVEMDL: END OF RUN: CLOSE FDB'
            CALL FLUSH(IU06)
          ENDIF
        ENDIF
      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('WAVEMDL',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE WAVEMDL
