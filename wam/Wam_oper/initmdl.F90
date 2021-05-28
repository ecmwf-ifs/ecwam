      SUBROUTINE INITMDL (NADV,                                         &
     &                    IREAD,                                        &
     &                    NFIELDS, NGPTOTG, NC, NR,                     &
     &                    FIELDS, LWCUR, MASK_IN, PRPLRADI)

! ----------------------------------------------------------------------

!**** *INITMDL* - INITIALIZES THE WAM MODEL.

!     L. ZAMBRESKY   GKSS/ECMWF    JULY 1988

!     MODIFIED BY:   H. GUNTHER    NOVEMBER 1989
!                    J. BIDLOT     FEBRUARY 1996-1997
!                    J. DOYLE      SEPTEMBER 1996
!                    J. BIDLOT     APRIL 1997
!                    J. BIDLOT     SEPTEMBER : MODIFY PARALLEL INPUT 
!                    S. ABDALLA    OCTOBER 2001: INCLUSION OF AIR 
!                                                DENSITY AND Zi/L
!                    J. BIDLOT     AUGUST 2008. ADD CALL TO PREWIND
!                                  TO GET CURRENTS EARLIER
!                                  (i.e BEFORE CALL TO GETSPEC).

!SHALLOW
!          DIFFERENCES FOR SHALLOW WATER RUNS TO DEEP WATER RUNS
!          ARE ENCLOSED IN COMMENT LINES : 'CSHALLOW'.
!SHALLOW
!NEST
!          DIFFERENCES FOR NESTED GRID RUNS TO NORMAL RUNS
!          ARE ENCLOSED IN COMMENT LINES : 'CNEST'.
!NEST
!REFRA
!          DIFFERENCES FOR REFRACTION RUNS TO NORMAL RUNS
!          ARE ENCLOSED IN COMMENT LINES : 'CREFRA'.
!REFRA

!*    PURPOSE.
!     --------

!       INITIALIZE THE WAM MODEL.

!**   INTERFACE.
!     ----------

!          ---- FORMAL PARAMETERS ----

!    *CALL* *INITMDL (NADV,
!    &                IREAD,
!    &                NFIELDS, NGPTOTG, NC, NR,
!    &                FIELDS, LWCUR, MASK_IN)*

!      *NADV*      INTEGER   NUMBER OF ADVECTION ITERATIONS
!      *IREAD*     INTEGER   PROCESSOR WHICH WILL ACCESS THE FILE ON DISK
!                            (IF NEEDED).
!      *NFIELDS*   INTEGER   NUMBER OF FIELDS HOLDING ATMOSPHERIC DATA
!      *NGPTOTG*   INTEGER   NUMBER OF ATMOSPHERIC GRID POINTS
!      *NC*        INTEGER   NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR
!      *NR*        INTEGER   NUMBER OF ATM. ROWS OF LATITUDES
!      *FIELDS*    REAL      ATM. FIELDS (U10, V10, AIR DENSITY, Zi/L, U and V CURRENTS)
!      *LWCUR*    LOGICAL   INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS
!      *MASK_IN*   INTEGER   MASK TO INDICATE WHICH PART OF FIELDS IS RELEVANT.
!      *PRPLRADI*  REAL      EARTH RADIUS REDUCTION FACTOR FOR SMALL PLANET

!          ---- INPUT/OUTPUT UNITS ---

!          THE NAMES ARE DEFINED IN SECTION 1. OF THIS PROGRAM,
!          IF IT IS NOT MENTIONED OTHERWISE.

!           *IU01*   - INPUT  UNIT UNBLOCKED WIND FILE.
!                      (SEE SUB READWIND).
!           *IU02*   - INPUT  UNIT OF BOUDARY VALUES FROM A PREVIOUS
!                      COARSE GRID IF THIS A FINE GRID RUN.
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'FBI'
!                      (OUTPUT OF BOUINT).
!           *IU06*   - PRINTER OUTPUT.
!           *IU07*   - INPUT  UNIT OF PRECOMPUTED GRID PARAMETERS.
!                      (OUTPUT OF PREPROC).
!           *IU08*   - INPUT  UNIT OF MODULE YOWUBUF.
!                      (OUTPUT OF PREPROC).
!NEST
!           *IU09*   - INPUT  UNIT MODULE YOWCPBO (OUTPUT OF PREPROC).
!           *IU10*   - INPUT  UNIT MODULE YOWFPBO (OUTPUT OF PREPROC).
!NEST
!           *IU11*   - INPUT UNIT OF SPECTRA AT ALL GRID POINTS.
!                      EACH PROPAGATION STEP THE FILES CONNECTED
!                      TO IU11 AND IU12 ARE INTERCHANGED.
!           *IU12*   - OUTPUT UNIT BLOCKS OF SPECTRA (SEE IU11).
!           *IU13*   - INPUT UNIT OF SPECTRA ON LAST LATITUDE
!                      OF A BLOCK. SPECTRA ARE SAVED FROM THE
!                      SECOND LATITUDE OF NEXT BLOCK.
!                      EACH PROPAGATION STEP THE FILES CONNECTED
!                      TO IU13 AND IU14 ARE INTERCHANGED.
!           *IU14*   - OUTPUT UNIT SECOND LATITUDES (SEE IU13).
!           *IU15*   - OUTPUT UNIT LAST WINDFIELDS.
!REFRA
!           *IU16*   - INPUT/OUTPUT UNIT OF MODULE YOWREFD.
!REFRA
!           *IU17*   - INPUT  UNIT OF BLOCKED WINDS.
!                      THIS FILE IS DYNAMICALLY ASSIGNED IN SUB
!                      IMPLSCH.
!           *IU18*   - INPUT  UNIT OF BLOCKED WINDS.
!                      THIS FILE IS DYNAMICALLY ASSIGNED IN SUB
!                      IMPLSCH.
!NEST
!           *IU19*   - OUTPUT UNIT OF BOUNDARY VALUES IF
!                      THIS IS A FINE GRID RUN.
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'CBO'
!NEST
!           *IU20*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
!                      OF THE TOTAL SPECTRUM (FIRST GUESS).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'MAP'
!ASSI
!           *IU22*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
!                      OF THE TOTAL SPECTRUM (ANALYSED).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'AMP'
!           *IU23*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
!                      OF SWELL AND WIND WAVES (ANALYSED)
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'ASS'
!           *IU27*   - OUTPUT UNIT OF SPECTRA AT CERTAIN GRID POINTS
!                      (ANALYSED).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'AUT'
!           *IU28*   - OUTPUT UNIT OF SWELL SPECTRA AT CERTAIN POINTS
!                      (ANALYSED).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'ASW'
!ASSI
!           *IUVELO* - OUTPUT UNIT OF BLOCKED WINDFILEDS.
!                      FILES ARE DYNAMICALLY ASSIGNED IN SUB
!                      NOTIM OR TIMIN.

!___PACK  THESE UNITS ARE USED FOR FILES CONTAINING GRIB DATA.

!         (IN CASE OF SPECTRA ALL SPECTRA OF THE MODEL ARE WRITTEN)
!           *IU30*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
!                      OF THE TOTAL SPECTRUM (FIRST GUESS).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'MPP'
!           *IU31*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
!                      OF SWELL AND WIND WAVES (FIRST GUESS).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'SWP'
!ASSI
!           *IU32*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
!                      OF THE TOTAL SPECTRUM (ANALYSED).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'APP'
!           *IU33*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
!                      OF SWELL AND WIND WAVES (ANALYSED)
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'AWP'
!ASSI
!           *IU35*   - OUTPUT UNIT OF SPECTRA AT ALL GRID POINTS
!                      (FIRST GUESS).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'OUP'
!           *IU36*   - OUTPUT UNIT OF SWELL SPECTRA AT ALL POINTS
!                      (FIRST GUESS).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'WSP'
!ASSI
!           *IU37*   - OUTPUT UNIT OF SPECTRA AT ALL GRID POINTS
!                      (ANALYSED).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'AUP'
!           *IU38*   - OUTPUT UNIT OF SWELL SPECTRA AT ALL POINTS
!                      (ANALYSED).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'ASP'
!ASSI
!___PACK

!          FOR A START THE RESTART FILES WILL DYNAMICALLY BE ASSIGNED
!          AND COPIED TO THE *OUTPUT* UNITS (IU12, IU14, IU15).
!          IF IT IS REQUESTED TO SAVE RESTART FILES THESE WILL BE COPIED
!          IN REGULAR INTERVALS FROM THE UNIT ALIAS FILES OF IU11 OR
!          IU12, IU13 OR IU14, AND IU15 TO THE PREMANENT RESTART FILES.
!          FOR DETAILS OF THE FILE NAMES SEE SUB GSFILE.

!          SUB NOTIM OR TIMIN OPENS FILES AND ASSIGNS THEM TO UNIT
!          IUVELO FOR THE BLOCKED WINDS. THESE FILES ARE READ AND
!          DELETED IN SUB IMPLSCH (IU17 AND IU18).

!          THE PROGRAM CLOSES AND DELETES ALL WORK FILES.

!          ALL PARAMETERS HAVE TO BE THE VALUES GIVEN AT THE END
!          OF THE PREPROC OUTPUT IN COLUMN 'REQUIRED'.

!     METHOD.
!     -------

!          THIS ROUTINE INITIALISES THE WAVEMODEL:
!            -  DEFINES THE UNITS FOR INPUT/OUTPUT,
!            -  READS THE USER INPUT FILE,
!            -  INITIALIZES SOME MODEL CONSTANTS,
!            -  GETS THE RECOVERY FILES,
!            -  READS THE YOWMON BLOCKS PRECOMPUTED BY PROG PREPROC,
!            -  DOES SOME GENERAL BOOKEEPING REGARDING
!               DATES, INTEGRATION TIME STEPS AND OUTPUT TIME STEPS.
!            -  READ MODULE YOWUBUF AND SPECTRA IF ONE BLOCK VERSION.
!REFRA
!            -  PRECOMPUTES AND WRITES TO IU16 REFRACTION TERMS.
!REFRA
!            -  OPENS THE FIRST RESULT FILES.

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.
!       *DIFDATE*   - COMPUTES A TIME DIFFERENCE.
!REFRA
!       *GRADI*     - COMPUTES DEPTH AND CURRENT GRADIENTS.
!REFRA
!       *GSFILE*    - ROUTINE TO DYNAMICALLY FETCH OR DISPOSE FILES.
!       *GETSPEC*   - READS SPECTRA
!       *GETSTRESS  - READS RESTART STRESS/WIND FILE
!       *HEADBC*    - WRITE BOUNDARY OUTPUT FILE HEADER.
!       *INCDATE*   - UPDATE DATE TIME GROUP.
!       *INIWCST*   - SETS GLOBAL CONSTANTS.
!       *PREWIND*   - REFORMAT FORCING FIELDS.
!REFRA
!       *PROPDOT*   - PRECOMPUTE REFRACTION.
!REFRA
!NEST
!       *READBOU*   - READS PREPROC BOUNDARY FILES.
!NEST
!       *SETMARSTYPE* SETS VARIABLE MARSTYPE
!       *SETWGRIBHD*- SETS DEFAULT GRIB HEADERS  
!       *USERIN*    - READS USER INPUT.

!     REFERENCE
!     ---------

!          A MORE DETAILED  DISCUSSION MAY BE FOUND IN SUB WAMODEL.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCPBO  , ONLY : IBOUNC   ,NBOUNC   ,IJARC    ,IGARC,        &
     &            GBOUNC  , IPOGBO   ,CBCPREF
      USE YOWCOUP  , ONLY : LWCOU    ,KCOUSTEP ,LWFLUX   ,              &
     &                      LWNEMOCOU,LWNEMOCOURECV
      USE YOWCOUT  , ONLY : COUTT    ,COUTLST  ,                        &
     &            FFLAG20  ,GFLAG20  ,                                  &
     &            NGOUT    ,IGAR     ,IJAR     ,NOUTT    ,LOUTINT  ,    &
     &            LWFLUXOUT
      USE YOWCURR  , ONLY : CDTCUR   ,IDELCUR  ,CDATECURA,U        ,    &
     &            V
      USE YOWFPBO  , ONLY : IBOUNF
      USE YOWGRIB_HANDLES , ONLY :NGRIB_HANDLE_WAM_I,NGRIB_HANDLE_WAM_S
      USE YOWFRED  , ONLY : FR       ,TH       ,DELTH   ,FR5      ,     &
     &            ZPIFR    ,                                            &
     &            FRM5     ,COFRM4   ,COEF4    ,FRATIO  ,FLOGSPRDM1,    &
     &            COSTH    ,SINTH    ,FLMAX    ,RHOWG_DFIM,             &
     &            DFIM     ,DFIMOFR  ,DFIMFR  ,DFIMFR2   ,              &
     &            DFIM_SIM ,DFIMOFR_SIM ,DFIMFR_SIM ,DFIMFR2_SIM ,      &
     &            DFIM_END_L, DFIM_END_U
      USE YOWGRIBHD, ONLY : LGRHDIFS
      USE YOWGRID  , ONLY : DELPHI   ,DELLAM   ,IJS     ,IJL      ,COSPH
      USE YOWICE   , ONLY : CICOVER  ,CITHICK  ,CIWA
      USE YOWINDN  , ONLY : MLSTHG   ,ENH
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,AMOWEP   ,AMOSOP   ,    &
     &            AMOEAP   ,AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,    &
     &            KMNOP    ,KMSOP    ,IPER
      USE YOWMEAN  , ONLY : PHIEPS   ,PHIAW    ,TAUOC    ,              &
     &                      TAUXD    ,TAUYD    ,WSEMEAN  ,WSFMEAN  ,    &
     &                      TAUOCXD  ,TAUOCYD  ,PHIOCD
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP     ,    &
     &            KTAG
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_RED ,NFRE_ODD ,    & 
     &            NGX      ,NGY      ,                                  &
     &            NIBLO    ,NIBLD    ,NBLD     ,NIBLC    ,NBLC
      USE YOWPCONS , ONLY : G        ,CIRC     ,PI       ,ZPI      ,    &
     &            RAD      ,ROWATER  ,ZPI4GM2  ,FM2FP
      USE YOWPHYS  , ONLY : ALPHAPMAX, ALPHAPMINFAC, FLMINFAC
      USE YOWREFD  , ONLY : THDD     ,THDC     ,SDOT
      USE YOWSHAL  , ONLY : NDEPTH   ,DEPTH    ,DEPTHA   ,DEPTHD   ,    &
     &            INDEP    ,TCGOND   ,IODP     ,IOBND    ,TOOSHALLOW,   &
     &            CINV     ,TFAK     ,GAM_B_J  ,EMAXDPT
      USE YOWSPEC  , ONLY : NBLKS    ,NBLKE    ,KLENTOP  ,KLENBOT  ,    &
     &            U10OLD   ,THWOLD   ,USOLD    ,Z0OLD    ,TAUW     ,    &
     &            Z0B      ,TAUWDIR,  ROAIRO   ,ZIDLOLD  ,              &
     &            FL1
      USE YOWSTAT  , ONLY : CDATEE   ,CDATEF   ,CDTPRO   ,CDTRES   ,    &
     &            CDTINTT  ,CDTBC    ,                                  &
     &            IDELPRO  ,IDELT    ,IDELWI   ,IDELWO   ,IDELRES  ,    &
     &            IDELINT  ,                                            &
     &            ISHALLO  ,IREFRA   ,                                  &
     &            IPHYS    ,                                            &
     &            NPROMA_WAM,                                           &
     &            CDATEA   ,MARSTYPE ,LANAONLY ,ISNONLIN ,IPROPAGS ,    &
     &            IDELWI_LST,IDELWO_LST,CDTW_LST,NDELW_LST
      USE YOWTABL  , ONLY : FAC0     ,FAC1     ,FAC2     ,FAC3     ,    &
     &            FAK      ,FRHF     ,DFIMHF
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWTEXT  , ONLY : LRESTARTED
      USE YOWWNDG  , ONLY : ICODE
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,KCOR     ,              &
     &            KRLAT    ,KRLON    ,LUPDTWGHT,                        &
     &            JXO      ,JYO      ,KCR      ,KPM      ,MPM      ,    &
     &            LSAMEDEPTH
      USE YOWUNIT  , ONLY : IU02     ,IU11     ,IU12     ,              &
     &            IU13     ,IU14     ,IU15     ,IU16     ,IU17     ,    &
     &            IU18     ,IU19     ,IU20     ,IU22     ,              &
     &            IU23     ,IU27     ,IU28     ,                        &
     &            IUVELO   ,IU30     ,IU31     ,                        &
     &            IU32     ,IU33     ,IU35     ,IU36     ,IU37     ,    &
     &            IU38
      USE YOWWAMI  , ONLY : CBPLTDT
      USE YOWWIND  , ONLY : CDA      ,CDAWIFL  ,CDATEWO  ,CDATEFL  ,    &
     &            LLNEWCURR,LLWSWAVE ,LLWDWAVE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWUNPOOL, ONLY : LLUNSTR, OUT_METHOD
      USE UNSTRUCT_BOUND , ONLY : IOBP
      USE OUTPUT_STRUCT, ONLY : INITIAL_OUTPUT_INITS

! -------------------------------------------------------------------

      IMPLICIT NONE
#include "wam_nproma.intfb.h"
#include "abort1.intfb.h"
#include "cigetdeac.intfb.h"
#include "cireduce.intfb.h"
#include "getspec.intfb.h"
#include "getstress.intfb.h"
#include "headbc.intfb.h"
#include "incdate.intfb.h"
#include "inisnonlin.intfb.h"
#include "init_sdiss_ardh.intfb.h"
#include "init_x0tauhf.intfb.h"
#include "initnemocpl.intfb.h"
#include "iniwcst.intfb.h"
#include "preset_wgrib_template.intfb.h"
#include "prewind.intfb.h"
#include "propdot.intfb.h"
#include "readbou.intfb.h"
#include "setmarstype.intfb.h"
#include "tabu_swellft.intfb.h"
#include "userin.intfb.h"

      INTEGER(KIND=JWIM), INTENT(OUT) :: NADV
      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      INTEGER(KIND=JWIM), INTENT(IN) :: NFIELDS
      INTEGER(KIND=JWIM), INTENT(IN) :: NGPTOTG
      INTEGER(KIND=JWIM), INTENT(IN) :: NC
      INTEGER(KIND=JWIM), INTENT(IN) :: NR
      REAL(KIND=JWRB), INTENT(IN) :: PRPLRADI
      LOGICAL, INTENT(IN) :: LWCUR
      INTEGER(KIND=JWIM),DIMENSION(NGPTOTG), INTENT(INOUT) :: MASK_IN
      REAL(KIND=JWRB),DIMENSION(NGPTOTG,NFIELDS), INTENT(IN) :: FIELDS


      INTEGER(KIND=JWIM) :: IFORCA
      INTEGER(KIND=JWIM) :: IJ, I, II, K, M, IP, LFILE, IX, IY, KX, ID
      INTEGER(KIND=JWIM) :: IC, ICR
      INTEGER(KIND=JWIM) :: KM1, KP1
      INTEGER(KIND=JWIM) :: JD
      INTEGER(KIND=JWIM) :: IDELWH
      INTEGER(KIND=JWIM) :: IU05, IU09, IU10
      INTEGER(KIND=JWIM) :: IWAM_GET_UNIT
      INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA
      INTEGER(KIND=JWIM) :: NTOT, MTHREADS
!$    INTEGER,EXTERNAL :: OMP_GET_MAX_THREADS


      REAL(KIND=JWRB) :: FCRANGE, XD
      REAL(KIND=JWRB) :: GVE, DPH, CFLP, CFLL, DLH, DLH_KX
      REAL(KIND=JWRB) :: FAC, SCDF_L, SCDF_U
      REAL(KIND=JWRB) :: D, OM, XK
      REAL(KIND=JWRB) :: XLOGFRATIO
      REAL(KIND=JWRB) :: GAM
      REAL(KIND=JWRB), PARAMETER :: ENH_MAX=10.0_JWRB
      REAL(KIND=JWRB), PARAMETER :: ENH_MIN=0.1_JWRB   ! to prevent ENH to become too small
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: AKI, TRANSF

      REAL(KIND=JWRB) :: XLA, XLO 

      CHARACTER(LEN=14) :: ZERO, CDUM
      CHARACTER(LEN=24) :: FILNM
      CHARACTER(LEN=80) :: FILENAME

      LOGICAL :: LLEXIST
      LOGICAL :: LLINIT
      LOGICAL :: LLALLOC_FIELDG_ONLY

!----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('INITMDL',0,ZHOOK_HANDLE)

      CALL INIWCST(PRPLRADI)

!     ADJUST NPROMA_WAM
      NTOT=IJL-IJS+1
      MTHREADS=1
!$    MTHREADS=OMP_GET_MAX_THREADS()
      NPROMA=NPROMA_WAM
      CALL WAM_NPROMA(NTOT, MTHREADS, NPROMA)
      NPROMA_WAM=NPROMA

!*    1. DEFINITION OF MODEL PARAMETERS.
!        -------------------------------

      ZERO = ' '

!*    1.1  DEFINE UNIT NAMES.
!          ------------------

!NEST
      IU02 = 0
!NEST


! GRID AND UBUF FILES UNITS
! -------------------------

!    DEFINITION HAS  BEEN MOVED TO MPDECOMP

      IU11 = 11


! RESTART FILES UNITS
! -------------------

! NOTE IF MESSAGE PASSING
! NOTE IU12 WILL ONLY BE CONNECTED TO UNIT 12 FOR PURE BINARY INPUT
! NOTE OTHERWISE IT WILL CONNECTED TO THE FILENAME BLS WITH
! NOTE AN EXTENSION FUNCTION OF THE PE.

      IF (LWCOU) THEN
        IU12 = 36
        IU15 = 37
      ELSE
        IU12 = 12
        IU15 = 15
      ENDIF


      IU13 = 13
      IU14 = 14
!REFRA
      IU16 = 16
!REFRA
      IU17 = 17
      IU18 = 18

      IU20 = 20
!ASSI
      IU22 = 22
      IU23 = 23
      IU27 = 27
      IU28 = 28
!ASSI

!__PACK
      IU30 = 30
      IU31 = 31
!ASSI
      IU32 = 32
      IU33 = 33
!ASSI
      IU35 = 35
      IU36 = 36
!ASSI
      IU37 = 37
      IU38 = 38
!ASSI
!__PACK

      IUVELO = 39

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. INITMDL: UNITS DEFINED'
        CALL FLUSH (IU06)
      ENDIF

      LOUTINT=.FALSE.

      LUPDTWGHT=.TRUE.

      LLNEWCURR=.TRUE.

      ICODE=0

!*    1.2 INPUT OF USER PARAMETER.
!         ------------------------

      CALL USERIN (IFORCA, LWCUR)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '    SUB. INITMDL: USERIN DONE'
        CALL FLUSH(IU06)
      ENDIF

      IF (LWFLUXOUT) THEN
        IF (.NOT.ALLOCATED(PHIEPS)) THEN 
          ALLOCATE(PHIEPS(IJS:IJL))
          PHIEPS(:) = 0.0_JWRB
        ENDIF
        IF (.NOT.ALLOCATED(PHIAW)) THEN 
          ALLOCATE(PHIAW(IJS:IJL))
          PHIAW(:) = 0.0_JWRB
        ENDIF
        IF (.NOT.ALLOCATED(TAUOC)) THEN
          ALLOCATE(TAUOC(IJS:IJL))
          TAUOC(:) = 0.0_JWRB
        ENDIF

        IF(.NOT.ALLOCATED(TAUXD)) ALLOCATE(TAUXD(IJS:IJL))
        IF(.NOT.ALLOCATED(TAUYD)) ALLOCATE(TAUYD(IJS:IJL))
        IF(.NOT.ALLOCATED(WSEMEAN)) ALLOCATE(WSEMEAN(IJS:IJL))
        IF(.NOT.ALLOCATED(WSFMEAN)) ALLOCATE(WSFMEAN(IJS:IJL))
        IF(.NOT.ALLOCATED(TAUOCXD)) ALLOCATE(TAUOCXD(IJS:IJL))
        IF(.NOT.ALLOCATED(TAUOCYD)) ALLOCATE(TAUOCYD(IJS:IJL))
        IF(.NOT.ALLOCATED(PHIOCD)) ALLOCATE(PHIOCD(IJS:IJL))

      ENDIF

!     DEFINE COEFFICIENT FOR MEAN PERIODS CALCULATION
      IF (ALLOCATED(DFIMOFR)) DEALLOCATE(DFIMOFR)
      ALLOCATE(DFIMOFR(NFRE))
      IF(ALLOCATED(DFIMFR)) DEALLOCATE(DFIMFR)
      ALLOCATE(DFIMFR(NFRE))
      IF(ALLOCATED(DFIMFR2)) DEALLOCATE(DFIMFR2)
      ALLOCATE(DFIMFR2(NFRE))

      DO M=1,NFRE
        DFIMOFR(M) = DFIM(M)/FR(M)
        DFIMFR(M)  = DFIM(M)*FR(M)
        DFIMFR2(M) = DFIM(M)*FR(M)**2
      ENDDO

!     DEFINE A FEW CONSTANTS FOR USE IN IMPLSCH

      IF (.NOT.ALLOCATED(ZPIFR)) ALLOCATE(ZPIFR(NFRE))
      IF (.NOT.ALLOCATED(FR5)) ALLOCATE(FR5(NFRE))
      IF (.NOT.ALLOCATED(FRM5)) ALLOCATE(FRM5(NFRE))
      IF (.NOT.ALLOCATED(COFRM4)) ALLOCATE(COFRM4(NFRE))
      IF (.NOT.ALLOCATED(FLMAX)) ALLOCATE(FLMAX(NFRE))
      IF (.NOT.ALLOCATED(DFIM_END_L)) ALLOCATE(DFIM_END_L(NFRE))
      IF (.NOT.ALLOCATED(DFIM_END_U)) ALLOCATE(DFIM_END_U(NFRE))

      SCDF_L = 0.5_JWRB*DELTH*(FRATIO-1.0_JWRB)
      SCDF_U = 0.5_JWRB*DELTH*(1.0_JWRB-1.0_JWRB/FRATIO)

      DO M=1,NFRE
        ZPIFR(M) = ZPI*FR(M)
        FR5(M) = FR(M)**5
        FRM5(M) = 1.0_JWRB/FR5(M)
        COFRM4(M) = COEF4*G/FR(M)**4
        FLMAX(M) = (ALPHAPMAX/PI)/(ZPI4GM2*FR5(M))
        DFIM_END_L(M) = SCDF_L*FR(M)
        DFIM_END_U(M) = SCDF_U*FR(M)
      ENDDO

      FLMINFAC = ALPHAPMINFAC*FM2FP*G/(PI*ZPI**3*FR(NFRE)**5)

      FLOGSPRDM1=1.0_JWRB/LOG10(FRATIO)

      XLOGFRATIO=LOG(FRATIO)

      IF (.NOT.ALLOCATED(RHOWG_DFIM)) ALLOCATE(RHOWG_DFIM(NFRE))
      RHOWG_DFIM(1) = 0.5_JWRB*ROWATER*G*DELTH*XLOGFRATIO*FR(1)
      DO M=2,NFRE-1
        RHOWG_DFIM(M) = ROWATER*G*DELTH*XLOGFRATIO*FR(M)
      ENDDO
      RHOWG_DFIM(NFRE) = 0.5_JWRB*ROWATER*G*DELTH*XLOGFRATIO*FR(NFRE)

      IF (.NOT.ALLOCATED(DFIM_SIM)) ALLOCATE(DFIM_SIM(NFRE))
      NFRE_ODD=NFRE-1+MOD(NFRE,2)
      DFIM_SIM(NFRE)=0.0_JWRB
      DFIM_SIM(1)=DELTH*XLOGFRATIO*FR(1)/3.0_JWRB
      DO M=2,NFRE_ODD-1,2
        DFIM_SIM(M)=4.0_JWRB*DELTH*XLOGFRATIO*FR(M)/3.0_JWRB
        DFIM_SIM(M+1)=2.0_JWRB*DELTH*XLOGFRATIO*FR(M+1)/3.0_JWRB
      ENDDO
      DFIM_SIM(NFRE_ODD)=DELTH*XLOGFRATIO*FR(NFRE_ODD)/3.0_JWRB

      IF(.NOT.ALLOCATED(DFIMOFR_SIM)) ALLOCATE(DFIMOFR_SIM(NFRE))
      IF(.NOT.ALLOCATED(DFIMFR_SIM)) ALLOCATE(DFIMFR_SIM(NFRE))
      IF(.NOT.ALLOCATED(DFIMFR2_SIM)) ALLOCATE(DFIMFR2_SIM(NFRE))
      DO M=1,NFRE
        DFIMOFR_SIM(M) = DFIM_SIM(M)/FR(M)
        DFIMFR_SIM(M)  = DFIM_SIM(M)*FR(M)
        DFIMFR2_SIM(M) = DFIM_SIM(M)*FR(M)**2
      ENDDO


      CALL TABU_SWELLFT

      CALL INIT_X0TAUHF

      IF (IPHYS.EQ.1) CALL INIT_SDISS_ARDH


      IF (.NOT. LLUNSTR) THEN

      IF (LMESSPASS) THEN

        KTAG=100

!      1.3 OUTPUT MODEL DECOMPOSITION DETAILS 
!          ----------------------------------

        IF (ITEST.GE.1) THEN
          WRITE(IU06,*)  ' MODEL DOMAIN DECOMPOSITION  : '
          WRITE(IU06,*)  ' =========================='
          DO IP=1,NPROC
            WRITE(IU06,*)
            WRITE(IU06,*) ' PROCESS NUMBER : ',IP
            WRITE(IU06,*) ' NBLKS  : ',NBLKS(IP)
            WRITE(IU06,*) ' NBLKE  : ',NBLKE(IP)
            WRITE(IU06,*) ' N      : ',NBLKE(IP)-NBLKS(IP)+1
            WRITE(IU06,*) ' KLENBOT: ',KLENBOT(IP)
            WRITE(IU06,*) ' KLENTOP: ',KLENTOP(IP)
            WRITE(IU06,*) ' ----------------------- '
          ENDDO
          CALL FLUSH(IU06)
        ENDIF
      ENDIF
      ENDIF ! .NOT. LLUNSTR
!     1.5 DETERMINE LAST OUTPUT DATE
!         --------------------------

      IF (NOUTT.GT.0) THEN
        COUTLST=COUTT(NOUTT)
      ELSE
        COUTLST=CDATEA
        IF (GFLAG20) THEN
          DO WHILE (COUTLST.LE.CDATEE.AND.IDELINT.GT.0)
            CALL INCDATE (COUTLST,IDELINT)
          ENDDO
          CALL INCDATE (COUTLST,-IDELINT)
        ENDIF

        IF (COUTLST.GT.CDATEE) COUTLST=CDATEE
      ENDIF

! ----------------------------------------------------------------------

!*    2. INPUT FROM PREPROCESSING PROGRAMS.
!        ----------------------------------

!     THE ACTUAL READING HAS BEEN MOVED TO MPDECOMP.

!*    2.1 READ MODULE YOWCPBO AND YOWFPBO.
!     ------------------------------------

      GBOUNC = 1

      IF (.NOT. LLUNSTR) THEN
 
      IF (IBOUNC.EQ.1 .OR. IBOUNF.EQ.1) THEN
        IF (IBOUNC.EQ.1) THEN
!         READ INFORMATION ABOUT WHERE THE FINE GRID(S) ARE
          FILENAME='wam_nested_grids_info'
          LFILE=0
          LLEXIST=.FALSE.
          IF (FILENAME.NE. ' ') LFILE=LEN_TRIM(FILENAME)
          INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)
          IF (.NOT. LLEXIST) THEN
            WRITE(IU06,*) '************************************'
            WRITE(IU06,*) '*                                  *'
            WRITE(IU06,*) '*  FATAL ERROR IN SUB. INITMDL     *'
            WRITE(IU06,*) '*  =============================   *'
            WRITE(IU06,*) '*  WITH OPTION IBOUNC = 1          *' 
            WRITE(IU06,*) '*  YOU MUST PROVIDE THE FILE:      *' 
            WRITE(IU06,*) '* ',FILENAME(1:LFILE)
            WRITE(*,*) '*  WAVE MODEL INPUT FILE ',FILENAME(1:LFILE),   &
     &        ' IS MISSING !!!!'
            WRITE(IU06,*) '*                                  *'
            WRITE(IU06,*) '************************************'
            CALL ABORT1
          ENDIF
          IU09 = IWAM_GET_UNIT(IU06, FILENAME(1:LFILE) , 'r', 'u', 0)
        ELSE
          IU09 = 9
        ENDIF

        IF (IBOUNF.EQ.1) THEN
!         READ INFORMATION ABOUT WHERE THE BOUNDARY VALUES ARE 
          FILENAME='wam_boundary_grid_info'
          LFILE=0
          LLEXIST=.FALSE.
          IF (FILENAME.NE. ' ') LFILE=LEN_TRIM(FILENAME)
          INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)
          IF (.NOT. LLEXIST) THEN
            WRITE(IU06,*) '************************************'
            WRITE(IU06,*) '*                                  *'
            WRITE(IU06,*) '*  FATAL ERROR IN SUB. INITMDL     *'
            WRITE(IU06,*) '*  =============================   *'
            WRITE(IU06,*) '*  WITH OPTION IBOUNF = 1          *' 
            WRITE(IU06,*) '*  YOU MUST PROVIDE THE FILE:      *' 
            WRITE(IU06,*) '* ',FILENAME(1:LFILE)
            WRITE(*,*) '*  WAVE MODEL INPUT FILE ',FILENAME(1:LFILE),   &
     &        ' IS MISSING !!!!'
            WRITE(IU06,*) '*                                  *'
            WRITE(IU06,*) '************************************'
            CALL ABORT1
          ENDIF
          IU10 = IWAM_GET_UNIT(IU06, FILENAME(1:LFILE) , 'r', 'u', 0)
        ELSE
          IU10 = 10 
        ENDIF

        CALL READBOU (IU09, IU10, IU06)

        CLOSE (UNIT=IU09, STATUS='KEEP')
        CLOSE (UNIT=IU10, STATUS='KEEP')
        IF (ITEST.GE.2) THEN
          WRITE(IU06,*)'    SUB. INITMDL: BOUNDARY POINTS READ ',       &
     &     ' AND FILES CLOSED '
          CALL FLUSH (IU06)
        ENDIF
      ENDIF

      ENDIF ! .NOT. LLUNSTR

      IF (.NOT.ALLOCATED(IU19)) ALLOCATE(IU19(GBOUNC))


!*    2.2.* SET GRIB HEADERS FOR INPUTS/OUTPUTS
!          ------------------------------------
      LANAONLY=.FALSE.
      IF ((CDATEA.EQ.CDATEE).AND.(CDATEA.EQ.CDATEF)) LANAONLY=.TRUE.

      CALL SETMARSTYPE

      IF (.NOT. LGRHDIFS) THEN
!       FOR INTEGRATED PARAMETERS
        CALL PRESET_WGRIB_TEMPLATE("I",NGRIB_HANDLE_WAM_I)
!       FOR SPECTRA 
        CALL PRESET_WGRIB_TEMPLATE("S",NGRIB_HANDLE_WAM_S)
      ENDIF

      IF (MARSTYPE.EQ.'cf' .OR. MARSTYPE.EQ.'pf' .OR.                   &
     &                          MARSTYPE.EQ.'fc'     ) THEN
        IF (ALLOCATED(FAC0))   DEALLOCATE(FAC0)
        IF (ALLOCATED(FAC1))   DEALLOCATE(FAC1)
        IF (ALLOCATED(FAC2))   DEALLOCATE(FAC2)
        IF (ALLOCATED(FAC3))   DEALLOCATE(FAC3)
        IF (ALLOCATED(FAK))    DEALLOCATE(FAK)
        IF (ALLOCATED(FRHF))   DEALLOCATE(FRHF)
        IF (ALLOCATED(DFIMHF)) DEALLOCATE(DFIMHF)
      ENDIF

! ----------------------------------------------------------------------

!*    3. PRINT INITIAL CONDITIONS AS READ FROM PERPROCESSING.
!        ----------------------------------------------------

      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' WAVE MODEL GRID ORGANISATION:'
      WRITE(IU06,3002) ' SOUTHERNMOST LATITUDE IN GRID IS .......: ',   &
     & AMOSOP, ' DEGREE'
      WRITE(IU06,3002) ' NORTHERNMOST LATITUDE IN GRID IS .......: ',   &
     & AMONOP, ' DEGREE'
      WRITE(IU06,3002) ' WESTERNMOST LONGITUDE IN GRID IS .......: ',   &
     & AMOWEP, ' DEGREE'
      WRITE(IU06,3002) ' EASTERNMOST LONGITUDE IN GRID IS .......: ',   &
     & AMOEAP, ' DEGREE'
      WRITE(IU06,3002) ' LATITUDE INCREMENT IS ..................: ',   &
     & XDELLA, ' DEGREE'
      WRITE(IU06,3002) ' LONGITUDE INCREMENT IS .................: ',   &
     & XDELLO, ' DEGREE'
      WRITE(IU06,*) '  '
      WRITE(IU06,3003) ' TOTAL LENGTH OF EACH BLOCK .............: ',   &
     & NIBLO
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' SPECTRAL RESOLUTION:'
      WRITE(IU06,3003) ' TOTAL NUMBER OF DIRECTIONS .............: ', NANG 
      WRITE(IU06,3003) ' TOTAL NUMBER OF FREQUENCIES I/O & PROPAG: ', NFRE_RED
      WRITE(IU06,3003) ' TOTAL NUMBER OF FREQUENCIES FOR  PHYSICS: ', NFRE

      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' PARALLEL ORGANISATION : '
      WRITE(IU06,*) ' NPROC      : ', NPROC
      WRITE(IU06,*) ' MTHREADS   : ', MTHREADS
      WRITE(IU06,*) ' NPROMA_WAM : ', NPROMA_WAM
      WRITE(IU06,*) '  '
      CALL FLUSH (IU06)

 3002 FORMAT(3x,a,f9.3,a)
 3003 FORMAT(3x,a,i8,  a)

!NEST
      IF (IBOUNC .EQ. 1) THEN
        WRITE(IU06,*)
        WRITE(IU06,*) ' COARSE GRID: BOUNDARY OUTPUT POINTS :'
        WRITE(IU06,*) ' TOTAL NUMBER OF BOUNDARY POINTS IS: ',NBOUNC
        IF (ITEST.GE.2) THEN
          WRITE(IU06,'(/,4X,''BLOCK NO'',6X,''INDEX NO'',               &
     &     8X,''LONGITUDE'',6X,''LATITUDE'')')
          DO I=1,NBOUNC
            IX  = IXLG(IJARC(I))
            KX  = KXLT(IJARC(I))
            XLO = AMOWEP+REAL(IX-1,JWRB)*ZDELLO(KX)
            XLA = AMOSOP+REAL(KX-1,JWRB)*XDELLA

            WRITE(IU06,'((6X,I3,10X,I5,7X,F10.3,4X,F10.3))')            &
     &       IGARC(I), IJARC(I), XLO, XLA
          ENDDO
        ENDIF
      ENDIF
!NEST

!*    3.1 COMPUTE SHALLOW WATER TABLE INDICES.
!         ------------------------------------

!     BUILD DEPTH POINTER AND RESET DEPTH TO ALL POSITIVE
      IF(.NOT.ALLOCATED(IODP)) ALLOCATE(IODP(IJS:IJL))
      DO IJ=IJS,IJL
        IF(DEPTH(IJ).LE.TOOSHALLOW) THEN
          IODP(IJ) = 0
        ELSE
          IODP(IJ) = 1
        ENDIF
      ENDDO

      IF (.NOT.ALLOCATED(IOBND)) ALLOCATE(IOBND(IJS:IJL))
      IF (.NOT. LLUNSTR) THEN
         IOBND(:)=1
      ELSE
        DO IJ = IJS, IJL
          IF (IOBP(IJ) .NE. 0) THEN
            IOBND(IJ)=0
          ELSE
            IOBND(IJ)=1
          ENDIF
        ENDDO
      ENDIF

      DO IJ=NINF,NSUP
        DEPTH(IJ) = MAX(DEPTH(IJ),DEPTHA)
      ENDDO

!     COMPUTE THE MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH
      IF(.NOT.ALLOCATED(EMAXDPT)) ALLOCATE(EMAXDPT(IJS:IJL))
      DO IJ = IJS, IJL
!       REDUCE GAMMA FOR SMALL DEPTH ( < 4m)
!       (might need to be revisted when grid is fine resolution)
        IF (DEPTH(IJ).LT.4.0_JWRB) THEN
          GAM=GAM_B_J*DEPTH(IJ)/4.0_JWRB
        ELSE
          GAM=GAM_B_J
        ENDIF
        EMAXDPT(IJ)=0.0625_JWRB*(GAM*DEPTH(IJ))**2
      ENDDO


!     COMPUTE SHALLOW WATER TABLE INDICES AND THE RECIPROCAL PHASE VELOCITY.
!     INDEP must contain the halo points.
      IF (.NOT.ALLOCATED(INDEP)) ALLOCATE(INDEP(NINF-1:NSUP))
      IF (.NOT.ALLOCATED(CINV)) ALLOCATE(CINV(NDEPTH,NFRE))
      INDEP(NINF-1)=NDEPTH
      IF (ISHALLO.NE.1) THEN
        DO IJ=NINF,NSUP
          XD = LOG(DEPTH(IJ)/DEPTHA)/LOG(DEPTHD)+1.0_JWRB
          ID = NINT(XD)
          ID = MAX(ID,1)
          INDEP(IJ) = MIN(ID,NDEPTH)
        ENDDO

        DO M=1,NFRE
          FAC = ZPI*FR(M)
          DO JD=1,NDEPTH
            CINV(JD,M)=TFAK(JD,M)/FAC
          ENDDO
        ENDDO

      ELSE

        DO IJ=NINF,NSUP
          INDEP(IJ) = NDEPTH
        ENDDO

        DO M=1,NFRE
          FAC = ZPI*FR(M)
          DO JD=1,NDEPTH
            CINV(JD,M)=FAC/G
          ENDDO
        ENDDO

      ENDIF


! ----------------------------------------------------------------------

!*    4. CONNECT RESTART FIELDS TO OUTPUT UNITS (IF PBIO SOFTWARE NOT
!        USED). AND READ IN LAST WINDFIELDS FROM RESTARTFILE.
!        ---------------------------------------------------------------

      IF (ITEST.GE.2)CALL FLUSH(IU06)

      IF ( (LWCOU .AND. LWCUR ) .OR.                                    &
     &        IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN
        IF (.NOT.ALLOCATED(U)) ALLOCATE(U(NINF-1:NSUP))
        IF (.NOT.ALLOCATED(V)) ALLOCATE(V(NINF-1:NSUP))
      ENDIF

      Z0B(:) = 0.0_JWRB

      CALL GETSTRESS(U10OLD,THWOLD,USOLD,TAUW,TAUWDIR,Z0OLD,            &
     &               ROAIRO,ZIDLOLD,CICOVER,CITHICK,                    &
     &               NBLKS,NBLKE,IREAD)

      CDA = CDTPRO

      IF (LWCOU) LLWSWAVE = .FALSE.
      IF (LWCOU) LLWDWAVE = .FALSE.

      IF (CDTPRO.NE.ZERO .OR. LRESTARTED) THEN

!*    4.1 MODEL STARTS FROM FILES OUT OF A PREVIOUS MODEL RUN.
!         ----------------------------------------------------

        IF (CDTPRO.LT.CDATEA .OR. CDTPRO.GT.CDATEE .OR.                 &
     &      (IFORCA.EQ.1 .AND. CDTPRO.GT.CDATEF) .OR.                   &
     &      (IFORCA.NE.1 .AND. CDTPRO.LE.CDATEF)       ) THEN
          WRITE(IU06,*) ' *******************************************'
          WRITE(IU06,*) ' *    FATAL ERROR IN SUB. INITMDL          *'
          WRITE(IU06,*) ' *    ===========================          *'
          WRITE(IU06,*) ' * START DATE FROM RESTART FIELD IS NOT    *'
          WRITE(IU06,*) ' * MODEL PERIOD.                           *'
          IF (IFORCA.EQ.1) THEN
            WRITE(IU06,*) ' *  IN ANALYSIS PERIOD AS REQUESTED.       *'
          ELSE
            WRITE(IU06,*) ' *  IN FORECAST PERIOD AS REQUESTED.       *'
          ENDIF
          WRITE(IU06,*) ' * START DATE OF RUN       IS CDATEA = ',      &
     &     CDATEA
          WRITE(IU06,*) ' * START DATE OF FORECAST  IS CDATEF = ',      &
     &     CDATEF
          WRITE(IU06,*) ' * END   DATE OF RUN       IS CDATEE = ',      &
     &     CDATEE
          WRITE(IU06,*) ' * START DATE FROM RESTART IS CDTPRO = ',      &
     &     CDTPRO
          WRITE(IU06,*) ' *                                         *'
          WRITE(IU06,*) ' * PROGRAM ABORTS     PROGRAM ABORTS       *'
          WRITE(IU06,*) ' *                                         *'
          WRITE(IU06,*) ' *******************************************'
          CALL ABORT1
        ENDIF
      ELSE

!*    4.2 MODEL STARTS FROM FIELDS CREATED BY PRESET.
!         -------------------------------------------

        CDTPRO = CDATEA
      ENDIF

      IF (LWCOU) THEN
        IDELWO = KCOUSTEP 
        IDELWI = KCOUSTEP 
        IDELCUR= KCOUSTEP
      ELSE
        IF (NDELW_LST.GT.0) THEN
          DO IC=1,NDELW_LST
            IF (CDTPRO.LT.CDTW_LST(IC)) THEN
              IDELWI=IDELWI_LST(IC)
              IDELWO=IDELWO_LST(IC)
              EXIT
            ENDIF
          ENDDO
        ENDIF
      ENDIF

      CDATEWO = CDTPRO
      CDTBC= CDATEA
      IF (IDELT.LT.IDELWO) CALL INCDATE(CDATEWO,IDELWO/2)
      CDAWIFL = CDTPRO
      IDELWH = MAX(IDELWI,IDELPRO)
      CALL INCDATE(CDAWIFL,IDELWH)
      CDATEFL = CDATEWO

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '    SUB. INITMDL: WIND FIELD AND ',              &
     &   ' COUNTER INITIALIZED'
        WRITE(IU06,*) '      NEXT WINDFIELD WILL BE READ AT     ',      &
     &   'CDATEWO = ',CDATEWO
        WRITE(IU06,*) '      NEXT WIND FILE WILL BE ACCESSED AT ',      &
     &   'CDATEFL = ',CDATEFL
        WRITE(IU06,*) '      NEXT WIND FILE NAME IS FROM        ',      &
     &   'CDAWIFL = ',CDAWIFL
      ENDIF

! ----------------------------------------------------------------------

!*    5. INITIALIZE MODEL TIME VARIABLES
!        -------------------------------

!*    5.1 OUTPUT TIME VARIABLES.
!         ----------------------

      IF (FFLAG20 .OR. GFLAG20) THEN
        CDTINTT= CBPLTDT
        IF (LANAONLY) THEN
          CDTINTT=CDATEA
        ELSE
          DO WHILE (CDTINTT.LE.CDTPRO)
            CALL INCDATE (CDTINTT, IDELINT)
          ENDDO
        ENDIF
      ELSE
        CDTINTT = ZERO
      ENDIF

!*    5.2 FILE DISPOSE TIME AND RESTART TIME.
!         -----------------------------------

      CDTRES = CBPLTDT
      CDTBC  = CDATEA
      CALL INCDATE(CDTRES, IDELRES)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '    SUB. INITMDL: TIME COUNTER INITIALIZED'
        CALL FLUSH(IU06)
      ENDIF

! ----------------------------------------------------------------------

!*    6. CONSISTENCY CHECK ACCORDING TO CFL CRITERION.
!        ---------------------------------------------
      IF (.NOT. LLUNSTR) THEN

      GVE = G/(ZPI*FR(1)*2.0_JWRB)
      DPH = DELPHI
      IF (DPH.EQ.0.0_JWRB) THEN
        CFLP= 0.0_JWRB
      ELSE
        CFLP= IDELPRO*GVE/DPH
      ENDIF
      DLH = CIRC 
!     FIND THE SMALLEST DLH
      DO KX=KMSOP,KMNOP
        DLH_KX = DELLAM(KX)*COSPH(KX)
        DLH = MIN(DLH_KX,DLH)
      ENDDO

      IF (DLH.EQ.0.0_JWRB) THEN
        CFLL= 0.0_JWRB
      ELSE
        CFLL= IDELPRO*GVE/DLH
      ENDIF
      IF (CFLP.GT.1.0_JWRB .OR. CFLL.GT.1.0_JWRB) THEN
        WRITE(IU06,*) ' **********************************************'
        WRITE(IU06,*) ' *                                            *'
        WRITE(IU06,*) ' *       FATAL ERROR IN SUB. INITMDL          *'
        WRITE(IU06,*) ' *       ===========================          *'
        WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED.               *'
        WRITE(IU06,*) ' * CFLP: ',CFLP,'  GROUP VELOCITY: ',GVE
        WRITE(IU06,*) ' * CFLL: ',CFLL,'  GROUP VELOCITY: ',GVE
        WRITE(IU06,*) ' * PROPAGATION TIME: ',IDELPRO
        WRITE(IU06,*) ' * GRID Y-DISTANCE AT MOST NORTH LATITUDE: ',DPH
        WRITE(IU06,*) ' * GRID X-DISTANCE AT MOST NORTH LATITUDE: ',DLH
        WRITE(IU06,*) ' *     PROGRAM ABORTS   PROGRAM ABORTS        *'
        WRITE(IU06,*) ' *                                            *'
        WRITE(IU06,*) ' **********************************************'
        CALL ABORT1
      ENDIF

      ENDIF ! .NOT. LLUNSTR

! ----------------------------------------------------------------------

!*    7. NUMBER OF PROPAGATION TIME STEPS PER CALL.
!        ------------------------------------------

      NADV = IDELWI/IDELPRO
      NADV = MAX(NADV,1)
      IF (LANAONLY) THEN
        NADV = 0
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' **** NOTE ****'
        WRITE(IU06,*) ' NO FORWARD TIME INTEGRATION WILL BE PERFORMED'
        WRITE(IU06,*) ' ONLY THE DATA ASSIMILATION.'
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' '
      ENDIF

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '    SUB. INITMDL: NUMBER OF PROPAGATION STEPS'
        WRITE(IU06,*) '                  IN ONE CALL OF SUB WAVEMDL'
        WRITE(IU06,*) '                  WILL BE NADV = ', NADV
        CALL FLUSH(IU06)
      ENDIF

! ----------------------------------------------------------------------


!*    8.2  PRECOMPUTE BOTTOM REFRACTION TERMS.
!          -----------------------------------
      IF (IREFRA.NE.0) THEN
        NIBLD=NIBLO
        NBLD=1
        NIBLC=NIBLO
        NBLC=1
      ELSE
        NIBLD=0
        NBLD=0
        NIBLC=0
        NBLC=0
      ENDIF

!     INITIALISE CDTCUR
      CDTCUR=CDATECURA
      IF (.NOT.LWCOU .AND. .NOT.LRESTARTED) THEN
        IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN
          CALL INCDATE(CDTCUR,-IDELCUR)
        ENDIF
      ENDIF

!     COMPUTE BOTTOM REFRACTION TERMS
      IF (IREFRA .NE. 0) THEN
        IF (.NOT. LLUNSTR) THEN
!         ARRAY TO KEEP DEPTH AND CURRENT REFRACTION FOR THETA DOT
!         AND SIGMA DOT
          IF (.NOT.ALLOCATED(THDC)) ALLOCATE(THDC(IJS:IJL,NANG))
          IF (.NOT.ALLOCATED(THDD)) ALLOCATE(THDD(IJS:IJL,NANG))
          IF (.NOT.ALLOCATED(SDOT)) ALLOCATE(SDOT(IJS:IJL,NANG,NFRE))

          NPROMA=NPROMA_WAM
!$OMP     PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
          DO JKGLO = IJS, IJL, NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,IJL)
            CALL PROPDOT(KIJS, KIJL, THDC(KIJS:KIJL,:), THDD(KIJS:KIJL,:), SDOT(KIJS:KIJL,:,:))
          ENDDO
!$OMP     END PARALLEL DO

          IF (ITEST.GE.2) THEN
            WRITE(IU06,*) ' SUB. INITMDL: REFRACTION TERMS INITIALIZED '
            CALL FLUSH(IU06)
          END IF
        END IF
      ENDIF


!     INITIALIZE THE NEMO COUPLING
      IF (LWNEMOCOU) CALL INITNEMOCPL(LWNEMOCOURECV)


      IF (LLUNSTR) THEN
        IF (OUT_METHOD .eq. 1) THEN
          CALL INITIAL_OUTPUT_INITS
        END IF
      ENDIF

!     NOTE THAT CURRENT REFRACTION TERMS ARE NOW COMPUTED WHEN
!     *GETCURR* IS CALLED (SEE *PREWIND*).
!     A CALL TO PREWIND IS NEEDED TO GET THE SURFACE CURRENTS
!     BEFORE A CALL TO GETSPEC DUE TO THE TRANSFORMATION FROM
!     ABSOLUTE TO RELATIVE FRAME OF REFERENCE. 

      LLINIT=.NOT.LRESTARTED
      LLALLOC_FIELDG_ONLY=.FALSE.

      CALL PREWIND (U10OLD,THWOLD,USOLD,Z0OLD,                          &
     &              ROAIRO, ZIDLOLD,                                    &
     &              CICOVER, CITHICK,                                   &
     &              LLINIT, LLALLOC_FIELDG_ONLY,                        &
     &              IREAD,                                              &
     &              NFIELDS, NGPTOTG, NC, NR,                           &
     &              FIELDS, LWCUR, MASK_IN)

      WRITE(IU06,*) ' SUB. INITMDL: PREWIND DONE'                   
      CALL FLUSH (IU06)


!    GET SEA ICE DIMENSIONLESS ENERGY ATTENUATION COEFFICIENT
!!!! might need to restrict call when needed !!!
      CALL CIGETDEAC
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '  SUB. INITMDL:  CIGETDEAC CALLED '
        WRITE(IU06,*) ''
      ENDIF

!     DETERMINE THE SEA ICE REDUCTION FACTOR
      CALL CIREDUCE (CICOVER,CITHICK,CIWA)

      WRITE(IU06,*) ' SUB. INITMDL: CIREDUCE DONE'                   
      CALL FLUSH (IU06)


! ----------------------------------------------------------------------

!*    9.1 READ SPECTRA
!         ------------

      CALL GETSPEC(FL1,IJS,IJL,NBLKS,NBLKE,IREAD)

      WRITE(IU06,*) '    SUB. INITMDL: SPECTRA READ IN'
      CALL FLUSH (IU06)


!     9.2 COMPUTE FREQUENCY DEPENDENT INDICES AND COEFFICIENTS FOR SNONLIN
!         AND THE FREQUENCY FRONT TAIl REDUCTION COEFFICIENTS.
!         ------------------------------------------------------------ 

      CALL INISNONLIN


!     9.2 COMPUTE THE NONLINEAR TRANSFER FUNCTION COEFFIxCIENTS FOR SNL
!         VALID FOR THE NEW FORMULATION AND OTHER STUFFS
!         ------------------------------------------------------------ 

      IF (.NOT.ALLOCATED(ENH)) ALLOCATE(ENH(IJS:IJL,MLSTHG))

      IF (ISNONLIN.EQ.1) THEN
        IF (ISHALLO.NE.1) THEN
            DO M=1,NFRE
               DO IJ = IJS, IJL 
                 D = DEPTH(IJ)
                 OM = ZPI*FR(M)
                 XK = AKI(OM,D)
                 ENH(IJ,M) = MAX(MIN(ENH_MAX,TRANSF(XK,D)),ENH_MIN)
               ENDDO
            ENDDO
            DO M=NFRE+1,MLSTHG
               DO IJ = IJS, IJL 
                 D = DEPTH(IJ)
                 OM = ZPI*FR(NFRE)*FRATIO**(M-NFRE)
!                NOTE THAT TFAK IS NOT DEFINED BEYOND M=NFRE
!                HENCE THE USE OF FUNCTIOn AKI.
                 XK = AKI(OM,D)
                 ENH(IJ,M) = MAX(MIN(ENH_MAX,TRANSF(XK,D)),ENH_MIN)
               ENDDO
            ENDDO
        ELSE
          DO M=1,MLSTHG
             DO IJ = IJS, IJL 
               ENH(IJ,M) = 1.0_JWRB
             ENDDO
          ENDDO
        ENDIF
      ENDIF

!     9.3 DETERMINE WHETHER A GRID POINT IS SURROUNDED BY
!         POINTS WITH THE SAME DEPTH INDEX (EXCLUDING LAND POINTS).
!         ---------------------------------------------------------
      IF (.NOT. LLUNSTR) THEN

      IF (.NOT.ALLOCATED(LSAMEDEPTH))                                   &
     &    ALLOCATE(LSAMEDEPTH(IJS:IJL))

      IF (IPROPAGS.EQ.2) THEN
         DO IJ = IJS, IJL 
           IF (INDEP(IJ).EQ.INDEP(KLON(IJ,1))   .AND.                   &
     &         INDEP(IJ).EQ.INDEP(KLON(IJ,2))   .AND.                   &
     &         INDEP(IJ).EQ.INDEP(KLAT(IJ,1,1)) .AND.                   &
     &         INDEP(IJ).EQ.INDEP(KLAT(IJ,2,1)) .AND.                   &
     &         INDEP(IJ).EQ.INDEP(KLAT(IJ,1,2)) .AND.                   &
     &         INDEP(IJ).EQ.INDEP(KLAT(IJ,2,2))      ) THEN
             LSAMEDEPTH(IJ) = .TRUE.
           ELSE
             LSAMEDEPTH(IJ) = .FALSE.
           ENDIF
           IF (LSAMEDEPTH(IJ)) THEN
             OUTER : DO IC=1,2
               DO ICR=1,4
                 IF (INDEP(IJ).NE.INDEP(KCOR(IJ,ICR,IC))) THEN
                   LSAMEDEPTH(IJ) = .FALSE.
                   EXIT OUTER
                 ENDIF 
               ENDDO
             ENDDO OUTER
           ENDIF

         ENDDO
      ELSEIF (IPROPAGS.EQ.1) THEN
         DO IJ = IJS, IJL 
           IF (INDEP(IJ) .EQ. INDEP(KLON (IJ,1))   .AND.                &
     &         INDEP(IJ) .EQ. INDEP(KLON (IJ,2))   .AND.                &
     &         INDEP(IJ) .EQ. INDEP(KLAT (IJ,1,1)) .AND.                &
     &         INDEP(IJ) .EQ. INDEP(KLAT (IJ,2,1)) .AND.                &
     &         INDEP(IJ) .EQ. INDEP(KLAT (IJ,1,2)) .AND.                &
     &         INDEP(IJ) .EQ. INDEP(KLAT (IJ,2,2)) .AND.                &
     &         INDEP(IJ) .EQ. INDEP(KRLON(IJ,1,1)) .AND.                &
     &         INDEP(IJ) .EQ. INDEP(KRLON(IJ,2,1)) .AND.                &
     &         INDEP(IJ) .EQ. INDEP(KRLAT(IJ,1,1)) .AND.                &
     &         INDEP(IJ) .EQ. INDEP(KRLAT(IJ,2,1)) .AND.                &
     &         INDEP(IJ) .EQ. INDEP(KRLAT(IJ,1,2)) .AND.                &
     &         INDEP(IJ) .EQ. INDEP(KRLAT(IJ,2,2)) .AND.                &
     &         INDEP(IJ) .EQ. INDEP(KRLON(IJ,1,2)) .AND.                &
     &         INDEP(IJ) .EQ. INDEP(KRLON(IJ,2,2))      ) THEN
             LSAMEDEPTH(IJ) = .TRUE.
           ELSE
             LSAMEDEPTH(IJ) = .FALSE.
           ENDIF
         ENDDO
      ELSE
         DO IJ = IJS, IJL 
           IF (INDEP(IJ) .EQ. INDEP(KLON(IJ,1))   .AND.                 &
     &         INDEP(IJ) .EQ. INDEP(KLON(IJ,2))   .AND.                 &
     &         INDEP(IJ) .EQ. INDEP(KLAT(IJ,1,1)) .AND.                 &
     &         INDEP(IJ) .EQ. INDEP(KLAT(IJ,2,1)) .AND.                 &
     &         INDEP(IJ) .EQ. INDEP(KLAT(IJ,1,2)) .AND.                 &
     &         INDEP(IJ) .EQ. INDEP(KLAT(IJ,2,2))      ) THEN
             LSAMEDEPTH(IJ) = .TRUE.
           ELSE
             LSAMEDEPTH(IJ) = .FALSE.
           ENDIF
         ENDDO
      ENDIF

!     9.4 DEFINE JXO, JYO, KCR IF NEEDED
!         ------------------------------
      IF (IPROPAGS.EQ.2) THEN

        IF (.NOT. ALLOCATED(MPM)) ALLOCATE(MPM(NFRE,-1:1))
        DO M=1,NFRE
          MPM(M,-1)= MAX(1,M-1)
          MPM(M,0) = M
          MPM(M,1) = MIN(NFRE,M+1)
        ENDDO

        IF (.NOT. ALLOCATED(KPM)) ALLOCATE(KPM(NANG,-1:1))

        IF (.NOT. ALLOCATED(JXO)) ALLOCATE(JXO(NANG,2))
        IF (.NOT. ALLOCATED(JYO)) ALLOCATE(JYO(NANG,2))
        IF (.NOT. ALLOCATED(KCR)) ALLOCATE(KCR(NANG,4))
        DO K=1,NANG

          KM1 = K-1
          IF (KM1.LT.1) KM1 = NANG
          KPM(K,-1)=KM1

          KPM(K,0)=K

          KP1 = K+1
          IF (KP1.GT.NANG) KP1 = 1
          KPM(K,1)=KP1


          IF (COSTH(K).GE.0.0_JWRB) THEN
            JYO(K,1)=1
            JYO(K,2)=2
            IF (SINTH(K).GE.0.0_JWRB) THEN
              JXO(K,1)=1
              JXO(K,2)=2
              KCR(K,1)=3
              KCR(K,2)=2
              KCR(K,3)=4
              KCR(K,4)=1
            ELSE
              JXO(K,1)=2
              JXO(K,2)=1
              KCR(K,1)=2
              KCR(K,2)=3
              KCR(K,3)=1
              KCR(K,4)=4
            ENDIF
          ELSE
            JYO(K,1)=2
            JYO(K,2)=1
            IF (SINTH(K).GE.0.0_JWRB) THEN
              JXO(K,1)=1
              JXO(K,2)=2
              KCR(K,1)=4
              KCR(K,2)=1
              KCR(K,3)=3
              KCR(K,4)=2
            ELSE
              JXO(K,1)=2
              JXO(K,2)=1
              KCR(K,1)=1
              KCR(K,2)=4
              KCR(K,3)=2
              KCR(K,4)=3
            ENDIF
          ENDIF
        ENDDO
      ENDIF

      IF (ITEST.GE.2)                                                   &
     &   WRITE(IU06,*) '    SUB. INITMDL: MODULE YOWUBUF ',             &
     &   'READ AND INITIALIZED '


! ----------------------------------------------------------------------
!NEST
!     10. WRITE BOUNDARY VALUE FILE HEADER.
!         ------------------------------
      IF (IBOUNC.EQ.1) THEN
        IF ((LMESSPASS.AND.IRANK.EQ.1).OR..NOT.LMESSPASS) THEN
          DO II=1,GBOUNC
            IU19(II)=IWAM_GET_UNIT(IU06, CBCPREF(II), 'w', 'u', 0)
!           make the unit available for a silly fort.unit output
!           we will need to recode this a bit better !!!

            CLOSE(IU19(II))
            CALL HEADBC (IPOGBO(II)-IPOGBO(II-1),IDELPRO,TH(1),FR(1),   &
     &                   IU19(II), IU06)
            IF (ITEST.GE.2)                                             &
     &       WRITE(IU06,'(''    SUB. INITMDL: HEADER FOR  '',           &
     &       ''COARSE GRID WAS WRITTEN OF UNIT = '',I4)')               &
     &        IU19(II)
          ENDDO
        ENDIF
      ENDIF
!NEST

      ENDIF ! .NOT. LLUNSTR

      IF (ITEST.GE.2)                                                   &
     &   WRITE(IU06,*) '   SUB. INITMDL: DEFAULT  GRIB HEADERS SET'


      IF (LHOOK) CALL DR_HOOK('INITMDL',1,ZHOOK_HANDLE)

      END SUBROUTINE INITMDL
