SUBROUTINE INITMDL (NADV,                                         &
     &              IREAD,                                        &
     &              NFIELDS, NGPTOTG, NC, NR,                     &
     &              FIELDS, LWCUR, MASK_IN, PRPLRADI)

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
!                      NOTIM

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

!          SUB NOTIM OPENS FILES AND ASSIGNS THEM TO UNIT
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
     &            NGOUT    ,IJAR     ,NOUTT    ,LOUTINT
      USE YOWCURR  , ONLY : CDTCUR   ,IDELCUR  ,CDATECURA
      USE YOWFPBO  , ONLY : IBOUNF
      USE YOWGRIB_HANDLES , ONLY :NGRIB_HANDLE_WAM_I,NGRIB_HANDLE_WAM_S
      USE YOWFRED  , ONLY : FR       ,TH       ,DELTH   ,FR5      ,     &
     &            ZPIFR    ,                                            &
     &            FRM5     ,COFRM4   ,COEF4    ,FRATIO  ,FLOGSPRDM1,    &
     &            FLMAX    ,RHOWG_DFIM,                                 &
     &            DFIM     ,DFIMOFR  ,DFIMFR  ,DFIMFR2   ,              &
     &            DFIM_SIM ,DFIMOFR_SIM ,DFIMFR_SIM ,DFIMFR2_SIM ,      &
     &            DFIM_END_L, DFIM_END_U
      USE YOWGRIBHD, ONLY : LGRHDIFS
      USE YOWGRID  , ONLY : DELPHI   ,DELLAM   ,IJS     ,IJL      ,COSPH
      USE YOWICE   , ONLY : CICOVER  ,CITHICK  ,CIWA
      USE YOWMAP   , ONLY : AMOWEP   ,AMOSOP   ,                        &
     &            AMOEAP   ,AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,    &
     &            KMNOP    ,KMSOP    ,IPER
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,KTAG
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_RED ,NFRE_ODD ,    & 
     &            NGX      ,NGY      ,                                  &
     &            NIBLO    ,NIBLD    ,NIBLC
      USE YOWPCONS , ONLY : G        ,CIRC     ,PI       ,ZPI      ,    &
     &            RAD      ,ROWATER  ,ZPI4GM2  ,FM2FP
      USE YOWPHYS  , ONLY : ALPHAPMAX, ALPHAPMINFAC, FLMINFAC
      USE YOWREFD  , ONLY : LLUPDTTD

      USE YOWSHAL  , ONLY : NDEPTH   ,DEPTH    ,DEPTHA   ,DEPTHD   ,    &
!                                   !!!!!!!!!!

     &            INDEP    ,CGROUP   ,IODP     ,IOBND    ,TOOSHALLOW
      USE YOWSPEC  , ONLY : NBLKS    ,NBLKE    ,KLENTOP  ,KLENBOT  ,    &
     &            WSWAVE   ,WDWAVE   ,UFRIC    ,Z0M    ,TAUW     ,    &
     &            Z0B      ,TAUWDIR,  AIRD   ,WSTAR  ,             &
     &            FL1
      USE YOWSTAT  , ONLY : CDATEE   ,CDATEF   ,CDTPRO   ,CDTRES   ,    &
     &            CDTINTT  ,CDTBC    ,                                  &
     &            IDELPRO  ,IDELT    ,IDELWI   ,IDELWO   ,IDELRES  ,    &
     &            IDELINT  ,                                            &
     &            IREFRA   ,                                            &
     &            IPHYS    ,                                            &
     &            NPROMA_WAM,                                           &
     &            CDATEA   ,MARSTYPE ,LANAONLY ,ISNONLIN ,IPROPAGS ,    &
     &            IDELWI_LST,IDELWO_LST,CDTW_LST,NDELW_LST
      USE YOWTABL  , ONLY : FAC0     ,FAC1     ,FAC2     ,FAC3     ,    &
     &            FAK      ,FRHF     ,DFIMHF
      USE YOWTEST  , ONLY : IU06
      USE YOWTEXT  , ONLY : LRESTARTED
      USE YOWWNDG  , ONLY : ICODE
      USE YOWUBUF  , ONLY : LUPDTWGHT
      USE YOWUNIT  , ONLY : IU02     ,IU11     ,IU12     ,              &
     &            IU13     ,IU14     ,IU15     ,IU16     ,IU17     ,    &
     &            IU18     ,IU19     ,IU20     ,IU22     ,              &
     &            IU23     ,IU27     ,IU28     ,                        &
     &            IUVELO   ,IU30     ,IU31     ,                        &
     &            IU32     ,IU33     ,IU35     ,IU36     ,IU37     ,    &
     &            IU38
      USE YOWWAMI  , ONLY : CBPLTDT
      USE YOWWIND  , ONLY : CDATEWL  ,CDAWIFL  ,CDATEWO  ,CDATEFL  ,    &
     &            LLNEWCURR,LLWSWAVE ,LLWDWAVE ,FF_NEXT
      USE YOWUNPOOL, ONLY : LLUNSTR, OUT_METHOD
      USE UNSTRUCT_BOUND , ONLY : IOBP
      USE OUTPUT_STRUCT, ONLY : INITIAL_OUTPUT_INITS

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
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
#include "initdpthflds.intfb.h"
#include "initnemocpl.intfb.h"
#include "iniwcst.intfb.h"
#include "preset_wgrib_template.intfb.h"
#include "prewind.intfb.h"
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
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

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

      LOUTINT=.FALSE.

      LUPDTWGHT=.TRUE.

      LLNEWCURR=.TRUE.

      ICODE=0

!*    1.2 INPUT OF USER PARAMETER.
!         ------------------------

      CALL USERIN (IFORCA, LWCUR)

!     DEFINE COEFFICIENT FOR MEAN PERIODS CALCULATION
      IF (ALLOCATED(DFIMOFR)) DEALLOCATE(DFIMOFR)
      ALLOCATE(DFIMOFR(NFRE))
      IF (ALLOCATED(DFIMFR)) DEALLOCATE(DFIMFR)
      ALLOCATE(DFIMFR(NFRE))
      IF (ALLOCATED(DFIMFR2)) DEALLOCATE(DFIMFR2)
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

      IF (.NOT.ALLOCATED(DFIMOFR_SIM)) ALLOCATE(DFIMOFR_SIM(NFRE))
      IF (.NOT.ALLOCATED(DFIMFR_SIM)) ALLOCATE(DFIMFR_SIM(NFRE))
      IF (.NOT.ALLOCATED(DFIMFR2_SIM)) ALLOCATE(DFIMFR2_SIM(NFRE))
      DO M=1,NFRE
        DFIMOFR_SIM(M) = DFIM_SIM(M)/FR(M)
        DFIMFR_SIM(M)  = DFIM_SIM(M)*FR(M)
        DFIMFR2_SIM(M) = DFIM_SIM(M)*FR(M)**2
      ENDDO


      CALL TABU_SWELLFT

      CALL INIT_X0TAUHF

      KTAG=100

!     1.5 DETERMINE LAST OUTPUT DATE
!         --------------------------

      IF (NOUTT > 0) THEN
        COUTLST=COUTT(NOUTT)
      ELSE
        COUTLST=CDATEA
        IF (GFLAG20) THEN
          DO WHILE (COUTLST <= CDATEE.AND.IDELINT > 0)
            CALL INCDATE (COUTLST,IDELINT)
          ENDDO
          CALL INCDATE (COUTLST,-IDELINT)
        ENDIF

        IF (COUTLST > CDATEE) COUTLST=CDATEE
      ENDIF

! ----------------------------------------------------------------------

!*    2. INPUT FROM PREPROCESSING PROGRAMS.
!        ----------------------------------

!     THE ACTUAL READING HAS BEEN MOVED TO MPDECOMP.

!*    2.1 READ MODULE YOWCPBO AND YOWFPBO.
!     ------------------------------------

      GBOUNC = 1

      IF (.NOT. LLUNSTR) THEN
 
      IF (IBOUNC == 1 .OR. IBOUNF == 1) THEN
        IF (IBOUNC == 1) THEN
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

        IF (IBOUNF == 1) THEN
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
      ENDIF

      ENDIF ! .NOT. LLUNSTR

      IF (.NOT.ALLOCATED(IU19)) ALLOCATE(IU19(GBOUNC))


!*    2.2.* SET GRIB HEADERS FOR INPUTS/OUTPUTS
!          ------------------------------------
      LANAONLY=.FALSE.
      IF ((CDATEA == CDATEE).AND.(CDATEA == CDATEF)) LANAONLY=.TRUE.

      CALL SETMARSTYPE

      IF (.NOT. LGRHDIFS) THEN
!       FOR INTEGRATED PARAMETERS
        CALL PRESET_WGRIB_TEMPLATE("I",NGRIB_HANDLE_WAM_I)
!       FOR SPECTRA 
        CALL PRESET_WGRIB_TEMPLATE("S",NGRIB_HANDLE_WAM_S)
      ENDIF

      IF (MARSTYPE == 'cf' .OR. MARSTYPE == 'pf' .OR. MARSTYPE == 'fc' ) THEN
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
      IF (IBOUNC == 1) THEN
        WRITE(IU06,*)
        WRITE(IU06,*) ' COARSE GRID: BOUNDARY OUTPUT POINTS :'
        WRITE(IU06,*) ' TOTAL NUMBER OF BOUNDARY POINTS IS: ',NBOUNC
      ENDIF
!NEST

!*    3.1 COMPUTE SHALLOW WATER TABLE INDICES.
!         ------------------------------------

!     BUILD DEPTH POINTER AND RESET DEPTH TO ALL POSITIVE
      IF (.NOT.ALLOCATED(IODP)) ALLOCATE(IODP(IJS:IJL))
      DO IJ=IJS,IJL
        IF (DEPTH(IJ) <= TOOSHALLOW) THEN
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
          IF (IOBP(IJ) /= 0) THEN
            IOBND(IJ)=0
          ELSE
            IOBND(IJ)=1
          ENDIF
        ENDDO
      ENDIF

!     SET DEPTH MINIMUM
      DO IJ = IJS, IJL
        DEPTH(IJ) = MAX(DEPTH(IJ),DEPTHA)
      ENDDO

!     COMPUTE INDEP
      IF (.NOT.ALLOCATED(INDEP)) ALLOCATE(INDEP(IJS:IJL))
      DO IJ=IJS,IJL
        XD = LOG(DEPTH(IJ)/DEPTHA)/LOG(DEPTHD)+1.0_JWRB
        ID = NINT(XD)
        ID = MAX(ID,1)
        INDEP(IJ) = MIN(ID,NDEPTH)
      ENDDO


!     INITIALISE GRID POINT FIELDS DEPENDENT ON WATER DEPTH AND FREQUENCY
      CALL INITDPTHFLDS


!     INITIALISATION FOR SDISS_ARD
      IF (IPHYS == 1) CALL INIT_SDISS_ARDH


! ----------------------------------------------------------------------

!*    4. CONNECT RESTART FIELDS TO OUTPUT UNITS (IF PBIO SOFTWARE NOT
!        USED). AND READ IN LAST WINDFIELDS FROM RESTARTFILE.
!        ---------------------------------------------------------------

      Z0B(:) = 0.0_JWRB

      CALL GETSTRESS(WSWAVE, WDWAVE, UFRIC, TAUW, TAUWDIR, Z0M,       &
     &               AIRD, WSTAR, CICOVER, CITHICK,                &
     &               NBLKS, NBLKE, IREAD)

      CDATEWL = CDTPRO

      IF (LWCOU) LLWSWAVE = .FALSE.
      IF (LWCOU) LLWDWAVE = .FALSE.

      IF (CDTPRO /= ZERO .OR. LRESTARTED) THEN

!*    4.1 MODEL STARTS FROM FILES OUT OF A PREVIOUS MODEL RUN.
!         ----------------------------------------------------

        IF (CDTPRO < CDATEA .OR. CDTPRO > CDATEE .OR.                 &
     &      (IFORCA == 1 .AND. CDTPRO > CDATEF) .OR.                  &
     &      (IFORCA /= 1 .AND. CDTPRO <= CDATEF)       ) THEN
          WRITE(IU06,*) ' *******************************************'
          WRITE(IU06,*) ' *    FATAL ERROR IN SUB. INITMDL          *'
          WRITE(IU06,*) ' *    ===========================          *'
          WRITE(IU06,*) ' * START DATE FROM RESTART FIELD IS NOT    *'
          WRITE(IU06,*) ' * MODEL PERIOD.                           *'
          IF (IFORCA == 1) THEN
            WRITE(IU06,*) ' *  IN ANALYSIS PERIOD AS REQUESTED.       *'
          ELSE
            WRITE(IU06,*) ' *  IN FORECAST PERIOD AS REQUESTED.       *'
          ENDIF
          WRITE(IU06,*) ' * START DATE OF RUN       IS CDATEA = ', CDATEA
          WRITE(IU06,*) ' * START DATE OF FORECAST  IS CDATEF = ', CDATEF
          WRITE(IU06,*) ' * END   DATE OF RUN       IS CDATEE = ', CDATEE
          WRITE(IU06,*) ' * START DATE FROM RESTART IS CDTPRO = ', CDTPRO
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
        IF (NDELW_LST > 0) THEN
          DO IC=1,NDELW_LST
            IF (CDTPRO < CDTW_LST(IC)) THEN
              IDELWI=IDELWI_LST(IC)
              IDELWO=IDELWO_LST(IC)
              EXIT
            ENDIF
          ENDDO
        ENDIF
      ENDIF

      CDATEWO = CDTPRO
      CDTBC= CDATEA
      IF (IDELT < IDELWO) CALL INCDATE(CDATEWO,IDELWO/2)
      CDAWIFL = CDTPRO
      IDELWH = MAX(IDELWI,IDELPRO)
      CALL INCDATE(CDAWIFL,IDELWH)
      CDATEFL = CDATEWO

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

! ----------------------------------------------------------------------

!*    6. CONSISTENCY CHECK ACCORDING TO CFL CRITERION.
!        ---------------------------------------------
      IF (.NOT. LLUNSTR) THEN

      GVE = G/(ZPI*FR(1)*2.0_JWRB)
      DPH = DELPHI
      IF (DPH == 0.0_JWRB) THEN
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

      IF (DLH == 0.0_JWRB) THEN
        CFLL= 0.0_JWRB
      ELSE
        CFLL= IDELPRO*GVE/DLH
      ENDIF
      IF (CFLP > 1.0_JWRB .OR. CFLL > 1.0_JWRB) THEN
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

! ----------------------------------------------------------------------


!*    8.2  PRECOMPUTE BOTTOM REFRACTION TERMS.
!          -----------------------------------
      IF (IREFRA /= 0) THEN
        NIBLD=NIBLO
        NIBLC=NIBLO
      ELSE
        NIBLD=0
        NIBLC=0
      ENDIF

!     INITIALISE CDTCUR
      CDTCUR=CDATECURA
      IF (.NOT.LWCOU .AND. .NOT.LRESTARTED) THEN
        IF (IREFRA == 2 .OR. IREFRA == 3) THEN
          CALL INCDATE(CDTCUR,-IDELCUR)
        ENDIF
      ENDIF

!     COMPUTE BOTTOM REFRACTION TERMS (SEE *PROPAGS_WAM*)
      IF (IREFRA /= 0) THEN
        LLUPDTTD = .TRUE.
      ELSE
        LLUPDTTD = .FALSE.
      ENDIF

!     INITIALIZE THE NEMO COUPLING
      IF (LWNEMOCOU) CALL INITNEMOCPL(LWNEMOCOURECV)


      IF (LLUNSTR) THEN
        IF (OUT_METHOD == 1) THEN
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

      CALL PREWIND (WSWAVE, WDWAVE, UFRIC, Z0M,                       &
     &              AIRD, WSTAR,                                      &
     &              CICOVER, CITHICK,                                 &
     &              FF_NEXT,                                          &
     &              LLINIT, LLALLOC_FIELDG_ONLY,                      &
     &              IREAD,                                            &
     &              NFIELDS, NGPTOTG, NC, NR,                         &
     &              FIELDS, LWCUR, MASK_IN)

      WRITE(IU06,*) ' SUB. INITMDL: PREWIND DONE'                   
      CALL FLUSH (IU06)


!    GET SEA ICE DIMENSIONLESS ENERGY ATTENUATION COEFFICIENT
!!!! might need to restrict call when needed !!!
      CALL CIGETDEAC

!     DETERMINE THE SEA ICE REDUCTION FACTOR
      CALL CIREDUCE (IJS, IJL, CGROUP, CICOVER, CITHICK, CIWA)

      WRITE(IU06,*) ' SUB. INITMDL: CIREDUCE DONE'                   
      CALL FLUSH (IU06)

! ----------------------------------------------------------------------

!*    9.1 READ SPECTRA
!         ------------

      CALL GETSPEC(FL1, IJS, IJL, NBLKS, NBLKE, IREAD)

      WRITE(IU06,*) '    SUB. INITMDL: SPECTRA READ IN'
      CALL FLUSH (IU06)


!     9.2 COMPUTE FREQUENCY DEPENDENT INDICES AND COEFFICIENTS FOR SNONLIN
!         AND THE FREQUENCY FRONT TAIl REDUCTION COEFFICIENTS.
!         ------------------------------------------------------------ 

      CALL INISNONLIN

! ----------------------------------------------------------------------
!NEST
!     10. WRITE BOUNDARY VALUE FILE HEADER.
!         ------------------------------
      IF (IBOUNC == 1 .AND. .NOT. LLUNSTR) THEN
        IF (IRANK == 1) THEN
          DO II=1,GBOUNC
            IU19(II)=IWAM_GET_UNIT(IU06, CBCPREF(II), 'w', 'u', 0)
!           make the unit available for a silly fort.unit output
!           we will need to recode this a bit better !!!

            CLOSE(IU19(II))
            CALL HEADBC (IPOGBO(II)-IPOGBO(II-1),IDELPRO,TH(1),FR(1),   &
     &                   IU19(II), IU06)
          ENDDO
        ENDIF
      ENDIF
!NEST


      IF (LHOOK) CALL DR_HOOK('INITMDL',1,ZHOOK_HANDLE)

END SUBROUTINE INITMDL
