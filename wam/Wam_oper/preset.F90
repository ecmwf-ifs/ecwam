      PROGRAM preset

! ----------------------------------------------------------------------

!**** *PRESET* - GENERATES ALL BINARY FILES REQUIRED FOR A WAMODEL START
!             OR GENERATES SPECTRAL GRIB FILE REQUIRED FOR A WAMODEL
!             START. IN THAT CASE NO DRAG COEFFICIENT FIELD IS REQUIRED
!             TO RUN WAMODEL (SEE LNOCDIN IN WAMODEL).
!             OR CONVERTS AND INTERPOLATES GRIB SPECTRA PARAMETER 250 
!             INTO PARAMETER 251. THE CORRESPONDING DRAG COEFFICIENT 
!             FIELD SHOULD BE PROVIDED WITH AN APPROPRIATE MARS REQUEST.

!     SUSANNE HASSELMANN  MPI     JULY 1986.
!     ANNEGRET SPEIDEL    MPI      MAY 1988 PARAMETER STATEMENTS.
!     ANNEGRET SPEIDEL    MPI NOVEMBER 1988 CRAY-2 VERSION.
!     CYCLE_3 MODICIFATIONS:
!     ----------------------
!     RENATE PORTZ       MPI      JUNE 1990 COMPUTATION OF INITIAL
!                                           JONSWAP SPECTRA FROM
!                                           INITIAL WIND FIELD.
!     CYCLE_4 MODICIFATIONS:
!     ----------------------
!     H. GUNTHER  GKSS/ECMWF  DECEMBER 1990
!     J. BIDLOT    ECMWF   FEBRUARY 1996  MESSAGE PASSING
!     J. BIDLOT    ECMWF   MARCH 1997  MODIFY ROUTINES FOR OUTPUT OF
!                          RESTART FILE 
!     B. HANSEN    ECMWF   APRIL    1997  RESTART FACILITY.
!     B. HANSEN    ECMWF   JANUARY  1998  NAMELIST INPUT.
!     B. HANSEN    ECMWF   FEBRUARY 1998  WRITE 2DSP TO FDB.
!     J. BIDLOT    ECMWF   MARCH 1998 PRODUCE GRIB SPECTRA (215) AND
!                                     INTERPOLATE TO OUTPUT GRID.
!                                     NO LAW FILE IS PRODUCED THEN.
!     J. BIDLOT    ECMWF   OCTOBER 1998   MODULES.

!*    PURPOSE.
!     --------
!       TO INITIALISE ALL FILES REQUESTED BY THE WAMODEL.
!**   INTERFACE.
!     ----------
!       *IU01*   INTEGER    INPUT  UNIT UNBLOCKED WIND FILE.
!                           (SEE SUB READWIND).
!       *IU05*   INTEGER    USER INPUT UNIT.
!       *IU06*   INTEGER    PRINTER OUTPUT.
!       *IU07*   INTEGER    INPUT  UNIT PREPROC GRID OUTPUT.
!       *IU12*   INTEGER    OUTPUT UNIT BLOCKS OF SPECTRA.
!       *IU14*   INTEGER    OUTPUT UNIT SECOND LAT OF BLOCKS.
!       *IU15*   INTEGER    OUTPUT UNIT LAST WINDFIELDS.
!                 !OR!      OUTPUT UNIT LAST DRAG COEFFICIENT.
!     METHOD.
!     -------
!       IT IS USED TO PRODUCED A COLD START:

!       A JONSWAP SPECTRAL SHAPE IS ASSUMED.
!       JONSWAP PARAMETERS ARE DEFINED EITHER BY USER INPUT (IOPTI=0) OR
!       BY FETCH LAWS (IOPTI=1, or 2). IOPTI=3 CAN BE USED TO SPECIFY
!       INITIAL LOCALISED SWELL SYSTEMS (SEE MSWELL) 
!       THE 2-D SPECTRA ARE COMPUTED FOR
!       EACH POINT IN A BLOCK, THE WAMODEL BLOCKS ARE INITIALISED BY
!       THESE SPECTRA AND ALL BLOCKS AND OVERLAPPING LATITUDES ARE SAVED
!       IF FETCH LAWS ARE USED TO DEFINE PARAMETERS THE FIRST WIND
!       FIELD IS GENERATED OTHERWISE THE WIND FIELD IS INITIALISED
!       WITH ZEROS. THE MAIN MODEL WILL RECONSTRUCT THE WIND ANY HOW.
!       THE FILE HANDLING OF THE RESTART FILES IS COMPUTER DEPENDENT.
!       SUB GSFILE HAS TO BE MODIFIED, TO COPY THE UNIT ALIAS FILES
!       (UNITS IU12,IU14, AND IU15) TO PERMANENT FILES.


!     EXTERNALS.
!     ----------
!       *ABORT1*    - TERMINATES PROCESSING.
!       *BUILDSTRESS- RECONSTRUCT THE LAW FILE FROM WINDS AND CD. 
!       *GRSTNAME*  - BUILD NAME FOR INPUT/OUTPUT FILES 
!       *GSFILE*    - GETS OR SAVES FILES (COMPUTER DEPENDENT).
!       *INCDATE*   - UPDATES A DATE TINE GROUP.
!       *INIWCST*   - SETS GLOBAL CONSTANTS.
!       *LOCINT*    - INTERPOLATES TO MODEL GRID.
!       *MSTART*    - GENERATES THE RESTART FILES.
!       *NOTIM*     - CONTROLS WIND GENERATION (NO TIME INTERPOLATION).
!       *PEAK*      - COMPUTES PARAMETERS BY FETCH LAWS.
!       *PREWIND*   - PREPARES WINDS.
!       *SAVSPEC*   - SAVES THE RESTART SPECTRA.
!       *SAVSTRESS* - SAVES THE RESTART WIND AND STRESS FIELDS.
!       *SPECTRA*   - COMPUTES SPECTRA FROM PARAMETERS.
!       *SPR*       - DIRECTIONAL DISTRIBUTION.
!       *READPRE*   - READS PREPROC OUTPUT.
!       *READWIND*   - READS A WIND FIELD.

!     REFERENCES
!     ----------
!       NONE.
! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRIB_HANDLES , ONLY :NGRIB_HANDLE_WAM_I,NGRIB_HANDLE_WAM_S
      USE YOWGRIBHD, ONLY : PPMISS   ,PPEPS    ,PPREC    ,NTENCODE ,    &
     &            NGRBRESS ,HOPERS   ,PPRESOL  ,LGRHDIFS ,LNEWLVTP
      USE YOWGRID  , ONLY : DELPHI   ,NLONRGG  ,IJS      , IJL,         &
     &            IJSLOC   ,IJLLOC   ,IJGLOBAL_OFFSET
      USE YOWICE   , ONLY : CICOVER  ,CITHICK
      USE YOWJONS  , ONLY : FM       ,ALFA     ,GAMMA    ,SA       ,    &
     &            SB       ,THETAQ
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,IRGG     ,AMOWEP   ,    &
     &            AMOSOP   ,AMOEAP   ,AMONOP   ,XDELLA   ,XDELLO   ,    &
     &            IFROMIJ  ,JFROMIJ
      USE YOWMESPAS, ONLY : LMESSPASS,LFDBIOOUT,LGRIBOUT
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP     ,    &
     &            KTAG     ,NPRECR   ,NPRECI
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      ,    &
     &            NIBLO    ,SWAMPWIND,CLDOMAIN ,LL1D
      USE YOWPCONS , ONLY : G        ,RAD      ,DEG      ,ZMISS    ,    &
     &            ROAIR
      USE YOWSTAT  , ONLY : MARSTYPE ,YCLASS   ,YEXPVER  ,CDATEA   ,    &
     &            CDATEE   ,CDATEF   ,CDTPRO   ,CDATER   ,CDATES   ,    &
     &            IDELPRO  ,IDELWI   ,IDELWO   ,                        &
     &            NENSFNB  ,NTOTENS  ,NSYSNB   ,NMETNB   ,NPROMA_WAM,   &
     &            IREFDATE ,ISTREAM  ,NLOCGRB  ,IREFRA
      USE YOWSPEC  , ONLY : NSTART   ,NEND     ,U10OLD   ,THWOLD   ,    &
     &            USOLD    ,Z0OLD    ,TAUW     ,TAUWDIR  ,FL1      ,    &
     &            ROAIRO   ,WSTAROLD  ,NBLKS    ,NBLKE
      USE YOWTABL  , ONLY :  FAC0     ,FAC1     ,FAC2     ,FAC3    ,    &
     &            FAK      ,FRHF      ,DFIMHF    , OMEGA   ,THH     ,   &
     &            DFDTH    ,IM_P      ,IM_M     ,TA       ,TB      ,    &
     &            TC_QL    ,TT_4M     ,TT_4P    ,TFAKH

      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWTEXT  , ONLY : ICPLEN   ,USERID   ,RUNID    ,PATH     ,    &
     &            CPATH
      USE YOWUNPOOL ,ONLY : LLUNSTR  ,LPREPROC
      USE YOWUNIT  , ONLY : IU12     ,IU14     ,IU15
      USE YOWWIND  , ONLY : CDATEWL  ,CDAWIFL  ,CDATEWO  ,CDATEFL  ,    &
     &            LLNEWCURR,NXFF     ,NYFF     ,WSPMIN   ,FF_NEXT
      USE UNSTRUCT_BOUND, ONLY : IOBPD

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
      USE MPL_MODULE,ONLY : MPL_INIT, MPL_END

! -------------------------------------------------------------------

      IMPLICIT NONE 
#include "abort1.intfb.h"
#include "cigetdeac.intfb.h"
#include "iniwcst.intfb.h"
#include "mstart.intfb.h"
#include "mswell.intfb.h"
#include "outspec.intfb.h"
#include "preset_wgrib_template.intfb.h"
#include "prewind.intfb.h"
#include "readpre.intfb.h"
#include "savspec.intfb.h"
#include "savstress.intfb.h"

      INTEGER(KIND=JWIM), PARAMETER :: NC=1
      INTEGER(KIND=JWIM), PARAMETER :: NR=1
      INTEGER(KIND=JWIM), PARAMETER :: NGPTOTG=NC*NR
      INTEGER(KIND=JWIM), PARAMETER :: NFIELDS=1
      INTEGER(KIND=JWIM) :: ILEN, IREAD, IOPTI
      INTEGER(KIND=JWIM) :: IJ, K, M
      INTEGER(KIND=JWIM) :: IU05, IU07 
      INTEGER(KIND=JWIM) :: I4(2)
      INTEGER(KIND=JWIM) :: MASK_IN(NGPTOTG)
      INTEGER(KIND=JWIM) :: IWAM_GET_UNIT

      REAL(KIND=JWRB) :: PRPLRADI
      REAL(KIND=JWRB) :: THETA, FETCH, FRMAX 
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: X4(2)
      REAL(KIND=JWRB) :: FIELDS(NGPTOTG,NFIELDS)

      CHARACTER(LEN=1) :: CLTUNIT
      CHARACTER(LEN=14) :: ZERO, CDUM
      CHARACTER(LEN=70) :: HEADER
      CHARACTER(LEN=120) :: SFILENAME, ISFILENAME

      PARAMETER (ZERO=' ', CDUM='00000000000000')

      LOGICAL :: LLINIT
      LOGICAL :: LLALLOC_FIELDG_ONLY
      LOGICAL :: LWCUR

! ----------------------------------------------------------------------

      NAMELIST /NALINE/ HEADER,                                         &
     &          IOPTI, ITEST, ITESTB,                                   &
     &          ALFA, FM, GAMMA, SA, SB, THETA, FETCH, SWAMPWIND ,      &
     &          USERID, RUNID, PATH, CPATH,                             &
     &          CDATEA, IDELWI, CLTUNIT,                                &
     &          LLUNSTR, LPREPROC,                                      &
     &          LGRIBOUT,                                               &
     &          MARSTYPE, YCLASS, YEXPVER, NPROMA_WAM

!     IOPTI : IT SELECTS COLD START SPECTRAL FORM
!             (0, 1, or 2 see MSTART)
!             (3 see MSWELL)
!     ITEST, ITESTB : TEST OUTPUT LEVEL SEE MODULE YOWTEST
!     ALFA, FM, GAMMA, SA, SB, THETA, FETCH : JONSWAP PARAMETERS
!     SWAMPWIND : CONSTANT WIND SPEED USED BY SWAMP CASE.
!     USERID, RUNID, PATH : OUT OF DATE, SHOULD ONLY BE USED WHENEVER 
!                           USER WANTS TO WRITE OUTPUT DIRECTLY TO ECFS
!                           WHICH IS HIGHLY NOT RECOMMENDATED ON THE VPP
!     CPATH : PATH FOR OUTPUT TO DISK
!     CDATEA : INPUT DATE.
!     IDELWI : INPUT WIND TIME STEP IN UNIT DEFINED BY CLTUNIT
!     CLTUNIT : INPUT WIND TIME STEP UNIT (S : seconds or H : hours)
!     LGRIBOUT : IF TRUE THE WAVE SPECTRA IS OUTPUT IN GRIB FORMAT, ELSE
!                THE BINARY RESTART FILES ARE PRODUCED.
!     MARSTYPE : DATA TYPE USED TO CODE DATA IN GRIB
!     YCLASS   : DATA CLASS USED TO CODE DATA IN GRIB.
!     YEXPVER : EXPERIMENT VERSION USED TO CODE DATA IN GRIB.
!     NPROMA_WAM: NUMBER OF THE GRID POINTS PER LOOP WHEN CUT INTO CHUNKS


! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('PRESET',0,ZHOOK_HANDLE)

      CALL MPL_INIT(KOUTPUT=1)

      PRPLRADI=1.0_JWRB
      CALL INIWCST(PRPLRADI)

!*    0. SET DEFAULT VALUES FOR THE NAMELIST ELEMENTS.
!        ---------------------------------------------

      HEADER = ZERO
      IOPTI  =    1
      ITEST  =   -9
      ITESTB =   -9
      ALFA   =    0.0_JWRB
      FM     =    0.0_JWRB
      GAMMA  =    0.0_JWRB
      SA     =    0.0_JWRB
      SB     =    0.0_JWRB
      THETA  =    0.0_JWRB
      FETCH  =    0.0_JWRB
      SWAMPWIND = 18.45_JWRB
      USERID = ZERO
      RUNID  = ZERO
      PATH   = ZERO
      CPATH  = ZERO
      CDATEA = ZERO
      CLTUNIT= 'H' 
      IDELWI =    0

      LLUNSTR  =.FALSE.
      LPREPROC =.FALSE.

      LMESSPASS = .FALSE.

      LGRIBOUT = .TRUE.

      MARSTYPE = 'an'
      YCLASS   = 'rd'
      YEXPVER  = USERID//'a'

      NPROMA_WAM = 1
      NPROMA_WAM = HUGE(NPROMA_WAM)/2

!*    1. DEFINE UNIT NAMES.
!        ------------------

      IU05 = 5
      IU06 = 6

      IU07 = IWAM_GET_UNIT(IU06, 'wam_grid_tables', 'r', 'u', 0)

      IU12 = 12
      IU14 = 14
      IU15 = 15

!     1.1 INITIALISATION OF MPP VARIABLES WITH DEFAULT VALUES
!         ---------------------------------------------------

      LFDBIOOUT = .FALSE. 
      IRANK=1
      NPROC=1
      LL1D=.FALSE.
      X4(:)=1.0_JWRB
      NPRECR = KIND(X4)
      I4(:)=1
      NPRECI = KIND(I4)
      KTAG=100

      LLNEWCURR=.TRUE.

! ----------------------------------------------------------------------

! ALLOCATE NECESSARY ARRAYS

      ALLOCATE (NSTART(NPROC),NEND(NPROC))
      ALLOCATE (NBLKS(NPROC),NBLKE(NPROC))

! ----------------------------------------------------------------------

!*    2. READ NAMELIST NALINE.
!        ---------------------

      CDATEWO = ' '
      CDAWIFL = ' '
      CDATEFL = ' '

      READ (IU05, NALINE)

      IF (CLTUNIT == 'H') IDELWI = IDELWI*3600
      CDATEF = CDATEA
      CDATER = '000000000000'
      CDATES = '000000000000'
      ICPLEN=LEN_TRIM(CPATH)

      NTOTENS = 0
      NENSFNB = 0  
      ISTREAM =1045 !! if changed to an ifs stream also change LNEWLVTP
      NLOCGRB   = 1
      NSYSNB  = -1
      NMETNB  = -1
      IREFDATE= 0

      LGRHDIFS=.FALSE.
      LNEWLVTP=.FALSE.

      LWCUR=.FALSE.

! ----------------------------------------------------------------------

!*    3. READ PREPROC OUTPUT.
!        --------------------

      CALL READPRE (IU07)

      NINF=1
      NSUP=NIBLO

      IF (LLUNSTR) THEN

        IJS = 1
        IJL = NIBLO
        NXFF=NIBLO
        NYFF=1
        NSTART=1
        NEND=IJL
        IJSLOC=1
        IJLLOC=IJL
        IJGLOBAL_OFFSET=0
        NBLKS=NSTART
        NBLKE=NEND
        IF (ALLOCATED(IFROMIJ)) DEALLOCATE(IFROMIJ)
        ALLOCATE(IFROMIJ(IJSLOC:IJLLOC))
        IF (ALLOCATED(JFROMIJ)) DEALLOCATE(JFROMIJ)
        ALLOCATE(JFROMIJ(IJSLOC:IJLLOC))
        DO IJ = IJSLOC, IJLLOC 
          IFROMIJ(IJ)=IJ
          JFROMIJ(IJ)=1
        ENDDO
      ELSE
        NXFF=NGX
        NYFF=NGY
        NSTART=1
        NEND=IJL
        IJSLOC=1
        IJLLOC=IJL
        IJGLOBAL_OFFSET=0
        NBLKS=NSTART
        NBLKE=NEND
        IF (ALLOCATED(IFROMIJ)) DEALLOCATE(IFROMIJ)
        ALLOCATE(IFROMIJ(IJSLOC:IJLLOC))
        IF (ALLOCATED(JFROMIJ)) DEALLOCATE(JFROMIJ)
        ALLOCATE(JFROMIJ(IJSLOC:IJLLOC))
        DO IJ = IJSLOC, IJLLOC 
          IFROMIJ(IJ)=IXLG(IJ)
          JFROMIJ(IJ)=NGY-KXLT(IJ)+1
        ENDDO
      ENDIF



!!!   deallocate big arrays that were read in with READPRE

      IF (ALLOCATED(FAC0)) DEALLOCATE(FAC0)
      IF (ALLOCATED(FAC1)) DEALLOCATE(FAC1)
      IF (ALLOCATED(FAC2)) DEALLOCATE(FAC2)
      IF (ALLOCATED(FAC3)) DEALLOCATE(FAC3)
      IF (ALLOCATED(FAK)) DEALLOCATE(FAK)
      IF (ALLOCATED(FRHF)) DEALLOCATE(FRHF)
      IF (ALLOCATED(DFIMHF)) DEALLOCATE(DFIMHF)

      IF (ALLOCATED(OMEGA)) DEALLOCATE(OMEGA)
      IF (ALLOCATED(THH))   DEALLOCATE(THH)
      IF (ALLOCATED(DFDTH)) DEALLOCATE(DFDTH)
      IF (ALLOCATED(IM_P)) DEALLOCATE(IM_P)
      IF (ALLOCATED(IM_M)) DEALLOCATE(IM_M) 
      IF (ALLOCATED(TA)) DEALLOCATE(TA)
      IF (ALLOCATED(TB)) DEALLOCATE(TB)
      IF (ALLOCATED(TC_QL)) DEALLOCATE(TC_QL)
      IF (ALLOCATED(TT_4M)) DEALLOCATE(TT_4M)
      IF (ALLOCATED(TT_4P)) DEALLOCATE(TT_4P)
      IF (ALLOCATED(TFAKH)) DEALLOCATE(TFAKH)


!*    3.* SET GRIB HEADERS FOR INPUTS/OUTPUTS
!         -----------------------------------
      IF (.NOT. LGRHDIFS) THEN
!       FOR INTEGRATED PARAMETERS
        CALL PRESET_WGRIB_TEMPLATE("I",NGRIB_HANDLE_WAM_I)
!       FOR SPECTRA 
        CALL PRESET_WGRIB_TEMPLATE("S",NGRIB_HANDLE_WAM_S)
      ENDIF

! ----------------------------------------------------------------------

!*    4. PRINTER PROTOCOL OF INPUT.
!        --------------------------

      WRITE (IU06,'(A70)') HEADER
      WRITE (IU06,'('' MODEL OPTIONS  :'',/)')

      WRITE (IU06,*)'  '
      WRITE (IU06,*)' ******************************************'
      IF (LGRIBOUT) THEN
        WRITE (IU06,*)' THE OUTPUT SPECTRA WILL BE GRIBBED '
        WRITE (IU06,*)' INTO ',NFRE*NANG,' FIELDS'
      ELSE
        WRITE (IU06,*)' THE OUTPUT OF SPECTRA WILL BE BINARY  '
      ENDIF
      WRITE (IU06,*)' ******************************************'
      WRITE (IU06,*)'  '

      IF (IOPTI == 0) THEN
        WRITE (IU06,'('' INITIAL VALUES ARE COMPUTED FROM'',            &
     &   '' INPUT PARAMETERS.'')')
      ELSEIF (IOPTI == 1) THEN
        WRITE (IU06,'('' INITIAL VALUES ARE COMPUTED FROM'',            &
     &   '' LOCAL WIND.'')')
        WRITE (IU06,'('' WAVE ENERGY IS ZERO IN CALM WIND AREAS.'')')
      ELSEIF (IOPTI == 2) THEN
        WRITE (IU06,'('' INITIAL VALUES ARE COMPUTED FROM'',            &
     &   '' LOCAL WIND.'')')
        WRITE (IU06,'('' PARAMETERS USED IN CALM WIND AREAS.'')')
      ELSEIF (IOPTI == 3) THEN
        WRITE (IU06,'('' INITIAL VALUES ARE COMPUTED FROM'',            &
     &   '' IMPOSED SWELL SYSTEMS.'')')
      ELSE
        WRITE (IU06,'('' INVALID INPUT OPTION. PROGRAM WILL ABORT '')')
        CALL ABORT1
      ENDIF

      WRITE(IU06,*) '  '

      WRITE (IU06,*) ' TEST OUTPUT LEVEL IS .......... ITEST = ', ITEST
      WRITE (IU06,*) ' TEST OUTPUT IN BLOCK LOOP UPTO ITESTB = ', ITESTB

      WRITE (IU06,'('' JONSWAP PARAMETERS  :'',/)')
      WRITE (IU06,'('' ALFA : '',F10.5,'' FM : '',F10.5,'' GAMMA : '',  &
     &              F10.5,'' SA : '',F10.5,'' SB : '',F10.5)')          &
     &              ALFA, FM, GAMMA, SA, SB
      WRITE (IU06,'('' MEAN WAVE DIRECTION :  THETA = '',F10.5,         &
     &              '' DEGREE'')')  THETA
      WRITE (IU06,*) '  '
      WRITE (IU06,*) ' WIND INPUT TIMESTEP (SECONDS)      : ',IDELWI
      WRITE (IU06,*) '  '
      WRITE (IU06,*) ' END OF USER INPUT PROTOCOLL'
      WRITE (IU06,'(''  NUMBER OF DIRECTION BINS  NANG = '',I4)') NANG
      WRITE (IU06,'(''  NUMBER OF FREQUENCY BINS  NFRE = '',I4)') NFRE

! ----------------------------------------------------------------------

!*    5. PREPARE WINDFIELD.
!        ------------------

      CDATEWL = ZERO
      CDTPRO  = CDATEA
      CDATEE  = CDATEA
      CDAWIFL = CDATEA
      IDELPRO = IDELWI
      IDELWO  = IDELWI

      LWCOU=.FALSE.
      IREFRA=0
      ILEN=1

      IREAD=1

      IF (.NOT.ALLOCATED(FL1)) ALLOCATE (FL1(NIBLO,NANG,NFRE))
      IF (.NOT.ALLOCATED(U10OLD)) ALLOCATE(U10OLD(NIBLO))
      IF (.NOT.ALLOCATED(THWOLD)) ALLOCATE(THWOLD(NIBLO))
      IF (.NOT.ALLOCATED(USOLD)) ALLOCATE(USOLD(NIBLO))
      IF (.NOT.ALLOCATED(Z0OLD)) ALLOCATE(Z0OLD(NIBLO))
      IF (.NOT.ALLOCATED(TAUW)) ALLOCATE(TAUW(NIBLO))
      IF (.NOT.ALLOCATED(TAUWDIR)) ALLOCATE(TAUWDIR(NIBLO))
      IF (.NOT.ALLOCATED(ROAIRO)) ALLOCATE(ROAIRO(NIBLO))
      IF (.NOT.ALLOCATED(WSTAROLD)) ALLOCATE(WSTAROLD(NIBLO))
      IF (.NOT.ALLOCATED(CICOVER)) ALLOCATE(CICOVER(NIBLO))
      IF (.NOT.ALLOCATED(CITHICK)) ALLOCATE(CITHICK(NIBLO))
      IF (.NOT.ALLOCATED(FF_NEXT)) ALLOCATE(FF_NEXT(NIBLO))

      U10OLD(:) = WSPMIN
      THWOLD(:) = 0.0_JWRB
      USOLD(:) =  U10OLD(:)*0.035847_JWRB
      TAUW(:) = 0.1_JWRB*USOLD(:)
      TAUWDIR(:) = THWOLD(:)
      Z0OLD(:) = 0.00001_JWRB
      ROAIRO(:) = ROAIR      
      WSTAROLD(:) = 0.0_JWRB
      CICOVER(:) = 0.0_JWRB
      CITHICK(:) = 0.0_JWRB


      IF (IOPTI > 0 .AND. IOPTI /= 3) THEN

!!!! might need to restict call when needed !!!
!!! remove that call in 40R3
      CALL CIGETDEAC

        LLINIT=.FALSE.
        LLALLOC_FIELDG_ONLY=.FALSE.

        CALL PREWIND (U10OLD,THWOLD,USOLD,Z0OLD,            &
     &                ROAIRO, WSTAROLD,                     &
     &                CICOVER, CITHICK,                     &
     &                FF_NEXT,                              &
     &                LLINIT, LLALLOC_FIELDG_ONLY,          &
     &                IREAD,                                &
     &                NFIELDS, NGPTOTG, NC, NR,             &
     &                FIELDS, LWCUR, MASK_IN)

      TAUW(:) = 0.1_JWRB * USOLD(:)**2

      ENDIF

! ----------------------------------------------------------------------

!*    6. DEFINE FETCH AND MAXIMUM PEAK FREQUENCY.
!        ----------------------------------------


      IF (FETCH < 0.1E-5_JWRB) FETCH = 0.5_JWRB*DELPHI
      FRMAX = FM
      IF (IOPTI /= 0 .AND. IOPTI /= 3) THEN
        WRITE (IU06,*) ' FETCH USED (METRES)       : ', FETCH
        WRITE (IU06,*) ' MAXIMUM PEAK FREQUENCY IS : ', FRMAX
      ENDIF

! ----------------------------------------------------------------------

!*    7. GENERATE AND WRITE START FILES.
!        -------------------------------

      IF (IOPTI /= 3) THEN
        THETAQ = THETA * RAD
        CALL MSTART (IU12, IU14, IU15, IOPTI, FETCH, FRMAX,             &
     &              IJS, IJL, FL1, U10OLD, THWOLD)
      ELSE
        CALL MSWELL (IJS, IJL, FL1)

        IF (LLUNSTR) THEN
!         reset points with no flux out of the boundary to 0
          DO M=1,NFRE
            DO K=1,NANG
              DO IJ=1,NIBLO
                FL1(IJ,K,M) = FL1(IJ,K,M)*IOBPD(K,IJ)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ENDIF

! ----------------------------------------------------------------------

!*    8. DISPOSE START FILES.
!        --------------------

!     STRESS RELATED FIELDS :
      WRITE(IU06,*) 'MINVAL OF U10OLD = ',MINVAL(U10OLD)
      WRITE(IU06,*) 'MAXVAL OF U10OLD = ',MAXVAL(U10OLD)
      IF (.NOT.LGRIBOUT) THEN
        CALL SAVSTRESS(U10OLD, THWOLD, USOLD, TAUW, TAUWDIR, Z0OLD,     &
     &                 ROAIRO, WSTAROLD, CICOVER, CITHICK,              &
     &                 NBLKS, NBLKE, CDATEA, CDATEA)
      ENDIF

      IF (ALLOCATED(U10OLD)) DEALLOCATE(U10OLD)
      IF (ALLOCATED(THWOLD)) DEALLOCATE(THWOLD)
      IF (ALLOCATED(USOLD)) DEALLOCATE(USOLD)
      IF (ALLOCATED(Z0OLD)) DEALLOCATE(Z0OLD)
      IF (ALLOCATED(TAUW)) DEALLOCATE(TAUW)
      IF (ALLOCATED(TAUWDIR)) DEALLOCATE(TAUWDIR)
      IF (ALLOCATED(ROAIRO)) DEALLOCATE(ROAIRO)
      IF (ALLOCATED(WSTAROLD)) DEALLOCATE(WSTAROLD)
      IF (ALLOCATED(CICOVER)) DEALLOCATE(CICOVER)
      IF (ALLOCATED(CITHICK)) DEALLOCATE(CITHICK)

      WRITE (IU06,*) ' SAVING SPECTRA '

      CALL FLUSH(IU06)
!     SPECTRA :
      IF (LGRIBOUT) THEN
!       THE COLD START SPECTRA WILL BE SAVED AS GRIB FILES.
        CDTPRO  = CDATEA
        CALL OUTSPEC(IJS, IJL, FL1, CICOVER)
      ELSE
        CALL SAVSPEC(IJS, IJL, FL1, NBLKS, NBLKE, CDATEA, CDATEA, CDUM)
      ENDIF

! ----------------------------------------------------------------------

!*    9. END OF JOB: DELETE WORK FILES.
!        ------------------------------

      WRITE (IU06,*) ' '
      WRITE (IU06,*) ' PROGRAM PRESET: ALL DONE'

      CALL MPL_END()

      IF (LHOOK) CALL DR_HOOK('PRESET',1,ZHOOK_HANDLE)

      END
