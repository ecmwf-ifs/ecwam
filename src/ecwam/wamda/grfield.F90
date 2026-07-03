SUBROUTINE GRFIELD(BLK2GLO, FL1, INTFLDS, MINIJS, MAXIJL, WHMOD,   &
 &                 CICVR, SIGMOD, DISTMAX)

!--------------------------------------------------------------------

!**** *GRFIELD* - READS GRIDDED MEASUREMENTS AND DISTRIBUTES THEM.
!                 ALSO COMPUTES THE ALTIMETER DATA CORRECTIONS.

!     B.HANSEN       ECMWF       JANUARY 1998
!     P.JANSSEN      ECMWF       FEBRUARY 1999 ADD SEA STATE DEPENDENT
!                                              CORRECTION TO ALTIMETER
!                                              WAVE HEIGHT.
!     S. ABDALLA     ECMWF       NOVEMBER 2011 ADD ODB;  CONVERT TO F90

!     PURPOSE.
!     --------

!       PREPARE THE DATA FOR OPTIMAL INTERPOLATION.

!**   INTERFACE.
!     ----------

!       *CALL* *GRFIELD(BLK2GLO,FL1,INTFLDS, MINIJS,MAXIJL,WHMOD,CICVR,SIGMOD,DISTMAX)*

!          *BLK2GLO*  BLOCK TO GRID TRANSFORMATION
!          *FL1*      2D SPECTRUM
!          *INTFLDS*  INTEGRATED/DERIVED PARAMETERS
!          *MINIJS*   MINUMUM INDEX OF WHMOD AND CICVR 
!          *MAXIJL*   MAXIMUM INDEX OF WHMOD AND CICVR
!          *WHMOD*    MODEL WAVE HEIGHT (FIRST GUESS) (global block indexing)
!          *CICVR*    MODEL SEA ICE (global block indexing)
!          *SIGMOD*   MODEL ERROR. 
!          *DISTMAX*  MAXIMUM SPREADING DISTANCE IN RADIAN !!


!     METHOD.
!     -------

!        READ PRE PREPARED SATELLITE DATA FROM FILE OR ODB AND DISTRIBUTE TO
!        OTHER PES.

!      EXTERNALS.
!      ----------

!         SKEWNESS       - COMPUTES PARAMETERS OF THE NEARLY-GAUSSIAN
!                          DISTRIBUTION OF OCEAN WAVES

!      MODIFICATONS.
!      -------------
!         NONE.

!-----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO, INTGT_PARAM_FIELDS

      USE YOWABORT , ONLY : WAM_ABORT
      USE YOWALTAS , ONLY : NUMALT   ,NIJALT   ,NALTDT   ,NALTAVLB ,    &
     &            IJALT    ,INTLMAX  ,KMINLMAX ,KMAXLMAX ,NOBSPE   ,    &
     &            ALTSDTHRSH, ALTBGTHRSH, ALTGRTHRSH, HSALTCUT,         &
     &            ALTDATA ,IBUFRSAT ,NALTUDT  ,ALTUNDATA,               &
     &            XKAPPA2  ,HSCOEFCOR,HSCONSCOR ,LALTCOR ,LALTLRGR ,    &
     &            LALTGRDOUT, LODBRALT, CSATNAME,LALTPASSIV,            &
     &            LRALTPREPROC
      USE YOWGRID  , ONLY : IJSLOC   ,IJLLOC   ,IJGLOBAL_OFFSET    ,    &
     &                      SINPH    ,COSPH,                            &
     &                      NPROMA_WAM, NCHNK, ICHNKFROMIJ, IPRMFROMIJ 
      USE YOWICE   , ONLY : CITHRSH_SAT  ,FLMIN
      USE YOWMAP   , ONLY : XDELLA   ,ZDELLO   ,NIBLO
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,KTAG
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,LL1D     ,LLUNSTR
      USE YOWPCONS , ONLY : RAD      ,ZMISS
      USE YOWSTAT  , ONLY : CDTPRO   ,YEXPVER  ,IDELALT  ,ISTREAM
      USE YOWSPEC  , ONLY : NSTART   ,NEND     ,IJ2NEWIJ
      USE YOWTEST  , ONLY : IU06
      USE YOWASSI  , ONLY : GETODBRALT

#ifdef WAM_HAVE_UNWAM
      USE YOWPD,     ONLY : NODES=>NODES_GLOBAL, RANK
#endif
      USE YOWSPHERE, ONLY : SPHERICAL_COORDINATE_DISTANCE

      USE EC_LUN   , ONLY : NULERR
      USE MPL_MODULE

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

#include "abort1.intfb.h"
#include "confile.intfb.h"

#include "grdata.intfb.h"
#include "iwam_get_unit.intfb.h"
#include "mpbcastintfld.intfb.h"
#include "secondhh.intfb.h"
#include "skewness.intfb.h"
#include "wstream_strg.intfb.h"

      TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(IN) :: FL1
      TYPE(INTGT_PARAM_FIELDS), INTENT(INOUT) :: INTFLDS
      INTEGER(KIND=JWIM), INTENT(IN) :: MINIJS, MAXIJL
      REAL(KIND=JWRB), INTENT(IN) :: DISTMAX, SIGMOD
      REAL(KIND=JWRB), DIMENSION(MINIJS:MAXIJL), INTENT(IN) :: WHMOD, CICVR

      INTEGER(KIND=JWIM) :: IOS,IREAD, IDUM, NOBS, NOBSPEOLD
      INTEGER(KIND=JWIM) :: ISAT
      INTEGER(KIND=JWIM) :: IOBS, JOBS, ITAG, IFIELD, IFAIL, IUME, NOBSDUM, IOBSPE
      INTEGER(KIND=JWIM) :: IOBS_S, KFROM, IOBS_EXTRA
      INTEGER(KIND=JWIM) :: NIJALTP1, IIJALT
      INTEGER(KIND=JWIM) :: NFPAS, NBLKL, NREJP, NREJICE, NREJSD, NREJBG, NSAT
      INTEGER(KIND=JWIM) :: NTOTFPAS, NTOTBLKL, NTOTREJP, NTOTREJICE, NTOTREJSD, NTOTREJBG, NTOTSAT
      INTEGER(KIND=JWIM) :: IJ, IJG, M, K, IR, JR, IALTDT, IJJ, KI, KJ
      INTEGER(KIND=JWIM) :: IP, IK, ISTART
      INTEGER(KIND=JWIM) :: KRCOUNT, KRTAG, MPLENGTH
      INTEGER(KIND=JWIM) :: IZCOMLEN, ICOUNT
      INTEGER(KIND=JWIM) :: IUMAIL
      INTEGER(KIND=JWIM) :: NERS
      INTEGER(KIND=JWIM) :: IZCOMLENMAX, NDISPE, IDISPE, IPR, NOBSIN
      INTEGER(KIND=JWIM) :: NOBSPEMAX
      INTEGER(KIND=JWIM) :: NFRPEALT, NOBSPELOCMAX
      INTEGER(KIND=JWIM) :: ISEND(1)
      INTEGER(KIND=JWIM) :: NOBSPSAT(NUMALT)
      INTEGER(KIND=JWIM) :: NSATOBS(NUMALT)
      INTEGER(KIND=JWIM),DIMENSION(7) :: ISUMBUF
      INTEGER(KIND=JWIM),DIMENSION(NPROC) :: ISENDLEN
      INTEGER(KIND=JWIM),ALLOCATABLE :: ISENDREQ(:)
      INTEGER(KIND=JWIM),ALLOCATABLE,DIMENSION(:,:) :: IJALTPE, ITEMP
      INTEGER(KIND=JWIM),ALLOCATABLE,DIMENSION(:) :: IJALT_LOC

      REAL(KIND=JWRB) :: COSDISTMAX
      REAL(KIND=JWRB) :: DELLON, COSLON, DIS2, XLONGKJ
      REAL(KIND=JWRB) :: THRSHLD, THRSHLDGR, HSCUT, REJECT, REJRATIO
      REAL(KIND=JWRB) :: RCOR, HS, XN, HS_NEW
      REAL(KIND=JWRB) :: WGT
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      REAL(KIND=JWRB),ALLOCATABLE :: XLONG(:)
      REAL(KIND=JWRB),ALLOCATABLE :: XKAPPA1(:), DELH_ALT(:)
      REAL(KIND=JWRB),ALLOCATABLE :: FALT(:,:,:)
      REAL(KIND=JWRB),ALLOCATABLE :: ZCOMBUFROD(:)
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:) :: WHGTAGB, CWHTAGB, RANGAGB, AVGWGT
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:,:) :: ALTDATAPE, ZTEMP

      REAL(KIND=JWRU)  :: XLONJD, XLATJD, DISTD, DISTMAXD

      REAL(KIND=JWRU),ALLOCATABLE,DIMENSION(:) :: XLONID, XLATID
      REAL(KIND=JWRU),ALLOCATABLE,DIMENSION(:,:) :: ALTUNDATAPE, DTEMP
      REAL(KIND=JWRU),ALLOCATABLE :: ZCOMBUFS(:,:)
      REAL(KIND=JWRU),ALLOCATABLE :: ZCOMBUFR(:)

      CHARACTER(LEN=2) :: CDUM
      CHARACTER(LEN=3) :: CLFID
      CHARACTER(LEN=4) :: CSTRM
      CHARACTER(LEN=14) :: CTIME

      LOGICAL :: LLOPENED, LLTEST
      LOGICAL :: LAVERAGE
      LOGICAL :: LALTCORRECTION
      LOGICAL :: LLNOTFOUND
      LOGICAL :: LDUM
      LOGICAL, SAVE :: LLFRSTSCOR
      LOGICAL :: LLCOMMON(NPROC)
      LOGICAL :: LLALT(NUMALT)
      LOGICAL,ALLOCATABLE, DIMENSION(:) :: LLWITHIN
      LOGICAL,ALLOCATABLE, DIMENSION(:,:) :: LLIN

      DATA LLFRSTSCOR / .TRUE. /


!*     VARIABLE     TYPE      PURPOSE.
!      --------     ----      --------

!     *IREAD*       INTGER    PROCESSOR WHICH WILL READ THE INPUT FILES.
!     *CLFID        CHARACTER THREE CHARACTER FILE IDENTIFIER.

      DATA CLFID, IREAD /'RFL', 1/


!     -------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GRFIELD',0,ZHOOK_HANDLE)

!the size of the 2nd dimension of IJALT is increased by 1
!to keep track of the input order (see below).
NIJALTP1=NIJALT+1

!*    1. GET THE ALTIMETER DATA
!        ----------------------

IDUM=0
NOBS=0
NSATOBS(:)=0
NOBSPSAT(:)=0

IF (LODBRALT) THEN
!* INPUT DATA FROM ODB (rfl4wam needs to be run as part of preprocessing of the data
  CALL GETODBRALT(IREAD, NIJALTP1, NOBS, NSATOBS)

ELSE
!* IF LODBRALT is false, read observations from file (either preprocessed or raw data)

  IF (LRALTPREPROC) THEN
!*  READ RALT OBSERVATIONS FROM A FILE AS PREPROCESSED BY RFL4WAM
    WRITE(IU06,*) '     GRFIELD :  ALTIMETER DATA ARE READ FROM A FILE AS PREPRARED BY RFL4WAM'
    WRITE(IU06,*) ''
    IF (IRANK == IREAD) THEN
        IFIELD = 0
        IFAIL = 1

!*       FIND SPARE UNIT TO ATTACH INPUT FILE.

        UNIT: DO IUME=33,99,1
          INQUIRE ( UNIT=IUME, OPENED=LLOPENED)
          IF ( .NOT. LLOPENED) THEN
            EXIT UNIT
          ENDIF
        ENDDO UNIT

!*       IF FOUND: SEEK INPUT FILE, READ HEADER AND FIELDS.
!*                 IFIELD MAY BE MISSING INDICATING THERE ARE NO FIELDS.

        IF (IUME < 100) THEN
          LLTEST=.FALSE.
          CALL CONFILE (IU06, IUME, CDTPRO, CLFID, IFAIL, LLTEST)
          IF (IFAIL == 0) THEN
            READ (IUME,ERR=9000,IOSTAT=IOS) CTIME, NOBS
            IF (NOBS > 0) THEN
              ALLOCATE(IJALT(NOBS,NIJALTP1))
              ALLOCATE(ALTDATA(NOBS,NALTDT))
              IF (LLUNSTR) ALLOCATE(ALTUNDATA(NOBS,NALTUDT))

              READ (IUME,ERR=9000,IOSTAT=IOS) ((IJALT(IOBS,IIJALT),IOBS = 1,NOBS),IIJALT = 1,NIJALT)
              READ (IUME,ERR=9000,IOSTAT=IOS) ALTDATA
              IF (LLUNSTR) READ (IUME,ERR=9000,IOSTAT=IOS) ALTUNDATA
            ENDIF
            CLOSE (IUME)
          ELSE
            WRITE(IU06,'("WARNING IN GRFIELD:  W A R N I N G  !!!!")')
            WRITE(IU06,'("NO FILES FOUND"                          )')
            WRITE(IU06,'("WARNING IN GRFIELD:  W A R N I N G  !!!!")')
            CALL FLUSH(IU06)
            NOBS=0
          ENDIF
        ELSE
!*        IF NOT FOUND: ISSUE WARNING MESSAGE.
          WRITE(IU06,'("WARNING IN GRFIELD:  W A R N I N G  !!!!!!")')
          WRITE(IU06,'("NO SPARE UNIT TO READ SATELLITE DATA FOUND")')
          WRITE(IU06,'("WARNING IN GRFIELD:  W A R N I N G  !!!!!!")')
          CALL FLUSH(IU06)
          NOBS=0
        ENDIF
    ENDIF

    GOTO 9001
9000 WRITE(IU06,'("READING ERROR IN GRFIELD: ")')
     WRITE(IU06,'(25X,"I/O STATUS IS",I8)') IOS
     WRITE(IU06,'(25X," ABORT !!!!!!!!")')
     CALL ABORT1
9001 CONTINUE

  ELSE
    WRITE(IU06,*) '     GRFIELD :  ALTIMETER DATA ARE READ FROM A FILE AS PREPRARED BY URAPRE'
    WRITE(IU06,*) ''
!!  Currently hardcoded 
    LAVERAGE=.TRUE.
!!
    IF (IRANK == IREAD) THEN
      CALL GRDATA (NOBS, IDELALT, LAVERAGE, NIJALTP1)
      WRITE(IU06,*) ''
    ENDIF
  ENDIF

ENDIF   ! end of if LODBRALT block

IF (IRANK == IREAD) THEN
  IF (NOBS == 0) THEN
    WRITE(IU06,*) '     GRFIELD: NO RADAR ALTIMETER (RALT) AVAILABLE'
  ELSE
    WRITE(IU06,*) '     GRFIELD: RADAR ALTIMETER, TOTAL NUMBER OF ENTRIES :',  NOBS
    DO ISAT=1,NUMALT
      IF (NSATOBS(ISAT) > 0 ) THEN
!       note that NSATOBS is currently only meaningful if LODBRALT is true
        WRITE(IU06,*) '  GRFIELD : NUMBER OF ENTRIES FOR ALTIMETER ',IBUFRSAT(ISAT),' IS ', NSATOBS(ISAT)
      ENDIF
    ENDDO
  ENDIF
ENDIF


IF (NOBS == 0) THEN
  DO IJ = IJSLOC, IJLLOC 
    IK = ICHNKFROMIJ(IJ)
    IP = IPRMFROMIJ(IJ)
    INTFLDS%ALTWH(IP, IK) = ZMISS
    INTFLDS%CALTWH(IP, IK) = ZMISS
    INTFLDS%RALTCOR(IP, IK) = ZMISS
  ENDDO
ENDIF

IF (IRANK == IREAD) THEN

      IF (NOBS > 0) THEN
!       REORGANISE THE IJ INDEX WHEN 2-D DECOMPOSITION
!!      For unstructured grid this conversion for local to global can only be done
!!      once  IJALT(IOBS,1) has been distributed (see below)
        IF (.NOT.LL1D .AND. NPROC > 1 .AND. LRALTPREPROC .AND. .NOT.LLUNSTR) THEN
          DO IOBS = 1, NOBS
            IJALT(IOBS,1) = IJ2NEWIJ(IJALT(IOBS,1))
          ENDDO
        ENDIF

!       TAG THE DATA SO THAT THE DATA CAN BE REORDERED
!       (to insure that the results are the same on any number of PE's)
        IF (LODBRALT) THEN
          ITAG=0
          DO ISAT=1,NUMALT
            IF (NSATOBS(ISAT) > 0 ) THEN
              DO IJ = 1, NIBLO
                DO IOBS = 1, NOBS
                  IF (IJALT(IOBS,2) == IBUFRSAT(ISAT)) THEN
                    IF (IJALT(IOBS,1) == IJ) THEN
                      ITAG=ITAG+1
                      IJALT(IOBS,NIJALTP1) = ITAG
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ELSE
          DO IOBS = 1, NOBS
            IJALT(IOBS,NIJALTP1) = IOBS
          ENDDO
        ENDIF


        DO ISAT=1,NUMALT
          DO IOBS = 1, NOBS
!       take only points that belong to the satellite
            IF (IJALT(IOBS,2) == IBUFRSAT(ISAT)) NOBSPSAT(ISAT)=NOBSPSAT(ISAT)+1
          ENDDO
        ENDDO
!       check whether sum(NOBSPSAT)=NOBS
        NOBSDUM=0
        DO ISAT=1,NUMALT
          NOBSDUM=NOBSDUM+NOBSPSAT(ISAT)
        ENDDO
        IF (NOBSDUM /= NOBS) THEN
          WRITE(IU06,*) ' '
          WRITE(IU06,'("IN GRFIELD:  W A R N I N G  !!!!!!")')
          WRITE(IU06,*) 'THE TOTAL NUMBER OF ALT DATA DOES NOT '
          WRITE(IU06,*) 'CORRESPOND TO THE SUM OF THE DIFFERENT '
          WRITE(IU06,*) 'CONTRIBUTION '
          WRITE(IU06,*) 'NOBS = ',NOBS
          DO ISAT=1,NUMALT
            WRITE(IU06,*) 'NOBSPSAT = ',NOBSPSAT(ISAT),ISAT
          ENDDO
          WRITE(IU06,*) ' '
          WRITE(IU06,*) 'NO DATA ASSIMILATION WILL TAKE PLACE !!!'
          WRITE(IU06,*) ' '
          WRITE(IU06,'("IN GRFIELD:  W A R N I N G  !!!!!!")')
          WRITE(IU06,*) ' '
          DO ISAT=1,NUMALT
            NOBSPSAT(ISAT)=0
          ENDDO
          CALL FLUSH(IU06)
        ENDIF
      ENDIF
ENDIF


!     SEND NOBS TO THE OTHER PE'S TO LET THEM KNOW WHETHER THERE ARE
!     DATA OR NOT.
      CALL MPBCASTINTFLD(IREAD,KTAG,NUMALT,1,NOBSPSAT)
      KTAG=KTAG+1
      NOBS=0
      DO ISAT=1,NUMALT
        NOBS=NOBS+NOBSPSAT(ISAT)
      ENDDO

!     IF NO DATA, NOTHING TO BE DONE
      IF (NOBS <= 0) THEN
        NOBSPE(:) = 0
!     -------------------------------------------------------------------------

       IF (LHOOK) CALL DR_HOOK('GRFIELD',1,ZHOOK_HANDLE)

!     -------------------------------------------------------------------------
        RETURN
      ENDIF

!     Continue processing the data

!     add a small number to DISTMAX as a safety precaution
!     i.e each PE will get slightly more data that they should
!     strictly have.
      COSDISTMAX=COS(DISTMAX+0.1_JWRB)
      DISTMAXD=REAL(DISTMAX+0.1_JWRB,JWRU)

! ----------------------------------------------------------------------

!*    2. SEND FIELDS TO ALL OTHER PES.
!        -----------------------------

      IF (NPROC > 1) THEN
!       GET THE NUMBER OF DATA PER PE.
        IF (IRANK == IREAD) THEN
          NOBSPE(:) = 0
          DO IR = 1, NPROC
            IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
              DO IOBS = 1, NOBS
                DO IP=1,RANK(IR)%NP
                  IF ( IJALT(IOBS,1) == RANK(IR)%IPLG(IP) ) THEN
                    NOBSPE(IR) = NOBSPE(IR) + 1
                    EXIT
                  ENDIF
                ENDDO
              ENDDO
#else
          CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
            ELSE
              DO IOBS = 1, NOBS
                IF ( IJALT(IOBS,1) >= NSTART(IR) .AND. &
     &               IJALT(IOBS,1) <= NEND(IR) ) THEN
                  NOBSPE(IR) = NOBSPE(IR) + 1
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDIF
        CALL MPBCASTINTFLD(IREAD,KTAG,NPROC,1,NOBSPE)
        KTAG=KTAG+1
        NOBSPEMAX=MAXVAL(NOBSPE(:))

!       DETERMINE THE NUMBER OF PE'S EACH PE WILL NEED TO GET DATA FROM
        NFRPEALT = 0
        NOBSPELOCMAX = 0
        DO JR = 1, IRANK-1
          IF (INTLMAX(JR) == 1 .AND. NOBSPE(JR) > 0) THEN
            NFRPEALT = NFRPEALT + 1
            NOBSPELOCMAX=MAX(NOBSPELOCMAX,NOBSPE(JR))
          ENDIF
        ENDDO
        DO JR = IRANK+1,NPROC
          IF (INTLMAX(JR) == 1 .AND. NOBSPE(JR) > 0) THEN
            NFRPEALT = NFRPEALT + 1
            NOBSPELOCMAX=MAX(NOBSPELOCMAX,NOBSPE(JR))
          ENDIF
        ENDDO


!       DISTRIBUTE THE DATA FOR EACH PE.
        IF (IRANK == IREAD) THEN

!         EXTRACT THE DATA THAT ARE ON EACH PE.
          ALLOCATE(LLIN(NOBS,NPROC))

          NDISPE=0
          DO IR = 1, NPROC
            IF (NOBSPE(IR) > 0 ) THEN
              NDISPE=NDISPE+1
              IOBSPE=0
              IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
                DO IOBS = 1, NOBS
                  LLIN(IOBS,IR)=.FALSE.
                  DO IP=1,RANK(IR)%NP
                    IF ( IJALT(IOBS,1) == RANK(IR)%IPLG(IP) ) THEN
                      LLIN(IOBS,IR)=.TRUE.
                      IOBSPE=IOBSPE+1
                      EXIT
                    ENDIF
                  ENDDO
                ENDDO
#else
          CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
              ELSE
                DO IOBS = 1, NOBS
                  IF ( IJALT(IOBS,1) >= NSTART(IR) .AND. IJALT(IOBS,1) <= NEND(IR) ) THEN
                    LLIN(IOBS,IR)=.TRUE.
                    IOBSPE=IOBSPE+1
                  ELSE
                    LLIN(IOBS,IR)=.FALSE.
                  ENDIF
                ENDDO
              ENDIF

              IF (IOBSPE /= NOBSPE(IR)) THEN
                WRITE(IU06,*)'ERROR IN GRFIELD: IOBSPE.NE.NOBSPE(IR) !'
                WRITE(IU06,*)'IOBSPE =', IOBSPE
                WRITE(IU06,*)'NOBSPE(IR) =',NOBSPE(IR)
                CALL FLUSH(IU06)
                WRITE(NULERR,*)'ERROR IN GRFIELD: IOBSPE.NE.NOBSPE(IR) !'
                WRITE(NULERR,*)'IOBSPE =', IOBSPE
                WRITE(NULERR,*)'NOBSPE(IR) =',NOBSPE(IR)
                CALL ABORT1
              ENDIF
            ENDIF
          ENDDO

          IF (LLUNSTR) THEN
            IZCOMLENMAX = (NIJALTP1+NALTDT+NALTUDT) * NOBSPEMAX
          ELSE
            IZCOMLENMAX = (NIJALTP1+NALTDT) * NOBSPEMAX
          ENDIF
          ALLOCATE(ZCOMBUFS(IZCOMLENMAX,NDISPE))

          IDISPE=0
          DO IR = 1, NPROC
            IF (NOBSPE(IR) > 0 ) THEN
              ICOUNT = 0
              IDISPE=IDISPE+1
              DO IIJALT = 1, NIJALTP1
                DO IOBS = 1, NOBS
                  IF (LLIN(IOBS,IR)) THEN
                    ICOUNT = ICOUNT + 1
                    ZCOMBUFS(ICOUNT,IDISPE) = REAL(IJALT(IOBS,IIJALT),JWRU)
                  ENDIF
                ENDDO
              ENDDO
              DO IALTDT = 1, NALTDT
                DO IOBS = 1, NOBS
                  IF (LLIN(IOBS,IR)) THEN
                    ICOUNT = ICOUNT + 1
                    ZCOMBUFS(ICOUNT,IDISPE) = REAL(ALTDATA(IOBS,IALTDT),JWRU)
                  ENDIF
                ENDDO
              ENDDO

              IF (LLUNSTR) THEN
                DO IALTDT = 1, NALTUDT
                  DO IOBS = 1, NOBS
                    IF (LLIN(IOBS,IR)) THEN
                      ICOUNT = ICOUNT + 1
                      ZCOMBUFS(ICOUNT,IDISPE) = ALTUNDATA(IOBS,IALTDT)
                    ENDIF
                  ENDDO
                ENDDO
              ENDIF

              ISENDLEN(IR)=ICOUNT
            ELSE
              ISENDLEN(IR)=0
            ENDIF
          ENDDO

          ALLOCATE(ISENDREQ(NDISPE))
          IDISPE=0
          DO IR = 1, NPROC
            IF (NOBSPE(IR) > 0 ) THEN
              IDISPE=IDISPE+1
              ICOUNT = ISENDLEN(IR)
              CALL GSTATS(618,0)
              CALL MPL_SEND(ZCOMBUFS(1:ICOUNT,IDISPE),KDEST=IR,KTAG=KTAG, &
     &                      KMP_TYPE=JP_NON_BLOCKING_STANDARD, &
     &                      KREQUEST=ISENDREQ(IDISPE), &
     &                      CDSTRING='GRFIELD 1:')
              CALL GSTATS(618,1)
            ENDIF
          ENDDO

          DEALLOCATE(LLIN)
          DEALLOCATE(IJALT)
          DEALLOCATE(ALTDATA)
          IF (LLUNSTR) DEALLOCATE(ALTUNDATA)
        ENDIF

        IF (LLUNSTR) THEN
          IZCOMLEN = (NIJALTP1+NALTDT+NALTUDT) * NOBSPEMAX
        ELSE
          IZCOMLEN = (NIJALTP1+NALTDT) * NOBSPEMAX
        ENDIF
        ALLOCATE(ZCOMBUFR(IZCOMLEN))

        IF (NOBSPE(IRANK) > 0 ) THEN
          CALL GSTATS(618,0)
!!CHECK!! CALL MPL_RECV(ZCOMBUFR(1:IZCOMLEN),KSOURCE=IREAD,KTAG=KTAG, &
!!CHECK!!                                    *******
!!OLD!!!! CALL MPL_RECV(ZCOMBUFR(1:IZCOMLEN),KFROM  =IREAD,KTAG=KTAG, &

          CALL MPL_RECV(ZCOMBUFR(1:IZCOMLEN),KSOURCE=IREAD,KTAG=KTAG, &
     &       KOUNT=KRCOUNT,KRECVTAG=KRTAG,CDSTRING='GRFIELD 1:')
          IF (KRTAG /= KTAG) CALL MPL_ABORT &
     &      ('MPL_RECV ERROR in GRFIELD 1 : MISMATCHED TAGS' )
          CALL GSTATS(618,1)

!         KEEP THE DATA ON EACH PE.
          IF (ALLOCATED(IJALT)) DEALLOCATE(IJALT)
          ALLOCATE(IJALT(NOBSPE(IRANK),NIJALTP1))
          IF (ALLOCATED(ALTDATA)) DEALLOCATE(ALTDATA)
          ALLOCATE(ALTDATA(NOBSPE(IRANK),NALTDT))
          IF (LLUNSTR) THEN
            IF (ALLOCATED(ALTUNDATA)) DEALLOCATE(ALTUNDATA)
            ALLOCATE(ALTUNDATA(NOBSPE(IRANK),NALTUDT))
          ENDIF
          ICOUNT = 0
          DO IIJALT = 1, NIJALTP1
            DO IOBS = 1, NOBSPE(IRANK)
              ICOUNT = ICOUNT + 1
              IJALT(IOBS,IIJALT) = NINT(ZCOMBUFR(ICOUNT))
            ENDDO
          ENDDO
          DO IALTDT = 1, NALTDT
            DO IOBS = 1, NOBSPE(IRANK)
              ICOUNT = ICOUNT + 1
              ALTDATA(IOBS,IALTDT) = REAL(ZCOMBUFR(ICOUNT),JWRB)
            ENDDO
          ENDDO
          IF (LLUNSTR) THEN
            DO IALTDT = 1, NALTUDT
              DO IOBS = 1, NOBSPE(IRANK)
                ICOUNT = ICOUNT + 1
                ALTUNDATA(IOBS,IALTDT) = ZCOMBUFR(ICOUNT)
              ENDDO
            ENDDO
          ENDIF

        ENDIF

! Wait any outstanding sends to complete

        CALL GSTATS(618,0)
        IF (IRANK == IREAD) CALL MPL_WAIT(KREQUEST=ISENDREQ,CDSTRING='GRFIELD 1')
        KTAG = KTAG + 1
        CALL GSTATS(618,1)
        IF (ALLOCATED(ISENDREQ)) DEALLOCATE(ISENDREQ)
        IF (ALLOCATED(ZCOMBUFS)) DEALLOCATE(ZCOMBUFS)
        IF (ALLOCATED(ZCOMBUFR)) DEALLOCATE(ZCOMBUFR)


      ELSE
        NOBSPE = NOBS
      ENDIF

!     FIND OUT WHICH SATELLITE IS PRESENT
      NALTAVLB=0
      DO ISAT=1,NUMALT
        IF (NOBSPSAT(ISAT) > 0) THEN
          LLALT(ISAT)=.TRUE.
          NALTAVLB=NALTAVLB+1
        ELSE
          LLALT(ISAT)=.FALSE.
        ENDIF
        WRITE(IU06,'(5X,I2,A,A25,A,I7,A,L7)') &
     &        ISAT,'. SATELLITE: ',CSATNAME(ISAT),&
     &     '  ... NOBSPSAT=',NOBSPSAT(ISAT), &
     &     '  ... LLALT=',LLALT(ISAT)
      ENDDO

! ----------------------------------------------------------------------

!*    3. APPLY CORRECTION TO WAVE HEIGHT DATA ACCORDING TO
!        SUPPLIED LINEAR REGRESSION COEFFICIENTS (obsolete !)
!        --------------------------------------------------

      DO ISAT=1,NUMALT
        IF (LALTLRGR(ISAT) .AND. LLALT(ISAT)) THEN
          DO IOBS = 1, NOBSPE(IRANK)
            IF (ALTDATA(IOBS,1) > 0._JWRB .AND. IJALT(IOBS,3) == 1 .AND. &
     &          IJALT(IOBS,2) == IBUFRSAT(ISAT) ) THEN
              ALTDATA(IOBS,1) = HSCOEFCOR(ISAT)*ALTDATA(IOBS,1)+HSCONSCOR(ISAT)
            ENDIF
          ENDDO
        ENDIF
      ENDDO

! ----------------------------------------------------------------------

!     4. APPLY SCREENING OF DATA BASED ON MODEL FIELD
!        --------------------------------------------

      ALLOCATE(IJALT_LOC(MAX(NOBSPE(IRANK),1)))
      DO IOBS = 1, NOBSPE(IRANK)
        IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
!         if unstructured grid IJALT(IOBS,1) is the global index, but IJ is not
          LLNOTFOUND=.TRUE.
          ISTART=RANK(IRANK)%ISTART
          DO IP=1,RANK(IRANK)%NP
            IF (IJALT(IOBS,1) == RANK(IRANK)%IPLG(IP)) THEN
              IJALT_LOC(IOBS)=IP
              LLNOTFOUND=.FALSE.
              EXIT
            ENDIF
          ENDDO
          IF (LLNOTFOUND) THEN
!           it should not happen !!!
            WRITE(IU06,*)'GRFIELD: GLOBAl TO LOCAL INDEXING PROBLEM !'
            WRITE(IU06,*)'IRANK =', IRANK
            CALL FLUSH(IU06)
            WRITE(NULERR,*)'GRFIELD: GLOBAl TO LOCAL INDEXING PROBLEM !'
            WRITE(NULERR,*)'IRANK =', IRANK
            CALL ABORT1
          ENDIF
#else
          CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
        ELSE
          IJALT_LOC(IOBS)=IJALT(IOBS,1)
        ENDIF
      ENDDO

      ALLOCATE(LLIN(MAX(NOBSPE(IRANK),1),IRANK:IRANK))
      DO IOBS = 1, NOBSPE(IRANK)
        IF ( IJALT_LOC(IOBS) >= IJSLOC .AND. IJALT_LOC(IOBS) <= IJLLOC ) THEN
          LLIN(IOBS,IRANK)=.TRUE.
        ELSE
          LLIN(IOBS,IRANK)=.FALSE.
        ENDIF
      ENDDO

      DO ISAT=1,NUMALT

        WRITE(IU06,*)' '
        WRITE(IU06,*)'   GRFIELD:'
        WRITE(IU06,*)'   ALTIMETER ',CSATNAME(ISAT)

        NFPAS=0
        NBLKL=0
        NREJP=0
        NREJICE=0
        NREJSD=0
        NREJBG=0
        NSAT=0

        IF (ALTSDTHRSH(ISAT) < 0._JWRB) THEN
          IF (XDELLA < 0.3_JWRB) THEN
            ALTSDTHRSH(ISAT) = 0.75_JWRB
          ELSE
            ALTSDTHRSH(ISAT) = 0.6_JWRB
          ENDIF
          WRITE(IU06,'(A,F6.2,A)') &
     &        '    SUSPICIOUS DATA THRESHOLD NOT GIVEN. IT IS SET TO ', &
     &         ALTSDTHRSH(ISAT)*100.,'%.'
        ENDIF

        IF (LLALT(ISAT)) THEN

          DO IOBS = 1, NOBSPE(IRANK)

!           take only points that belong to the satellite
            IF (IJALT(IOBS,2) == IBUFRSAT(ISAT)) THEN

              IF (LALTPASSIV(IJALT(IOBS,2))) NFPAS=NFPAS+1

              IF (LLIN(IOBS,IRANK)) THEN 
                IF (IJALT(IOBS,3) == 0) NBLKL=NBLKL+1
                IF (IJALT(IOBS,3) < 0) NREJP=NREJP+1
              ENDIF

!             EXCLUDE POINTS WHICH ARE ON MODEL SEA ICE

              IJG=IJALT(IOBS,1)
              IF (CICVR(IJG) > CITHRSH_SAT .AND. IJALT(IOBS,3) >= 0) THEN
                IJALT(IOBS,3)=-5
                IF (LLIN(IOBS,IRANK)) NREJICE=NREJICE+1
              ENDIF

!             LOOK FOR SUSPICIOUS DATA (if too many then the whole dataset will be discarded)
!             ON VALID DATA ONLY.

              IF (ALTDATA(IOBS,1) > 0._JWRB .AND. IJALT(IOBS,3) >= 0 .AND. WHMOD(IJG) > 0._JWRB) THEN
                THRSHLD=ABS(ALTDATA(IOBS,1)-WHMOD(IJG))/WHMOD(IJG)
                IF (THRSHLD >= ALTSDTHRSH(ISAT)) THEN
                  IF (LLIN(IOBS,IRANK)) NREJSD=NREJSD+1
                ENDIF
              ENDIF

!              BACKGROUND CHECK (remove data that deviate too much from the first guess)
!              ON VALID DATA ONLY.

              IF (ALTDATA(IOBS,1) > 0._JWRB .AND. IJALT(IOBS,3) >=0 .AND. WHMOD(IJG) > 0._JWRB) THEN
                THRSHLD=(ALTDATA(IOBS,2)/SIGMOD)* &
     &                  (ALTDATA(IOBS,1)-WHMOD(IJG))**2/(ALTDATA(IOBS,1)*WHMOD(IJG))
                THRSHLDGR=ABS(ALTDATA(IOBS,1)-WHMOD(IJG))
                HSCUT=MIN(HSALTCUT(ISAT),ALTDATA(IOBS,2))
                IF (THRSHLD >= ALTBGTHRSH(ISAT) .OR. THRSHLDGR >= ALTGRTHRSH(ISAT) ) THEN
                  IJALT(IOBS,3)=-6
                  IF (LLIN(IOBS,IRANK)) NREJBG=NREJBG+1
                ELSE IF (ALTDATA(IOBS,1) <= HSCUT .OR. WHMOD(IJG) <= HSCUT ) THEN 
                  IJALT(IOBS,3)=-6
                  IF (LLIN(IOBS,IRANK)) NREJBG=NREJBG+1
                ENDIF
              ENDIF
            ENDIF
          ENDDO
          IF (LLALT(ISAT)) THEN
            NSAT=1
          ELSE
            NSAT=0
          ENDIF
        ENDIF

        ISUMBUF(1)=NFPAS
        ISUMBUF(2)=NBLKL
        ISUMBUF(3)=NREJP
        ISUMBUF(4)=NREJICE
        ISUMBUF(5)=NREJSD
        ISUMBUF(6)=NREJBG
        ISUMBUF(7)=NSAT

        IF (NPROC > 1) THEN
          CALL GSTATS(618,0)
          CALL MPL_ALLREDUCE(ISUMBUF,'SUM',CDSTRING='GRFIELD:')
          CALL GSTATS(618,1)
        ENDIF

        NTOTFPAS=ISUMBUF(1)
        NTOTBLKL=ISUMBUF(2)
        NTOTREJP=ISUMBUF(3)
        NTOTREJICE=ISUMBUF(4)
        NTOTREJSD=ISUMBUF(5)
        NTOTREJBG=ISUMBUF(6)
        NTOTSAT=ISUMBUF(7)

        WRITE(IU06,*)'   TOTAL NUMBER OF GRIDDED WAVE HEIGHTS= ', NOBSPSAT(ISAT)
        IF (NOBSPSAT(ISAT) > 0) THEN
          WRITE(IU06,*)'   NUMBER OF SUBAREA WITH DATA= ',NTOTSAT
          WRITE(IU06,*)''
          WRITE(IU06,*)'   NUMBER PASSIVE = ', NTOTFPAS
          WRITE(IU06,*)'   NUMBER BLACKLISTED = ', NTOTBLKL
          WRITE(IU06,*)'   NUMBER REJECTED BY PRE-PROCESSING = ', NTOTREJP
          WRITE(IU06,*)'   NUMBER OVER MODEL SEA ICE= ',NTOTREJICE
          WRITE(IU06,*)'   NUMBER REJECTED BY THE BACKGROUND CHECK= ',NTOTREJBG
          WRITE(IU06,*)'   NUMBER OF SUSPICIOUS DATA= ',NTOTREJSD
        ENDIF
        WRITE(IU06,*)''
        WRITE(IU06,*)'-------------------------------------------------------'
        CALL FLUSH(IU06)

        IF (NTOTSAT > 0) THEN

          IF (XDELLA < 0.3_JWRB) THEN
            REJECT=0.75_JWRB
          ELSE
            REJECT=0.5_JWRB
          ENDIF

          REJRATIO=FLOAT(NTOTREJSD)/FLOAT(NOBS)

          IF (REJRATIO > REJECT) THEN
            WRITE(IU06,*)'*********************************************'
            WRITE(IU06,*)' '
            WRITE(IU06,*)'   ALTIMETER ',IBUFRSAT(ISAT)
            WRITE(IU06,*)'  THE NUMBER OF SUSPICIOUS ALT DATA IS OVER ', &
     &                   REJECT*100,' %'
            WRITE(IU06,*)' '
            WRITE(IU06,*)'  THEY ARE ALL REJECTED AS A PRECAUTION !'
            WRITE(IU06,*)' '
            WRITE(IU06,*)'*********************************************'
            DO IOBS = 1, NOBSPE(IRANK)
              IF (IJALT(IOBS,2) == IBUFRSAT(ISAT)) IJALT(IOBS,3)=-6
            ENDDO
            IF (IRANK == 1) THEN
              IUMAIL=IWAM_GET_UNIT(IU06,'MM','S','F',0,'READWRITE')
              WRITE(IUMAIL,*)'          !!!!!! WARNING !!!!!! '
              WRITE(IUMAIL,*)' NUMBER OF SUSPICIOUS DATA IS OVER', &
     &                       REJECT*100,' % ',REJRATIO*100,' %'
              WRITE(IUMAIL,*)' IN EXPERIMENT ',YEXPVER
              CALL WSTREAM_STRG(ISTREAM,CSTRM,IDUM,IDUM,CDUM,IDUM,LDUM)
              WRITE(IUMAIL,*)' STREAM        ',CSTRM
              WRITE(IUMAIL,*)' AT TIME       ',CDTPRO
              WRITE(IUMAIL,*)'          !!!!!! WARNING !!!!!! '
              CLOSE(IUMAIL)
              CALL SYSTEM (' mail wab < MM ' )
            ENDIF
          ENDIF

        ENDIF

      ENDDO ! END LOOP ON SATELLITE
      DEALLOCATE(LLIN)
      WRITE(IU06,*)' '
      CALL FLUSH(IU06)


! ----------------------------------------------------------------------

!*    5. APPLY CORRECTION TO WAVE HEIGHT DATA ACCORDING TO SROKOSZ.
!        ---------------------------------------------------------
!        SEE ECMWF TECH MEMO 269 or JANSSEN 2000: ECMWF wave modelling
!        and satellite altimeter wave data, Satellites, Oceanography and
!        Society, ed. D Halpern, Elsevier: Amsterdam. pp 35-56.

!     5.1 COLLOCATE WAMSPECTRA WITH ALT-HS OBSERVATIONS (if needed).
!         ---------------------------------------------

      LALTCORRECTION=.FALSE.
      DO ISAT=1,NUMALT
        IF (LLALT(ISAT).AND.LALTCOR(ISAT)) LALTCORRECTION=.TRUE.
      ENDDO


!     The second order correction has not been used since the days of ERS, so it is a bit obsolete....
      IF (LALTCORRECTION .AND. LLFRSTSCOR) THEN
        LLFRSTSCOR=.FALSE.
        CALL SECONDHH
      ENDIF


      ALLOCATE(WHGTAGB(IJSLOC:IJLLOC))
      WHGTAGB(:)=0._JWRB
      ALLOCATE(CWHTAGB(IJSLOC:IJLLOC))
      CWHTAGB(:)=0._JWRB
      ALLOCATE(RANGAGB(IJSLOC:IJLLOC))
      RANGAGB(:)=0._JWRB
      ALLOCATE(AVGWGT(IJSLOC:IJLLOC))
      AVGWGT(:)=0._JWRB


      IF (NOBSPE(IRANK) > 0) THEN
        NERS=NOBSPE(IRANK)
        ALLOCATE(XKAPPA1(NERS))
        ALLOCATE(DELH_ALT(NERS))

        IF (LALTCORRECTION) THEN
          ALLOCATE(FALT(NERS,NANG,NFRE))
          DO IOBS = 1, NOBSPE(IRANK)
            IF (ALTDATA(IOBS,1) > 0._JWRB .AND. IJALT(IOBS,3) == 1) THEN
              IJ=IJALT_LOC(IOBS)
              IK = ICHNKFROMIJ(IJ)
              IP = IPRMFROMIJ(IJ)
              DO M=1,NFRE
                DO K=1,NANG
                  FALT(IOBS,K,M) = FL1(IP, K, M, IK)
                ENDDO
              ENDDO
            ELSE
              DO M=1,NFRE
                DO K=1,NANG
                 FALT(IOBS,K,M)=FLMIN
                ENDDO
             ENDDO
            ENDIF
          ENDDO

          CALL SKEWNESS(IU06,FALT,NERS,XKAPPA1,DELH_ALT)
        ELSE
          XKAPPA1(:)=1._JWRB
          DELH_ALT(:)=0._JWRB
        ENDIF


!       5.3 COMPUTE GRIB BOX AVERAGE VALUES FOR GRIB OUPUT AND APPLY SROKOSZ CORRECTION (if LALTCOR(ISAT)).
!           -----------------------------------------------------------------------------------------------
     
        DO ISAT=1,NUMALT
          IF (LLALT(ISAT)) THEN

            DO IOBS = 1, NOBSPE(IRANK)
!             TAKE ONLY VALID OBSERVATIONS
              IF (ALTDATA(IOBS,1) > 0._JWRB .AND. IJALT(IOBS,3) == 1  .AND. &
     &            IJALT(IOBS,2) == IBUFRSAT(ISAT)) THEN
                RCOR = DELH_ALT(IOBS)
                HS = ALTDATA(IOBS,1)
                IF (LALTCOR(ISAT)) THEN
                  XN = HS**2/16.0_JWRB + XKAPPA2(ISAT)
                  HS_NEW = 16.0_JWRB*(XKAPPA1(IOBS)*XN - XKAPPA2(ISAT))
                  IF ( HS_NEW >  0.0_JWRB ) THEN
                    HS_NEW = SQRT(HS_NEW)
                    ALTDATA(IOBS,1) = HS_NEW
                  ELSE
                    HS_NEW = 0.0_JWRB
                  ENDIF
                ELSE
                  HS_NEW = HS
                ENDIF

!               SAVE ALTIMETER DATA TO GRIB
                IF (LALTGRDOUT(ISAT)) THEN
                  IJ=IJALT_LOC(IOBS)
                  IF (IJ >= IJSLOC .AND. IJ <= IJLLOC) THEN
!                   WHGTAGB CONTAINS THE ORIGINAL DATA
!                   A WEIGHTED AVERAGE IS COMPUTED IF MORE THAN ONE SATELLITE
                    IF (HS_NEW > 0._JWRB) THEN
                      WGT=1._JWRB/ALTDATA(IOBS,2)**2
                      WHGTAGB(IJ) = WHGTAGB(IJ) + WGT*ALTDATA(IOBS,3)
                      CWHTAGB(IJ) = CWHTAGB(IJ) + WGT*HS_NEW
                      RANGAGB(IJ) = RANGAGB(IJ) + WGT*HS_NEW*RCOR
                      AVGWGT(IJ) = AVGWGT(IJ) + WGT
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDDO

        DO IJ = IJSLOC, IJLLOC
          IF (WHGTAGB(IJ) > 0._JWRB) THEN
            WHGTAGB(IJ) = WHGTAGB(IJ)/AVGWGT(IJ)
            CWHTAGB(IJ) = CWHTAGB(IJ)/AVGWGT(IJ)
            RANGAGB(IJ) = RANGAGB(IJ)/AVGWGT(IJ)
            RANGAGB(IJ) = RANGAGB(IJ)/CWHTAGB(IJ)
          ELSE
            WHGTAGB(IJ) = ZMISS
            CWHTAGB(IJ) = ZMISS
            RANGAGB(IJ) = ZMISS
          ENDIF
        ENDDO

        DEALLOCATE(XKAPPA1)
        DEALLOCATE(DELH_ALT)
        IF (LALTCORRECTION) DEALLOCATE(FALT)
      ELSE
        DO IJ = IJSLOC, IJLLOC
          WHGTAGB(IJ) = ZMISS
          CWHTAGB(IJ) = ZMISS
          RANGAGB(IJ) = ZMISS
        ENDDO
      ENDIF ! ON NOBSPE


!     COLLECT THE OUTPUT OF THE ALTIMETER CORRECTIONS
      DO IJ = IJSLOC, IJLLOC 
        IK = ICHNKFROMIJ(IJ)
        IP = IPRMFROMIJ(IJ)
        INTFLDS%ALTWH(IP, IK) = WHGTAGB(IJ)
        INTFLDS%CALTWH(IP, IK) = CWHTAGB(IJ)
        INTFLDS%RALTCOR(IP, IK) = RANGAGB(IJ)
      ENDDO

      DEALLOCATE(WHGTAGB)
      DEALLOCATE(CWHTAGB)
      DEALLOCATE(RANGAGB)
      DEALLOCATE(AVGWGT)

      DEALLOCATE(IJALT_LOC)

!     GET ALL RELEVANT DATA ON EACH PE.
!     BY DISTRIBUTING THEM ONTO THE DIFFERENT PE'S

      IF (NPROC > 1) THEN
        IF (LLUNSTR) THEN
          IZCOMLENMAX = (NIJALTP1+NALTDT+NALTUDT) * NOBSPE(IRANK)
        ELSE
          IZCOMLENMAX=(NIJALTP1+NALTDT) * NOBSPE(IRANK)
        ENDIF
        NDISPE=0
        DO IR = 1, NPROC
          IF (INTLMAX(IR) == 1 .AND. IR /= IRANK) THEN
            IF (NOBSPE(IRANK) > 0 ) THEN
              LLCOMMON(IR)=.TRUE.
              NDISPE=NDISPE+1
            ELSE
              LLCOMMON(IR)=.FALSE.
            ENDIF
          ELSE
            LLCOMMON(IR)=.FALSE.
          ENDIF
        ENDDO

        IF (IZCOMLENMAX > 0) THEN
!         PACKING THE CONTRIBUTION WHICH MIGHT BE COMMON
          ALLOCATE(ZCOMBUFS(IZCOMLENMAX,NDISPE))
          ALLOCATE(ISENDREQ(NDISPE))
          IDISPE=0
          DO IR = 1, NPROC
            IF (LLCOMMON(IR)) THEN
              IDISPE=IDISPE+1
              ICOUNT = 0
              DO IIJALT = 1, NIJALTP1
                DO IOBS = 1, NOBSPE(IRANK)
                  ICOUNT = ICOUNT + 1
                  ZCOMBUFS(ICOUNT,IDISPE) = REAL(IJALT(IOBS,IIJALT),JWRU)
                ENDDO
              ENDDO
              DO IALTDT = 1, NALTDT
                DO IOBS = 1, NOBSPE(IRANK)
                  ICOUNT = ICOUNT + 1
                  ZCOMBUFS(ICOUNT,IDISPE) = REAL(ALTDATA(IOBS,IALTDT),JWRU)
                ENDDO
              ENDDO
              IF (LLUNSTR) THEN
                DO IALTDT = 1, NALTUDT
                  DO IOBS = 1, NOBSPE(IRANK)
                    ICOUNT = ICOUNT + 1
                    ZCOMBUFS(ICOUNT,IDISPE) = ALTUNDATA(IOBS,IALTDT)
                  ENDDO
                ENDDO
              ENDIF
              ISENDLEN(IR)=ICOUNT
            ELSE
              ISENDLEN(IR)=0
            ENDIF
          ENDDO

!         SENDING THE CONTRIBUTION WHICH MIGHT BE COMMON
!         SEND NON BLOCKING THE BUFFERS

          IDISPE=0
          CALL GSTATS(618,0)
          DO IR = 1, NPROC
            IF (LLCOMMON(IR)) THEN
              IDISPE=IDISPE+1
              ICOUNT = ISENDLEN(IR)
              CALL MPL_SEND(ZCOMBUFS(1:ICOUNT,IDISPE), &
     &                      KDEST=IR,KTAG=KTAG, &
     &                      KMP_TYPE=JP_NON_BLOCKING_STANDARD, &
     &                      KREQUEST=ISENDREQ(IDISPE), &
     &                      CDSTRING='GRFIELD 3:')
            ENDIF
          ENDDO
          CALL GSTATS(618,1)
        ENDIF

        IF (LLUNSTR) THEN
          ALLOCATE(XLONID(IJSLOC:IJLLOC))
          ALLOCATE(XLATID(IJSLOC:IJLLOC))
#ifdef WAM_HAVE_UNWAM
!!! should be replaced by permanent array
          DO IP=1,RANK(IRANK)%NP
            XLONID(IP) = NODES(RANK(IRANK)%IPLG(IP))%X
            XLATID(IP) = NODES(RANK(IRANK)%IPLG(IP))%Y
          ENDDO
#else
          CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
        ELSE
          ALLOCATE(XLONG(IJSLOC:IJLLOC))
          DO IJ=IJSLOC,IJLLOC
            XLONG(IJ) = REAL(BLK2GLO%IXLG(IJ)-1,JWRB)*ZDELLO(BLK2GLO%KXLT(IJ))
          ENDDO
        ENDIF

!       RECEIVING THE CONTRIBUTIONS WHICH MIGHT BE COMMON

        DO IPR = 1, NFRPEALT

          IF (LLUNSTR) THEN
            IZCOMLEN = (NIJALTP1+NALTDT+NALTUDT) * NOBSPELOCMAX
          ELSE
            IZCOMLEN = (NIJALTP1+NALTDT) * NOBSPELOCMAX
          ENDIF
          ALLOCATE(ZCOMBUFR(IZCOMLEN))

          CALL GSTATS(618,0)
          CALL MPL_RECV(ZCOMBUFR(1:IZCOMLEN),KFROM=IR,KTAG=KTAG, &
     &           KOUNT=KRCOUNT,KRECVTAG=KRTAG,CDSTRING='GRFIELD 4:')
          CALL GSTATS(618,1)
          IF (KRTAG /= KTAG) CALL MPL_ABORT &
     &    ('MPL_RECV ERROR in GRFIELD 3 : MISMATCHED TAGS' )


          ALLOCATE(IJALTPE(NOBSPE(IR),NIJALTP1))
          ICOUNT = 0
          DO IIJALT = 1, NIJALTP1
            DO IOBS = 1, NOBSPE(IR)
              ICOUNT = ICOUNT + 1
              IJALTPE(IOBS,IIJALT) = NINT(ZCOMBUFR(ICOUNT))
            ENDDO
          ENDDO

          ALLOCATE(ALTDATAPE(NOBSPE(IR),NALTDT))
          DO IALTDT = 1, NALTDT
            DO IOBS = 1, NOBSPE(IR)
              ICOUNT = ICOUNT + 1
              ALTDATAPE(IOBS,IALTDT) = REAL(ZCOMBUFR(ICOUNT),JWRB)
            ENDDO
          ENDDO

          IF (LLUNSTR) THEN
            ALLOCATE(ALTUNDATAPE(NOBSPE(IR),NALTUDT))
            DO IALTDT = 1, NALTUDT
              DO IOBS = 1, NOBSPE(IR)
                ICOUNT = ICOUNT + 1
                ALTUNDATAPE(IOBS,IALTDT) = ZCOMBUFR(ICOUNT)
              ENDDO
            ENDDO
          ENDIF

          DEALLOCATE(ZCOMBUFR)

!         how many observations are actually relevant to PE IRANK
!         Keep only the relevant information
          NOBSIN=0
          ALLOCATE(LLWITHIN(NOBSPE(IR)))

          IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
            DO IOBS = 1, NOBSPE(IR)
              XLATJD = ALTUNDATAPE(IOBS,1)
              XLONJD = ALTUNDATAPE(IOBS,2)

!!! limit of lat separation ???
              LLWITHIN(IOBS) = .FALSE.

                DO IJ=IJSLOC,IJLLOC
                  CALL SPHERICAL_COORDINATE_DISTANCE(XLONID(IJ),XLONJD,XLATID(IJ),XLATJD,DISTD)
                  IF (DISTD <= DISTMAXD) THEN
                    NOBSIN=NOBSIN+1
                    LLWITHIN(IOBS) = .TRUE.
                    EXIT
                  ENDIF
                ENDDO

            ENDDO
#else
          CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif

          ELSE
            DO IOBS = 1, NOBSPE(IR)
              IJJ=IJALTPE(IOBS,1)
              KJ = BLK2GLO%KXLT(IJJ)
              XLONGKJ = REAL(BLK2GLO%IXLG(IJJ)-1,JWRB)*ZDELLO(KJ)
              LLWITHIN(IOBS) = .FALSE.
              IF (KJ >= KMINLMAX(IRANK) .AND.  KJ <= KMAXLMAX(IRANK) ) THEN
                 DO IJ=IJSLOC,IJLLOC
                     KI = BLK2GLO%KXLT(IJ)
                     DELLON = XLONGKJ-XLONG(IJ)
                     COSLON = COS(DELLON*RAD)
                     DIS2 = COSLON*COSPH(KJ)*COSPH(KI)+ &
     &                    SINPH(KJ)*SINPH(KI)
                     IF (DIS2 >= COSDISTMAX) THEN
                       NOBSIN=NOBSIN+1
                       LLWITHIN(IOBS) = .TRUE.
                       EXIT
                    ENDIF
                ENDDO
              ENDIF
            ENDDO
          ENDIF

          IF (NOBSIN > 0) THEN
            NOBSPEOLD = NOBSPE(IRANK)
            NOBSPE(IRANK) = NOBSPE(IRANK) + NOBSIN

            IF (NOBSPEOLD > 0) THEN
              ALLOCATE(ITEMP(NOBSPEOLD,NIJALTP1))
              ITEMP(:,:) = IJALT(:,:)
              DEALLOCATE(IJALT)
              ALLOCATE(IJALT(NOBSPE(IRANK),NIJALTP1))
              DO IIJALT = 1, NIJALTP1
                DO IOBS = 1, NOBSPEOLD
                  IJALT(IOBS,IIJALT) = ITEMP(IOBS,IIJALT)
                ENDDO
              ENDDO
              DEALLOCATE(ITEMP)

              ALLOCATE(ZTEMP(NOBSPEOLD,NALTDT))
              ZTEMP(:,:) = ALTDATA(:,:)
              DEALLOCATE(ALTDATA)
              ALLOCATE(ALTDATA(NOBSPE(IRANK),NALTDT))
              DO IALTDT = 1, NALTDT
                DO IOBS = 1, NOBSPEOLD
                  ALTDATA(IOBS,IALTDT) = ZTEMP(IOBS,IALTDT)
                ENDDO
              ENDDO
              DEALLOCATE(ZTEMP)

              IF (LLUNSTR) THEN
                ALLOCATE(DTEMP(NOBSPEOLD,NALTUDT))
                DTEMP(:,:) = ALTUNDATA(:,:)
                DEALLOCATE(ALTUNDATA)
                ALLOCATE(ALTUNDATA(NOBSPE(IRANK),NALTUDT))
                DO IALTDT = 1, NALTUDT
                  DO IOBS = 1, NOBSPEOLD
                    ALTUNDATA(IOBS,IALTDT) = DTEMP(IOBS,IALTDT)
                  ENDDO
                ENDDO
                DEALLOCATE(DTEMP)
              ENDIF
            ELSE
              IF (ALLOCATED(IJALT)) DEALLOCATE(IJALT)
              ALLOCATE(IJALT(NOBSPE(IRANK),NIJALTP1))
              IF (ALLOCATED(ALTDATA)) DEALLOCATE(ALTDATA)
              ALLOCATE(ALTDATA(NOBSPE(IRANK),NALTDT))
              IF (LLUNSTR) THEN
                IF (ALLOCATED(ALTUNDATA)) DEALLOCATE(ALTUNDATA)
                ALLOCATE(ALTUNDATA(NOBSPE(IRANK),NALTUDT))
              ENDIF
            ENDIF

            DO IIJALT = 1, NIJALTP1
              ICOUNT = 0
              DO IOBS = 1, NOBSPE(IR)
                IF (LLWITHIN(IOBS)) THEN
                  ICOUNT = ICOUNT + 1
                  IJALT(NOBSPEOLD+ICOUNT,IIJALT) = IJALTPE(IOBS,IIJALT)
                ENDIF
              ENDDO
            ENDDO

            DO IALTDT = 1, NALTDT
              ICOUNT = 0
              DO IOBS = 1, NOBSPE(IR)
                IF (LLWITHIN(IOBS)) THEN
                  ICOUNT = ICOUNT + 1
                  ALTDATA(NOBSPEOLD+ICOUNT,IALTDT) = ALTDATAPE(IOBS,IALTDT)
                ENDIF
              ENDDO
            ENDDO

            IF (LLUNSTR) THEN
              DO IALTDT = 1, NALTUDT
                ICOUNT = 0
                DO IOBS = 1, NOBSPE(IR)
                  IF (LLWITHIN(IOBS)) THEN
                    ICOUNT = ICOUNT + 1
                    ALTUNDATA(NOBSPEOLD+ICOUNT,IALTDT) = ALTUNDATAPE(IOBS,IALTDT)
                  ENDIF
                ENDDO
              ENDDO
            ENDIF

          ENDIF

          DEALLOCATE(LLWITHIN)
          DEALLOCATE(IJALTPE)
          DEALLOCATE(ALTDATAPE)
          IF (LLUNSTR) DEALLOCATE(ALTUNDATAPE)

        ENDDO

        IF (LLUNSTR) THEN
          DEALLOCATE(XLONID)
          DEALLOCATE(XLATID)
        ELSE
          DEALLOCATE(XLONG)
        ENDIF

! Wait any outstanding sends to complete

        CALL GSTATS(618,0)
        IF (IZCOMLENMAX > 0) CALL MPL_WAIT(KREQUEST=ISENDREQ,CDSTRING='GRFIELD')
        CALL GSTATS(618,1)

        IF (ALLOCATED(ISENDREQ)) DEALLOCATE(ISENDREQ)
        KTAG = KTAG + 1
        IF (ALLOCATED(ZCOMBUFS)) DEALLOCATE(ZCOMBUFS)


!       MAKE ALL PE's AWARE OF THE NEW NOBSPE
        ISENDLEN(:)=1
        CALL GSTATS(618,0)
        ISEND=NOBSPE(IRANK:IRANK)
        CALL MPL_ALLGATHERV(ISEND,NOBSPE,ISENDLEN,CDSTRING='GRFIELD:')
        CALL GSTATS(618,1)

        KTAG = KTAG + 1

      ENDIF

!     REORDER THE DATA (to insure that the results are the same
!                       on any number of PE's)
      IF (NOBSPE(IRANK) > 0) THEN
        ALLOCATE(ITEMP(NOBS,1))
        ITEMP(:,1)=0
        DO IOBS = 1, NOBSPE(IRANK)
            ITEMP(IJALT(IOBS,NIJALTP1),1)=IJALT(IOBS,NIJALTP1)
        ENDDO
        ICOUNT=0
        DO IOBS = 1, NOBS
          IF (ITEMP(IOBS,1) > 0) THEN
            ICOUNT=ICOUNT+1
            ITEMP(IOBS,1)=ICOUNT
          ENDIF
        ENDDO
        DO IOBS = 1, NOBSPE(IRANK)
          IJALT(IOBS,NIJALTP1)=ITEMP(IJALT(IOBS,NIJALTP1),1)
        ENDDO
        DEALLOCATE(ITEMP)

        ALLOCATE(ITEMP(NOBSPE(IRANK),NIJALT))
        DO IIJALT = 1, NIJALT
            DO IOBS = 1, NOBSPE(IRANK)
              ITEMP(IJALT(IOBS,NIJALTP1),IIJALT)=IJALT(IOBS,IIJALT)
            ENDDO
        ENDDO
        ALLOCATE(ZTEMP(NOBSPE(IRANK),NALTDT))
        DO IALTDT = 1, NALTDT
          DO IOBS = 1, NOBSPE(IRANK)
            ZTEMP(IJALT(IOBS,NIJALTP1),IALTDT) = ALTDATA(IOBS,IALTDT)
          ENDDO
        ENDDO

        IF (LLUNSTR) THEN
          ALLOCATE(DTEMP(NOBSPE(IRANK),NALTUDT))
          DO IALTDT = 1, NALTUDT
            DO IOBS = 1, NOBSPE(IRANK)
              DTEMP(IJALT(IOBS,NIJALTP1),IALTDT) = ALTUNDATA(IOBS,IALTDT)
            ENDDO
          ENDDO
        ENDIF

        DEALLOCATE(IJALT)
        ALLOCATE(IJALT(NOBSPE(IRANK),NIJALT))
        DO IIJALT = 1, NIJALT
          DO IOBS = 1, NOBSPE(IRANK)
            IJALT(IOBS,IIJALT) = ITEMP(IOBS,IIJALT)
          ENDDO
        ENDDO
        DEALLOCATE(ITEMP)

        DO IALTDT = 1, NALTDT
          DO IOBS = 1, NOBSPE(IRANK)
            ALTDATA(IOBS,IALTDT) = ZTEMP(IOBS,IALTDT)
          ENDDO
        ENDDO
        DEALLOCATE(ZTEMP)

        IF (LLUNSTR) THEN
          DO IALTDT = 1, NALTUDT
            DO IOBS = 1, NOBSPE(IRANK)
              ALTUNDATA(IOBS,IALTDT) = DTEMP(IOBS,IALTDT)
            ENDDO
          ENDDO
          DEALLOCATE(DTEMP)
        ENDIF

      ENDIF

IF (LHOOK) CALL DR_HOOK('GRFIELD',1,ZHOOK_HANDLE)

END SUBROUTINE GRFIELD
