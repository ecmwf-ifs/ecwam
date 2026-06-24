      SUBROUTINE GRDATA(NWH, IDATAWL, LAVERAGE)

!--------------------------------------------------------------------

!**** *GRDATA* - TRANSFERS MEASUREMENTS TO MODEL GLOBAL GRID,
!****            DISCARDING UNRELIABLE DATA.

!     P.LIONELLO     ECMWF       APRIL 1990

!     PURPOSE.
!     --------

!       PREPARE THE DATA FOR OPTIMAL INTERPOLATION.

!**   INTERFACE.
!     ----------

!       *CALL* *GRDATA(NWH, IDATAWL, LAVERAGE )*
!        *NWH*  NUMBER OF RELEVANT GRIDDED OBSERVATIONS.
!        *LAVERAGE* LOGICAL TRUE IF SUPER-OBSERVATIONS ARE TO BE 
!                                   ASSIMILATED

!     METHOD.
!     -------

!        A BOX OF SIZE DELLA*DELLO IS CENTERED IN EACH GRID POINT.
!        THE AVERAGE VALUE OF THE SATELLITE MEASUREMENTS INSIDE THE
!        BOX IS TAKEN TO UPDATE THE WAM MODEL. THE VALUE IS DISCARDED
!        IF THE ROOT MEAN SQUARE ERROR IS TOO BIG.
!        MOREOVER THE VALUE IS DISCARDED IF THERE ARE LESS THEN FOUR
!        MEASUREMENTS INSIDE THE BOX.

!      EXTERNALS.
!      ----------

!         READSAT  - READ A MEASUREMENT.

!      MODIFICATONS.
!      -------------
!         FEB 92  Bjoern Hansen :
!                 Interface to readsat extended in order to transfer
!                 the flags added to the satellite information by the
!                 quality control preura.
!                 The as unreliable flagged satellite measurements
!                 are rejected.
!                 The input unit of the measurements is automatically
!                 looked up each time grdata is called. It is assumed
!                 that the needed input file is present in the working
!                 directory, contains data of the same period of time
!                 as the output time step and is named RFLyyyymmddhhmm,
!                 where yyyymmddhhmm is a DATE/TIME group.

!          Aug 97 Jean Bidlot : in coupled mode IDELWO is not the
!                 appropriate time span for the selection of the data
!                 to be used but rather IDATAWL, A NEW VARIABLE WHICH
!                 IS SUPPLIED IN WAMINPUT

!          Oct 97 Jean Bidlot : use of the new ice mask instead of 
!                 SST field to determine whether an observation which
!                 falls on a model sea point.

!          May 01 Jean Bidlot : do not use wind speed observations

!                 Jean Bidlot : only use grid information to thin
!                               and smooth the data.

!          J. BIDLOT       *ECMWF*       SEP 2001
!            INTRODUCE THE CONCEPT OF A BLACKLIST TO REMOVE DATA FROM
!            BEING PROCESSED. IT IS TRIGGERED IF THE FILE
!            wave_data_blacklist IS PRESENT.
!            THE BLACKLIST WILL CONTAIN THE FOLLOWING INFORMATION
!            STATID : THE SATELLITE IDENTIFICATION NUMBER (AS IN BUFR)
!            SENSOR : THE INSTRUMENT IDENTIFICATION NUMBER (AS IN BUFR)
!            THE DATA WILL BE BLACKLISTED BETWEEN THE DATE CDATESTART
!            AND CDATEEND AND BETWEEN LATMI<= LATITUDE <= LATMA
!            AND BETWEEN LONGMI <= LONGITUDE <= LONGMA
!            FOR the OBSERVED PARAMETER PRESCRIBED BY PARAM (AS IN BUFR)
!            PROVIDED IT IS BETWEEN THE VALUES VALUE_MIN AND VALUE_MAX

!-----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWABORT , ONLY : WAM_ABORT
      USE YOWALTAS , ONLY : NUMALT   ,NIJALT   ,NALTDT   ,IJALT    ,    &
     &                      ALTDATA  ,IBUFRSAT ,NALTEDT  ,NALTUDT  ,    &
     &                      ALTEXDATA,ALTUNDATA,CDATEOBS
      USE YOWGRID  , ONLY : NCHNK  ,KIJL4CHNK, IJFROMCHNK
      USE YOWMAP   , ONLY : BLK2GLO  ,IPER     ,AMOWEP   ,              &
     &                      AMOSOP   ,XDELLA   ,ZDELLO   ,NLONRGG  ,    &
     &                      NGX      ,NGY     ,NIBLO
      USE YOWMPP   , ONLY : NPROC
      USE YOWPARAM , ONLY : LLUNSTR
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWTEST  , ONLY : IU06
#ifdef WAM_HAVE_UNWAM
      USE OUTPUT_STRUCT, ONLY : IPELEMENT_CLOSEST
#endif
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "difdate.intfb.h"
#include "incdate.intfb.h"
#include "iwam_get_unit.intfb.h"
#include "readsat.intfb.h"

      INTEGER(KIND=JWIM), INTENT(OUT) :: NWH
      INTEGER(KIND=JWIM), INTENT(IN) :: IDATAWL
      LOGICAL, INTENT(IN) :: LAVERAGE

      INTEGER(KIND=JWIM), PARAMETER :: MAXSATBC=25
      INTEGER(KIND=JWIM) :: NNMIN
      INTEGER(KIND=JWIM) :: IJ, IJ2, I, J, IC, JSN, J1, IP1
      INTEGER(KIND=JWIM) :: ICHNK, IPRM
      INTEGER(KIND=JWIM) :: IUME, IOSZX
      INTEGER(KIND=JWIM) :: NBL, IOBL , IBL
      INTEGER(KIND=JWIM) :: NN, NUWH, IOBS
      INTEGER(KIND=JWIM) :: IOBCT, NSATBC, MAXNBCINC 
      INTEGER(KIND=JWIM) :: IBCT, IBCS, MAXBCT
      INTEGER(KIND=JWIM) :: ISAT, KSAT, IDUM, NBCINC
      INTEGER(KIND=JWIM) :: IBLKLIST_SWH
      INTEGER(KIND=JWIM) :: IDENTI, ISENSOR, IDES_HS, IDES_WS
      INTEGER(KIND=JWIM) :: NRGG
      INTEGER(KIND=JWIM) :: ISHIFT
      INTEGER(KIND=JWIM) :: IE, IE2
      INTEGER(KIND=JWIM), DIMENSION(NUMALT) :: ICOUNT, ICOUNTT, ICOUNTD
      INTEGER(KIND=JWIM), DIMENSION(NUMALT) :: ICOUNTB, ICOUNTQCF, ICOUNTBL, IJOLD
      INTEGER(KIND=JWIM), DIMENSION(MAXSATBC) :: ISATBC, NBCT
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:) :: STATID, SENSOR, PARAM
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:) :: NBCINCMAX
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:) :: I_J_TOIJ
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:) :: IFLAG 
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:) :: NUMBWH, NR

      REAL(KIND=JWRB) :: HSMIN
      REAL(KIND=JWRB) :: WSMIN
      REAL(KIND=JWRB) :: BCINC0
      REAL(KIND=JWRB) :: DX, XI
      REAL(KIND=JWRB) :: RLAT, RLON, SWH, WS, WHSE2, RFWM
      REAL(KIND=JWRB) :: XX, ZI, BCIV
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      REAL(KIND=JWRU) :: XLA, XLO
      REAL(KIND=JWRU) :: XLAC, XLOC, DISTMIN
      REAL(KIND=JWRU) :: XLAC2, XLOC2, DISTMIN2

      REAL(KIND=JWRB), DIMENSION(NUMALT) :: SWHO
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:) :: LATMI, LONGMI, LATMA, LONGMA
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:) :: VALUE_MIN, VALUE_MAX, DLON
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:) :: BCINC
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:,:) :: BCTABLE
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:) :: WHME, WHSE
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:) :: XWSTMP
      REAL(KIND=JWRU), ALLOCATABLE, DIMENSION(:,:) :: XLATTMP, XLONTMP

      CHARACTER (LEN=3) :: YOFID
      CHARACTER (LEN=14) :: CBEGINDT, CENDDT, CDATE
      CHARACTER (LEN=19) :: YOFILE
      CHARACTER (LEN=27) :: CDATES
      CHARACTER (LEN=200) :: CDUM
      CHARACTER (LEN=14) :: CDATEO(NUMALT)
      CHARACTER (LEN=12), ALLOCATABLE, DIMENSION(:) :: CDATESTART
      CHARACTER (LEN=12), ALLOCATABLE, DIMENSION(:) :: CDATEEND
      CHARACTER (LEN=12), ALLOCATABLE, DIMENSION(:,:) :: CBCDATESTART
      CHARACTER (LEN=12), ALLOCATABLE, DIMENSION(:,:) :: CBCDATEEND
      CHARACTER (LEN=14), ALLOCATABLE, DIMENSION(:,:) :: CDATETMP

      LOGICAL :: LEOFD
      LOGICAL :: LOPENED
      LOGICAL :: LLBLEXIST
      LOGICAL :: LLBCTEXIST, LLDOBC, LNEW
      LOGICAL :: LBLKLIST_SWH

!*     VARIABLE     TYPE     PURPOSE.
!      --------     ----     --------

!     *CBEGINDT*    CHAR*14  BEGINNING DATE OF THE TIME WINDOW
!                            DURING WHICH DATA ARE USED FOR ANALYSIS.
!     *CENDDT*      CHAR*14  END DATE OF THE TIME WINDOW
!                            DURING WHICH DATA ARE SAVE FOR ANALYSIS.
!     *WHME*        REAL     GRID WAVE HEIGHT.
!     *WHSE*        REAL     WAVE HEIGHT ROOT MEAN SQUARE ERROR.
!     *NUMBWH*      INTEGER  NUMBER OF HS MEASUREMENTS INSIDE A BOX.

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('GRDATA',0,ZHOOK_HANDLE)

!*    1. FIX TIME WINDOW EQUAL TO ONE WIND TIME STEP.
!        --------------------------------------------

      CBEGINDT = CDTPRO
      IF (IDATAWL <= 10800) THEN
        CALL INCDATE (CBEGINDT,-IDATAWL)
      ELSE 
        CALL INCDATE (CBEGINDT,-IDATAWL/2)
      ENDIF
      CENDDT = CBEGINDT
      CALL INCDATE (CENDDT,IDATAWL)

      HSMIN=0.05_JWRB
      WSMIN=0.01_JWRB

      WRITE(IU06,*) ' GRDATA : SELECTING ALTIMETER MEASUREMENTS '
      WRITE(IU06,*) '          PERIOD FROM ', CBEGINDT,' TO ', CENDDT

!     FIX THE MINIMUM NUMBER OF OBSERVATIONS PER GRID BOX REQUIRED
!     IN ORDER TO ACCEPT THE DATA AT THE GRIB POINTS

      IF (LAVERAGE) THEN
!       if superobservations fall within the same grid box, 
!       they will be averaged.
        NNMIN=1
      ELSEIF (LLUNSTR) THEN
        NNMIN=1
      ELSEIF (XDELLA >= 0.5_JWRB) THEN
        NNMIN=4
      ELSEIF (XDELLA > 0.25_JWRB) THEN
        NNMIN=3
      ELSE
        NNMIN=2
      ENDIF
      WRITE(IU06,*) '   GRIDDED DATA WILL BE KEPT WHEN THE NUMBER OF'
      WRITE(IU06,*) '   OBSERVATIONS PER GRID BOX IS LARGER THAN ',NNMIN-1
      WRITE(IU06,*) ' '

      IF (LAVERAGE) THEN
        YOFID='AFA'
      ELSE 
        YOFID='QCP'
      ENDIF
! ----------------------------------------------------------------------

!*    2. INITIALIZE ARRAYS, counters and units.
!        --------------------------------------

      IF (.NOT. LLUNSTR) THEN
        ALLOCATE(I_J_TOIJ(NGX,NGY))
        I_J_TOIJ(:,:) = 0 

        DO IJ = 1, NIBLO 
          I = BLK2GLO%IXLG(IJ)
          J = NGY-BLK2GLO%KXLT(IJ)+1
          I_J_TOIJ(I,J) = IJ
        ENDDO

      ENDIF

      ALLOCATE(IFLAG(NIBLO,NUMALT))
      ALLOCATE(WHME(NIBLO,NUMALT))
      ALLOCATE(WHSE(NIBLO,NUMALT))
      ALLOCATE(NUMBWH(NIBLO,NUMALT))
      ALLOCATE(NR(NIBLO,NUMALT))

      ! XLATTMP(real8), XLONTMP(real8), XWSTMP, CDATETMP   ARE
      ! TEMPORARY ARRAYS TO HOLD INDIVIDUAL TIME AND LOCATION 
      ! OF THE LAST OBSERVATION AT A GIVEN MODEL GRID POINT.
      ! THIS IS OK IF THERE IS NO AVERAGING INVOLVED.
      ! HOWEVER, IT SHOULD BE CORRECTED LATER TO TAKE ORIGINAL 
      ! OBSERVATIONS INTO ACCOUNT.
      ALLOCATE(XLATTMP (NIBLO,NUMALT))
      ALLOCATE(XLONTMP (NIBLO,NUMALT))
      ALLOCATE(XWSTMP  (NIBLO,NUMALT))
      ALLOCATE(CDATETMP(NIBLO,NUMALT))

      DO ISAT=1,NUMALT
        DO IJ=1,NIBLO
          IFLAG(IJ,ISAT) = -1 
          WHME(IJ,ISAT) = 0._JWRB
          WHSE(IJ,ISAT) = 0._JWRB
          NUMBWH(IJ,ISAT) = 0
          NR(IJ,ISAT) = 0
        ENDDO
        ICOUNT(ISAT)=0
        ICOUNTT(ISAT)=0
        ICOUNTD(ISAT)=0
        ICOUNTB(ISAT)=0
        ICOUNTQCF(ISAT)=0
        ICOUNTBL(ISAT)=0
        IJOLD(ISAT) = 0
        IBUFRSAT(ISAT)=0 
        CDATEO(ISAT) = CBEGINDT
      ENDDO


      DO IUME=33,99,1
        INQUIRE ( UNIT=IUME, OPENED=LOPENED)
        IF ( .NOT. LOPENED) THEN
          WRITE(IU06,*) ' GRDATA : TAKES UNIT ',IUME,' AS ALTIMETER DATA INPUT UNIT'
          WRITE(IU06,*) ' '
          GOTO 2000
        ENDIF
      ENDDO

!     CHECK IF BLACKLISTING FILE IS PRESENT
2000  CONTINUE
      NBL=0
      YOFILE='wave_data_blacklist'
      INQUIRE (FILE=YOFILE,EXIST=LLBLEXIST)
      IF (LLBLEXIST) THEN
        IOBL = IWAM_GET_UNIT(IU06, YOFILE, 'r', 'f', 0, 'READ')
        WRITE(IU06,*) '   GRDATA : READING BLACKLISTING FILE ', YOFILE
        WRITE(IU06,*) ' '

        DO IC=1,21
          READ(IOBL,'(A200)',END=200,ERR=200) CDUM
        ENDDO

        DO WHILE(.TRUE.)
          READ(IOBL,'(A200)',END=100,ERR=100) CDUM
          NBL=NBL+1
        ENDDO
100     WRITE(IU06,*) '            TOTAL NUMBER OF ENTRIES :',NBL
        WRITE(IU06,*) ' '
        REWIND(IOBL)

        IF (NBL > 0) THEN
          ALLOCATE(STATID(NBL))
          ALLOCATE(SENSOR(NBL))
          ALLOCATE(CDATESTART(NBL))
          ALLOCATE(CDATEEND(NBL))
          ALLOCATE(LATMI(NBL))
          ALLOCATE(LONGMI(NBL))
          ALLOCATE(LATMA(NBL))
          ALLOCATE(LONGMA(NBL))
          ALLOCATE(DLON(NBL))
          ALLOCATE(PARAM(NBL))
          ALLOCATE(VALUE_MIN(NBL))
          ALLOCATE(VALUE_MAX(NBL))
        ENDIF

        DO IC=1,21
          READ(IOBL,'(A200)',END=200,ERR=200) CDUM
        ENDDO

        DO IBL=1,NBL
          READ(IOBL,111,ERR=200) STATID(IBL), SENSOR(IBL),              &
     &         CDATESTART(IBL), CDATEEND(IBL), LATMI(IBL), LONGMI(IBL), &
     &         LATMA(IBL), LONGMA(IBL), PARAM(IBL), VALUE_MIN(IBL),     &
     &         VALUE_MAX(IBL)
        ENDDO
        GOTO 300
200     CONTINUE
        WRITE(IU06,*) '************************************************'
        WRITE(IU06,*) '   GRDATA : ERROR READING ', YOFILE
        WRITE(IU06,*) '   THE PROGRAM WILL ABORT !!!!'
        WRITE(IU06,*) '************************************************'
        CALL ABORT

300     CLOSE(IOBL)
        DO IBL=1,NBL
          IF (ABS(LONGMI(IBL)-LONGMA(IBL)) /= 360._JWRB) THEN
            LONGMI(IBL) = MOD(LONGMI(IBL)+720._JWRB,360._JWRB)
            LONGMA(IBL) = MOD(LONGMA(IBL)+720._JWRB,360._JWRB)
          ENDIF
          IF (LONGMI(IBL) >= LONGMA(IBL)) LONGMI(IBL)=LONGMI(IBL)-360._JWRB
          DLON(IBL) =  LONGMA(IBL) - LONGMI(IBL)

        ENDDO
      ELSE
        WRITE(IU06,*) '************************************************'
        WRITE(IU06,*) '   GRDATA : THE BLACK LISTING FILE ',YOFILE
        WRITE(IU06,*) '   WAS NOT FOUND !!!'
        WRITE(IU06,*) '   THE PROCESSING WILL CONTINUE WITHOUT IT !!!!'
        WRITE(IU06,*) '************************************************'
      ENDIF
111   FORMAT(2(I8,1x),2(A12,1x),2(f5.1,1x,f6.1,1x),I6,1x,2(f9.3,1x))

! ----------------------------------------------------------------------
!     CHECK IF THE FILE WITH "BIAS CORRECTION TABLES" IS PRESENT
2500  CONTINUE
      NBCT=0
      YOFILE='biascorrection.swh'
      INQUIRE (FILE=YOFILE,EXIST=LLBCTEXIST)
      IF (LLBCTEXIST) THEN
        IOBCT = IWAM_GET_UNIT(IU06, YOFILE, 'r', 'f', 0, 'READ')
        WRITE(IU06,*) '   GRDATA : READING FILE WITH BIAS CORRECTION '  &
     &                    //'TABLES  ', YOFILE

        WRITE(IU06,*) ' '

        DO IC=1,21
          READ(IOBCT,'(A)',END=2800,ERR=2800) CDUM
        ENDDO

        NSATBC=0
        MAXNBCINC=0
        DO WHILE(.TRUE.)
          READ(IOBCT,'(I4,1X,I4,A27,I4,1X,F6.3)',END=2600,ERR=2600)     &
     &             ISAT, IDUM, CDATES, NBCINC, BCINC0
          MAXNBCINC=MAX(MAXNBCINC,NBCINC)
          LNEW=.TRUE.
          DO IBCS=1, NSATBC
            IF (ISAT == ISATBC(IBCS)) THEN
              LNEW=.FALSE.
              NBCT(IBCS) = NBCT(IBCS) + 1
              EXIT
            ENDIF
          ENDDO
          IF (LNEW) THEN
            NSATBC=NSATBC+1
            IF (NSATBC > MAXSATBC) THEN
              WRITE(IU06,*) '******************************************'
              WRITE(IU06,*) '   GRDATA : ERROR:  NUMBER OF SATELLITES ' &
     &             //'IN BIAS CORRECTION FILE EXCEEDS MAXSATBC =',      &
     &                                                MAXSATBC 
              WRITE(IU06,*) '       INCREASE MAXSATBC IN GRDATA'
              WRITE(IU06,*) '   THE PROGRAM WILL ABORT !!!!'
              WRITE(IU06,*) '******************************************'
              CALL ABORT
            ENDIF
            ISATBC(NSATBC)=ISAT
            NBCT(NSATBC) = 1
          ENDIF
          READ(IOBCT,'(A)',END=2800,ERR=2800) CDUM
        ENDDO
2600    WRITE(IU06,*) '            NUMBER OF SATELLITE ENTRIES :',NSATBC
        MAXBCT=NBCT(1)
        DO IBCS=1, NSATBC
          WRITE(IU06,*) '              SATELLITE ',ISATBC(IBCS),        &
     &                                  ' HAS ', NBCT(IBCS), ' TABLES'
          MAXBCT=MAX(MAXBCT,NBCT(IBCS))
        ENDDO
        REWIND(IOBCT)

        IF (NSATBC > 0) THEN
          ALLOCATE(NBCINCMAX(NSATBC,MAXBCT))
          ALLOCATE(BCINC(NSATBC,MAXBCT))
          ALLOCATE(BCTABLE(NSATBC,MAXBCT,MAXNBCINC))
          ALLOCATE(CBCDATESTART(NSATBC,MAXBCT))
          ALLOCATE(CBCDATEEND(NSATBC,MAXBCT))
        ENDIF

        DO IC=1,21
          READ(IOBCT,'(A)',END=2800,ERR=2800) CDUM
        ENDDO

        DO IBCS=1, NSATBC
          NBCT(IBCS) = 0
        ENDDO
        DO WHILE(.TRUE.)
          READ(IOBCT,'(I4,1X,I4,A27,I4,1X,F6.3)',END=2900,ERR=2900)     &
     &          ISAT, IDUM, CDATES, NBCINC, BCINC0
          DO IBCS=1, NSATBC
            IF (ISAT == ISATBC(IBCS)) THEN
              NBCT(IBCS) = NBCT(IBCS) + 1
              NBCINCMAX(IBCS,NBCT(IBCS))=NBCINC
              BCINC(IBCS,NBCT(IBCS))=BCINC0
              CBCDATESTART(IBCS,NBCT(IBCS))=CDATES(2:13)
              CBCDATEEND  (IBCS,NBCT(IBCS))=CDATES(15:26)
              READ(IOBCT,'(9999F6.3)',END=2800,ERR=2800)                &
     &             (BCTABLE(IBCS,NBCT(IBCS),I),I=1,NBCINC)
              EXIT
            ENDIF
          ENDDO
        ENDDO

2800    CONTINUE
        WRITE(IU06,*) '************************************************'
        WRITE(IU06,*) '   GRDATA : ERROR READING ', YOFILE
        WRITE(IU06,*) '   THE PROGRAM WILL ABORT !!!!'
        WRITE(IU06,*) '************************************************'
        CALL ABORT

2900    CONTINUE
        CLOSE (IOBCT)
        WRITE(IU06,*) '            BIAS CORRECTION TABLES READ IN'
        WRITE(IU06,*) ' '

      ELSE
        WRITE(IU06,*) '************************************************'
        WRITE(IU06,*) '   GRDATA : FILE WITH BIAS CORRECTION TABLES ', YOFILE
        WRITE(IU06,*) '   WAS NOT FOUND !!!'
        WRITE(IU06,*) '   THE PROCESSING WILL CONTINUE WITHOUT IT !!!!'
        WRITE(IU06,*) '************************************************'
      ENDIF
! ----------------------------------------------------------------------

!*    3. READ THE MEASUREMENT TAPE AND CUMULATE THE DATA ON GRID POINTS.
!        ---------------------------------------------------------------

 3000 CONTINUE
      CALL READSAT (IU06, IUME, CDTPRO, IDENTI, ISENSOR,                &
     &              CDATE, RLAT, RLON, IDES_HS, SWH,                    &
     &              IDES_WS, WS, LEOFD, YOFID)

      IF (LEOFD) GOTO 4000

!     WHICH SATELLITE IS IT
      DO ISAT=1,NUMALT
        IF (IBUFRSAT(ISAT) == 0) THEN
          IBUFRSAT(ISAT)=IDENTI
          KSAT=ISAT
          EXIT
        ENDIF
        IF (IBUFRSAT(ISAT) == IDENTI) THEN
          KSAT=ISAT
          EXIT
        ENDIF
      ENDDO
      IF (ISAT > NUMALT) THEN
        WRITE(IU06,*) ' GRDATA : WARNING MORE SATELLITES IN INPUT THAN ALLOWED BY NUMALT !!!'
        WRITE(IU06,*) '          THE SATELLITE IS ',IDENTI
        WRITE(IU06,*) '          NEW DATA WILL NOT BE ALLOWED IN ' 
        GOTO 3000
      ENDIF

      ICOUNTD(KSAT) = ICOUNTD(KSAT)+1

!     CHECK DATE.

      IF (CDATE < CBEGINDT)  then
        ICOUNTB(KSAT) = ICOUNTB(KSAT) + 1

      ELSE IF (CDATE <= CENDDT) THEN
!       DATE IS INSIDE THE TIME WINDOW AND RECORD FLAGGED AS RELIABLE.
!       THE INDICES ARE COMPUTED

        ICOUNTT(KSAT) = ICOUNTT(KSAT)+1

!       FIND THE MODEL INDEX CLOSEST TO THE OBSERVATION THAT IS WITHIN AN ELEMENT (IE)
        IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
!       UNSTRUCTURED GRID:
!       We will limit to observations that are within the elements
          XLA = REAL(RLAT,JWRU)
          XLO = REAL(RLON,JWRU)
          CALL IPELEMENT_CLOSEST(XLO,XLA,IJ,IE,XLOC,XLAC,DISTMIN)

!!!    the find_element does not do periodicity very well (will need to change !!!!)
          IF (IJ <= 0) THEN
            XLO = REAL(RLON,JWRU)-360._JWRU
            CALL IPELEMENT_CLOSEST(XLO,XLA,IJ,IE,XLOC,XLAC,DISTMIN)
          ENDIF

          IF (IJ <= 0) THEN
            XLO = REAL(RLON,JWRU)+360._JWRU
            CALL IPELEMENT_CLOSEST(XLO,XLA,IJ,IE,XLOC,XLAC,DISTMIN)
          ENDIF
#else
          CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif

        ELSE
!       STRUCTURED GRID:
          J       = NGY-NINT((RLAT-AMOSOP)/XDELLA)
          XI      = MOD(RLON-AMOWEP+720._JWRB,360._JWRB)
          JSN     = NGY-J+1
          IF (J <= 0 .OR. J > NGY) JSN=1
          NRGG    = NLONRGG(JSN)
          I       = NINT(XI/ZDELLO(JSN)) + 1
          IF (IPER == 1 .AND. I == NRGG+1) I = 1

!         IF THE INDICES ARE INSIDE THE GRID THE MEASUREMENT IS CUMULATED
          IF (I <= NRGG .AND. I >= 1 .AND. J <= NGY .AND. J > 1) THEN
            IJ=I_J_TOIJ(I,J)
          ELSE
            IJ=0
          ENDIF
        ENDIF

!       DATA ON THE MODEL GRID
        IF (IJ > 0) THEN
          ICOUNT(KSAT)=ICOUNT(KSAT)+1
          IF (LLUNSTR) THEN
            XLATTMP (IJ,KSAT)=XLAC
            XLONTMP (IJ,KSAT)=XLOC
          ELSE
            XLATTMP (IJ,KSAT)=REAL(RLAT,JWRU)
            XLONTMP (IJ,KSAT)=REAL(RLON,JWRU)
          ENDIF
          XWSTMP  (IJ,KSAT)=WS
          CDATETMP(IJ,KSAT)=CDATE

!         BLACKLISTED ?
          LBLKLIST_SWH=.FALSE.
          IF (LLBLEXIST) THEN
            DO IBL=1,NBL
              IF (STATID(IBL) == IDENTI .AND.                           &
     &           SENSOR(IBL) == ISENSOR .AND.                           &
     &           CDATESTART(IBL) <= CDATE(1:12) .AND.                   &
     &           CDATEEND(IBL) >= CDATE(1:12) .AND.                     &
!             BLACKLISTING OF WAVE HEIGHT ONLY !!
     &           PARAM(IBL) == IDES_HS .AND.                            &
     &           LATMI(IBL) <= RLAT .AND.                               &
     &           LATMA(IBL) >= RLAT .AND.                               &
     &           VALUE_MIN(IBL) <= SWH .AND.                            &
     &           VALUE_MAX(IBL) >= SWH                                  &
     &          ) THEN
                DX = MOD(RLON-LONGMI(IBL)+720._JWRB,360._JWRB)
                IF (DX <= DLON(IBL)) THEN
                  ICOUNTBL(KSAT)=ICOUNTBL(KSAT)+1
                  LBLKLIST_SWH=.TRUE.
                ENDIF
              ENDIF
            ENDDO
          ENDIF


          IF (LAVERAGE) THEN
!           the superobservation is assigned to the closest grid point.
            IF (SWH >= 0.0_JWRB) THEN
!             if blacklisted, mark the grid box for that satellite as blacklisted
              IF (LBLKLIST_SWH) IFLAG(IJ,ISAT) = 0
              WHME(IJ,KSAT) = WHME(IJ,KSAT) + SWH
              NUMBWH(IJ,KSAT) = NUMBWH(IJ,KSAT) + 1
            ENDIF
          ELSE
!           grid box average
            IF (IJ /= IJOLD(KSAT)) THEN
!               NEW BOX
              CALL DIFDATE(CDATEO(KSAT),CDATE,ISHIFT)
              IF (ISHIFT < 3 .AND. ICOUNT(KSAT) > 1) THEN
                IF (SWH-SWHO(KSAT) > 2.0_JWRB) THEN
                    NR(IJ,KSAT)=NR(IJ,KSAT)+1
                  GO TO 3000
                ENDIF
              ENDIF
              IJOLD(KSAT) = IJ
              SWHO(KSAT) = SWH
            ELSE
!             SUDDEN POSITVE STEPS ARE REFUSED.
              IF (SWH-SWHO(KSAT) > 1.0_JWRB) THEN
                NR(IJ,KSAT) = NR(IJ,KSAT)+1
                GO TO 3000
              ENDIF
            ENDIF

            IF (SWH >= 0._JWRB) THEN
!             if blacklisted, mark the grid box for that satellite as blacklisted
              IF (LBLKLIST_SWH) IFLAG(IJ,ISAT) = 0
              WHME(IJ,KSAT) = WHME(IJ,KSAT) + SWH
              WHSE (IJ,KSAT) = WHSE(IJ,KSAT) + SWH**2
              NUMBWH(IJ,KSAT) = NUMBWH(IJ,KSAT) + 1
            ENDIF

!           A PREVIOUS SUDDEN NEGATIVE STEP IS COMPENSATED.
            IF (SWHO(KSAT)-SWH > 2.0_JWRB) THEN
              NR(IJ,KSAT) = NR(IJ,KSAT)+1
              WHME(IJ,KSAT) = WHME(IJ,KSAT) - SWHO(KSAT)
              WHSE (IJ,KSAT) = WHSE(IJ,KSAT) - SWHO(KSAT)**2
              NUMBWH(IJ,KSAT) = NUMBWH(IJ,KSAT) - 1
            ENDIF

            SWHO(KSAT) = SWH
            CDATEO(KSAT) = CDATE
          ENDIF

        ENDIF ! VALID IJ

      ENDIF  ! end test whether obs in analsis window
      GOTO 3000

! ----------------------------------------------------------------------

!*    4. DETERMINE AVERAGE VALUE AND ERROR OF THE MEASUREMENTS IN THE
!*       BOX AROUND GRIDPOINT I,J.
!        ------------------------------------------------------------

 4000 CONTINUE

      NWH = 0

      DO ISAT=1,NUMALT
        IF (IBUFRSAT(ISAT) /= 0) THEN
          NUWH=0
          DO IJ = 1, NIBLO
            NUWH=NUWH+NUMBWH(IJ,ISAT)
          ENDDO

          WRITE (IU06,'(" GRDATA: INPUT STATISTICS:")')
          WRITE (IU06,'("         THE SATELLITE IS",                    &
     &                  " ......:", I8)') IBUFRSAT(ISAT)
          WRITE (IU06,'("         NUMBER OF RECORDS READ FROM DATA ",   &
     &                  "FILE IS ......:", I8)') ICOUNTD(ISAT)
          WRITE (IU06,'("         THERE OF NUMBER OF BLACKLISTED ",      &
     &                  "RECORDS IS .....:", I8)') ICOUNTBL(ISAT)
          WRITE (IU06,'("         NUMBER OF DATA IN SELECTED PERIOD ",  &
     &                  "IS ..........:", I8)') ICOUNTT(ISAT)
          WRITE (IU06,'("         NUMBER OF DATA BEFORE SELECTED ",     &
     &                  "PERIOD IS ......:", I8)') ICOUNTB(ISAT)
          WRITE (IU06,'("         NUMBER OF DATA IN BOXES AROUND ",     &
     &                  "OVER SEA IS ..:", I8)') ICOUNT(ISAT)
          WRITE (IU06,'("         NUMBER OF GOOD DATA ABOVE LAND IS ",  &
     &                  ".............:", I8)') ICOUNTQCF(ISAT)
          WRITE (IU06,'("         NUMBER OF DATA RETAINED BY ",         &
     &                  "''GRDATA'' IS FOR HS .:", I8)') NUWH
          WRITE (IU06,*) ' '

          DO IJ = 1, NIBLO
            NN = NUMBWH(IJ,ISAT)
            IF (NN >= NNMIN) THEN
!             MEAN VALUES
              WHME(IJ,ISAT) = WHME(IJ,ISAT)/REAL(NN,JWRB)
              IF ( IFLAG(IJ,ISAT) /= 0 ) IFLAG(IJ,ISAT) = 1 
!             COMPUTE THE STANDARD DEVIATION.
              IF (.NOT. LAVERAGE) THEN
                WHSE2 =  (WHSE(IJ,ISAT)-WHME(IJ,ISAT)**2*REAL(NN,JWRB)) &
     &                   / (REAL(NN,JWRB)-1._JWRB)
                WHSE(IJ,ISAT) = SQRT( MAX(WHSE2,0._JWRB))
              ENDIF
            ENDIF
          ENDDO

! ----------------------------------------------------------------------

!*    5. FLAG GRID POINTS, IF THERE ARE FEW OR STRONGLY SCATTERED DATA.
!        --------------------------------------------------------------

!         ONLY KEEP THE RELEVANT INFORMATION
!         DATA THAT HAVE FAILED THE BASIC QC WILL BE KEPT AND FLAGGED WITH
!         IJALT(:,3) < 0
!         BLACKLISTED DATA WILL BE FLAGGED WITH IJALT(:,3) = 0

          IF (LAVERAGE) THEN
            DO IJ = 1, NIBLO
              NN = NUMBWH(IJ,ISAT)
              IF (NN == 0) THEN
!               IF THERE ARE NO MEASUREMENTS.
                IFLAG(IJ,ISAT) = -1
              ELSEIF (WHME(IJ,ISAT) < HSMIN) THEN
!               HS BELOW MINIMUM
                IFLAG(IJ,ISAT) = -2
                NWH = NWH+1
              ELSE
                IF ( IFLAG(IJ,ISAT) /= 0 ) IFLAG(IJ,ISAT) = 1 
                NWH = NWH+1
              ENDIF
            ENDDO
          ELSE

            DO IJ = 1, NIBLO
              NN = NUMBWH(IJ,ISAT)
              IF (NN == 0) THEN
!               IF THERE ARE NO MEASUREMENTS.
                IFLAG(IJ,ISAT) = -1 
              ELSEIF (NN > 0 .AND. NN < NNMIN) THEN
!               IF THERE ARE FEW MEASUREMENTS.
                IFLAG(IJ,ISAT) = -3
                NWH = NWH+1
              ELSEIF (NN >= NNMIN) THEN
!               IF THE VARIANCE IS TOO LARGE.
!               OR SWH IS TOO SMALL  OR THERE ARE TOO MANY SPIKES
                XX = REAL(NR(IJ,ISAT),JWRB)/REAL(NN,JWRB)
                RFWM = MAX(0.5_JWRB,0.25_JWRB*WHME(IJ,ISAT))
                IF (WHSE(IJ,ISAT) > RFWM .OR.                          &
     &              WHME(IJ,ISAT) < HSMIN .OR.                         &
     &              XX > 0.1_JWRB) THEN
                  IFLAG(IJ,ISAT) = -4
                  NWH = NWH+1
                ELSE
                  NWH = NWH+1
                ENDIF
              ENDIF
            ENDDO
          ENDIF ! LAVERAGE

        ENDIF  ! IBUFRSAT 
      ENDDO  ! ON ISAT

      DEALLOCATE(WHSE)
      DEALLOCATE(NUMBWH)
      DEALLOCATE(NR)

      IF (NWH > 0) THEN
        ALLOCATE(IJALT(NWH,NIJALT))
        IJALT(:,:)=0
        ALLOCATE(ALTDATA(NWH,NALTDT))
        ALTDATA(:,:)=0._JWRB
        ALLOCATE(ALTEXDATA(NWH,NALTEDT))
        ALTEXDATA(:,:)=0._JWRB
        IF (LLUNSTR) THEN
          ALLOCATE(ALTUNDATA(NWH,NALTUDT))
          ALTUNDATA(:,:)=0._JWRU
        ENDIF
        ALLOCATE(CDATEOBS (NWH))

        IOBS=0
        DO ISAT=1,NUMALT
          LLDOBC=.FALSE.
          IF (IBUFRSAT(ISAT) /= 0) THEN
            IF (LLBCTEXIST) THEN
!             loop over all satellites
              DO IBCS=1, NSATBC
                IF (IBUFRSAT(ISAT) == ISATBC(IBCS)) THEN
!                 loop over number of table for a given satellite
                  DO IBCT=1,NBCT(IBCS)
                    IF (CDTPRO(1:12) >= CBCDATESTART(IBCS,IBCT) .AND.   &
     &                  CDTPRO(1:12) <= CBCDATEEND  (IBCS,IBCT)) THEN
                      LLDOBC=.TRUE.
                      EXIT
                    ENDIF
                  ENDDO
                  EXIT
                ENDIF
              ENDDO
            ENDIF
            IF (LLDOBC) THEN
              WRITE(IU06,*) '  BIAS CORRECTION FOR SATELLITE ',IBUFRSAT(ISAT)
            ENDIF

            DO IJ = 1, NIBLO
              IF (WHME(IJ,ISAT) > 0._JWRB) THEN
                IOBS=IOBS+1
                ALTDATA(IOBS,3) = WHME(IJ,ISAT)  !  1 REPLACED BY 3
                IJALT(IOBS,1) = IJ
                IJALT(IOBS,2) = IBUFRSAT(ISAT)
                IJALT(IOBS,3) = IFLAG(IJ,ISAT) 

                ! NOTE THAT WE NEED TO PASS ALL DATA TO ODB EVEN THOSE REJECTED.
                ! NEW FLAGS NEED TO BE DESIGNED FOR THAT.
                ! FOR THE TIME BEING ONLY ACCEPTED OBS ARE PASSED TO ODB.
                ALTEXDATA(IOBS,1) = REAL(XLATTMP(IJ,ISAT),JWRB)
                ALTEXDATA(IOBS,2) = REAL(XLONTMP(IJ,ISAT),JWRB)
                ALTEXDATA(IOBS,3) = XWSTMP(IJ,ISAT)
                CDATEOBS (IOBS)   = CDATETMP(IJ,ISAT)

                IF (LLUNSTR) THEN
                  ALTUNDATA(IOBS,1) = XLATTMP(IJ,ISAT)
                  ALTUNDATA(IOBS,2) = XLONTMP(IJ,ISAT)
                ENDIF 

                !! NO REAl QC ON WIND SPEED YET !!!
                !! only that there is a wave height and wind speed > WSMIN
                IF (XWSTMP(IJ,ISAT) > WSMIN) THEN
                  !! wind data are passive
                  IJALT(IOBS,4) = 0 
                ELSE
                  IJALT(IOBS,4) = -2 
                ENDIF

                ! BIAS CORRECTION IS DONE IF NEEDED
                IF (LLDOBC) THEN
                  ZI=MIN(REAL(NBCINCMAX(IBCS,IBCT),JWRB), WHME(IJ,ISAT)/BCINC(IBCS,IBCT))
                  I=INT(ZI)+1
                  I=MIN(NBCINCMAX(IBCS,IBCT), I)
                  IP1=MIN(NBCINCMAX(IBCS,IBCT), I+1)
                  IF (I == NBCINCMAX(IBCS,IBCT)) THEN
                    BCIV=BCTABLE(IBCS,IBCT,I)
                  ELSE
                    BCIV=(I-ZI)*BCTABLE(IBCS,IBCT,I)+(ZI-I+1)  *BCTABLE(IBCS,IBCT,IP1)
                  ENDIF
!                 Removing the bias
                  ALTDATA(IOBS,1)=WHME(IJ,ISAT)-BCIV
                ELSE
                  ALTDATA(IOBS,1)=WHME(IJ,ISAT)
                ENDIF

              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF

      DEALLOCATE(IFLAG)
      DEALLOCATE(WHME)

      DEALLOCATE(XLATTMP)
      DEALLOCATE(XLONTMP)
      DEALLOCATE(XWSTMP)
      DEALLOCATE(CDATETMP)

      IF (NBL > 0) THEN
        DEALLOCATE(STATID)
        DEALLOCATE(SENSOR)
        DEALLOCATE(CDATESTART)
        DEALLOCATE(CDATEEND)
        DEALLOCATE(LATMI)
        DEALLOCATE(LONGMI)
        DEALLOCATE(LATMA)
        DEALLOCATE(LONGMA)
        DEALLOCATE(DLON)
        DEALLOCATE(PARAM)
        DEALLOCATE(VALUE_MIN)
        DEALLOCATE(VALUE_MAX)
      ENDIF

      IF (ALLOCATED(I_J_TOIJ)) DEALLOCATE(I_J_TOIJ)

      WRITE (IU06,*) ' GRDATA: OVERALL GRIDDED MEASUREMENTS ARE SELECTED:'
      WRITE (IU06,*) '         WAVE HEIGHT AT ', NWH,' GRID POINTS'
      WRITE (IU06,*) ' '

      CLOSE (IUME,IOSTAT=IOSZX)

      IF (LHOOK) CALL DR_HOOK('GRDATA',1,ZHOOK_HANDLE)

      END SUBROUTINE GRDATA
