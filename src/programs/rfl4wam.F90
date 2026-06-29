! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

PROGRAM RFL4WAM

! ----------------------------------------------------------------------

!**** *RFL4WAM* - PREPARES QUALITY CHECKED SATELLITE DATA FOR WAMASSI.

!     B. HANSEN      ECMWF         DECEMBER 1997
!     S. ABDALLA     ECMWF         JUNE     2011 TEXT FILES ARE PRODUCED 
!                                                TO BE USED BY simulobs2odb

!*    PURPOSE.
!     --------

!       SPEED UP ASSIMILATION BY TUNING THE I/O OF OBSERVATIONS.
!       AND PEFORMS BLACKLISTING AND CORRECTION OF THE DATA.

!**   INTERFACE.
!     ----------

!          ---- FORMAL PARAMETERS ----

!    *rfl4wam*  -d date [ -d date ... ] [ -i itest ] [ -a idatawl ]
!                   [ -e satID:weight [ -e ... ] ]
!                   [ -t type ]
!      *date*      STRING   DATE/TIME GROUP YYYYMMDDHHMM.
!      *itest*     INTEGER  TEST OUTPUT LEVEL (*COMMON* *TESTO*).
!      *idatawl*   INTEGER  LENGTH OF THE DATA WINDOW (SECONDS).
!      *satID*     INTEGER  SATELLITE BUFR ID.
!      *weight*    REAL     ERROR WEIGHT FOR CORRESPONDING ALTIMETER
!      *type*      STRING*1 OBSERVATION TYPE (i: INDIVIDUAL, a: AVERAGED)
!      if the file wave_data_blacklist is present, information from it
!      will be used to blacklist data.

!          ---- INPUT/OUTPUT UNITS ---

!           *IU00*   - ERROR STANDARD OUTPUT.
!           *IU05*   - USER INPUT UNIT.
!           *IU06*   - PRINTER OUTPUT.
!           *wam_namelist* with the following information
!             &NALINE
!             NANG=
!             NFRE=
!             NFRE_RED=
!             LLUNSTR =
!             /
!           *wam_grid_tables*   - INPUT  UNIT OF PRECOMPUTED COMMON BLOCKS.
!                      (OUTPUT OF PREPROC).
!           *wam_subgrid*   - INPUT  UNIT OF COMMON UBUF.
!                      (OUTPUT OF PREPROC).
!           *IUME*   - ALTIMETER DATA INPUT UNIT (FILES: RFL$DATE).
!           *IO"     - GRIDDED ALTIMETER DATA OUTPUT UNIT
!                       (FILES: rfl$DATE).
!           *IHDR"   - OUTPUT UNIT - ASCII TXT FILE WITH ODB HEADER FOR RALT OBS
!                       (FILES: hdr$DATE.txt).
!           *IBDY"   - OUTPUT UNIT - ASCII TXT FILE WITH ODB BODY   FOR RALT OBS
!                       (FILES: body$DATE.txt).

!          MESSAGE PASSING NOT CATERED FOR YET.
!          IF MESSAGE PASSING ENABLED THE UNITS WILL BE REASSIGNED.

!     METHOD.
!     -------
!          READ SATELLITE DATA AND DISTRIBUTE ONTO MODEL GRID.
!          FOR EACH COMMAND LINE OPTION '-d date' A CORRESPONDING FILE
!          NAMED RFL$DATE IS EXPECTED TO BE ON CURRENT WORKING
!          DIRECTORY. THE TWO
!          OUTPUT FILES FROM PREPROC WILL BE ATTACHED TO IU07 AND IU08.
!          NO FURTHER USER INPUT IS REQUIRED. FOR EACH FILE RFL$DATE ONE
!          CORRESPONDING FILE rfl$DATE WILL BE GENERATED CONTAINING ONLY
!          RECORDS:
!          A) THE CHARACTER STRING "CDTPRO"
!          B) THE ARRAYS IJALT and ALTDATA

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWALTAS , ONLY : NIJALT   ,IJALT    ,ALTDATA, ALTEXDATA,  ALTUNDATA, CDATEOBS
      USE YOWCMDLO , ONLY : WAM_GETCLO, WAM_GETCLA
      USE YOWCOUP  , ONLY : LWCOU
      USE YOWMAP   , ONLY : AMOWEP   ,AMOSOP   ,AMOEAP   ,AMONOP   , &
     &            XDELLA   ,XDELLO
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPREVIOUS,NNEXT    , &
     &            MPMAXLENGTH,KTAG     ,NPRECR   , &
     &            NPRECI
      USE YOWPARAM , ONLY : LL1D     , LLUNSTR
      USE YOWSTAT  , ONLY : CDTPRO   ,YCLASS   ,IPROPAGS ,LSUBGRID
      USE YOWTEST  , ONLY : IU06     ,ITEST

      USE MPL_MODULE,ONLY : MPL_INIT, MPL_END
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

#ifdef WITH_ODB
      USE PARKIND1  ,ONLY : JPIM, JPRB, JPRD
      use, intrinsic :: iso_c_binding
      use odc_c_binding, only : odb_write_new, odb_start
      use odb2_flag_definitions
      use odbmap_reportype, only : find_obstype_codetype, find_reportype
#endif

! -------------------------------------------------------------------

      IMPLICIT NONE
#include "grdata.intfb.h"
#include "iwam_get_unit.intfb.h"
#include "iniwcst.intfb.h"
#include "mpdecomp.intfb.h"

#ifdef WITH_ODB
#include "rfl2odb.intfb.h"
#endif

#include "wvwaminit.intfb.h"

      INTEGER(KIND=JWIM) :: IDATAWL
      INTEGER(KIND=JWIM) :: IDATES, IOBSER, IMO, IOEMAX
      INTEGER(KIND=JWIM) :: I1, JL, LEN1, IND1, NOBS, IOBS
      INTEGER(KIND=JWIM) :: IOPTVAL, MORARG 
      INTEGER(KIND=JWIM) :: IO
      INTEGER(KIND=JWIM) :: NPR, MAXLEN, IREAD, JLD
      INTEGER(KIND=JWIM) :: NGAUSSW, NLON, NLAT
      INTEGER(KIND=JWIM) :: I4(2)
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:) :: IOESATID, IOERDUM

      REAL(KIND=JWRB) :: PRPLRADI, RSOUTW, RNORTW
      REAL(KIND=JWRU) :: Z8
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: Z4(2)
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:) :: OERR, OERRDUM

      CHARACTER(LEN=1) :: CLOPTLET
      CHARACTER(LEN=3) :: CLL1
      CHARACTER(LEN=10) :: CLOPTS
      CHARACTER(LEN=12) :: CLFMT
      CHARACTER(LEN=128) :: CLARG
      CHARACTER (LEN=12),ALLOCATABLE :: CLDATES(:), CLMD(:)

      LOGICAL :: LSQRT, LLIRANK, LLWVENVI
      LOGICAL :: LAVERAGE 
      LOGICAL :: LLRNL
      LOGICAL :: LNOT_ACTIVE, LNOT_ACTIVE_B
      LOGICAL :: LLCLOSE

      DATA CLOPTS/'d;i;a;e;t;'/
      DATA LAVERAGE /.FALSE./

#ifdef WITH_ODB
!     ODB-2 definitions
      TYPE(C_PTR)                        :: odb_handle, odb_it
      INTEGER(KIND=C_INT)                :: cerr
      CHARACTER(KIND=C_CHAR, LEN=256)    :: config = C_NULL_CHAR
      CHARACTER(KIND=C_CHAR, LEN=256)    :: odb_file
      INTEGER(KIND=JWIM), PARAMETER      :: mdi = 2147483647
      CHARACTER(LEN=256)                 :: one_colname, ctype
      CHARACTER(LEN=256)                 :: bitfield_names, bitfield_sizes
      INTEGER(KIND=C_INT)                :: c_ncolumns
      REAL(KIND=C_DOUBLE), ALLOCATABLE   :: odb_values(:), odb_values_b(:)
      INTEGER(KIND=JWIM)                 :: idate, itime, iseqno_offset
      CHARACTER(LEN=256)                 :: cenv_tmp
      INTEGER(KIND=JWIM)                 :: igroupid, isatid
      INTEGER(KIND=JWIM)                 :: ibufrtype,isubtype
      INTEGER(KIND=JWIM)                 :: iobstype,icodetype
      INTEGER(KIND=JWIM)                 :: isensor,ireportype,isatinst
      INTEGER(KIND=JWIM)                 :: NumVar
#endif

! ----------------------------------------------------------------------

#ifdef WITH_WAMASSI

      CALL MPL_INIT(KOUTPUT=1)

      LHOOK = .TRUE.
      IF (LHOOK) CALL DR_HOOK('RFL4WAM',0,ZHOOK_HANDLE)


!*    1. INITIAL VALUES SET AND CRACK COMMAND LINE.
!        -----------------------------------------

      PRPLRADI=1.0_JWRB
      CALL INIWCST(PRPLRADI)

      ITEST   =     0
      IDATES  =     0
      IOBSER  =     0
      IMO     =    10
      IOEMAX  =    10
      IDATAWL = 21600 ! by default the input ralt data will cover a maximum period of 21600 s = 6 hours
      IPROPAGS = 0
      LSUBGRID = .FALSE.
      IU06 = 6

      ALLOCATE (CLDATES(IMO))
      ALLOCATE (IOESATID(IOEMAX), OERR(IOEMAX))

      CMDLINE: DO
        IOPTVAL=WAM_GETCLO(CLOPTS,CLARG)
        IF (IOPTVAL <= 0 )  THEN
          IF (ITEST >= 1)WRITE (*,'("IOPTVAL = ", I3)') IOPTVAL
          EXIT CMDLINE
        ENDIF

        CLOPTLET=CHAR(IOPTVAL)
        IF (ITEST >= 1) WRITE(*,*) ' OPTION = ', CLOPTLET

!       GETS VARIABLE ARGUMENT FOR OPTION
        MORARG=WAM_GETCLA(CLARG)
        IF (ITEST >= 1) WRITE(*,*) ' MORARG : ', MORARG
        IF (MORARG /= 0) THEN
          IF ( CLOPTLET == 'i' ) THEN

! ITEST.
! ------
            I1=LEN_TRIM(CLARG)
            WRITE (CLL1,'(I3)') I1
            CLFMT = '(I'//CLL1//')'
            READ (CLARG(1:I1),FMT=CLFMT) ITEST
          ELSEIF ( CLOPTLET == 'a' ) THEN

! IDATAWL.
! --------
            I1=LEN_TRIM(CLARG)
            WRITE (CLL1,'(I3)') I1
            CLFMT = '(I'//CLL1//')'
            READ (CLARG(1:I1),FMT=CLFMT) IDATAWL
          ELSEIF ( CLOPTLET == 't' ) THEN

! TYPE: INDIVIDUAL OR AVERAGED.
! -----------------------------
            IF (CLARG(1:1) == 'a' .OR. CLARG(1:1) == 'A') LAVERAGE=.TRUE. 

          ELSEIF ( CLOPTLET == 'e' ) THEN

! OBSERVATION ERROR ARRAY.
! ------------------------
            IOBSER = IOBSER + 1
            IF (IOBSER > IOEMAX) THEN
              ALLOCATE(IOERDUM(IOEMAX), OERRDUM(IOEMAX))
              IOERDUM = IOESATID
              OERRDUM = OERR
              DEALLOCATE(IOESATID, OERR)
              IOEMAX = IOEMAX + 10
              ALLOCATE (IOESATID(IOEMAX), OERR(IOEMAX))
              DO JL = 1, IOEMAX-10
                OERR(JL) = OERRDUM(JL)
                IOESATID(JL) = IOERDUM(JL)
              ENDDO
              DEALLOCATE(IOERDUM, OERRDUM)
            ENDIF
            LEN1=LEN_TRIM(CLARG)
            IND1=INDEX(CLARG, ':')
            IF (IND1 <= 1) THEN
              IDATES = 0
              EXIT CMDLINE
            ENDIF
            WRITE (CLL1,'(I3)') IND1-1
            CLFMT = '(I'//CLL1//')'
            READ (CLARG(1:IND1-1),FMT=CLFMT) IOESATID(IOBSER)
            WRITE (CLL1,'(I3)') LEN1-IND1
            CLFMT = '(F'//CLL1//'.0)'
            READ (CLARG(IND1+1:LEN1),FMT=CLFMT) OERR(IOBSER)
            IF (ITEST >= 1) THEN
              WRITE(IU06,*) 'OBSERVATION ERROR WEIGHT FOR SAT. ', &
     &            IOESATID(IOBSER), ' IS ', OERR(IOBSER)
            ENDIF

          ELSEIF ( CLOPTLET == 'd' ) THEN

! DATES ARRAY.
! ------------
            IDATES = IDATES + 1
            IF (IDATES > IMO) THEN
              ALLOCATE(CLMD(IMO))
              CLMD = CLDATES
              DEALLOCATE(CLDATES)
              IMO = IMO + 10
              ALLOCATE(CLDATES(IMO))
              DO JL = 1, IMO - 10
                CLDATES(JL) = CLMD(JL)
              ENDDO
              DEALLOCATE(CLMD)
            ENDIF
            CLDATES(IDATES) = CLARG
          ENDIF
        ENDIF
      ENDDO CMDLINE

      IF (ITEST >= 1) THEN
        WRITE(*,'("NUMBER OF FILES TO PROCESS: ",I2)') IDATES
        WRITE(*,'(5a13)') (CLDATES(JL),JL=1,IDATES)
        IF (LAVERAGE) THEN
          WRITE(*,*) '(AVERAGED) SUPER-OBSERVATIONS ARE USED.'
        ELSE
          WRITE(*,*) 'INDIVIDUAL OBSERVATIONS ARE USED.'
        ENDIF
      ENDIF

      IF (IDATES == 0 ) THEN
        WRITE(*,*) 'NOTHING TO DO?'
        WRITE(*,*) 'YOU MAY HAVE GIVEN THE WRONG COMMAND LINE OPTIONS:'
        WRITE(*,*) 'USAGE: rfl4wam -d date [ -d date [ ... ] ] ', &
     &             '[ -i itest ] [ -t type ] [ -a idatawl ] ' &
     &           , '[ -e satID:weight [ -e ... ] ]'
        WRITE(*,*) '       type is either: i  for individual (default)'
        WRITE(*,*) '                   or: a  for averaged super-obs.'
        CALL EXIT (1)
      ENDIF

! ---------------------------------------------------------------------

!     2. GLOBAL SETTINGS.
!        ----------------

! Initialize ODB-2
!
#ifdef WITH_ODB
      CALL ODB_START()
      odb_handle  = odb_write_new(config, cerr)
#endif

      Z4=1._JWRB
      NPRECR = KIND(Z4)
      I4=1
      NPRECI = KIND(I4)

      IRANK = 1
      NPROC = 1

      NPREVIOUS=IRANK-1
      IF (IRANK == NPROC) THEN
        NNEXT=0
      ELSE
        NNEXT=IRANK+1
      ENDIF

      LWCOU = .FALSE.
      NLON=1
      NLAT=1
      LLRNL=.TRUE.

      CALL WVWAMINIT (LWCOU,IU06,LLRNL,NGAUSSW,NLON,NLAT,RSOUTW,RNORTW)


      NPR=NPROC
      LL1D = .TRUE.
      LLIRANK = .FALSE.
      LLWVENVI = .FALSE.

      CALL MPDECOMP(NPR, MAXLEN, LLIRANK, LLWVENVI)

      MPMAXLENGTH=MAXLEN

! SET IREAD TO THE PROCESSOR WHICH WILL READ THE INPUT FILES.
! -----------------------------------------------------------
      IREAD=1

! ----------------------------------------------------------------------


!*    5. PRINT INITIAL CONDITIONS AS READ FROM PREPROCESSING.
!        ----------------------------------------------------

      IF (ITEST >= 1) THEN
        WRITE(IU06,*) '  '
        WRITE(IU06,*) ' WAVE MODEL GRID ORGANISATION:'
        WRITE(IU06,3002) ' SOUTHERNMOST LATITUDE IN GRID IS .......: ', &
     &                 AMOSOP, ' DEGREE'
        WRITE(IU06,3002) ' NORTHERNMOST LATITUDE IN GRID IS .......: ', &
     &                 AMONOP, ' DEGREE'
        WRITE(IU06,3002) ' WESTERNMOST LONGITUDE IN GRID IS .......: ', &
     &               AMOWEP, ' DEGREE'
        WRITE(IU06,3002) ' EASTERNMOST LONGITUDE IN GRID IS .......: ', &
     &               AMOEAP, ' DEGREE'
        WRITE(IU06,3002) ' LATITUDE INCREMENT IS ..................: ', &
     &                 XDELLA, ' DEGREE'
        WRITE(IU06,3002) ' LONGITUDE INCREMENT IS .................: ', &
     &               XDELLO, ' DEGREE'
        WRITE(IU06,*) '  '
        WRITE(IU06,*) '  '
 3002   FORMAT(3x,a,f8.2,a)
 3003   FORMAT(3x,a,i6,  a)
      ENDIF

#ifdef WITH_ODB
      iseqno_offset = 0
      CALL GET_ENVIRONMENT_VARIABLE('ODB_SEQNO_OFFSET', cenv_tmp)
      if (cenv_tmp /= '') then
        read(cenv_tmp,*) iseqno_offset
      endif
#endif

      DO JLD = 1, IDATES
        CDTPRO(1:12) = CLDATES(JLD)
        CDTPRO(13:14) = '00'
        io = IWAM_GET_UNIT(IU06, 'rfl'//CLDATES(JLD), 'w', 'u', 0, 'READWRITE')

        IF (IOBSER > 0) THEN
          CALL GRDATA (NOBS, IDATAWL, LAVERAGE, NIJALT, IOBSERRSATID=IOESATID(1:IOBSER), OBSERR=OERR(1:IOBSER))
        ELSE
          CALL GRDATA (NOBS, IDATAWL, LAVERAGE, NIJALT)
        ENDIF

        WRITE(io) CDTPRO, NOBS

        IF (NOBS > 0) THEN
          WRITE (io) IJALT 
          WRITE (io) ALTDATA
          IF (LLUNSTR) WRITE (io) ALTUNDATA
        ENDIF
        CLOSE (io)

!       TRANSFER THE CONTENT OF THE ALTIMETER RFL FILED TO ODB
        LLCLOSE = (JLD == IDATES)
#ifdef WITH_ODB
        CALL RFL2ODB(NOBS, CLDATES(JLD), LLCLOSE) 
#endif

        IF (NOBS > 0) THEN
          DEALLOCATE(IJALT)
          DEALLOCATE(ALTDATA)
          DEALLOCATE(ALTEXDATA)
          DEALLOCATE(CDATEOBS)
        ENDIF

      ENDDO

      IF (LHOOK) CALL DR_HOOK('RFL4WAM',1,ZHOOK_HANDLE)

      CALL MPL_END()

! end WITH_WAMASSI
#endif

END PROGRAM rfl4wam
