SUBROUTINE RFL2ODB(NOBS, CLDATES, LLCLOSE)

!--------------------------------------------------------------------

!**** *RFL2ODB* - TRANSFERS THE CONTENT OF THE ALTIMETER RFL FILED TO ODB

! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWALTAS , ONLY : IJALT, ALTDATA, ALTEXDATA, CDATEOBS
USE YOWTEST  , ONLY : IU06

USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK

#ifdef WITH_ODB
USE PARKIND1  ,ONLY : JPIM, JPRB, JPRD
use, intrinsic :: iso_c_binding
use odc_c_binding
use odb2_flag_definitions
use odbmap_reportype, only : find_obstype_codetype, find_reportype
#endif

! -------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JWIM), INTENT(IN) :: NOBS
CHARACTER(LEN=12), INTENT(IN) :: CLDATES
LOGICAL, INTENT(IN) :: LLCLOSE 

#ifdef WITH_ODB

INTEGER(KIND=JWIM) :: IOBS
INTEGER(KIND=JWIM), DIMENSION(999) :: ICODE
INTEGER(KIND=JWIM), DIMENSION(999) :: ISENSORNB

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JWRU) :: Z8

LOGICAL, SAVE :: LFRSTIME

DATA LFRSTIME / .TRUE. /
DATA ICODE /999*999/ ! 999 is for Missing
DATA ISENSORNB /999*0/ ! 0 is for Missing

!     ODB-2 definitions
TYPE(C_PTR), SAVE                  :: odb_handle
INTEGER(KIND=JWIM), SAVE           :: iseqno_offset

TYPE(C_PTR), SAVE                  :: odb_it
INTEGER(KIND=C_INT)                :: cerr
CHARACTER(KIND=C_CHAR, LEN=256)    :: config = C_NULL_CHAR
CHARACTER(KIND=C_CHAR, LEN=256)    :: odb_file
INTEGER(KIND=JWIM), PARAMETER      :: mdi = 2147483647
CHARACTER(LEN=256)                 :: one_colname, ctype
CHARACTER(LEN=256)                 :: bitfield_names, bitfield_sizes
INTEGER(KIND=C_INT), PARAMETER     :: c_ncolumns=32
REAL(KIND=C_DOUBLE)                :: odb_values(c_ncolumns), odb_values_b(c_ncolumns)
INTEGER(KIND=JWIM)                 :: idate, itime
CHARACTER(LEN=256)                 :: cenv_tmp
INTEGER(KIND=JWIM)                 :: igroupid, isatid
INTEGER(KIND=JWIM)                 :: ibufrtype,isubtype
INTEGER(KIND=JWIM)                 :: iobstype,icodetype
INTEGER(KIND=JWIM)                 :: isensor,ireportype,isatinst
INTEGER(KIND=JWIM)                 :: NumVar
LOGICAL                            :: lnot_active, lnot_active_b

DATA iseqno_offset /0/

#include "wam_fcobs.h"

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('RFL2ODB',0,ZHOOK_HANDLE)


! Initialize ODB-2
!
IF (LFRSTIME) THEN
  CALL ODB_START()
  odb_handle  = odb_write_new(config, cerr)

  iseqno_offset = 0
  CALL GET_ENVIRONMENT_VARIABLE('ODB_SEQNO_OFFSET', cenv_tmp)
  IF (cenv_tmp /= '') THEN
    READ(cenv_tmp,*) iseqno_offset
  ENDIF
  LFRSTIME=.FALSE.
ENDIF

!! CREATE ASCII FILES TO NEEDED TO GENERATE THE RALT ODB (USING THE 
!!    UTILITY  simulobs2odb )

        ! TO CONVERT FROM SAT. ID TO BUFR SUBTYPE/ODB CODETYPE
        ! SPECIFY THE SENSOR NUMBER FOR ODB
        ! NOT THE BEST WAY TO DO THIS, BUT THE SIMPLEST
        ! NOTE THAT ICODE() AS DEFINED BELOW, REFERS TO SUBTYPE (NOT CODETYPE)
        !           CODETYPE IS 123 FOR ALL ALTIMETER DATA

        ! ERS-1
        ICODE(1)=123
        ISENSORNB(1)=64

        ! ERS-2
        ICODE(2)=123
        ISENSORNB(2)=64

        ! Cryosat-2
        ICODE(47)=217
        ISENSORNB(47)=177

        ! ENVISAT
        ICODE(60)=213
        ISENSORNB(60)=147

        ! Sentinel-3A
        ICODE(61)=220
        ISENSORNB(61)=178

        ! Sentinel-3B
        ICODE(65)=220
        ISENSORNB(65)=178

        ! Sentinel-6A
        ICODE(66)=233
        ISENSORNB(66)=57      ! Code is defined in ODB but not the sensor! (Saleh)

        ! Sentinel-6B
        ICODE(67)=233
        ISENSORNB(67)=57      ! Code is defined in ODB but not the sensor! (Saleh)

        ! Jason-1
        ICODE(260)=214
        ISENSORNB(260)=9

        ! Jason-2
        ICODE(261)=216
        ISENSORNB(261)=367

        ! Jason-3
        ICODE(262)=216
        ISENSORNB(262)=367

        ! SARAL/AltiKa
        ICODE(441)=218
        ISENSORNB(441)=368

        ! CFOSat
        ICODE(802)=220        ! We need to define this in ODB (Saleh)
        ISENSORNB(802)=368    ! sensor number is not correct for SWIM, WMO needs to provide the correct one (Saleh)

        NumVar=2

        odb_file = 'wam_ralt'//CLDATES//'.odb'//C_NULL_CHAR

        odb_it = odb_write_iterator_new(odb_handle, odb_file, cerr)

! Creation of ODB-2 header
!!!!!!!! it seems that the order of the column is essential when it comes to convert
!!!!!!!! this odb2 file into odb1 (see odb/tools/Odb2_to_Odb1_ralt.F90)
!!!!!!!! all hdr entries must come first, followed by body, errsat and then sat. 
        cerr = odb_write_set_no_of_columns(odb_it, c_ncolumns)
! Seqno
        one_colname='seqno@hdr'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 0, ODB_INTEGER, one_colname)
! obstype
        one_colname='obstype@hdr'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 1, ODB_INTEGER, one_colname)
! codetype
        one_colname='codetype@hdr'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 2, ODB_INTEGER, one_colname)

! bufrtype
        one_colname='bufrtype@hdr'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 3, ODB_INTEGER, one_colname)
! subtype
        one_colname='subtype@hdr'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 4, ODB_INTEGER, one_colname)

! groupid
        one_colname='groupid@hdr'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 5, ODB_INTEGER, one_colname) 
! reportype
        one_colname='reportype@hdr'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 6, ODB_INTEGER, one_colname)

! date
        one_colname='date@hdr'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 7, ODB_INTEGER, one_colname) 
! time
        one_colname='time@hdr'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 8, ODB_INTEGER, one_colname) 
! lat
        one_colname='lat@hdr'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 9, ODB_REAL, one_colname) 
! lon
        one_colname='lon@hdr'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it,10, ODB_REAL, one_colname) 
! gp_number
        one_colname='gp_number@hdr'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it,11, ODB_INTEGER, one_colname) 
! sensor 
        one_colname='sensor@hdr'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it,12, ODB_INTEGER, one_colname) 
! report_status
        ctype='status_t'
        call get_odb_flags(ctype,bitfield_names, bitfield_sizes)
        bitfield_names=trim(bitfield_names)//C_NULL_CHAR
        bitfield_sizes=trim(bitfield_sizes)//C_NULL_CHAR
        one_colname='report_status@hdr'//C_NULL_CHAR
        cerr = odb_write_set_bitfield(odb_it, 13, ODB_INTEGER, one_colname, bitfield_names, bitfield_sizes)
! report_event1
        ctype='report_event1_t'
        call get_odb_flags(ctype,bitfield_names, bitfield_sizes)
        bitfield_names=trim(bitfield_names)//C_NULL_CHAR
        bitfield_sizes=trim(bitfield_sizes)//C_NULL_CHAR
        one_colname='report_event1@hdr'//C_NULL_CHAR
        cerr = odb_write_set_bitfield(odb_it, 14, ODB_INTEGER, one_colname, bitfield_names, bitfield_sizes)
! report_event2
        one_colname='report_event2@hdr'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 15, ODB_INTEGER, one_colname)
! distribtype
        one_colname='distribtype@hdr'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 16, ODB_INTEGER, one_colname) 
! body.len
        one_colname='LINKLEN(body)@hdr'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it,17, ODB_INTEGER, one_colname) 
! sat.len
        one_colname='LINKLEN(sat)@hdr'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it,18, ODB_INTEGER, one_colname) 
! errstat.len
        one_colname='LINKLEN(errstat)@hdr'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 19, ODB_INTEGER, one_colname) 

! entryno
        one_colname='entryno@body'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 20, ODB_INTEGER, one_colname) 
! varno
        one_colname='varno@body'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 21, ODB_INTEGER, one_colname) 
! obsvalue
        one_colname='obsvalue@body'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 22, ODB_REAL, one_colname) 
! biascorr
        one_colname='biascorr@body'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 23, ODB_REAL,  one_colname)
! fg_depar
        one_colname='fg_depar@body'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 24, ODB_REAL, one_colname) 
! an_depar
        one_colname='an_depar@body'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 25, ODB_REAL, one_colname) 

! vertco_reference_1
        one_colname='vertco_reference_1@body'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 26, ODB_REAL, one_colname) 

! datum_status
        ctype='status_t'
        call get_odb_flags(ctype,bitfield_names, bitfield_sizes)
        bitfield_names=trim(bitfield_names)//C_NULL_CHAR
        bitfield_sizes=trim(bitfield_sizes)//C_NULL_CHAR
        one_colname='datum_status@body'//C_NULL_CHAR
        cerr = odb_write_set_bitfield(odb_it, 27, ODB_INTEGER, one_colname, bitfield_names, bitfield_sizes)
! datum_event1
        ctype='datum_event1_t'
        call get_odb_flags(ctype,bitfield_names, bitfield_sizes)
        bitfield_names=trim(bitfield_names)//C_NULL_CHAR
        bitfield_sizes=trim(bitfield_sizes)//C_NULL_CHAR
        one_colname='datum_event1@body'//C_NULL_CHAR
        cerr = odb_write_set_bitfield(odb_it, 28, ODB_INTEGER, one_colname, bitfield_names, bitfield_sizes)
! datum_event2
        one_colname='datum_event2@body'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 29, ODB_INTEGER, one_colname) 
! obs_error
        one_colname='obs_error@errstat'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 30, ODB_REAL, one_colname) 
! satellite_identifier
        one_colname='satellite_identifier@sat'//C_NULL_CHAR
        cerr = odb_write_set_column(odb_it, 31, ODB_INTEGER, one_colname) 

        cerr = odb_write_header(odb_it)
! End ODB-2 header

        IF (NOBS > 0) THEN
          IF (NOBS > 9999999) THEN
            WRITE(IU06,*)'WARNING: FORMAT I7 IS NOT ENOUGH IN rfl2odb.F'
            WRITE(IU06,*)'         FIRST  9999999  ARE OUTPUT'
          ENDIF

          DO IOBS=1, NOBS

!           Flag missing data
            IF (ALTDATA(IOBS,1) > mdi .OR. IJALT(IOBS,3) == -1 ) THEN
              ALTDATA(IOBS,1) = -mdi
              ALTDATA(IOBS,3) = -mdi
              IJALT(IOBS,3) = -1
            ENDIF
            IF (ALTEXDATA(IOBS,3) > mdi .OR. IJALT(IOBS,4) == -1) THEN
              ALTEXDATA(IOBS,3) = -mdi
              IJALT(IOBS,4) = -1
            ENDIF

            read(CDATEOBS(IOBS)(1:8), '(i8)') idate
            read(CDATEOBS(IOBS)(9:14), '(i6)') itime
! Columns that are common to the 2 varno
            iobstype = 12
            icodetype = 123   ! RALT data always have code 123
            ibufrtype = 12
            isubtype = ICODE(IJALT(IOBS,2)) 
            if (isubtype  == 999 ) then
               WRITE(IU06,*)'WARNING IN rfl2odb: subtype ',IJALT(IOBS,2), 'does not have an odb equivalent !'
               CALL EXIT (1)
            endif
            isatid =  IJALT(IOBS,2) 
            isensor = ISENSORNB(IJALT(IOBS,2)) ! No sensor needed for RALT data
            isatinst= ISENSORNB(IJALT(IOBS,2)) ! No need of satellite instrument for RALT data 
            odb_values(1) = IOBS + iseqno_offset ! seqno@hdr
            odb_values(2) = iobstype             ! obstype@hdr
            odb_values(3) = icodetype            ! codetype@hdr
            odb_values(4) = ibufrtype            ! bufrtype@hdr
            odb_values(5) = isubtype             ! subtype@hdr

            igroupid = mdi
            ireportype = mdi

            ! If there is an error like "Error: No reportype mapping", then make sure that there is a report type for the same satellite/sensor exists in file
            !    ${XDATA}/${IFS_CYCLE}/an/odb_code_mappings.dat  (e.g. /home/rdx/data/48r1/an/odb_code_mappings.dat)
            call find_reportype(igroupid, isatid, ibufrtype,isubtype, &
               & iobstype,icodetype,isensor,ireportype,isatinst)

            odb_values(6) = igroupid             ! groupid@hdr
            odb_values(7) = ireportype           ! reportype@hdr
            odb_values(8) = idate                ! date@hdr
            odb_values(9) = itime                ! time@hdr
            odb_values(10) = ALTEXDATA(IOBS,1)   ! lon@hdr
            odb_values(11) = ALTEXDATA(IOBS,2)   ! lon@hdr
            odb_values(12) = IJALT(IOBS,1)       ! gp_number@hdr WAM block index
            odb_values(13) = isensor             ! sensor@hdr

            Z8 = 0.0_JWRU
            odb_values(14) =  ZCHSTAT_ACTIVE(Z8) ! report_status@hdr
            odb_values(15) =  0                  ! report_event1@hdr
            odb_values(16) =  0                  ! report_event2@hdr
            odb_values(17) =  0                  ! distribtype@hdr ! Needed for CCMA
                                                 ! 0 to preserve gp_number@hdr (2 leads to overwriting it by IFS)

            odb_values(18) = NumVar              ! body.len@hdr (=2 because 2 varno)
            odb_values(19) = 1                   ! sat.len@hdr
            odb_values(20) = NumVar              ! errstat.len@hdr (=body.len)

            odb_values_b(:) = odb_values(:)

! columns depending on the varno Var. No (220=SWH)
            odb_values(21) = 1                   ! entryno@body
            odb_values(22) = 220                 ! varno@body
            odb_values(23) = ALTDATA(IOBS,1)     ! obsvalue@body
            odb_values(24) = ALTDATA(IOBS,1)-ALTDATA(IOBS,3)   ! biascorr@body
            odb_values(25) = -mdi                 ! fg_depar@body
            odb_values(26) = -mdi                 ! an_depar@body 
            odb_values(27) = 0.0_JWRB             ! vertco_reference_1@body

            Z8 = 0.0_JWRU
            lnot_active = .true.
            IF (IJALT(IOBS,3) > 0) THEN
              odb_values(28) = ZCHSTAT_ACTIVE(Z8)   ! datum_status@body
              lnot_active = .false.
            ELSEIF (IJALT(IOBS,3) == 0) THEN
              odb_values(28) = ZCHSTAT_BLACKLST(Z8) ! datum_status@body
            ELSE
              odb_values(28) = ZCHSTAT_REJECT(Z8)   ! datum_status@body
            ENDIF

            odb_values(29) = 0                   ! datum_event1@body

            odb_values(30) =  IJALT(IOBS,3)      ! datum_event2@body


            IF ( ALTDATA(IOBS,2) <= 0.0_JWRB ) THEN
              odb_values(31) = -mdi              ! obs_error@errstat
            ELSE
              odb_values(31) = ALTDATA(IOBS,2)   ! obs_error@errstat
            ENDIF

            odb_values(32) =  IJALT(IOBS,2)      ! satellite_identifier@sat


! Var. No (221=Wnd Spd)
            odb_values_b(21) = 2                   ! entryno@body
            odb_values_b(22) = 221                 ! varno@body
            odb_values_b(23) = ALTEXDATA(IOBS,3)   ! obsvalue@body
            odb_values_b(24) = 0.0_JWRB            ! biascorr@body
            odb_values_b(25) = -mdi                ! fg_depar@body
            odb_values_b(26) = -mdi                ! an_depar@body
            odb_values_b(27) = 0.0_JWRB            ! vertco_reference_1@body

            Z8 = 0.0_JWRU
            lnot_active_b = .true.
            IF (IJALT(IOBS,4) > 0) THEN
              odb_values_b(28) = ZCHSTAT_ACTIVE(Z8)  ! datum_status@body
              lnot_active_b = .false.
            ELSEIF (IJALT(IOBS,4) == 0) THEN
              odb_values_b(28) = ZCHSTAT_PASSIVE(Z8) ! datum_status@body
            ELSE
              odb_values_b(28) = ZCHSTAT_REJECT(Z8)  ! datum_status@body
            ENDIF

            odb_values_b(29) = 0                   ! datum_event1@body

            odb_values_b(30) =  IJALT(IOBS,4)      ! datum_event2@body

            odb_values_b(31) = -mdi                ! obs_error@errstat

            odb_values_b(32) =  IJALT(IOBS,2)      ! satellite_identifier@sat


            ! Adjust report_status if all datum_status are inactive
            IF ( (lnot_active).and.(lnot_active_b))THEN
                odb_values(14) = 0.0_JWRB
                odb_values_b(14) = 0.0_JWRB
            ENDIF

            ! Write both rows 
            cerr = odb_write_set_next_row(odb_it,odb_values,c_ncolumns)
            cerr = odb_write_set_next_row(odb_it,odb_values_b,c_ncolumns)

          ENDDO

        ENDIF

        cerr = odb_write_iterator_delete(odb_it)
        iseqno_offset = NOBS

IF (LLCLOSE) cerr = odb_write_delete(odb_handle)

IF (LHOOK) CALL DR_HOOK('RFL2ODB',1,ZHOOK_HANDLE)

#endif
END SUBROUTINE RFL2ODB
