! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE WRITENEMOCPLFLDS( NEMO2WAM, KFIELDS )

! ----------------------------------------------------------------------

!**** *WRITENEMOCPLFLDS* -
!     ----------

!  OUTPUT COUPLING FIELDS EITHER SENT OR RECEIVED FROM NEMO

!  THE GRIB DATA WILL BE SAVED TO FILE (see OUTFILE)

!  It will temporarilly reset the output set up to only output what is necessary
!  All is reset at the end.

! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
USE YOWDRVTYPE  , ONLY : OCEAN2WAVE
USE YOWCOUT  , ONLY : JPPFLAG, FFLAG, GFLAG, NFLAG, NIPRMOUT, LFDB
USE YOWGRIB_HANDLES , ONLY : NGRIB_HANDLE_WAM_I
USE YOWGRIBHD, ONLY : PPMISS, PPEPS, PPREC, PPRESOL, PPMIN_RESET, NTENCODE, NGRBRESS, LNEWLVTP
USE YOWGRID  , ONLY : IJS, IJL, NPROMA_WAM, NCHNK, KIJL4CHNK, IJFROMCHNK, IPRMFROMIJ
USE YOWINTP  , ONLY : GOUT
USE YOWMAP   , ONLY : IRGG, AMONOP, AMOSOP, XDELLA, NLONRGG, BLK2GLO, NIBLO, NGX, NGY, CLDOMAIN 
USE YOWMPP   , ONLY : IRANK, NPRECI 
USE YOWPARAM , ONLY : LL1D
USE YOWPCONS , ONLY : ZMISS
USE YOWSPEC  , ONLY : IJ2NEWIJ
USE YOWTEST  , ONLY : IU06, ITEST
USE YOWTEXT  , ONLY : ICPLEN, CPATH
USE YOWSTAT  , ONLY : CDATEF, CDTPRO, IDELPRO


USE YOMHOOK  , ONLY : LHOOK ,DR_HOOK, JPHOOK
USE YOWGRIB  , ONLY : JPKSIZE_T, &
                    & IGRIB_OPEN_FILE, &
                    & IGRIB_CLOSE_FILE, &
                    & IGRIB_WRITE_BYTES, &
                    & IGRIB_GET_MESSAGE, &
                    & IGRIB_GET_MESSAGE_SIZE,&
                    & IGRIB_CLONE, &
                    & IGRIB_RELEASE

! ----------------------------------------------------------------------

IMPLICIT NONE

#include "abort1.intfb.h"
#include "mpcrtbl.intfb.h"
#include "outgrid.intfb.h"
#include "preset_wgrib_template.intfb.h"
#include "wgribencode.intfb.h"
#include "difdate.intfb.h"

TYPE(OCEAN2WAVE), INTENT(IN) :: NEMO2WAM
INTEGER, INTENT(IN) :: KFIELDS

INTEGER(KIND=JWIM), PARAMETER :: NBITSPERVALUE = 24  !! need a higher precision to code the grid point index 

INTEGER(KIND=JWIM) :: IJ, IPRM, IJGLO, ITMIN, ITMAX, IK, IM, IR, IFLAG, IDT, KSTEP, IX, JSN, I,J,K ,IK1, IK2
INTEGER(KIND=JWIM) :: LFILE, IUOUT, ICOUNT
INTEGER(KIND=JWIM) :: IFCST, IT, IPARAM, ITABLE, IZLEV
INTEGER(KIND=JWIM) :: ICHNK, KIJS, KIJL, IJSB, IJLB
INTEGER(KIND=JWIM) :: IGRIB_HANDLE
INTEGER(KIND=JWIM) :: ISIZE
INTEGER(KIND=JPKSIZE_T) :: KBYTES
INTEGER(KIND=JWIM), ALLOCATABLE :: KGRIB_BUFR(:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

CHARACTER(LEN= 2) :: MARSTYPE_DUM
CHARACTER(LEN=14) :: CDATE
CHARACTER(LEN= 8) :: CDI8
CHARACTER(LEN=296) :: OUTFILE

LOGICAL :: LLCREATE, LFDBBAK
LOGICAL :: LGRHDIFS_DUM, LPADPOLES_DUM, LRSTST0_DUM
LOGICAL, DIMENSION(JPPFLAG) :: FFLAGBAK, GFLAGBAK, NFLAGBAK
! TEMPORARY ARRAYS FOR LOCAL MASK AND GLOBAL INDICES
INTEGER(KIND=JWIM), DIMENSION(NGX,NGY) :: IGLOBAL
INTEGER(KIND=JWIM), DIMENSION(IJS:IJL) :: NLOCMSK, NGLOIND
REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:, :, :) :: BOUT

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WRITENEMOCPLFLDS',0,ZHOOK_HANDLE)

CALL DIFDATE(CDATEF, CDTPRO, IFCST)
CDATE = CDATEF

! Output will be save in
WRITE(CDI8,'(I8.8)') IFCST/IDELPRO
OUTFILE='wam_nemorecvfields_'//CDI8//'.grib'
IF ( ICPLEN > 0 ) THEN
  OUTFILE = CPATH(1:ICPLEN)//'/'//TRIM(OUTFILE)
ENDIF
LFILE=LEN_TRIM(OUTFILE)

WRITE(IU06,*) ' WRITENEMOCPLFLDS : The coupling fields are written to ', OUTFILE(1:LFILE)

! save output parameter selection (it will be overwritten and then reset)
FFLAGBAK(:) = FFLAG(:)
GFLAGBAK(:) = GFLAG(:)
NFLAGBAK(:) = NFLAG(:)
LFDBBAK = LFDB

! Create a new output parameter selection that will only output
! parameters relevant for the description of the model decomposition
FFLAG(:) = .FALSE.
GFLAG(:) = .FALSE.
NFLAG(:) = .FALSE.
LFDB = .FALSE. 

GFLAG(JPPFLAG-4) = .TRUE.
GFLAG(JPPFLAG-3) = .TRUE.
GFLAG(JPPFLAG-2) = .TRUE.
GFLAG(JPPFLAG-1) = .TRUE.
GFLAG(JPPFLAG) = .TRUE.


! Set output parameter mapping
CALL MPCRTBL

ALLOCATE ( BOUT(NPROMA_WAM, MAX(NIPRMOUT,1), NCHNK) )

! Defining the ouput fields:
! -------------------------

IK1 = 0
IK2 = 0

DO ICHNK = 1, NCHNK

   KIJS = 1
   IJSB = IJFROMCHNK(KIJS,ICHNK)
   KIJL = KIJL4CHNK(ICHNK)
   IJLB = IJFROMCHNK(KIJL,ICHNK) 

  IR = 0

  ! SSTL
  IR = IR + 1
  BOUT(KIJS:KIJL, IR, ICHNK) = NEMO2WAM%NEMOSST(1:KIJL4CHNK(ICHNK),ICHNK)

  ! CICOVER:
  IR = IR + 1
  BOUT(KIJS:KIJL, IR, ICHNK) = NEMO2WAM%NEMOSST(1:KIJL4CHNK(ICHNK),ICHNK)

  ! CICHICK:
  IR = IR + 1
  BOUT(KIJS:KIJL, IR, ICHNK) = NEMO2WAM%NEMOCITHICK(1:KIJL4CHNK(ICHNK),ICHNK)

  ! UCUR:
  IR = IR + 1
  BOUT(KIJS:KIJL, IR, ICHNK) = NEMO2WAM%NEMOUCUR(1:KIJL4CHNK(ICHNK),ICHNK)

  ! VCUR:
  IR = IR + 1
  BOUT(KIJS:KIJL, IR, ICHNK) = NEMO2WAM%NEMOVCUR(1:KIJL4CHNK(ICHNK),ICHNK)


  IF ( IR /= NIPRMOUT ) THEN
    WRITE(IU06,*) '******************************************'
    WRITE(IU06,*) ' WRITENEMOCPLFLDS : ABORTING' 
    WRITE(IU06,*) ' IR /= NIPRMOUT ', IR, NIPRMOUT 
    WRITE(IU06,*) '******************************************'
    CALL ABORT1
  ENDIF

ENDDO


! Gather data for output (to IRANK = 1)
CALL OUTGRID(BOUT)

IF(IRANK == 1) THEN
  ! Grib output to file:
  CALL IGRIB_OPEN_FILE(IUOUT, OUTFILE(1:LFILE), 'w')

  ! Prepare grib template
  LLCREATE = .TRUE.
  CALL PRESET_WGRIB_TEMPLATE("I", NGRIB_HANDLE_WAM_I, LLCREATE=LLCREATE, NBITSPERVALUE=NBITSPERVALUE )


  ! keep looping over all posible output variables (as in outint)
  ICOUNT=0
  DO IFLAG = JPPFLAG-4, JPPFLAG

    ITABLE=140

    IF ( GFLAG(IFLAG) ) THEN

      ICOUNT = ICOUNT+1

      SELECT CASE (IFLAG)
        CASE (JPPFLAG-4)
          ITABLE=151
          IPARAM=159
          IZLEV=0
        CASE (JPPFLAG-3)
          ITABLE=3
          IPARAM=91
          IZLEV=0
        CASE (JPPFLAG-2)
          ITABLE=3
          IPARAM=92
          IZLEV=0
        CASE (JPPFLAG-1)
          ITABLE=3
          IPARAM=49
          IZLEV=0
        CASE (JPPFLAG)
          ITABLE=3
          IPARAM=50
          IZLEV=0
      END SELECT

      ITMIN = 0
      ITMAX = 0
      IK = 0
      IM = 0

      MARSTYPE_DUM = 'an'

      LGRHDIFS_DUM = .FALSE.
      IDT = 0
      KSTEP = 0
      LPADPOLES_DUM = .FALSE.
      LRSTST0_DUM = .FALSE.

      IGRIB_HANDLE=-99
      CALL IGRIB_CLONE(NGRIB_HANDLE_WAM_I, IGRIB_HANDLE)

      ! grib coding
      CALL WGRIBENCODE( IU06, ITEST, &
 &                      NGX, NGY, &
 &                      GOUT(ICOUNT,:,:),  &
 &                      ITABLE, IPARAM, &
 &                      IZLEV, &
 &                      ITMIN, ITMAX, &
 &                      IK, IM, &
 &                      CDATE, IFCST, MARSTYPE_DUM, &
 &                      PPMISS, PPEPS, PPREC, PPRESOL, PPMIN_RESET, NTENCODE, &
 &                      LGRHDIFS_DUM, &
 &                      IDT, &
 &                      NGRBRESS, LNEWLVTP, LPADPOLES_DUM, &
 &                      SIZE(NLONRGG), NLONRGG, IRGG, &
 &                      AMONOP, AMOSOP, XDELLA, CLDOMAIN, &
 &                      KSTEP, LRSTST0_DUM, &
 &                      ZMISS, &
 &                      IGRIB_HANDLE )

      ! Save grib data to file

      CALL IGRIB_GET_MESSAGE_SIZE(IGRIB_HANDLE, KBYTES)
      ISIZE = (KBYTES+NPRECI-1)/NPRECI
      ALLOCATE(KGRIB_BUFR(ISIZE))

      CALL IGRIB_GET_MESSAGE(IGRIB_HANDLE, KGRIB_BUFR)

      CALL IGRIB_WRITE_BYTES(IUOUT, KGRIB_BUFR, KBYTES)

      CALL IGRIB_RELEASE(IGRIB_HANDLE)

      DEALLOCATE(KGRIB_BUFR)

    ENDIF
  ENDDO

  CALL IGRIB_CLOSE_FILE(IUOUT)
  CALL IGRIB_RELEASE(NGRIB_HANDLE_WAM_I)

  IF (ALLOCATED(GOUT)) DEALLOCATE(GOUT)

ENDIF


DEALLOCATE ( BOUT )

! Restore output configuration: 
FFLAG(:) = FFLAGBAK(:)
GFLAG(:) = GFLAGBAK(:)
NFLAG(:) = NFLAGBAK(:)
LFDB = LFDBBAK
!  reset output field mapping
CALL MPCRTBL

WRITE(IU06,*) ' '

IF (LHOOK) CALL DR_HOOK('WRITENEMOCPLFLDS',1,ZHOOK_HANDLE)

END SUBROUTINE WRITENEMOCPLFLDS
