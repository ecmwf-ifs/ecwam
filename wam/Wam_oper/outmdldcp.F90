SUBROUTINE OUTMDLDCP (IJS, IJL)

! ----------------------------------------------------------------------

!**** *OUTMDLDCP* -
!     ----------

!  OUPUT THE MODEL DECOMPOSITON 
!    MPI RNAK (IRANK) IN GRIB AS EXTRA PARAMETER 80
!    GRID POINT INDEX (IJ) AS EXTRA PARAMETER 81  

!  THE GRIB DATA WILL BE SAVED TO FILE (see OUTFILE)

!  It will temporarilly reset the output set up to only output what is necessary
!  All is reset at the end.

! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWCOUT  , ONLY : LFDB, JPPFLAG, FFLAG, GFLAG, BOUT, ITOBOUT, INFOBOUT
USE YOWINTP  , ONLY : GOUT
USE YOWMPP   , ONLY : IRANK 
USE YOWPARAM , ONLY : NGX, NGY
USE YOWSTAT  , ONLY : MARSTYPE, CDATEA
USE YOWTEST  , ONLY : IU06, ITEST
USE YOWTEXT  , ONLY : ICPLEN, CPATH

USE YOMHOOK  ,ONLY  : LHOOK ,DR_HOOK
USE GRIB_API_INTERFACE

! ----------------------------------------------------------------------

IMPLICIT NONE
#include "mpcrtbl.intfb.h"
#include "outgrid.intfb.h"
#include "wgribenout.intfb.h"

INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL


INTEGER(KIND=JWIM) :: IJ, IRMPIRNK, IRGRDPT, IFLAG
INTEGER(KIND=JWIM) :: LFILE, IUOUT, ICOUNT
INTEGER(KIND=JWIM) :: IFCST, IT, IPARAM, ITABLE, IZLEV

REAL(KIND=JWRB) :: ZHOOK_HANDLE

CHARACTER(LEN= 2) :: MARSTYPEBAK
CHARACTER(LEN=14) :: CDATE
CHARACTER(LEN=296) :: OUTFILE

LOGICAL :: LFDBBAK
LOGICAL :: FFLAGBAK(JPPFLAG), GFLAGBAK(JPPFLAG)

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('OUTMDLDCP',0,ZHOOK_HANDLE)

! Output will be save in 
OUTFILE = CPATH(1:ICPLEN)//'/wam_model_mpi_decomposition.grb'
LFILE=LEN_TRIM(OUTFILE)

WRITE(IU06,*) '  OUTMDLDCP : '
WRITE(IU06,*) '  The MPI decomposition will be written to ', OUTFILE(1:LFILE)

! save output parameter selection (it will be overwritten and then reset)
LFDBBAK = LFDB
FFLAGBAK(:) = FFLAG(:)
GFLAGBAK(:) = GFLAG(:)
MARSTYPEBAK = MARSTYPE

! create a new output parameter selection that will only output
! parameters relevant for the description of the model decomposition
LFDB = .FALSE.  ! data will be written to file
FFLAG(:) = .FALSE.
GFLAG(:) = .FALSE.
MARSTYPE = 'an'
IFCST = 0

! use the extra field codes (currently the last 5 fields of the list of potential output parameters:
! MPI RANK:
IRMPIRNK = JPPFLAG-4
GFLAG(IRMPIRNK) = .TRUE.
BOUT(IJS:IJL, IRMPIRNK) = IRANK

! GRID POINT INDEX:
IRGRDPT = JPPFLAG-3
GFLAG(IRGRDPT) = .TRUE.
DO IJ = IJS, IJL
  BOUT(IJ, IRGRDPT) = IJ 
ENDDO

! Set output parameter mapping 
CALL MPCRTBL


! Gather data for output (to IRANK = 1)
CALL OUTGRID

! Grib output to file:
IF(IRANK == 1) THEN
  CALL IGRIB_OPEN_FILE(IUOUT,OUTFILE(1:LFILE),'w')

  ! keep looping over all posible output varaibles (as in outint)
  ICOUNT=0
  DO IFLAG = 1, JPPFLAG
    IF ( GFLAG(IFLAG) ) THEN
      ICOUNT = ICOUNT+1

      IT = ITOBOUT(IFLAG)
      ITABLE = INFOBOUT(IT,1)
      IPARAM = INFOBOUT(IT,2)
      IZLEV = INFOBOUT(IT,3)

      CDATE = CDATEA  ! set date to start of the run

      CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, GOUT(ICOUNT,:,:),  &
 &                    ITABLE, IPARAM, IZLEV, 0 , 0,             &
 &                    CDATE, IFCST, MARSTYPE, LFDB, IUOUT)
    ENDIF
  ENDDO

  CALL IGRIB_CLOSE_FILE(IUOUT)

ENDIF


! Restore output configuration: 
LFDB = LFDBBAK
FFLAG(:)=FFLAGBAK(:)
GFLAG(:)=GFLAGBAK(:)
MARSTYPE=MARSTYPEBAK
!  reset output field mapping
CALL MPCRTBL

WRITE(IU06,*) ' '
WRITE(IU06,*) '  OUTMDLDCP ALL DONE '

IF (LHOOK) CALL DR_HOOK('OUTMDLDCP',1,ZHOOK_HANDLE)

END SUBROUTINE OUTMDLDCP
