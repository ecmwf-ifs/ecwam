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
USE YOWGRIB_HANDLES , ONLY : NGRIB_HANDLE_WAM_I
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
#include "preset_wgrib_template.intfb.h"
#include "wgribenout.intfb.h"

INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL


INTEGER(KIND=JWIM), PARAMETER :: NBITSPERVALUE = 32  !! need a higher precision to code the grid point index 

INTEGER(KIND=JWIM) :: IJ, IRMPIRNK, IRGRDPT, IFLAG
INTEGER(KIND=JWIM) :: LFILE, IUOUT, ICOUNT
INTEGER(KIND=JWIM) :: IFCST, IT, IPARAM, ITABLE, IZLEV

REAL(KIND=JWRB) :: ZHOOK_HANDLE

CHARACTER(LEN= 2) :: MARSTYPEBAK
CHARACTER(LEN=14) :: CDATE
CHARACTER(LEN=296) :: OUTFILE

LOGICAL :: LFDBBAK, LLCREATE
LOGICAL :: FFLAGBAK(JPPFLAG), GFLAGBAK(JPPFLAG)

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('OUTMDLDCP',0,ZHOOK_HANDLE)

! Output will be save in 
IF ( ICPLEN > 0 ) THEN
  OUTFILE = CPATH(1:ICPLEN)//'/wam_model_mpi_decomposition.grb'
ELSE
  OUTFILE = 'wam_model_mpi_decomposition.grb'
ENDIF
LFILE=LEN_TRIM(OUTFILE)

WRITE(IU06,*) ' OUTMDLDCP : The MPI decomposition is written to ', OUTFILE(1:LFILE)

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

! GRID POINT INDEX:
IRGRDPT = JPPFLAG-3
GFLAG(IRGRDPT) = .TRUE.

! Set output parameter mapping (and allocate BOUT) 
CALL MPCRTBL


! Defining the ouput fields:
! -------------------------
BOUT(IJS:IJL, IRMPIRNK) = REAL(IRANK, JWRB)
DO IJ = IJS, IJL
  BOUT(IJ, IRGRDPT) = REAL(IJ, JWRB)
ENDDO


! Gather data for output (to IRANK = 1)
CALL OUTGRID

IF(IRANK == 1) THEN
  ! Grib output to file:
  CALL IGRIB_OPEN_FILE(IUOUT,OUTFILE(1:LFILE),'w')

  ! Prepare grib template
  LLCREATE = .TRUE.
write(0,*) 'debile before PRESET_WGRIB_TEMPLATE '
  CALL PRESET_WGRIB_TEMPLATE("I", NGRIB_HANDLE_WAM_I, LLCREATE=LLCREATE, NBITSPERVALUE=NBITSPERVALUE )
write(0,*) 'debile after PRESET_WGRIB_TEMPLATE '


  ! keep looping over all posible output varaibles (as in outint)
  ICOUNT=0
  DO IFLAG = 1, JPPFLAG
    IF ( GFLAG(IFLAG) ) THEN
      ICOUNT = ICOUNT+1

      IT = ITOBOUT(IFLAG)
      ITABLE = INFOBOUT(IT,1)
      IPARAM = INFOBOUT(IT,2)
      IZLEV  = INFOBOUT(IT,3)

      CDATE = CDATEA  ! set date to start of the run

      CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, GOUT(ICOUNT,:,:),  &
 &                    ITABLE, IPARAM, IZLEV, 0 , 0,             &
 &                    CDATE, IFCST, MARSTYPE, LFDB, IUOUT)
    ENDIF
  ENDDO

  CALL IGRIB_CLOSE_FILE(IUOUT)
  CALL IGRIB_RELEASE(NGRIB_HANDLE_WAM_I)

  IF(ALLOCATED(GOUT)) DEALLOCATE(GOUT)

ENDIF


! Restore output configuration: 
LFDB = LFDBBAK
FFLAG(:) = FFLAGBAK(:)
GFLAG(:) = GFLAGBAK(:)
MARSTYPE = MARSTYPEBAK
!  reset output field mapping
CALL MPCRTBL

WRITE(IU06,*) ' '

IF (LHOOK) CALL DR_HOOK('OUTMDLDCP',1,ZHOOK_HANDLE)

END SUBROUTINE OUTMDLDCP
