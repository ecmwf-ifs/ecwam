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

USE YOWCOUT  , ONLY : JPPFLAG, FFLAG, GFLAG, NFLAG, NIPRMOUT, BOUT, ITOBOUT, INFOBOUT
USE YOWGRIB_HANDLES , ONLY : NGRIB_HANDLE_WAM_I
USE YOWGRIBHD, ONLY : PPMISS, PPEPS, PPREC, PPRESOL, PPMIN_RESET, NTENCODE, NGRBRESS, LNEWLVTP
USE YOWGRID  , ONLY : NLONRGG
USE YOWINTP  , ONLY : GOUT
USE YOWMAP   , ONLY : IRGG, AMONOP, AMOSOP, XDELLA
USE YOWMPP   , ONLY : IRANK, NPRECI 
USE YOWPARAM , ONLY : NIBLO, NGX, NGY, LL1D, CLDOMAIN
USE YOWPCONS , ONLY : ZMISS
USE YOWSPEC  , ONLY : IJ2NEWIJ
USE YOWTEST  , ONLY : IU06, ITEST
USE YOWTEXT  , ONLY : ICPLEN, CPATH



USE YOMHOOK  ,ONLY  : LHOOK ,DR_HOOK
USE GRIB_API_INTERFACE

! ----------------------------------------------------------------------

IMPLICIT NONE
#include "abort1.intfb.h"
#include "mpcrtbl.intfb.h"
#include "outgrid.intfb.h"
#include "preset_wgrib_template.intfb.h"
#include "wgribencode.intfb.h"

INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL


INTEGER(KIND=JWIM), PARAMETER :: NBITSPERVALUE = 24  !! need a higher precision to code the grid point index 

INTEGER(KIND=JWIM) :: IJ, IJGLO, IK, IM, IR, IFLAG, IDT, KSTEP
INTEGER(KIND=JWIM) :: LFILE, IUOUT, ICOUNT
INTEGER(KIND=JWIM) :: IFCST, IT, IPARAM, ITABLE, IZLEV
INTEGER(KIND=JWIM) :: IGRIB_HANDLE
INTEGER(KIND=JWIM) :: ISIZE
INTEGER(KIND=JPKSIZE_T) :: KBYTES
INTEGER(KIND=JWIM), ALLOCATABLE :: KGRIB_BUFR(:)

REAL(KIND=JWRB) :: ZHOOK_HANDLE

CHARACTER(LEN= 2) :: MARSTYPE_DUM
CHARACTER(LEN=14) :: CDATE
CHARACTER(LEN=296) :: OUTFILE

LOGICAL :: LLCREATE
LOGICAL :: LGRHDIFS_DUM, LPADPOLES_DUM, LRSTST0_DUM
LOGICAL, DIMENSION(JPPFLAG) :: FFLAGBAK, GFLAGBAK, NFLAGBAK

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
FFLAGBAK(:) = FFLAG(:)
GFLAGBAK(:) = GFLAG(:)
NFLAGBAK(:) = NFLAG(:)

! Create a new output parameter selection that will only output
! parameters relevant for the description of the model decomposition
FFLAG(:) = .FALSE.
GFLAG(:) = .FALSE.
NFLAG(:) = .FALSE.

! Use the extra field codes (currently the last 5 fields of the list of potential output parameters:
! MPI RANK:
GFLAG(JPPFLAG-4) = .TRUE.

! GRID POINT INDEX (as used inside the wave model):
GFLAG(JPPFLAG-3) = .TRUE.

! GLOBAL POINT INDEX (as used for non parallel binary input/output)
! It is different than the local grid point index if LL1D is false
GFLAG(JPPFLAG-2) = .TRUE.

!!  GFLAG(JPPFLAG-1) and GFLAG(JPPFLAG) are still available if more entries are needed.


! Set output parameter mapping (and allocate BOUT) 
CALL MPCRTBL


! Defining the ouput fields:
! -------------------------
IR = 0

IR = IR + 1
BOUT(IJS:IJL, IR) = REAL(IRANK, JWRB)


IR = IR + 1
DO IJ = IJS, IJL
  BOUT(IJ, IR) = REAL(IJ, JWRB)
ENDDO


IR = IR + 1
IF ( .NOT. LL1D) THEN
  DO IJGLO = 1, NIBLO
    IJ = IJ2NEWIJ(IJGLO) 
    IF ( IJ >= IJS .AND. IJ <= IJL ) THEN 
       BOUT(IJ, IR) = REAL(IJGLO, JWRB)
    ENDIF
  ENDDO
ELSE
  DO IJ = IJS, IJL
    BOUT(IJ, IR) = REAL(IJ, JWRB)
  ENDDO
ENDIF


IF ( IR /= NIPRMOUT ) THEN
  WRITE(IU06,*) '******************************************'
  WRITE(IU06,*) ' OUTMDLDCP : ABORTING' 
  WRITE(IU06,*) ' IR /= NIPRMOUT ', IR, NIPRMOUT 
  WRITE(IU06,*) '******************************************'
  CALL ABORT1
ENDIF


! Gather data for output (to IRANK = 1)
CALL OUTGRID



IF(IRANK == 1) THEN
  ! Grib output to file:
  CALL IGRIB_OPEN_FILE(IUOUT,OUTFILE(1:LFILE),'w')

  ! Prepare grib template
  LLCREATE = .TRUE.
  CALL PRESET_WGRIB_TEMPLATE("I", NGRIB_HANDLE_WAM_I, LLCREATE=LLCREATE, NBITSPERVALUE=NBITSPERVALUE )


  ! keep looping over all posible output varaibles (as in outint)
  ICOUNT=0
  DO IFLAG = 1, JPPFLAG
    IF ( GFLAG(IFLAG) ) THEN
      ICOUNT = ICOUNT+1

      IT = ITOBOUT(IFLAG)
      ITABLE = INFOBOUT(IT,1)
      IPARAM = INFOBOUT(IT,2)
      IZLEV  = INFOBOUT(IT,3)

      IK = 0
      IM = 0

      CDATE = '20220220000000' ! set any date and time as the start date of the run might still be unknown
      IFCST = 0
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
      ISIZE=(KBYTES+NPRECI-1)/NPRECI
      ALLOCATE(KGRIB_BUFR(ISIZE))

      CALL IGRIB_GET_MESSAGE(IGRIB_HANDLE, KGRIB_BUFR)

      CALL IGRIB_WRITE_BYTES(IUOUT, KGRIB_BUFR, KBYTES)

      CALL IGRIB_RELEASE(IGRIB_HANDLE)

      DEALLOCATE(KGRIB_BUFR)

    ENDIF
  ENDDO

  CALL IGRIB_CLOSE_FILE(IUOUT)
  CALL IGRIB_RELEASE(NGRIB_HANDLE_WAM_I)

  IF(ALLOCATED(GOUT)) DEALLOCATE(GOUT)

ENDIF


! Restore output configuration: 
FFLAG(:) = FFLAGBAK(:)
GFLAG(:) = GFLAGBAK(:)
NFLAG(:) = NFLAGBAK(:)
!  reset output field mapping
CALL MPCRTBL

WRITE(IU06,*) ' '

IF (LHOOK) CALL DR_HOOK('OUTMDLDCP',1,ZHOOK_HANDLE)

END SUBROUTINE OUTMDLDCP
