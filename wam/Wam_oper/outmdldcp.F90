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

USE YOWCOUT  , ONLY : LFDB, JPPFLAG, FFLAG, GFLAG, NFLAG, NIPRMOUT, BOUT, ITOBOUT, INFOBOUT
USE YOWGRIB_HANDLES , ONLY : NGRIB_HANDLE_WAM_I, NGRIB_HANDLE_IFS
USE YOWGRIBHD, ONLY : LPADPOLES
USE YOWINTP  , ONLY : GOUT
USE YOWMPP   , ONLY : IRANK 
USE YOWPARAM , ONLY : NIBLO, NGX, NGY, LL1D
USE YOWSPEC  , ONLY : IJ2NEWIJ
USE YOWSTAT  , ONLY : MARSTYPE, CDATEA
USE YOWTEST  , ONLY : IU06, ITEST
USE YOWTEXT  , ONLY : ICPLEN, CPATH
USE YOWMAP   , ONLY : IXLG, KXLT
USE YOWGRID  , ONLY : NLONRGG

USE YOMHOOK  ,ONLY  : LHOOK ,DR_HOOK
USE GRIB_API_INTERFACE

! ----------------------------------------------------------------------

IMPLICIT NONE
#include "abort1.intfb.h"
#include "mpcrtbl.intfb.h"
#include "outgrid.intfb.h"
#include "preset_wgrib_template.intfb.h"
#include "wgribenout.intfb.h"

INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL


INTEGER(KIND=JWIM), PARAMETER :: NBITSPERVALUE = 24  !! need a higher precision to code the grid point index 

INTEGER(KIND=JWIM) :: IJ, IJGLO, IR, IFLAG, IX, JSN, I,J,K
INTEGER(KIND=JWIM) :: LFILE, IUOUT, ICOUNT
INTEGER(KIND=JWIM) :: IFCST, IT, IPARAM, ITABLE, IZLEV

REAL(KIND=JWRB) :: ZHOOK_HANDLE

CHARACTER(LEN= 2) :: MARSTYPEBAK
CHARACTER(LEN=14) :: CDATE
CHARACTER(LEN=296) :: OUTFILE

LOGICAL :: LFDBBAK, LPADPOLESBAK, LLCREATE
LOGICAL, DIMENSION(JPPFLAG) :: FFLAGBAK, GFLAGBAK, NFLAGBAK
! TEMPORARY ARRAYS FOR LOCAL MASK AND GLOBAL INDICES
INTEGER(KIND=JWIM) :: IGLOBAL(NGX,NGY)
INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:) :: NLOCMSK, NGLOIND

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
LPADPOLESBAK = LPADPOLES
FFLAGBAK(:) = FFLAG(:)
GFLAGBAK(:) = GFLAG(:)
NFLAGBAK(:) = NFLAG(:)
MARSTYPEBAK = MARSTYPE

! Create a new output parameter selection that will only output
! parameters relevant for the description of the model decomposition
LFDB = .FALSE.  ! data will be written to file
LPADPOLES = .FALSE.  ! Do not pad the poles
FFLAG(:) = .FALSE.
GFLAG(:) = .FALSE.
NFLAG(:) = .FALSE.
MARSTYPE = 'an'
IFCST = 0

! Use the extra field codes (currently the last 5 fields of the list of potential output parameters:
! MPI RANK:
GFLAG(JPPFLAG-8) = .TRUE.
GFLAG(JPPFLAG-7) = .TRUE.

! LAND sea mask
GFLAG(JPPFLAG-6) = .TRUE.

! GRID POINT INDEX (as used inside the wave model):
GFLAG(JPPFLAG-5) = .TRUE.
GFLAG(JPPFLAG-4) = .TRUE.

! GRID POINT INDEX (as used inside the wave model):
GFLAG(JPPFLAG-3) = .TRUE.
GFLAG(JPPFLAG-2) = .TRUE.

! GLOBAL POINT INDEX (as used for non parallel binary input/output)
! It is different than the local grid point index if LL1D is false
GFLAG(JPPFLAG-1) = .TRUE.
GFLAG(JPPFLAG) = .TRUE.

!!  GFLAG(JPPFLAG-1) and GFLAG(JPPFLAG) are still available if more entries are needed.


! Set output parameter mapping (and allocate BOUT) 
CALL MPCRTBL


! SETUP GLOBAL INDICES 1D ARRAY + LOCAL MASK
! SINCE WE ONLY HAVE OCEAN POINTS THE LOCAL MASK IS ONE FOR ALL POINTS
IGLOBAL(:,:) = -1
K=0
DO J = 1, NGY
  DO I = 1, NLONRGG(NGY-J+1)
    K= K + 1
    IGLOBAL(I,NGY-J+1) = K
  ENDDO
ENDDO
!
ALLOCATE( NLOCMSK(IJS:IJL), NGLOIND(IJS:IJL) )
DO IJ=IJS,IJL
   IX          = IXLG(IJ,1)
   JSN         = KXLT(IJ,1)
   NLOCMSK(IJ) = 1
   NGLOIND(IJ) = IGLOBAL(IX,JSN)
ENDDO

! Defining the ouput fields:
! -------------------------
IR = 0

IR = IR + 1
BOUT(IJS:IJL, IR) = REAL(MOD(IRANK,65536), JWRB)

IR = IR + 1
BOUT(IJS:IJL, IR) = REAL(IRANK/65536, JWRB)

IR = IR + 1
BOUT(IJS:IJL, IR) = REAL(NLOCMSK(IJS:IJL), JWRB)

IR = IR + 1
BOUT(IJS:IJL, IR) = REAL(MOD(NGLOIND(IJS:IJL),65536), JWRB)

IR = IR + 1
BOUT(IJS:IJL, IR) = REAL(NGLOIND(IJS:IJL)/65536, JWRB)

IR = IR + 1
DO IJ = IJS, IJL
  BOUT(IJ, IR) = REAL(MOD(IJ,65536), JWRB)
ENDDO

IR = IR + 1
DO IJ = IJS, IJL
  BOUT(IJ, IR) = REAL(IJ/65536, JWRB)
ENDDO


IR = IR + 1
IF ( .NOT. LL1D) THEN
  DO IJGLO = 1, NIBLO
    IJ = IJ2NEWIJ(IJGLO) 
    IF ( IJ >= IJS .AND. IJ <= IJL ) THEN 
       BOUT(IJ, IR) = REAL(MOD(IJGLO,65536), JWRB)
    ENDIF
  ENDDO
ELSE
  DO IJ = IJS, IJL
    BOUT(IJ, IR) = REAL(MOD(IJ,65536), JWRB)
  ENDDO
ENDIF

IR = IR + 1
IF ( .NOT. LL1D) THEN
  DO IJGLO = 1, NIBLO
    IJ = IJ2NEWIJ(IJGLO) 
    IF ( IJ >= IJS .AND. IJ <= IJL ) THEN 
       BOUT(IJ, IR) = REAL(IJGLO/65536, JWRB)
    ENDIF
  ENDDO
ELSE
  DO IJ = IJS, IJL
    BOUT(IJ, IR) = REAL(IJ/65536, JWRB)
  ENDDO
ENDIF

DEALLOCATE( NLOCMSK, NGLOIND )

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
  IF ( NGRIB_HANDLE_IFS < 0 ) NGRIB_HANDLE_IFS = NGRIB_HANDLE_WAM_I


  ! keep looping over all posible output varaibles (as in outint)
  ITABLE=140
  ICOUNT=0
  DO IFLAG = JPPFLAG-8, JPPFLAG


    IF ( GFLAG(IFLAG) ) THEN
      ICOUNT = ICOUNT+1

      SELECT CASE (IFLAG)
        CASE (JPPFLAG-8)
          IPARAM=80
          IZLEV=0
        CASE (JPPFLAG-7)
          IPARAM=80
          IZLEV=1
        CASE (JPPFLAG-6)
          IPARAM=81
          IZLEV=0
        CASE (JPPFLAG-5)
          IPARAM=82
          IZLEV=0
        CASE (JPPFLAG-4)
          IPARAM=82
          IZLEV=1
        CASE (JPPFLAG-3)
          IPARAM=83
          IZLEV=0
        CASE (JPPFLAG-2)
          IPARAM=83
          IZLEV=1
        CASE (JPPFLAG-1)
          IPARAM=84
          IZLEV=0
        CASE (JPPFLAG)
          IPARAM=84
          IZLEV=1
      END SELECT

      CDATE = '20220220000000' ! set any date and time as the start date of the run might still be unknown

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
LPADPOLES = LPADPOLESBAK
FFLAG(:) = FFLAGBAK(:)
GFLAG(:) = GFLAGBAK(:)
NFLAG(:) = NFLAGBAK(:)
MARSTYPE = MARSTYPEBAK
!  reset output field mapping
CALL MPCRTBL

WRITE(IU06,*) ' '

IF (LHOOK) CALL DR_HOOK('OUTMDLDCP',1,ZHOOK_HANDLE)

END SUBROUTINE OUTMDLDCP
