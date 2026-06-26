      SUBROUTINE FLD2WAM (IPARAM, FILNM, IJS, IJL, IFROMIJ, JFROMIJ, VALMISS, &
     &                    NXS, NXE, NYS, NYE, FIELDG,                         &
     &                    BLOCK, CDATEIN)

!--------------------------------------------------------------------

!**** *FLD2WAM* 

!     PURPOSE.
!     --------

!     TO READ !ONE! GRIB FIELD IN AND TO TRANSFER IT ONTO THE
!     WAVE MODEL LOCAL SEA POINTS.

!**   INTERFACE
!     ---------
!     *CALL* *FLD2WAM(IPARAM, FILNM, IJS, IJL, VALMISS,
!     &               NXS, NXE, NYS, NYE, FIELDG,
!    &                BLOCK, CDATEIN)


!     *IPARAM*    GRIB PARAMETER EXPECTED.
!     *FILNM*     DATA INPUT FILENAME.
!            !!!  FILNM WILL NOT BE KEPT OPENED !!!!
!     *IJS*       INDEX OF FIRST GRIDPOINT
!     *IJL*       INDEX OF LAST GRIDPOINT
!     *VALMISS*   VALUE GIVEN TO BLOCK IF MISSING DATA ARE DECODED
!     *NXS:NXE*  FIRST DIMENSION OF FIELDG
!     *NYS:NYE*  SECOND DIMENSION OF FIELDG
!     *FIELDG* - INPUT FORCING FIELDS ON THE WAVE MODEL GRID
!     *BLOCK*     LOCAL MODEL SEA POINT VALUES OF INPUT FIELD
!     *CDATEIN*   DATE OF THE DECODED DATA. 


!     METHOD.
!     -------


!    EXTERNAL.
!    ---------

!---------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : FORCING_FIELDS

      USE YOWGRIBHD, ONLY : PPEPS    ,PPREC
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWMAP   , ONLY : IRGG     ,NGY      ,NLONRGG
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPRECI
      USE YOWPARAM , ONLY : LLUNSTR
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSPEC  , ONLY : NSTART   ,NEND
      USE YOWTEST  , ONLY : IU06
      USE YOWPD, ONLY : MNP => npa

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE MPL_MODULE
      USE YOWGRIB

!-----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "grib2wgrid.intfb.h"
#include "kgribsize.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IPARAM
      CHARACTER(LEN=24), INTENT(IN) :: FILNM
      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL), INTENT(IN) :: IFROMIJ  ,JFROMIJ
      REAL(KIND=JWRB), INTENT(IN) :: VALMISS 
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), INTENT(IN) :: FIELDG
      REAL(KIND=JWRB), INTENT(OUT) :: BLOCK(IJS:IJL)
      CHARACTER(LEN=14), INTENT(OUT) :: CDATEIN


      INTEGER(KIND=JWIM) :: NBIT = 1600000

      INTEGER(KIND=JWIM) :: IREAD
      INTEGER(KIND=JWIM) :: KFILE_HANDLE1
      INTEGER(KIND=JWIM) :: LFILE, KGRIB_HANDLE
      INTEGER(KIND=JWIM) :: IRET, ISIZE
      INTEGER(KIND=JWIM) :: KK, MM
      INTEGER(KIND=JWIM) :: IFORP, IPARAM_IN, KZLEV, IJ, IC, IX, JY
      INTEGER(KIND=JWIM) :: IBUF(1)
      INTEGER(KIND=JWIM), ALLOCATABLE :: INGRIB(:)
      INTEGER(KIND=JWIM) :: NLONRGG_LOC(NGY)
      INTEGER(KIND=JPKSIZE_T) :: KBYTES

      REAL(KIND=JWRB) :: ZDUM
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NXS:NXE, NYS:NYE) :: FIELD 

      LOGICAL :: LLCHKINT
! --------------------------------------------------------------------  

      IF (LHOOK) CALL DR_HOOK('FLD2WAM',0,ZHOOK_HANDLE)

      IREAD=1

      DO IJ = IJS, IJL
        IX = IFROMIJ(IJ)
        JY = JFROMIJ(IJ)
        IF (IX < NXS .OR. IX > NXE .OR. JY < NYS .OR. JY > NYE) THEN
          WRITE(IU06,*) '**************************************'
          WRITE(IU06,*) '* FLD2WAM: IJS to IJL SPAN TOO MUCH !'
          WRITE(IU06,*) '* IX, JY = ', IX, JY
          WRITE(IU06,*) '* NXS, NXE  = ', NXS, NXE
          WRITE(IU06,*) '* NYS, NYE  = ', NYS, NYE
          WRITE(IU06,*) '**************************************'
          CALL ABORT1
        ENDIF
      ENDDO

      IF (LLUNSTR) THEN
        NLONRGG_LOC(:)=MNP
      ELSE
        NLONRGG_LOC(:)=NLONRGG(:)
      ENDIF

      LFILE = LEN_TRIM(FILNM)

      IF (IRANK == IREAD) THEN

        KFILE_HANDLE1=-99
        CALL IGRIB_OPEN_FILE(KFILE_HANDLE1,FILNM(1:LFILE),'r')
        ISIZE=NBIT

1021    ISIZE=NBIT
        KBYTES=ISIZE*NPRECI
        IF (.NOT.ALLOCATED(INGRIB)) ALLOCATE(INGRIB(ISIZE))
        CALL IGRIB_READ_FROM_FILE(KFILE_HANDLE1,INGRIB,KBYTES,IRET)
        IF (IRET == JPGRIB_BUFFER_TOO_SMALL) THEN
!!!        *IGRIB_READ_FROM_FILE* does not read through the file if
!!!         the size is too small, so figure out the size and read again.
          CALL KGRIBSIZE(IU06, KBYTES, NBIT, 'FLD2WAM')
          DEALLOCATE(INGRIB)
          GOTO 1021
        ELSEIF (IRET == JPGRIB_END_OF_FILE) THEN
          WRITE(IU06,*) '**************************************'
          WRITE(IU06,*) '* FLD2WAM: END OF FILE ENCOUNTED *'
          WRITE(IU06,*) '**************************************'
          CALL ABORT1
        ELSEIF (IRET /= JPGRIB_SUCCESS) THEN
          WRITE(IU06,*) '************************************'
          WRITE(IU06,*) '* FLD2WAM: FILE HANDLING ERROR *'
          WRITE(IU06,*) '************************************'
          CALL ABORT1
        ENDIF

        WRITE(IU06,*) ' SUB. FLD2WAM - READ FROM ',FILNM
      ENDIF


!     BROADCAST GRIB DATA TO OTHER PE'S
      CALL GSTATS(622,0)
      IF (NPROC > 1) THEN

        CALL MPL_BARRIER(CDSTRING='FLD2WAM: INGRIB ')

        IF (IRANK == IREAD) THEN
          IBUF(1)=ISIZE
        ENDIF
        CALL MPL_BROADCAST(IBUF(1:1),KROOT=IREAD,KTAG=1,                &
     &                       CDSTRING='FLD2WAM IBUF:')
        IF (IRANK /= IREAD) THEN
          ISIZE=IBUF(1)
          ALLOCATE(INGRIB(ISIZE))
        ENDIF

        CALL MPL_BROADCAST(INGRIB(1:ISIZE),KROOT=IREAD,KTAG=2,          &
     &                       CDSTRING='FLD2WAM INGRIB:')

      ENDIF
      CALL GSTATS(622,1)
! ----------------------------------------------------------------------

!*    UNPACK GRIB FIELDS.
!     -------------------
      CALL IGRIB_NEW_FROM_MESSAGE(KGRIB_HANDLE,INGRIB)

      LLCHKINT = .TRUE.
      CALL GRIB2WGRID (IU06, NPROMA_WAM,                              &
     &                 KGRIB_HANDLE, INGRIB, ISIZE,                   &
     &                 LLUNSTR, LLCHKINT,                             &
     &                 NGY, IRGG, NLONRGG_LOC,                        &
     &                 NXS, NXE, NYS, NYE,                            &
     &                 FIELDG%XLON, FIELDG%YLAT,                      &
     &                 ZMISS, PPREC, PPEPS,                           &
     &                 CDATEIN, IFORP, IPARAM_IN,KZLEV, KK, MM, FIELD)

      CALL IGRIB_RELEASE(KGRIB_HANDLE)

      IF (IPARAM_IN == IPARAM) THEN
        DO IJ = IJS, IJL
          IX = IFROMIJ(IJ)
          JY = JFROMIJ(IJ)
          BLOCK(IJ) = FIELD(IX,JY)
!         SOME WAM MODEL GRID POINTS MAY HAVE A MISSING DATA FROM
!         OCEAN MODEL. THEY ARE SET TO VALMISS 
          IF (BLOCK(IJ) /= ZMISS) BLOCK(IJ)=VALMISS
        ENDDO
      ELSE
        WRITE(IU06,*) ' *****************************************'
        WRITE(IU06,*) ' *                                       *'
        WRITE(IU06,*) ' *    FATAL ERROR IN SUB. FLD2WAM        *'
        WRITE(IU06,*) ' *      ============================     *'
        WRITE(IU06,*) ' * THE INPUT PARAMETER ',IPARAM_IN 
        WRITE(IU06,*) ' * IS NOT THE REQUESTED ONE: ',IPARAM 
        WRITE(IU06,*) ' *                                       *'
        WRITE(IU06,*) ' *  PROGRAM ABORTS     PROGRAM ABORTS    *'
        WRITE(IU06,*) ' *                                       *'
        WRITE(IU06,*) ' *****************************************'
        CALL ABORT1
      ENDIF

      DEALLOCATE(INGRIB)

      IF (LHOOK) CALL DR_HOOK('FLD2WAM',1,ZHOOK_HANDLE)

      END SUBROUTINE FLD2WAM
