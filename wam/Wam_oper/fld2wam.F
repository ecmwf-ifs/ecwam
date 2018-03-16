      SUBROUTINE FLD2WAM (IPARAM, FILNM, IJS, IJL, IG, VALMISS,
     &                    BLOCK, CDATEIN)

!--------------------------------------------------------------------

!**** *FLD2WAM* 

!     PURPOSE.
!     --------

!     TO READ !ONE! GRIB FIELD IN AND TO TRANSFER IT ONTO THE
!     WAVE MODEL LOCAL SEA POINTS.

!**   INTERFACE
!     ---------
!     *CALL* *FLD2WAM(IPARAM, FILNM, IJS, IJL, IG, VALMISS,
!    &                    BLOCK, CDATEIN)


!     *IPARAM*    GRIB PARAMETER EXPECTED.
!     *FILNM*     DATA INPUT FILENAME.
!            !!!  FILNM WILL NOT BE KEPT OPENED !!!!
!     *IJS*       INDEX OF FIRST GRIDPOINT
!     *IJL*       INDEX OF LAST GRIDPOINT
!     *IG*        BLOCK NUMBER
!     *VALMISS*   VALUE GIVEN TO BLOCK IF MISSING DATA ARE DECODED
!     *BLOCK*     LOCAL MODEL SEA POINT VALUES OF INPUT FIELD
!     *CDATEIN*   DATE OF THE DECODED DATA. 


!     METHOD.
!     -------


!    EXTERNAL.
!    ---------

!---------------------------------------------------------------------

      USE YOWCURR  , ONLY : U        ,V
      USE YOWGRIBHD, ONLY : PPEPS    ,PPREC
      USE YOWGRID  , ONLY : NLONRGG
      USE YOWMAP   , ONLY : IRGG     ,XDELLA   ,XDELLO   ,ZDELLO   ,
     &            IFROMIJ  ,JFROMIJ
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPRECI
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSTAT  , ONLY : NPROMA_WAM 
      USE YOWSPEC  , ONLY : NSTART   ,NEND
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWPD, ONLY : MNP => npa
      USE YOWUNPOOL ,ONLY : LLUNSTR
      USE YOWWIND  , ONLY : NXFF     ,NYFF     ,FIELDG

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE MPL_MODULE
      USE GRIB_API_INTERFACE

!-----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "grib2wgrid.intfb.h"
#include "kgribsize.intfb.h"

      INTEGER, INTENT(IN) :: IPARAM, IJS, IJL, IG

      REAL, INTENT(IN) :: VALMISS 
      REAL, INTENT(OUT) :: BLOCK(IJS:IJL)

      CHARACTER(LEN=24), INTENT(IN) :: FILNM
      CHARACTER(LEN=14), INTENT(OUT) :: CDATEIN

      INTEGER :: NBIT = 550000

      INTEGER :: IREAD
      INTEGER :: KFILE_HANDLE1
      INTEGER :: LFILE, KGRIB_HANDLE
      INTEGER :: KRET, KPLENG, KLEN, KLENG, IDM
      INTEGER :: IRET, ISIZE
      INTEGER :: KK, MM
      INTEGER :: IFORP, IPARAM_IN, KZLEV, IJ, IC, IX, JY
      INTEGER :: IBUF(2)
      INTEGER, ALLOCATABLE :: INGRIB(:)
      INTEGER :: NLONRGG_LOC(NYFF)
      INTEGER(KIND=JPKSIZE_T) :: KBYTES

      REAL :: ZDUM
      REAL :: ZHOOK_HANDLE
      REAL :: FIELD(NXFF,NYFF)


! --------------------------------------------------------------------  
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('FLD2WAM',0,ZHOOK_HANDLE)
#endif

      IREAD=1

      IF(LLUNSTR) THEN
        NLONRGG_LOC(:)=MNP
      ELSE
        NLONRGG_LOC(:)=NLONRGG(:)
      ENDIF

      LFILE = LEN_TRIM(FILNM)

      IF ((LMESSPASS .AND. (IRANK.EQ.IREAD)) .OR.
     &    .NOT.LMESSPASS) THEN

        KFILE_HANDLE1=-99
        CALL IGRIB_OPEN_FILE(KFILE_HANDLE1,FILNM(1:LFILE),'r')
        ISIZE=NBIT
      ENDIF


        IF ((LMESSPASS .AND. (IRANK.EQ.IREAD)) .OR.
     &      .NOT.LMESSPASS) THEN

1021      ISIZE=NBIT
          KBYTES=ISIZE*NPRECI
          IF(.NOT.ALLOCATED(INGRIB)) ALLOCATE(INGRIB(ISIZE))
          CALL IGRIB_READ_FROM_FILE(KFILE_HANDLE1,INGRIB,KBYTES,IRET)
          IF(IRET.EQ.JPGRIB_BUFFER_TOO_SMALL) THEN
!!!         *IGRIB_READ_FROM_FILE* does not read through the file if
!!!          the size is too small, so figure out the size and read again.
            CALL KGRIBSIZE(IU06, KBYTES, NBIT, 'FLD2WAM')
            DEALLOCATE(INGRIB)
            GOTO 1021
          ELSEIF(IRET.EQ.JPGRIB_END_OF_FILE) THEN
             WRITE(IU06,*) '**************************************'
             WRITE(IU06,*) '* FLD2WAM: END OF FILE ENCOUNTED *'
             WRITE(IU06,*) '**************************************'
             CALL ABORT1
          ELSEIF(IRET.NE.JPGRIB_SUCCESS) THEN
             WRITE(IU06,*) '************************************'
             WRITE(IU06,*) '* FLD2WAM: FILE HANDLING ERROR *'
             WRITE(IU06,*) '************************************'
             CALL ABORT1
          ENDIF

          WRITE(IU06,*) ' SUB. FLD2WAM - READ FROM ',FILNM
        ENDIF

!       BROADCAST GRIB DATA TO OTHER PE'S
        CALL GSTATS(622,0)
        IF(LMESSPASS .AND. NPROC.GT.1) THEN

          CALL MPL_BARRIER(CDSTRING='FLD2WAM: INGRIB ')

          IF(IRANK.EQ.IREAD) THEN
            IBUF(1)=ISIZE
            IBUF(2)=KLEN
          ENDIF
          CALL MPL_BROADCAST(IBUF(1:2),KROOT=IREAD,KTAG=1,
     &                       CDSTRING='FLD2WAM IBUF:')
          IF(IRANK.NE.IREAD) THEN
            ISIZE=IBUF(1)
            KLEN=IBUF(2)
            ALLOCATE(INGRIB(ISIZE))
          ENDIF

          CALL MPL_BROADCAST(INGRIB(1:ISIZE),KROOT=IREAD,KTAG=2,
     &                       CDSTRING='FLD2WAM INGRIB:')

        ENDIF
        CALL GSTATS(622,1)
! ----------------------------------------------------------------------

!*      UNPACK GRIB FIELDS.
!       -------------------
        CALL IGRIB_NEW_FROM_MESSAGE(KGRIB_HANDLE,INGRIB)

        CALL GRIB2WGRID (IU06, ITEST, NPROMA_WAM, 
     &                   KGRIB_HANDLE, INGRIB, ISIZE,
     &                   LLUNSTR,
     &                   NXFF, NYFF, NLONRGG_LOC,
     &                   IRGG, XDELLA, ZDELLO,
     &                   FIELDG%XLON, FIELDG%YLAT,
     &                   ZMISS, PPREC, PPEPS,
     &                   CDATEIN, IFORP, IPARAM_IN,KZLEV,KK,MM,FIELD)

        CALL IGRIB_RELEASE(KGRIB_HANDLE)

        IF(IPARAM_IN.EQ.IPARAM) THEN
           DO IJ = IJS, IJL 
             IX = IFROMIJ(IJ,IG)
             JY = JFROMIJ(IJ,IG)
             BLOCK(IJ) = FIELD(IX,JY)
!            SOME WAM MODEL GRID POINTS MAY HAVE A MISSING DATA FROM
!            OCEAN MODEL. THEY ARE SET TO VALMISS 
             IF(BLOCK(IJ).LE.ZMISS) BLOCK(IJ)=VALMISS
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


#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('FLD2WAM',1,ZHOOK_HANDLE)
#endif

      RETURN
      END SUBROUTINE FLD2WAM
