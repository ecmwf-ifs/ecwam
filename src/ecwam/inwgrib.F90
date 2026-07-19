! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE INWGRIB (FILNM, IREAD, CDATE, IPARAMID, KZLEV,           &
     &                    NXS, NXE, NYS, NYE, FIELDG, FIELD,              &
     &                    LLCHKINT, LLFIXEDSIZE, NPR,                     &
     &                    KANGNB, KFRENB, NFILE_HANDLE)

! -----------------------------------------------------------------     

!***  *INWGRIB* - READS IN AND UNPACKS WAVE MODEL GRIB FIELDS
!                 IT WILL ALSO INTERPOLATE THE DATA TO THE MODEL
!                 GRID (IF NECESSARY).
!                 (Except for wave spectra. See *GETSPEC*)
!                 MESSAGE PASSING WILL BE USE TO DISTRIBUTED THE
!                 CODED DATA TO ALL PE'S BEFORE DECODING.

!      J. BIDLOT    ECMWF    APRIL 2010. 

!**   INTERFACE.                                                        
!     ----------                                                        

!      *CALL INWGRIB* (FILNM, IREAD, CDATE, IPARAMID, KZLEV,
!    &                 NXS, NXE, NYS, NYE, FIELDG, FIELD,
!    &                 LLCHKINT, LLFIXEDSIZE, NPR,
!    &                 KANGNB, KFRENB, NFILE_HANDLE)

!        *FILNM*  - DATA INPUT FILENAME. 
!        *IREAD*  - ACCESS TO FILE ONLY FOR PE=IREAD.
!        *CDATE*  - DATE/TIME OF THE DATA READ.         
!        *IPARAMID* PARAMETER IDENTIFIER NUMBER. 
!        *KZLEV*  - REFERENCE LEVEL IN FULL METER
!                   SHOULD BE 0 EXCEPT FOR 140233, 140245 AND 140249 WHERE IT
!                   MIGHT BE DIFFERENT THAN ZERO. 
!        *NXS:NXE*  FIRST DIMENSION OF FIELDG
!        *NYS:NYE*  SECOND DIMENSION OF FIELDG
!        *FIELDG* - INPUT FORCING FIELDS ON THE WAVE MODEL GRID
!        *FIELD*  - UNPACKED DATA.

!        OPTIONAL INPUT:
!        *LLCHKINT*    - IF TRUE CHECK ON WHETHER OR NOT THE INPUT DATA WILL NEED TO BE INTERPOLATED TO THE MODEL GRID
!                        OTHERWISE INPUT AND TARGET GRID ARE ASSUMED TO BE THE SAME
!        *LLFIXEDSIZE* - IF TRUE THE SIZE OF THE GRIB RECORD WILL BE FIXEDa TO ITS FIRST GUESS (NBIT)
!                        RATHER THAN TRYING TO ADJUST IT.
!                        IF IT IS TOO SMALL, THE RUN WILL ABORT !!
!        *NPR*    - NUMBER OF SUBDOMAINS (USUALLY THE NUMBER OF PE'S )

!        OPTIONAL OUTPUT
!        *KANGNB*  INDEX FOR SPECTRAL DIRECTION
!        *KFRENB*  INDEX FOR SPECTRAL FREQUENCY

!     EXTERNALS.                                                        
!     ----------                                                        

!     GRIBAPI 

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : FORCING_FIELDS

      USE YOWABORT,  ONLY : WAM_ABORT
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWGRIBHD, ONLY : PPEPS    ,PPREC
      USE YOWPARAM , ONLY : LLUNSTR
      USE YOWMAP   , ONLY : NGY      ,NIBLO    ,IRGG     ,NLONRGG
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPRECI 
      USE YOWPCONS , ONLY : ZMISS
      USE YOWTEST  , ONLY : IU06
#ifdef WAM_HAVE_UNWAM
      USE YOWPD, ONLY : MNP => npa
#endif

      USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE EC_LUN   ,  ONLY : NULERR
      USE MPL_MODULE, ONLY : MPL_BARRIER, MPL_BROADCAST
      USE YOWGRIB, ONLY    : IGRIB_OPEN_FILE, IGRIB_CLOSE_FILE, IGRIB_RELEASE, &
                           & IGRIB_NEW_FROM_MESSAGE, IGRIB_READ_FROM_FILE, &
                           & JPKSIZE_T, JPGRIB_BUFFER_TOO_SMALL, &
                           & JPGRIB_END_OF_FILE, JPGRIB_SUCCESS

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "grib2wgrid.intfb.h"
#include "kgribsize.intfb.h"

      CHARACTER(LEN=*), INTENT(IN)     :: FILNM
      INTEGER(KIND=JWIM), INTENT(IN)   :: IREAD
      CHARACTER(LEN=14), INTENT(OUT)   :: CDATE
      INTEGER(KIND=JWIM), INTENT(OUT)  :: IPARAMID, KZLEV
      INTEGER(KIND=JWIM), INTENT(IN)   :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), INTENT(IN) :: FIELDG
      REAL(KIND=JWRB), INTENT(OUT)     :: FIELD(NXS:NXE, NYS:NYE)

      LOGICAL, INTENT(IN), OPTIONAL               :: LLCHKINT 
      LOGICAL, INTENT(IN), OPTIONAL               :: LLFIXEDSIZE
      INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL    :: NPR 
      INTEGER(KIND=JWIM), INTENT(INOUT), OPTIONAL :: NFILE_HANDLE

      INTEGER(KIND=JWIM), INTENT(OUT), OPTIONAL   :: KANGNB
      INTEGER(KIND=JWIM), INTENT(OUT), OPTIONAL   :: KFRENB 


      INTEGER(KIND=JWIM), SAVE :: NBIT

      INTEGER(KIND=JWIM) :: IPARAM 
      INTEGER(KIND=JWIM) :: NPRC 
      INTEGER(KIND=JWIM) :: IFORP
      INTEGER(KIND=JWIM) :: LFILE, KFILE_HANDLE, KGRIB_HANDLE 
      INTEGER(KIND=JWIM) :: IRET, ISIZE
      INTEGER(KIND=JWIM) :: KK, MM
      INTEGER(KIND=JWIM) :: IBUF(2) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: INGRIB(:)
      INTEGER(KIND=JWIM) :: NLONRGG_LOC(NGY)

      INTEGER(KIND=JPKSIZE_T) :: KBYTES

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      LOGICAL :: LLEXIST
      LOGICAL :: LLOPENFILE, LLHANDEL, LLFIX, LLCHK

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('INWGRIB',0,ZHOOK_HANDLE)


      IF( PRESENT(NFILE_HANDLE) ) THEN
        LLHANDEL = .TRUE.
        IF ( NFILE_HANDLE < 0 .AND. IRANK == IREAD ) THEN
          LLOPENFILE = .TRUE.
        ELSE
          LLOPENFILE = .FALSE.
        ENDIF
      ELSE
        LLHANDEL = .FALSE.
        LLOPENFILE = .TRUE.
      ENDIF

      IF( PRESENT(LLCHKINT) ) THEN
        LLCHK = LLCHKINT 
      ELSE
        LLCHK = .TRUE.
      ENDIF

      IF( PRESENT(LLFIXEDSIZE) ) THEN
        LLFIX = LLFIXEDSIZE
      ELSE
        LLFIX = .FALSE.
      ENDIF

      IF( PRESENT(NPR) ) THEN
        NPRC = NPR
      ELSE
        NPRC = NPROC
      ENDIF

      NBIT=NIBLO
      ISIZE=NBIT

      IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
        NLONRGG_LOC(:)=MNP
#else
        CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      ELSE
        NLONRGG_LOC(:)=NLONRGG(:)
      ENDIF

!     READ DATA ON PE IREAD
      IF (IRANK == IREAD) THEN
        IF (LLOPENFILE) THEN
          !! INPUT FROM FILE in FILNM
          LLEXIST=.FALSE.
          LFILE = LEN_TRIM(FILNM)
          INQUIRE(FILE=FILNM(1:LFILE),EXIST=LLEXIST)
          IF (.NOT. LLEXIST .OR. LFILE == 0) THEN
            WRITE (IU06,*) '*************************************'
            WRITE (IU06,*) '*                                   *'
            WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO INQUIRE  *'
            WRITE (IU06,*) '*  IN INWGRIB:                      *'
            WRITE (IU06,*) '*  COULD NOT FIND FILE ',FILNM
            WRITE (IU06,*) '*                                   *'
            WRITE (IU06,*) '*************************************'
            WRITE (NULERR,*) '*************************************'
            WRITE (NULERR,*) '*                                   *'
            WRITE (NULERR,*) '*  ERROR FOLLOWING CALL TO INQUIRE  *'
            WRITE (NULERR,*) '*  IN INWGRIB:                      *'
            WRITE (NULERR,*) '*  COULD NOT FIND FILE ',FILNM
            WRITE (NULERR,*) '*                                   *'
            WRITE (NULERR,*) '*************************************'
            CALL ABORT1
          ENDIF

          CALL IGRIB_OPEN_FILE(KFILE_HANDLE,FILNM(1:LFILE),'r')

          IF( LLHANDEL ) NFILE_HANDLE = KFILE_HANDLE
        ELSE
          IF( LLHANDEL ) KFILE_HANDLE = NFILE_HANDLE
        ENDIF

1021    ISIZE=NBIT
        KBYTES=ISIZE*NPRECI
        IF (.NOT.ALLOCATED(INGRIB)) ALLOCATE(INGRIB(ISIZE))
          CALL IGRIB_READ_FROM_FILE(KFILE_HANDLE,INGRIB,KBYTES,IRET)

        IF ( IRET == JPGRIB_BUFFER_TOO_SMALL .AND. LLFIX ) THEN
          WRITE(IU06,*) '****************************************************'
          WRITE(IU06,*) '* INWGRIB: BUFFER TOO SMALL BUT LLFIXEDSIZE IS TRUE'
          WRITE(IU06,*) '* FILE: ',FILNM(1:LFILE)
          WRITE(NULERR,*) '* INWGRIB: BUFFER TOO SMALL BUT LLFIXEDSIZE IS TRUE'
          WRITE(NULERR,*) '* FILE: ',FILNM(1:LFILE)
          WRITE(IU06,*) '****************************************************'
          CALL ABORT1
        ELSEIF ( IRET == JPGRIB_BUFFER_TOO_SMALL .AND. .NOT. LLFIX ) THEN
!!!       *IGRIB_READ_FROM_FILE* does not read through the file if
!!!       the size is too small, so figure out the size and read again.
          CALL KGRIBSIZE(IU06, KBYTES, NBIT, 'INWGRIB')
          DEALLOCATE(INGRIB)
          GOTO 1021
        ELSEIF ( IRET == JPGRIB_END_OF_FILE ) THEN
          WRITE(IU06,*) '**********************************'
          WRITE(IU06,*) '* INWGRIB: END OF FILE ENCOUNTED'
          WRITE(IU06,*) '* FILE: ',FILNM(1:LFILE)
          WRITE(NULERR,*) '* INWGRIB: END OF FILE ENCOUNTED'
          WRITE(NULERR,*) '* FILE: ',FILNM(1:LFILE)
          WRITE(IU06,*) '**********************************'
          CALL ABORT1
        ELSEIF ( IRET /= JPGRIB_SUCCESS ) THEN
          WRITE(IU06,*) '**********************************'
          WRITE(IU06,*) '* INWGRIB: FILE HANDLING ERROR'
          WRITE(IU06,*) '* FILE: ',FILNM(1:LFILE)
          WRITE(NULERR,*) '* INWGRIB: FILE HANDLING ERROR'
          WRITE(NULERR,*) '* FILE: ',FILNM(1:LFILE)
          WRITE(IU06,*) '**********************************'
          CALL ABORT1
        ENDIF

        IF (LLOPENFILE) THEN
          WRITE(IU06,*) ' SUB. INWGRIB - READ FROM ',FILNM(1:LFILE)
        ENDIF

      ENDIF


!     SEND GRIB DATA TO THE OTHER PE'S
      IF (NPRC > 1) THEN
        CALL GSTATS(619,0)

        IF ( .NOT. LLFIX ) THEN
           IF (IRANK == IREAD) THEN
            IBUF(1)=ISIZE
            IBUF(2)=KBYTES
          ENDIF
          CALL MPL_BROADCAST(IBUF(1:2),KROOT=IREAD, KTAG=1, CDSTRING='INWGRIB IBUF:') 
          IF (IRANK /= IREAD) THEN
            ISIZE=IBUF(1)
            KBYTES=IBUF(2)
          ENDIF
        ENDIF

        IF (IRANK /= IREAD) ALLOCATE(INGRIB(ISIZE))
        CALL MPL_BROADCAST(INGRIB(1:ISIZE), KROOT=IREAD, KTAG=2, CDSTRING='INWGRIB: INGRIB')
        CALL GSTATS(619,1)
      ENDIF

      CALL MPL_BARRIER(CDSTRING='INWGRIB: DATA INPUT IN')

!     DECODE THE GRIB DATA
!     (and interpolate to model grid if necessary)

      KGRIB_HANDLE=-99
      CALL IGRIB_NEW_FROM_MESSAGE(KGRIB_HANDLE,INGRIB)

      CALL GRIB2WGRID (IU06, NPROMA_WAM,                                &
     &                 KGRIB_HANDLE, INGRIB, ISIZE,                     &
     &                 LLUNSTR, LLCHK,                                  &
     &                 NGY, IRGG, NLONRGG_LOC,                          &
     &                 NXS, NXE, NYS, NYE,                              &
     &                 FIELDG%XLON, FIELDG%YLAT,                        &
     &                 ZMISS, PPREC, PPEPS,                             &
     &                 CDATE, IFORP, IPARAM, KZLEV, KK, MM, FIELD,      &
     &                 IPARAMID=IPARAMID)

      IF( PRESENT(KANGNB) ) KANGNB = KK 
      IF( PRESENT(KFRENB) ) KFRENB = MM 

      CALL IGRIB_RELEASE(KGRIB_HANDLE)

      IF (ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)

      IF (LHOOK) CALL DR_HOOK('INWGRIB',1,ZHOOK_HANDLE)

      END SUBROUTINE INWGRIB
