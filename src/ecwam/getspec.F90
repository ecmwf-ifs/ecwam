! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE GETSPEC(FL1, BLK2GLO, BLK2LOC, WVENVI, NBLKS, NBLKE, IREAD)
! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      SEPTEMBER 1997 
!     J. BIDLOT    ECMWF      MARCH 2010: modified to use gribapi 

!*    PURPOSE.
!     --------
!     READ THE SPECTRA FROM DISK.

!**   INTERFACE.
!     ----------
!     *CALL* *GETSPEC(FL1, BLK2GLO, WVENVI,NBLKS,NBLKE,IREAD)
!     *FL1       ARRAY CONTAINING THE SPECTRA CONTRIBUTION ON EACH PE

!     *BLK2GLO*  BLOCK TO GRID TRANSFORMATION
!     *BLK2LOC*  POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!     *WVENVI*   WAVE ENVIRONMENT FIELDS
!     *NBLKS*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!     *NBLKE*    INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!     *IREAD*    PROCESSOR WHICH WILL ACCESS THE FILE ON DISK 


!     METHOD.
!     -------

!     IN CASE THE INPUT SPECTRA ARE IN GRIB THEN EITHER THEY ARE READ
!     FROM A FILE ON DISK.
!     GRIB SPECTRA WILL BE DECODED. IN CASE OF BINARY
!     DATA, USE READFL TO READ IN THE SPECTRA DEPENDING ON THE USE OF
!     THE PBIO SOFTWARE OR NOT. PBIO WILL LIMIT THE SIZE OF THE ARRAY
!     NECESSARY TO READ THE INPUT SPECTRA. THE READING IS ONLY DONE ON
!     PE 1, THEREFORE THE RELEVANT INFORMATION IS SENT TO THE OTHER
!     PE'S USING MPDISTRIBFL

!     EXTERNALS.
!     ----------
!     GETENV
!     GRSTNAME
!     MPDISTRIBFL
!     MPDISTRIBSCFLD
!     MPL_BARRIER
!     READFL

!     REFERENCE.
!     ----------
!     NONE

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO, WVGRIDLOC, ENVIRONMENT, FORCING_FIELDS

      USE YOWCOUT  , ONLY : KDEL     ,MDEL     ,LRSTPARALR
      USE YOWFRED  , ONLY : FR       ,TH       ,FR5      ,FRM5
      USE YOWGRIBHD, ONLY : PPEPS    ,PPREC
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, KIJL4CHNK, IJFROMCHNK
      USE YOWMAP   , ONLY : IRGG     ,NLONRGG
      USE YOWMESPAS, ONLY : LGRIBIN
      USE YOWMAP   , ONLY : NGY      ,NIBLO
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,                        &
     &                      KTAG     ,NPRECR   ,NPRECI
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_RED ,LLUNSTR
      USE YOWPCONS , ONLY : G        ,DEG      ,R        ,ZMISS    ,    &
     &                      EPSMIN
      USE YOWSTAT  , ONLY : CDATEF   ,CDTPRO   ,IREFRA  ,LNSESTART
      USE YOWTEST  , ONLY : IU06
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH    ,LRESTARTED
#ifdef WAM_HAVE_UNWAM
      USE YOWPD, ONLY : MNP => npa
#endif
      USE YOWWIND  , ONLY : NXFFS    ,NXFFE    ,NYFFS    ,NYFFE

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE MPL_MODULE, ONLY : MPL_RECV, MPL_SEND, MPL_WAIT, MPL_ABORT, &
                           & MPL_BARRIER, MPL_WAITANY, &
                           & JP_NON_BLOCKING_STANDARD, JP_BLOCKING_STANDARD
      USE YOWGRIB, ONLY : IGRIB_OPEN_FILE, IGRIB_CLOSE_FILE, IGRIB_RELEASE, &
                        & IGRIB_READ_FROM_FILE, IGRIB_NEW_FROM_MESSAGE, &
                        & JPGRIB_SUCCESS, JPGRIB_BUFFER_TOO_SMALL, &
                        & JPGRIB_END_OF_FILE, JPKSIZE_T
      USE YOWABORT, ONLY : WAM_ABORT
      USE EC_LUN   , ONLY : NULERR

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "sdepthlim.intfb.h"
#include "abort1.intfb.h"
#include "expand_string.intfb.h"
#include "grib2wgrid.intfb.h"
#include "grstname.intfb.h"
#include "init_fieldg.intfb.h"
#include "kgribsize.intfb.h"
#include "mpdistribfl.intfb.h"
#include "readfl.intfb.h"

      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(OUT) :: FL1
      TYPE(WVGRIDGLO), INTENT(IN)                                            :: BLK2GLO
      TYPE(WVGRIDLOC), INTENT(IN)                                            :: BLK2LOC
      TYPE(ENVIRONMENT), INTENT(IN)                                          :: WVENVI
      INTEGER(KIND=JWIM), DIMENSION(NPROC), INTENT(IN)                       :: NBLKS, NBLKE
      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD


      INTEGER(KIND=JWIM) :: NBIT

      INTEGER(KIND=JWIM) :: ISEND
      INTEGER(KIND=JWIM) :: IJ, K, M, IC, ICC, JSN, IDUM, IX, IY, ID, MR, KR, IP
      INTEGER(KIND=JWIM) :: KLOOP, MLOOP, KINF, KSUP, MINF, MSUP
      INTEGER(KIND=JWIM) :: IERR, KRCOUNT, KRTAG
      INTEGER(KIND=JWIM) :: IRA
      INTEGER(KIND=JWIM) :: IUNIT
      INTEGER(KIND=JWIM) :: IFCST
      INTEGER(KIND=JWIM) :: LNAME
      INTEGER(KIND=JWIM) :: IBREAD, NBREAD, NBREAD_AGAIN 
      INTEGER(KIND=JWIM) :: IPARAM, KZLEV, KK, MM
      INTEGER(KIND=JWIM) :: IYYYY, JCONS, IFORP, KDEXN
      INTEGER(KIND=JWIM) :: ISTEP, ISTEP_LOCAL
      INTEGER(KIND=JWIM) :: KRET, IPLENG, ISIZE, KLEN, ILENG, ILENB, KWORD
      INTEGER(KIND=JWIM) :: IRET
      INTEGER(KIND=JWIM) :: LFILE, KFILE_HANDLE, KGRIB_HANDLE
      INTEGER(KIND=JWIM) :: JKGLO
      INTEGER(KIND=JWIM) :: IJSG, IJLG, IJSB, IJLB, KIJS, KIJL, IPRM, ICHNK
      INTEGER(KIND=JWIM) :: IPROC, ITAG, ISREQ, IRREQ, JNR, INR, IST, IEND, KSEND 
      INTEGER(KIND=JWIM) :: ISENDREQ(NPROC), IRECVREQ(NPROC), KFROM(NPROC)
      INTEGER(KIND=JWIM) :: NLONRGG_LOC(NGY)
      INTEGER(KIND=JWIM), ALLOCATABLE :: INGRIB(:), INTMP(:)
      INTEGER(KIND=JPKSIZE_T) :: KBYTES

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: RMONOP
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:) :: EM
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:) :: WORK
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:) :: ZRECVBUF
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:,:) :: RFL
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:) :: FIELD
      TYPE(FORCING_FIELDS) :: FIELDG

      CHARACTER(LEN= 14) :: CDATE 
      CHARACTER(LEN=296) :: FILENAME

      LOGICAL :: LLINIALL, LLOCAL
      LOGICAL :: LFRSDECODE, LOUNIT, LCUNIT, LLEXIST
      LOGICAL :: LLRESIZING=.FALSE.
      LOGICAL :: LLEPSMIN=.TRUE.
      LOGICAL :: LLCHKINT

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GETSPEC',0,ZHOOK_HANDLE)

      WRITE(IU06,*)' GETSPEC :'
      CALL FLUSH(IU06)

      LFRSDECODE=.TRUE.

      NBIT=NIBLO+200

      LOUNIT = .TRUE.
      LCUNIT = .TRUE.
      ISEND=IREAD

      IF (LNSESTART .AND. .NOT.LRESTARTED) THEN
!     BY-PASSED INPUT BY STARTING WITH SPECTRA AT NOISE LEVEL
!     =======================================================

!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, M, K, IJ)
        DO ICHNK = 1, NCHNK
          DO M = 1, NFRE
            DO K = 1, NANG
              DO IJ = 1, NPROMA_WAM
                FL1(IJ, K, M, ICHNK) = EPSMIN 
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO

      ELSEIF (LGRIBIN .AND. .NOT.LRESTARTED) THEN
!     INPUT SPECTRA ARE IN GRIB FORMAT
!     ================================

        IF (IRANK == IREAD) THEN
          FILENAME='specwavein'
          LFILE=LEN_TRIM(FILENAME)

          INQUIRE(FILE=FILENAME(1:LFILE), EXIST=LLEXIST)
          IF (.NOT.LLEXIST) THEN
            WRITE(IU06,*)'**************************************'
            WRITE(IU06,*)'*                                    *'
            WRITE(IU06,*)'*GETSPEC : GRIB SPECTRA NOT FOUND IN *'
            WRITE(IU06,*)  FILENAME
            WRITE(IU06,*)'*PROGRAM WILL ABORT                  *'
            WRITE(IU06,*)'*                                    *'
            WRITE(IU06,*)'**************************************'
            WRITE(NULERR,*)'**************************************'
            WRITE(NULERR,*)'*                                    *'
            WRITE(NULERR,*)'*GETSPEC : GRIB SPECTRA NOT FOUND IN *'
            WRITE(NULERR,*)  FILENAME
            WRITE(NULERR,*)'*PROGRAM WILL ABORT                  *'
            WRITE(NULERR,*)'*                                    *'
            WRITE(NULERR,*)'**************************************'
            CALL ABORT1
          ENDIF
        ENDIF

        NBREAD=0
        NBREAD_AGAIN=0
1121    CONTINUE

!         CONNECT INPUT PE (IREAD) WITH INPUT FILE
        IF (IRANK == IREAD) THEN
          CALL IGRIB_OPEN_FILE(KFILE_HANDLE, FILENAME(1:LFILE),'r')
        ENDIF
        IF (.NOT.ALLOCATED(WORK)) ALLOCATE(WORK(NIBLO))

!       GET GRIB DATA FROM (NFRE_RED*NANG) FIELDS

!       READ NPROC-1 FIELDS AND SEND THEM SUCCESSIVELY TO ALL OTHER PE'S
!       FOR DECODING (IF IN MESSAGE PASSING MODE)

        ISIZE=NBIT
        ISTEP=MAX(NPROC-1,1)

        IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
          NLONRGG_LOC(:)=MNP
#else
          CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
        ELSE
          NLONRGG_LOC(:)=NLONRGG(:)
        ENDIF

!       GRIB2WGRID REQUIRES FIELDG !

        LLINIALL = .FALSE.
        LLOCAL = .FALSE.

        IF(.NOT. FIELDG%LALLOC) CALL FIELDG%ALLOC(LBOUNDS=[NXFFS, NYFFS], UBOUNDS=[NXFFE, NYFFE])
        CALL INIT_FIELDG(BLK2LOC, LLINIALL, LLOCAL,         &
     &                   NXFFS, NXFFE, NYFFS, NYFFE, FIELDG)



        ALL_FILE: DO IC=1,NFRE_RED*NANG,ISTEP 

          ISTEP_LOCAL=ISTEP
          IF (IC+ISTEP > NFRE_RED*NANG ) ISTEP_LOCAL=NFRE_RED*NANG-IC+1

          ALL_DECODE_PE : DO IDUM=1,ISTEP_LOCAL
            IF (NPROC == 1) THEN
              KSEND=1
            ELSEIF (IDUM < IREAD) THEN
              KSEND=IDUM
            ELSE
              KSEND=IDUM+1
            ENDIF
            M = (((IC-1)+IDUM-1)/NANG)+1
            K = (IC-1)+IDUM-(M-1)*NANG

!           DATA ARE READ IN ON PE IREAD

!           LOAD THE DATA
            IF (IRANK == IREAD) THEN
1021          ISIZE=NBIT
              KBYTES=ISIZE*NPRECI
              IF (.NOT.ALLOCATED(INGRIB)) ALLOCATE(INGRIB(ISIZE))
              NBREAD=NBREAD+1
              CALL IGRIB_READ_FROM_FILE(KFILE_HANDLE,INGRIB,KBYTES,IRET)
              IF (IRET == JPGRIB_BUFFER_TOO_SMALL) THEN
                IF (.NOT.LLRESIZING) NBREAD_AGAIN=NBREAD
                CALL KGRIBSIZE(IU06, KBYTES, NBIT, 'GETSPEC')
                DEALLOCATE(INGRIB)
                LLRESIZING=.TRUE.
                GOTO 1021
              ELSEIF (LLRESIZING .AND. IRET /= JPGRIB_END_OF_FILE) THEN
!               LOOP UNTIL YOU HAVE EXPLORE THE SIZE FOR THE WHOLE FILE.
                DEALLOCATE(INGRIB)
                GOTO 1021
              ELSEIF (LLRESIZING .AND. IRET == JPGRIB_END_OF_FILE) THEN
!               WE SHOULD HAVE THE MAXIMUM SIZE NECESSARY, START ALL OVER.
                WRITE(IU06,*) ''
                WRITE(IU06,*) '* GETSPEC: WE SHOULD HAVE THE MAXIMUM SIZE NECESSARY, START ALL OVER.'
                WRITE(IU06,*) ''
                CALL FLUSH(IU06)
                DEALLOCATE(INGRIB)
                LLRESIZING=.FALSE.
                CALL IGRIB_CLOSE_FILE(KFILE_HANDLE)
                CALL IGRIB_OPEN_FILE(KFILE_HANDLE,FILENAME(1:LFILE),'r')
                ISIZE=NBIT
                IF (.NOT.ALLOCATED(INGRIB)) ALLOCATE(INGRIB(ISIZE))
!               READ AGAIN UNTIL THE FIRST TIME WE ENCOUNTERED JPGRIB_BUFFER_TOO_SMALL 
                DO IBREAD=1,NBREAD_AGAIN
                  KBYTES=ISIZE*NPRECI
                  CALL IGRIB_READ_FROM_FILE(KFILE_HANDLE,INGRIB,KBYTES,IRET)
                  IF (IRET == JPGRIB_BUFFER_TOO_SMALL) THEN
                    WRITE(IU06,*) '****************************************************'
                    WRITE(IU06,*) '* GETSPEC: JPGRIB_BUFFER_TOO_SMALL SHOULD NOT HAPPEN'
                    WRITE(NULERR,*) '* GETSPEC: JPGRIB_BUFFER_TOO_SMALL SHOULD NOT HAPPEN'
                    WRITE(IU06,*) '****************************************************'
                    CALL ABORT1
                  ELSEIF (IRET == JPGRIB_END_OF_FILE) THEN
                    WRITE(IU06,*) '**********************************'
                    WRITE(IU06,*) '* GETSPEC: END OF FILE ENCOUNTED'
                    WRITE(NULERR,*) '* GETSPEC: END OF FILE ENCOUNTED'
                    WRITE(IU06,*) '**********************************'
                    CALL ABORT1
                  ELSEIF (IRET /= JPGRIB_SUCCESS) THEN
                    WRITE(IU06,*) '**********************************'
                    WRITE(IU06,*) '* GETSPEC: FILE HANDLING ERROR'
                    WRITE(NULERR,*) '* GETSPEC: FILE HANDLING ERROR'
                    WRITE(IU06,*) '**********************************'
                    CALL ABORT1
                  ENDIF
                ENDDO
                NBREAD=IBREAD-1
                NBREAD_AGAIN=0

              ENDIF

            ENDIF

!           IN CASE OF MESSAGE PASSING THE DECODING WILL OCCUR ON KSEND
!           SEND DATA TO KEND: 
            CALL GSTATS(623,0)
            IF (IRANK == KSEND .AND. NPROC /= 1) THEN
!             RECEIVE GRIB DATA SIZE FROM IREAD
              ITAG=(M-1)*NANG+K
              CALL MPL_RECV(ISIZE,KSOURCE=IREAD,KTAG=ITAG,              &
     &          KOUNT=KRCOUNT,KRECVTAG=KRTAG,KERROR=IERR,               &
     &          CDSTRING='GETSPEC 0:')
              IF (IERR < 0) CALL MPL_ABORT('MPL_RECV ERROR AT 1 in GETSPEC' )
              IF (KRTAG /= ITAG) CALL MPL_ABORT                          &
     &          ('MPL_RECV ERROR AT 1 in GETSPEC:  MISMATCHED TAGS' )

              IF (ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)
              ALLOCATE(INGRIB(ISIZE))
            ELSE IF (IRANK == IREAD .AND. NPROC /= 1) THEN
!             SEND GRIB DATA SIZE TO PE KSEND
              ITAG=(M-1)*NANG+K
              CALL MPL_SEND(ISIZE,KDEST=KSEND,KTAG=ITAG,KERROR=IERR,CDSTRING='GETSPEC 0:')
              IF (IERR < 0) CALL MPL_ABORT('MPL_SEND ERROR AT 1 in GETSPEC' )
            ENDIF

            IF (IRANK == KSEND .AND. NPROC /= 1) THEN
!             RECEIVED GRIB DATA FROM PE IREAD
              ITAG=NFRE_RED*NANG+(M-1)*NANG+K
              ALLOCATE(INTMP(1:ISIZE))
              CALL MPL_RECV(INTMP(1:ISIZE),KSOURCE=IREAD,KTAG=ITAG,     &
     &             KOUNT=KRCOUNT,KRECVTAG=KRTAG,KERROR=IERR,            &
     &             CDSTRING='GETSPEC 1:')
              IF (IERR < 0) CALL MPL_ABORT                              &
     &                     ('MPL_RECV ERROR AT 2 in GETSPEC ' )
              IF (KRCOUNT /= ISIZE) CALL MPL_ABORT                       &
     &        ('MPL_RECV ERROR in 2 in GETSPEC:MISMATCHED MSG LENGTH')
              IF (KRTAG /= ITAG) CALL MPL_ABORT                          &
     &        ('MPL_RECV ERROR in 2 in GETSPEC:MISMATCHED TAGS' )

            ELSE IF (IRANK == IREAD .AND. NPROC /= 1) THEN
!             SEND GRIB DATA TO PE KSEND
              ITAG=NFRE_RED*NANG+(M-1)*NANG+K
              CALL MPL_SEND(INGRIB(1:ISIZE),KDEST=KSEND,KTAG=ITAG,      &
     &         KMP_TYPE=JP_BLOCKING_STANDARD,  &
     &         KERROR=IERR,CDSTRING='GETSPEC 1:')
              IF (IERR < 0) CALL MPL_ABORT('MPL_SEND ERROR AT 2 in GETSPEC' )
              IF (ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)
            ENDIF

            IF (ALLOCATED(INTMP)) THEN
              INGRIB(1:ISIZE) = INTMP(1:ISIZE)
              DEALLOCATE(INTMP)
            ENDIF
            CALL GSTATS(623,1)

!           DECODE THE GRIB DATA ON PE KSEND 
            IF (IRANK == KSEND .OR. NPROC == 1) THEN

              KGRIB_HANDLE=-99
              CALL IGRIB_NEW_FROM_MESSAGE(KGRIB_HANDLE,INGRIB)

              IF (.NOT.ALLOCATED(FIELD)) ALLOCATE(FIELD(NXFFS:NXFFE, NYFFS:NYFFE))

              LLCHKINT = .TRUE.

              CALL GRIB2WGRID (IU06, NPROMA_WAM,                           &
     &                         KGRIB_HANDLE, INGRIB, ISIZE,                &
     &                         LLUNSTR, LLCHKINT,                          &
     &                         NGY, IRGG, NLONRGG_LOC,                     &
     &                         NXFFS, NXFFE, NYFFS, NYFFE,                 &
     &                         FIELDG%XLON, FIELDG%YLAT,                   &
     &                         ZMISS, PPREC, PPEPS,                        &
     &                         CDATE, IFORP, IPARAM, KZLEV, KK, MM, FIELD)

              CALL IGRIB_RELEASE(KGRIB_HANDLE)

              IF (CDATE /= CDTPRO) THEN
                WRITE(NULERR,*)'**********************************'
                WRITE(NULERR,*)'*                                *'
                WRITE(NULERR,*)'* FATAL ERROR IN SUB GETSPEC     *'
                WRITE(NULERR,*)'* ===========================    *'
                WRITE(NULERR,*)'*                                *'
                WRITE(NULERR,*)'* REQUESTED DATE IS NOT EQUAL TO *'
                WRITE(NULERR,*)'* RETRIEVED DATE.                *'
                WRITE(NULERR,*)'* IN FILE: ',FILENAME
                WRITE(NULERR,*)'* CDATE = ',CDATE
                WRITE(NULERR,*)'* CDTPRO = ',CDTPRO
                WRITE(NULERR,*)'*                                *'
                WRITE(NULERR,*)'**********************************'
                CALL ABORT1
              ENDIF
              IF (K /= KK) THEN
                WRITE(NULERR,*) '************************************'
                WRITE(NULERR,*) '* FATAL ERROR IN SUB. GETSPEC      *'
                WRITE(NULERR,*) '* REQUESTED AND DECODED DIRECTIONAL*'
                WRITE(NULERR,*) '* INDEX ARE DIFFERENT :            *'
                WRITE(NULERR,*) '* REQUESTED : ',K 
                WRITE(NULERR,*) '* DECODED   : ',KK
                WRITE(NULERR,*) '*                                  *'
                WRITE(NULERR,*) '************************************'
                CALL ABORT1
              ENDIF
              IF (M /= MM) THEN
                WRITE(NULERR,*) '************************************'
                WRITE(NULERR,*) '* FATAL ERROR IN SUB. GETSPEC      *'
                WRITE(NULERR,*) '* REQUESTED AND DECODED FREQUENCY  *'
                WRITE(NULERR,*) '* INDEX ARE DIFFERENT :            *'
                WRITE(NULERR,*) '* REQUESTED : ',M 
                WRITE(NULERR,*) '* DECODED   : ',MM
                WRITE(NULERR,*) '*                                  *'
                WRITE(NULERR,*) '************************************'
                CALL ABORT1
              ENDIF

!$OMP         PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO, KIJS, KIJL, IJ, IX, IY)
              DO JKGLO = 1, NIBLO, NPROMA_WAM
                KIJS=JKGLO
                KIJL=MIN(KIJS+NPROMA_WAM-1, NIBLO)
                DO IJ = KIJS, KIJL
                  IX = BLK2GLO%IXLG(IJ)
                  IY = NGY- BLK2GLO%KXLT(IJ) +1
                  WORK(IJ) = FIELD(IX,IY)
                ENDDO
              ENDDO
!$OMP         END PARALLEL DO


              DEALLOCATE(FIELD)

            ENDIF ! end decode on KSEND

          ENDDO ALL_DECODE_PE

          IF (LLUNSTR) THEN
              WRITE(*,*) 'In GETSPEC : NOT YET READY!'
              CALL ABORT1
          ELSE
            CALL GSTATS(623,0)
!           RECEIVE THE RESPECTIVE CONTRIBUTIONS OF WORK TO EACH PE.
            IPROC=IRANK
            IST=NBLKS(IPROC)
            IEND=NBLKE(IPROC)
            ILENB=(IEND-IST)+1
            IF( ALLOCATED(ZRECVBUF)) THEN
               IF( SIZE(ZRECVBUF) /= (ILENB*ISTEP_LOCAL) ) THEN
                  DEALLOCATE(ZRECVBUF)
                  ALLOCATE(ZRECVBUF(IST:IEND,ISTEP_LOCAL))
               ENDIF
            ELSE
              ALLOCATE(ZRECVBUF(IST:IEND,ISTEP_LOCAL))
            ENDIF
            IRREQ=0
            DO IDUM = 1, ISTEP_LOCAL
              IF (NPROC == 1) THEN
                KSEND=1
              ELSEIF (IDUM < IREAD) THEN
                KSEND=IDUM
              ELSE
                KSEND=IDUM+1
              ENDIF

              M=(((IC-1)+IDUM-1)/NANG)+1
              K=(IC-1)+IDUM-(M-1)*NANG
              
              IF ( IRANK /= KSEND ) THEN
!               RECEIVE INFORMATION FROM KSEND (that sets M and K)
                IRREQ=IRREQ+1
                KFROM(IRREQ)=KSEND
                IF (KFROM(IRREQ) < IREAD) THEN
                  ID=KFROM(IRREQ)
                ELSE
                  ID=KFROM(IRREQ)-1
                ENDIF
                MR = (((IC-1)+ID-1)/NANG)+1
                KR = (IC-1)+ID-(MR-1)*NANG

                ITAG = 2*NFRE_RED*NANG+(MR-1)*NANG+KR

                CALL MPL_RECV(ZRECVBUF(IST:IEND,IRREQ),                 &
     &                    KSOURCE=KSEND,KREQUEST=IRECVREQ(IRREQ),      &
     &                    KMP_TYPE=JP_NON_BLOCKING_STANDARD,KTAG=ITAG, &
     &                    CDSTRING='GETSPEC: RECEIVING WORK' )
              ELSE
!               SAVE LOCAL CONTRIBUTION

                IF (NBLKS(IRANK) /= IJFROMCHNK(1,1) .OR. NBLKE(IRANK) /= IJFROMCHNK(KIJL4CHNK(NCHNK), NCHNK) ) THEN
                  WRITE(IU06,*)'* GETSPEC : SERIOUS ISSUE WITH THE MODEL DECOMPOSITION FOR THE LOCAL PTS *'
                  WRITE(NULERR,*)'*************************************************************************'
                  WRITE(NULERR,*)'* IRANK = ',IRANK
                  WRITE(NULERR,*)'* GETSPEC : SERIOUS ISSUE WITH THE MODEL DECOMPOSITION FOR THE LOCAL PTS *'
                  WRITE(NULERR,*)'* THE FOLLOWING TWO NUMBERS SHOULD BE EQUAL !!!'
                  WRITE(NULERR,*)'* NBLKS(IRANK) = ', NBLKS(IRANK)
                  WRITE(NULERR,*)'* IJFROMCHNK(1,1) = ',IJFROMCHNK(1,1)
                  WRITE(NULERR,*)'* AND OR THE FOLLOWING TWO NUMBERS SHOULD BE EQUAL !!!'
                  WRITE(NULERR,*)'* NBLKE(IRANK) = ', NBLKE(IRANK)
                  WRITE(NULERR,*)'* IJFROMCHNK(KIJL4CHNK(NCHNK), NCHNK) = ', IJFROMCHNK(KIJL4CHNK(NCHNK), NCHNK)
                  WRITE(NULERR,*)'*                                                                       *'
                  WRITE(NULERR,*)'*************************************************************************'
                  CALL ABORT1
                ENDIF

!$OMP           PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, IJSB, KIJL, IJLB, IJ) 
                DO ICHNK = 1, NCHNK
                  KIJS = 1
                  IJSB = IJFROMCHNK(KIJS, ICHNK)
                  KIJL = KIJL4CHNK(ICHNK)
                  IJLB = IJFROMCHNK(KIJL, ICHNK)

                  FL1(KIJS:KIJL, K, M, ICHNK) = WORK(IJSB:IJLB)

                  DO IJ = KIJS, KIJL
                    IF (FL1(IJ, K, M, ICHNK) ==  ZMISS) FL1(IJ, K, M, ICHNK) = EPSMIN 
                  ENDDO

                  IF (KIJL < NPROMA_WAM) THEN
                    FL1(KIJL+1:NPROMA_WAM, K, M, ICHNK) = FL1(1, K, M, ICHNK)
                  ENDIF
                ENDDO
!$OMP           END PARALLEL DO
              ENDIF
            ENDDO
            
            ISREQ=0
            DO IDUM=1,ISTEP_LOCAL
              IF (NPROC == 1) THEN
                KSEND=1
              ELSEIF (IDUM < IREAD) THEN
                KSEND=IDUM
              ELSE
                KSEND=IDUM+1
              ENDIF
              M=(((IC-1)+IDUM-1)/NANG)+1
              K=(IC-1)+IDUM-(M-1)*NANG

              ITAG=2*NFRE_RED*NANG+(M-1)*NANG+K

              IF (IRANK == KSEND) THEN
!               SEND TO ALL OTHER TASKS
                DO IP=1,NPROC-1
                  IPROC=MOD(IRANK+IP-1,NPROC)+1
                  ISREQ=ISREQ+1
                  IST=NBLKS(IPROC)
                  IEND=NBLKE(IPROC)
                  CALL MPL_SEND(WORK(IST:IEND),                           &
     &                          KDEST=IPROC,KTAG=ITAG,                    &
     &                          KMP_TYPE=JP_NON_BLOCKING_STANDARD,        &
     &                          KREQUEST=ISENDREQ(ISREQ),                 &
     &                          CDSTRING='GETSPEC: SENDING WORK' )
                ENDDO
              ENDIF
            ENDDO

            DO JNR=1,IRREQ
              CALL MPL_WAITANY(KREQUEST=IRECVREQ(1:IRREQ),KINDEX=INR,CDSTRING='GETSPEC: WAIT FOR ANY RECEIVES')
              IF (KFROM(INR) < IREAD) THEN
                ID=KFROM(INR)
              ELSE
                ID=KFROM(INR)-1
              ENDIF
              MR = (((IC-1)+ID-1)/NANG)+1
              KR = (IC-1)+ID-(MR-1)*NANG

!$OMP         PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, IJSB, KIJL, IJLB) 
              DO ICHNK = 1, NCHNK
                KIJS = 1
                IJSB = IJFROMCHNK(KIJS, ICHNK)
                KIJL = KIJL4CHNK(ICHNK)
                IJLB = IJFROMCHNK(KIJL, ICHNK)

                FL1(KIJS:KIJL, KR, MR, ICHNK) = ZRECVBUF(IJSB:IJLB,INR)

                DO IJ = KIJS, KIJL
                  IF (FL1(IJ, KR, MR, ICHNK) ==  ZMISS) FL1(IJ, KR, MR, ICHNK) = EPSMIN 
                ENDDO

                IF (KIJL < NPROMA_WAM) THEN
                  FL1(KIJL+1:NPROMA_WAM, KR, MR, ICHNK) = FL1(1, KR, MR, ICHNK)
                ENDIF
              ENDDO
!$OMP         END PARALLEL DO
            ENDDO   

!           ENSURE ALL SENDS ARE FINISHED.
            IF (ISREQ > 0) THEN
              CALL MPL_WAIT(KREQUEST=ISENDREQ(1:ISREQ), CDSTRING='GETSPEC: WAIT FOR SENDS')
            ENDIF
            CALL GSTATS(623,1)


!           MAKE SURE THAT ALL RECEIVE ARE FINISHED ON ALL PE'S
!           BEFORE PROCESSING ANOTHER BATCH.
            CALL MPL_BARRIER(CDSTRING='GETSPEC:')
          ENDIF  ! end llunstr not ready

        ENDDO ALL_FILE

        IF(ALLOCATED(ZRECVBUF)) DEALLOCATE(ZRECVBUF)
        CALL FIELDG%DEALLOC()

        IF (IRANK == IREAD) THEN
          CALL IGRIB_CLOSE_FILE(KFILE_HANDLE)
        ENDIF

        IF (ALLOCATED(WORK)) DEALLOCATE(WORK)


!       FILL MISSING PART (if any) AND
!       CHECK THAT INPUT SPECTRA ARE CONSISTENT WITH MODEL DEPTH
!       RESCALE IF NOT
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, M, K, IJ)
        DO ICHNK = 1, NCHNK
          DO M = NFRE_RED+1, NFRE
            DO K = 1, NANG
              DO IJ = 1, NPROMA_WAM
                FL1(IJ, K, M, ICHNK) = FL1(IJ, K, NFRE_RED, ICHNK) * FR5(NFRE_RED)*FRM5(M) 
              ENDDO
            ENDDO
          ENDDO

          CALL SDEPTHLIM(1, NPROMA_WAM, WVENVI%EMAXDPT(:,ICHNK), FL1(:,:,:,ICHNK))
        ENDDO
!$OMP   END PARALLEL DO


      ELSE
 
!     BINARY INPUT:
!     =============

         IFCST = 0
         CALL GRSTNAME(CDTPRO, CDATEF, IFCST, 'BLS', ICPLEN, CPATH, FILENAME)

         IUNIT=0

         IF (LRSTPARALR) THEN
!          RESTART FILES FROM ALL PE's
           LNAME = LEN_TRIM(FILENAME)
           FILENAME = FILENAME(1:LNAME)//'.%p_%n'
           CALL EXPAND_STRING(IRANK,NPROC, 0, 0, FILENAME, 1)

           IJSG = IJFROMCHNK(1,1)
           IJLG = IJSG + SUM(KIJL4CHNK) - 1
           ALLOCATE(RFL(IJSG:IJLG, NANG, NFRE))

           CALL READFL(RFL, IJSG, IJLG, 1, NANG, 1, NFRE,              &
     &                 FILENAME, IUNIT, LOUNIT, LCUNIT, LRSTPARALR)

!$OMP      PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, IJSB, KIJL, IJLB, M, K)
           DO ICHNK = 1, NCHNK
             KIJS = 1
             IJSB = IJFROMCHNK(KIJS, ICHNK)
             KIJL = KIJL4CHNK(ICHNK)
             IJLB = IJFROMCHNK(KIJL, ICHNK)

             FL1(KIJS:KIJL, :, :, ICHNK) = RFL(IJSB:IJLB, :, :)

             IF (KIJL < NPROMA_WAM) THEN
                DO M = 1, NFRE 
                  DO K = 1, NANG 
                    FL1(KIJL+1:NPROMA_WAM, K, M, ICHNK) = FL1(1, K, M, ICHNK)
                  ENDDO
                ENDDO
             ENDIF
           ENDDO
!$OMP      END PARALLEL DO

           DEALLOCATE(RFL)

         ELSE

           DO MLOOP= 1, NFRE, MDEL
             MINF=MLOOP
             MSUP=MIN(MLOOP+MDEL-1,NFRE)
             DO KLOOP=1,NANG,KDEL
               KINF=KLOOP
               KSUP=MIN(KLOOP+KDEL-1, NANG)

               ALLOCATE(RFL(1:NIBLO, KINF:KSUP, MINF:MSUP))
!              READ RESTART SPECTRA FROM PE ISEND (IREAD) 
               IF (IRANK == ISEND) THEN
                 LOUNIT = .FALSE.
                 LCUNIT = .FALSE.
                 IF (MINF == 1 .AND. KINF == 1) LOUNIT = .TRUE.
                 IF (MSUP == NFRE .AND. KSUP == NANG) LCUNIT = .TRUE.

                 CALL READFL(RFL, 1, NIBLO, KINF, KSUP, MINF, MSUP,       &
     &                       FILENAME, IUNIT, LOUNIT, LCUNIT, LRSTPARALR)
               ENDIF

               CALL MPDISTRIBFL(ISEND, KTAG, NBLKS, NBLKE, KINF, KSUP, MINF, MSUP, RFL)
               KTAG=KTAG+1

!             KEEP CORRESPONDING CONTRIBUTION TO FL1
!$OMP         PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, IJSB, KIJL, IJLB, K, M)
              DO ICHNK = 1, NCHNK
                KIJS = 1
                IJSB = IJFROMCHNK(KIJS, ICHNK)
                KIJL = KIJL4CHNK(ICHNK)
                IJLB = IJFROMCHNK(KIJL, ICHNK)

                DO M = MINF, MSUP
                  DO K = KINF, KSUP
                    FL1(KIJS:KIJL, K, M, ICHNK) = RFL(IJSB:IJLB, K, M)
                  ENDDO
                ENDDO

                IF (KIJL < NPROMA_WAM) THEN
                  DO M = MINF, MSUP
                    DO K = KINF, KSUP
                      FL1(KIJL+1:NPROMA_WAM, K, M, ICHNK) = FL1(1, K, M, ICHNK)
                    ENDDO
                  ENDDO
                ENDIF

              ENDDO
!$OMP         END PARALLEL DO

               DEALLOCATE(RFL)
             ENDDO
           ENDDO

         ENDIF

      ENDIF

      WRITE(IU06,*) ' SPECTRUM FILE READ IN............... CDTPRO  = ', CDTPRO
      WRITE(IU06,*) ' '
      CALL FLUSH (IU06)

IF (LHOOK) CALL DR_HOOK('GETSPEC',1,ZHOOK_HANDLE)

END SUBROUTINE GETSPEC
