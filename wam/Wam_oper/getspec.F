      SUBROUTINE GETSPEC(FL,SL,NSTART,NEND,IREAD,U10OLD,THWOLD)
! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      SEPTEMBER 1997 
!     J. BIDLOT    ECMWF      MARCH 2010: modified to use gribapi 

!*    PURPOSE.
!     --------
!     READ THE SPECTRA FROM DISK.

!**   INTERFACE.
!     ----------
!     *CALL* *GETSPEC(FL,SL,NSTART,NEND,IREAD,U10OLD,THWOLD)
!     *FL*        ARRAY CONTAINING THE SPECTRA CONTRIBUTION ON EACH PE
!     *SL*        SPARE ARRAY TO MANIPUATE SPECTRA
!     *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!     *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!     *IREAD*     PROCESSOR WHICH WILL ACCESS THE FILE ON DISK 
!     *U10OLD*    WIND SPEED. (used with fetch law to fill empty 
!                 sea points)
!     *THWOLD*    WIND DIRECTION (RADIANS).


!     METHOD.
!     -------

!     IN CASE THE INPUT SPECTRA ARE IN GRIB THEN EITHER THEY ARE READ
!     FROM A FILE ON DISK USING THE PBGRIB SOFTWARE.
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

      USE YOWCOUT  , ONLY : KDEL     ,MDEL
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRIBHD, ONLY : PPEPS    ,PPREC
      USE YOWGRID  , ONLY : DELPHI   ,NLONRGG
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,IRGG     ,AMOSOP   ,
     &            AMONOP   ,AMOWEP   ,AMOEAP   ,XDELLA   ,ZDELLO

      USE YOWMESPAS, ONLY : LMESSPASS,LPBIOIN  ,LFDBIOIN ,LGRIBIN
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP     ,
     &            KTAG     ,NPRECR   ,NPRECI
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      ,
     &            NBLO     ,NIBLO    ,CLDOMAIN
      USE YOWPCONS , ONLY : G        ,DEG      ,R        ,ZMISS    ,
     &            EPSMIN
      USE YOWSTAT  , ONLY : CDATEF   ,CDTPRO   ,NWAM_BLKS,IREFRA   ,
     &            LNSESTART
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH    ,LRESTARTED
      USE MPL_MODULE
      USE FDBSUBS_MOD
      USE GRIB_API_INTERFACE

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: NBIT

      INTEGER :: ISEND, IREAD
      INTEGER :: IJ, K, M, IC, ICC, JSN, IDUM, IX, IY, ID, MR, KR, IP
      INTEGER :: KLOOP, MLOOP, KINF, KSUP, MINF, MSUP
      INTEGER :: IERR, KRCOUNT, KFROM, KRTAG
      INTEGER :: IRA, IG
      INTEGER :: IUNIT
      INTEGER :: IBREAD, NBREAD, NBREAD_AGAIN 
      INTEGER :: IPARAM, KZLEV,KK,MM
      INTEGER :: IYYYY, JCONS, IFORP, KDEXN
      INTEGER :: ISTEP, ISTEP_LOCAL
      INTEGER :: KRET, IPLENG, ISIZE, KLEN, ILENG, KWORD
      INTEGER :: IRET
      INTEGER :: LFILE, KFILE_HANDLE, KGRIB_HANDLE
      INTEGER :: JKGLO,KIJS,KIJL,NPROMA
      INTEGER :: IPROC, ITAG, IREQ, IST, IEND, KSEND
      INTEGER :: ISENDREQ(NPROC)
      INTEGER, DIMENSION(NPROC) :: NSTART, NEND
      INTEGER, ALLOCATABLE :: INGRIB(:), INTMP(:)
      INTEGER(KIND=JPKSIZE_T) :: KBYTES

      REAL :: ZTHETA,ZFRE
      REAL :: RMONOP
      REAL, DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL, SL
      REAL, DIMENSION(NINF:NSUP,NBLO) :: U10OLD,THWOLD
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK
      REAL, ALLOCATABLE, DIMENSION(:) :: ZRECVBUF
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: RFL
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FIELD

      CHARACTER(LEN=14) :: CDATE 
      CHARACTER(LEN=120) :: FILENAME

      LOGICAL :: LLEPSMIN
      LOGICAL :: LFRSDECODE, LPBOPEN, LPBCLOSE, LLEXIST
      LOGICAL :: LLRESIZING=.FALSE.

! ----------------------------------------------------------------------

      LFRSDECODE=.TRUE.

      NBIT=NIBLO+200

      IF (ITEST > 3 ) THEN
        WRITE(IU06,*) ' SUB. GETSPEC:'
        WRITE(IU06,*) ' ABOUT TO READ SPECTRUM IN........'
        WRITE(IU06,*) ' CDTPRO    =', CDTPRO
        WRITE(IU06,*) ' LGRIBIN   =', LGRIBIN
        WRITE(IU06,*) ' LPBIOIN   =', LPBIOIN
        WRITE(IU06,*) ' LFDBIOIN  =', LFDBIOIN
        WRITE(IU06,*) ' LRESTARTED=', LRESTARTED
        WRITE(IU06,*) ' IRANK=', IRANK
        WRITE(IU06,*) ' IREAD=', IREAD
        CALL FLUSH(IU06)
      ENDIF
C
      LPBOPEN = .TRUE.
      LPBCLOSE = .TRUE.
      ISEND=IREAD

      IF(LNSESTART .AND. .NOT.LRESTARTED) THEN
!     BY-PASSED INPUT BY STARTING WITH SPECTRA AT NOISE LEVEL
!     =======================================================

        NPROMA=(NEND(IRANK)-NSTART(IRANK))/NWAM_BLKS+1
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ,K,M)
        DO JKGLO=NSTART(IRANK),NEND(IRANK),NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,NEND(IRANK))
          DO IJ=KIJS,KIJL
            DO K=1,NANG
              DO M=1,NFRE
                FL(IJ,K,M) = EPSMIN 
              ENDDO
           ENDDO
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO

         DO K=1,NANG
           DO M=1,NFRE
             FL(NINF-1,K,M) = 0. 
           ENDDO
         ENDDO

      ELSEIF(LGRIBIN.AND..NOT.LRESTARTED) THEN
!     INPUT SPECTRA ARE IN GRIB FORMAT
!     ================================

        IF (LFDBIOIN) THEN
          WRITE(IU06,*)'**************************************'
          WRITE(IU06,*)'*                                    *'
          WRITE(IU06,*)'*GETSPEC : THIS OPTION IS OBSOLETE !!*'
          WRITE(IU06,*)'*PROGRAM WILL ABORT                  *'
          WRITE(IU06,*)'*                                    *'
          WRITE(IU06,*)'**************************************'
          CALL ABORT1
        ENDIF

        FILENAME='specwavein'
        LFILE=LEN_TRIM(FILENAME)

        INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)
        IF(.NOT.LLEXIST) THEN
          WRITE(IU06,*)'**************************************'
          WRITE(IU06,*)'*                                    *'
          WRITE(IU06,*)'*GETSPEC : GRIB SPECTRA NOT FOUND IN *'
          WRITE(IU06,*)  FILENAME
          WRITE(IU06,*)'*PROGRAM WILL ABORT                  *'
          WRITE(IU06,*)'*                                    *'
          WRITE(IU06,*)'**************************************'
          CALL ABORT1
        ENDIF

        NBREAD=0
        NBREAD_AGAIN=0
1121    CONTINUE

!         CONNECT INPUT PE (IREAD) WITH INPUT FILE
        IF (IRANK.EQ.IREAD) THEN
          CALL IGRIB_OPEN_FILE(KFILE_HANDLE,FILENAME(1:LFILE),'r')
        ENDIF
        IF (.NOT.ALLOCATED(WORK)) ALLOCATE(WORK(NIBLO))

!       GET GRIB DATA FROM (NFRE*NANG) FIELDS

!       READ NPROC-1 FIELDS AND SEND THEM SUCCESSIVELY TO ALL OTHER PE'S
!       FOR DECODING (IF IN MESSAGE PASSING MODE)

        ISIZE=NBIT
        ISTEP=MAX(NPROC-1,1)

        ALL_FILE : DO IC=1,NFRE*NANG,ISTEP 

          ISTEP_LOCAL=ISTEP
          IF(IC+ISTEP.GT.NFRE*NANG ) ISTEP_LOCAL=NFRE*NANG-IC+1

          ALL_DECODE_PE : DO IDUM=1,ISTEP_LOCAL
            IF(NPROC.EQ.1) THEN
              KSEND=1
            ELSEIF(IDUM.LT.IREAD) THEN
              KSEND=IDUM
            ELSE
              KSEND=IDUM+1
            ENDIF
            M=(((IC-1)+IDUM-1)/NANG)+1
            K=(IC-1)+IDUM-(M-1)*NANG

!           DATA ARE READ IN ON PE IREAD

!           LOAD THE DATA
            IF (IRANK.EQ.IREAD) THEN
1021          ISIZE=NBIT
              KBYTES=ISIZE*NPRECI
              IF(.NOT.ALLOCATED(INGRIB)) ALLOCATE(INGRIB(ISIZE))
              NBREAD=NBREAD+1
              CALL IGRIB_READ_FROM_FILE(KFILE_HANDLE,INGRIB,
     &                                  KBYTES,IRET)
              IF(IRET.EQ.JPGRIB_BUFFER_TOO_SMALL) THEN
                IF(.NOT.LLRESIZING) NBREAD_AGAIN=NBREAD
                CALL KGRIBSIZE(IU06, KBYTES, NBIT, 'GETSPEC')
                DEALLOCATE(INGRIB)
                LLRESIZING=.TRUE.
                GOTO 1021
              ELSEIF(LLRESIZING .AND. IRET.NE.JPGRIB_END_OF_FILE) THEN
!               LOOP UNTIL YOU HAVE EXPLORE THE SIZE FOR THE WHOLE FILE.
                DEALLOCATE(INGRIB)
                GOTO 1021
              ELSEIF(LLRESIZING .AND. IRET.EQ.JPGRIB_END_OF_FILE) THEN
!               WE SHOULD HAVE THE MAXIMUM SIZE NECESSARY, START ALL OVER.
                DEALLOCATE(INGRIB)
                LLRESIZING=.FALSE.
                CALL IGRIB_CLOSE_FILE(KFILE_HANDLE)
                CALL IGRIB_OPEN_FILE(KFILE_HANDLE,FILENAME(1:LFILE),'r')
                ISIZE=NBIT
                IF(.NOT.ALLOCATED(INGRIB)) ALLOCATE(INGRIB(ISIZE))

                NBREAD=NBREAD+1

                CALL IGRIB_READ_FROM_FILE(KFILE_HANDLE,INGRIB,
     &                                    KBYTES,IRET)
                IF(IRET.EQ.JPGRIB_BUFFER_TOO_SMALL) THEN
                  IF(.NOT.LLRESIZING) NBREAD_AGAIN=NBREAD
                  CALL KGRIBSIZE(IU06, KBYTES, NBIT, 'GETSPEC')
                  DEALLOCATE(INGRIB)
                  LLRESIZING=.TRUE.
                  GOTO 1021
                ELSEIF(LLRESIZING .AND. IRET.NE.JPGRIB_END_OF_FILE) THEN
!                 LOOP UNTIL YOU HAVE EXPLORE THE SIZE FOR THE WHOLE FILE.
                  DEALLOCATE(INGRIB)
                  GOTO 1021
                ELSEIF(LLRESIZING .AND. IRET.EQ.JPGRIB_END_OF_FILE) THEN
!                 WE SHOULD HAVE THE MAXIMUM SIZE NECESSARY, START ALL OVER.
                  DEALLOCATE(INGRIB)
                  LLRESIZING=.FALSE.
                  CALL IGRIB_CLOSE_FILE(KFILE_HANDLE)
                  CALL IGRIB_OPEN_FILE(KFILE_HANDLE,
     &                                 FILENAME(1:LFILE),'r')
                  ISIZE=NBIT
                  IF(.NOT.ALLOCATED(INGRIB)) ALLOCATE(INGRIB(ISIZE))
                  DO IBREAD=1,NBREAD_AGAIN
                    KBYTES=ISIZE*NPRECI
                    CALL IGRIB_READ_FROM_FILE(KFILE_HANDLE,INGRIB,
     &                                        KBYTES,IRET)
                  ENDDO
                  NBREAD=IBREAD-1
                  NBREAD_AGAIN=0

                ELSEIF(IRET.EQ.JPGRIB_END_OF_FILE) THEN
                  WRITE(IU06,*) '**********************************'
                  WRITE(IU06,*) '* GETSPEC: END OF FILE ENCOUNTED'
                  WRITE(IU06,*) '**********************************'
                  CALL ABORT1
                ELSEIF(IRET.NE.JPGRIB_SUCCESS) THEN
                  WRITE(IU06,*) '**********************************'
                  WRITE(IU06,*) '* GETSPEC: FILE HANDLING ERROR'
                  WRITE(IU06,*) '**********************************'
                  CALL ABORT1
                ENDIF
              ENDIF

            ENDIF

            IF(LMESSPASS) THEN
!             IN CASE OF MESSAGE PASSING THE DECODING WILL OCCUR ON KSEND
!             SEND DATA TO KEND: 
              CALL GSTATS(623,0)
              IF(IRANK.EQ.IREAD.AND.NPROC.NE.1) THEN
!               SEND GRIB DATA SIZE TO PE KSEND
                ITAG=(M-1)*NANG+K
                CALL MPL_SEND(ISIZE,KDEST=KSEND,KTAG=ITAG,KERROR=IERR,
     &            CDSTRING='GETSPEC 0:')
                IF(IERR.LT.0) CALL MPL_ABORT
     &                    ('MPL_SEND ERROR AT 1 in GETSPEC' )
              ENDIF

              IF (IRANK.EQ.KSEND.AND.NPROC.NE.1) THEN
!               RECEIVE GRIB DATA SIZE FROM IREAD
                ITAG=(M-1)*NANG+K
                CALL MPL_RECV(ISIZE,KSOURCE=IREAD,KTAG=ITAG,
     &            KOUNT=KRCOUNT,KRECVTAG=KRTAG,KERROR=IERR,
     &            CDSTRING='GETSPEC 0:')
                IF(IERR.LT.0) CALL MPL_ABORT
     &                           ('MPL_RECV ERROR AT 1 in GETSPEC' )
                IF(KRTAG.NE.ITAG) CALL MPL_ABORT
     &          ('MPL_RECV ERROR AT 1 in GETSPEC:  MISMATCHED TAGS' )

                IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)
                ALLOCATE(INGRIB(ISIZE))
              ENDIF

              IF(IRANK.EQ.IREAD.AND.NPROC.NE.1) THEN
!               SEND GRIB DATA TO PE KSEND
                ITAG=NFRE*NANG+(M-1)*NANG+K
                CALL MPL_SEND(INGRIB(1:ISIZE),KDEST=KSEND,KTAG=ITAG,
     &           KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(1),
     &           KERROR=IERR,CDSTRING='GETSPEC 1:')
                IF(IERR.LT.0) CALL MPL_ABORT
     &                      ('MPL_SEND ERROR AT 2 in GETSPEC' )
              ENDIF

              IF (IRANK.EQ.KSEND.AND.NPROC.NE.1) THEN
!               RECEIVED GRIB DATA FROM PE IREAD
                ITAG=NFRE*NANG+(M-1)*NANG+K
                ALLOCATE(INTMP(1:ISIZE))
                CALL MPL_RECV(INTMP(1:ISIZE),KSOURCE=IREAD,KTAG=ITAG,
     &               KOUNT=KRCOUNT,KRECVTAG=KRTAG,KERROR=IERR,
     &               CDSTRING='GETSPEC 1:')
                IF(IERR.LT.0) CALL MPL_ABORT
     &                       ('MPL_RECV ERROR AT 2 in GETSPEC ' )
                IF(KRCOUNT.NE.ISIZE) CALL MPL_ABORT
     &          ('MPL_RECV ERROR in 2 in GETSPEC:MISMATCHED MSG LENGTH')
                IF(KRTAG.NE.ITAG) CALL MPL_ABORT
     &          ('MPL_RECV ERROR in 2 in GETSPEC:MISMATCHED TAGS' )
              ENDIF

              IF(IRANK.EQ.IREAD.AND.NPROC.NE.1) THEN
                CALL MPL_WAIT(INGRIB(1:ISIZE),KREQUEST=ISENDREQ(1), 
     &           CDSTRING='GETSPEC: WAIT FOR SEND')
                IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)
              ENDIF

              IF (ALLOCATED(INTMP)) THEN
                INGRIB(1:ISIZE) = INTMP(1:ISIZE)
                DEALLOCATE(INTMP)
              ENDIF

              CALL GSTATS(623,1)

            ENDIF


!           DECODE THE GRIB DATA ON PE KSEND 
            IF (IRANK.EQ.KSEND.OR.NPROC.EQ.1) THEN

              KGRIB_HANDLE=-99
              CALL IGRIB_NEW_FROM_MESSAGE(KGRIB_HANDLE,INGRIB)

              IF(.NOT.ALLOCATED(FIELD)) ALLOCATE(FIELD(NGX,NGY))
              CALL GRIB2WGRID (IU06, ITEST, NWAM_BLKS, 
     &                         KGRIB_HANDLE, INGRIB, ISIZE,
     &                         NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                         AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &                         ZMISS, PPREC, PPEPS,
     &                         CDATE, IFORP, IPARAM, KZLEV,KK,MM,FIELD)

              CALL IGRIB_RELEASE(KGRIB_HANDLE)

              IF (CDATE.NE.CDTPRO) THEN
                WRITE(IU06,*)'**********************************'
                WRITE(IU06,*)'*                                *'
                WRITE(IU06,*)'* FATAL ERROR IN SUB GETSPEC     *'
                WRITE(IU06,*)'* ===========================    *'
                WRITE(IU06,*)'*                                *'
                WRITE(IU06,*)'* REQUESTED DATE IS NOT EQUAL TO *'
                WRITE(IU06,*)'* RETRIEVED DATE.                *'
                WRITE(IU06,*)'* IN FILE: ',FILENAME
                WRITE(IU06,*)'* CDATE = ',CDATE
                WRITE(IU06,*)'* CDTPRO = ',CDTPRO
                WRITE(IU06,*)'*                                *'
                WRITE(IU06,*)'**********************************'
                CALL ABORT1
              ENDIF
              IF(K.NE.KK) THEN
                WRITE(IU06,*) '************************************'
                WRITE(IU06,*) '* FATAL ERROR IN SUB. GETSPEC      *'
                WRITE(IU06,*) '* REQUESTED AND DECODED DIRECTIONAL*'
                WRITE(IU06,*) '* INDEX ARE DIFFERENT :            *'
                WRITE(IU06,*) '* REQUESTED : ',K 
                WRITE(IU06,*) '* DECODED   : ',KK
                WRITE(IU06,*) '*                                  *'
                WRITE(IU06,*) '************************************'
                CALL ABORT1
              ENDIF
              IF(M.NE.MM) THEN
                WRITE(IU06,*) '************************************'
                WRITE(IU06,*) '* FATAL ERROR IN SUB. GETSPEC      *'
                WRITE(IU06,*) '* REQUESTED AND DECODED FREQUENCY  *'
                WRITE(IU06,*) '* INDEX ARE DIFFERENT :            *'
                WRITE(IU06,*) '* REQUESTED : ',M 
                WRITE(IU06,*) '* DECODED   : ',MM
                WRITE(IU06,*) '*                                  *'
                WRITE(IU06,*) '************************************'
                CALL ABORT1
              ENDIF

              NPROMA=NEND(NPROC)/NWAM_BLKS+1
!$OMP         PARALLEL DO SCHEDULE(STATIC)
!$OMP+        PRIVATE(JKGLO,KIJS,KIJL,IJ,IX,IY)
              DO JKGLO=1,NEND(NPROC),NPROMA
                KIJS=JKGLO
                KIJL=MIN(KIJS+NPROMA-1,NEND(NPROC))
                DO IJ=KIJS,KIJL
                  IX = IXLG(IJ,1)
                  IY = NGY- KXLT(IJ,1) +1
                  IF (FIELD(IX,IY) .NE. ZMISS) THEN
                    WORK(IJ) = FIELD(IX,IY)
                  ELSE
                    WORK(IJ) = EPSMIN
                  ENDIF
                ENDDO
              ENDDO
!$OMP         END PARALLEL DO
              DEALLOCATE(FIELD)

              ELSE

            ENDIF ! end decode on KSEND

          ENDDO ALL_DECODE_PE

!         SEND THE RESPECTIVE CONTRIBUTIONS OF WORK TO EACH PE.
          CALL GSTATS(623,0)
          IREQ=0
          DO IDUM=1,ISTEP_LOCAL
            IF(NPROC.EQ.1) THEN
              KSEND=1
            ELSEIF(IDUM.LT.IREAD) THEN
              KSEND=IDUM
            ELSE
              KSEND=IDUM+1
            ENDIF
            M=(((IC-1)+IDUM-1)/NANG)+1
            K=(IC-1)+IDUM-(M-1)*NANG

            ITAG=2*NFRE*NANG+(M-1)*NANG+K

            IF(IRANK.EQ.KSEND) THEN
!             SEND TO ALL OTHER TASKS
              DO IP=1,NPROC-1
                IPROC=MOD(IRANK+IP-1,NPROC)+1
                IREQ=IREQ+1
                IST=NSTART(IPROC)
                IEND=NEND(IPROC)
                CALL MPL_SEND(WORK(IST:IEND),
     &                        KDEST=IPROC,KTAG=ITAG,
     &                        KMP_TYPE=JP_NON_BLOCKING_STANDARD,
     &                        KREQUEST=ISENDREQ(IREQ),
     &                        CDSTRING='GETSPEC: SENDING WORK' )
              ENDDO

            ENDIF
          ENDDO

!         RECEIVE THE RESPECTIVE CONTRIBUTIONS OF WORK TO EACH PE.
          IPROC=IRANK
          IST=NSTART(IPROC)
          IEND=NEND(IPROC)
          ALLOCATE(ZRECVBUF(IST:IEND))
          DO IDUM=1,ISTEP_LOCAL
            IF(NPROC.EQ.1) THEN
              KSEND=1
            ELSEIF(IDUM.LT.IREAD) THEN
              KSEND=IDUM
            ELSE
              KSEND=IDUM+1
            ENDIF

            M=(((IC-1)+IDUM-1)/NANG)+1
            K=(IC-1)+IDUM-(M-1)*NANG

            IF(IRANK.EQ.KSEND) THEN
!             SAVE LOCAL CONTRIBUTION
              DO IJ = NSTART(IRANK),NEND(IRANK) 
                FL(IJ,K,M)=WORK(IJ)
              ENDDO
            ELSE
!             RECEIVE INFORMATION FROM KSEND (that sets M and K)
              CALL MPL_RECV(ZRECVBUF(IST:IEND),
     &                      KFROM=KFROM,KRECVTAG=KRTAG,
     &                      KMP_TYPE=JP_BLOCKING_STANDARD,
     &                      CDSTRING='GETSPEC: RECEIVING WORK' )

              IF(KFROM.LT.IREAD) THEN
                ID=KFROM
              ELSE
                ID=KFROM-1
              ENDIF
              MR=(((IC-1)+ID-1)/NANG)+1
              KR=(IC-1)+ID-(MR-1)*NANG
              ITAG=2*NFRE*NANG+(MR-1)*NANG+KR
              IF(KRTAG.NE.ITAG) THEN
                WRITE(0,*)'MPL_RECV ERROR in GETSPEC: MISMATCHED TAGS'
                WRITE(0,*)'IRANK = ',IRANK
                WRITE(0,*)'KFROM = ',KFROM
                WRITE(0,*)'KRTAG = ',KRTAG
                WRITE(0,*)'ITAG  = ',ITAG
                WRITE(0,*)' ABORTING !!!!'
                CALL ABORT1
              ENDIF

              DO IJ = IST,IEND 
                FL(IJ,KR,MR)=ZRECVBUF(IJ) 
              ENDDO

            ENDIF

            FL(NINF-1,K,M)=0.
          ENDDO
          DEALLOCATE(ZRECVBUF)

!         ENSURE ALL SENDS ARE FINISHED.
          IF(IREQ.GT.0) THEN
            CALL MPL_WAIT(WORK(:),KREQUEST=ISENDREQ(1:IREQ),
     &                    CDSTRING='GETSPEC: WAIT SENDING WORK')
          ENDIF
          CALL GSTATS(623,1)


!         MAKE SURE THAT ALL RECEIVE ARE FINISHED ON ALL PE'S
!         BEFORE PROCESSING ANOTHER BATCH.
          CALL MPL_BARRIER(CDSTRING='GETSPEC:')

        ENDDO ALL_FILE

        IF (IRANK.EQ.IREAD) THEN
          CALL IGRIB_CLOSE_FILE(KFILE_HANDLE)
        ENDIF

        IF(ALLOCATED(WORK)) DEALLOCATE(WORK)

      ELSE
 
!     BINARY INPUT:
!     =============
         IF (ITEST.GE.2) THEN
           WRITE(IU06,*)''
           WRITE(IU06,*)' IN GETSPEC :'
           IF(MDEL.EQ.1) THEN
             WRITE(IU06,*)'  ',MDEL,
     &                  ' FREQUENCY IS DISTRIBUTED AT ONE TIME'
           ELSE
             WRITE(IU06,*)'  ',MDEL,
     &                    ' FREQUENCIES ARE DISTRIBUTED EACH TIME'
           ENDIF
         ENDIF

         IF (IRANK.EQ.ISEND)
     &       CALL GRSTNAME(CDTPRO,CDATEF,'BLS',ICPLEN,CPATH,FILENAME)

         IUNIT=0

         DO MLOOP=1,NFRE,MDEL
           MINF=MLOOP
           MSUP=MIN(MLOOP+MDEL-1,NFRE)
           DO KLOOP=1,NANG,KDEL
             KINF=KLOOP
             KSUP=MIN(KLOOP+KDEL-1,NANG)

             ALLOCATE(RFL(0:NIBLO,KINF:KSUP,MINF:MSUP))
!            READ RESTART SPECTRA FROM PE ISEND (IREAD) 
             IF (IRANK.EQ.ISEND) THEN
               LPBOPEN = .FALSE.
               LPBCLOSE = .FALSE.
               IF(MINF.EQ.1 .AND. KINF.EQ.1) LPBOPEN = .TRUE.
               IF(MSUP.EQ.NFRE .AND. KSUP.EQ.NANG) LPBCLOSE = .TRUE.

               CALL READFL(RFL, KINF, KSUP, MINF, MSUP,
     &                     FILENAME, IUNIT, LPBOPEN, LPBCLOSE)
             ENDIF

             CALL MPDISTRIBFL(ISEND,KTAG,NSTART,NEND,KINF,KSUP,
     &                        MINF,MSUP,RFL)
             KTAG=KTAG+1

             IF (ITEST.GE.2)
     &        WRITE(IU06,*)
     &        'SUB. GETSPEC: RESTART SPECTRUM COLLECTED, PART:',MLOOP

!            KEEP CORRESPONDING CONTRIBUTION TO FL
             DO M=MINF,MSUP
               DO K=KINF,KSUP
                 DO IJ=NSTART(IRANK),NEND(IRANK)
                   FL(IJ,K,M) = RFL(IJ,K,M)
                 ENDDO
               ENDDO
             ENDDO

             DO M=MINF,MSUP
               DO K=KINF,KSUP
                 FL(NINF-1,K,M) = 0. 
               ENDDO
             ENDDO

             DEALLOCATE(RFL)
           ENDDO
         ENDDO

      ENDIF

      WRITE(IU06,*) ' SPECTRUM FILE READ IN...............',
     & ' CDTPRO  = ', CDTPRO
      WRITE(IU06,*) ' '
      CALL FLUSH (IU06)

      RETURN
      END SUBROUTINE GETSPEC
