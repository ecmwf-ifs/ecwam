      SUBROUTINE OUTWSPEC (SPEC)

!----------------------------------------------------------------------

!**** *OUTWSPEC*  ENCODES SPECTRA AS PARAMETER 251 USING GRIB API
!                 AND WRITES TO FILE OR TO FDB.

!     J. BIDLOT   ECMWF  APRIL 2010 

!*    PURPOSE.
!     --------

!       ENCODES SPECTRA AND WRITES TO FDB OR TO A FILE.

!**   INTERFACE.
!     ----------

!     SUBROUTINE OUTWSPEC (SPEC)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *SPEC*     REAL      LOCAL SPECTRA OF CURRENT PE.

!     METHOD.
!     -------

!           ENCODE SPECTRA PER FREQUENCY AND DIRECTION
!      INTO GRIB AND WRITE TO FDB OR TO A SINGLE FILE.

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!       NONE.

!-------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWCOUT  , ONLY : NWRTOUTWAM
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,IRGG     ,AMOWEP   ,
     &            AMOSOP   ,AMOEAP   ,AMONOP   ,XDELLA   ,XDELLO
      USE YOWMESPAS, ONLY : LFDBIOOUT
      USE YOWMPP   , ONLY : NPRECI   ,IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSTAT  , ONLY : CDATEE   ,CDATEF   ,CDTPRO   ,CDATEA   ,
     &            MARSTYPE ,NWAM_BLKS,NWFDBREF ,CFDB2DSP ,LFDBOPEN
      USE YOWSPEC, ONLY   : NSTART   ,NEND
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH
      USE MPL_MODULE
      USE FDBSUBS_MOD
      USE GRIB_API_INTERFACE
      USE YOMHOOK  , ONLY : LHOOK, DR_HOOK

!-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: IFCST
      INTEGER :: IC, IST, ISTEP, IM, IK, J1, ICNT, IP, IR
      INTEGER :: IJ, IX, IY, I, J, K
      INTEGER :: JKGLO, KIJS, KIJL, NPROMA
      INTEGER :: IPARAM, ITABLE
      INTEGER :: IGLOBAL, ILOCAL
      INTEGER :: ITAG, LFILE, IUOUT
      INTEGER :: NRCV
      INTEGER :: IGRIB_HANDLE, IGRBHNDL
      INTEGER :: MAXMSGSIZE
      INTEGER :: IERR, KRCOUNT, KRTAG, KFROM
      INTEGER(KIND=JPKSIZE_T) :: KBYTES
      INTEGER, DIMENSION(NPROC) :: IPR, ISENDCOUNTS, IRECVCOUNTS
      INTEGER, DIMENSION(NPROC) :: MSGSIZE
      INTEGER, DIMENSION(NPROC) :: ISENDBUF, IRECVBUF
      INTEGER, DIMENSION(NPROC) :: ISENDREQ
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ISENDMSG, IRECVMSG

      REAL, DIMENSION(NINF-1:NSUP, NANG, NFRE), INTENT(IN) :: SPEC
      REAL :: ZHOOK_HANDLE
      REAL :: FIELD(NGX, NGY)
      REAL, ALLOCATABLE, DIMENSION(:) :: ZSENDBUF,ZRECVBUF

      CHARACTER(LEN=14) :: CDATE 
      CHARACTER(LEN=120) :: OUTFILEN

!-----------------------------------------------------------------------

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTWSPEC',0,ZHOOK_HANDLE)
#endif

      ILOCAL=0

      IF (ITEST.GT.1) THEN
        WRITE(IU06,*) '*      THIS IS OUTWSPEC         *'
        CALL FLUSH (IU06)
      ENDIF

      IF(CDTPRO.LE.CDATEF) THEN
!*    0.1.  THIS IS AN ANALYSIS DATE.
        IF (LWCOU .AND. MARSTYPE .EQ. 'fg') THEN
          CDATE=CDATEA
          CALL DIFDATE (CDATEA, CDTPRO, IFCST)
          IFCST = IFCST/3600
        ELSEIF (LWCOU .AND. MARSTYPE .EQ. '4v') THEN
          CDATE=CDATEA
          CALL DIFDATE (CDATEA, CDTPRO, IFCST)
          IFCST = IFCST/3600
        ELSE
          CDATE=CDTPRO
          IFCST = 0
        ENDIF
      ELSE
!*    0.2.  THIS IS A  FORECAST DATE.
        CDATE=CDATEF
        CALL DIFDATE (CDATEF, CDTPRO, IFCST)
        IFCST = IFCST/3600
      ENDIF
!-----------------------------------------------------------------------

!*    3. OUTPUT GRIB DATA (USING A MAXIMUM OF NPROC/NWRTOUTWAM TASKS)
!     -------------------

      DO J = 1,NGY
        DO I = 1,NGX
          FIELD(I,J)=ZMISS
        ENDDO
      ENDDO

!     IPR: for a given PE, it gives the PE to which data should be sent
!     in order to do the writing to FDB or to file.
      IF (LFDBIOOUT) THEN
        DO IP=1,NPROC
          IPR(IP)=1+(IP-1)/NWRTOUTWAM
        ENDDO
      ELSE
        IPR(:)=1
      ENDIF

      IUOUT=0
      IF (.NOT.LFDBIOOUT .AND. IPR(IRANK).EQ.IRANK) THEN
!       output to file should only take place on one PE
        CALL GRSTNAME(CDTPRO,CDATEF,'SGS',ICPLEN,CPATH,OUTFILEN)
        LFILE=LEN_TRIM(OUTFILEN)
        CALL IGRIB_OPEN_FILE(IUOUT,OUTFILEN(1:LFILE),'w')
      ENDIF

      ALLOCATE(ZSENDBUF(NPROC*(NEND(IRANK)-NSTART(IRANK)+1)))
      ALLOCATE(ZRECVBUF(NEND(NPROC)))

      ISTEP=NPROC

      DO IC=1,NFRE*NANG,ISTEP

        IST=ISTEP
        IF(IC+ISTEP.GT.NFRE*NANG ) IST=NFRE*NANG-IC+1

!       COLLECT CONTRIBUTIONS ON ALL PE's (or a subset if less fields
!       ---------------------------------- than PE's)
        ICNT = 0
        DO IP=1,NPROC
          IF (IP <= IST) THEN
            IM=(((IC-1)+IP-1)/NANG)+1
            IK=(IC-1)+IP-(IM-1)*NANG

            DO J1=NSTART(IRANK),NEND(IRANK)
              ICNT = ICNT + 1
              ZSENDBUF(ICNT) = SPEC(J1,IK,IM)
            ENDDO
                          
            ISENDCOUNTS(IP) = NEND(IRANK) - NSTART(IRANK) + 1
          ELSE
            ISENDCOUNTS(IP) = 0
          ENDIF

          IF (IRANK <= IST) THEN
            IRECVCOUNTS(IP) = NEND(IP) - NSTART(IP) + 1
          ELSE
             IRECVCOUNTS(IP) = 0 
          ENDIF
        ENDDO

        CALL GSTATS(674,0)
        CALL MPL_ALLTOALLV(ZSENDBUF,ISENDCOUNTS,
     &                     ZRECVBUF,IRECVCOUNTS,
     &                      CDSTRING='OUTWSPEC:')
        CALL GSTATS(674,1)
        IF (ITEST.GT.1)
     &    WRITE(IU06,*) 'SUB OUTWSPEC: MPL_ALLTOALLV DONE' 

        CALL GSTATS(1496,0)
        NPROMA=NEND(NPROC)/NWAM_BLKS+1
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ,IX,IY)
        DO JKGLO=1,NEND(NPROC),NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,NEND(NPROC))
          DO IJ=KIJS,KIJL
            IX = IXLG(IJ,1)
            IY = NGY- KXLT(IJ,1) +1
            FIELD(IX,IY) = ZRECVBUF(IJ)
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1496,1)

!-----------------------------------------------------------------------

!       OUTPUT FROM ALL PE's WHICH CONTAIN A CONTRIBUTION
!       -------------------------------------------------

        IPARAM=251
        ITABLE=140

!       GRIB ENCODING
        MSGSIZE(:)=0
        IF(IRANK.LE.IST) THEN
          IM=(((IC-1)+IRANK-1)/NANG)+1
          IK=(IC-1)+IRANK-(IM-1)*NANG

          CALL WGRIBENCODE(IU06, ITEST, NGX, NGY, FIELD,
     &                     ITABLE, IPARAM, 0, IK , IM, 
     &                     CDATE, IFCST, MARSTYPE,
     &                     IGRIB_HANDLE)

          CALL IGRIB_GET_MESSAGE_SIZE(IGRIB_HANDLE,KBYTES)
          MSGSIZE(IRANK)=(KBYTES+NPRECI-1)/NPRECI
        ENDIF

!       SHARE MESSAGE SIZES
        ISENDBUF(:)=MSGSIZE(IRANK)
        ISENDCOUNTS(:)=1
        IRECVCOUNTS(:)=1
        CALL GSTATS(674,0)
        CALL MPL_ALLTOALLV(ISENDBUF,ISENDCOUNTS,
     &                     IRECVBUF,IRECVCOUNTS,
     &                     CDSTRING='OUTWSPEC IGRIB_HANDLE:')
        CALL GSTATS(674,1)
        MSGSIZE(:)=IRECVBUF(:)


        ITAG=IC
        IF (MSGSIZE(IRANK).NE.0) THEN 
          ALLOCATE(ISENDMSG(MSGSIZE(IRANK)))
          CALL IGRIB_GET_MESSAGE(IGRIB_HANDLE,ISENDMSG)

!         SEND GRIB MESSAGE TO RELEVANT OUTPUT PE (IF IT IS NOT ONE)
          IF (IPR(IRANK).NE.IRANK) THEN
            CALL MPL_SEND(ISENDMSG(1:MSGSIZE(IRANK)),
     &                    KDEST=IPR(IRANK),KTAG=ITAG,
     &                    KMP_TYPE=JP_NON_BLOCKING_STANDARD,
     &                    KREQUEST=ISENDREQ(1),
     &                    KERROR=IERR,CDSTRING='OUTWSPEC:')
            IF(IERR.LT.0) CALL MPL_ABORT('MPL_SEND ERROR IN OUTWSPEC')

          ELSE
!           OUTPUT GRIB DATA THAT ARE ALREADY THERE
            CALL WGRIBOUT(IU06, ITEST,
     &                    LFDBIOOUT, CFDB2DSP, NWFDBREF, LFDBOPEN,
     &                    IUOUT,
     &                    IGRIB_HANDLE,MSGSIZE(IRANK),ISENDMSG(1))

            ILOCAL=ILOCAL+1
            DEALLOCATE(ISENDMSG)

          ENDIF
          CALL IGRIB_RELEASE(IGRIB_HANDLE)
        ENDIF

!       ANYTHING TO RECEIVE ON PE IRANK FROM OTHER PEs ?
        NRCV=0
        MAXMSGSIZE=0
        DO IR=1,IRANK-1
          IF (IPR(IR).EQ.IRANK.AND.MSGSIZE(IR).NE.0) THEN
            NRCV=NRCV+1
            MAXMSGSIZE=MAX(MAXMSGSIZE,MSGSIZE(IR))
          ENDIF
        ENDDO
        DO IR=IRANK+1,NPROC
          IF (IPR(IR).EQ.IRANK.AND.MSGSIZE(IR).NE.0) THEN
            NRCV=NRCV+1
            MAXMSGSIZE=MAX(MAXMSGSIZE,MSGSIZE(IR))
          ENDIF
        ENDDO

        IF(MAXMSGSIZE.GT.0) THEN
          ALLOCATE(IRECVMSG(MAXMSGSIZE))
!         RECEIVE FROM RELEVANT PE's
          DO IR=1,NRCV
            CALL MPL_RECV(IRECVMSG,KFROM=KFROM,KTAG=ITAG,
     &                    KOUNT=KRCOUNT,KRECVTAG=KRTAG,KERROR=IERR,
     &                    CDSTRING='OUTWSPEC :')
            IF(IERR.LT.0) CALL MPL_ABORT('MPL_RECV ERROR IN OUTWSPEC ')
            IF(KRTAG.NE.ITAG) CALL MPL_ABORT
     &          ('MPL_RECV ERROR OUTWSPEC:MISMATCHED TAGS' )
            IF(KRCOUNT.NE.MSGSIZE(KFROM)) CALL MPL_ABORT
     &          ('MPL_RECV ERROR OUTWSPEC:MISMATCHED MSG LENGTH')

!           OUTPUT GRIB DATA (IF ANY)
            CALL IGRIB_NEW_FROM_MESSAGE(IGRBHNDL,
     &                                IRECVMSG(1:MSGSIZE(KFROM)))
            CALL WGRIBOUT(IU06, ITEST,
     &                    LFDBIOOUT, CFDB2DSP, NWFDBREF, LFDBOPEN,
     &                    IUOUT,
     &                    IGRBHNDL,MSGSIZE(KFROM),IRECVMSG(1))

            ILOCAL=ILOCAL+1

            CALL IGRIB_RELEASE(IGRBHNDL)
          ENDDO
        ENDIF

!       IS THE SEND FINISHED?
        IF (MSGSIZE(IRANK).NE.0 .AND. IPR(IRANK).NE.IRANK) THEN
          CALL MPL_WAIT(ISENDMSG(1:MSGSIZE(IRANK)),KREQUEST=ISENDREQ(1),
     &           CDSTRING='OUTWSPEC: WAIT FOR SEND')
            DEALLOCATE(ISENDMSG)
        ENDIF

        IF(ALLOCATED(IRECVMSG)) DEALLOCATE(IRECVMSG)

      ENDDO ! IC

      DEALLOCATE(ZSENDBUF)
      DEALLOCATE(ZRECVBUF)

      IF(LFDBIOOUT.AND.ILOCAL.GT.0) THEN
        IGLOBAL = NANG*NFRE
        IERR = ISETFIELDCOUNTFDBSUBS(NWFDBREF,IGLOBAL,ILOCAL)
          IF(IERR.NE.0)THEN
            WRITE(IU06,*) ' ------------------------'
            WRITE(IU06,*) ' ERROR setting fdb field count '
            WRITE(IU06,*) ' in routine OUTWSPEC '
            WRITE(IU06,*) ' FDB ERROR CODE IS ',IERR
            WRITE(IU06,*) ' IGLOBAL = ',IGLOBAL
            WRITE(IU06,*) ' ILOCAL = ',ILOCAL
            WRITE(IU06,*) ' ------------------------'
            CALL FLUSH(IU06)
            CALL ABORT1
          ENDIF
      ELSEIF (.NOT.LFDBIOOUT .AND. IPR(IRANK).EQ.IRANK) THEN
!       done writing to file
        CALL IGRIB_CLOSE_FILE(IUOUT)
      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('OUTWSPEC',1,ZHOOK_HANDLE)
#endif

      RETURN

      END SUBROUTINE OUTWSPEC
