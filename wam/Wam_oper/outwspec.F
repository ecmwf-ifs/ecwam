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
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,IRGG     ,AMOWEP   ,
     &            AMOSOP   ,AMOEAP   ,AMONOP   ,XDELLA   ,XDELLO
      USE YOWMESPAS, ONLY : LFDBIOOUT
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
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
      INTEGER :: IC, IST, ISTEP, IM, IK, J1, ICNT, IP
      INTEGER :: IJ, IX, IY, I, J, K
      INTEGER :: JKGLO, KIJS, KIJL, NPROMA
      INTEGER :: IPARAM, ITABLE
      INTEGER :: IGLOBAL, ILOCAL
      INTEGER :: IERR, LFILE, IUOUT
      INTEGER, DIMENSION(NPROC) :: ISENDCOUNTS,IRECVCOUNTS

      REAL :: ZHOOK_HANDLE
      REAL :: FIELD(NGX, NGY)
      REAL :: SPEC (NINF-1:NSUP, NANG, NFRE)
      REAL, ALLOCATABLE, DIMENSION(:) :: ZSENDBUF,ZRECVBUF

      CHARACTER(LEN=14) :: CDATE 
      CHARACTER(LEN=120) :: OUTFILEN

      IF (LHOOK) CALL DR_HOOK('OUTWSPEC',0,ZHOOK_HANDLE)

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

!*    3. OUTPUT GRIB DATA 
!     -------------------

      DO J = 1,NGY
        DO I = 1,NGX
          FIELD(I,J)=ZMISS
        ENDDO
      ENDDO
      
      ISTEP=NPROC
      IF (.NOT.LFDBIOOUT) ISTEP=1

      ALLOCATE(ZSENDBUF(NPROC*(NEND(IRANK)-NSTART(IRANK)+1)))
      ALLOCATE(ZRECVBUF(NEND(NPROC)))

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

        CALL GSTATS(1433,0)
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
        CALL GSTATS(1433,1)

!-----------------------------------------------------------------------

!       OUTPUT FROM ALL PE's WHICH CONTAIN A CONTRIBUTION
!       -------------------------------------------------

        IPARAM=251
        ITABLE=140

        IF(IRANK.LE.IST) THEN
          IM=(((IC-1)+IRANK-1)/NANG)+1
          IK=(IC-1)+IRANK-(IM-1)*NANG

          ILOCAL=ILOCAL+1

          IF (.NOT.LFDBIOOUT) THEN
!         output to file should only take place on one PE, thus istep=1
            IF(IM.EQ.1.AND.IK.EQ.1) THEN
              CALL GRSTNAME(CDTPRO,CDATEF,'SGS',ICPLEN,CPATH,OUTFILEN)
              LFILE=LEN_TRIM(OUTFILEN)
              CALL IGRIB_OPEN_FILE(IUOUT,OUTFILEN(1:LFILE),'w')
            ENDIF
          ENDIF

          CALL GSTATS(1704,0)

          CALL WGRIBOUT(IU06, ITEST, NGX, NGY, FIELD,
     &                  ITABLE, IPARAM, 0, IK , IM, 
     &                  CDATE, IFCST, MARSTYPE,
     &                  LFDBIOOUT, CFDB2DSP, NWFDBREF, LFDBOPEN,IUOUT)

          CALL GSTATS(1704,1)

          IF (.NOT.LFDBIOOUT) THEN
            IF(IM.EQ.NFRE.AND.IK.EQ.NANG) THEN
              CALL IGRIB_CLOSE_FILE(IUOUT)
            ENDIF
          ENDIF

        ENDIF

      ENDDO

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
      ENDIF

      IF (LHOOK) CALL DR_HOOK('OUTWSPEC',1,ZHOOK_HANDLE)
      RETURN

      END SUBROUTINE OUTWSPEC
