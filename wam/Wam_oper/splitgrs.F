!----------------------------------------------------------------------

      SUBROUTINE SPLITGRS(IU06, INFILEN, KSPLIT, CDATE, ITEST,
     &                    LFDBIOOUT,
     &                    CDCLASS, CDEXPVER, CDFDB2DSP,
     &                    CDDOMAIN)

!----------------------------------------------------------------------

!**** *SPLITGRS*  PERFORMS THE SPLITTING OF GRIB WAVE MODEL SPECTRA 
!                 PARAMETER 250 INTO GRIB PARAMETER 251.  

!     J. BIDLOT   ECMWF   MARCH 1998 

!*    PURPOSE.
!     --------

!       PERFORMS THE SPLITTING OF GRIB WAVE SPECTRA PARAMETER 250
!       (1 or 2 FIELDS) INTO GRIB PARAMETER 251 NFRE*NANG FIELDS.
!       IT ALSO PERFORMS THE OUTPUT OF THESE NEW FIELDS.

!**   INTERFACE.
!     ----------

!       *CALL* *SPLITGRS*(IU06,INFILEN,KSPLIT,CDATE,ITEST,
!                         LFDBIOOUT,
!                         CDCLASS, CDEXPVER, CDFDB2DSP,
!                         CDDOMAIN)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *IU06*      INTEGER   OUTPUT UNIT FOR STANDARD OUTPUT.
!      *INFILEN*   CHARACTER FILENAME OF INPUT GRIBBED SPECTRA
!      *KSPLIT*    INTEGER   NUMBER OF SUBDOMAINS USED TO SAVE THE FULL
!                            SPECTRA IN GRIB (1 or 2)
!                  !!!       IF THERE ARE 2 DOMAINS, THE FIRST REQUEST
!                            MUST BE THE NORTHERN HEMISPHERE INCLUDING
!                            THE EQUATOR, AND THE SECOND DOMAIN THE
!                            SOUTHERN HEMISPHERE WIHTOUT THE EQUATOR.
!                            THE WHOLE REQUEST BEING SAVED IN FILENAME.
!      *CDATE*     CHARACTER DATE OF THE REQUESTED FIELD
!      *ITEST*     INTEGER   TEST OUTPUT LEVEL
!      *LFDBIOOUT* LOGICAL   SELECT TYPE OF OUTPUT OF SPLIT FILES 
!                            IF TRUE THE FDB SOFTWARE WILL BE USED FOR
!                            OUTPUT.
!      *CDCLASS*   CHARACTER CLASS OF CURRENT RUN.
!      *CDEXPVER*  CHARACTER EXPERIMENT VERSION NUMBER.
!      *CDDOMAIN*  CHARACTER FDB DOMAIN IDENTIFIER.

!     METHOD.
!     -------

!      GET THE GRIB SPECTRA, DECODE IT AND SPLIT IT INTO NFRE*NANG FILES
!      BY GRIBBING. ALSO ADJUST HEADERS TO BE CONSISTENT WITH WHAT IS
!      EXPECTED AS OUTPUT. ALSO PERFORM THE OUTPUT OF THESE NEW FIELDS.

!     EXTERNALS.
!     ----------

!      *ABORT1*
!      *GRIBEX* 
!      *GRSTNAME* 
!      *PBCLOSE*
!      *PBGRIB*
!      *PBOPEN*
!      *PBSIZE*
!      *PBWRITE*

!     REFERENCE.
!     ----------

!       NONE.

!-------------------------------------------------------------------

      USE YOWFRED  , ONLY : FR       ,C        ,TH
      USE YOWGRIBHD, ONLY : NKSEK0   ,NKSEK1   ,NKSEK2   ,NKSEK3,
     &            NKSEK4   ,NPSEK2   ,NPSEK3   ,KSEK0    ,KSEK1,
     &            KSEK2    ,KSEK3    ,KSEK4    ,PSEK2    ,PSEK3
      USE YOWGRID  , ONLY : NLONRGG
      USE YOWMAP   , ONLY : IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,XDELLA
      USE YOWMPP   , ONLY : NPRECI
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY
      USE YOWPCONS , ONLY : DEG      ,R        ,ZMISS
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH
      USE FDBSUBS

!-----------------------------------------------------------------------

      INTEGER :: KSEC0(NKSEK0), KSEC1(NKSEK1), KSEC2(NKSEK2),
     &           KSEC3(NKSEK3), KSEC4(NKSEK4)
      INTEGER :: ISEC0(NKSEK0), ISEC1(NKSEK1), ISEC2(NKSEK2),
     &           ISEC3(NKSEK3), ISEC4(NKSEK4)

      INTEGER, ALLOCATABLE :: INGRIBN(:),INGRIBS(:),OUTGRIB(:)

      REAL, ALLOCATABLE :: PSEC4(:),PSEC4OUT(:)
      REAL :: PSEC2(NPSEK2), PSEC3(NPSEK3)
      REAL :: ZSEC2(NPSEK2), ZSEC3(NPSEK3)

      CHARACTER*1 MODE
      CHARACTER*12 CDATE 
      CHARACTER*40 CLMSG
      CHARACTER*120 INFILEN,OUTFILEN
      LOGICAL LLEXIST, LFDBIOOUT, LFDBOPEN

      CHARACTER *(*) CDFDB2DSP, CDCLASS, CDEXPVER, CDDOMAIN
!-----------------------------------------------------------------------

!*    0. INITIALISATION 
!     ------------------

      IF(KSPLIT.EQ.2) THEN
        SPLITLAT=0.0
      ELSE
        SPLITLAT=-90.0
      ENDIF

      KSEC0 = KSEK0
      KSEC1 = KSEK1
      KSEC3 = KSEK3
      PSEC3 = PSEK3

!*    1. INPUT OF GRIB DATA 
!     ----------------------

      LFILE=0
      LLEXIST=.FALSE.
      IF (INFILEN.NE. ' ') LFILE=LEN_TRIM(INFILEN)
      INQUIRE(FILE=INFILEN(1:LFILE),EXIST=LLEXIST)
      IF(LLEXIST) THEN
        MODE='r'
        CALL PBOPEN(IUGRS,INFILEN(1:LFILE),MODE,KRET)
        IF(KRET.LT.0) THEN
          WRITE (IU06,*) '****************************************'
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO PBOPEN     *'
          WRITE (IU06,*) '*   IN SPLITGRS                        *'
          IF(KRET.EQ.-1)
     &      WRITE (IU06,*) 'COULD NOT OPEN FILE ',INFILEN
          IF(KRET.EQ.-2)
     &      WRITE (IU06,*) 'INVALID FILENAME ',INFILEN
          IF(KRET.EQ.-3) WRITE (IU06,*) 'INVALID OPEN MODE SPECIFIED'
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '****************************************'
          CALL ABORT1
        ENDIF
      ELSE
        WRITE(IU06,*)'***********************************************'
        WRITE(IU06,*)'*                                             *'
        WRITE(IU06,*)'*SUB. SPLITGRS : GRIBBED SPECTRA NOT FOUND IN *'
        WRITE(IU06,*)  INFILEN
        WRITE(IU06,*)'*PROGRAM WILL ABORT                           *'
        WRITE(IU06,*)'*                                             *'
        WRITE(IU06,*)'***********************************************'
        CALL ABORT1
      ENDIF


      DO ISPLIT=1,KSPLIT
        IF(ISPLIT.EQ.1) THEN
          CALL PBSIZE (IUGRS, IPLENGN)
          ISIZEN=(IPLENGN+NPRECI-1)/NPRECI
          IPLENGN=ISIZEN*NPRECI
          ALLOCATE(INGRIBN(ISIZEN))
          CALL PBGRIB(IUGRS,INGRIBN,IPLENGN,ILENGN,IRET)

        ELSE
          CALL PBSIZE (IUGRS, IPLENGS)
          ISIZES=(IPLENGS+NPRECI-1)/NPRECI
          IPLENGS=ISIZES*NPRECI
          ALLOCATE(INGRIBS(ISIZES))
          CALL PBGRIB(IUGRS,INGRIBS,IPLENGS,ILENGS,IRET)
        ENDIF
        IF (IRET.EQ.-1) THEN
          WRITE (IU06,*) ' SUB. SPLITGRS : '
          WRITE (IU06,*) ' REACHED EOF ON UNIT ',IUGRS
          CALL ABORT1
        ELSEIF (iret .EQ. -2) THEN
          WRITE (IU06,*) ' SUB. SPLITGRS : '
          WRITE (IU06,*) ' ERROR IN FILE HANDLING ON UNIT ',IUGRS
          CALL ABORT1
        ELSEIF (iret .EQ. -3) THEN
          WRITE (IU06,*) ' SUB. SPLITGRS : '
          WRITE (IU06,*) ' INGRIB IS NOT BIG ENOUGH.'
          CALL ABORT1
        ENDIF

      ENDDO

      K=1
      DO WHILE (AMONOP-(K-1)*XDELLA.GT.SPLITLAT)
        K=K+1
      ENDDO
      IKSPLIT=K

      IF (IRGG.EQ.1) THEN
        KLENP=0
        DO K=1,NGY
          JSN=NGY-K+1
          KLENP = NLONRGG(JSN) + KLENP 
        ENDDO
        KLENP = KLENP * NANG * NFRE
      ELSE
        KLENP = NGY * NGX * NFRE * NANG
      ENDIF

      ALLOCATE(PSEC4(KLENP))
      PSEC4=ZMISS
      IF(ITEST.GT.3) CALL GRSDBG(1)
      CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &           PSEC4, KLENP, INGRIBN, ISIZEN , KWORD, 'D', KRET)
      CLMSG='DECODING INGRIBN'
      IF(KRET.GT.0) GOTO 3500 

      IF(KSPLIT.EQ.2) THEN
        IF (IRGG.EQ.1) THEN
          NOFF=0
          DO K=1,IKSPLIT
            JSN=NGY-K+1
            NOFF=NOFF+NLONRGG(JSN)
          ENDDO
          NOFF=NOFF* NANG * NFRE

          KLENPS=0
          DO K=IKSPLIT+1,NGY
            JSN=NGY-K+1
            KLENPS = NLONRGG(JSN) + KLENPS 
          ENDDO
          KLENPS = KLENPS * NANG * NFRE
        ELSE
          NOFF= IKSPLIT *  NGX * NFRE * NANG
          KLENPS = (NGY-IKSPLIT) * NGX * NFRE * NANG
        ENDIF
        CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &  PSEC4(NOFF+1), KLENPS, INGRIBS, ISIZES , KWORD, 'D', KRET)
        CLMSG='DECODING INGRIBS'
        IF(KRET.GT.0) GOTO 3500 
      ENDIF

      IF(ALLOCATED(INGRIBN)) DEALLOCATE(INGRIBN)
      IF(ALLOCATED(INGRIBS)) DEALLOCATE(INGRIBS)
      IF(KSEC1(38).EQ.1) THEN
        WRITE(IU06,*) "CDCLASS='OD'"
      ELSEIF(KSEC1(38).EQ.2) THEN
        WRITE(IU06,*) "CDCLASS='RD'"
      ELSEIF(KSEC1(38).EQ.3) THEN
        WRITE(IU06,*) "CDCLASS='ER'"
      ELSEIF(KSEC1(38).EQ.5) THEN
        WRITE(IU06,*) "CDCLASS='E4'"
      ELSE
        WRITE(IU06,*)' CLASS ', KSEC1(38) ,' NOT SUPPORTED '
        WRITE(IU06,*)' PROGRAM ABORTS '
        CALL ABORT1
      ENDIF
      IF (ITEST .GT. 0) CALL FLUSH (IU06)

!    2.0 SPLITTING 
!        ---------
      IF (IRGG.EQ.1) THEN
        NTOT=0
        DO K=1,NGY
          JSN=NGY-K+1
          NTOT = NLONRGG(JSN) + NTOT 
        ENDDO
        NTOT = NTOT
        NGRIBTOT = NTOT
      ELSE
        NTOT = NGY * NGX
        NGRIBTOT = NGX * (NINT(180./XDELLA) + 1 )
      ENDIF

      ALLOCATE(PSEC4OUT(NGRIBTOT))

      NFRANG=NANG*NFRE
      ILENGOUT=NGX*NGY/NPRECI
      ALLOCATE(OUTGRIB(ILENGOUT))

      ISEC0 = KSEC0
      ISEC1 = KSEC1
      ISEC2 = KSEC2
      ISEC3 = KSEC3
      ISEC4 = KSEC4
      ZSEC2 = PSEC2
      ZSEC3 = PSEC3

      ISEC1( 6) =  251 

      ISEC1(37) = 13

      IF ( CDCLASS .EQ. 'OD' .OR. CDCLASS .EQ. 'od' ) THEN
        ISEC1(38) = 1
      ELSEIF (CDCLASS .EQ. 'ER' .OR. CDCLASS .EQ. 'er' ) THEN
        ISEC1(38) = 3
      ELSEIF (CDCLASS .EQ. 'E4' .OR. CDCLASS .EQ. 'e4' ) THEN
        ISEC1(38) = 5
      ELSE
        ISEC1(38) = 2
      ENDIF

      ISEC1(41) = KSEK1(41) 

      ISEC1(44) = K 
      ISEC1(45) = M 

      ISEC1(46) = NANG
      ISEC1(47) = NFRE

      ISEC1(48) = 1000
      ISEC1(49) = 1000000

      DO IC=1,NANG
        ISEC1(49+IC) = INT(ISEC1(48)*TH(IC)*DEG)
      ENDDO 

      DO IC=1,NFRE
        ISEC1(49+NANG+IC) = INT(ISEC1(49)*FR(IC))
      ENDDO 

      ISEC2( 3) = KSEK2( 3) 
      ISEC2( 4) = KSEK2( 4) 
      ISEC2( 5) = KSEK2( 5) 
      ISEC2( 7) = KSEK2( 7) 
      ISEC2( 8) = KSEK2( 8) 

      IF (IRGG.EQ.1) THEN
        IF( CDDOMAIN == 'm' ) THEN
          DO IC=1,NGY
            JSN=NGY-IC+1
            ISEC2(22+IC) = NLONRGG(JSN)
          ENDDO
        ELSE
          ISTART = NINT((90. - AMONOP ) / XDELLA)
          DO IC=1,ISTART
            JSN=NGY-IC+1
            ISEC2(22+IC) = 0
          ENDDO
          DO IC=1,NGY
            JSN=NGY-IC+1
            ISEC2(22+IC+ISTART) = NLONRGG(JSN)
          ENDDO
        ENDIF
      ENDIF

      ISEC4(1) = NGRIBTOT

!   SET HEADER TO SINGLE DATUM AT EACH GRID POINT AND
!   NO SECONDARY BITMAPS !!
      ISEC4( 8) = 0
      ISEC4( 9) = 0
      ISEC4(50) = 1
      ISEC4(51) = 1

      IF (LFDBIOOUT) THEN
        CONTINUE
      ELSE
        CALL GRSTNAME(CDATE,CDATE,'SGS',ICPLEN,CPATH,OUTFILEN)
        LFILE=LEN_TRIM(OUTFILEN)
        CALL PBOPEN(IUOUT,OUTFILEN(1:LFILE),'w',KRET)
        IF(KRET.LT.0) GOTO 3200
      ENDIF

      DO M=1,NFRE
        DO K=1,NANG

          ISEC1(44) = K
          ISEC1(45) = M

          IF(IRGG.EQ.1 .OR. CDDOMAIN == 'm' ) THEN
            INDEX=0
          ELSE
            INDEX = (NINT((90. - AMONOP ) / XDELLA))*NGX
            PSEC4OUT = ZMISS
          ENDIF

          DO IJ=1,NTOT
            PSEC4OUT(IJ+INDEX)=PSEC4(K+(M-1)*NANG+(IJ-1)*NFRANG)
          ENDDO

          IF (ITEST.GT.0)
     &      WRITE(IU06,*) ' CALL TO GRIBEX, ILENGOUT=', ILENGOUT
          CALL GRIBEX (ISEC0, ISEC1, ISEC2, ZSEC2, ISEC3, ZSEC3, ISEC4,
     &                 PSEC4OUT, NGRIBTOT,OUTGRIB, ILENGOUT ,KWORD,'C',
     &                 KRET)
          CLMSG='ENCODING PSEC4OUT'
          IF(KRET.GT.0) GOTO 3500

          IF (LFDBIOOUT) THEN
            CALL FLD2FDB (IU06, itest, isec1, outgrib, kword, 'write',
     &                    ierr,
     &                    cdclass,cdexpver,cdfdb2dsp,cddomain,
     &                    k, m, kfdb,lfdbopen)
            ISTAT = iflushfdbsubs (kfdb)
            CALL FLD2FDB (IU06, itest, isec1, outgrib, kount, 'close',
     &                    kerr,
     &                    cdclass,cdexpver,cdfdb2dsp,cddomain,
     &                    k, m, kfdb,lfdbopen)
          ELSE
            KOUNT=KWORD*NPRECI
            CALL PBWRITE(IUOUT,OUTGRIB,KOUNT,KRET)
            IF(KRET.LT.0) GOTO 3300
          ENDIF

        ENDDO
      ENDDO

      IF (.NOT.LFDBIOOUT) THEN
        CALL PBCLOSE(IUOUT,KRET)
        IF(KRET.LT.0) GOTO 3400 
      ENDIF

      IF(ALLOCATED(OUTGRIB)) DEALLOCATE(OUTGRIB)
      IF(ALLOCATED(PSEC4OUT)) DEALLOCATE(PSEC4OUT)

!      4.0 DEALLOCATION
!          ------------
 
      IF(ALLOCATED(PSEC4)) DEALLOCATE(PSEC4)

      CALL PBCLOSE(IUGRS,KRET)
      IF(KRET.LT.0) THEN
        WRITE (IU06,*) '************************************'
        WRITE (IU06,*) '*                                  *'
        WRITE (IU06,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
        WRITE (IU06,*) '* IN SPLITGRS                      *'
        WRITE (IU06,*) '* FILE ',INFILEN
        WRITE (IU06,*) '*                                  *'
        WRITE (IU06,*) '************************************'
        CALL ABORT1
      ENDIF


      RETURN

!*    ERROR MESSAGES

3200  CONTINUE
      WRITE (IU06,*) '*************************************'
      WRITE (IU06,*) '*                                   *'
      WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO PBOPEN  *'
      WRITE (IU06,*) '*   IN SPLITGRS                     *'
      IF(KRET.EQ.-1)
     &        WRITE (IU06,*) ' COULD NOT OPEN FILE ',OUTFILEN
      IF(KRET.EQ.-2)
     &        WRITE (IU06,*) ' INVALID FILENAME ',OUTFILEN
      IF(KRET.EQ.-3) 
     &        WRITE (IU06,*)' INVALID OPEN MODE SPECIFIED'
      WRITE (IU06,*) '*                                   *'
      WRITE (IU06,*) '*************************************'
      CALL ABORT1
3300  CONTINUE
      WRITE (IU06,*) '**************************************'
      WRITE (IU06,*) '*                                    *'
      WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO PBWRITE   *'
      WRITE (IU06,*) '*  IN SPLITGRS                       *'
      WRITE (IU06,*) '*  FILE ',OUTFILEN
      WRITE (IU06,*) '*                                    *'
      WRITE (IU06,*) '**************************************'
      CALL ABORT1
3400  CONTINUE
      WRITE (IU06,*) '************************************'
      WRITE (IU06,*) '*                                  *'
      WRITE (IU06,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
      WRITE (IU06,*) '* IN SPLITGRS                      *'
      WRITE (IU06,*) '* FILE ',OUTFILEN
      WRITE (IU06,*) '*                                  *'
      WRITE (IU06,*) '************************************'
      CALL ABORT1
3500  CONTINUE
      WRITE(IU06,*) '*********************************************'
      WRITE(IU06,*) '*                                           *'
      WRITE(IU06,*) '*         FATAL ERROR IN SUB. SPLITGRS      *'
      WRITE(IU06,*) '*         ============================      *'
      WRITE(IU06,*) '* GRIBEX ERROR - KRET = ',KRET
      WRITE(IU06,*) '* WHILE ',CLMSG
      WRITE(IU06,*) '*                                           *'
      WRITE(IU06,*) '*    PROGRAM ABORTS     PROGRAM ABORTS      *'
      WRITE(IU06,*) '*                                           *'
      WRITE(IU06,*) '*********************************************'
      CALL ABORT1

      END SUBROUTINE SPLITGRS
