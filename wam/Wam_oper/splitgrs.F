C----------------------------------------------------------------------
C
      SUBROUTINE SPLITGRS(IU06, INFILEN, KSPLIT, CDATE, ITEST,
     &                    LFDBIOOUT,
     &                    CDTYPE, CDCLASS, CDEXPVER, CDFDB2DSP,
     &                    CDDOMAIN)
C
C----------------------------------------------------------------------
C
C**** *SPLITGRS*  PERFORMS THE SPLITTING OF GRIBBED WAVE MODEL SPECTRA 
C                 INTO NFRE*NANG FILES. 
C
C     J. BIDLOT   ECMWF   MARCH 1998 
C
C*    PURPOSE.
C     --------
C
C       PERFORMS THE SPLITTING OF GRIBBED WAVE SPECTRA INTO NFRE*NANG FILES
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *SPLITGRS*(IU06,INFILEN,KSPLIT,CDATE,ITEST,
C                         LFDBIOOUT,
C                         CDTYPE, CDCLASS, CDEXPVER, CDFDB2DSP,
C                         CDDOMAIN)
C
C*     VARIABLE.   TYPE.     PURPOSE.
C      ---------   -------   --------
C      *IU06*      INTEGER   OUTPUT UNIT FOR STANDARD OUTPUT.
C      *INFILEN*   CHARACTER FILENAME OF INPUT GRIBBED SPECTRA
C      *KSPLIT*    INTEGER   NUMBER OF SUBDOMAINS USED TO SAVE THE FULL
C                            SPECTRA IN GRIB (1 or 2)
C                  !!!       IF THERE ARE 2 DOMAINS, THE FIRST REQUEST
C                            MUST BE THE NORTHERN HEMISPHERE INCLUDING
C                            THE EQUATOR, AND THE SECOND DOMAIN THE
C                            SOUTHERN HEMISPHERE WIHTOUT THE EQUATOR.
C                            THE WHOLE REQUEST BEING SAVED IN FILENAME.
C      *CDATE*     CHARACTER DATE OF THE REQUESTED FIELD
C      *ITEST*     INTEGER   TEST OUTPUT LEVEL
C      *LFDBIOOUT* LOGICAL   SELECT TYPE OF OUTPUT OF SPLIT FILES 
C                            IF TRUE THE FDB SOFTWARE WILL BE USED FOR
C                            OUTPUT.
C      *CDTYPE*    CHARACTER TYPE OF FIELD ( FG, AN OR FC).
C      *CDCLASS*   CHARACTER CLASS OF CURRENT RUN.
C      *CDEXPVER*  CHARACTER EXPERIMENT VERSION NUMBER.
C      *CDDOMAIN*  CHARACTER FDB DOMAIN IDENTIFIER.
C
C     METHOD.
C     -------
C
C	   GET THE GRIB SPECTRA, DECODE IT AND SPLIT IT INTO NFRE*NANG FILES 
C      BY GRIBBING. ALSO ADJUST HEADERS TO BE CONSISTENT WITH WHAT IS
C      EXPECTED AS OUTPUT
C
C     EXTERNALS.
C     ----------
C
C      *ABORT1*
C      *GRIBEX* 
C      *GRSTNAME* 
C      *PBCLOSE*
C      *PBGRIB*
C      *PBOPEN*
C      *PBSIZE*
C      *PBWRITE*
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C-------------------------------------------------------------------
C
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.
C
#include "param.h"
C
#include "parcons.h"
C
#include "comfred.h"
C
#include "comgrid.h"
C
#include "commap.h"
C
#include "commpp.h"
C
C-----------------------------------------------------------------------
C
      INTEGER :: OUTGRIBSIZE
      INTEGER :: KSEC0(2),KSEC1(49+NANG+NFRE),KSEC2(384),
     &           KSEC3(2),KSEC4(128)
      INTEGER :: ISEC0(2),ISEC1(49+NANG+NFRE),ISEC2(384),
     &           ISEC3(2),ISEC4(128)
      INTEGER, ALLOCATABLE :: INGRIBN(:),INGRIBS(:),OUTGRIB(:)

      REAL, ALLOCATABLE :: PSEC4(:),PSEC4OUT(:)
      REAL :: PSEC2(96),PSEC3(2)
      REAL :: ZSEC2(96),ZSEC3(2)
C
      CHARACTER*1 MODE
      CHARACTER*12 CDATE 
      CHARACTER*40 CLMSG
      CHARACTER*80 INFILEN,OUTFILEN
      LOGICAL LLEXIST,LDFDB, LFDBIOOUT
C
      CHARACTER *(*) CDFDB2DSP, CDTYPE, CDCLASS, CDEXPVER, CDDOMAIN
C-----------------------------------------------------------------------
C
C*    0. INITIALISATION 
C     ------------------
C
      IF(KSPLIT.EQ.2) THEN
        SPLITLAT=0.0
      ELSE
        SPLITLAT=-90.0
      ENDIF

      KSEC1=0
      KSEC3(1)=0
      KSEC3(2)=0
      PSEC3(2)=0.
      PSEC3(2)=ZMISS

C
C*    1. INPUT OF GRIB DATA 
C     ----------------------
C
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
          ALLOCATE(INGRIBN(ISIZEN))
          CALL PBGRIB(IUGRS,INGRIBN,IPLENGN,ILENGN,IRET)

        ELSE
          CALL PBSIZE (IUGRS, IPLENGS)
          ISIZES=(IPLENGS+NPRECI-1)/NPRECI
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
C
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
C
C    2.0 SPLITTING 
C        ---------
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
C
      ALLOCATE(PSEC4OUT(NGRIBTOT))
C
      NFRANG=NANG*NFRE
      ILENGOUT=NGX*NGY/NPRECI
      ALLOCATE(OUTGRIB(ILENGOUT))
C
      ISEC0 = KSEC0
      ISEC1 = KSEC1
      ISEC2 = KSEC2
      ISEC3 = KSEC3
      ISEC4 = KSEC4
      ZSEC2 = PSEC2
      ZSEC3 = PSEC3
C
      ISEC1( 6) =  251 
C
      ISEC1(37) = 13
C 
      IF ( CDCLASS .EQ. 'OD' .OR. CDCLASS .EQ. 'od' ) THEN
        ISEC1(38) = 1
      ELSEIF (CDCLASS .EQ. 'ER' .OR. CDCLASS .EQ. 'er' ) THEN
        ISEC1(38) = 3
      ELSEIF (CDCLASS .EQ. 'E4' .OR. CDCLASS .EQ. 'e4' ) THEN
        ISEC1(38) = 5
      ELSE
        ISEC1(38) = 2
      ENDIF
C
      ISEC1(41) = 0
      DO JI = 1, 4
        ISEC1(41) = ISEC1(41)*256 + ICHAR(CDEXPVER(JI:JI))
      ENDDO
C
      ISEC1(44) = K 
      ISEC1(45) = M 
C
      ISEC1(46) = NANG
      ISEC1(47) = NFRE
C
      ISEC1(48) = 1000
      ISEC1(49) = 1000000
C
      DO IC=1,NANG
         ISEC1(49+IC) = INT(ISEC1(48)*TH(IC)*DEG)
      ENDDO 

      DO IC=1,NFRE
         ISEC1(49+NANG+IC) = INT(ISEC1(49)*FR(IC))
      ENDDO 
C
      IF( CDDOMAIN == 'm' ) THEN
        ISEC2( 3) =  NGY
        ISEC2( 4) =  NINT(AMONOP*1000)
        ISEC2( 7) =  NINT(AMOSOP*1000)
      ELSE
        ISEC2( 4) =  90000 
        ISEC2( 7) = -90000 
        ISEC2( 3) =  NINT(180./XDELLA) + 1 
      ENDIF
      ISEC2( 5) =  NINT(AMOWEP*1000)
      ISEC2( 8) =  NINT(AMOEAP*1000)
C
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
C
      ISEC4(1) = NGRIBTOT
C  
C   SET HEADER TO SINGLE DATUM AT EACH GRID POINT AND NO SECONDARY BITMAPS !!
      ISEC4( 8) = 0
      ISEC4( 9) = 0
      ISEC4(50) = 1
      ISEC4(51) = 1
C
      IF (LFDBIOOUT) THEN
        CONTINUE
      ELSE
        CALL GRSTNAME(IU06,CDATE,CDATE,1,'SGS',IRANK,NPROC,OUTFILEN)
        LFILE=LEN_TRIM(OUTFILEN)
        CALL PBOPEN(IUOUT,OUTFILEN(1:LFILE),'w',KRET)
        IF(KRET.LT.0) GOTO 3200
      ENDIF
C
      DO M=1,NFRE
        DO K=1,NANG
C
          ISEC1(44) = K
          ISEC1(45) = M
C
          IF(IRGG.EQ.1 .OR. CDDOMAIN == 'm' ) THEN
            INDEX=0
          ELSE
            INDEX = (NINT((90. - AMONOP ) / XDELLA))*NGX
            PSEC4OUT = ZMISS
          ENDIF
C
          DO IJ=1,NTOT
            PSEC4OUT(IJ+INDEX)=PSEC4(K+(M-1)*NANG+(IJ-1)*NFRANG)
          ENDDO
C
          IF (ITEST.GT.0)
     .      WRITE(IU06,*) ' CALL TO GRIBEX, ILENGOUT=', ILENGOUT
          CALL GRIBEX (ISEC0, ISEC1, ISEC2, ZSEC2, ISEC3, ZSEC3, ISEC4,
     &                 PSEC4OUT, NGRIBTOT,OUTGRIB, ILENGOUT ,KWORD,'C',
     &                 KRET)
          CLMSG='ENCODING PSEC4OUT'
          IF(KRET.GT.0) GOTO 3500
C
          IF (LFDBIOOUT) THEN
            CALL FLD2FDB (IU06, itest, isec1, outgrib, kword, 'write',
     .                    ierr,
     .                    cdtype, cdclass,cdexpver,cdfdb2dsp,cddomain,
     .                    k, m, kfdb)
            ISTAT = iflushfdb (kfdb)
            CALL FLD2FDB (IU06, itest, isec1, outgrib, kount, 'close',
     .                    kerr,
     .                    cdtype, cdclass,cdexpver,cdfdb2dsp,cddomain,
     .                    k, m, kfdb)
          ELSE
            KOUNT=KWORD*NPRECI
            CALL PBWRITE(IUOUT,OUTGRIB,KOUNT,KRET)
            IF(KRET.LT.0) GOTO 3300
          ENDIF

        ENDDO
      ENDDO
C
      IF (.NOT.LFDBIOOUT) THEN
        CALL PBCLOSE(IUOUT,KRET)
        IF(KRET.LT.0) GOTO 3400 
      ENDIF
C
      IF(ALLOCATED(OUTGRIB)) DEALLOCATE(OUTGRIB)
      IF(ALLOCATED(PSEC4OUT)) DEALLOCATE(PSEC4OUT)
C
C      4.0 DEALLOCATION
C          ------------
 
      IF(ALLOCATED(PSEC4)) DEALLOCATE(PSEC4)
C     
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
C
C*    ERROR MESSAGES
C
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

      END
