C
C ----------------------------------------------------------------------
C
      PROGRAM CONVERT_GRBSPEC
C
C     J BIDLOT ECMWF      MAY 1998
C
C     PURPOSE.
C     --------
C     CONVERTS GRIB WAVE SPECTRA FROM PARAMETER DEFINITION 250 TO 251
C     AND VICE VERSA.
C  
C     USAGE : convert_grbspec [-i input_filename] [-o output_filename]
C     ------                  [-t itest]
C    
C             default value for input_filename : input_spectra 
C                               output_filename : output_spectra 
C                               itest : 0 (no diagnostics)
C
C            itest = 1 : prints built in diagnostics
C            itest = 2 : prints built in diagnostics + gribex debugger 
C                        messages (!!! it can be big)
C
C
C     INPUT FILE REQUIREMENT :
C     -----------------------
C
C     The input file can only contain one parameter (250 or 251) and
C     one time and date. 
C
C     If the input file contains parameter 250, it must have been obtained
C     for all areas it was originally created. Namelly, for global model 
C     fields DOMAIN=G prior to 1996/12/5 12Z and DOMAIN=N and DOMAIN=S
C     afterwards.  
C
C     If the input file contains parameter 251, it must have been obtained
C     for all directions and frequencies in the default order. Namelly, the
C     mars request should be done with DIRECTION=1/TO/nang, and
C     FREQUENCY=1/TO/nfre. At this time nfre=25 and nang=12 for global model
C     data abd 24 for mediterranean data.
C
C
C     IMPORTANT NOTE :
C     --------------
C
C     !!!! converting parameter 251 to 250 and back to 251 will NOT
C     !!!! necessarily reproduce the same data as the initial file.
C     !!!! It is due to the fact that the coding of each field 251
C     !!!! for a given direction and frequency uses the maximum
C     !!!! of that particular field. However, when coding the data as
C     !!!! 250, the maximum of all field (for the coded area) is considered
C     !!!! and since only 9 bits are used to code, the truncation will be
C     !!!! different. 
C
C     RUNNING RUN :
C     ------------
C
C     If you are dealing with global 0.5 spectra, you will need at least
C     250 Mb with you are using a 32bit version and around 450 for a 64 bit.
C 
C     LIBRARY : EMOSLIB
C     -------
C
C ----------------------------------------------------------------------
C
      PARAMETER(NANGMAX=24,NFREMAX=25)
      INTEGER I4
      INTEGER GETCLO, GETCLA, IOPTVAL
      INTEGER :: KSEC0(2),KSEC1(49+NANGMAX+NFREMAX),KSEC2(384),
     &           KSEC3(2),KSEC4(128)
      INTEGER :: ISEC0(2),ISEC1(49+NANGMAX+NFREMAX),ISEC2(384),
     &           ISEC3(2),ISEC4(128)
      INTEGER, ALLOCATABLE :: INGRIB(:),OUTGRIB(:),KDOMRGG(:)
      REAL :: ONE,ZTHETA,ZFRE
      REAL :: PSEC2(96),PSEC3(2)
      REAL :: ZSEC2(96),ZSEC3(2)
      REAL, ALLOCATABLE :: PSEC4(:),PSEC4OUT(:)
      CHARACTER     CLOPTLET
      CHARACTER*  3 CLL1
      CHARACTER*  6 CLOPTS
      CHARACTER* 12 CLFMT
      CHARACTER*40 ERRMSG(-4:6), MSG
      CHARACTER*128 CLARG, FNAMEIN, FNAMEOUT
      LOGICAL LLEXIST
C
      DATA CLOPTS/'i;o;t;'/
      DATA ERRMSG/
     4            ' DECODED WITH BIT MAP                  *',
     3            '                                       *',
     2            '                                       *',
     1            '                                       *',
     1            ' NO ERROR                              *',
     1            ' END OF FILE ENCOUNTED                 *',
     2            ' DECODING ERROR SEE GRIBEX DESCRIPTION *',
     3            ' SUSPICIOUS TIME UNIT IN BLOCK 1       *',
     4            ' DIMENSION 1   IS TOO SMALL            *',
     5            ' NO 2D SPECTRA IN INPUT USE "INMARSB"  *',
     6            ' ENCODING ERROR SEE GRIBEX DESCRIPTION *'/
C
C ----------------------------------------------------------------------
C
C*    INITIAL VALUES SET AND CRACK COMMAND LINE.
C     -----------------------------------------
C
      ZMISS=-999.
      I4=1
      NPRECI = KIND(I4) 
      FNAMEIN='input_spectra'
      FNAMEOUT='output_spectra'
      ITEST=0
C
      CMDLINE: DO
        IOPTVAL=GETCLO(CLOPTS,CLARG)
        IF (IOPTVAL .LE. 0 )  THEN
          EXIT CMDLINE
        ENDIF
        CLOPTLET=CHAR(IOPTVAL)
C       GETS VARIABLE ARGUMENT FOR OPTION
        MORARG=GETCLA(CLARG)
        IF (MORARG.NE.0) THEN
          IF ( CLOPTLET .EQ. 'i' ) THEN
            FNAMEIN=CLARG
          ELSE IF ( CLOPTLET .EQ. 'o' ) THEN
            FNAMEOUT=CLARG
          ELSE IF ( CLOPTLET .EQ. 't' ) THEN
            I1=LEN_TRIM(CLARG)
            WRITE (CLL1,'(I3)') I1
            CLFMT = '(I'//CLL1//')'
            READ (CLARG(1:I1),FMT=CLFMT) ITEST
          ENDIF
        ENDIF
      ENDDO CMDLINE
C
      KSEC1=0
      KSEC2=0
      KSEC4=0
      KSEC3(1)=0
      KSEC3(2)=0
      PSEC2=0.
      PSEC3(2)=ZMISS
C
C*    INPUT FIRST GRIB DATA FILE
C     --------------------------
C
      LFILE=0
      LLEXIST=.FALSE.
      IF (FNAMEIN.NE. ' ') LFILE=LEN_TRIM(FNAMEIN)
      INQUIRE(FILE=FNAMEIN(1:LFILE),EXIST=LLEXIST)
      IF(LLEXIST) THEN
        CALL PBOPEN(IUGRS,FNAMEIN(1:LFILE),'r',KRET)
        IF(KRET.LT.0) THEN
          WRITE (*,*) '****************************************'
          WRITE (*,*) '*                                      *'
          WRITE (*,*) '*   ERROR FOLLOWING CALL TO PBOPEN     *'
          IF(KRET.EQ.-1)
     &      WRITE (*,*) 'COULD NOT OPEN FILE ',FNAMEIN
          IF(KRET.EQ.-2)
     &      WRITE (*,*) 'INVALID FILENAME ',FNAMEIN
          IF(KRET.EQ.-3) WRITE (*,*) 'INVALID OPEN MODE SPECIFIED'
          WRITE (*,*) '*                                      *'
          WRITE (*,*) '****************************************'
          CALL ABORT
        ENDIF
      ELSE
        WRITE(*,*)'****************************'
        WRITE(*,*)'*                          *'
        WRITE(*,*)'*GRIB SPECTRA NOT FOUND IN *'
        WRITE(*,*)  FNAMEIN 
        WRITE(*,*)'*PROGRAM WILL ABORT        *'
        WRITE(*,*)'*                          *'
        WRITE(*,*)'****************************'
        CALL ABORT
      ENDIF
C
C     GET SIZE OF INGRIB
C
      CALL PBSIZE (IUGRS, IPLENG)
      ISIZE=(IPLENG+NPRECI-1)/NPRECI
      ALLOCATE(INGRIB(ISIZE))
C
C     GET FIRST DATA FILE
C
      CALL PBGRIB(IUGRS,INGRIB,IPLENG,ILENG,KRET)
      IF     (KRET.EQ.-1) THEN
        WRITE (*,*) ' REACHED EOF IN ',FNAMEIN 
      ELSEIF (KRET .EQ. -2) THEN
        WRITE (*,*) ' ERROR IN FILE HANDLING IN ',FNAMEIN 
      ELSEIF (KRET .EQ. -3) THEN
        WRITE (*,*) ' INGRIB IS NOT BIG ENOUGH. CURRENT SIZE: ',IPLENG
        WRITE (*,*) ' REQUIRED SIZE SHOULD BE: ',ILENG
      ENDIF
C
C*    GET GRIB HEADERS 
C
      KRET  = 1
      ILENP = 1
      ALLOCATE(PSEC4(ILENP))
      IF (ITEST.gt.0)
     &     WRITE(*,*)' GETTING GRIB HEADER OF 1st INPUT FIELD' 
      IF (ITEST.GT.1) CALL GRSDBG (1)
      CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &             PSEC4, ILENP, INGRIB, ILENG , KWORD, 'J', KRET)
      IF (ITEST.GT.0) WRITE(*,*)' GRIBEX DONE status=' , KRET 
      IF(KRET.GT.0) THEN
        MSG = ERRMSG(2)
        CALL GRPRS0 (KSEC0)
        CALL GRPRS1 (KSEC0, KSEC1)
        WRITE (*,*) MSG
        CALL ABORT
      ENDIF
      IF(ALLOCATED(PSEC4)) DEALLOCATE(PSEC4)
C
C*    DETERMINE DATA FIELD CHARACTERISTICS 
C
      IPARAM = KSEC1(6)
      IRGG = KSEC2(17)
      AMONOP = FLOAT(KSEC2(4)/1000)+0.1*(MOD(KSEC2(4),1000)/100)+
     &         0.01*(MOD(KSEC2(4),100)/10)+0.001*MOD(KSEC2(4),10)
      AMOSOP = FLOAT(KSEC2(7)/1000)+0.1*(MOD(KSEC2(7),1000)/100)+
     &         0.01*(MOD(KSEC2(7),100)/10)+0.001*MOD(KSEC2(7),10)
      XDELLA = FLOAT(KSEC2(10))/1000
C
      NGY = KSEC2(3)
      NGYN = NGY 
      IF (ITEST.gt.0) WRITE(*,*) ' THE INPUT PARAMETER IS ',IPARAM 
      IF(IPARAM.EQ.250) THEN
        IF(IRGG.EQ.0) THEN
          NGYMAX = NGY
        ELSE
          IF(AMONOP.GT.0..AND.AMOSOP.LE.0.) THEN
            NGYMAX = 2*KSEC2(3)-1
          ELSE IF(AMONOP.GT.0..AND.AMOSOP.GT.0.) THEN
            NGYMAX = NGY 
          ELSE
            WRITE (*,*) 'THE GRIB DATA SHOULD BE INPUT SUCH THAT' 
            WRITE (*,*) 'THE NORTHERN HEMISPHERE WAS FIRST RETRIEVED'
            WRITE (*,*) 'AND THEN THE SOUTHERN PART.'
            WRITE (*,*) 'THE PROGRAM WILL ABORT'
            CALL ABORT
          ENDIF
        ENDIF
        NANG = KSEC4(53) 
        NFRE = KSEC4(55) 
      ELSE IF(IPARAM.EQ.251) THEN
        NGYMAX = NGY
        NANG = KSEC1(46)
        NFRE = KSEC1(47)
      ELSE
        WRITE(*,*) 'THE INPUT GRIB PARAMETER IS NOT 250 OR 251 BUT',
     &               IPARAM 
        WRITE(*,*) 'WHICH IS NOT A WAVE SPECTRUM PARAMETER !!!'
        WRITE(*,*) 'PROGRAM WILL ABORT'
        CALL ABORT
      ENDIF
C
      ALLOCATE(KDOMRGG(NGYMAX))
      NFRANG=NANG*NFRE
      IF(IRGG.EQ.0) THEN
        NGX = KSEC2(2)
        KDOMRGG = NGX
      ELSE
        ISTART=0
        DO WHILE(KSEC2(23+ISTART).EQ.0)
           ISTART=ISTART+1
        ENDDO
        NGX = 0
        DO J=1,NGY-ISTART
           KDOMRGG(J) = KSEC2(22+J+ISTART)
           NGX = MAX(NGX,KDOMRGG(J))
        ENDDO
      ENDIF
C
C     DECODE INPUT GRIB DATA
C     ----------------------
C
      IF(IPARAM.EQ.250) THEN
C
        IF(IRGG.EQ.0) THEN
          ILENP = NGY*NGX*NFRANG
          ILEN1 = NGY*NGX
        ELSE
          ILENP=0
          DO J=1,NGY-1
            ILENP = ILENP + KDOMRGG(J) 
          ENDDO
          ILENPS = ILENP
          IOFF = ILENP + NGX
          ILENP = 2*ILENP
          ILENP = ILENP + NGX
          ILEN1 = ILENP
          ILENP = ILENP*NFRANG
          ILENPS = ILENPS*NFRANG
          IOFF = IOFF*NFRANG
        ENDIF
C
        ALLOCATE(PSEC4(ILENP))
        PSEC4=ZMISS
        KSEC3(2)=0
        PSEC3(2)=ZMISS
C
        IF (ITEST.gt.0) WRITE(*,*)' DECODING 1st INPUT FIELD' 
        CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &               PSEC4, ILENP, INGRIB, ILENG , KWORD, 'D', KRET)
        IF (ITEST.gt.0) THEN 
          WRITE(*,*)' GRIBEX DONE status=' , KRET 
          WRITE(*,*) ' KSEC0 : '
          WRITE(*,*) KSEC0
          WRITE(*,*) ' KSEC1 : '
          WRITE(*,*) KSEC1
          WRITE(*,*) ' KSEC2 : '
          WRITE(*,*) KSEC2
          WRITE(*,*) ' KSEC3 : '
          WRITE(*,*) KSEC3
          WRITE(*,*) ' KSEC4 : '
          WRITE(*,*) KSEC4
          WRITE(*,*) ' PSEC2 : '
          WRITE(*,*) PSEC2
        ENDIF
        IF(KRET.GT.0) THEN
          MSG = ERRMSG(2)
          CALL GRPRS0 (KSEC0)
          CALL GRPRS1 (KSEC0, KSEC1)
          WRITE (*,*) MSG
          CALL ABORT
        ENDIF
C
        IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)
C
C       GET SECOND FIELD IF IT EXISTS
C
        IF(IRGG.EQ.1.AND.AMONOP.EQ.81..AND.AMOSOP.EQ.0.) THEN
          CALL PBSIZE (IUGRS, IPLENG)
          ISIZE=(IPLENG+NPRECI-1)/NPRECI
          ALLOCATE(INGRIB(ISIZE))
          CALL PBGRIB(IUGRS,INGRIB,IPLENG,ILENG,KRET)
          IF     (KRET.EQ.-1) THEN
            WRITE (*,*) ' REACHED END OF FILE IN ',FNAMEIN 
          ELSEIF (KRET .EQ. -2) THEN
            WRITE (*,*) ' ERROR IN FILE HANDLING IN ',FNAMEIN 
          ELSEIF (KRET .EQ. -3) THEN
            WRITE (*,*)'INGRIB IS NOT BIG ENOUGH. CURRENT SIZE: ',IPLENG
            WRITE (*,*) ' REQUIRED SIZE SHOULD BE: ',ILENG
          ENDIF
C
          KSEC3(2)=0
          PSEC3(2)=ZMISS
          IF (ITEST.gt.0) WRITE(*,*)' DECODING 2nd INPUT FIELD' 
          CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &         PSEC4(IOFF+1), ILENPS, INGRIB, ILENG , KWORD, 'D', KRET)
          IF (ITEST.gt.0) THEN 
            WRITE(*,*)' GRIBEX DONE status=' , KRET 
            WRITE(*,*) ' KSEC0 : '
            WRITE(*,*) KSEC0
            WRITE(*,*) ' KSEC1 : '
            WRITE(*,*) KSEC1
            WRITE(*,*) ' KSEC2 : '
            WRITE(*,*) KSEC2
            WRITE(*,*) ' KSEC3 : '
            WRITE(*,*) KSEC3
            WRITE(*,*) ' KSEC4 : '
            WRITE(*,*) KSEC4
            WRITE(*,*) ' PSEC2 : '
            WRITE(*,*) PSEC2
          ENDIF
          IF(KRET.GT.0) THEN
            MSG = ERRMSG(2)
            CALL GRPRS0 (KSEC0)
            CALL GRPRS1 (KSEC0, KSEC1)
            WRITE (*,*) MSG
            CALL ABORT
          ENDIF
C
          IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)
C
          NGY = KSEC2(3)
          DO J=1,NGY
             KDOMRGG(J+NGYN) = KSEC2(22+J)
          ENDDO
C
        ENDIF
C
      ELSE IF(IPARAM.EQ.251) THEN
C
        IF(IRGG.EQ.0) THEN
          ILENP = NGY*NGX*NFRANG
          ILEN1 = NGY*NGX
        ELSE
          ILENP=0
          DO J=1,NGY
            ILENP = ILENP + KDOMRGG(J) 
          ENDDO
          ILEN1 = ILENP
          ILENP = ILENP*NFRANG
        ENDIF
C
        ALLOCATE(PSEC4(ILEN1))
        ALLOCATE(PSEC4OUT(ILENP))
        PSEC4=ZMISS
        PSEC4OUT=ZMISS
C
        IF(IRGG.EQ.1) THEN
          INDEXN = 0 
          INDEXS = 0 
        ELSE
          IF(AMONOP.EQ.90..AND.AMOSOP.EQ.-90.) THEN
            INDEXN = (NINT((90. - 81. ) / XDELLA))*NGX
            INDEXS = (NINT((90. - 81. ) / XDELLA))*NGX
          ELSE
            INDEXN = 0 
            INDEXS = 0 
          ENDIF
        ENDIF
C
        DO M=1,NFRE
          DO K=1,NANG
C
            KSEC3(2)=0
            PSEC3(2)=ZMISS
            IF (ITEST.gt.0) 
     &          WRITE(*,*)' DECODING INPUT FIELD ',(M-1)*NANG+K 
            CALL GRIBEX(KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &               PSEC4, ILEN1, INGRIB, ILENG , KWORD, 'D', KRET)
            IF (ITEST.gt.0) THEN 
              WRITE(*,*)' GRIBEX DONE status=' , KRET 
              WRITE(*,*) ' KSEC0 : '
              WRITE(*,*) KSEC0
              WRITE(*,*) ' KSEC1 : '
              WRITE(*,*) KSEC1
              WRITE(*,*) ' KSEC2 : '
              WRITE(*,*) KSEC2
              WRITE(*,*) ' KSEC3 : '
              WRITE(*,*) KSEC3
              WRITE(*,*) ' KSEC4 : '
              WRITE(*,*) KSEC4
              WRITE(*,*) ' PSEC2 : '
              WRITE(*,*) PSEC2
            ENDIF
            IF(KRET.GT.0) THEN
              MSG = ERRMSG(2)
              CALL GRPRS0 (KSEC0)
              CALL GRPRS1 (KSEC0, KSEC1)
              WRITE (*,*) MSG
              CALL ABORT
            ENDIF
C
            IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)
C
            KK = KSEC1(44)
            MM = KSEC1(45)
            IF(KK.NE.K.AND.MM.NE.M) THEN
              WRITE (*,*) '************************************'
              WRITE (*,*) '* DIRECTION AND FREQUENCY INDEX ARE'
              WRITE (*,*) '* OUT OF ORDER. '
              WRITE (*,*) '* IN FILE ',FNAMEIN
              WRITE (*,*) '* WE EXPECTED DIRECTION : ', K 
              WRITE (*,*) '* AND GOT ', KK
              WRITE (*,*) '* WE EXPECTED FREQUENCY : ', M 
              WRITE (*,*) '* AND GOT ', MM 
              WRITE (*,*) '* VERIFY THAT YOUR MARS REQUEST WAS'
              WRITE (*,*) '* DONE FOR ALL DIRECTIONS AND FREQUENCIES'
              WRITE (*,*) '* IN THE DEFAULT ORDER.'
              WRITE (*,*) '************************************'
              CALL ABORT
            ENDIF
C
            DO IJ=1+INDEXN,ILEN1-INDEXS
             PSEC4OUT(K+(M-1)*NANG+(IJ-1)*NFRANG) =  PSEC4(IJ)
            ENDDO
C
C           GET NEXT FIELD
C
            IF(.NOT.(K.EQ.NANG.AND.M.EQ.NFRE)) THEN
              CALL PBSIZE (IUGRS, IPLENG)
              ISIZE=(IPLENG+NPRECI-1)/NPRECI
              ALLOCATE(INGRIB(ISIZE))
C
              CALL PBGRIB(IUGRS,INGRIB,IPLENG,ILENG,KRET)
              IF     (KRET.EQ.-1) THEN
                WRITE (*,*) ' REACHED EOF IN ',FNAMEIN 
                WRITE (*,*) ' CHECK WHETHER YOU HAVE PROVIDED THE INPUT' 
                WRITE (*,*) ' SPECTRA FOR ALL DIRECTIONS ie: ',NANG
                WRITE (*,*) ' AND FOR ALL FREQUENCIES ie: ',NFRE 
              ELSEIF (KRET .EQ. -2) THEN
                WRITE (*,*) ' ERROR IN FILE HANDLING IN ',FNAMEIN 
              ELSEIF (KRET .EQ. -3) THEN
               WRITE (*,*) ' INGRIB IS NOT BIG ENOUGH. CURRENT SIZE: ',
     &                       IPLENG
               WRITE (*,*) ' REQUIRED SIZE SHOULD BE: ',ILENG
              ENDIF
            ENDIF
C
          ENDDO
        ENDDO
C
        IF(ALLOCATED(PSEC4)) DEALLOCATE(PSEC4)
C
      ENDIF
C
      CALL PBCLOSE(IUGRS,KRET)
      IF(KRET.LT.0) THEN
        WRITE (*,*) '************************************'
        WRITE (*,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
        WRITE (*,*) '* FILE ',FNAMEIN
        WRITE (*,*) '************************************'
        CALL ABORT
      ENDIF
C
C     CONVERT INTO the OTHER FORMAT 
C     *****************************
C
      LFILE=0
      IF (FNAMEOUT.NE. ' ') LFILE=LEN_TRIM(FNAMEOUT)
      CALL PBOPEN(IUOUT,FNAMEOUT(1:LFILE),'w',KRET)
      IF(KRET.LT.0) THEN
        WRITE (*,*) '****************************************'
        WRITE (*,*) '*                                      *'
        WRITE (*,*) '*   ERROR FOLLOWING CALL TO PBOPEN     *'
        IF(KRET.EQ.-1)
     &    WRITE (*,*) 'COULD NOT OPEN FILE ',FNAMEOUT
        IF(KRET.EQ.-2)
     &    WRITE (*,*) 'INVALID FILENAME ',FNAMEOUT
        IF(KRET.EQ.-3) WRITE (*,*) 'INVALID OPEN MODE SPECIFIED'
        WRITE (*,*) '*                                      *'
        WRITE (*,*) '****************************************'
        CALL ABORT
      ENDIF
C
C     CONVERSION FROM 250 TO 251
C     --------------------------
C
      IF(IPARAM.EQ.250) THEN
C
C
        ILENGOUT=(1000+NGX*NGYMAX)/NPRECI
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
        ISEC1(46) = NANG
        ISEC1(47) = NFRE
C
        ISEC1(48) = 1000
        ISEC1(49) = 1000000
C
        DO IC=1,NANG
           ZTHETA = TRANSFER(KSEC4(59+IC),ONE) 
           ISEC1(49+IC) = NINT(ISEC1(48)*ZTHETA)
        ENDDO
C
        DO IC=1,NFRE
           ZFRE = TRANSFER(KSEC4(59+NANG+IC),ONE) 
           ISEC1(49+NANG+IC) = NINT(ISEC1(49)*ZFRE)
        ENDDO
C
        ISTART = 0
        INDEXN = 0 
        INDEXS = 0 
        IF((AMONOP.EQ.81..AND.AMOSOP.EQ.0.).OR.
     &     (AMONOP.EQ.81..AND.AMOSOP.EQ.-81.).OR.
     &     (AMONOP.EQ.72..AND.AMOSOP.EQ.-63.)) THEN
           ISEC2( 4) =  90000
           ISEC2( 7) = -90000
           ISEC2( 3) =  NINT(180./XDELLA) + 1
           ISTART = NINT((90. - AMONOP ) / XDELLA)
           IF(IRGG.EQ.1) THEN
             INDEXN=0
             INDEXS=0
           ELSE
             INDEXN = (NINT((90. - AMONOP ) / XDELLA))*NGX
             IF(AMOSOP.NE.0.) THEN
              INDEXS = (NINT((AMOSOP + 90. ) / XDELLA))*NGX
             ELSE
              INDEXS = INDEXN 
             ENDIF
           ENDIF
        ENDIF
C
        ILENO = ILEN1+INDEXN+INDEXS
        ALLOCATE(PSEC4OUT(ILENO))
C
        IF(IRGG.EQ.1) THEN
          DO IC=1,ISTART
            ISEC2(22+IC) = 0
          ENDDO
          DO IC=1,NGYMAX
            ISEC2(22+IC+ISTART) = KDOMRGG(IC)
          ENDDO
        ENDIF
C
        ISEC4(1) = ILENO 
C
        ISEC4( 8) = 0
        ISEC4( 9) = 0
        ISEC4(50) = 1
        ISEC4(51) = 1
C
        DO M=1,NFRE
          DO K=1,NANG
C
            ISEC1(44) = K
            ISEC1(45) = M
C
            ISEC3(2)=0
            ZSEC3(2)=ZMISS
C
            PSEC4OUT = ZMISS
            DO IJ=1,ILEN1
              PSEC4OUT(IJ+INDEXN)=PSEC4(K+(M-1)*NANG+(IJ-1)*NFRANG)
            ENDDO
C
            IF (ITEST.gt.0) THEN 
              WRITE(*,*)' CODING OUTPUT FIELD ',(M-1)*NANG+K 
              WRITE(*,*) ' ISEC0 : '
              WRITE(*,*) ISEC0
              WRITE(*,*) ' ISEC1 : '
              WRITE(*,*) ISEC1
              WRITE(*,*) ' ISEC2 : '
              WRITE(*,*) ISEC2
              WRITE(*,*) ' ISEC3 : '
              WRITE(*,*) ISEC3
              WRITE(*,*) ' ISEC4 : '
              WRITE(*,*) ISEC4
              WRITE(*,*) ' ZSEC2 : '
              WRITE(*,*) ZSEC2
            ENDIF
            CALL GRIBEX(ISEC0, ISEC1, ISEC2, ZSEC2, ISEC3, ZSEC3, ISEC4,
     &           PSEC4OUT, ILENO,OUTGRIB, ILENGOUT ,KWORD,'C',KRET)
            IF (ITEST.GT.0) WRITE(*,*)' GRIBEX DONE status=' , KRET 
            IF(KRET.GT.0) THEN
              MSG = ERRMSG(6)
              CALL GRPRS0 (KSEC0)
              CALL GRPRS1 (KSEC0, KSEC1)
              WRITE (*,*) MSG
              CALL ABORT
            ENDIF
C
            KOUNT=KWORD*NPRECI
            CALL PBWRITE(IUOUT,OUTGRIB,KOUNT,KRET)
            IF(KRET.LT.0) THEN
              WRITE (*,*) '**************************************'
              WRITE (*,*) '*  ERROR FOLLOWING CALL TO PBWRITE   *'
              WRITE (*,*) '**************************************'
              CALL ABORT
            ENDIF
          ENDDO
        ENDDO
C
        IF(ALLOCATED(OUTGRIB)) DEALLOCATE(OUTGRIB)
        IF(ALLOCATED(PSEC4OUT)) DEALLOCATE(PSEC4OUT)
C
C     CONVERSION FROM 251 TO 250
C     --------------------------
C
      ELSE IF(IPARAM.EQ.251) THEN
C
        ISEC0 = KSEC0
        ISEC1 = KSEC1
        ISEC2 = KSEC2
        ISEC3 = KSEC3
        ISEC4 = KSEC4
        ZSEC2 = PSEC2
        ZSEC3 = PSEC3
C
        ISEC1( 6) =  250
C
        ISEC1(37) = 1 
C
        ILENO = ILENP
        ILENGOUT=(1000+NGX*NGYMAX*NFRANG)/NPRECI
C
        DO IC=23,49+NANGMAX+NFREMAX
          ISEC2(IC) = 0 
        ENDDO
C
        IF(AMONOP.EQ.90..AND.AMOSOP.EQ.-90..AND.XDELLA.EQ.0.5) THEN
          ISEC2( 4) =  81000
          ISEC2( 7) =  0
          ISEC2( 3) = 163 
          ISHIFTN = ISEC2(3)
          ILENGOUT = (1000+NGX*(ISEC2(3)+1)*NFRANG)/NPRECI
        ELSE IF (AMONOP.EQ.90..AND.AMOSOP.EQ.-90.AND.XDELLA.EQ.1.5) THEN
          ISEC2( 4) =  81000
          ISEC2( 7) = -81000
          ISEC2( 3) = 109 
        ELSE IF (AMONOP.EQ.90..AND.AMOSOP.EQ.-90.AND.XDELLA.EQ.3.0) THEN
          ISEC2( 4) =  72000
          ISEC2( 7) = -63000
          ISEC2( 3) =  48 
        ENDIF
C
        IF(IRGG.EQ.0) THEN
          ILENO = ISEC2(3)*NGX*NFRANG
        ELSE
          ILENO = 0 
          DO IC=1,ISEC2(3)
            ISEC2(22+IC) = KDOMRGG(IC)
            ILENO = ILENO + KDOMRGG(IC)
          ENDDO
          ILENO = ILENO*NFRANG
        ENDIF
        ILENON = ILENO
C
        ISEC3(2)=0
        ZSEC3(2)=ZMISS
C
        ISEC4(1) = ILENO 
C
        ISEC4( 8) = 64 
        ISEC4( 9) = 32 
        ISEC4(50) = NANG 
        ISEC4(51) = NFRE 
C
        DO IC=1,NANG
           ZTHETA = FLOAT(KSEC1(49+IC))/KSEC1(48) 
           ISEC4(59+IC) = TRANSFER(ZTHETA,1) 
        ENDDO
C
        DO IC=1,NFRE
           ZFRE = FLOAT(KSEC1(49+NANG+IC))/KSEC1(49) 
           ISEC4(59+NANG+IC) = TRANSFER(ZFRE,1) 
        ENDDO
C
        ALLOCATE(OUTGRIB(ILENGOUT))
C
        IF (ITEST.gt.0) THEN 
          WRITE(*,*)' CODING 1st OUTPUT FIELD '
          WRITE(*,*) ' ISEC0 : '
          WRITE(*,*) ISEC0
          WRITE(*,*) ' ISEC1 : '
          WRITE(*,*) ISEC1
          WRITE(*,*) ' ISEC2 : '
          WRITE(*,*) ISEC2
          WRITE(*,*) ' ISEC3 : '
          WRITE(*,*) ISEC3
          WRITE(*,*) ' ISEC4 : '
          WRITE(*,*) ISEC4
          WRITE(*,*) ' ZSEC2 : '
          WRITE(*,*) ZSEC2
        ENDIF
        CALL GRIBEX(ISEC0, ISEC1, ISEC2, ZSEC2, ISEC3, ZSEC3, ISEC4,
     &         PSEC4OUT(1), ILENO, OUTGRIB, ILENGOUT ,KWORD,'C',KRET)
        IF (ITEST.gt.0) WRITE(*,*)' GRIBEX DONE status=' , KRET 
        IF(KRET.GT.0) THEN
          MSG = ERRMSG(6)
          CALL GRPRS0 (KSEC0)
          CALL GRPRS1 (KSEC0, KSEC1)
          WRITE (*,*) MSG
          CALL ABORT
        ENDIF
C
        KOUNT=KWORD*NPRECI
        CALL PBWRITE(IUOUT,OUTGRIB,KOUNT,KRET)
        IF(KRET.LT.0) THEN
        WRITE (*,*) '**************************************'
        WRITE (*,*) '*  ERROR FOLLOWING CALL TO PBWRITE   *'
        WRITE (*,*) '**************************************'
        CALL ABORT
        ENDIF
C
C       FOR SPLIT CONTRIBUTION
C
        IF(AMONOP.EQ.90..AND.AMOSOP.EQ.-90..AND.XDELLA.EQ.0.5) THEN
          ISEC2( 4) = -500 
          ISEC2( 7) = -81000 
          ISEC2( 3) = 162 
C 
          DO IC=23,49+NANGMAX+NFREMAX
            ISEC2(IC) = 0 
          ENDDO
          IF(IRGG.EQ.0) THEN
            ILENO = ISEC2(3)*NGX*NFRANG
          ELSE
            ILENO = 0 
            DO IC=1,ISEC2(3)
              ISEC2(22+IC) = KDOMRGG(ISHIFTN+IC)
              ILENO = ILENO + KDOMRGG(ISHIFTN+IC)
            ENDDO
            ILENO = ILENO*NFRANG
          ENDIF
C
          ISEC4(1) = ILENO 
C
          IF (ITEST.gt.0) THEN 
            WRITE(*,*)' CODING 2nd OUTPUT FIELD '
            WRITE(*,*) ' ISEC0 : '
            WRITE(*,*) ISEC0
            WRITE(*,*) ' ISEC1 : '
            WRITE(*,*) ISEC1
            WRITE(*,*) ' ISEC2 : '
            WRITE(*,*) ISEC2
            WRITE(*,*) ' ISEC3 : '
            WRITE(*,*) ISEC3
            WRITE(*,*) ' ISEC4 : '
            WRITE(*,*) ISEC4
            WRITE(*,*) ' ZSEC2 : '
            WRITE(*,*) ZSEC2
          ENDIF
          CALL GRIBEX(ISEC0, ISEC1, ISEC2, ZSEC2, ISEC3, ZSEC3, ISEC4,
     &    PSEC4OUT(ILENON+1), ILENO, OUTGRIB, ILENGOUT ,KWORD,'C',KRET)
          IF (ITEST.gt.0) WRITE(*,*)' GRIBEX DONE status=' , KRET 
          IF(KRET.GT.0) THEN
            MSG = ERRMSG(6)
            CALL GRPRS0 (KSEC0)
            CALL GRPRS1 (KSEC0, KSEC1)
            WRITE (*,*) MSG
            CALL ABORT
          ENDIF
C
          KOUNT=KWORD*NPRECI
          CALL PBWRITE(IUOUT,OUTGRIB,KOUNT,KRET)
          IF(KRET.LT.0) THEN
          WRITE (*,*) '**************************************'
          WRITE (*,*) '*  ERROR FOLLOWING CALL TO PBWRITE   *'
          WRITE (*,*) '**************************************'
          CALL ABORT
          ENDIF
        ENDIF
        IF(ALLOCATED(OUTGRIB)) DEALLOCATE(OUTGRIB)
        IF(ALLOCATED(PSEC4OUT)) DEALLOCATE(PSEC4OUT)
C
      ENDIF
C
      CALL PBCLOSE(IUOUT,KRET)
      IF(KRET.LT.0) THEN 
        WRITE (*,*) '************************************'
        WRITE (*,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
        WRITE (*,*) '************************************'
      ENDIF
C
      END
C
C#######################################################################
C      FUNCTION getclo(yaoptions, yaargument)
C      INTEGER getclo, getcla
C      CHARACTER*   1 yolastarg
C      CHARACTER* (*) yaoptions, yaargument
C      CHARACTER* 120 arg
CC
C      INTEGER here, imorearg, ivarg
C      DATA here, imorearg, ivarg, arg / 1, 0, 0, "  " /
C      DATA yolastarg / " " /
CC
C      arg=' '
C      CALL getarg(here,arg)
CC-->  PRINT*,'-------------getclo---------------'
CC-->  PRINT*,' ###',arg,'###'
CC
C      iol=len_trim(arg)
CC-->  PRINT*,' iol: ', iol
C      IF (iol .EQ. 2 .AND. arg(1:1) .EQ. '-' .AND. ivarg .EQ. 0 ) THEN
C        iol = len_trim(yaoptions)
CC-->    PRINT*,' iol 2',iol,' options: ', yaoptions
C        DO jl=1,iol
C          getclo = 0
C          IF ( yaoptions(jl:jl) .EQ. arg(2:2) ) THEN
C            getclo = ichar(arg(2:2))
CC-->        PRINT *,' FOUND ',yaoptions(jl:jl), ' IN THE COMMAND LINE',
CC--> .               yaoptions(jl+1:jl+1)
C            IF (yaoptions(jl+1:jl+1) .EQ. ':' ) THEN
CC-->          PRINT*, yaoptions(jl:jl),' requires arguments'
C              yolastarg=yaoptions(jl:jl)
C              ivarg=1
C            ENDIF
C            EXIT
C          ENDIF
C        ENDDO
C      ELSEIF ( ivarg .EQ. 1 ) THEN
C         WRITE(*,*) ' option -', yolastarg, ' requires arguments'
C         getclo=-1
C      ELSEIF (iol .EQ. 0) THEN
C        getclo=0
C      ELSE
C         WRITE(*,*) 'illegal option: ',arg(1:iol)
C         getclo=-1
C      ENDIF
CC-->  PRINT *,' HERE ins getclo', here
C      here = here + 1
C      RETURN
C
C      ENTRY getcla(yaargument)
CC-->  PRINT*,'-------------getcla--------------'
CC-->  PRINT*, 'HERE ins getcla :', here,' options: ', yaoptions
C
C      getcla = 1
C      CALL getarg(here,arg)
CC-->  PRINT*,' arg in getcla ', arg
C      IF ( arg (1:1) .NE. '-' ) THEN
C        here = here + 1
C        yaargument=arg
C      ELSE
C        IF (ivarg.EQ.1) THEN
C          WRITE(*,*)' refused to take ', arg (1:2) ,' as argument for',
C     .    ' the option -',yolastarg
C          getcla = -1
C        ELSE
C          getcla = 0
C        ENDIF
C      ENDIF
C      ivarg=0
CC-->  PRINT*,' getcla in getcla ', getcla
CC
C      RETURN
C      END
