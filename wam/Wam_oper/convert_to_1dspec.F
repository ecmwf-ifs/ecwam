
! ----------------------------------------------------------------------

      PROGRAM CONVERT_TO_1DSPEC

!     J BIDLOT ECMWF   SEPTEMBER 1998

!     PURPOSE.
!     --------
!     CONVERTS GRIB WAVE SPECTRA FROM PARAMETER DEFINITION 251
!     TO 1-D SPECTRA (252 : this is my own definition for now).
!     NOTE: IT CAN ALSO DO 250 IF IT'S PROVIDED AS ARCHIVED.
!     FOR EACH INDIVIDUAL FREQUENCY BIN OR INTEGRATED BETWEEN TWO
!     FREQUENCIES. IF THE  1-D SPECTRA IS INTEGRATED BETWEEN TWO
!     FREQUENCIES, THEN A WAVE HEIGHT TYPE OF INFORMATION CAN BE COMPUTED
!     BY TAKING 4*SQRT OF THE INTEGRATED QUATITY.

!     The grib headers are not entirely consistent with this new definition

!     !!! if you want to use this converter to plot 1-D spectra field
!         remember to provide data on a regular grid as metview will
!         not plot irregular grid data.
  
!     USAGE : convert_to_1dspec [-i input_filename] [-o output_filename]
!     ------                    [-f first_frequency] [-l last_frequency]
!                               [-s or -h] [-t itest]
!    
!             default value for input_filename : input_spectra 
!                               output_filename : output_spectra 
!                               first_frequency : 1
!                               last_frequency : last frequency available
!                               -s : not specified. No integration from
!                                    first to last frequency. (last_frequency-first_frequency+1)
!                                    fields will be produced.
!                               -h : not specified. No integration from
!                                    first to last frequency. (last_frequency-first_frequency+1)
!                                    fields will be produced and no conversion to a wave height
!                                    type of infomation.
!                               itest : 0 (no diagnostics)
!                               itest : 0 (no diagnostics)
!                               itest : 1 (some diagnostics)
!                               itest : 2 (some more diagnostics)
!                               itest : 3 (all diagnostics)

!            itest = 1 : prints built in diagnostics (only indication to calls)
!            itest = 2 : prints built in diagnostics
!            itest = 3 : prints built in diagnostics + gribex debugger 
!                        messages (!!! it can be big)


!     INPUT FILE REQUIREMENT :
!     -----------------------

!     The input file can only contain one parameter (250 or 251) and
!     one time and date.

!     If the input file contains parameter 250, it must have been obtained
!     for all areas it was originally created. Namelly, for global model
!     fields DOMAIN=G prior to 1996/12/5 12Z and DOMAIN=N and DOMAIN=S
!     afterwards.

!     If the input file contains parameter 251, it must have been obtained
!     for all directions and frequencies in the default order. Namelly, the
!     mars request should be done with DIRECTION=1/TO/nang, and
!     FREQUENCY=1/TO/nfre. At this time nfre=25,30 or 36 and nang=12 or 24
!     for global model data and 24 for mediterranean data.

!     RUNNING RUN :
!     ------------

!     LIBRARY : EMOSLIB
!     -------

! ----------------------------------------------------------------------

      PARAMETER(NANGMAX=24,NFREMAX=36)
      PARAMETER (G = 9.806, PI = 3.1415927, CIRC = 40000000.,
     1           ZPI = 2.*PI, RAD = PI/180., DEG = 180./PI,
     2           R = CIRC/ZPI, EPSMIN=0.1E-32, ZMISS=-999.0)
      PARAMETER (CO = 1.1)

      INTEGER I4
      INTEGER GETCLO, GETCLA, IOPTVAL
      INTEGER :: KSEC0(2),KSEC1(54+NANGMAX+NFREMAX),KSEC2(384),
     &           KSEC3(2),KSEC4(128)
      INTEGER :: ISEC0(2),ISEC1(54+NANGMAX+NFREMAX),ISEC2(384),
     &           ISEC3(2),ISEC4(128)
      INTEGER, ALLOCATABLE :: INGRIB(:),OUTGRIB(:),KDOMRGG(:)
      REAL :: ONE,ZTHETA,ZFRE
      REAL :: PSEC2(96),PSEC3(2)
      REAL :: ZSEC2(96),ZSEC3(2)
      REAL, ALLOCATABLE :: PSEC4(:),PSEC4OUT(:),PSEC41D(:,:)
      REAL, ALLOCATABLE :: FR(:), DFIM(:)
      CHARACTER     CLOPTLET
      CHARACTER*  3 CLL1
      CHARACTER*  14 CLOPTS
      CHARACTER* 12 CLFMT
      CHARACTER*40 ERRMSG(-4:6), MSG
      CHARACTER*128 CLARG, FNAMEIN, FNAMEOUT
      LOGICAL LLEXIST, LINTEGRATE, LWAVEHGT, LLASTFRE

      DATA CLOPTS/'i;o;f:l;s;h;t;'/
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

! ----------------------------------------------------------------------

!*    INITIAL VALUES SET AND CRACK COMMAND LINE.
!     -----------------------------------------

      I4=1
      NPRECI = KIND(I4) 
      FNAMEIN='input_spectra'
      FNAMEOUT='output_spectra'
      IFRE_FIRST=1
      IFRE_LAST=NFREMAX
      LINTEGRATE=.FALSE.
      LWAVEHGT=.FALSE.
      LLASTFRE=.FALSE.
      ITEST=0

      CMDLINE: DO
        IOPTVAL=GETCLO(CLOPTS,CLARG)
        IF (IOPTVAL .LE. 0 )  THEN
          EXIT CMDLINE
        ENDIF
        CLOPTLET=CHAR(IOPTVAL)
!       GETS VARIABLE ARGUMENT FOR OPTION
        MORARG=GETCLA(CLARG)
        IF (MORARG.NE.0) THEN
          IF ( CLOPTLET .EQ. 'i' ) THEN
            FNAMEIN=CLARG
          ELSE IF ( CLOPTLET .EQ. 'o' ) THEN
            FNAMEOUT=CLARG
          ELSE IF ( CLOPTLET .EQ. 'f' ) THEN
            I1=LEN_TRIM(CLARG)
            WRITE (CLL1,'(I3)') I1
            CLFMT = '(I'//CLL1//')'
            READ (CLARG(1:I1),FMT=CLFMT) IFRE_FIRST
          ELSE IF ( CLOPTLET .EQ. 'l' ) THEN
            I1=LEN_TRIM(CLARG)
            WRITE (CLL1,'(I3)') I1
            CLFMT = '(I'//CLL1//')'
            READ (CLARG(1:I1),FMT=CLFMT) IFRE_LAST
            LLASTFRE=.TRUE.
          ELSE IF ( CLOPTLET .EQ. 's' ) THEN
            LINTEGRATE=.TRUE.
          ELSE IF ( CLOPTLET .EQ. 'h' ) THEN
            LWAVEHGT=.TRUE.
          ELSE IF ( CLOPTLET .EQ. 't' ) THEN
            I1=LEN_TRIM(CLARG)
            WRITE (CLL1,'(I3)') I1
            CLFMT = '(I'//CLL1//')'
            READ (CLARG(1:I1),FMT=CLFMT) ITEST
          ENDIF
        ELSE IF ( CLOPTLET .EQ. 's' ) THEN
          LINTEGRATE=.TRUE.
        ELSE IF ( CLOPTLET .EQ. 'h' ) THEN
          LWAVEHGT=.TRUE.
        ENDIF
      ENDDO CMDLINE

      KSEC1=0
      KSEC2=0
      KSEC4=0
      KSEC3(1)=0
      KSEC3(2)=0
      PSEC2=0.
      PSEC3(2)=ZMISS

!*    INPUT FIRST GRIB DATA FILE
!     --------------------------

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

!     GET SIZE OF INGRIB

      CALL PBSIZE (IUGRS, IPLENG)
      ISIZE=(IPLENG+NPRECI-1)/NPRECI
      ISIZEMAX=ISIZE
      ALLOCATE(INGRIB(ISIZE))

!     GET FIRST DATA FILE

      CALL PBGRIB(IUGRS,INGRIB,IPLENG,ILENG,KRET)
      IF     (KRET.EQ.-1) THEN
        WRITE (*,*) ' REACHED EOF IN ',FNAMEIN 
      ELSEIF (KRET .EQ. -2) THEN
        WRITE (*,*) ' ERROR IN FILE HANDLING IN ',FNAMEIN 
      ELSEIF (KRET .EQ. -3) THEN
        WRITE (*,*) ' INGRIB IS NOT BIG ENOUGH. CURRENT SIZE: ',IPLENG
        WRITE (*,*) ' REQUIRED SIZE SHOULD BE: ',ILENG
      ENDIF

!*    GET GRIB HEADERS 

      KRET  = 1
      ILENP = 1
      ALLOCATE(PSEC4(ILENP))
      IF (ITEST.gt.0)
     &     WRITE(*,*)' GETTING GRIB HEADER OF 1st INPUT FIELD' 
      IF (ITEST.GT.2) CALL GRSDBG (1)
      CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &             PSEC4, ILENP, INGRIB, ILENG , KWORD, 'J', KRET)
      IF (ITEST.GT.1) WRITE(*,*)' GRIBEX DONE status=' , KRET 
!     the error code 811: 'Cannot handle 2ndary bitmaps for J option'
!     is actually a warning message and should not cause a failure.
      IF(KRET.GT.0.and.kret.ne.811) THEN
        MSG = ERRMSG(2)
        CALL GRPRS0 (KSEC0)
        CALL GRPRS1 (KSEC0, KSEC1)
        WRITE (*,*) MSG
        CALL ABORT
      ENDIF
      IF(ALLOCATED(PSEC4)) DEALLOCATE(PSEC4)

!*    DETERMINE DATA FIELD CHARACTERISTICS 

      IPARAM = KSEC1(6)
      IRGG = KSEC2(17)
      AMONOP = FLOAT(KSEC2(4)/1000)+0.1*(MOD(KSEC2(4),1000)/100)+
     &         0.01*(MOD(KSEC2(4),100)/10)+0.001*MOD(KSEC2(4),10)
      AMOSOP = FLOAT(KSEC2(7)/1000)+0.1*(MOD(KSEC2(7),1000)/100)+
     &         0.01*(MOD(KSEC2(7),100)/10)+0.001*MOD(KSEC2(7),10)
      XDELLA = FLOAT(KSEC2(10))/1000

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

      IF(.NOT.LLASTFRE) IFRE_LAST=NFRE

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

      ALLOCATE(FR(NFRE))

      DELTH=ZPI/NANG

!     DECODE INPUT GRIB DATA
!     ----------------------

      IF(IPARAM.EQ.250) THEN

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

        ALLOCATE(PSEC4(ILENP))
        PSEC4=ZMISS
        KSEC3(2)=0
        PSEC3(2)=ZMISS

        IF (ITEST.gt.0) WRITE(*,*)' DECODING 1st INPUT FIELD' 
        CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &               PSEC4, ILENP, INGRIB, ILENG , KWORD, 'D', KRET)
        IF (ITEST.gt.1) THEN 
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

        DO IC=1,NFRE
          FR(IC) = TRANSFER(KSEC4(59+NANG+IC),ONE) 
        ENDDO

        IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)

!       GET SECOND FIELD IF IT EXISTS

        IF(IRGG.EQ.1.AND.AMONOP.EQ.81..AND.AMOSOP.EQ.0.) THEN
          CALL PBSIZE (IUGRS, IPLENG)
          ISIZE=(IPLENG+NPRECI-1)/NPRECI
          ISIZEMAX=MAX(ISIZEMAX,ISIZE)
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

          KSEC3(2)=0
          PSEC3(2)=ZMISS
          IF (ITEST.gt.0) WRITE(*,*)' DECODING 2nd INPUT FIELD' 
          CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &         PSEC4(IOFF+1), ILENPS, INGRIB, ILENG , KWORD, 'D', KRET)
          IF (ITEST.gt.1) THEN 
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

          IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)

          NGY = KSEC2(3)
          DO J=1,NGY
             KDOMRGG(J+NGYN) = KSEC2(22+J)
          ENDDO

        ENDIF

        PPREC=0.0
        PPEPS=1.0e-10
        DO IJ=1,ILENP
          IF (PSEC4(IJ) .NE. ZMISS) THEN
            PSEC4(IJ) = 10.**(PSEC4(IJ)-ABS(PPREC))- PPEPS
          ENDIF
        ENDDO

      ELSE IF(IPARAM.EQ.251) THEN

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

        ALLOCATE(PSEC4OUT(ILEN1))
        PSEC4OUT=ZMISS

        ALLOCATE(PSEC41D(ILEN1,NFRE))
        PSEC41D=0

        DO M=1,NFRE
          DO K=1,NANG

            KSEC3(2)=0
            PSEC3(2)=ZMISS
            IF (ITEST.gt.0) 
     &          WRITE(*,*)' DECODING INPUT FIELD ',(M-1)*NANG+K 
            CALL GRIBEX(KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &               PSEC4OUT, ILEN1, INGRIB, ILENG , KWORD, 'D', KRET)
            IF (ITEST.gt.1) THEN 
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

            IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)

            KK = KSEC1(44)
            MM = KSEC1(45)
            FR(MM) =  FLOAT(KSEC1(49+NANG+MM))/KSEC1(49)
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

            PPREC=0.0
            PPEPS=1.0e-10
            DO IJ=1,ILEN1
              IF (PSEC4OUT(IJ) .NE. ZMISS) THEN
                PSEC4OUT(IJ) = 10.**(PSEC4OUT(IJ)-ABS(PPREC))- PPEPS
              ELSE
                PSEC4OUT(IJ) = 0.0
              ENDIF
              PSEC41D(IJ,M) = PSEC41D(IJ,M) + PSEC4OUT(IJ)
            ENDDO

!           GET NEXT FIELD

            IF(.NOT.(K.EQ.NANG.AND.M.EQ.NFRE)) THEN
              CALL PBSIZE (IUGRS, IPLENG)
              ISIZE=(IPLENG+NPRECI-1)/NPRECI
              ISIZEMAX=MAX(ISIZEMAX,ISIZE)
              ALLOCATE(INGRIB(ISIZE))

              IF (ITEST.gt.1) WRITE (*,*) ' CALLING PBGRIB ' 
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

          ENDDO
        ENDDO

        IF(ALLOCATED(PSEC4OUT)) DEALLOCATE(PSEC4OUT)

      ENDIF

      CALL PBCLOSE(IUGRS,KRET)
      IF(KRET.LT.0) THEN
        WRITE (*,*) '************************************'
        WRITE (*,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
        WRITE (*,*) '* FILE ',FNAMEIN
        WRITE (*,*) '************************************'
        CALL ABORT
      ENDIF

!     DETERMINE THE 1-D SPECTRA 
!     *************************

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

!     COMPUTE FREQUENCY SPACING IF NEEDED
      IF(LINTEGRATE.OR.LWAVEHGT) THEN
        ALLOCATE(DFIM(NFRE))
        CO1 = 0.5*(CO-1.)
        DFIM(1)= CO1*FR(1)
        DO M=2,NFRE-1
          DFIM(M)=CO1 * (FR(M)+FR(M-1))
        ENDDO
        DFIM(NFRE)=CO1*FR(NFRE-1)
      ENDIF

!     CONVERSION FROM 250 to 1-D spectra
!     -------------------

      IF(IPARAM.EQ.250) THEN

        ILENGOUT=2*(1000+NGX*NGYMAX)/NPRECI
        ALLOCATE(OUTGRIB(ILENGOUT))

        ISEC0 = KSEC0
        ISEC1 = KSEC1
        ISEC2 = KSEC2
        ISEC3 = KSEC3
        ISEC4 = KSEC4
        ZSEC2 = PSEC2
        ZSEC3 = PSEC3

        ISEC1( 6) = 252

        ISEC1(37) = 13

        ISEC1(46) = 0 
        ISEC1(47) = NFRE 

        ISEC1(48) = 1000
        ISEC1(49) = 1000000

        DO IC=1,NFRE
          ZFRE = TRANSFER(KSEC4(59+NANG+IC),ONE) 
          ISEC1(49+0+IC) = NINT(ISEC1(49)*ZFRE)
        ENDDO

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

        ILENO = ILEN1+INDEXN+INDEXS

        ALLOCATE(PSEC4OUT(ILENO))

        IF(IRGG.EQ.1) THEN
          DO IC=1,ISTART
            ISEC2(22+IC) = 0
          ENDDO
          DO IC=1,NGYMAX
            ISEC2(22+IC+ISTART) = KDOMRGG(IC)
          ENDDO
        ENDIF

        ISEC4(1) = ILENO 

        ISEC4( 8) = 0
        ISEC4( 9) = 0
        ISEC4(50) = 1
        ISEC4(51) = 1

        IF(LINTEGRATE.OR.LWAVEHGT) THEN

        IF (ITEST.gt.0)
     &      WRITE (*,*) ' COMPUTING integrated 1-D spectra from 250'

          ISEC1(44) = 0
          ISEC1(45) = IFRE_FIRST 
          ISEC3(2)=0
          ZSEC3(2)=ZMISS

          PSEC4OUT=0.

          DO M=IFRE_FIRST,IFRE_LAST
            DO K=1,NANG
              DO IJ=1,ILEN1
                IC=K+(M-1)*NANG+(IJ-1)*NFRANG 
                IF(PSEC4(IC).NE.ZMISS) THEN
                  PSEC4OUT(IJ+INDEXN)=PSEC4OUT(IJ+INDEXN)
     &                               +PSEC4(IC)*DFIM(M)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
          DO IJ=1,ILEN1
            PSEC4OUT(IJ+INDEXN)=PSEC4OUT(IJ+INDEXN)*DELTH
            IF(PSEC4OUT(IJ+INDEXN).LE.0.) PSEC4OUT(IJ+INDEXN)=ZMISS
          ENDDO

!         ADD TAIL (IF LAST FREQUENCY IS CONSIDERED)
 
          IF(IFRE_LAST.EQ.NFRE) THEN
            IF (ITEST.gt.0) WRITE (*,*) ' ADDING TAIL TO SPECTRUM' 
            DELT25 = FR(NFRE)/4.*DELTH
            M=IFRE_LAST
            DO K=1,NANG
              DO IJ=1,ILEN1
                IC=K+(M-1)*NANG+(IJ-1)*NFRANG
                IF(PSEC4(IC).NE.ZMISS.AND.
     &             PSEC4OUT(IJ+INDEXN).NE.ZMISS) THEN
                    PSEC4OUT(IJ+INDEXN)=PSEC4OUT(IJ+INDEXN)
     &                                  +DELT25*PSEC4(IC)
                ENDIF
              ENDDO
            ENDDO
          ENDIF

          IF(LWAVEHGT) THEN
            IF (ITEST.gt.0)
     &      WRITE (*,*)' integrated spectra is converted to wave height'
            DO IJ=1,ILEN1
              IF(PSEC4OUT(IJ+INDEXN).NE.ZMISS)
     &           PSEC4OUT(IJ+INDEXN)=4*SQRT(PSEC4OUT(IJ+INDEXN))
            ENDDO
          ENDIF

          DO IJ=1,INDEXN
             PSEC4OUT(IJ)=ZMISS
          ENDDO

          DO IJ=1,INDEXS
             PSEC4OUT(IJ+INDEXN+ILEN1)=ZMISS
          ENDDO

          IF (ITEST.gt.0) THEN 
            WRITE(*,*)' CODING OUTPUT FIELD ' 
            WRITE(*,*)' CORRESPONDING TO FREQUENCY FROM ',
     &      FLOAT(ISEC1(49+IFRE_FIRST))/ISEC1(49)
            WRITE(*,*)' TO',
     &      FLOAT(ISEC1(49+IFRE_LAST))/ISEC1(49)
          ENDIF
          IF (ITEST.gt.1) THEN 
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
     &         PSEC4OUT, ILENO,OUTGRIB, ILENGOUT ,KWORD,'C',KRET)
          IF (ITEST.GT.1) WRITE(*,*)' GRIBEX DONE status=' , KRET 
          IF(KRET.GT.0) THEN
            MSG = ERRMSG(6)
            CALL GRPRS0 (KSEC0)
            CALL GRPRS1 (KSEC0, KSEC1)
            WRITE (*,*) MSG
            CALL ABORT
          ENDIF

          KOUNT=KWORD*NPRECI
          CALL PBWRITE(IUOUT,OUTGRIB,KOUNT,KRET)
          IF(KRET.LT.0) THEN
            WRITE (*,*) '**************************************'
            WRITE (*,*) '*  ERROR FOLLOWING CALL TO PBWRITE   *'
            WRITE (*,*) '**************************************'
            CALL ABORT
          ENDIF

        ELSE

        IF (ITEST.gt.0) WRITE (*,*) ' COMPUTING 1-D spectra from 250'

          DO M=IFRE_FIRST,IFRE_LAST

            ISEC1(44) = 0
            ISEC1(45) = M


            ISEC3(2)=0
            ZSEC3(2)=ZMISS

            PSEC4OUT=0.
            DO K=1,NANG
              DO IJ=1,ILEN1
                IC=K+(M-1)*NANG+(IJ-1)*NFRANG 
                IF(PSEC4(IC).NE.ZMISS) THEN
                  PSEC4OUT(IJ+INDEXN)=PSEC4OUT(IJ+INDEXN)+PSEC4(IC)
                ENDIF
              ENDDO
            ENDDO

            DO IJ=1,ILEN1
              PSEC4OUT(IJ+INDEXN)=PSEC4OUT(IJ+INDEXN)*DELTH
              IF(PSEC4OUT(IJ+INDEXN).LE.0.) PSEC4OUT(IJ+INDEXN)=ZMISS
            ENDDO

            DO IJ=1,INDEXN
               PSEC4OUT(IJ)=ZMISS
            ENDDO

            DO IJ=1,INDEXS
               PSEC4OUT(IJ+INDEXN+ILEN1)=ZMISS
            ENDDO

            IF (ITEST.gt.0) THEN 
              WRITE(*,*)' CODING OUTPUT FIELD ',M 
              WRITE(*,*)' CORRESPONDING TO FREQUENCY ',
     &        FLOAT(ISEC1(49+M))/ISEC1(49)
            ENDIF
            IF (ITEST.gt.1) THEN 
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
            IF (ITEST.GT.1) WRITE(*,*)' GRIBEX DONE status=' , KRET 
            IF(KRET.GT.0) THEN
              MSG = ERRMSG(6)
              CALL GRPRS0 (KSEC0)
              CALL GRPRS1 (KSEC0, KSEC1)
              WRITE (*,*) MSG
              CALL ABORT
            ENDIF

            KOUNT=KWORD*NPRECI
            CALL PBWRITE(IUOUT,OUTGRIB,KOUNT,KRET)
            IF(KRET.LT.0) THEN
              WRITE (*,*) '**************************************'
              WRITE (*,*) '*  ERROR FOLLOWING CALL TO PBWRITE   *'
              WRITE (*,*) '**************************************'
              CALL ABORT
            ENDIF
          ENDDO
        ENDIF

        IF(ALLOCATED(OUTGRIB)) DEALLOCATE(OUTGRIB)
        IF(ALLOCATED(PSEC4OUT)) DEALLOCATE(PSEC4OUT)

!     CONVERSION FROM 251 to 1-D spectra
!     -------------------

      ELSE IF(IPARAM.EQ.251) THEN

        ILENGOUT=2*(1000+NGX*NGYMAX)/NPRECI
        ALLOCATE(OUTGRIB(ILENGOUT))

        ISEC0 = KSEC0
        ISEC1 = KSEC1
        ISEC2 = KSEC2
        ISEC3 = KSEC3
        ISEC4 = KSEC4
        ZSEC2 = PSEC2
        ZSEC3 = PSEC3

        ISEC1( 6) =  252

        ISEC1(46) = 0 

        DO IC=1,NFRE
          ISEC1(49+0+IC) = ISEC1(49+NANG+IC)
        ENDDO

        ILENO = ILEN1

        ALLOCATE(PSEC4OUT(ILENO))

        ISEC4(1) = ILENO 

        IF(LINTEGRATE.OR.LWAVEHGT) THEN

        IF (ITEST.gt.0)
     &      WRITE (*,*) ' COMPUTING integrated 1-D spectra from 251'

          ISEC1(44) = 0
          ISEC1(45) = IFRE_FIRST 

          ISEC3(2)=0
          ZSEC3(2)=ZMISS

          PSEC4OUT=0.
          DO M=IFRE_FIRST,IFRE_LAST
            DO IJ=1,ILENO
              PSEC4OUT(IJ)=PSEC4OUT(IJ)+DFIM(M)*PSEC41D(IJ,M)*DELTH
            ENDDO
          ENDDO

          DO IJ=1,ILENO
            IF(PSEC4OUT(IJ).LE.0.) PSEC4OUT(IJ)=ZMISS
          ENDDO

!         ADD TAIL (IF LAST FREQUENCY IS CONSIDERED)
 
          IF(IFRE_LAST.EQ.NFRE) THEN
            IF (ITEST.gt.0) WRITE (*,*) ' ADDING TAIL TO SPECTRUM' 
            DELT25 = FR(NFRE)/4.*DELTH
            DO IJ=1,ILENO
              IF(PSEC4OUT(IJ).NE.ZMISS)
     &        PSEC4OUT(IJ)=PSEC4OUT(IJ)+DELT25*PSEC41D(IJ,IFRE_LAST)
            ENDDO
          ENDIF

          IF(LWAVEHGT) THEN
            IF (ITEST.gt.0)
     &      WRITE (*,*)' integrated spectra is converted to wave height'
            DO IJ=1,ILENO
              IF(PSEC4OUT(IJ).NE.ZMISS)
     &        PSEC4OUT(IJ)=4*SQRT(PSEC4OUT(IJ))
            ENDDO
          ENDIF

          IF (ITEST.gt.0) THEN 
            WRITE(*,*)' CODING OUTPUT FIELD ' 
            WRITE(*,*)' CORRESPONDING TO FREQUENCY FROM ',
     &      FLOAT(ISEC1(49+IFRE_FIRST))/ISEC1(49)
            WRITE(*,*)' TO',
     &      FLOAT(ISEC1(49+IFRE_LAST))/ISEC1(49)
          ENDIF
          IF (ITEST.gt.1) THEN 
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
     &         PSEC4OUT, ILENO,OUTGRIB, ILENGOUT ,KWORD,'C',KRET)
          IF (ITEST.GT.1) WRITE(*,*)' GRIBEX DONE status=' , KRET 
          IF(KRET.GT.0) THEN
            MSG = ERRMSG(6)
            CALL GRPRS0 (KSEC0)
            CALL GRPRS1 (KSEC0, KSEC1)
            WRITE (*,*) MSG
            CALL ABORT
          ENDIF

          KOUNT=KWORD*NPRECI
          CALL PBWRITE(IUOUT,OUTGRIB,KOUNT,KRET)
          IF(KRET.LT.0) THEN
            WRITE (*,*) '**************************************'
            WRITE (*,*) '*  ERROR FOLLOWING CALL TO PBWRITE   *'
            WRITE (*,*) '**************************************'
            CALL ABORT
          ENDIF

        ELSE

        IF (ITEST.gt.0) WRITE (*,*) ' COMPUTING 1-D spectra from 251'

          DO M=IFRE_FIRST,IFRE_LAST

            ISEC1(44) = 0
            ISEC1(45) = M

            ISEC3(2)=0
            ZSEC3(2)=ZMISS

            DO IJ=1,ILENO
              PSEC4OUT(IJ)=PSEC41D(IJ,M)*DELTH
              IF(PSEC4OUT(IJ).LE.0.) PSEC4OUT(IJ)=ZMISS
            ENDDO

            IF (ITEST.gt.0) THEN 
              WRITE(*,*)' CODING OUTPUT FIELD ',M 
              WRITE(*,*)' CORRESPONDING TO FREQUENCY ',
     &        FLOAT(ISEC1(49+M))/ISEC1(49)
            ENDIF
            IF (ITEST.gt.1) THEN 
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
            IF (ITEST.GT.1) WRITE(*,*)' GRIBEX DONE status=' , KRET 
            IF(KRET.GT.0) THEN
              MSG = ERRMSG(6)
              CALL GRPRS0 (KSEC0)
              CALL GRPRS1 (KSEC0, KSEC1)
              WRITE (*,*) MSG
              CALL ABORT
            ENDIF

            KOUNT=KWORD*NPRECI
            CALL PBWRITE(IUOUT,OUTGRIB,KOUNT,KRET)
            IF(KRET.LT.0) THEN
              WRITE (*,*) '**************************************'
              WRITE (*,*) '*  ERROR FOLLOWING CALL TO PBWRITE   *'
              WRITE (*,*) '**************************************'
              CALL ABORT
            ENDIF
          ENDDO
        ENDIF

        IF(ALLOCATED(OUTGRIB)) DEALLOCATE(OUTGRIB)
        IF(ALLOCATED(PSEC4OUT)) DEALLOCATE(PSEC4OUT)

      ENDIF

      CALL PBCLOSE(IUOUT,KRET)
      IF(KRET.LT.0) THEN 
        WRITE (*,*) '************************************'
        WRITE (*,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
        WRITE (*,*) '************************************'
      ENDIF

      END

!#######################################################################
!     FUNCTION getclo(yaoptions, yaargument)
!     INTEGER getclo, getcla, rtb
!     CHARACTER*   1 yolastarg
!     CHARACTER* (*) yaoptions, yaargument
!     CHARACTER* 120 arg
!
!     INTEGER here, imorearg, ivarg
!     DATA here, imorearg, ivarg, arg / 1, 0, 0, "  " /
!     DATA yolastarg / " " /
!     arg=' '
!     CALL getarg(here,arg)
!
!     iol=rtb(arg)
!     IF (iol .EQ. 2 .AND. arg(1:1) .EQ. '-' .AND. ivarg .EQ. 0 ) THEN
!       iol = rtb(yaoptions)
!       DO jl=1,iol
!         getclo = 0
!         IF ( yaoptions(jl:jl) .EQ. arg(2:2) ) THEN
!           getclo = ichar(arg(2:2))
!           IF (yaoptions(jl+1:jl+1) .EQ. ':' ) THEN
!             yolastarg=yaoptions(jl:jl)
!             ivarg=1
!           ENDIF
!           EXIT
!         ENDIF
!       ENDDO
!     ELSEIF ( ivarg .EQ. 1 ) THEN
!        WRITE(*,*) ' option -', yolastarg, ' requires arguments'
!        getclo=-1
!     ELSEIF (iol .EQ. 0) THEN
!       getclo=0
!     ELSE
!        WRITE(*,*) 'illegal option: ',arg(1:iol)
!        getclo=-1
!     ENDIF
!     here = here + 1
!     RETURN
!     ENTRY getcla(yaargument)
!     getcla = 1
!     CALL getarg(here,arg)
!     IF ( arg (1:1) .NE. '-' ) THEN
!       here = here + 1
!       yaargument=arg
!     ELSE
!       IF (ivarg.EQ.1) THEN
!         WRITE(*,*)' refused to take ', arg (1:2) ,' as argument for',
!    .    ' the option -',yolastarg
!         getcla = -1
!       ELSE
!         getcla = 0
!       ENDIF
!     ENDIF
!     ivarg=0
!     RETURN
!     END
