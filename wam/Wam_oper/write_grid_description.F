!
! ----------------------------------------------------------------------
!
      PROGRAM WRITE_GRIB_DESCRIPTION 
!
!     J BIDLOT ECMWF
!
!     PURPOSE.
!     --------
!     writes out the information about the grid found in the input file

!!! this program was adapted from a more complicated program, it might
!!! still contains too much !!!!

!     USAGE : write_grid_description [-i input_filename]
!                                    [-t itest] 
!    
!             default value for input_filename : input 
!                               output_filename : output 
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
!     The input file can only contain one parameter
 
!     LIBRARY : FORTRAN 90, and EMOSLIB
!     -------

!     OUTPUT FILE FORMAT :
!     --------------------


! ----------------------------------------------------------------------

      PARAMETER(NANGMAX=24,NFREMAX=25)
      PARAMETER (G = 9.806, PI = 3.1415927, CIRC = 40000000.,
     1           ZPI = 2.*PI, RAD = PI/180., DEG = 180./PI,
     2           R = CIRC/ZPI, EPSMIN=0.1E-32, ZMISS=-999.0)

      INTEGER I4
      INTEGER GETCLO, GETCLA, IOPTVAL
      INTEGER :: KSEC0(2),KSEC1(51+NANGMAX+NFREMAX),KSEC2(1024),
     &           KSEC3(2),KSEC4(128)
      INTEGER, ALLOCATABLE :: INGRIB(:),OUTGRIB(:),KDOMRGG(:)
      REAL :: ONE,ZTHETA,ZFRE
      REAL :: PSEC2(96),PSEC3(2)
      REAL, ALLOCATABLE :: PSEC4OUT(:)
      REAL, ALLOCATABLE :: FIELD(:,:)
      CHARACTER     CLOPTLET
      CHARACTER*  3 CLL1
      CHARACTER*  6 CLOPTS
      CHARACTER* 12 CLFMT
      CHARACTER*40 ERRMSG(-4:6), MSG
      CHARACTER*128 CLARG, FNAMEIN, FNAMEOUT
      LOGICAL LLEXIST, LINTEGRATE, LWAVEHGT, LLASTFRE, llfrst

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
     5            ' error  *',
     6            ' ENCODING ERROR SEE GRIBEX DESCRIPTION *'/

! ----------------------------------------------------------------------

!*    INITIAL VALUES SET AND CRACK COMMAND LINE.
!     -----------------------------------------

      I4=1
      NPRECI = KIND(I4) 
      FNAMEIN='input'
      FNAMEOUT='output'
      IFRE_FIRST=1
      IFRE_LAST=NFREMAX
      LINTEGRATE=.FALSE.
      LWAVEHGT=.FALSE.
      LLASTFRE=.FALSE.
      ITEST=0
      llfrst=.true.

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
          ELSE IF ( CLOPTLET .EQ. 't' ) THEN
            I1=LEN_TRIM(CLARG)
            WRITE (CLL1,'(I3)') I1
            CLFMT = '(I'//CLL1//')'
            READ (CLARG(1:I1),FMT=CLFMT) ITEST
          ENDIF
        ENDIF
      ENDDO CMDLINE

      KSEC1=0
      KSEC2=0
      KSEC4=0
      KSEC3(1)=0
      KSEC3(2)=0
      PSEC2=0.
      PSEC3(2)=ZMISS

!*    INPUT GRIB DATA FILE
!     --------------------
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
        WRITE(*,*)'*GRIB DATA NOT FOUND IN *'
        WRITE(*,*)  FNAMEIN 
        WRITE(*,*)'*PROGRAM WILL ABORT        *'
        WRITE(*,*)'*                          *'
        WRITE(*,*)'****************************'
        CALL ABORT
      ENDIF
 
!     GET SIZE OF INGRIB


      CALL PBSIZE (IUGRS, IPLENG)
      ISIZE=(IPLENG+NPRECI-1)/NPRECI
1111  IPLENG=ISIZE*NPRECI
      IF(.NOT.ALLOCATED(INGRIB)) ALLOCATE(INGRIB(ISIZE))

!     GET DATA FILE

      CALL PBGRIB(IUGRS,INGRIB,IPLENG,ILENG,KRET)
      IF     (KRET.EQ.-1) THEN
        WRITE (*,*) ' REACHED EOF IN ',FNAMEIN 
        GOTO 2222
      ELSEIF (KRET .EQ. -2) THEN
        WRITE (*,*) ' ERROR IN FILE HANDLING IN ',FNAMEIN 
      ELSEIF (KRET .EQ. -3) THEN
        DEALLOCATE(INGRIB)

        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' SIZE  OF KGRIB IS NOT BIG ENOUGH.'
        WRITE(IU06,*) ' IT WAS ', ISIZE

        KKRET=0

        ISIZE=(ILENG+KIND(ISIZE)-1)/KIND(ISIZE)
        WRITE(IU06,*) ' IT SHOULD AT LEAST BE ', ISIZE
        WRITE(IU06,*) ' THE SIZE WAS RESET AUTOMATICALLY'
        WRITE(IU06,*) ' '
        CALL FLUSH(IU06)

!       RESET THE FILE POINTER TO READ FIELD AGAIN
        KOFFSET=-ILENG
        CALL PBSEEK(IUGRS,KOFFSET,1,KKRET)
        IF(KKRET.EQ.-1) THEN
          WRITE(IU06,*) '***********************************'
          WRITE(IU06,*) '*  PBSEEK : END OF FILE ENCOUNTED'
          WRITE(IU06,*) '***********************************'
          CALL ABORT
        ENDIF
        IF(KKRET.EQ.-2) THEN
          WRITE(IU06,*) '***********************************'
          WRITE(IU06,*) '*  PBSEEK : FILE HANDLING ERROR'
          WRITE(IU06,*) '***********************************'
          CALL ABORT
        ENDIF

        GOTO 1111

      ENDIF

!*    GET GRIB HEADERS 

      KRET  = 1
      ILEN1 = 1
      ALLOCATE(PSEC4OUT(ILEN1))
      IF (ITEST.gt.0)
     &     WRITE(*,*)' GETTING GRIB HEADER OF 1st INPUT FIELD' 
      IF (ITEST.GT.2) CALL GRSDBG (1)
      CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &             PSEC4OUT, ILEN1, INGRIB, ILENG , KWORD, 'J', KRET)
      IF (ITEST.GT.1) WRITE(*,*)' GRIBEX DONE status=' , KRET 
      IF(KRET.GT.0) THEN
        MSG = ERRMSG(2)
        CALL GRPRS0 (KSEC0)
        CALL GRPRS1 (KSEC0, KSEC1)
        WRITE (*,*) MSG
        CALL ABORT
      ENDIF
      IF(ALLOCATED(PSEC4OUT)) DEALLOCATE(PSEC4OUT)

!*    DETERMINE DATA FIELD CHARACTERISTICS 

      IPARAM = KSEC1(6)
      IREPR = KSEC2(1)
      IRGG = KSEC2(17)
      AMONOP = FLOAT(KSEC2(4)/1000)+0.1*(MOD(KSEC2(4),1000)/100)+
     &         0.01*(MOD(KSEC2(4),100)/10)+0.001*MOD(KSEC2(4),10)
      AMOSOP = FLOAT(KSEC2(7)/1000)+0.1*(MOD(KSEC2(7),1000)/100)+
     &         0.01*(MOD(KSEC2(7),100)/10)+0.001*MOD(KSEC2(7),10)
      AMOWEP = FLOAT(KSEC2(5)/1000)+0.1*(MOD(KSEC2(5),1000)/100)+
     &         0.01*(MOD(KSEC2(5),100)/10)+0.001*MOD(KSEC2(5),10)

      AMOEAP = FLOAT(KSEC2(8)/1000)+0.1*(MOD(KSEC2(8),1000)/100)+
     &         0.01*(MOD(KSEC2(8),100)/10)+0.001*MOD(KSEC2(8),10)
      CALL ADJUST (AMOWEP, AMOEAP)

      IYYYY=(KSEC1(21)-1)*100+KSEC1(10)
      IMM=KSEC1(11)
      IDD=KSEC1(12)
      IHH=KSEC1(13)
      IMI=KSEC1(14)

      NGY = KSEC2(3)
      IF (ITEST.gt.0) WRITE(*,*) ' THE INPUT PARAMETER IS ',IPARAM 

      ALLOCATE(KDOMRGG(NGY))
      IF(IRGG.EQ.0) THEN
        NGX = KSEC2(2)
        KDOMRGG = NGX
        XDELLO = FLOAT(KSEC2(9))/1000
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
        XDELLO=(AMOEAP-AMOWEP)/MAX(1,NGX-1)
      ENDIF

      IF(IREPR.EQ.4) THEN
!       for Gaussian grid XDELLA is the averaged one !
        XDELLA = (AMONOP-AMOSOP)/(NGY-1)
!       Gaussian grid are assumed to be periodic !!!
        IPER = 1
      ELSE
        XDELLA = FLOAT(KSEC2(10))/1000
        IPER = 0
        IF (AMOEAP-AMOWEP+1.5*XDELLO.GE.360.) IPER = 1
      ENDIF




!     output
!     *********
      write(*,*) 'AMONOP = ',AMONOP
      write(*,*) 'AMOSOP = ',AMOSOP
      write(*,*) 'AMOWEP =',AMOWEP
      write(*,*) 'AMOEAP =',AMOEAP
      write(*,*) 'XDELLA =',XDELLA
      write(*,*) 'NGY =',NGY

      open(21,file='grid_description',form='formatted')
      write(21,'(i4)') NGY-1
      write(21,'(f8.3)') AMONOP
      write(21,'(f8.3)') AMOSOP
      write(21,'(f8.3)') AMOWEP
      write(21,'(f8.3)') AMOEAP
      write(21,'(i4)') IPER 
      write(21,'(i4)') IRGG 
      write(21,'(i4)') NGY
      DO J=1,NGY
        write(21,'(i4)') KDOMRGG(J)
      ENDDO

      DEALLOCATE(KDOMRGG)
!!! only read the first field      GOTO 1111


2222  CONTINUE
      CALL PBCLOSE(IUGRS,KRET)
      IF(KRET.LT.0) THEN
        WRITE (*,*) '************************************'
        WRITE (*,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
        WRITE (*,*) '* FILE ',FNAMEIN
        WRITE (*,*) '************************************'
        CALL ABORT
      ENDIF

      CLOSE(10)

      END

!#######################################################################
      FUNCTION getclo(yaoptions, yaargument)
      INTEGER getclo, getcla, rtb
      CHARACTER*   1 yolastarg
      CHARACTER* (*) yaoptions, yaargument
      CHARACTER* 120 arg

      INTEGER here, imorearg, ivarg
      DATA here, imorearg, ivarg, arg / 1, 0, 0, "  " /
      DATA yolastarg / " " /

      arg=' '
      CALL getarg(here,arg)
!-->  PRINT*,'-------------getclo---------------'
!-->  PRINT*,' ###',arg,'###'

      iol=rtb(arg)
!-->  PRINT*,' iol: ', iol
      IF (iol .EQ. 2 .AND. arg(1:1) .EQ. '-' .AND. ivarg .EQ. 0 ) THEN
        iol = rtb(yaoptions)
!-->    PRINT*,' iol 2',iol,' options: ', yaoptions
        DO jl=1,iol
          getclo = 0
          IF ( yaoptions(jl:jl) .EQ. arg(2:2) ) THEN
            getclo = ichar(arg(2:2))
!-->        PRINT *,' FOUND ',yaoptions(jl:jl), ' IN THE COMMAND LINE',
!--> .               yaoptions(jl+1:jl+1)
            IF (yaoptions(jl+1:jl+1) .EQ. ':' ) THEN
!-->          PRINT*, yaoptions(jl:jl),' requires arguments'
              yolastarg=yaoptions(jl:jl)
              ivarg=1
            ENDIF
            EXIT
          ENDIF
        ENDDO
      ELSEIF ( ivarg .EQ. 1 ) THEN
         WRITE(*,*) ' option -', yolastarg, ' requires arguments'
         getclo=-1
      ELSEIF (iol .EQ. 0) THEN
        getclo=0
      ELSE
         WRITE(*,*) 'illegal option: ',arg(1:iol)
         getclo=-1
      ENDIF
!-->  PRINT *,' HERE ins getclo', here
      here = here + 1
      RETURN

      ENTRY getcla(yaargument)
!-->  PRINT*,'-------------getcla--------------'
!-->  PRINT*, 'HERE ins getcla :', here,' options: ', yaoptions

      getcla = 1
      CALL getarg(here,arg)
!-->  PRINT*,' arg in getcla ', arg
      IF ( arg (1:1) .NE. '-' ) THEN
        here = here + 1
        yaargument=arg
      ELSE
        IF (ivarg.EQ.1) THEN
          WRITE(*,*)' refused to take ', arg (1:2) ,' as argument for',
     .    ' the option -',yolastarg
          getcla = -1
        ELSE
          getcla = 0
        ENDIF
      ENDIF
      ivarg=0
!-->  PRINT*,' getcla in getcla ', getcla

      RETURN
      END
      SUBROUTINE ADJUST(WEST, EAST)

! ----------------------------------------------------------------------

!**** *ADJUST* - ROUTINE TO CORRECT BORDERS OF INTERVALS.

!     H.GUNTHER            ECMWF       04/04/1990

!*    PURPOSE.
!     -------

!       ADJUSTS INTERVAL BORDERS GIVEN IN DEGREE.

!**   INTERFACE.
!     ----------

!       *CALL* *ADJUST (WEST, EAST)*
!          *WEST*    - LEFT INTERVAL BORDER IN DEGREE.
!          *EAST*    - RIGHT INTERVAL BORDER IN DEGREE.

!     METHOD.
!     -------

!       THE INTERVAL BORDERS ARE CHANGED TO FULLFILL:
!         0. .LE. EAST  .AND. EAST .LT. 360. .AND. WEST .LE. EAST

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

!* 1. CORRECT BORDERS.
!     ----------------

      WEST = MOD(WEST+720.,360.)
      EAST = MOD(EAST+720.,360.)
      IF (WEST.GT.EAST) WEST = WEST-360.

      RETURN
      END SUBROUTINE ADJUST






