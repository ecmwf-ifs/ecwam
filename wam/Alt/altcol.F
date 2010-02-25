!
! this simple program is intended to produce a collocation between
! model and altimeter gridded observations.
! The fields should be retrieved on their initial resolution
! and therefore have the same resolution for all of them.
!
      INTEGER, PARAMETER :: NBIT = 430000
      INTEGER, PARAMETER :: NKSEC0=2
      INTEGER, PARAMETER :: NKSEC1=2048
      INTEGER, PARAMETER :: NKSEC2=2070
      INTEGER, PARAMETER :: NKSEC3=2
      INTEGER, PARAMETER :: NKSEC4=252

      INTEGER, PARAMETER :: NPSEC2=96
      INTEGER, PARAMETER :: NPSEC3=2

      INTEGER,ALLOCATABLE :: KGRIB(:)

      INTEGER, PARAMETER :: NPARAM=12
      INTEGER :: IUNIT(nparam),jparam(nparam),jrgg(nparam),kr(nparam)
      CHARACTER*3 :: CPARAM
      CHARACTER*24 :: FILE(nparam)
      CHARACTER*12 :: CDATE(nparam)
      LOGICAL :: lfrt
      LOGICAL :: lpresent(nparam)

      INTEGER :: KSEC0(NKSEC0), KSEC1(NKSEC1), KSEC2(NKSEC2),
     &           KSEC3(NKSEC3), KSEC4(NKSEC4)
      REAL ::  PSEC2(NPSEC2),PSEC3(NPSEC3)
      INTEGER,ALLOCATABLE :: RLONRGG(:)
      REAL, ALLOCATABLE :: PSEC4(:), RDELLO(:), FIELD(:,:,:)

      iu06=6      

!     list of the collocated parameters (only 229,246,247 are
!     always required)
      jparam(1)=246   ! uncorrected alt Hs (must be provided)
      jparam(2)=247   ! corrected alt Hs (must be provided)
      jparam(3)=248   ! alt relative sea state bias correction
      jparam(4)=229   ! Hs (must be prodived)
      jparam(5)=221   ! Tz
      jparam(6)=231   ! Tp
      jparam(7)=230   ! mean wave direction
      jparam(8)=244   ! mean square slope
      jparam(9)=245   ! U10
      jparam(10)=234  ! Hs windsea
      jparam(11)=236  ! Tm windsea
      jparam(12)=235  ! mean wind sea direction


      open(11,file='collocation_file',form='formatted')

      lfrt=.true.
      ZMISS=-999.
      KSEC3(1) = 0
      PSEC3(2) = ZMISS
      lpresent = .false.


!     open files

      do ifile=1,nparam
        write(cparam,'(i3.3)') jparam(ifile)
        file(ifile)='param_'//cparam 
        INQUIRE (FILE=file(ifile), EXIST=lpresent(ifile))
        if(lpresent(ifile)) then
          CALL PBOPEN(IUNIT(ifile),file(ifile),'r',KRET)
          IF(KRET.LT.0) THEN
            IF(KRET.EQ.-1) THEN
               WRITE(IU06,*) '*  could not open '//FILE(ifile)
               if(ifile.eq.1 .or. ifile.eq.2 .or. ifile.eq.4) then
                 WRITE(IU06,*) '*  this file is always needed !!!'
                 WRITE(IU06,*) '*  the program must abort !!!'
                 call abort
               else
                 WRITE(IU06,*) '*  this file is optional !!!'
                 WRITE(IU06,*) '*  parameter ',cparam,
     &                         ' will not be collocated'
                 lpresent(ifile)=.false.
               endif
            ENDIF
            IF(KRET.EQ.-2) THEN
              WRITE(IU06,*) '***********************************'
              WRITE(IU06,*) '*                                 *'
              WRITE(IU06,*) '*  INVALID FILENAME'//FILE(ifile)
              WRITE(IU06,*) '*                                 *'
              WRITE(IU06,*) '***********************************'
              CALL ABORT
             ENDIF
          ENDIF
        else
          WRITE(IU06,*) '*  file ',file(ifile),' is optional !!!'
          WRITE(IU06,*) '*  parameter ',cparam,
     &                  ' will not be collocated'
        endif
      enddo

!     read in and decode files

      do while (.true.)

      do ifile=1,nparam
        if(lpresent(ifile)) then
        ISIZE=NBIT
1001    KPLENG=ISIZE*KIND(ISIZE)
        IF(.NOT.ALLOCATED(KGRIB)) ALLOCATE(KGRIB(ISIZE))
        CALL PBGRIB(IUNIT(ifile),KGRIB,KPLENG,KLENPB,KRET)
        IF(KRET.EQ.-1) THEN
!         end of file 1. It means that all has been read
!         or nothing. In either case you must move on
          goto 2000 
        ENDIF
        IF(KRET.EQ.-2) THEN
          WRITE(IU06,*) '***********************************'
          WRITE(IU06,*) '*                                 *'
          WRITE(IU06,*) '*  PBGRIB : FILE HANDLING ERROR'
          WRITE(IU06,*) '*                                 *'
          WRITE(IU06,*) '***********************************'
          CALL ABORT
        ENDIF
        IF(KRET.EQ.-3) THEN
          WRITE(IU06,*) ' '
          WRITE(IU06,*) ' ***** WARNING IN MODEL FOR WAVES *****'
          WRITE(IU06,*) ' SIZE  OF KGRIB IS NOT BIG ENOUGH.'
          WRITE(IU06,*) ' IT WAS ', ISIZE

          KRET=0
          DEALLOCATE(KGRIB)
          ISIZE=(KLENPB+KIND(ISIZE)-1)/KIND(ISIZE)
          WRITE(IU06,*) ' IT SHOULD AT LEAST BE ', ISIZE
          WRITE(IU06,*) ' THE SIZE WAS RESET AUTOMATICALLY'
          WRITE(IU06,*) ' AND THE FIELD READ WITH THE NEW SIZE'
          WRITE(IU06,*) ' IF THIS PROBLEM OCCURS TOO OFTEN'
          WRITE(IU06,*) ' MODIFY THE VALUE OF NBIT IN SOURCE'
          WRITE(IU06,*) ' ***** WARNING ****** WARNING *****'
          WRITE(IU06,*) ' '

!         RESET THE FILE POINTER TO READ FIELD AGAIN
          KOFFSET=-KLENPB
          CALL PBSEEK(IUNIT(ifile),KOFFSET,1,KRET)
          IF(KRET.EQ.-1) THEN
            WRITE(IU06,*) '***********************************'
            WRITE(IU06,*) '*                                 *'
            WRITE(IU06,*) '*  PBSEEK : END OF FILE ENCOUNTED'
            WRITE(IU06,*) '*                                 *'
            WRITE(IU06,*) '***********************************'
            CALL ABORT
          ENDIF
          IF(KRET.EQ.-2) THEN
            WRITE(IU06,*) '***********************************'
            WRITE(IU06,*) '*                                 *'
            WRITE(IU06,*) '*  PBSEEK : FILE HANDLING ERROR'
            WRITE(IU06,*) '*                                 *'
            WRITE(IU06,*) '***********************************'
            CALL ABORT
          ENDIF

          GOTO 1001
        ENDIF

!       DECODE DATA

        if(lfrt) then
!         GET SIZE OF THE UNPACKED DATA
          ALLOCATE(PSEC4(1))
          KLENG = ISIZE
          KLENP = 1
          KRET  = 1
!         CALL GRSDBG(1)
          CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &                 PSEC4, KLENP, KGRIB, KLENG , KWORD, 'I', KRET)
          IF(KRET.GT.0) THEN
            WRITE(IU06,*) '***************************************'
            WRITE(IU06,*) '*                                     *'
            WRITE(IU06,*) '*  GRIBEX FAILED WITH ERROR CODE ',KRET
            WRITE(IU06,*) '*                                     *'
            WRITE(IU06,*) '*     THE PROGRAM ABORTS              *'
            WRITE(IU06,*) '***************************************'
            CALL ABORT
          ENDIF                

          NR=KSEC2(3)
          IF(NR.GT.NKSEC2-22) THEN
            WRITE(IU06,*) '***************************************'
            WRITE(IU06,*) '*                                     *'
            WRITE(IU06,*) '*  THE DIMENSION FOR KSEC2 ',NR+22
            WRITE(IU06,*) '*  IS LARGER THAN THE PRESCRIBED ONE ',NKSEC2
            WRITE(IU06,*) '*                                     *'
            WRITE(IU06,*) '*     THE PROGRAM ABORTS              *'
            WRITE(IU06,*) '***************************************'
            CALL ABORT
          ENDIF
          IRGG = KSEC2(17)
          JRGG(ifile) = IRGG

          IF(irgg.EQ.1) THEN
            NC=0
            NTOT=0
            DO J=1,NR
            NC = MAX(NC,KSEC2(22+J))
            NTOT=NTOT+KSEC2(22+J)
            ENDDO

            IR=0
            DO J=1,NR
              IF(KSEC2(22+J).NE.0) IR=IR+1
            ENDDO
            NR=IR
            KLEN=NTOT
          ELSEIF(irgg.EQ.0) THEN
            NC=KSEC2(2)
            KLEN = NC*NR
          ELSE
            WRITE(IU06,*)
     &  '   REPRESENTATION OF THE FIELD NOT KNOWN'
            CALL ABORT
          ENDIF

          IREPR=KSEC2(1)

          IF(IREPR.NE.0.AND.IREPR.NE.4) THEN
            WRITE(IU06,*) '***************************************'
            WRITE(IU06,*) '*                                     *'
            WRITE(IU06,*) '*  UNKNOWN GRID REPRESENTATION = ',IREPR
            WRITE(IU06,*) '*  CAN ONLY DEAL WITH                 *'
            WRITE(IU06,*) '*  LATITUDE/LONGITUDE GRID (IREPR=0)  *'
            WRITE(IU06,*) '*   OR GAUSSIAN (IREPR=4)             *'
            WRITE(IU06,*) '*                                     *'
            WRITE(IU06,*) '*     THE PROGRAM ABORTS              *'
            WRITE(IU06,*) '***************************************'
            CALL ABORT
          ENDIF

          DEALLOCATE(PSEC4)
          lfrt=.false.

        endif

!       GET THE DATA
        KLENG = ISIZE
        KLENP = KLEN
        IF(.NOT.ALLOCATED(PSEC4)) ALLOCATE(PSEC4(KLENP))
        KRET  = 1

        CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &               PSEC4, KLENP, KGRIB, KLENG, KWORD, 'D', KRET)
        IF(KRET.GT.0) THEN
          WRITE(IU06,*) '***************************************'
          WRITE(IU06,*) '*                                     *'
          WRITE(IU06,*) '*                                     *'
          WRITE(IU06,*) '*  GRIBEX FAILED WITH ERROR CODE ',KRET
          WRITE(IU06,*) '*                                     *'
          WRITE(IU06,*) '*     THE PROGRAM ABORTS              *'
          WRITE(IU06,*) '***************************************'
          CALL ABORT
        ENDIF

!*      DETERMINE INFORMATION ABOUT THE DECODED DATA

!       START DATE.
        IYYYY=(KSEC1(21)-1)*100+KSEC1(10)
        WRITE(CDATE(ifile),'(I4.4,4I2.2)') IYYYY,KSEC1(11),KSEC1(12),
     &                             KSEC1(13),KSEC1(14)

!*      DETERMINE TIME UNIT.
        IF (KSEC1(15).EQ.0) THEN
          JCONS = 60
        ELSEIF (KSEC1(15).EQ.1) THEN
          JCONS = 3600
        ELSEIF (KSEC1(15).EQ.2) THEN
          JCONS = 86400
        ELSE
          WRITE(IU06,*) '***********************************'
          WRITE(IU06,*) '*  UNABLE TO DECODE THE TIME UNIT *'
          WRITE(IU06,*) '***********************************'
          CALL ABORT
        ENDIF

!       DETERMINE FORECAST PERIOD.
        IFORP = JCONS*KSEC1(16)

!       GET DATE OF THE FORECAST INSTEAD OF STARTING DATE
        IF(KSEC1(39).NE.1) THEN
          CALL INCDATE (CDATE(ifile),IFORP)
        ENDIF

!*      DETERMINE CODE FOR DATA FIELD TYPE.
        IPARAM = KSEC1(6)
        if(IPARAM.ne.jparam(ifile)) then
          WRITE(IU06,*) '***********************************'
          WRITE(IU06,*) '*  WRONG PARAMETER ORDER. ABORT !!!*'
          WRITE(IU06,*) '***********************************'
          CALL ABORT
        endif

        IF(KSEC1(7).EQ.105) THEN
          KZLEV=KSEC1(8)
        ELSE
          KZLEV=0
        ENDIF

!*      DETERMINE GRID PARAMETERS.
        NR=KSEC2(3)
        kr(ifile)=NR
        IRGG = KSEC2(17)
        JRGG(ifile) = IRGG

        IF(.NOT.ALLOCATED(RLONRGG)) ALLOCATE(RLONRGG(NR))
        RLONRGG=0
        IF(.NOT.ALLOCATED(RDELLO)) ALLOCATE(RDELLO(NR))

        RMONOP = REAL(KSEC2(4)/1000)+0.1*(MOD(KSEC2(4),1000)/100)+
     &           0.01*(MOD(KSEC2(4),100)/10)+0.001*MOD(KSEC2(4),10)

        RMOWEP = REAL(KSEC2(5)/1000)+0.1*(MOD(KSEC2(5),1000)/100)+
     &           0.01*(MOD(KSEC2(5),100)/10)+0.001*MOD(KSEC2(5),10)

        RMOSOP = REAL(KSEC2(7)/1000)+0.1*(MOD(KSEC2(7),1000)/100)+
     &           0.01*(MOD(KSEC2(7),100)/10)+0.001*MOD(KSEC2(7),10)

        RMOEAP = REAL(KSEC2(8)/1000)+0.1*(MOD(KSEC2(8),1000)/100)+
     &           0.01*(MOD(KSEC2(8),100)/10)+0.001*MOD(KSEC2(8),10)


        IF(irgg.EQ.1) THEN

          ISTART=0
          DO WHILE(KSEC2(23+ISTART).EQ.0)
            ISTART=ISTART+1
          ENDDO

          ISTOP=0
          DO WHILE(KSEC2(23+NR-1-ISTOP).EQ.0)
            ISTOP=ISTOP+1
          ENDDO

          IR=0
          DO J=1,NR
            IF(KSEC2(22+J).NE.0) IR=IR+1
          ENDDO
          NR=IR

          DO J=1,NR
            JSN=NR-J+1
            RLONRGG(JSN) = KSEC2(22+J+ISTART)
          ENDDO

          IDUM = KSEC2(4)-ISTART*KSEC2(10)
          RMONOP = REAL(IDUM/1000)+0.1*(MOD(IDUM,1000)/100)+
     &             0.01*(MOD(IDUM,100)/10)+0.001*MOD(IDUM,10)

          IDUM = KSEC2(7)+ISTOP*KSEC2(10)
          RMOSOP = REAL(IDUM/1000)+0.1*(MOD(IDUM,1000)/100)+
     &           0.01*(MOD(IDUM,100)/10)+0.001*MOD(IDUM,10)

        ELSEIF(irgg.EQ.0) THEN
          RLONRGG=NC
        ELSE
          WRITE(IU06,*)
     &    ' REPRESENTATION OF THE FIELD NOT KNOWN'
          CALL ABORT
        ENDIF

        PMISS=PSEC3(2)

        IPERIODIC = 0
        DELLO=(RMOEAP-RMOWEP)/MAX(1,NC-1)
        IF (RMOEAP-RMOWEP+DELLO.GE.360.) IPERIODIC = 1

        DELLA=(RMONOP-RMOSOP)/MAX(1,NR-1)

        IF(IPERIODIC.EQ.1) THEN
          DELLO=(RMOEAP-RMOWEP)/MAX(1,NC)
          DO J=1,NR
            JSN=NR-J+1
            RDELLO(JSN) = (RMOEAP-RMOWEP+DELLO)/MAX(1,RLONRGG(JSN))
          ENDDO
        ELSE
          DELLO=(RMOEAP-RMOWEP)/MAX(1,NC-1)
          DO J=1,NR
            JSN=NR-J+1
            RDELLO(JSN) = (RMOEAP-RMOWEP)/MAX(1,RLONRGG(JSN)-1)
          ENDDO
        ENDIF


        IF(.NOT.ALLOCATED(FIELD)) ALLOCATE(FIELD(NC,NR,nparam))

        L = 0
        DO K=1,NR
          JSN=NR-K+1
          DO I=1,RLONRGG(JSN)
            L = L+1
            FIELD(I,K,IFILE) = PSEC4(L)
          ENDDO
        ENDDO
        endif
      enddo

!     flag missing parameters
      do ifile=1,nparam
        if(.not.lpresent(ifile)) then
          DO K=1,NR
            JSN=NR-K+1
            DO I=1,RLONRGG(JSN)
              FIELD(I,K,IFILE) = ZMISS 
            ENDDO
          ENDDO
        endif
      enddo

!     COLLOCATIONS

      do ip=1,nparam
        if(lpresent(ip)) then
        do jp=ip+1,nparam
          if(lpresent(jp)) then
            if(kr(jp).ne.kr(ip)) then
            WRITE(IU06,*) '*************************************'
            WRITE(IU06,*) '* NR MISMATCHED !!!'
            WRITE(IU06,*) ip, KR(ip),'for parameter ',jparam(ip) 
            WRITE(IU06,*) jp, KR(jp),'for parameter ',jparam(jp) 
            WRITE(IU06,*) '*************************************'
            CALL ABORT
          endif
          endif
        enddo
        endif
      enddo

      do ip=1,nparam
        if(lpresent(ip)) then
        do jp=ip+1,nparam
          if(lpresent(jp)) then
            if(cdate(ip).ne.cdate(jp)) then
              WRITE(IU06,*) '***********************************'
              WRITE(IU06,*) '* DATE AND TIME MISMATCHED !!!'
              WRITE(IU06,*) cdate(1) 
              WRITE(IU06,*) cdate(2) 
              WRITE(IU06,*) '***********************************'
              CALL ABORT
            endif
          endif
        enddo
        endif
      enddo

      do ip=1,nparam
        if(lpresent(ip)) then
        do jp=ip+1,nparam
          if(lpresent(jp)) then
            if(jrgg(ip).ne.jrgg(jp)) then
              WRITE(IU06,*) '***********************************'
              WRITE(IU06,*) '* GRID TYPE MISMATCHED !!!'
              WRITE(IU06,*) jrgg(1) 
              WRITE(IU06,*) jrgg(2) 
              WRITE(IU06,*) '***********************************'
              CALL ABORT
            endif
          endif
        enddo
        endif
      enddo

      write(*,*) 'COLLOCATION FOR ', CDATE(1) 

        DO K=1,NR
          JSN=NR-K+1
          DO I=1,RLONRGG(JSN)
            if(FIELD(I,K,1).ne.PMISS .and. FIELD(I,K,1).gt.0 .and.
     &         FIELD(I,K,2).ne.PMISS .and. FIELD(I,K,2).gt.0 .and.
     &         FIELD(I,K,4).ne.PMISS .and. FIELD(I,K,4).gt.0 ) then
              xlon=RMOWEP+(I-1)*RDELLO(JSN)
              ylat=RMOSOP + (JSN-1)*DELLA
              write(11,'(a12,1x,24(f11.6,1x))') cdate(1),xlon,ylat,
     &                  (FIELD(I,K,iparam),iparam=1,nparam)
            endif
          ENDDO
        ENDDO

      enddo

2000  continue

      close(11)

      end
