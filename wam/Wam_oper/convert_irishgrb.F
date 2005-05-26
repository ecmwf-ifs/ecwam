
! ----------------------------------------------------------------------

      PROGRAM CONVERT_IRISHGRB

!     J BIDLOT      ECMWF      JUNE 1998

!     PURPOSE.
!     --------
!     CONVERTS THE IRISH MET SERVICE GRIB WAVE SPECTRA 
!     FROM PARAMETER DEFINITION 251 TO 250

!     USAGE :
!     ------  
!            input_filename : input.dat
!            output_filename : output.dat

!     LIBRARY : EMOSLIB
!     -------

!     NOTE: IT IS ASSUMED THAT THE PARAMETER 251 FIELDS WERE PRODUCED
!     ----  ON 3 SUBAREAS.
!           IF OTHER PARAMETERS ARE PRESENT IN THE INPUT FILE, THEY
!           WILL BE UNTOUCHED AND WILL BE WRITTEN OUT TO THE OUTPUT
!           FILE FIRST. 

! ----------------------------------------------------------------------

      PARAMETER(NANGMAX=24,NFREMAX=36)
!     PARAMETER(Nxymax=100*50)
      PARAMETER(Nxymax=300*80)
      PARAMETER(IPLENG=10000)

      INTEGER I4
      INTEGER IAREA(6)
      INTEGER KSEC0(2),KSEC1(255),KSEC2(2060),KSEC3(2),KSEC4(659)
      INTEGER ISEC0(2),ISEC1(255),ISEC2(2060),ISEC3(2),ISEC4(659)
      INTEGER INGRIB(IPLENG),OUTGRIB(IPLENG*NFREMAX*NANGMAX)

      REAL ZTHETA,ZFRE
      REAL PSEC2(96),PSEC3(2)
      REAL ZSEC2(96),ZSEC3(2)
      REAL PSEC4(Nxymax), PSEC4OUT(Nxymax*NFREMAX*NANGMAX)

! ----------------------------------------------------------------------

      ZMISS=-999.
      I4=1
      NPRECI = KIND(I4) 
      IAREA=-999
      ICOMPTOT=1

      CALL PBOPEN(IUOUT,'output.dat','w',KRET)
      IF(KRET.LT.0) THEN
        IF(KRET.EQ.-1)
     &   WRITE (*,*) ' ERROR; COULD NOT OPEN FILE ','output.dat'
        IF(KRET.EQ.-2)
     &   WRITE (*,*) ' ERROR; INVALID FILENAME ','output.dat'
        IF(KRET.EQ.-3)
     &   WRITE (*,*) ' ERROR; INVALID OPEN MODE SPECIFIED'
        CALL ABORT
      ENDIF

!    LOOP 3 times to extract the 3 subdomains

  200 CONTINUE

      KSEC1=0
      KSEC2=0
      KSEC4=0
      KSEC3(1)=0
      KSEC3(2)=0
      PSEC2=0.
      PSEC3(2)=ZMISS
      PSEC4OUT=ZMISS

      CALL PBOPEN(IUGRS,'input.dat','r',KRET)
      IF(KRET.LT.0) THEN
        IF(KRET.EQ.-1)
     &   WRITE (*,*) ' ERROR; COULD NOT OPEN FILE ','input.dat'
        IF(KRET.EQ.-2)
     &   WRITE (*,*) ' ERROR; INVALID FILENAME ','input.dat'
        IF(KRET.EQ.-3)
     &   WRITE (*,*) ' ERROR; INVALID OPEN MODE SPECIFIED'
        CALL ABORT
      ENDIF

      ICOMPE=0

!     Loop on all entry from input file 

  100 CONTINUE

      CALL PBGRIB(IUGRS,INGRIB,IPLENG,ILENG,KRET)

      IF (KRET.EQ.-1) THEN
        WRITE (*,*) ' REACHED EOF IN ','input.dat' 
        CALL PBCLOSE(IUGRS,KRET)

        IF(ICOMPE.NE.NANG*NFRE) THEN
          WRITE (*,*) 'ERROR  THE NUMBER OF PARAMETER 251 FIELDS' 
          WRITE (*,*) 'ICOMPE = ',ICOMPE
          WRITE (*,*) 'IS NOT EQUAL TO THE EXPECTED NUMBER'
          WRITE (*,*) 'NANG*NFRE = ',NANG*NFRE
          WRITE (*,*) 'PROGRAM WILL ABORT '
          CALL ABORT()
        ENDIF

        WRITE (*,*) ' parameter 250 will be created'

        ISEC1( 6) =  250
        ISEC1(37) = 1 
        ISEC1(24) = 0 ! No local ext.
        ISEC3(2) = 0
        ZSEC3(2) = ZMISS
        ISEC4( 8) = 64
        ISEC4( 9) = 32
        ISEC4(50) = NANG
        ISEC4(51) = NFRE

        DO IC=1,NANG
          ZTHETA = REAL(ISEC1(49+IC))/ISEC1(48)
          ISEC4(59+IC) = TRANSFER(ZTHETA,1)
        ENDDO

        DO IC=1,NFRE
          ZFRE = REAL(ISEC1(49+NANG+IC))/ISEC1(49)
          ISEC4(59+NANG+IC) = TRANSFER(ZFRE,1)
        ENDDO

        ILENGOUT=IPLENG*NFREMAX*NANGMAX

        ILENO = NXY*NANG*NFRE
        ISEC4(1) = ILENO 

        IF( ILENO.GT.Nxymax*NFREMAX*NANGMAX) THEN
          WRITE(*,*) ' ERROR; PSEC4OUT TO SMALL',ILENO,
     &     Nxymax*NFREMAX*NANGMAX
          CALL ABORT()
        ENDIF

        CALL GRIBEX(ISEC0, ISEC1, ISEC2, ZSEC2, ISEC3, ZSEC3, ISEC4,
     &         PSEC4OUT, ILENO, OUTGRIB, ILENGOUT ,KWORD,'C',KRET)
        IF(KRET.GT.0) THEN
          CALL GRPRS0 (ISEC0)
          CALL GRPRS1 (ISEC0, ISEC1)
          WRITE (*,*) ' ENCODING ERROR SEE GRIBEX DESCRIPTION *' 
          CALL ABORT()
        ENDIF

        WRITE(*,*) ' writing the results to output file '
        KOUNT=KWORD*NPRECI
        CALL PBWRITE(IUOUT,OUTGRIB,KOUNT,KRET)
        IF(KRET.LT.0) THEN
          WRITE (*,*) '**************************************'
          WRITE (*,*) '*  ERROR FOLLOWING CALL TO PBWRITE   *'
          WRITE (*,*) '**************************************'
          CALL ABORT
        ENDIF

        IF (ICOMPTOT.EQ.1) THEN
          IAREA(3)=IAREA(1)
          IAREA(4)=IAREA(2)
        ELSEIF (ICOMPTOT.EQ.2) THEN
          IAREA(5)=IAREA(1)
          IAREA(6)=IAREA(2)
        ENDIF

        ICOMPTOT=ICOMPTOT+1
        IF(ICOMPTOT.LE.2) GOTO 200

        CALL PBCLOSE(IUOUT,KRET)
        STOP

      ELSEIF (KRET .EQ. -2) THEN
        WRITE (*,*) ' ERROR IN FILE HANDLING IN ','input.dat' 
      ELSEIF (KRET .EQ. -3) THEN
        WRITE (*,*) ' INGRIB IS NOT BIG ENOUGH. CURRENT SIZE: ',IPLENG
        WRITE (*,*) ' REQUIRED SIZE SHOULD BE: ',ILENG
      ENDIF

      KSEC3(2)=0
      PSEC3(2)=ZMISS
      ILEN1=Nxymax
      CALL GRIBEX(KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &            PSEC4, ILEN1, INGRIB, IPLENG , KWORD, 'D', KRET)
      IF (KRET.GT.0) THEN
        CALL GRPRS0 (KSEC0)
        CALL GRPRS1 (KSEC0, KSEC1)
        WRITE (*,*) ' DECODING ERROR SEE GRIBEX DESCRIPTION *' 
        CALL ABORT
      ENDIF

      IF (KSEC1(6).ne.251) THEN
        IF (ICOMPTOT.EQ.1) THEN
          KOUNT=ILENG
          CALL PBWRITE(IUOUT,INGRIB,KOUNT,KRET)
          IF(KRET.LT.0) THEN
            WRITE (*,*) '**************************************'
            WRITE (*,*) '*  ERROR FOLLOWING CALL TO PBWRITE   *'
            WRITE (*,*) '*  FOR PARAMETER ',KSEC1(6)
            WRITE (*,*) '**************************************'
            CALL ABORT
          ENDIF
        ENDIF
        GO TO 100
      ENDIF

      IF (ICOMPTOT.EQ.1.AND.ICOMPE.EQ.0) THEN
        print*,'SUBDOMAIN :'
        print*,' KSEC2(4) = ',KSEC2(4),' KSEC2(5) = ',KSEC2(5)
        print*,' KSEC2(7) = ',KSEC2(7),' KSEC2(8) = ',KSEC2(8)
        IAREA(1)=KSEC2(4)
        IAREA(2)=KSEC2(5)
        NXY=KSEC2(2)*KSEC2(3)
        If (NXY.GT.NXYMAX) THEN
          WRITE(*,*) 'DATA ARRAY TOO SMALL.'
          CALL ABORT()
        ENDIF
      ELSEIF (ICOMPTOT.EQ.2.AND.ICOMPE.EQ.0) THEN
        IF(KSEC2(4).EQ.IAREA(3).AND.KSEC2(5).EQ.IAREA(4)) GOTO 100
        print*,'SUBDOMAIN :'
        print*,' KSEC2(4) = ',KSEC2(4),' KSEC2(5) = ',KSEC2(5)
        print*,' KSEC2(7) = ',KSEC2(7),' KSEC2(8) = ',KSEC2(8)
        IAREA(1)=KSEC2(4)
        IAREA(2)=KSEC2(5)
        NXY=KSEC2(2)*KSEC2(3)
        If (NXY.GT.NXYMAX) THEN
          WRITE(*,*) 'DATA ARRAY TOO SMALL.'
          CALL ABORT()
        ENDIF
      ELSEIF (ICOMPTOT.EQ.3.AND.ICOMPE.EQ.0) THEN
        IF(KSEC2(4).EQ.IAREA(3).AND.KSEC2(5).EQ.IAREA(4)) GOTO 100
        IF(KSEC2(4).EQ.IAREA(5).AND.KSEC2(5).EQ.IAREA(6)) GOTO 100
        print*,'SUBDOMAIN :'
        print*,' KSEC2(4) = ',KSEC2(4),' KSEC2(5) = ',KSEC2(5)
        print*,' KSEC2(7) = ',KSEC2(7),' KSEC2(8) = ',KSEC2(8)
        IAREA(1)=KSEC2(4)
        IAREA(2)=KSEC2(5)
        NXY=KSEC2(2)*KSEC2(3)
        If (NXY.GT.NXYMAX) THEN
          WRITE(*,*) 'DATA ARRAY TOO SMALL.'
          CALL ABORT()
        ENDIF
      ENDIF

      IF(KSEC2(4).EQ.IAREA(1).AND.KSEC2(5).EQ.IAREA(2)) THEN
        ICOMPE=ICOMPE+1
        ISEC0 = KSEC0
        ISEC1 = KSEC1
        ISEC2 = KSEC2
        ISEC3 = KSEC3
        ISEC4 = KSEC4
        ZSEC2 = PSEC2
        ZSEC3 = PSEC3
        KK = KSEC1(44)
        MM = KSEC1(45)
        NANG = KSEC1(46)
        NFRE = KSEC1(47)
        NFRANG = NANG*NFRE
        DO IJ=1,Nxy
          PSEC4OUT(KK+(MM-1)*NANG+(IJ-1)*NFRANG) =  PSEC4(IJ)
        ENDDO
      ENDIF

      GO TO 100

      END


