C
C ----------------------------------------------------------------------
C
      PROGRAM CONVERT_IRISHGRB
C
C     J BIDLOT      ECMWF      JUNE 1998
C
C     PURPOSE.
C     --------
C     CONVERTS THE IRISH MET SERVICE GRIB WAVE SPECTRA 
C     FROM PARAMETER DEFINITION 251 TO 250
C  
C     USAGE :
C     ------  
C            input_filename : input.dat
C            output_filename : output.dat
C
C     LIBRARY : EMOSLIB
C     -------
C
C     NOTE: IT IS ASSUMED THAT THE PARAMETER 251 FIELDS WERE PRODUCED
C     ----  ON 3 SUBAREAS.
C           IF OTHER PARAMETERS ARE PRESENT IN THE INPUT FILE, THEY
C           WILL BE UNTOUCHED AND WILL BE WRITTEN OUT TO THE OUTPUT
C           FILE FIRST. 
C
C ----------------------------------------------------------------------
C
      PARAMETER(NANGMAX=12,NFREMAX=25)
      PARAMETER(Nxymax=100*50)
      PARAMETER(IPLENG=6500)
C
      INTEGER I4
      INTEGER IAREA(6)
      INTEGER KSEC0(2),KSEC1(128),KSEC2(384),KSEC3(2),KSEC4(659)
      INTEGER ISEC0(2),ISEC1(128),ISEC2(384),ISEC3(2),ISEC4(659)
      INTEGER INGRIB(IPLENG),OUTGRIB(IPLENG*NFREMAX*NANGMAX)
C
      REAL ZTHETA,ZFRE
      REAL PSEC2(96),PSEC3(2)
      REAL ZSEC2(96),ZSEC3(2)
      REAL PSEC4(Nxymax), PSEC4OUT(Nxymax*NFREMAX*NANGMAX)
C
C ----------------------------------------------------------------------
C
      ZMISS=-999.
      I4=1
      NPRECI = KIND(I4) 
      IAREA=-999
      ICOMPTOT=1
C
      CALL PBOPEN(IUOUT,'output.dat','w',KRET)
      IF(KRET.LT.0) THEN
          IF(KRET.EQ.-1)
     &      WRITE (*,*) ' ERROR; COULD NOT OPEN FILE ','output.dat'
          IF(KRET.EQ.-2)
     &      WRITE (*,*) ' ERROR; INVALID FILENAME ','output.dat'
          IF(KRET.EQ.-3)
     &       WRITE (*,*) ' ERROR; INVALID OPEN MODE SPECIFIED'
          CALL ABORT
      ENDIF
C
C    LOOP 3 times to extract the 3 subdomains
c
  200 CONTINUE
C
      KSEC1=0
      KSEC2=0
      KSEC4=0
      KSEC3(1)=0
      KSEC3(2)=0
      PSEC2=0.
      PSEC3(2)=ZMISS
      PSEC4OUT=ZMISS
C
      CALL PBOPEN(IUGRS,'input.dat','r',KRET)
      IF(KRET.LT.0) THEN
          IF(KRET.EQ.-1)
     &      WRITE (*,*) ' ERROR; COULD NOT OPEN FILE ','input.dat'
          IF(KRET.EQ.-2)
     &      WRITE (*,*) ' ERROR; INVALID FILENAME ','input.dat'
          IF(KRET.EQ.-3)
     &       WRITE (*,*) ' ERROR; INVALID OPEN MODE SPECIFIED'
          CALL ABORT
      ENDIF
C
      ICOMPE=0
C
C     Loop on all entry from input file 
C
  100 CONTINUE
C
      CALL PBGRIB(IUGRS,INGRIB,IPLENG,ILENG,KRET)
C
      IF (KRET.EQ.-1) THEN
        WRITE (*,*) ' REACHED EOF IN ','input.dat' 
        CALL PBCLOSE(IUGRS,KRET)
C
        IF(ICOMPE.NE.NANG*NFRE) THEN
          WRITE (*,*) 'ERROR  THE NUMBER OF PARAMETER 251 FIELDS' 
          WRITE (*,*) 'ICOMPE = ',ICOMPE
          WRITE (*,*) 'IS NOT EQUAL TO THE EXPECTED NUMBER'
          WRITE (*,*) 'NANG*NFRE = ',NANG*NFRE
          WRITE (*,*) 'PROGRAM WILL ABORT '
          CALL ABORT()
        ENDIF
C
        WRITE (*,*) ' parameter 250 will be created'
C
        ISEC1( 6) =  250
        ISEC1(37) = 1 
        ISEC1(24) = 0 ! No local ext.
        ISEC3(2) = 0
        ZSEC3(2) = ZMISS
        ISEC4( 8) = 64
        ISEC4( 9) = 32
        ISEC4(50) = NANG
        ISEC4(51) = NFRE
C
        DO IC=1,NANG
           ZTHETA = FLOAT(ISEC1(49+IC))/ISEC1(48)
           ISEC4(59+IC) = TRANSFER(ZTHETA,1)
        ENDDO
C
        DO IC=1,NFRE
           ZFRE = FLOAT(ISEC1(49+NANG+IC))/ISEC1(49)
           ISEC4(59+NANG+IC) = TRANSFER(ZFRE,1)
        ENDDO
C
        ILENGOUT=IPLENG*NFREMAX*NANGMAX
C
        ILENO = NXY*NANG*NFRE
        ISEC4(1) = ILENO 
C
        IF( ILENO.GT.Nxymax*NFREMAX*NANGMAX) THEN
           WRITE(*,*) ' ERROR; PSEC4OUT TO SMALL',ILENO,
     &                  Nxymax*NFREMAX*NANGMAX
           CALL ABORT()
        ENDIF
C
        CALL GRIBEX(ISEC0, ISEC1, ISEC2, ZSEC2, ISEC3, ZSEC3, ISEC4,
     &         PSEC4OUT, ILENO, OUTGRIB, ILENGOUT ,KWORD,'C',KRET)
        IF(KRET.GT.0) THEN
          CALL GRPRS0 (ISEC0)
          CALL GRPRS1 (ISEC0, ISEC1)
          WRITE (*,*) ' ENCODING ERROR SEE GRIBEX DESCRIPTION *' 
          CALL ABORT()
        ENDIF
C
        WRITE(*,*) ' writing the results to output file '
        KOUNT=KWORD*NPRECI
        CALL PBWRITE(IUOUT,OUTGRIB,KOUNT,KRET)
        IF(KRET.LT.0) THEN
          WRITE (*,*) '**************************************'
          WRITE (*,*) '*  ERROR FOLLOWING CALL TO PBWRITE   *'
          WRITE (*,*) '**************************************'
          CALL ABORT
        ENDIF
C 
        IF (ICOMPTOT.EQ.1) THEN
          IAREA(3)=IAREA(1)
          IAREA(4)=IAREA(2)
        ELSE IF (ICOMPTOT.EQ.2) THEN
          IAREA(5)=IAREA(1)
          IAREA(6)=IAREA(2)
        ENDIF
C
        ICOMPTOT=ICOMPTOT+1
        IF(ICOMPTOT.LE.3) GOTO 200
C
        CALL PBCLOSE(IUOUT,KRET)
        STOP
C
      ELSEIF (KRET .EQ. -2) THEN
        WRITE (*,*) ' ERROR IN FILE HANDLING IN ','input.dat' 
      ELSEIF (KRET .EQ. -3) THEN
        WRITE (*,*) ' INGRIB IS NOT BIG ENOUGH. CURRENT SIZE: ',IPLENG
        WRITE (*,*) ' REQUIRED SIZE SHOULD BE: ',ILENG
      ENDIF
C
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
C
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
C
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
C
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
C
      GO TO 100
C
      END

