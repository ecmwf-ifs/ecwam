      SUBROUTINE INITIALINT (IU06, ITEST, NCA, NRA, KSEC2,
     &                       NGX, NGY, KRGG, KLONRGG, XDELLA, ZDELLO,
     &                       AMOWEP, AMOSOP, AMOEAP, AMONOP, IPERIODIC,
     &                       RLONRGG, 
     &                       LLINTERPOL,DK1, DII1, DIIP1, KK, II, IIP)
! ----------------------------------------------------------------------    

!***  *INITIALINT* - INITIALIZES THE INTERPOLATION PARAMETERS TO BE USED
!                    IN INTERPOLATING FIELDS PASSED FROM THE ATMOSPHIRIC 
!                    SIDE TO THE WAVE SIDE OF THE MODEL.

!      S. ABDALLA  ECMWF  DECEMBER 2000. 

!     PURPOSE.                                                          
!     --------                                                          

!     IT COMPUTES THE WEIGHT VALUES AND THE 4 GRID POINTS TO BE USED
!     FOR INTERPOLATING THE FIELDS PASSED FROM THE ATMOSPHIRIC SIDE TO  
!     THE WAVE SIDE OF THE MODEL. "INITIALINT" MUST BE CALLED BEFORE
!     ANY CALL TO "FLDINTER".

!     LIMITATIONS.
!     ------------
!     1. IT IS ASSUMED THAT THE FIELDS ARE NOT ANGULAR (DIRECTIONAL).
!     2. IT IS ASSUMED THAT THERE IS NO MISSING DATA POINTS.
!     3. IT ONLY WORKS IF THE INPUT AND OUTPUT GRIDS ARE LATITUDE/LONGITUDE
!        (REGULAR OR IRREGULAR) OR GAUSSIAN GRIDS (FULL OR REDUCED) !!!!

!**   INTERFACE.                                                        
!     ----------                                                        

!      *CALL* *INITIALINT* (IU06, ITEST, NCA, NRA, KSEC2,
!    &                      NGX, NGY, KRGG, KLONRGG, XDELLA, ZDELLO,
!    &                      AMOWEP, AMOSOP, AMOEAP, AMONOP, IPERIODIC,
!    &                      RLONRGG,
!    &                      LLINTERPOL,DK1, DII1, DIIP1, KK, II, IIP)

!        *IU06*   - OUTPUT UNIT.
!        *ITEST*  - TEST MESSAGE LEVEL.
!       ATMOSPHERIC MODEL GRID:
!        *NCA*    - NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR
!        *NRA*    - NUMBER OF ATM. ROWS OF LATITUDES
!        *KSEC2*  - GRIB SECTION 2 TO EXTRACT ATM. GRID CHARACTERISTICS.
!        *RLONRGG - NUMBER OF GRID POINTS FOR EACH LATITUDE (ATM. GRID)
!       WAVE MODEL GRID SPECIFICATION:
!        *NGX  *  - NUMBER OF COLUMNS IN ARRAY FIELD USED.              
!        *NGY  *  - NUMBER OF ROWS    IN ARRAY FIELD USED.              
!        *KRGG*   - GRID DEFINITION PARAMETER (0=REGULAR, 1=IRREGULAR)
!        *KLONRGG - NUMBER OF GRID POINTS FOR EACH LATITUDE
!        *XDELLA* - GRID POINT SPACING BETWEEN LATITUDES.
!        *ZDELLO* - GRID POINT SPACING PER LATITUDES. 
!        *AMOWEP* - MOST WESTERN LONGITUDE IN GRID (  1, ? ).
!        *AMOSOP* - MOST SOUTHERN LATITUDE IN GRID.( ? ,NGY).
!        *AMOEAP* - MOST EASTERN LONGITUDE IN GRID (NGX, ? ).
!        *AMONOP* - MOST NORTHERN LATITUDE IN GRID ( ? , 1 ).
!       OUTPUT: INTERPOLATION COEFFICIENTS & INDECES
!        *IPERIODIC* - SPECIFIES IF ATM. GRID IS PERIODIC (=1) OR NOT (=0)
!        *LLINTERPOL* - FLAG (TRUE=DO INTERPOLATION, FALSE=ASSIGN ONLY)
!        *DK1M*   - COEFFICIENT
!        *DII1M*  - COEFFICIENT
!        *DIIP1M* - COEFFICIENT
!        *KKM*    - INDEX
!        *IIM*    - INDEX
!        *IIPM*   - INDEX

!     EXTERNALS.                                                        
!     ----------                                                        

!        none
! ----------------------------------------------------------------------

      INTEGER NCA, NRA
      INTEGER KSEC2(22+NRA)

      INTEGER KLONRGG(NGY)
      INTEGER RLONRGG(NRA)

      REAL ZDELLO(NGY)
      REAL, ALLOCATABLE :: RDELLO(:)

      INTEGER  KK(NGY), II(NGX,NGY), IIP(NGX,NGY)
      REAL DK1(NGY), DII1(NGX,NGY), DIIP1(NGX,NGY)

      LOGICAL LLINTERPOL, LLNEAREST
                                                                        
! ----------------------------------------------------------------------

      NR=KSEC2(3)
      JRGG = KSEC2(17)

      IF(JRGG.EQ.1) THEN
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
      ELSEIF(JRGG.EQ.0) THEN
        NC=KSEC2(2)
        KLEN = NC*NR
      ELSE
        WRITE(IU06,*)
     &  ' SUB INITIALINT : REPRESENTATION OF THE FIELD NOT KNOWN'
        CALL ABORT1
      ENDIF

      IF (NC.NE.NCA .OR. NR.NE.NRA) THEN
        WRITE(IU06,*) '***************************************'
        WRITE(IU06,*) '***************************************'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*    WARNING   IN SUB. INITIALINT     *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '* NC and NCA are NOT equal:           *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,'(A,I6,A,I6,A)')
     &             '*  NC = ', NC, '   and    NCA = ',NCA,'  *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*              OR                     *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '* NR and NRA are NOT equal:           *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,'(A,I6,A,I6,A)')
     &             '*  NR = ', NR, '   and    NRA = ',NRA,'  *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '***************************************'
      ENDIF

      LLNEAREST=.FALSE.

      IREPR=KSEC2(1)

      IF(IREPR.NE.0.AND.IREPR.NE.4) THEN
        WRITE(IU06,*) '***************************************'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*  FATAL ERROR IN SUB. INITIALINT     *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*  UNKNOWN GRID REPRESENTATION = ',IREPR
        WRITE(IU06,*) '*  GRB2WGRD CAN ONLY DEAL WITH        *'
        WRITE(IU06,*) '*  LATITUDE/LONGITUDE GRID (IREPR=0)  *'
        WRITE(IU06,*) '*   OR GAUSSIAN (IREPR=4)              *' 
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*     THE PROGRAM ABORTS              *'
        WRITE(IU06,*) '***************************************'
        CALL ABORT1
      ENDIF


! ----------------------------------------------------------------------

!*    DETERMINE INFORMATION ABOUT THE DECODED DATA 
!     --------------------------------------------
!*    DETERMINE GRID PARAMETERS.                                    
!     --------------------------                                    

      RLONRGG=0
      ALLOCATE(RDELLO(NR))

      RMONOP = REAL(KSEC2(4)/1000)+0.1*(MOD(KSEC2(4),1000)/100)+
     &         0.01*(MOD(KSEC2(4),100)/10)+0.001*MOD(KSEC2(4),10)

      RMOWEP = REAL(KSEC2(5)/1000)+0.1*(MOD(KSEC2(5),1000)/100)+
     &         0.01*(MOD(KSEC2(5),100)/10)+0.001*MOD(KSEC2(5),10)

      RMOSOP = REAL(KSEC2(7)/1000)+0.1*(MOD(KSEC2(7),1000)/100)+
     &         0.01*(MOD(KSEC2(7),100)/10)+0.001*MOD(KSEC2(7),10)

!!! THERE IS A DANGER THAT THE DEFINITON FOR RMOEAP MIGHT VARY DUE TO
!!! THE AMBIGOUS DEFINITION FOR IRREGULAR GRIDS. FOR NON WAVE FIELDS,
!!! A GAUSSIAN GRID IMPLIES THAT THE GRID IS GLOBAL, THEREFORE
!!! RMOEAP IS IMPLICITLY KNOWN.
      IF(IREPR.EQ.4 .AND. .NOT. LLNEAREST) THEN
        DELLO = 360./MAX(1,NC)
        RMOEAP = RMOWEP+360. - DELLO
        IPERIODIC = 1
      ELSE
        RMOEAP = REAL(KSEC2(8)/1000)+0.1*(MOD(KSEC2(8),1000)/100)+
     &           0.01*(MOD(KSEC2(8),100)/10)+0.001*MOD(KSEC2(8),10)
        IF(JRGG.EQ.1) CALL ADJUST (RMOWEP, RMOEAP)
        IPERIODIC = 0
        DELLO=(RMOEAP-RMOWEP)/MAX(1,NC-1)
        IF (RMOEAP-RMOWEP+1.5*DELLO.GE.360.) IPERIODIC = 1
      ENDIF


      IF(JRGG.EQ.1) THEN

        ISTART=0
        DO WHILE(KSEC2(23+ISTART).EQ.0)
          ISTART=ISTART+1 
        ENDDO

        ISTOP=0
        DO WHILE(KSEC2(23+NR-1-ISTOP).EQ.0)
          ISTOP=ISTOP+1 
        ENDDO

        DO J=1,NR
          JSN=NR-J+1
          RLONRGG(JSN) = KSEC2(22+J+ISTART) 
        ENDDO

        IDUM = KSEC2(4)-ISTART*KSEC2(10)
        RMONOP = REAL(IDUM/1000)+0.1*(MOD(IDUM,1000)/100)+
     &           0.01*(MOD(IDUM,100)/10)+0.001*MOD(IDUM,10)

        IDUM = KSEC2(7)+ISTOP*KSEC2(10)
        RMOSOP = REAL(IDUM/1000)+0.1*(MOD(IDUM,1000)/100)+
     &         0.01*(MOD(IDUM,100)/10)+0.001*MOD(IDUM,10)

      ELSEIF(JRGG.EQ.0) THEN
        RLONRGG=NC
      ELSE
        WRITE(IU06,*)
     &  ' SUB INITIALINT: REPRESENTATION OF THE FIELD NOT KNOWN'
        CALL ABORT1
      ENDIF


!     FIND WHETHER INTERPOLATION IS NEEDED

      DELLA=(RMONOP-RMOSOP)/MAX(1,NR-1)

      KAMOWEP=NINT(AMOWEP*100.)
      KAMOEAP=NINT(AMOEAP*100.)
      KAMONOP=NINT(AMONOP*100.)
      KAMOSOP=NINT(AMOSOP*100.)
      KRMOWEP=NINT(RMOWEP*100.)
      KRMOEAP=NINT(RMOEAP*100.)
      KRMONOP=NINT(RMONOP*100.)
      KRMOSOP=NINT(RMOSOP*100.)

      IF(IPERIODIC.EQ.1) THEN
        DELLO=360./MAX(1,NC)
        DO J=1,NR
          JSN=NR-J+1
          RDELLO(JSN) = 360./MAX(1,RLONRGG(JSN)) 
        ENDDO
      ELSE
        DELLO=(RMOEAP-RMOWEP)/MAX(1,NC-1)
        DO J=1,NR
          JSN=NR-J+1
          RDELLO(JSN) = (RMOEAP-RMOWEP)/MAX(1,RLONRGG(JSN)-1) 
        ENDDO
      ENDIF

      IF(KAMONOP.GT.KRMONOP.OR.KAMONOP.LT.KRMOSOP.OR.
     &   KAMOSOP.LT.KRMOSOP.OR.KAMOSOP.GT.KRMONOP.OR.
     &   (JRGG.EQ.0.AND.IREPR.NE.4.AND.IPERIODIC.NE.1
     &    .AND.KAMOWEP.LT.KRMOWEP).OR.
     &   (JRGG.EQ.0.AND.IREPR.NE.4.AND.IPERIODIC.NE.1
     &    .AND.KAMOEAP.GT.NINT((RMOEAP+DELLO)*100.)) ) THEN

         WRITE(IU06,*) '***********************************'
         WRITE(IU06,*) '*                                 *'
         WRITE(IU06,*) '*  FATAL ERROR IN SUB. INITIALINT *'
         WRITE(IU06,*) '*  ============================== *'
         WRITE(IU06,*) '*                                 *'
         WRITE(IU06,*) '*  THE MODEL DOMAIN IS OUTSIDE    *' 
         WRITE(IU06,*) '*  THE INPUT DOMAIN FOR FIELDS    *'
         WRITE(IU06,*) '*                                 *'
         WRITE(IU06,*) '*  AMOSOP: ', AMOSOP, 'RMOSOP : ',RMOSOP
         WRITE(IU06,*) '*  AMONOP: ', AMONOP, 'RMONOP : ',RMONOP
         WRITE(IU06,*) '*  AMOWEP: ', AMOWEP, 'RMOWEP : ',RMOWEP
         WRITE(IU06,*) '*  AMOEAP: ', AMOEAP, 'RMOEAP : ',RMOEAP
         WRITE(IU06,*) '*  DELLO: ', DELLO
         WRITE(IU06,*) '*                                 *'
         WRITE(IU06,*) '*     THE PROGRAM ABORTS          *'
         WRITE(IU06,*) '***********************************'
         CALL ABORT1
      ENDIF

      LLINTERPOL=.TRUE.

      IF(KAMONOP.EQ.KRMONOP.AND.KAMOSOP.EQ.KRMOSOP.AND.
     &   KAMOWEP.EQ.KRMOWEP.AND.KAMOEAP.EQ.KRMOEAP ) THEN
         IF(JRGG.EQ.KRGG.AND.
     &      NC.EQ.NGX.AND.NR.EQ.NGY) THEN

            LLINTERPOL=.FALSE.

            IF(KRGG.EQ.1) THEN
              DO J=1,NGY
                IF(RLONRGG(J).NE.KLONRGG(J)) THEN
                  LLINTERPOL=.TRUE.
                  EXIT
                ENDIF
              ENDDO
            ENDIF

         ENDIF
      ENDIF

      IF(LLINTERPOL) THEN

!       INTERPOLATE TO WAVE MODEL GRID
!       ------------------------------

        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' THE FIELDS FROM ATM. MODEL'
        IF(IREPR.EQ.0.AND.JRGG.EQ.0) THEN
          WRITE(IU06,*) ' ON A REGULAR LATITUDE/LONGITUDE GRID '
        ELSEIF(IREPR.EQ.0.AND.JRGG.EQ.1) THEN
          WRITE(IU06,*) ' ON A REDUCED LATITUDE/LONGITUDE GRID '
        ELSEIF(IREPR.EQ.4.AND.JRGG.EQ.0) THEN
          WRITE(IU06,*) ' ON A REGULAR GAUSSIAN GRID '
        ELSE
          WRITE(IU06,*) ' ON A REDUCED GAUSSIAN GRID '
        ENDIF
        WRITE(IU06,*)' ARE TO BE INTERPOLATED ONTO WAVEMODEL GRID'
        WRITE(IU06,*) ' '

        DO K=1,NGY                                                

          JSN=NGY-K+1

          XK = RMONOP - (AMOSOP + (JSN-1)*XDELLA)
          XK = XK/DELLA+1.
          KK(K) = MAX(1,INT(XK))
          KSN = NR-KK(K)+1
          KK1 = MIN(KK(K)+1,NR)
          KSN1 = NR-KK1+1
          DK1(K)=XK-REAL(KK(K))

          DO I=1,KLONRGG(JSN)

            XI = AMOWEP+(I-1)*ZDELLO(JSN) - RMOWEP
            XI = MOD(XI+720.,360.)

            XII = XI/RDELLO(KSN)+1.
            III = MAX(1-IPERIODIC,INT(XII))
            II(I,K) = MIN(III,RLONRGG(KSN))
            DII1(I,K) = XII-REAL(II(I,K))

            XIIP = XI/RDELLO(KSN1)+1.
            IIIP = MAX(1-IPERIODIC,INT(XIIP))
            IIP(I,K) = MIN(IIIP,RLONRGG(KSN1)) 
            DIIP1(I,K)=XIIP-REAL(IIP(I,K))

          ENDDO
        ENDDO

      ELSE
        WRITE(IU06,*) ' '
        WRITE(IU06,*)
     &   ' THE FIELDS FROM ATM. MODEL MATCH THE WAVEMODEL GRID '
        WRITE(IU06,*) ' '
      ENDIF

      DEALLOCATE (RDELLO)

      RETURN
      END SUBROUTINE INITIALINT

