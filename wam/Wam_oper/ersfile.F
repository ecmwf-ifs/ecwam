      SUBROUTINE ERSFILE (IU91, IU92, IG, IGL)                          
                                                                        
! -------------------------------------------------------------------   

!**** *ERSFILE* - FILE MANAGEMENT OF ERS FILES.                      

!     H. GUNTHER               ECMWF         JULY 1991                  

!     PURPOSE.                                                          
!     --------                                                          

!       CONTROL INPUT AND OUTPUT FILES OF ERS SPECTRAL DATA.          


!*** INTERFACE.                                                         
!    ----------                                                         
!      *CALL* *ERSFILE (IU91, IU92 IG, IGL)*                            
!         *IU91*   - INPUT UNIT OF LCATIOB FILE.                        
!         *IU92*   - OUTPUT UNIT FOR SPECTRA.                           
!         *IG*     - BLOCK NUMBER.                                      
!         *IGL*    - NUMBER OF BLOCKS.                                  

!     EXTERNALS.                                                        
!     ----------                                                        

!        NONE.                                                          


!     METHOD.                                                           
!     -------                                                           

!       NONE.                                                           

!     REFERENCES.                                                       
!     -----------                                                       

!        NONE.                                                          

!     MODIFICATIONS.                                                    
!     --------------                                                    

!        B. HANSEN    *ECMWF*      FEB 93'                              
!           ACCESS TO ECFILE DELETED. INPUT FILES COLYYMMDDhhmm HAVE    
!           TO BE ONLINE ALREADY OTHERWISE NO OUTPUT OF WAVE SPECTRA    
!           CLOSE TO POSITIONS OF ERS1-SAR-SPECTRA.                     
!           WAVE SPECTRA WRITTEN TO UNIT 71. FILES ARE OPEND FOLLOWING  
!           THE KNOWN NAMING KONVENTION COSYYMMDDhhmm. THESE FILES      
!           ARE CLOSED (STATUS=KEEP) AND ARE TO BE SAVED AFTER THE      
!           MODEL RUN.                                                  

!     J. BIDLOT       *ECMWF*      JUNE 1996 MIGRATION TO FUJITSU

!---------------------------------------------------------------------- 

      USE YOWCOER  , ONLY : NERS     ,CDTERS   ,IDELERS  ,IERS     ,
     &            IJERS    ,IGERS
      USE YOWSTAT  , ONLY : CDATEF   ,CDTPRO   ,CDATEA   ,MARSTYPE
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK

! ----------------------------------------------------------------------

      INTEGER, ALLOCATABLE :: IDUM(:)
      CHARACTER * 2 CLH       ! HELPER STRING TO STORE HOUR.            
      CHARACTER *15 FNAME                                               
      CHARACTER *12 CL_DATE
      CHARACTER *14 CL_TIME
      LOGICAL FRSTIME                                                   
      DATA FRSTIME / .TRUE./                                            

! ----------------------------------------------------------------------

!*    1.  SAVE SPECTRA FILE.                                            
!         ------------------                                            

      IF (IG.EQ.-9 .AND.IGL.EQ.-9) THEN

        IF(ALLOCATED(IJERS)) DEALLOCATE(IJERS)
        IF(ALLOCATED(IGERS)) DEALLOCATE(IGERS)

!       DO NOT DEAL WITH IU91 AND IU92 IF NOT ON PE 1
        IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
          CLOSE (IU92,STATUS='KEEP')
          IF(ITEST.GT.0) THEN
            WRITE (IU06,*) '      SUB. ERSFILE: UNIT ',
     &       IU92, '  CLOSED    '
            WRITE (IU06,*) '                                      '
            CALL FLUSH (IU06)
          ENDIF
        ENDIF
        IERS = 0
        RETURN
      ENDIF

! ----------------------------------------------------------------------

!*    2. INITIAL DATES IF FIRST CALL.                                   
!        ----------------------------                                   

      IF (FRSTIME) THEN                                                 
        NERS = 1000                                                       
        IERS = 0                                                       
        IDELERS = 21600                                                
        CDTERS  = CDATEA                                               
        CLH=CDATEA(9:10)
        IF (CLH.EQ.'09'.OR.CLH.EQ.'15'.OR.
     &   CLH.EQ.'21'.OR.CLH.EQ.'03') THEN
          CALL INCDATE (CDTERS, -10800)
        ENDIF
        FRSTIME = .FALSE.                                              
        IF (ITEST.GE.3) THEN
          WRITE(IU06,*) " SUB: ERSFILE (for the first time) "//
     &     " CDATEA=", CDATEA, " CDTERS=", CDTERS, " CDATEF=", CDATEF
          CALL FLUSH (IU06)
        ENDIF
      ENDIF                                                             

      DO WHILE (CDTERS.LT.CDTPRO)
        CALL INCDATE (CDTERS, IDELERS)                                  
      ENDDO

! ----------------------------------------------------------------------

!*    3. FETCH LOCATION FILE.                                           
!        --------------------                                           

!     DO NOT DEAL WITH IU91 AND IU92 IF NOT ON PE 1
      IF(LMESSPASS.AND.(IRANK.NE.1)) RETURN

      IERS = 0                                                          
      IF (CDTERS.LE.CDATEF.OR.MARSTYPE.EQ.'an'.OR.MARSTYPE.EQ.'fg') THEN
        WRITE(FNAME,'(''COL'',A12)') CDTERS                            
        OPEN (UNIT=IU91, FILE=FNAME, FORM='UNFORMATTED', STATUS='OLD', 
     &   IOSTAT=IERR)                                             
        IF (IERR.NE.0) THEN                                            
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            IF(ITEST.GT.0) THEN
              WRITE (IU06,*) '                                     '
              WRITE (IU06,*) '      SUB. ERSFILE: FILE ',
     &         FNAME, ' NOT OPENED'
              WRITE (IU06,*) '                                     ' 
            ENDIF
          ENDIF
          RETURN                                                      
        ELSE                                                           
          READ(IU91)CL_DATE
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            IF(ITEST.GT.0) THEN
              WRITE (IU06,*) '                                     '
              WRITE (IU06,*) '      SUB. ERSFILE: FILE ',
     &         FNAME, ' OPENED    '
            ENDIF
          ENDIF

          WRITE(FNAME,'(''COS'',A12)') CDTERS                         
          OPEN (UNIT=IU92, FILE=FNAME(1:LEN_TRIM(FNAME))//MARSTYPE,   
     &     FORM='UNFORMATTED', STATUS='UNKNOWN')                 
          IF(ITEST.GT.0) THEN
            WRITE (IU06,*) '      SUB. ERSFILE: FILE ',
     &       FNAME, ' OPENED'
            CALL FLUSH (IU06)
          ENDIF
          WRITE(IU92) CDTPRO                                    
        ENDIF

        ALLOCATE(IJERS(NERS))
        ALLOCATE(IGERS(NERS))

 3001   CONTINUE                                                       
        IERS = IERS + 1                                                

        IF (IERS.GT.NERS) THEN                                         
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++'   
            WRITE (IU06,*) ' +                                     +'   
            WRITE (IU06,*) ' +        WARNING SUB ERSFILE.         +'   
            WRITE (IU06,*) ' +     ==========================      +'   
            WRITE (IU06,*) ' + NUMBER OF OUTPUT POINTS             +'   
            WRITE (IU06,*) ' + EXCEEDS DIMENSION NERS = ', NERS         
            WRITE (IU06,*) ' + THE DIMENSION WILL BE DOUBLED       +'   
            WRITE (IU06,*) ' + IN ORDER TO CONTNUE                 +'   
            WRITE (IU06,*) ' +                                     +'   
            WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++'   
          ENDIF
          NERSOLD=NERS
          NERS=NERS*2
          ALLOCATE(IDUM(NERSOLD))
          IDUM=IJERS
          DEALLOCATE(IJERS)
          ALLOCATE(IJERS(NERS))
          DO JC=1,NERSOLD
            IJERS(JC)=IDUM(JC)
          ENDDO

          IDUM=IGERS
          DEALLOCATE(IGERS)
          ALLOCATE(IGERS(NERS))
          DO JC=1,NERSOLD
            IGERS(JC)=IDUM(JC)
          ENDDO
        
          DEALLOCATE(IDUM)  

        ENDIF                                                          

        READ (IU91, END=3002, ERR=3002, IOSTAT=IOS) CL_TIME,
     &   ZLATIS, ZLONGS, ZLATIW, ZLONGW, IGERS(IERS), IJERS(IERS)
        IF (IGERS(IERS).EQ.0) THEN                                     
          IERS = IERS-1                                               
          GOTO 3001                                                   
        ENDIF                                                          
        IF (IERS.GT.1) THEN                                            
          DO I=1,IERS-1                                          
            IF (IGERS(I).EQ.IGERS(IERS) .AND.                        
     &       IJERS(I).EQ.IJERS(IERS)) THEN                        
              IERS = IERS-1                                         
              GOTO 3001                                             
            ENDIF                                                    
          ENDDO
        ENDIF                                                          
        GOTO 3001                                                      
 3002   CONTINUE                                                       
        IERS = IERS - 1                                                
        CLOSE (IU91, STATUS='KEEP')                                    
        IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
          IF(ITEST.GT.0) THEN
            WRITE(IU06,*) '      SUB. ERSFILE: NEW OUTPUT',
     &       ' POINTS FOR SATELLITE COLOCATION.  IERS = ', IERS,
     &       ' LAST READ STATUS = ', IOS
            CALL FLUSH (IU06)
          ENDIF
        ENDIF
      ENDIF                                                             

      RETURN                                                            
      END SUBROUTINE ERSFILE
