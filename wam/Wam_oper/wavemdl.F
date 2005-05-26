      SUBROUTINE WAVEMDL (CBEGDAT, KDELWI, KDURAT,
     &                    NFIELDS, NGPTOTG, NC, NR, FCRANGE,
     &                    ISEC1, ISEC2, FIELDS,
     &                    NWVFIELDS, WVFLDG,
     &                    NLONW, NLATW, LDSTOP, LDWRRE,
     &                    LDRESTARTED, ZDELATM, KQGAUSS,
     &                    LDWCOUNORMS, MASK_IN, MASK_OUT)

!****  *WAVEMDL* - SUPERVISES EXECUTION OF MAIN MODULES                 
!****              OF THE WAVE MODEL                                    

!      LIANA ZAMBRESKY    GKSS/ECMWF    OCTOBER 1988                    

!      MODIFICATION.
!      -------------
!         H. GUNTHER   ECMWF  MARCH    1990                      
!         P. LIONELLO  ECMWF  APRIL    1990  DATA ASSIMILATION   
!                                            MODULE WAMASSI ADDED.
!         J. BIDLOT    ECMWF  FEBRUARY 1996  MESSAGE PASSING.
!         J. DOYLE     ECMWF  OCTOBER  1996  ATMOSPHERIC COUPLING  .
!         J. BIDLOT    ECMWF  FEBRUARY 1997  MESSAGE PASSING.
!         B. HANSEN    ECMWF  MARCH    1997  SIGNAL HANDLING.
!          LDSTOP* - SET .TRUE. IF STOP SIGNAL RECEIVED.
!          LDWRRE* - SET .TRUE. IF RESTART FILE SIGNAL RECEIVED.
!          KDELWI* - WIND INPUT TIME STEP WHEN COUPLED TO THE IFS.
!         J. BIDLOT    ECMWF  April    1997  ADD ZDELATM IN PARAM. LIST 
!         S. ABDALLA   ECMWF  OCTOBER  2001  INCLUSION OF AIR DENSITY & Zi/L
!                                            GENERALIZE THE INTERFACE WITH
!                                            THE ATMOSPHERIC MODEL
!         G.Mozdzynski ECMWF  January  2005  OPTIMISE COUPLING COMMS
!                                            

!     PURPOSE.                                                          
!     --------                                                          

!          THIS SUBROUTINE SUPERVISES THE EXECUTION OF                  
!          MAIN MODULES FOR WAM MODEL INITIALIZATION,                   
!          WIND FIELD PREPROCESSING, WAM MODEL EXECUTION,               
!          AND WAVE DATA ASSIMILATION.                                  

!*    INTERFACE.                                                        
!     ----------                                                        

!          SEE MAIN MODULES SUB INITMDL, PREWIND, WAMODEL, WAMASSI.     

!     METHOD.                                                           
!     -------                                                           

!          THE FIRST TIME WAVEMDL IS CALLED, THE WAM MODEL IS           
!          INITIALIZED. THIS INITIALIZATION INCLUDES GETTING            
!          FROM ECFILE THE INITIAL SEA STATE FILES, FILLING             
!          COMMON BLOCKS DEFINING THE GRID AND SETTING GENERAL          
!          PARAMETERS. IN THE FIRST AND ALL SUBSEQUENT CALLS TO         
!          WAVEMDL  PREWIND REFORMATS THE WINDS INTO THE WAM            
!          MODEL BLOCKED STRUCTURE AND THE WAM MODEL IS EXECUTED.       
!          EACH CALL TO WAMODEL INTEGRATES THE WAVE SPECTRA FORWARD     
!          IN TIME BY ONE INPUT WIND TIME STEP OR PROPAGATION TIME      
!          STEP, WHAT EVER IS GREATER.                                  

!     EXTERNALS.                                                        
!     ----------                                                        

!          INITMDL  -  INITIALIZES THE WAM MODEL.                       
!                      GETS RECOVERY FILES OUT OF ECFILE,               
!                      SETS COMMON BLOCKS NECESSARY TO DEFINE           
!                      THE GRID AND BLOCKING STRUCTURE.                 
!                      DEFINES GENERAL PARAMETERS.                      

!          PREWIND  -  REFORMATS WINDS ON THE GAUSSIAN GRID             
!                      INTO THE WAM MODEL BLOCKED STRUCTURE.            

!          WAMODEL  -  INTEGRATES THE WAVE SPECTRA FORWARD IN TIME BY   
!                      ONE WIND INPUT TIME STEP OR ONE PROPAGATION      
!                      TIME STEP, WHATEVER IS GREATER.                  

!          WAMASSI  -  SUPERVISES DATA ASSIMILATION:                    
!                      PREPROCESSES DATA; PRODUCES ANALYSED             
!                      INTEGRATED QUANTITIES BY OPTIMAL INTERPOLATION;  
!                      ANALYSES WAVE SPECTRA ; SAVE ANALYSIS FOR        
!                      OUTPUT AND NEXT TIME STEP MODEL COMPUTATION.     

!     REFERENCES.                                                       
!     -----------                                                       

!          NONE                                                         

! -------------------------------------------------------------------   

      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWCOUT  , ONLY : CASS     ,NASS
      USE YOWCOUP  , ONLY : LWCOU    ,LWCOU2W  ,LWCOUNORMS
      USE YOWSTAT  , ONLY : MARSTYPE ,CDATEA   ,CDATEE   ,CDATEF   ,
     &            CDTPRO   ,IDELPRO  ,IDELWI   ,IDELWO   ,IASSI    ,
     &            MFDBSF   ,LSARINV
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWMAP   , ONLY : ZDELLO   ,IQGAUSS
      USE YOWTEXT  , ONLY : LRESTARTED
      USE YOWSPEC, ONLY : U10OLD   ,THWOLD   ,USOLD    ,Z0OLD    ,
     &            TAUW   ,ROAIRO   ,ZIDLOLD  ,USNEW    ,Z0NEW    ,
     &            NSTART ,NEND
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWWIND  , ONLY : CDAWIFL
      USE FDBSUBS

! --------------------------------------------------------------------- 


!     COUPLING ARRAYS
      REAL WVFLDG(NLONW,NLATW,NWVFIELDS), ZDELATM(NLATW)
      REAL, DIMENSION(NWVFIELDS) :: FAVG,FMIN,FMAX

      INTEGER MASK_IN(NGPTOTG)
      INTEGER MASK_OUT(NLONW,NLATW)

      INTEGER NFIELDS, NGPTOTG, NC, NR
      REAL FCRANGE
      INTEGER ISEC1(255), ISEC2(22+NR)
      REAL FIELDS(NGPTOTG,NFIELDS)
      
      CHARACTER CBEGDAT*12, CDTASS*12
      CHARACTER*12 :: FLABEL(NWVFIELDS)

      LOGICAL FRSTIME, LDSTOP, LDWRRE, LDRESTARTED, LDWCOUNORMS

      DATA FRSTIME, NADV /.TRUE., 0 /

!     NADV  NUMBER OF PROPAGATION TIME STEPS IN ONE CALL OF WAMODEL.    
! --------------------------------------------------------------------- 

!*    1.  THE FIRST CALL TO WAVEMDL PERFORMS INITIALIZATION.            
!         --------------------------------------------------            

 1010 FORMAT ('  **************************************************')
 1011 FORMAT ('  *                                                *')
 1012 FORMAT ('  *    WAVEMDL: INITMDL STARTS FOR (RE-)START.     *')
 1013 FORMAT ('  *    ========================================    *')
 1014 FORMAT ('  *    TWO-WAY INTERACTION WIND AND WAVES          *')
 1015 FORMAT ('  *    ONE-WAY INTERACTION WIND AND WAVES          *')
 1016 FORMAT ('  * START DATE OF RUN      ', A12,  ' (CDATEA)   *')
 1018 FORMAT ('  * END DATE OF RUN        ', A12,  ' (CDATEE)   *')
 1019 FORMAT ('  * TOTAL RUN IN HOURS     ', I10,  ' (KDURAT)     *')
 1020 FORMAT ('  * INTERACTION TIME STEP  ', I10,  ' (KDELWI)     *')

      IF (FRSTIME) THEN                                                 
!       if you change the -5 value, do the same in wamodel
        MFDBSF = -5
        IF(LWCOU) THEN
          IU06=20
          IL = LEN_TRIM(CBEGDAT)
          IF(IL.NE.12) THEN
            WRITE (IU06,*) ' NON-Y2K COMPLIANT DATE IN CALLING WAVEMDL'
            WRITE (IU06,*) ' CBEGDAT =  ',CBEGDAT
            WRITE (IU06,*) ' PROGRAM WILL ABORT '
            CALL FLUSH(IU06)
            CALL ABORT1
          ENDIF
          CDATEA = CBEGDAT
          CDATEF = CDATEA
          CDATEE = CDATEA
          IDURAT = KDURAT*60
          CALL INCDATE (CDATEE,IDURAT)
          WRITE (IU06,1010)
          WRITE (IU06,1011)
          WRITE (IU06,1012)
          WRITE (IU06,1013)
          IF (LWCOU2W) THEN
            WRITE (IU06,1014)
          ELSE
            WRITE (IU06,1015)
          ENDIF
          WRITE (IU06,1011)
          WRITE (IU06,1016) CDATEA
          WRITE (IU06,1018) CDATEE
          WRITE (IU06,1019) KDURAT/60
          WRITE (IU06,1020) KDELWI
          WRITE (IU06,1011)
          WRITE (IU06,1010)

        ENDIF
        CALL INITMDL (KDELWI,NADV)
        KQGAUSS=IQGAUSS
        FRSTIME = .FALSE.                                              

        IF (ITEST.GE.1) THEN
          WRITE(IU06,*) ' SUB. WAVEMDL: INITMDL DONE'                 
        ENDIF

        IF(LWCOU) THEN
          IF(NLONW.NE.NGX.OR.NLATW.NE.NGY) THEN
            WRITE (IU06,*) ' *********************************'
            WRITE (IU06,*) ' *                               *'
            WRITE (IU06,*) ' * PROBLEM IN WAVEMDL..........  *'
            WRITE (IU06,*) ' * PROBLEM WITH NLONW AND NLATW  *'
            WRITE (IU06,*) ' * NOT EQUAL TO NGX   AND NGY  : *'
            WRITE (IU06,*) ' * ============================= *'
            WRITE (IU06,*) ' *                               *'
            WRITE (IU06,*) ' * NLONW=',NLONW
            WRITE (IU06,*) ' * NLATW=',NLATW
            WRITE (IU06,*) ' * NGX=',NGX
            WRITE (IU06,*) ' * NGY=',NGY
            WRITE (IU06,*) ' *                               *'
            WRITE (IU06,*) ' *                               *'
            WRITE (IU06,*) ' *********************************'
            CALL ABORT1
          ENDIF
          IF (LRESTARTED)  THEN
            LDRESTARTED = .TRUE.
!            RETURN
          ELSE
            LDRESTARTED = .FALSE.
          ENDIF
        ENDIF

      ELSE
!       CHECK IF THE INTERACTION TIME STEP HAS NOT CHANGED
!       OTHERWISE RESET NADV
        IF(IDELWI.NE.KDELWI .AND. LWCOU) THEN
          IDELWO = KDELWI
          IDELWI = KDELWI
          NADV = IDELWI/IDELPRO
          CDAWIFL=CDTPRO
          CALL INCDATE(CDAWIFL,IDELWI)
          WRITE (IU06,*) ' ***********************************'
          WRITE (IU06,*) ' *                                 *'
          WRITE (IU06,*) ' * IN WAVEMDL :                    *'
          WRITE (IU06,*) ' * INTERACTION TIME STEP WAS RESET *' 
          WRITE (IU06,*) ' * KDELWI = ', KDELWI 
          WRITE (IU06,*) ' *                                 *'
          WRITE (IU06,*) ' ***********************************'
          CALL FLUSH(IU06)
          IF(IDELPRO*NADV.NE.IDELWI) THEN
            WRITE (IU06,*) ' ****************************************'
            WRITE (IU06,*) ' *                                      *'
            WRITE (IU06,*) ' * PROBLEM IN WAVEMDL :                 *'
            WRITE (IU06,*) ' * THE NEW INTERACTION TIME STEP IS NOT *'
            WRITE (IU06,*) ' * A MULTIPLE OF IDELPRO  !!!           *'
            WRITE (IU06,*) ' * KDELWI = ', KDELWI 
            WRITE (IU06,*) ' * IDELPRO = ', IDELPRO 
            WRITE (IU06,*) ' *                                      *'
            WRITE (IU06,*) ' ****************************************'
            CALL ABORT1
          ENDIF
        ENDIF

!       KEEP ATMOSPHERIC MODEL INFORMED ABOUT OUR GRID
        KQGAUSS=IQGAUSS

      ENDIF

! --------------------------------------------------------------------  

!*    2.0  INTEGRATE THE WAVE SPECTRA FORWARD IN TIME.                  
!          -------------------------------------------                  


!*    2.1  REFORMAT WINDS FROM GAUSSIAN TO BLOCKED.                     
!          ----------------------------------------                     

      CALL PREWIND (U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     &              ROAIRO, ZIDLOLD,
     &              NFIELDS, NGPTOTG, NC, NR, FCRANGE,
     &              ISEC1, ISEC2, FIELDS, MASK_IN)

!     WRITE(IU06,*) ' WAVEMDL: AFTER PREWIND, NGPTOTG=', NGPTOTG,
!    & ' SUM(MASK_IN)=', SUM(MASK_IN)

      IF (ITEST.GE.1) THEN
        WRITE(IU06,*) ' SUB. WAVEMDL: PREWIND DONE'                   
      ENDIF

!*    2.2  INTEGRATE THE WAVE SPECTRA FORWARD IN TIME.                  
!          -------------------------------------------                  

      CALL SETMARSTYPE

      CALL WAMODEL (NADV, LDSTOP, LDWRRE)

      IF (ITEST.GE.1) THEN
        WRITE(IU06,*) ' SUB. WAVEMDL: WAMODEL DONE'
      ENDIF

!*    2.3  DATA ASSIMILATION AND OR SAR INVERSION.
!*         IF REQUESTED AND MODEL IS IN ANALYSIS PERIOD.
!          ---------------------------------------------

      IF (IASSI.EQ.1) THEN

        MARSTYPE = 'an'

!       UPDATE ANALYSIS TIME
        DO J=1,NASS
          IF(CDTPRO.EQ.CASS(J)) THEN
            CDTASS=CDTPRO
            EXIT
          ENDIF
        ENDDO

        IF (NASS.GT.0 ) THEN
          IF ( CDTPRO.EQ.CDTASS ) THEN
            CALL WAMASSI (LDSTOP, LDWRRE)
            IF (ITEST.GE.1) THEN
              WRITE(IU06,*) ' SUB. WAVEMDL : WAMASSI DONE'
            ENDIF
          ENDIF
        ELSEIF ((.NOT.LWCOU .AND. CDTPRO.LE.CDATEF ) 
     &      .OR. 
     &      (LWCOU .AND. CDTPRO.EQ.CDATEF)) THEN
          CALL WAMASSI (LDSTOP, LDWRRE)
          IF (ITEST.GE.1) THEN
            WRITE(IU06,*) ' SUB. WAVEMDL: WAMASSI DONE'
          ENDIF
        ENDIF
      ELSE IF (IASSI.NE.1 .AND. LSARINV) THEN

!!!!
!!!!! we may have to introduce a list for sar inversion time !!!!
!       UPDATE ANALYSIS TIME (which is now the inversion time)
        DO J=1,NASS
          IF(CDTPRO.EQ.CASS(J)) THEN
            CDTASS=CDTPRO
            EXIT
          ENDIF
        ENDDO

        IF (NASS.GT.0 ) THEN
          IF ( CDTPRO.EQ.CDTASS ) THEN
            CALL SARINVERT
            IF(ITEST.GE.1)
     &         WRITE(IU06,*) '   SUB. WAVEMDL: SARINVERT DONE'
          ENDIF
        ELSEIF ((.NOT.LWCOU .AND. CDTPRO.LE.CDATEF ) 
     &      .OR. 
     &      (LWCOU .AND. CDTPRO.EQ.CDATEF)) THEN
          CALL SARINVERT
          IF(ITEST.GE.1) WRITE(IU06,*) '   SUB. WAVEMDL: SARINVERT done'
        ENDIF
      ENDIF

! PUT ZDELLO IN ZDELATM FOR THE ATMOSPHERIC COUPLING
      IF(LWCOU) THEN
        DO K=1,NLATW
          ZDELATM(K) = ZDELLO(K)
        ENDDO

! Fields to be passed to the atmospheric model are:
!  1. Charnock field
        FLABEL(1)='Charnock'
        DO IG=1,IGL
          DO IJ=IJS(IG),IJL(IG) 
            USNEW(IJ) = USOLD(IJ,IG)
            Z0NEW(IJ) = Z0OLD(IJ,IG)
          ENDDO
          CALL OUTBETA (IJS(IG), IJL(IG), IG, IGL,USNEW,Z0NEW,
     &                    NSTART, NEND, WVFLDG(1,1,1),
     &                    NLONW, NLATW)
        ENDDO
!  2. Test field (Hs).
!RMIT!  IF (NWVFIELDS.GE.2) THEN 
!RMIT!!  WVFLDG(:,:,2)=WHGTTG             ! Test ONLY !!!!!!!!!!!!!!!!!!!!!!!!!
!RMIT!!  WVFLDG(:,:,2)=10.*WVFLDG(:,:,1)  ! Test ONLY !!!!!!!!!!!!!!!!!!!!!!!!!
!RMIT!   DO IG=1,IGL
!RMIT!    CALL OUTTHEOTHER (IJS(IG), IJL(IG), IG, IGL,USNEW,Z0NEW,
!RMIT&                    NSTART, NEND, WVFLDG(1,1,2),
!RMIT&                    NLONW, NLATW)
!RMIT!   ENDDO 
!RMIT!  ENDIF
! Print TEST values:
! TEST ONLY !!!!!!!!!!!!!!!!!!!!!!!!!
!RMIT!  IF (NWVFIELDS.GE.2) THEN
!RMIT!   WRITE(6,'(A)') ' Sent second field (10.*Charnock):'
!RMIT!!  DO JLA=51,120
!RMIT!!   WRITE(6,'(60F6.3)') (WVFLDG(JLO,JLA,2),JLO=201,250)
!RMIT!!  ENDDO
!RMIT!!  CALL FLUSH(6)
!RMIT!   WRITE(6,'(A)') ' Summary Characteristics (W2):'
!RMIT!   IF (ANY(WVFLDG(:,:,2).GT.0.0)) THEN
!RMIT!    WRITE(6,'(A,G23.16)') '  * Abs. Minimum Value W2: ', 
!RMIT&      MINVAL(WVFLDG(:,:,2))
!RMIT!    WRITE(6,'(A,G23.16)') '  * Minimum Value W2: ', 
!RMIT&      MINVAL(WVFLDG(:,:,2),MASK=WVFLDG(:,:,2).GT.0.0)
!RMIT!    WRITE(6,'(A,G23.16)') '  * Maximum Value W2: ', 
!RMIT&      MAXVAL(WVFLDG(:,:,2),MASK=WVFLDG(:,:,2).GT.0.0)
!RMIT!    WRITE(6,'(A,G23.16)') '  * Mean Value W2: ', 
!RMIT&      SUM(WVFLDG(:,:,2),MASK=WVFLDG(:,:,2).GT.0.0) / 
!RMIT&      COUNT(WVFLDG(:,:,2).GT.0.0)
!RMIT!   ELSE
!RMIT!   WRITE(6,'(A)') ' WARNING: No values W2 were found GT 0.0!'
!RMIT!   ENDIF
!RMIT!   CALL FLUSH(6)
!RMIT!  ENDIF
!------------------------------------------------------

!       COMPUTATION OF THE NORMS OF OUTPUT FIELDS
        IF( LWCOUNORMS )THEN
          WRITE(IU06,*) ' '
          WRITE(IU06,*) ' NORM OF FIELDS RETURNED TO THE ATM. MODEL : '
          NTOT=NLATW*NLONW
          DO IFLD=1,NWVFIELDS
            FAVG(IFLD)=0.
            FMIN(IFLD)=999999999999.
            FMAX(IFLD)=-99999999999.
            DO J=1,NLATW
              DO I=1,NLONW
                FAVG(IFLD)=FAVG(IFLD)+WVFLDG(I,J,IFLD)
                FMIN(IFLD)=MIN(FMIN(IFLD),WVFLDG(I,J,IFLD))
                FMAX(IFLD)=MAX(FMAX(IFLD),WVFLDG(I,J,IFLD))
              ENDDO
            ENDDO
            FAVG(IFLD)=FAVG(IFLD)/NTOT
            WRITE(IU06,*) FLABEL(IFLD),
     *                    FAVG(IFLD),FMIN(IFLD),FMAX(IFLD),NTOT
          ENDDO
          WRITE(IU06,*) ' '
          IF(ITEST.GT.0) CALL FLUSH(IU06)
        ENDIF
!       TELL ATMOS MODEL THE WAM REQUIREMENT FOR GLOBAL NORMS
        LDWCOUNORMS=LWCOUNORMS
      ENDIF

      IF ( MFDBSF .GT. 0 .AND. CDATEE .EQ. CDTPRO) THEN
        CALL GSTATS(1787,0)
        ISTAT = ICLOSEFDBSUBS (MFDBSF)
        CALL GSTATS(1787,1)
        IF (ITEST.GE.1) THEN
          WRITE(IU06,*) ' SUB. WAVEMDL: END OF RUN: CLOSE FDB'
          CALL FLUSH(IU06)
        ENDIF
      ENDIF

      RETURN                                                            
      END SUBROUTINE WAVEMDL
