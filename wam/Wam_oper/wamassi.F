      SUBROUTINE WAMASSI(BETAG,NLONW,NLATW,LDSTOP,LDWRRE)

! ----------------------------------------------------------------------

!****  *WAMASSI* - SUPERVISES EXECUTION OF MAIN MODULES                 
!****              OF THE WAVE DATA ASSIMILATION.                       

!      P. LIONELLO     ECMWF     APRIL   1990
!      J. BIDLOT       ECMWF     FEBRARY 1996 MESSAGE PASSING.
!      J. BIDLOT       ECMWF     FEBRARY 1997 MODULE and ATMOSPHERIC
!                                COUPLING.
!      J. BIDLOT       ECMWF     MARCH   1997  ADD SAVSTRESS AND SAVSPEC
!      B. HANSEN       ECMWF     APRIL   1997  SIGNAL HANDLING.
!      B. HANSEN       ECMWF     JANUARY 1998  1.1 OUTSOURCED.

!     PURPOSE.                                                          
!     --------                                                          

!         TO ANALYSE THE WAVE FIELD SUBSTITUTING THE FIRST GUESS        
!         SPECTRA WITH ANALYSED SPECTRA AND FIRST GUESS VALUES          
!         WITH ANALYSED VALUES IN THE GLOBAL GRID                       

!*    INTERFACE.                                                        
!     ----------                                                        

!     *CALL* *WAMASSI(BETAG,NLONW,NLATW,LDSTOP,LDWRRE)
!        *BETAG*     REAL      CHARNOCK PARAMETER TO BE PASSED BACK TO
!                              THE ATMOSPHERIC MODEL
!        *NLONW*     INTEGER   FIRST DIMENSION OF BETAG
!        *NLATW*     INTEGER   SECOND DIMENSION OF BETAG
!        *LDSTOP*    LOGICAL   SET .TRUE. IF STOP SIGNAL RECEIVED.
!        *LDWRRE*    LOGICAL   SET .TRUE. IF RESTART SIGNAL RECEIVED.

!     METHOD.                                                           
!     -------                                                           

!     EXTERNALS.                                                        
!     ----------                                                        

!        *ANALYSE*  - UTILITY OF OIFIELD . IT PRODUCES THE OPTIMALLY    
!                     INTERPOLATED VALUE AT A SINGLE GRIDPOINT.         
!        *FDUR*     - COMPUTES WINDSEA DURATION.                        
!        *FEMEAN*   - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT. 
!        *FUSTAR*   - COMPUTES ANALYSED USTAR ACCORDING TO ANALYSED     
!                     WINDSEA ENERGY.                                   
!        *FWSEA*    - DETERMINES IF WINDSEA IS PRESENT IN THE SPECTRUM,,
!                     COMPUTING EVENTUALLY WINDSEA ENERGY AND MEAN FREQ.
!        *GRFIELD*  - READS GRIDDED MEASUREMENTS AND DISTRIBUTES THEM.
!        *GRDATA*   - TRANSFERS MEASUREMENTS TO THE MODEL GRID.         
!        *GSFILE*   - ROUTINE TO DYNAMICALLY FETCH OR DISPOSE FILES.    
!        *INCDATE*  - UPDATE DATE TIME GROUP.                           
!        *MAKEBLO*  - MAKES BLOCKS FROM GRIDDED FIEDS.                  
!        *MAKEGRID* - MAKES GRIDDED FIELDS.                             
!        *OIFIELD*  -  MERGES MODEL AND MEASUREMENTS BY O.I.            
!        *OUTAN*    - TRANSFERS INTEGRATED ANALYSED QUANTITIES TO FILE. 
!        *PEAKFR*   - COMUTE PEAK FREQUENCY.                            
!        *SAVSTRESS*- DISPOSE STRESS/WIND RESTART FILES.
!        *SAVSPEC   - DISPOSE SPECTRUM RESTART FILES.
!        *SEMEAN*   - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.   
!        *SERMAT*   - UTILITY OF OIFIELD TO SUM TWO MATRICES.           
!        *SPLITBL*  - STORING OVERLAPPING LATITUDES OF BLOCKS.          
!        *STHQ*     - COMPUTATION OF MEAN WAVE DIRECTION AT EACH        
!                     GRID POINT.                                       
!        *UPDATE*   - ANALYSES THE WAM MODEL SPECTRA.                   
!        *UPWSPEC*  - PRODUCES THE UPDATED SPECTRUM.                    
!        *WSMFEN*   - UTILITY OF FWSEA , COMPUTING  WINDSEA ENERGY AND  
!                     DURATION IN EACH GRID POINT.                      
!        *OUTWNORM* - COMPUTES A FEW NORMS OF GRIDDED FIELDS

!     REFERENCES.                                                       
!     -----------                                                       

!          NONE                                                         

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWCOUT  , ONLY : FFLAG20  ,FFLAG21  ,FFLAG25  ,FFLAG26  ,
     &            GFLAG20  ,LFDB     ,IPFGTBL
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRIBHD, ONLY : NKSEK0   ,NKSEK1   ,NKSEK2   ,NKSEK3   ,
     &            NKSEK4   ,NPSEK2   ,NPSEK3   ,KSEK0    ,KSEK1    ,
     &            KSEK2    ,KSEK3    ,KSEK4    ,PSEK2    ,PSEK3    ,
     &            PPMISS   ,PPEPS    ,PPREC    ,NTENCODE
      USE YOWGRID  , ONLY : NLONRGG  ,IGL      ,IJS      ,IJL2     ,IJL
      USE YOWINTP  , ONLY : WHGTTG   ,WDIRTG   ,WPKFTG   ,WMNFTG   ,
     &            USTARG   ,UDIRG    ,TAUWG    ,CDG      ,SMEANG   ,
     &            U10G
      USE YOWINTS  , ONLY : WHGTSG   ,WDIRSG   ,WMNFSG   ,WHGTWG   ,
     &            WDIRWG   ,WMNFWG
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,IRGG     ,AMOWEP   ,
     &            AMOSOP   ,AMOEAP   ,AMONOP   ,XDELLA   ,XDELLO
      USE YOWMESPAS, ONLY : LMESSPASS,LFDBIOOUT,LGRIBOUT
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP     ,KTAG
      USE YOWMLTSK , ONLY : NCHUNK
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      ,
     &            NBLO     ,NIBLO    ,CLDOMAIN
      USE YOWPCONS , ONLY : ZMISS    ,DEG
      USE YOWSPEC, ONLY : NSTART   ,NEND     ,U10NEW   ,U10OLD   ,
     &            THWNEW   ,THWOLD   ,USNEW    ,USOLD    ,Z0NEW    ,
     &            Z0OLD    ,TAUW     ,FL1      ,FL2      ,FL3
      USE YOWSTAT  , ONLY : CDATEE   ,CDATEF   ,CDTPRO   ,CDTRES   ,
     &            CDATER   ,CDATES   ,CDTINTT  ,CDTINTS  ,CDTSPT   ,
     &            CDTSPS   ,CFDBSF   ,CFDB2DSP ,MFDBSF   ,IDELWI   ,
     &            IREST    ,IASSI    ,NTASKS   ,NSIZE    ,NENSFNB  ,
     &            NTOTENS  ,CDATEA   ,MARSTYPE ,YCLASS   ,YEXPVER
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH    ,CWI
      USE YOWUNIT  , ONLY : IU11     ,IU12     ,IU13     ,IU14     ,
     &            IU22     ,IU23     ,IU27     ,IU28     ,IU04     ,IU32
      USE YOWWAMI  , ONLY : CBPLTDT  ,CEPLTDT  ,IANALPD  ,IFOREPD  ,
     &            IDELWIN  ,NFCST    ,ISTAT

! ----------------------------------------------------------------------

!*    *ADDITIONAL PARAMETERS FOR MULTITASKING

      INTEGER ITSKCA(9,NTASKS)
      EXTERNAL MT_ALLOT

! ----------------------------------------------------------------------
!     DEFINE ALLOCATABLE ARRAYS

      REAL,ALLOCATABLE,DIMENSION(:,:) :: USOI, USME, WHOI, WHME
      REAL,ALLOCATABLE,DIMENSION(:) :: WHOIB, USOIB, ETOIB 


!*    VARIABLE.     TYPE.        PURPOSE.                               
!     ---------     -----        --------                               
!      *USOI*        REAL        USTAR FIELD FROM O.I.                  
!      *USME*        REAL        USTAR FIELD FROM MEASUREMENTS          
!      *WHOI*        REAL        WAVE HEIGHT FIELD FROM O.I.            
!      *WHME*        REAL        WAVE HEIGHT FIELD FROM MEASUREMENTS    
!      *WHOIB*       REAL        WAVE HEIGHT FROM O.I. (BLOCKED)        
!      *USOIB*       REAL        USTAR FROM O.I. (BLOCKED)              
!      *ETOIB*       REAL        TOTAL ENERGY FROM O.I.                 

!-----------------------------------------------------------------------

      NAMELIST /OIP/ LMAX, SIGOBS, SIGMOD

!*    *NAMELIST* *OIP* PARAMETERS FOR THE OPTIMUM INTERPOLARTION.

      LOGICAL LOIP

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *LMAX*      INTEGER   MAXIMUM NUMBER OF GRID STEPS TO SPREAD
!                            MEASURED HS AND USTAR
!      *SIGOBS*    REAL      MEASUREMENT SCATTER.
!      *SIGMOD*    REAL      MODEL SCATTER.
!      *LOIP*      LOGICAL   TRUE IF INQUIRE ON FILE 'WTR' SUCCESSFULL.

! ----------------------------------------------------------------------
      REAL BETAG(NLONW,NLATW)
! ----------------------------------------------------------------------

      LOGICAL LDSTOP, LDWRRE, LLFLUSH
      WRITE (IU06,*) '  '                                               
      WRITE (IU06,*)                                                    
     &   ' START OF DATA ASSIMILATION: DATE IS CDTPRO: ', CDTPRO        

      LLFLUSH = .TRUE.
      ALLOCATE(WHME(NGX,NGY), USME(NGX,NGY))
      ALLOCATE(WHOI(NGX,NGY), USOI(NGX,NGY))
      ALLOCATE(WHOIB(NINF:NSUP),USOIB(NINF:NSUP),ETOIB(NINF:NSUP))

!    0. COLLECT GRIDDED FIELDS USED BY THE ASSIMILATION
!       -----------------------------------------------

      IF(LMESSPASS) THEN
        KTAG=KTAG+1
        CALL MPGATHEROIFLD(IPFGTBL(1),KTAG,WHGTTG)
        KTAG=KTAG+1
        CALL MPGATHEROIFLD(IPFGTBL(4),KTAG,USTARG)
        KTAG=KTAG+1
        CALL MPGATHEROIFLD(IPFGTBL(7),KTAG,CDG)
      ENDIF


!*    1. MEASUREMENT ARE TRANSFERRED TO THE MODEL GRID.                 
!        ----------------------------------------------                 

      CALL GRFIELD (WHME, USME, NSTART, NEND)
                                                                        
      IF (ITEST.GE.2)                                                   
     &   WRITE(IU06,*) '   SUB. WAMASSI: MEASUREMENTS ---> GRID'        

!*    2. MEASUREMENTS AND MODEL ARE MERGED BY OPTIMUM INTERPOLATION.    
!        -----------------------------------------------------------    

!*    2.1 GET PARAMETERS FROM NAMELIST OIP IF PRESENT.
!         --------------------------------------------
      LMAX=11
      IF(XDELLA.LT.0.5) LMAX=22
      SIGOBS=0.5
      SIGMOD=0.5

      LOIP=.FALSE.
      INQUIRE(FILE='oip',EXIST=LOIP)
      IF (LOIP) THEN 

        OPEN(UNIT=92,FILE='oip')
        READ (92,OIP)
        CLOSE(UNIT=92)

        IF ( LMAX .GT. 22 .OR. SIGOBS .GT. 1.0 .OR. SIGMOD .GT. 1.0
     &       .OR. LMAX .LT. 1 .OR. SIGOBS .LT. 0.1 .OR. SIGMOD .LT.
     &       0.1 ) THEN
          WRITE(IU06,*) '*******************************************'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*    FATAL ERROR IN SUB. WAMASSI          *'
          WRITE(IU06,*) '*    ===========================          *'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '* CONSTRAINTS:                            *'
          WRITE(IU06,*) '* 0 < LMAX < 23        CURRENT SETTING: ',
     &                     LMAX,'*'
          WRITE(IU06,*) '* 0.1 < SIGOBS <= 1.0  CURRENT SETTING: ',
     &                     SIGOBS,'*'
          WRITE(IU06,*) '* 0.1 < SIGMOD <= 1.0  CURRENT SETTING: ',
     &                     SIGMOD,'*'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.       *'
          WRITE(IU06,*) '* ---------------   --------------        *'
          WRITE(IU06,*) '*******************************************'
          CALL ABORT1
        ENDIF
      ENDIF
      NDIM1 = 2*LMAX+1
      NDIM2 = NDIM1*NDIM1
      IF(ITEST.GT.0) THEN
        WRITE(IU06,*) ' PARAMETER FOR OI: '
        WRITE(IU06,OIP)
        WRITE(IU06,*)
      ENDIF

!*    2.2 OI.
!         ---

      CALL OIFIELD (USME, USTARG, USOI, SIGOBS, SIGMOD, 
     &              LMAX, NDIM2)

      IF (ITEST.GE.2)
     &   WRITE(IU06,*) '   SUB. WAMASSI: OIFIELD DONE FOR WINDS'

      CALL OIFIELD (WHME, WHGTTG, WHOI, SIGOBS, SIGMOD, 
     &              LMAX, NDIM2)

      IF (ITEST.GE.2)
     &   WRITE(IU06,*) '   SUB. WAMASSI: OIFIELD DONE FOR WAVE HEIGHTS'


!*    3. ANALYSING THE SPECTRA.                                         
!        ----------------------                                         


!*    3.1 UNITS ARE PREPARED.                                           
!         -------------------                                           

      IF(.NOT.LMESSPASS) REWIND IU12
      IF (IGL.NE.1) THEN                                                
        REWIND IU11                                                    
        REWIND IU12                                                    
        REWIND IU13                                                    
        REWIND IU14                                                    

!*    3.2 CHANGE WORKING FILE UNITS.                                    
!         --------------------------                                    

        IFIL=IU11                                                      
        IU11=IU12                                                      
        IU12=IFIL                                                      
        IFIL=IU13                                                      
        IU13=IU14                                                      
        IU14=IFIL                                                      

!*    3.3 THE FIRST BLOCK IS READ IN.                                   
!         ---------------------------                                   

        READ(IU11) FL1             
        IF (ITEST.GE.2)                                                
     &   WRITE(IU06,*) '   SUB. WAMASSI: START OF READING',         
     &   ' FIRST BLOCK OF SPECTRA'                     
      ENDIF

!*    3.4 LOOP OVER THE BLOCKS.                                         
!         ---------------------                                         

      BLOCK : DO IG=1,IGL

!*    3.4.1 TRANSFERS RESULTS OF OPTIMAL INT TO BLOCKS.                 
!           -------------------------------------------                 

        CALL MAKEBLO (USOIB(IJS(IG)), USOI, IJS(IG), IJL(IG), IG)
        CALL MAKEBLO (WHOIB(IJS(IG)), WHOI, IJS(IG), IJL(IG), IG)

        DO IJ=IJS(IG),IJL(IG)                                     
          USNEW(IJ) = USOLD(IJ,IG)                                    
          THWNEW(IJ)   = THWOLD(IJ,IG)                                
        ENDDO
        DO IJ=IJS(IG),IJL(IG)                                     
          ETOIB(IJ) = -99.                                            
          IF (WHOIB(IJ).GT.0.) ETOIB(IJ) = WHOIB(IJ)**2/16.           
        ENDDO

        IF (ITEST.GE.2) THEN
          IF (IG.LE.ITESTB) WRITE(IU06,*)                             
     &     '   SUB. WAMASSI: FIELDS BLOCKED FOR IG = ', IG             
        ENDIF

!*    3.4.2 CHECK UNITS.                                                
!           ------------                                                

        IF (MOD(IG,2).EQ.1) THEN                                       

!*    3.4.2.1  ODD  BLOCKS.                                             
!              ------------                                             
          IF (LMESSPASS) THEN

            CALL UPDATE (FL1, ETOIB, USOIB, USNEW, THWNEW,
     &                   IJS(IG), IJL(IG))

            IF (ITEST.GE.2) THEN
              IF (IG.LE.ITESTB) WRITE(IU06,*)                          
     &         '   SUB. WAMASSI: SPECTRA UPDATED '   
            ENDIF

          ELSE

!*             START READING FOR NEXT BLOCK.                           

            IF (IG.NE.IGL)                                              
     &         READ(IU11) FL2       

!              UPDATE SPECTRA AND USTAR.                                
!              (PARALLEL PROCESSING - HENCE MT_ALLOT).

            NCHUNK = 0
            NLENGTH = (IJL(IG)-IJS(IG)+1)
            NWORK   = NLENGTH/NSIZE
            IF (MOD(NLENGTH,NSIZE).NE.0) NWORK = NWORK + 1

            NTSKCAP = 3
            NSTKSZ  = 1
            NSTKINC = 1
            DO JROC=1,NTASKS-1
              ITSKCA(1,JROC) = NTSKCAP
              ITSKCA(3,JROC) = JROC
              ITSKCA(4,JROC) = 1 !  JUST TO GET PAST COMP 'STACKSZW'
              ITSKCA(5,JROC) = NSTKSZ
              ITSKCA(6,JROC) = 1 !  JUST TO GET PAST COMP 'STACKINW'
              ITSKCA(7,JROC) = NSTKINC

              CALL TSKSTART(ITSKCA(1,JROC), MT_ALLOT, FL1, ETOIB,
     &         USOIB, USNEW, THWNEW, IJS(IG), IJL(IG),
     &         NSIZE,NWORK)
            ENDDO

            CALL MT_ALLOT (FL1, ETOIB, USOIB, USNEW, THWNEW,
     &       IJS(IG), IJL(IG), NSIZE,NWORK)

            DO JROC=1,NTASKS-1
              CALL TSKWAIT(ITSKCA(1,JROC))
            ENDDO


            IF (ITEST.GE.2) THEN
              IF (IG.LE.ITESTB) WRITE(IU06,*)                          
     &         '   SUB. WAMASSI: SPECTRA UPDATED FOR BLOCK IG = ', IG   
            ENDIF

!              OUTPUT OF OVERLAPPING LATITUDES AND                      
!              START OF OUTPUT OF ALL SPECTRA IN BLOCK.                 

            IF (IGL.NE.1) THEN                                          
              IF (IG.NE.1) CALL SPLITBL(FL1, IJS(IG), IJL2(IG), IU14)  
              WRITE(IU12) FL1      
              IF (ITEST.GE.2) THEN
                IF (ITESTB.GE.IG) WRITE(IU06,*)                       
     &           '   SUB. WAMASSI: OVERLAPPING LATITUDES',          
     &           ' SAVED AND START OF BLOCK OUTPUT.'                
              ENDIF
            ENDIF
          ENDIF

!           THE WIND FIELD CURRENTLY USED IS UPDATED.                   

          DO IJ=IJS(IG),IJL(IG)                                     
            U10OLD(IJ,IG) = USNEW(IJ)*U10OLD(IJ,IG)/USOLD(IJ,IG)        
            U10NEW(IJ) = U10OLD(IJ,IG)                                  
          ENDDO
          IF (ITEST.GE.2) THEN
            IF (IG.LE.ITESTB) WRITE(IU06,*)                             
     &      '   SUB. WAMASSI: WIND REPLACED FOR BLOCK IG = ', IG        
          ENDIF

!              MODEL OUTPUT INTEGRATED DATA ARE SAVED IN COMMON.        
!              SPECTRA ARE SAVED ON DISK AND PRINTED.                   

          IF (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO .OR.           
     &     CDTSPT.EQ.CDTPRO .OR. CDTSPS.EQ.CDTPRO ) THEN
            CALL OUTBS (FL1, FL3, IJS(IG), IJL(IG), IG, IGL,         
     &       IU27, IU28,THWOLD,USOLD,U10NEW,THWNEW,USNEW,TAUW,
     &       NSTART,NEND)
            IF (ITEST.GE.2) THEN
              IF (ITEST.GE.2) WRITE(IU06,*)                         
     &         '   SUB. WAMASSI: MODEL OUTPUT DONE IN SUB OUTBS'
            ENDIF
          ENDIF
        ELSE                                                           

!*    3.4.2.2  EVEN BLOCKS.                                             
!              ------------                                             

!              START READING FOR NEXT BLOCK.                           

          IF (IG.NE.IGL)                                              
     &     READ(IU11) FL1       

!              UPDATE SPECTRA AND USTAR.                                

!              (PARALLEL PROCESSING - HENCE MT_ALLOT).


          NCHUNK = 0
          NLENGTH = (IJL(IG)-IJS(IG)+1)
          NWORK   = NLENGTH/NSIZE
          IF (MOD(NLENGTH,NSIZE).NE.0) NWORK = NWORK + 1

          NTSKCAP = 3
          NSTKSZ  = 1
          NSTKINC = 1
          DO JROC=1,NTASKS-1
            ITSKCA(1,JROC) = NTSKCAP
            ITSKCA(3,JROC) = JROC
            ITSKCA(4,JROC) = 1 !  JUST TO GET PAST COMP 'STACKSZW'
            ITSKCA(5,JROC) = NSTKSZ
            ITSKCA(6,JROC) = 1 !  JUST TO GET PAST COMP 'STACKINW'
            ITSKCA(7,JROC) = NSTKINC

            CALL TSKSTART(ITSKCA(1,JROC), MT_ALLOT, FL2, ETOIB,
     &       USOIB, USNEW, THWNEW, IJS(IG), IJL(IG),
     &       NSIZE,NWORK)
          ENDDO

          CALL MT_ALLOT (FL2, ETOIB, USOIB, USNEW, THWNEW,
     &     IJS(IG), IJL(IG), NSIZE, NWORK)

          DO JROC=1,NTASKS-1
            CALL TSKWAIT(ITSKCA(1,JROC))
          ENDDO

          IF (ITEST.GE.2) THEN
            IF (IG.LE.ITESTB) WRITE(IU06,*)                          
     &       '   SUB. WAMASSI: SPECTRA UPDATED FOR BLOCK IG = ', IG   
          ENDIF

!           THE WIND FIELD CURRENTLY USED IS UPDATED.                   

          DO IJ=IJS(IG),IJL(IG)                                     
            U10OLD(IJ,IG) = USNEW(IJ)*U10OLD(IJ,IG)/USOLD(IJ,IG)        
            U10NEW(IJ) = U10OLD(IJ,IG)                                  
          ENDDO
          IF (ITEST.GE.2) THEN
            IF (IG.LE.ITESTB) WRITE(IU06,*)                             
     &      '   SUB. WAMASSI: WIND REPLACED FOR BLOCK IG = ', IG        
          ENDIF

!              OUTPUT OF OVERLAPPING LATITUDES AND                      
!              START OF OUTPUT OF ALL SPECTRA IN BLOCK.                 

          IF (IGL.NE.1) THEN                                          
            IF (IG.NE.1) CALL SPLITBL(FL2, IJS(IG), IJL2(IG), IU14)  
            WRITE(IU12) FL2      
            IF (ITEST.GE.2) THEN
              IF (ITESTB.GE.IG) WRITE(IU06,*)                       
     &         '   SUB. WAMASSI: OVERLAPPING LATITUDES',          
     &         ' SAVED AND START OF BLOCK OUTPUT.'                
            ENDIF
          ENDIF

!              MODEL OUTPUT INTEGRATED DATA ARE SAVED IN COMMON.        
!              SPECTRA ARE SAVED ON DISK AND PRINTED.                   

          IF (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO .OR.           
     &     CDTSPT.EQ.CDTPRO .OR. CDTSPS.EQ.CDTPRO ) THEN
            CALL OUTBS (FL2, FL3, IJS(IG), IJL(IG), IG, IGL,         
     &       IU27, IU28,THWOLD,USOLD,U10NEW,THWNEW,USNEW,TAUW,
     &       NSTART,NEND)

            IF (ITEST.GE.2) THEN
              IF (ITESTB.GE.IG) WRITE(IU06,*)                       
     &         '   SUB. WAMASSI: MODEL OUTPUT DONE IN SUB OUTBS'  
            ENDIF
          ENDIF
        ENDIF

!*    BRANCHING BACK TO 3.4 FOR NEXT BLOCK OF LATITUDES                 

        IF(LWCOU) THEN
          CALL OUTBETA ( IJS(IG), IJL(IG), IG, IGL,USNEW,Z0NEW,
     &     NSTART, NEND, BETAG, NLONW, NLATW)
        ENDIF


      ENDDO BLOCK

!     FREE UP SOME MEMORY

      DEALLOCATE(USOI,USME,WHOI,WHME,WHOIB,USOIB,ETOIB)

!*    3.5 SAVE FILES.                                                   
!         -----------                                                   

!     3.5.1 WRITE INTEGRATED DATA TO FILE AND/OR PRINTER                
!           DATA WERE COLLECTED INSIDE THE BLOCK LOOP.                  
!           --------------------------------------------                

      IF (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO) THEN                

!       PRINTS OUT NORMS
        CALL OUTWNORM

        CALL OUTINT (IU06, IU22, IU23, IU32, ITEST, MFDBSF)

        IF(LMESSPASS) THEN
          CALL MPE_BARRIER(IERR)
          IF(IERR.LT.0) CALL MPEI_ABORT('MPE_BARRIER ERROR in WAMASSI')
        ENDIF
        IF (ITEST.GE.2) WRITE(IU06,*) '   SUB. WAMASSI: OUTINT DONE'
      ENDIF

      IF (CDATEE.EQ.CDTPRO.OR.LDWRRE) THEN
C                                                                       

!*    3.5.2 STORE 2-D SPECTRA AT SELECTED GRID POINTS,                  
!           INTEGRATED PARAMETERS OF ENTIRE GRID.                       
!           ------------------------------------------                  

!           SAVE INTEGRATED DATA FILE.                                  

        IF (GFLAG20 .AND. .NOT. LFDB) THEN                             
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            CALL GSFILE (IU06, IU32, 0, CDTPRO, CDATEF, 'APP', 'S')
          ENDIF
        ENDIF
        IF (FFLAG20) THEN                                              
          IF((LMESSPASS.AND.(IRANK.EQ.NPROC)).OR.(.NOT.LMESSPASS)) THEN
            CALL GSFILE (IU06, IU22, 0, CDTPRO, CDATEF, 'AMP', 'S')
          ENDIF
        ENDIF

!           SAVE INTEGRATED SWELL AND SEA DATA FILE.                    
!           ----------------------------------------                    

        IF (FFLAG21) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.MIN(2,NPROC)))
     &     .OR.(.NOT.LMESSPASS)) THEN
            CALL GSFILE (IU06, IU23, 0, CDTPRO, CDATEF, 'ASS', 'S')
          ENDIF
        ENDIF

!           SAVE SPECTRA FILE.                                          

        IF (FFLAG25) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            CALL GSFILE (IU06, IU27, 0, CDTPRO, CDATEF, 'AUT', 'S')
          ENDIF
        ENDIF

!           SAVE SWELL SPECTRA FILE.                                    

        IF (FFLAG26) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            CALL GSFILE (IU06, IU28, 0, CDTPRO, CDATEF, 'ASW', 'S')
          ENDIF
        ENDIF
      ENDIF

      IF (CDTRES.EQ.CDTPRO.OR.CDATEE.EQ.CDTPRO.OR.CDTPRO.EQ.CDATER.OR.
     &    LDWRRE) THEN


!*    3.5.4 SAVE RESTART FIELDS.                                        
!           --------------------                                        

        IF ( (IREST.EQ.1.AND.(CDTPRO.EQ.CDATER .OR. CDTPRO.LE.CDATES))
     &       .OR.
     &        LDWRRE ) THEN 

!*          STRESS FIELD IS DUMPED FOR RESTART.
!           -----------------------------------
!           NOTE :
!           STRESS FIELD IS DUMPED FOR RESTART
!           ONLY IF CURRENT TIME CORRESPONDS TO CDATER

          IF((CDTPRO .EQ. CDATER  ).OR.
     &       (CDATER .EQ. '000000000000' ).OR.
     &        LDWRRE) THEN

            IF ( .NOT.LGRIBOUT ) THEN
              CALL SAVSTRESS(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,NSTART,NEND,
     &                       CDTPRO,CDATEF)
            ENDIF

          ENDIF

!           SAVE SPECTRUM
!           -------------

          IF ( LGRIBOUT ) THEN
            CALL SPEC2FDB (IU06, NPROC, IRANK, LWCOU, TH, FR,
     &                     FL1, NINF, NSUP, NANG, NFRE, NGX, NGY,
     &                     NLONRGG, KXLT, IXLG, NIBLO, NBLO,
     &                     AMONOP, XDELLA, NSTART, NEND,
     &                     NKSEK0, NKSEK1, NKSEK2, NKSEK3, NKSEK4,
     &                     NPSEK2, NPSEK3, KSEK0, KSEK1, KSEK2, KSEK3,
     &                     KSEK4, PSEK2, PSEK3, NTENCODE,
     &                     DEG, PPMISS, PPEPS, PPREC,
     &                     CDTPRO, CDATEF, CDATEA,
     &                     CLDOMAIN, MARSTYPE, YCLASS, YEXPVER,
     &                     CFDB2DSP, NENSFNB, NTOTENS, LFDBIOOUT,
     &                     ICPLEN, CPATH, ITEST, MFDBSF)
            ISTAT = IFLUSHFDB (MFDBSF)
            LLFLUSH = .FALSE.
          ELSE
            CALL SAVSPEC(FL1,NSTART,NEND,CDTPRO,CDATEF,CDATER,
     &                   MARSTYPE,YCLASS, YEXPVER, CFDB2DSP)
          ENDIF
          WRITE(IU06,*) ' '
          WRITE(IU06,*) '  WAVE SPECTRA FIELDS DISPOSED AT........',
     &                  ' CDTPRO  = ', CDTPRO
          WRITE(IU06,*) ' '


!*    3.5.5 UPDATE, WRITE AND SAVE WAMINFO FILE.
!           ------------------------------------

          IF (LDSTOP .AND. 
     &        (LMESSPASS .AND. IRANK .EQ. 1 .OR. .NOT. LMESSPASS )) THEN
            ICH = 7
            CALL DIFDATE (CDATEF,CDATEE,IFOREPD)
            CALL DIFDATE (CDTPRO,CDATEF,IANALPD)
            IF (CDTPRO.LE.CDATEF) THEN
              NFCST = 1
            ELSE
              NFCST = 0
              CALL DIFDATE (CDTPRO,CDATEE,IFOREPD)
            ENDIF
            ISTAT = 0
            IF (CDTPRO.EQ.CDATEE) ISTAT(1) = 1
            IDELWIN = IDELWI
            CBPLTDT = CDATEF
            CEPLTDT = CDATEE

            OPEN (IU04,FILE=CWI(1:ICPLEN+8),
     &            STATUS='UNKNOWN',FORM='FORMATTED')
            CALL WRITSTA (IU04, CDTPRO, CDATEE, IANALPD, IFOREPD,
     &                    IDELWIN, CDATER, CDATES, CBPLTDT, CEPLTDT,
     &                    IASSI, NFCST, ISTAT)

            CLOSE (IU04)
            WRITE(IU06,*) ' WAMINFO FILE WRITTEN FOR RESTART...',
     &                      ' CDTPRO  = ', CDTPRO
            WRITE(IU06,*) '                                    ',
     &                      ' CDATEF  = ', CDATEF
            WRITE(IU06,*) ' TO ', CWI(1:ICPLEN+8)
          ENDIF
        ENDIF

        IF(LMESSPASS) THEN
          CALL MPE_BARRIER(IERR)
          IF(IERR.LT.0) CALL MPEI_ABORT('MPE_BARRIER ERROR IN  WAMASSI')
        ENDIF
      ENDIF                                                            
      IF (LFDB .AND. LLFLUSH .AND. (
     &  (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO))) THEN
        ISTAT = IFLUSHFDB (MFDBSF)
        WRITE(IU06,*) ' DB ', MFDBSF , ' FLUSHED AT ',
     &    CDTPRO, ' FROM WAMASSI. '
        CALL FLUSH (IU06)
      ENDIF

!     DEALLOCATE OUTPUT ARRAYS

      IF(ALLOCATED(WHGTTG)) DEALLOCATE(WHGTTG)
      IF(ALLOCATED(WDIRTG)) DEALLOCATE(WDIRTG)
      IF(ALLOCATED(WPKFTG)) DEALLOCATE(WPKFTG)
      IF(ALLOCATED(WMNFTG)) DEALLOCATE(WMNFTG)
      IF(ALLOCATED(USTARG)) DEALLOCATE(USTARG)
      IF(ALLOCATED(UDIRG)) DEALLOCATE(UDIRG)
      IF(ALLOCATED(TAUWG)) DEALLOCATE(TAUWG)
      IF(ALLOCATED(CDG)) DEALLOCATE(CDG)
      IF(ALLOCATED(SMEANG)) DEALLOCATE(SMEANG)
      IF(ALLOCATED(U10G)) DEALLOCATE(U10G)

      IF(ALLOCATED(WHGTSG)) DEALLOCATE(WHGTSG)
      IF(ALLOCATED(WDIRSG)) DEALLOCATE(WDIRSG)
      IF(ALLOCATED(WMNFSG)) DEALLOCATE(WMNFSG)
      IF(ALLOCATED(WHGTWG)) DEALLOCATE(WHGTWG)
      IF(ALLOCATED(WDIRWG)) DEALLOCATE(WDIRWG)
      IF(ALLOCATED(WMNFWG)) DEALLOCATE(WMNFWG)

      RETURN                                                            
      END SUBROUTINE WAMASSI
