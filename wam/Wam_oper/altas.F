      SUBROUTINE ALTAS

! ----------------------------------------------------------------------

!****  *ALTAS* - WAVE DATA ASSIMILATION USING ALTIMETER DATA ALONE. 

!      J. BIDLOT       ECMWF     MAY     1999  SPLIT FROM WAMASSI 

!     PURPOSE.                                                          
!     --------                                                          

!         TO ANALYSE THE WAVE FIELD SUBSTITUTING THE FIRST GUESS        
!         SPECTRA WITH ANALYSED SPECTRA AND FIRST GUESS VALUES          
!         WITH ANALYSED VALUES IN THE GLOBAL GRID                       

!*    INTERFACE.                                                        
!     ----------                                                        

!     *CALL* *ALTAS

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
!        *MAKEBLO*  - MAKES BLOCKS FROM GRIDDED FIEDS.                  
!        *MPGATHEROIFLD - COLLECTS FIELD PRODUCED BY THE O.I. ON PE'S.
!        *OIFIELD*  -  MERGES MODEL AND MEASUREMENTS BY O.I.            
!        *PEAKFR*   - COMUTE PEAK FREQUENCY.                            
!        *SEMEAN*   - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.   
!        *SERMAT*   - UTILITY OF OIFIELD TO SUM TWO MATRICES.           
!        *SPLITBL*  - STORING OVERLAPPING LATITUDES OF BLOCKS.          
!        *STHQ*     - COMPUTATION OF MEAN WAVE DIRECTION AT EACH        
!                     GRID POINT.                                       
!        *UPDATE*   - ANALYSES THE WAM MODEL SPECTRA.                   
!        *UPWSPEC*  - PRODUCES THE UPDATED SPECTRUM.                    
!        *WSMFEN*   - UTILITY OF FWSEA , COMPUTING  WINDSEA ENERGY AND  
!                     DURATION IN EACH GRID POINT.                      

!     REFERENCES.                                                       
!     -----------                                                       

!          NONE                                                         

! ----------------------------------------------------------------------

      USE YOWCOUT  , ONLY : IPFGTBL
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL2     ,IJL
      USE YOWINTP  , ONLY : WHGTTG   ,USTARG   ,CDG
      USE YOWMAP   , ONLY : XDELLA
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : NINF     ,NSUP     ,KTAG
      USE YOWMLTSK , ONLY : NCHUNK
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWSPEC,   ONLY : U10OLD   ,THWNEW   ,THWOLD   ,USNEW    ,
     &            USOLD    ,FL1      ,FL2
      USE YOWSTAT  , ONLY : CDTPRO   ,NTASKS   ,NSIZE
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWUNIT  , ONLY : IU11     ,IU12     ,IU13     ,IU14

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
      LOGICAL LWAVE

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

      CALL GRFIELD (WHME, USME)

      IF (ITEST.GE.2)
     &   WRITE(IU06,*) '   SUB. ALTAS: MEASUREMENTS ---> GRID'

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
          WRITE(IU06,*) '*    FATAL ERROR IN SUB. ALTAS          *'
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

      LWAVE=.FALSE.
      CALL OIFIELD (USME, USTARG, USOI, SIGOBS, LWAVE, SIGMOD,
     &              LMAX, NDIM2)

      IF (ITEST.GE.2)
     &   WRITE(IU06,*) '   SUB. ALTAS: OIFIELD DONE FOR WINDS'

      LWAVE=.TRUE.
      CALL OIFIELD (WHME, WHGTTG, WHOI, SIGOBS, LWAVE, SIGMOD,
     &              LMAX, NDIM2)


      IF (ITEST.GE.2)
     &   WRITE(IU06,*) '   SUB. ALTAS: OIFIELD DONE FOR WAVE HEIGHTS'

!*    3. ANALYSING THE SPECTRA. 
!        ----------------------

      IF(.NOT.LMESSPASS) THEN
        REWIND IU12
        IF (IGL.NE.1) THEN
          REWIND IU11
          REWIND IU12
          REWIND IU13
          REWIND IU14
          IFIL=IU11
          IU11=IU12
          IU12=IFIL
          IFIL=IU13
          IU13=IU14
          IU14=IFIL
          READ(IU11) FL1
          IF (ITEST.GE.2)
     &     WRITE(IU06,*) '   SUB. ALTAS: START OF READING',
     &     ' FIRST BLOCK OF SPECTRA'
        ENDIF
      ENDIF

!*    3.1 LOOP OVER THE BLOCKS.
!         ---------------------
      BLOCK : DO IG=1,IGL

!*      3.1.1 TRANSFERS RESULTS OF OPTIMAL INT TO BLOCKS.
!             -------------------------------------------
        CALL MAKEBLO (USOIB(IJS(IG)), USOI, IJS(IG), IJL(IG), IG)
        CALL MAKEBLO (WHOIB(IJS(IG)), WHOI, IJS(IG), IJL(IG), IG)

        DO IJ=IJS(IG),IJL(IG)
          USNEW(IJ) = USOLD(IJ,IG)
          THWNEW(IJ)= THWOLD(IJ,IG)
        ENDDO
        DO IJ=IJS(IG),IJL(IG)
          IF (WHOIB(IJ).GT.0.) THEN
            ETOIB(IJ) = WHOIB(IJ)**2/16.
          ELSE
            ETOIB(IJ) = -99.
          ENDIF
        ENDDO

        IF (ITEST.GE.2) THEN
          IF (IG.LE.ITESTB) WRITE(IU06,*) 
     &     '   SUB. ALTAS: FIELDS BLOCKED FOR IG = ', IG
        ENDIF

        IF (MOD(IG,2).EQ.1) THEN
!*      3.1.2  ODD  BLOCKS.
!              ------------
          IF (LMESSPASS) THEN

            CALL UPDATE (FL1, ETOIB, USOIB, USNEW, THWNEW,
     &                   IJS(IG), IJL(IG))

            IF (ITEST.GE.2) THEN
              IF (IG.LE.ITESTB) WRITE(IU06,*)
     &         '   SUB. ALTAS: SPECTRA UPDATED '
            ENDIF

          ELSE
!*          START READING FOR NEXT BLOCK.
            IF (IG.NE.IGL) READ(IU11) FL2
!           UPDATE SPECTRA AND USTAR.
!           (PARALLEL PROCESSING - HENCE MT_ALLOT).
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
     &         '   SUB. ALTAS: SPECTRA UPDATED FOR BLOCK IG = ', IG
            ENDIF
!           OUTPUT OF OVERLAPPING LATITUDES AND
!           START OF OUTPUT OF ALL SPECTRA IN BLOCK.
            IF (IGL.NE.1) THEN
              IF (IG.NE.1) CALL SPLITBL(FL1, IJS(IG), IJL2(IG), IU14)
              WRITE(IU12) FL1
              IF (ITEST.GE.2) THEN
                IF (ITESTB.GE.IG) WRITE(IU06,*)
     &           '   SUB. ALTAS: OVERLAPPING LATITUDES',
     &           ' SAVED AND START OF BLOCK OUTPUT.'
              ENDIF
            ENDIF
          ENDIF

!         THE WIND FIELD CURRENTLY USED IS UPDATED.

          DO IJ=IJS(IG),IJL(IG)
            U10OLD(IJ,IG) = USNEW(IJ)*U10OLD(IJ,IG)/USOLD(IJ,IG)
            USOLD(IJ,IG) = USNEW(IJ)
          ENDDO
          IF (ITEST.GE.2) THEN
            IF (IG.LE.ITESTB) WRITE(IU06,*)
     &      '   SUB. ALTAS: WIND REPLACED FOR BLOCK IG = ', IG
          ENDIF

        ELSE 
!*      3.1.3  EVEN BLOCKS.
!              ------------
!         START READING FOR NEXT BLOCK.
          IF (IG.NE.IGL) READ(IU11) FL1
!         UPDATE SPECTRA AND USTAR.
!         (PARALLEL PROCESSING - HENCE MT_ALLOT).
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
     &       '   SUB. ALTAS: SPECTRA UPDATED FOR BLOCK IG = ', IG
          ENDIF

!         THE WIND FIELD CURRENTLY USED IS UPDATED.
          DO IJ=IJS(IG),IJL(IG)
            U10OLD(IJ,IG) = USNEW(IJ)*U10OLD(IJ,IG)/USOLD(IJ,IG)
            USOLD(IJ,IG) = USNEW(IJ)
          ENDDO
          IF (ITEST.GE.2) THEN
            IF (IG.LE.ITESTB) WRITE(IU06,*)
     &      '   SUB. ALTAS: WIND REPLACED FOR BLOCK IG = ', IG
          ENDIF

!         OUTPUT OF OVERLAPPING LATITUDES AND
!         START OF OUTPUT OF ALL SPECTRA IN BLOCK.
          IF (IGL.NE.1) THEN
            IF (IG.NE.1) CALL SPLITBL(FL2, IJS(IG), IJL2(IG), IU14)
            WRITE(IU12) FL2
            IF (ITEST.GE.2) THEN
              IF (ITESTB.GE.IG) WRITE(IU06,*)
     &         '   SUB. ALTAS: OVERLAPPING LATITUDES',
     &         ' SAVED AND START OF BLOCK OUTPUT.'
            ENDIF
          ENDIF

        ENDIF

!*    BRANCHING BACK TO 3.1 FOR NEXT BLOCK OF LATITUDES
      ENDDO BLOCK

      DEALLOCATE(USOI,USME,WHOI,WHME,WHOIB,USOIB,ETOIB)

      RETURN 
      END SUBROUTINE ALTAS
