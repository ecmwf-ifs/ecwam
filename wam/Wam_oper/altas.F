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

      USE YOWALTAS , ONLY : INTLMAX  ,KMINLMAX ,KMAXLMAX  ,NOBSPE
      USE YOWCOUT  , ONLY : IPFGTBL
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL2     ,IJL
      USE YOWMAP   , ONLY : XDELLA   ,KXLT
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,AKMEAN
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    , NINF     ,NSUP     ,
     &            KTAG
      USE YOWMLTSK , ONLY : NCHUNK
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NIBLO
      USE YOWPCONS , ONLY : DEG      ,R
      USE YOWSPEC,   ONLY : NSTART   ,NEND     ,U10NEW   ,U10OLD   ,
     &            THWNEW   ,THWOLD   ,USNEW    ,USOLD    ,Z0NEW    ,
     &            Z0OLD    ,TAUW     ,FL1      ,FL2      ,SL       ,
     &            ROAIRN   ,ZIDLNEW
      USE YOWSTAT  , ONLY : CDTPRO   ,NTASKS   ,NSIZE
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWUNIT  , ONLY : IU11     ,IU12     ,IU13     ,IU14
      USE YOWWIND  , ONLY : WSPMIN

! ----------------------------------------------------------------------

!*    *ADDITIONAL PARAMETERS FOR MULTITASKING

      INTEGER ITSKCA(9,NTASKS)
      EXTERNAL MT_ALLOT

! ----------------------------------------------------------------------
!     DEFINE ALLOCATABLE ARRAYS

      REAL,ALLOCATABLE,DIMENSION(:) :: WHMOD
      REAL,ALLOCATABLE,DIMENSION(:) :: ETOIB

!*    VARIABLE.     TYPE.        PURPOSE.
!     ---------     -----        --------
!      *WHMOD*      REAL         MODEL WAVE HEIGHT (FIRST GUESS).
!      *ETOIB*      REAL         WAVE HEIGHT FROM O.I. UNTIL
!                                CONVERTED TOTAL ENERGY BEFORE UPDATE.

!-----------------------------------------------------------------------

      LOGICAL :: LLAK

      WRITE(IU06,*)'  -----ALT DATA ASSIMILATION BEGINS--------'

!    1. GET MODEL WAVE HEIGHT 
!       ---------------------

      IF(LMESSPASS) THEN
        LLAK=.FALSE.
        CALL FEMEAN (FL1, NSTART(IRANK), NEND(IRANK), EMEAN, FMEAN,
     &               AKMEAN, LLAK)

        ALLOCATE(WHMOD(NIBLO))

        DO IJ = NSTART(IRANK), NEND(IRANK)
          IF(EMEAN(IJ).GT.0.) THEN
            WHMOD(IJ)=4.*SQRT(EMEAN(IJ)) 
          ELSE
            WHMOD(IJ)=0.
          ENDIF
        ENDDO

!       make WHMOD a global array
        KTAG=KTAG+1
        CALL MPBCASTSCFLD(KTAG,NSTART,NEND,WHMOD) 
      ELSE
        write(*,*) 'ALTAS DOES NOT WORK FOR THIS CONFIGURATION !!!'
        WRITE(IU06,*) 'ALTAS DOES NOT WORK FOR THIS CONFIGURATION !!!'
        CALL ABORT1
      ENDIF

!*    2. MEASUREMENTS AND MODEL ARE MERGED BY OPTIMUM INTERPOLATION.
!        -----------------------------------------------------------

!     specify the correlation distance (DIST) and the maximum spreading
!     distance (DISTMAX)
      IF(XDELLA.LT.0.5) THEN
        DIST  = 150000./R
      ELSE
        DIST  = 300000./R
      ENDIF

      DISTMAX = SQRT(2.) * 2.*DIST

!     specify the model error
!     note that the observation error is now specified in rfl4wam
      SIGMOD=0.5

!     DETERMINE THE TABLE WHICH TELLS WHICH PE POTENTIALLY SHARES 
!     ALTIMETER DATA WITH OTHER PE'S.
      IF(.NOT. ALLOCATED(INTLMAX) ) THEN
        ALLOCATE(INTLMAX(NPROC,NPROC))
        ALLOCATE(KMINLMAX(NPROC))
        ALLOCATE(KMAXLMAX(NPROC))

        CALL MPPEWITHINDIST(DISTMAX,INTLMAX,KMINLMAX,KMAXLMAX)

      ENDIF

!*    2.1 THINNED MEASUREMENTS ARE READ IN
!         --------------------------------

      ALLOCATE (NOBSPE(NPROC))
      CALL GRFIELD(WHMOD,SIGMOD,DISTMAX)

      IF (ITEST.GE.2)
     &   WRITE(IU06,*) '   SUB. ALTAS: MEASUREMENTS READ IN'


!*    2.2 OI.
!         ---

      ALLOCATE(ETOIB(NINF:NSUP))

      CALL OIFIELD (WHMOD, ETOIB, SIGMOD, DIST, DISTMAX)

      DEALLOCATE(WHMOD)

      IF (ITEST.GE.2)
     &   WRITE(IU06,*) '   SUB. ALTAS: OIFIELD DONE FOR WAVE HEIGHTS'

!*    3. ANALYSING THE SPECTRA. 
!        ----------------------

!*    3.1 LOOP OVER THE BLOCKS.
!         ---------------------
      IG=1

!*    3.1.1 TRANSFERS RESULTS OF OPTIMAL INT TO BLOCKS.
!           -------------------------------------------

      DO IJ=IJS(IG),IJL(IG)
        USNEW(IJ) = USOLD(IJ,IG)
        THWNEW(IJ)= THWOLD(IJ,IG)
        Z0NEW(IJ)= Z0OLD(IJ,IG)
      ENDDO


      CINV16=1./16.
      DO IJ=IJS(IG),IJL(IG)
        IF (ETOIB(IJ).GT.0.) THEN
          ETOIB(IJ) = ETOIB(IJ)**2*CINV16
        ELSE
          ETOIB(IJ) = -99.
        ENDIF
      ENDDO

      CALL UPDATE (FL1, ETOIB, USNEW, THWNEW, IJS(IG), IJL(IG))

      DEALLOCATE(ETOIB)

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. ALTAS: SPECTRA UPDATED '
      ENDIF

!     THE WIND FIELD CURRENTLY USED IS UPDATED.

      CALL UPDATEWD(FL1, SL, IJS(IG),IJL(IG),IG, THWNEW, USNEW,
     &              USOLD, TAUW, Z0NEW, U10NEW, ROAIRN, ZIDLNEW)

      DO IJ=IJS(IG),IJL(IG)
        U10OLD(IJ,IG) = U10NEW(IJ)
        Z0OLD(IJ,IG) = Z0NEW(IJ)
        USOLD(IJ,IG) = USNEW(IJ)
      ENDDO

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. ALTAS: WIND REPLACED'
      ENDIF

      DEALLOCATE (NOBSPE)

      WRITE(IU06,*) '  ALTAS ENDS NORMALLY'
      WRITE(IU06,*) ' '
      CALL FLUSH(IU06)

      RETURN 
      END SUBROUTINE ALTAS
