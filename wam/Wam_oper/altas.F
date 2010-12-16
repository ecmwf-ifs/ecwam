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
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL2     ,IJL      ,
     &            COSPH
      USE YOWMAP   , ONLY : XDELLA   ,KXLT
      USE YOWMAP   , ONLY : XDELLA   ,KXLT     ,AMOSOP
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    , NINF    ,NSUP     ,
     &            KTAG     ,MPMAXLENGTH
      USE YOWMLTSK , ONLY : NCHUNK
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NIBLO    ,NANG     ,
     &            NFRE
      USE YOWPCONS , ONLY : DEG      ,R
      USE YOWSPEC,   ONLY : NSTART   ,NEND     ,U10NEW   ,U10OLD   ,
     &            THWNEW   ,THWOLD   ,USNEW    ,USOLD    ,Z0NEW    ,
     &            Z0OLD    ,TAUW     ,FL1      ,FL2      ,SL       ,
     &            ROAIRN   ,ZIDLNEW
      USE YOWSTAT  , ONLY : CDTPRO   ,NTASKS   ,NSIZE    ,NWAM_BLKS
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWUNIT  , ONLY : IU11     ,IU12     ,IU13     ,IU14
      USE YOWWIND  , ONLY : WSPMIN
      USE MPL_MODULE

! ----------------------------------------------------------------------

!*    *ADDITIONAL PARAMETERS FOR MULTITASKING

      INTEGER ITSKCA(9,NTASKS)
      EXTERNAL MT_ALLOT

! Mod for OPENMP
      INTEGER JKGLO,KIJS,KIJL,NPROMA
! End Mod for OPENMP

      LOGICAL NEWREAD, NEWFIL
! ----------------------------------------------------------------------
!     DEFINE ALLOCATABLE ARRAYS

      INTEGER,DIMENSION(NPROC) :: ISENDLEN
      INTEGER,ALLOCATABLE :: ISENDREQ(:)

      REAL,ALLOCATABLE,DIMENSION(:) :: DIST 
      REAL,ALLOCATABLE,DIMENSION(:) :: WHMOD
      REAL,ALLOCATABLE,DIMENSION(:) :: ETOIB
      REAL,ALLOCATABLE :: ZCOMBUFS(:,:)
      REAL,ALLOCATABLE :: ZCOMBUFR(:)

!-----------------------------------------------------------------------

      LOGICAL :: LLAK

      WRITE(IU06,*)' '
      WRITE(IU06,*)'  -----ALT DATA ASSIMILATION BEGINS--------'

      IG=1
      IF(.NOT.LMESSPASS) THEN
        write(*,*) 'ALTAS DOES NOT WORK FOR THIS CONFIGURATION !!!'
        WRITE(IU06,*) 'ALTAS DOES NOT WORK FOR THIS CONFIGURATION !!!'
        CALL ABORT1
      ENDIF

!*    1. MEASUREMENTS AND MODEL ARE MERGED BY OPTIMUM INTERPOLATION.
!        -----------------------------------------------------------

!     specify the correlation distance (DIST) and the maximum spreading
!     factor  EXTENDMAX and get first guess wave height (WHMOD)

!!! if DIST were to be reset during the run, then MPPEWITHINDIST has to be called
!!! again !!!!
      ALLOCATE(DIST(NSTART(IRANK):NEND(IRANK)))

!!! for now, we will disactivate the geographically varying DIST

      IF(XDELLA.LT.0.3) THEN
        XLMIN=150000.
        XLMAX=0.
      ELSE
        XLMIN=300000.
        XLMAX=0.
      ENDIF

      SIGH=0.75

      ALLOCATE(WHMOD(NIBLO))

      NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1

      LLAK=.FALSE.
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL,IJ,KJ)
      DO JKGLO=IJS(IG),IJL(IG),NPROMA
        KIJS=JKGLO
        KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
        CALL FEMEAN (FL1(KIJS:KIJL,1:NANG,1:NFRE), KIJS, KIJL,
     &               EMEAN(KIJS), FMEAN(KIJS),LLAK)

        DO IJ=KIJS,KIJL
          IF(EMEAN(IJ).GT.0.) THEN
            WHMOD(IJ)=4.*SQRT(EMEAN(IJ)) 
          ELSE
            WHMOD(IJ)=0.
          ENDIF
        ENDDO

        DO IJ=KIJS,KIJL
          KJ=KXLT(IJ,IG)
          DIST(IJ) = (XLMIN+XLMAX*COSPH(KJ))/R
        ENDDO
      ENDDO
!$OMP       END PARALLEL DO

      EXTENDMAX = SQRT(2.) * 2.

      DISTMAX = 0.
      DO IJ = NSTART(IRANK), NEND(IRANK)
         DISTMAX = MAX(DISTMAX,DIST(IJ))
      ENDDO
!!    DMAX IS THE ABSOLUTE MAXIMUM OF DISTMAX
      DMAX = EXTENDMAX*(XLMIN+XLMAX)

!     DETERMINE THE TABLE WHICH TELLS WHICH PE POTENTIALLY SHARES 
!     ALTIMETER DATA WITH OTHER PE'S.
!     IT IS BASED ON THE VALUES OF DIST BEFORE REDUCTION FOR LOW WAVE HEIGHT
      IF(.NOT. ALLOCATED(INTLMAX) ) THEN
        ALLOCATE(INTLMAX(NPROC,NPROC))
        ALLOCATE(KMINLMAX(NPROC))
        ALLOCATE(KMAXLMAX(NPROC))

        CALL MPPEWITHINDIST(DMAX,INTLMAX,KMINLMAX,KMAXLMAX)
      ENDIF


!     SPECIFY THE MODEL ERROR:

!     note that the observation error is now specified in rfl4wam
      SIGMOD=0.5

!     REDUCE DIST FOR LOW WAVE HEIGHT:
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL,IJ,REDUC)
      DO JKGLO=IJS(IG),IJL(IG),NPROMA
        KIJS=JKGLO
        KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
        DO IJ=KIJS,KIJL
          REDUC=1.-EXP(-MIN((WHMOD(IJ)/SIGH)**2,50.)) 
          DIST(IJ) = MAX(0.1,REDUC)*DIST(IJ)
        ENDDO
      ENDDO
!$OMP       END PARALLEL DO


!     DISTRIBUTE WHMOD TO WHERE IT IS NEEDED

      IF(NPROC.GT. 1) THEN

        KTAG=KTAG+1
        NDISPE=0
        DO IR = 1, NPROC
          IF(INTLMAX(IR,IRANK).EQ.1) THEN
            NDISPE=NDISPE+1
            ISENDLEN(IR)=MPMAXLENGTH
          ELSE
            ISENDLEN(IR)=0
          ENDIF
        ENDDO
!       PACKING THE CONTRIBUTION WHICH MIGHT BE COMMON
        ALLOCATE(ZCOMBUFS(MPMAXLENGTH,NDISPE))
        ALLOCATE(ISENDREQ(NDISPE))
        IDISPE=0
        DO IR = 1, NPROC
          IF(INTLMAX(IR,IRANK).EQ.1) THEN
            IDISPE=IDISPE+1
            ICOUNT = 0
            DO IJ = NSTART(IRANK), NEND(IRANK)
              ICOUNT = ICOUNT + 1
              ZCOMBUFS(ICOUNT,IDISPE) = WHMOD(IJ)
           ENDDO
          ENDIF
        ENDDO

!       SENDING THE CONTRIBUTION WHICH MIGHT BE COMMON
!       SEND NON BLOCKING THE BUFFERS

        IDISPE=0
        CALL GSTATS(618,0)
        DO IR = 1, NPROC
          IF(INTLMAX(IR,IRANK).EQ.1) THEN
            IDISPE=IDISPE+1
            ICOUNT = ISENDLEN(IR) 
            CALL MPL_SEND(ZCOMBUFS(1:ICOUNT,IDISPE),
     &                    KDEST=IR,KTAG=KTAG,
     &                    KMP_TYPE=JP_NON_BLOCKING_STANDARD,
     &                    KREQUEST=ISENDREQ(IDISPE),
     &                    CDSTRING='ALTAS:')
          ENDIF
        ENDDO
        CALL GSTATS(618,1)

!       RECEIVING THE CONTRIBUTIONS WHICH MIGHT BE COMMON
!       (in whatever order they arrive.)
        ALLOCATE(ZCOMBUFR(MPMAXLENGTH))
        DO IPR =1,NDISPE
            CALL GSTATS(618,0)
            CALL MPL_RECV(ZCOMBUFR(1:MPMAXLENGTH),KFROM=IR,KTAG=KTAG,
     &             KOUNT=KRCOUNT,KRECVTAG=KRTAG,CDSTRING='ALTAS:')
            CALL GSTATS(618,1)
            IF(KRTAG.NE.KTAG) CALL MPL_ABORT
     &      ('MPL_RECV ERROR in ALTAS: MISMATCHED TAGS' )
             ICOUNT = 0
             DO IJ = NSTART(IR), NEND(IR)
               ICOUNT = ICOUNT + 1
               WHMOD(IJ)=ZCOMBUFR(ICOUNT)
             ENDDO
        ENDDO
        DEALLOCATE(ZCOMBUFR)

!       WAIT ANY OUTSTANDING SENDS TO COMPLETE
        CALL GSTATS(618,0)
        IF(MPMAXLENGTH.GT.0) THEN
          CALL MPL_WAIT(ZCOMBUFS(1:MPMAXLENGTH,1),KREQUEST=ISENDREQ,
     &                  CDSTRING='ALTAS')
        ENDIF
        CALL GSTATS(618,1)
        IF(ALLOCATED(ISENDREQ)) DEALLOCATE(ISENDREQ)
        KTAG = KTAG + 1
        IF(ALLOCATED(ZCOMBUFS)) DEALLOCATE(ZCOMBUFS)

      ENDIF


!*    2.2 THINNED MEASUREMENTS ARE READ IN
!         --------------------------------

      ALLOCATE (NOBSPE(NPROC))
      CALL GRFIELD(WHMOD,SIGMOD,DMAX)

      IF (ITEST.GE.2)
     &   WRITE(IU06,*) '   SUB. ALTAS: MEASUREMENTS READ IN'


!*    2.3 OI.
!         ---

      ALLOCATE(ETOIB(NSTART(IRANK):NEND(IRANK)))

      CALL OIFIELD (WHMOD, ETOIB, SIGMOD, DIST, DISTMAX, EXTENDMAX)

      DEALLOCATE(DIST)
      DEALLOCATE(WHMOD)

      IF (ITEST.GE.2)
     &   WRITE(IU06,*) '   SUB. ALTAS: OIFIELD DONE FOR WAVE HEIGHTS'

!*    3. ANALYSING THE SPECTRA. 
!        ----------------------

!*    3.1 LOOP OVER THE BLOCKS.
!         ---------------------

!*    3.1.1 TRANSFERS RESULTS OF OPTIMAL INT TO BLOCKS.
!           -------------------------------------------

      NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1

      CINV16=1./16.

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL,IJ)
      DO JKGLO=IJS(IG),IJL(IG),NPROMA
        KIJS=JKGLO
        KIJL=MIN(KIJS+NPROMA-1,IJL(IG))

        DO IJ=KIJS,KIJL
          USNEW(IJ) = USOLD(IJ,IG)
          THWNEW(IJ)= THWOLD(IJ,IG)
          Z0NEW(IJ)= Z0OLD(IJ,IG)
        ENDDO

        DO IJ=KIJS,KIJL
          IF (ETOIB(IJ).GT.0.) THEN
            ETOIB(IJ) = ETOIB(IJ)**2*CINV16
          ELSE
            ETOIB(IJ) = -99.
          ENDIF
        ENDDO

!       UPDATE SPECTRA
        CALL UPDATE (FL1(KIJS:KIJL,1:NANG,1:NFRE),
     &               ETOIB(KIJS), USNEW(KIJS), THWNEW(KIJS), KIJS, KIJL)

!       THE WIND FIELD CURRENTLY USED IS UPDATED.
        CALL UPDATEWD(FL1(KIJS:KIJL,1:NANG,1:NFRE),
     &                SL(KIJS:KIJL,1:NANG,1:NFRE),
     &                KIJS,KIJL,
     &                THWNEW(KIJS), USNEW(KIJS),
     &                USOLD(KIJS:KIJL,IG), TAUW(KIJS:KIJL,IG),
     &                Z0NEW(KIJS), U10NEW(KIJS),
     &                ROAIRN(KIJS), ZIDLNEW(KIJS) )

        DO IJ=KIJS,KIJL
          U10OLD(IJ,IG) = U10NEW(IJ)
          Z0OLD(IJ,IG) = Z0NEW(IJ)
          USOLD(IJ,IG) = USNEW(IJ)
        ENDDO
      ENDDO
!$OMP       END PARALLEL DO

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. ALTAS: SPECTRA AND WIND UPDATED '
      ENDIF

      DEALLOCATE(ETOIB)
      DEALLOCATE (NOBSPE)

      WRITE(IU06,*) '  ALTAS ENDS NORMALLY'
      WRITE(IU06,*) ' '
      CALL FLUSH(IU06)

      RETURN 
      END SUBROUTINE ALTAS
