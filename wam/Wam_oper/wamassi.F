      SUBROUTINE WAMASSI(BETAG,NLONW,NLATW,LDSTOP,LDWRRE)
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C****  *WAMASSI* - SUPERVISES EXECUTION OF MAIN MODULES                 
C****              OF THE WAVE DATA ASSIMILATION.                       
C                                                                       
C      P. LIONELLO     ECMWF     APRIL   1990
C      J. BIDLOT       ECMWF     FEBRARY 1996 MESSAGE PASSING.
C      J. BIDLOT       ECMWF     FEBRARY 1997 MODULE and ATMOSPHERIC
C                                COUPLING.
C      J. BIDLOT       ECMWF     MARCH   1997  ADD SAVSTRESS AND SAVSPEC
C      B. HANSEN       ECMWF     APRIL   1997  SIGNAL HANDLING.
C      B. HANSEN       ECMWF     JANUARY 1998  1.1 OUTSOURCED.
C                                                                       
C     PURPOSE.                                                          
C     --------                                                          
C                                                                       
C         TO ANALYSE THE WAVE FIELD SUBSTITUTING THE FIRST GUESS        
C         SPECTRA WITH ANALYSED SPECTRA AND FIRST GUESS VALUES          
C         WITH ANALYSED VALUES IN THE GLOBAL GRID                       
C                                                                       
C*    INTERFACE.                                                        
C     ----------                                                        
C                                                                       
C     *CALL* *WAMASSI(BETAG,NLONW,NLATW,LDSTOP,LDWRRE)
C        *BETAG*     REAL      CHARNOCK PARAMETER TO BE PASSED BACK TO
C                              THE ATMOSPHERIC MODEL
C        *NLONW*     INTEGER   FIRST DIMENSION OF BETAG
C        *NLATW*     INTEGER   SECOND DIMENSION OF BETAG
C        *LDSTOP*    LOGICAL   SET .TRUE. IF STOP SIGNAL RECEIVED.
C        *LDWRRE*    LOGICAL   SET .TRUE. IF RESTART SIGNAL RECEIVED.
C
C     METHOD.                                                           
C     -------                                                           
C                                                                       
C     EXTERNALS.                                                        
C     ----------                                                        
C                                                                       
C        *ANALYSE*  - UTILITY OF OIFIELD . IT PRODUCES THE OPTIMALLY    
C                     INTERPOLATED VALUE AT A SINGLE GRIDPOINT.         
C        *FDUR*     - COMPUTES WINDSEA DURATION.                        
C        *FEMEAN*   - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT. 
C        *FUSTAR*   - COMPUTES ANALYSED USTAR ACCORDING TO ANALYSED     
C                     WINDSEA ENERGY.                                   
C        *FWSEA*    - DETERMINES IF WINDSEA IS PRESENT IN THE SPECTRUM,,
C                     COMPUTING EVENTUALLY WINDSEA ENERGY AND MEAN FREQ.
C        *GRFIELD*  - READS GRIDDED MEASUREMENTS AND DISTRIBUTES THEM.
C        *GRDATA*   - TRANSFERS MEASUREMENTS TO THE MODEL GRID.         
C        *GSFILE*   - ROUTINE TO DYNAMICALLY FETCH OR DISPOSE FILES.    
C        *INCDATE*  - UPDATE DATE TIME GROUP.                           
C        *MAKEBLO*  - MAKES BLOCKS FROM GRIDDED FIEDS.                  
C        *MAKEGRID* - MAKES GRIDDED FIELDS.                             
C        *OIFIELD*  -  MERGES MODEL AND MEASUREMENTS BY O.I.            
C        *OUTAN*    - TRANSFERS INTEGRATED ANALYSED QUANTITIES TO FILE. 
C        *PEAKFR*   - COMUTE PEAK FREQUENCY.                            
C        *SAVSTRESS*- DISPOSE STRESS/WIND RESTART FILES.
C        *SAVSPEC   - DISPOSE SPECTRUM RESTART FILES.
C        *SEMEAN*   - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.   
C        *SERMAT*   - UTILITY OF OIFIELD TO SUM TWO MATRICES.           
C        *SPLITBL*  - STORING OVERLAPPING LATITUDES OF BLOCKS.          
C        *STHQ*     - COMPUTATION OF MEAN WAVE DIRECTION AT EACH        
C                     GRID POINT.                                       
C        *UPDATE*   - ANALYSES THE WAM MODEL SPECTRA.                   
C        *UPWSPEC*  - PRODUCES THE UPDATED SPECTRUM.                    
C        *WSMFEN*   - UTILITY OF FWSEA , COMPUTING  WINDSEA ENERGY AND  
C                     DURATION IN EACH GRID POINT.                      
C        *OUTWNORM* - COMPUTES A FEW NORMS OF GRIDDED FIELDS
C                                                                       
C     REFERENCES.                                                       
C     -----------                                                       
C                                                                       
C          NONE                                                         
C                                                                       
C ----------------------------------------------------------------------
C
      USE WAVE_MODULE
C
C ----------------------------------------------------------------------
C                                                                       
#include "param.h"
C                                                                       
#include "txtpara.h"
C                                                                       
#include "comcoup.h"
C
#include "txtcoup.h"
C
#include "comcout.h"
C                                                                       
#include "txtcout.h"
C                                                                       
#include "comfred.h"
C                                                                       
#include "txtfred.h"
C                                                                       
#include "comgrid.h"
C                                                                       
#include "txtgrid.h"
C                                                                       
#include "comintp.h"
C                                                                       
#include "txtintp.h"
C                                                                       
#include "commap.h"
C                                                                       
#include "txtmap.h"
C                                                                       
#include "commean.h"
C                                                                       
#include "txtmean.h"
C                                                                       
#include "comover.h"
C                                                                       
#include "txtover.h"
C                                                                       
#include "comstat.h"
C                                                                       
#include "txtstat.h"
C                                                                       
#include "comtest.h"
C                                                                       
#include "txttest.h"
C                                                                       
#include "comunit.h"
C                                                                       
#include "txtunit.h"
C                                                                       
#include "comwind.h"
C                                                                       
#include "txtwind.h"
C                                                                       
#include "parcons.h"
C
#include "commltsk.h"
C
C*    *ADDITIONAL PARAMETERS FOR MULTITASKING
C
      INTEGER ITSKCA(9,NTASKS)
      EXTERNAL MT_ALLOT
C
#include "commesspass.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
#include "comwami.h"
C
#include "comtext.h"
C
C ----------------------------------------------------------------------
C     define allocatable arrays

      REAL,ALLOCATABLE,DIMENSION(:,:) :: USOI, USME, WHOI, WHME
      REAL,ALLOCATABLE,DIMENSION(:) :: WHOIB, USOIB, ETOIB 

      CHARACTER TITL*100                                                
C                                                                       
C*    VARIABLE.     TYPE.        PURPOSE.                               
C     ---------     -----        --------                               
C      *USOI*        REAL        USTAR FIELD FROM O.I.                  
C      *USME*        REAL        USTAR FIELD FROM MEASUREMENTS          
C      *WHOI*        REAL        WAVE HEIGHT FIELD FROM O.I.            
C      *WHME*        REAL        WAVE HEIGHT FIELD FROM MEASUREMENTS    
C      *WHOIB*       REAL        WAVE HEIGHT FROM O.I. (BLOCKED)        
C      *USOIB*       REAL        USTAR FROM O.I. (BLOCKED)              
C      *ETOIB*       REAL        TOTAL ENERGY FROM O.I.                 
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      NAMELIST /OIP/ LMAX, SIGOBS, SIGMOD
C
C*    *NAMELIST* *OIP* PARAMETERS FOR THE OPTIMUM INTERPOLARTION.
C
      LOGICAL LOIP
C
C*     VARIABLE.   TYPE.     PURPOSE.
C      ---------   -------   --------
C      *LMAX*      INTEGER   MAXIMUM NUMBER OF GRID STEPS TO SPREAD
C                            MEASURED HS AND USTAR
C      *SIGOBS*    REAL      MEASUREMENT SCATTER.
C      *SIGMOD*    REAL      MODEL SCATTER.
C      *LOIP*      LOGICAL   TRUE IF INQUIRE ON FILE 'WTR' SUCCESSFULL.
C                                                                       
C ----------------------------------------------------------------------
      REAL BETAG(NLONW,NLATW)
C ----------------------------------------------------------------------
C
      CHARACTER*1 CMODE
      CHARACTER*80 FILENAME
      LOGICAL LOPENED, LDSTOP, LDWRRE, LLFAILED, LLFLUSH
      WRITE (IU06,*) '  '                                               
      WRITE (IU06,*)                                                    
     1   ' START OF DATA ASSIMILATION: DATE IS CDTPRO: ', CDTPRO        

      LLFLUSH = .TRUE.
      ALLOCATE(WHME(NGX,NGY), USME(NGX,NGY))
      ALLOCATE(WHOI(NGX,NGY), USOI(NGX,NGY))
      ALLOCATE(WHOIB(NINF:NSUP),USOIB(NINF:NSUP),ETOIB(NINF:NSUP))
C
C    0. COLLECT GRIDDED FIELDS USED BY THE ASSIMILATION
C       -----------------------------------------------

      IF(LMESSPASS) THEN
        KTAG=KTAG+1
        CALL MPGATHEROIFLD(IPFGTBL(1),KTAG,NSTART,NEND,WHGTTG)
        KTAG=KTAG+1
        CALL MPGATHEROIFLD(IPFGTBL(4),KTAG,NSTART,NEND,USTARG)
        KTAG=KTAG+1
        CALL MPGATHEROIFLD(IPFGTBL(7),KTAG,NSTART,NEND,CDG)
      ENDIF

C                                                                       
C*    1. MEASUREMENT ARE TRANSFERRED TO THE MODEL GRID.                 
C        ----------------------------------------------                 
C                                                                       
      CALL GRFIELD (WHME, USME, NSTART, NEND)
                                                                        
      IF (ITEST.GE.2)                                                   
     1   WRITE(IU06,*) '   SUB. WAMASSI: MEASUREMENTS ---> GRID'        
C                                                                       
C*    2. MEASUREMENTS AND MODEL ARE MERGED BY OPTIMUM INTERPOLATION.    
C        -----------------------------------------------------------    
C                                                                       
C*    2.1 GET PARAMETERS FROM NAMELIST OIP IF PRESENT.
C         --------------------------------------------
      LMAX=11
      IF(XDELLA.LT.0.5) LMAX=22
      SIGOBS=0.5
      SIGMOD=0.5
C
      LOIP=.FALSE.
      INQUIRE(FILE='oip',EXIST=LOIP)
      IF (LOIP) THEN 
C
        OPEN(UNIT=92,FILE='oip')
        READ (92,OIP)
        CLOSE(UNIT=92)
C
        IF ( LMAX .GT. 22 .OR. SIGOBS .GT. 1.0 .OR. SIGMOD .GT. 1.0
     .       .OR. LMAX .LT. 1 .OR. SIGOBS .LT. 0.1 .OR. SIGMOD .LT.
     .       0.1 ) THEN
          WRITE(IU06,*) '*******************************************'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '*    FATAL ERROR IN SUB. WAMASSI          *'
          WRITE(IU06,*) '*    ===========================          *'
          WRITE(IU06,*) '*                                         *'
          WRITE(IU06,*) '* CONSTRAINTS:                            *'
          WRITE(IU06,*) '* 0 < LMAX < 23        CURRENT SETTING: ',
     .                     LMAX,'*'
          WRITE(IU06,*) '* 0.1 < SIGOBS <= 1.0  CURRENT SETTING: ',
     .                     SIGOBS,'*'
          WRITE(IU06,*) '* 0.1 < SIGMOD <= 1.0  CURRENT SETTING: ',
     .                     SIGMOD,'*'
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
C
C*    2.2 OI.
C         ---
C
      CALL OIFIELD (USME, USTARG, USOI, SIGOBS, SIGMOD, 
     .              LMAX, NDIM1, NDIM2)
C
      IF (ITEST.GE.2)
     1   WRITE(IU06,*) '   SUB. WAMASSI: OIFIELD DONE FOR WINDS'
C
      CALL OIFIELD (WHME, WHGTTG, WHOI, SIGOBS, SIGMOD, 
     .              LMAX, NDIM1, NDIM2)
C
      IF (ITEST.GE.2)
     1   WRITE(IU06,*) '   SUB. WAMASSI: OIFIELD DONE FOR WAVE HEIGHTS'
C
CHG   TITL = ' OIFIELD:    WAVE HEIGHT ( 0.1 METRES )'                  
CHG   CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 10.,          
CHG  1            WHOI, AMOWEP, AMOSOP, AMOEAP, AMONOP)                 
CHG   TITL = ' OIFIELD: USTAR    ( 0.01 METRES/SECOND )'                
CHG   CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 100.,         
CHG  1            USOI, AMOWEP, AMOSOP, AMOEAP, AMONOP)                 
C                                                                       
C*    3. ANALYSING THE SPECTRA.                                         
C        ----------------------                                         
C                                                                       
 3000 CONTINUE                                                          
C                                                                       
C*    3.1 UNITS ARE PREPARED.                                           
C         -------------------                                           
C                                                                       
      IF(.NOT.LMESSPASS) REWIND IU12
      IF (IGL.NE.1) THEN                                                
         REWIND IU11                                                    
         REWIND IU12                                                    
         REWIND IU13                                                    
         REWIND IU14                                                    
C                                                                       
C*    3.2 CHANGE WORKING FILE UNITS.                                    
C         --------------------------                                    
C                                                                       
         IFIL=IU11                                                      
         IU11=IU12                                                      
         IU12=IFIL                                                      
         IFIL=IU13                                                      
         IU13=IU14                                                      
         IU14=IFIL                                                      
C                                                                       
C*    3.3 THE FIRST BLOCK IS READ IN.                                   
C         ---------------------------                                   
C                                                                       
         READ(IU11) FL1             
         IF (ITEST.GE.2)                                                
     1      WRITE(IU06,*) '   SUB. WAMASSI: START OF READING',         
     2                    ' FIRST BLOCK OF SPECTRA'                     
      ENDIF                                                             
C                                                                       
C*    3.4 LOOP OVER THE BLOCKS.                                         
C         ---------------------                                         
C                                                                       
 3400 CONTINUE                                                          
      DO 3401 IG=1,IGL                                                  
C                                                                       
C*    3.4.1 TRANSFERS RESULTS OF OPTIMAL INT TO BLOCKS.                 
C           -------------------------------------------                 
C                                                                       
         CALL MAKEBLO (USOIB(IJS(IG)), USOI, IJS(IG), IJL(IG), IG)
         CALL MAKEBLO (WHOIB(IJS(IG)), WHOI, IJS(IG), IJL(IG), IG)
C
         DO 3411 IJ=IJS(IG),IJL(IG)                                     
            USNEW(IJ) = USOLD(IJ,IG)                                    
            THWNEW(IJ)   = THWOLD(IJ,IG)                                
 3411    CONTINUE                                                       
         DO 3412 IJ=IJS(IG),IJL(IG)                                     
            ETOIB(IJ) = -99.                                            
            IF (WHOIB(IJ).GT.0.) ETOIB(IJ) = WHOIB(IJ)**2/16.           
 3412    CONTINUE                                                       
C
         IF (ITEST.GE.2) THEN                                           
            IF (IG.LE.ITESTB) WRITE(IU06,*)                             
     1      '   SUB. WAMASSI: FIELDS BLOCKED FOR IG = ', IG             
         ENDIF                                                          
C                                                                       
C*    3.4.2 CHECK UNITS.                                                
C           ------------                                                
C                                                                       
         IF (MOD(IG,2).EQ.1) THEN                                       
C                                                                       
C*    3.4.2.1  ODD  BLOCKS.                                             
C              ------------                                             
           IF (LMESSPASS) THEN
C
            CALL UPDATE (FL1, ETOIB, USOIB, USNEW, THWNEW,
     1                   IJS(IG), IJL(IG))
C
            IF (ITEST.GE.2) THEN                                        
               IF (IG.LE.ITESTB) WRITE(IU06,*)                          
     1         '   SUB. WAMASSI: SPECTRA UPDATED '   
            ENDIF
C
           ELSE
C
C*             START READING FOR NEXT BLOCK.                           
C                                                                       
            IF (IG.NE.IGL)                                              
     1         READ(IU11) FL2       
C                                                                       
C              UPDATE SPECTRA AND USTAR.                                
C              (PARALLEL PROCESSING - HENCE MT_ALLOT).
C
               NCHUNK = 0
               NLENGTH = (IJL(IG)-IJS(IG)+1)
               NWORK   = NLENGTH/NSIZE
               IF (MOD(NLENGTH,NSIZE).NE.0) NWORK = NWORK + 1
C
               NTSKCAP = 3
               NSTKSZ  = 1
               NSTKINC = 1
               DO 3423 JROC=1,NTASKS-1
                  ITSKCA(1,JROC) = NTSKCAP
                  ITSKCA(3,JROC) = JROC
                  ITSKCA(4,JROC) = 1 !  JUST TO GET PAST COMP 'STACKSZW'
                  ITSKCA(5,JROC) = NSTKSZ
                  ITSKCA(6,JROC) = 1 !  JUST TO GET PAST COMP 'STACKINW'
                  ITSKCA(7,JROC) = NSTKINC
C
                  CALL TSKSTART(ITSKCA(1,JROC), MT_ALLOT, FL1, ETOIB,
     1                          USOIB, USNEW, THWNEW, IJS(IG), IJL(IG),
     2                          NSIZE,NWORK)
 3423          CONTINUE
C
               CALL MT_ALLOT (FL1, ETOIB, USOIB, USNEW, THWNEW,
     1                        IJS(IG), IJL(IG), NSIZE,NWORK)
C
               DO 3424 JROC=1,NTASKS-1
                  CALL TSKWAIT(ITSKCA(1,JROC))
 3424          CONTINUE
C
C
            IF (ITEST.GE.2) THEN                                        
               IF (IG.LE.ITESTB) WRITE(IU06,*)                          
     1         '   SUB. WAMASSI: SPECTRA UPDATED FOR BLOCK IG = ', IG   
            ENDIF                                                       
C                                                                       
C              OUTPUT OF OVERLAPPING LATITUDES AND                      
C              START OF OUTPUT OF ALL SPECTRA IN BLOCK.                 
C                                                                       
            IF (IGL.NE.1) THEN                                          
               IF (IG.NE.1) CALL SPLITBL(FL1, IJS(IG), IJL2(IG), IU14)  
               WRITE(IU12) FL1      
               IF (ITEST.GE.2) THEN                                     
                  IF (ITESTB.GE.IG) WRITE(IU06,*)                       
     1               '   SUB. WAMASSI: OVERLAPPING LATITUDES',          
     2               ' SAVED AND START OF BLOCK OUTPUT.'                
               ENDIF                                                    
            ENDIF                                                       
           ENDIF
C                                                                       
C           THE WIND FIELD CURRENTLY USED IS UPDATED.                   
C                                                                       
         DO 3431 IJ=IJS(IG),IJL(IG)                                     
            U10OLD(IJ,IG) = USNEW(IJ)*U10OLD(IJ,IG)/USOLD(IJ,IG)        
            U10NEW(IJ) = U10OLD(IJ,IG)                                  
 3431    CONTINUE                                                       
         IF (ITEST.GE.2) THEN                                           
            IF (IG.LE.ITESTB) WRITE(IU06,*)                             
     1      '   SUB. WAMASSI: WIND REPLACED FOR BLOCK IG = ', IG        
         ENDIF                                                          
C                                                                       
C              MODEL OUTPUT INTEGRATED DATA ARE SAVED IN COMMON.        
C              SPECTRA ARE SAVED ON DISK AND PRINTED.                   
C                                                                       
            IF (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO .OR.           
     1          CDTSPT.EQ.CDTPRO .OR. CDTSPS.EQ.CDTPRO ) THEN
               CALL OUTBS (FL1, FL3, IJS(IG), IJL(IG), IG, IGL,         
     1            IU27, IU28,THWOLD,USOLD,U10NEW,THWNEW,USNEW,TAUW,
     2            NSTART,NEND)
               IF (ITEST.GE.2) THEN                                     
                  IF (ITEST.GE.2) WRITE(IU06,*)                         
     1               '   SUB. WAMASSI: MODEL OUTPUT DONE IN SUB OUTBS'
               ENDIF                                                    
            ENDIF                                                       
         ELSE                                                           
C                                                                       
C*    3.4.2.2  EVEN BLOCKS.                                             
C              ------------                                             
C                                                                       
C              START READING FOR NEXT BLOCK.                           
C                                                                       
            IF (IG.NE.IGL)                                              
     1         READ(IU11) FL1       
C                                                                       
C              UPDATE SPECTRA AND USTAR.                                
C                                                                       
C              (PARALLEL PROCESSING - HENCE MT_ALLOT).
C
C
               NCHUNK = 0
               NLENGTH = (IJL(IG)-IJS(IG)+1)
               NWORK   = NLENGTH/NSIZE
               IF (MOD(NLENGTH,NSIZE).NE.0) NWORK = NWORK + 1
C
               NTSKCAP = 3
               NSTKSZ  = 1
               NSTKINC = 1
               DO 3432 JROC=1,NTASKS-1
                  ITSKCA(1,JROC) = NTSKCAP
                  ITSKCA(3,JROC) = JROC
                  ITSKCA(4,JROC) = 1 !  JUST TO GET PAST COMP 'STACKSZW'
                  ITSKCA(5,JROC) = NSTKSZ
                  ITSKCA(6,JROC) = 1 !  JUST TO GET PAST COMP 'STACKINW'
                  ITSKCA(7,JROC) = NSTKINC
C
                  CALL TSKSTART(ITSKCA(1,JROC), MT_ALLOT, FL2, ETOIB,
     1                          USOIB, USNEW, THWNEW, IJS(IG), IJL(IG),
     2                          NSIZE,NWORK)
 3432          CONTINUE
C
               CALL MT_ALLOT (FL2, ETOIB, USOIB, USNEW, THWNEW,
     1                        IJS(IG), IJL(IG), NSIZE, NWORK)
C
               DO 3433 JROC=1,NTASKS-1
                  CALL TSKWAIT(ITSKCA(1,JROC))
 3433          CONTINUE
C
            IF (ITEST.GE.2) THEN                                        
               IF (IG.LE.ITESTB) WRITE(IU06,*)                          
     1         '   SUB. WAMASSI: SPECTRA UPDATED FOR BLOCK IG = ', IG   
            ENDIF                                                       
C                                                                       
C           THE WIND FIELD CURRENTLY USED IS UPDATED.                   
C                                                                       
         DO 3434 IJ=IJS(IG),IJL(IG)                                     
            U10OLD(IJ,IG) = USNEW(IJ)*U10OLD(IJ,IG)/USOLD(IJ,IG)        
            U10NEW(IJ) = U10OLD(IJ,IG)                                  
 3434    CONTINUE                                                       
         IF (ITEST.GE.2) THEN                                           
            IF (IG.LE.ITESTB) WRITE(IU06,*)                             
     1      '   SUB. WAMASSI: WIND REPLACED FOR BLOCK IG = ', IG        
         ENDIF                                                          
C                                                                       
C              OUTPUT OF OVERLAPPING LATITUDES AND                      
C              START OF OUTPUT OF ALL SPECTRA IN BLOCK.                 
C                                                                       
            IF (IGL.NE.1) THEN                                          
               IF (IG.NE.1) CALL SPLITBL(FL2, IJS(IG), IJL2(IG), IU14)  
               WRITE(IU12) FL2      
               IF (ITEST.GE.2) THEN                                     
                  IF (ITESTB.GE.IG) WRITE(IU06,*)                       
     1               '   SUB. WAMASSI: OVERLAPPING LATITUDES',          
     2               ' SAVED AND START OF BLOCK OUTPUT.'                
               ENDIF                                                    
            ENDIF                                                       
C                                                                       
C              MODEL OUTPUT INTEGRATED DATA ARE SAVED IN COMMON.        
C              SPECTRA ARE SAVED ON DISK AND PRINTED.                   
C                                                                       
            IF (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO .OR.           
     1          CDTSPT.EQ.CDTPRO .OR. CDTSPS.EQ.CDTPRO ) THEN
               CALL OUTBS (FL2, FL3, IJS(IG), IJL(IG), IG, IGL,         
     1            IU27, IU28,THWOLD,USOLD,U10NEW,THWNEW,USNEW,TAUW,
     2            NSTART,NEND)
C                                                                       
               IF (ITEST.GE.2) THEN                                     
                  IF (ITESTB.GE.IG) WRITE(IU06,*)                       
     1               '   SUB. WAMASSI: MODEL OUTPUT DONE IN SUB OUTBS'  
               ENDIF                                                    
            ENDIF                                                       
         ENDIF                                                          
C                                                                       
C*    BRANCHING BACK TO 3.4 FOR NEXT BLOCK OF LATITUDES                 
C                                                                       
      IF(LWCOU) THEN
        CALL OUTBETA ( IJS(IG), IJL(IG), IG, IGL,USNEW,Z0NEW,
     1                 NSTART, NEND, BETAG, NLONW, NLATW)
      ENDIF
C
C
 3401 CONTINUE                                                          
C                                                                       
C     FREE UP SOME MEMORY
C
      DEALLOCATE(USOI,USME,WHOI,WHME,WHOIB,USOIB,ETOIB)
C
C*    3.5 SAVE FILES.                                                   
C         -----------                                                   
C                                                                       
 3500 CONTINUE                                                          
C     3.5.1 WRITE INTEGRATED DATA TO FILE AND/OR PRINTER                
C           DATA WERE COLLECTED INSIDE THE BLOCK LOOP.                  
C           --------------------------------------------                
C                                                                       
      IF (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO) THEN                
C
C       PRINTS OUT NORMS
        CALL OUTWNORM
C
        CALL OUTINT (IU06, IU22, IU23, IU32, ITEST, MFDBSF)
C
        IF(LMESSPASS) THEN
          CALL MPE_BARRIER(IERR)
          IF(IERR.LT.0) CALL MPEI_ABORT('MPE_BARRIER ERROR in WAMASSI')
        ENDIF
         IF (ITEST.GE.2) WRITE(IU06,*) '   SUB. WAMASSI: OUTINT DONE'
      ENDIF
C                                                                       
      IF (CDATEE.EQ.CDTPRO.OR.LDWRRE) THEN
C                                                                       
C*    3.5.2 STORE 2-D SPECTRA AT SELECTED GRID POINTS,                  
C           INTEGRATED PARAMETERS OF ENTIRE GRID.                       
C           ------------------------------------------                  
C                                                                       
C           SAVE INTEGRATED DATA FILE.                                  
C                                                                       
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
C                                                                       
C           SAVE INTEGRATED SWELL AND SEA DATA FILE.                    
C           ----------------------------------------                    
C                                                                       
         IF (FFLAG21) THEN
            IF((LMESSPASS.AND.(IRANK.EQ.MIN(2,NPROC)))
     $         .OR.(.NOT.LMESSPASS)) THEN
            CALL GSFILE (IU06, IU23, 0, CDTPRO, CDATEF, 'ASS', 'S')
            ENDIF
          ENDIF
C                                                                       
C           SAVE SPECTRA FILE.                                          
C                                                                       
         IF (FFLAG25) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            CALL GSFILE (IU06, IU27, 0, CDTPRO, CDATEF, 'AUT', 'S')
          ENDIF
         ENDIF
C                                                                       
C           SAVE SWELL SPECTRA FILE.                                    
C                                                                       
         IF (FFLAG26) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            CALL GSFILE (IU06, IU28, 0, CDTPRO, CDATEF, 'ASW', 'S')
          ENDIF
         ENDIF
      ENDIF
C
      IF (CDTRES.EQ.CDTPRO.OR.CDATEE.EQ.CDTPRO.OR.CDTPRO.EQ.CDATER.OR.
     &    LDWRRE) THEN
C                                                                       
C                                                                       
C*    3.5.4 SAVE RESTART FIELDS.                                        
C           --------------------                                        
C                                                                       
        IF (  ( (IREST.GE.1 .AND. IREST.LE.4)
     1          .AND.
     2          (CDTPRO.EQ.CDATER .OR. CDTPRO.LE.CDATES) )
     3        .OR.
     4        LDWRRE  ) THEN
C
C*          STRESS FIELD IS DUMPED FOR RESTART.
C           -----------------------------------
C           NOTE :
C           STRESS FIELD IS DUMPED FOR RESTART
C           ONLY IF CURRENT TIME CORRESPONDS TO CDATER

          IF((CDTPRO .EQ. CDATER  ).OR.
     &       (CDATER .EQ. '000000000000' ).OR.
     &        LDWRRE) THEN
C
            IF ( .NOT.LGRIBIN ) THEN
              CALL SAVSTRESS(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,NSTART,NEND,
     1                       LDSTOP, LDWRRE, IREST,CDTPRO,CDATEF)
            ENDIF
C
          ENDIF
C
C           SAVE SPECTRUM
C           -------------
C
          IF ( LGRIBIN ) THEN
            CALL SPEC2FDB (IU06, NPROC, IRANK, LWCOU, TH, FR,
     1                     FL1, NINF, NSUP, NANG, NFRE, NGX, NGY,IRGG,
     2                     NLONRGG, KXLT, IXLG, NIBLO, NBLO, ZMISS,
     3                     AMONOP, AMOWEP, AMOEAP, AMOSOP, XDELLA,
     4                     XDELLO, NSTART, NEND,
     5                     CDTPRO, CDATEF, CDATER, CDATES, CDATEA,
     6                     MARSTYPE, YCLASS, YEXPVER, CFDBSF,
     7                     NENSFNB, NTOTENS, LFDBIOOUT, ITEST, MFDBSF)
            ISTAT = IFLUSHFDB (MFDBSF)
            LLFLUSH = .FALSE.
          ELSE
            CALL SAVSPEC(FL1,NSTART,NEND,LDSTOP,LDWRRE,IREST,
     1                   CDTPRO,CDATEF,CDATER,CDATES,MARSTYPE,
     2                   YCLASS, YEXPVER, CFDB2DSP)
          ENDIF
          WRITE(IU06,*) ' '
          WRITE(IU06,*) '  WAVE SPECTRA FILES DISPOSED AT........',
     1                  ' CDTPRO  = ', CDTPRO
          WRITE(IU06,*) ' '
C
C
C*    3.5.5 UPDATE, WRITE AND SAVE WAMINFO FILE.
C           ------------------------------------
C
          IF (LDSTOP .AND. 
     1        (LMESSPASS .AND. IRANK .EQ. 1 .OR. .NOT. LMESSPASS )) THEN
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
C
            OPEN (IU04,FILE=CWI(1:ICPLEN+8),
     1            STATUS='UNKNOWN',FORM='FORMATTED')
            CALL WRITSTA (IU04, CDTPRO, CDATEE, IANALPD, IFOREPD,
     1                    IDELWIN, CDATER, CDATES, CBPLTDT, CEPLTDT,
     2                    IASSI, NFCST, ISTAT)
C
            CLOSE (IU04)
            WRITE(IU06,*) ' WAMINFO FILE WRITTEN FOR RESTART...',
     1                      ' CDTPRO  = ', CDTPRO
            WRITE(IU06,*) '                                    ',
     1                      ' CDATEF  = ', CDATEF
            WRITE(IU06,*) ' TO ', CWI(1:ICPLEN+8)
          ENDIF
        ENDIF
C
        IF(LMESSPASS) THEN
          CALL MPE_BARRIER(IERR)
          IF(IERR.LT.0) CALL MPEI_ABORT('MPE_BARRIER ERROR IN  WAMASSI')
        ENDIF
      END IF                                                            
      IF (LFDB .AND. LLFLUSH .AND. (
     .  (CDTINTT.EQ.CDTPRO. OR. CDTINTS.EQ.CDTPRO))) THEN
        ISTAT = IFLUSHFDB (MFDBSF)
        WRITE(IU06,*) ' DB ', MFDBSF , ' FLUSHED AT ',
     .    CDTPRO, ' FROM WAMASSI. '
        CALL FLUSH (IU06)
      ENDIF
C
      RETURN                                                            
      END                                                               
