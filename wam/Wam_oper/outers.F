      SUBROUTINE OUTERS (FL3, IG, IGL, CDTPRO, NSTART, NEND,
     1                   U10NEW,THWNEW,USNEW)
                                                                        
C -----------------------------------------------------------------     
C                                                                       
C**** *OUTERS* -  OUTPUT OF ERS_1 COLOCATION SPECTRA.                   
C                                                                       
C     H. GUNTHER     GKSS/ECMWF            JULY 1991                    
!     J. BIDLOT      ECMWF           FEBRARY 1996  MESSAGE PASSING
!     J. BIDLOT      ECMWF           JANUARY 1997  MODIFY DEFINITION OF TH0 
C                                                                       
C*** INTERFACE.                                                         
C    ----------                                                         
C                                                                       
C       *CALL  OUTERS (FL3, IG, IGL, CDTPRO, NSTART, NEND,
C                      U10NEW,THWNEW,USNEW)* 
C          *FL3*    -   SPECTRUM.                                       
C          *IG*     -   BLOCK NUMBER.                                   
C          *IGL*    -   NUMBER OF BLOCKS.                               
C          *CDTPRO* -   MODEL PROPAGATION TIME.                         
C          *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
C          *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
C          *U10NEW*    NEW WIND SPEED IN M/S.
C          *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
C                      NOTATION (POINTING ANGLE OF WIND VECTOR,
C                      CLOCKWISE FROM NORTH).
C          *USNEW*     NEW FRICTION VELOCITY IN M/S.
C                                                                       
C    EXTERNALS.                                                         
C    ----------                                                         
C                                                                       
C      ERSFILE   - FETCHES COLOCATION INFORMATION DONE BY PREUWA.       
C                                                                       
C    METHOD.                                                            
C    -------                                                            
C                                                                       
C      NONE.                                                            
C                                                                       
C    REFERENCES.                                                        
C    -----------                                                        
C                                                                       
C      NONE.                                                            
C                                                                       
C---------------------------------------------------------------------- 
C                                                                       
#include "param.h"
C                                                                       
#include "comcoer.h"
C                                                                       
#include "comfred.h"
C                                                                       
#include "commap.h"
C                                                                       
#include "commean.h"
C                                                                       
#include "comtest.h"
C                                                                       
#include "comwind.h"
C                                                                       
#include "commesspass.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C ----------------------------------------------------------------------
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL3
      REAL,DIMENSION(NINF:NSUP) :: U10NEW,THWNEW,USNEW

C ----------------------------------------------------------------------
C                                                                       
#include "parcons.h"
C                                                                       
      PARAMETER ( FREI = 1.1)                                  
C                                                                       
C                FREI   REAL  FREQUENCY INCREMENT                       
C                                                                       
      CHARACTER NAME*20, CDTPRO*10                                                 
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C     WORK ARRAYS
      INTEGER, DIMENSION(2*NERS+1) :: IBUFF
      REAL,ALLOCATABLE,DIMENSION(:) :: U10,THW,US
      REAL,ALLOCATABLE :: FLPTS(:,:,:,:)

      IU91 = 91                                                         
      IU92 = 92                                                         

      IF (IG.EQ.1) THEN
         CALL ERSFILE (IU91, IU92, IG, IGL)

       IF(LMESSPASS) THEN
C      EXCHANGE THE VALUE OF IERS,IGERS, AND IJERS WITH THE OTHER PROCESS
          ITAG=92
        IF(IRANK.EQ.1) THEN
          KCOUNT=0
          DO NGOU=1,NERS
             KCOUNT=KCOUNT+1
             IBUFF(KCOUNT)=IGERS(NGOU)
          ENDDO
          DO NGOU=1,NERS
             KCOUNT=KCOUNT+1
             IBUFF(KCOUNT)=IJERS(NGOU)
          ENDDO
          KCOUNT=KCOUNT+1
          IBUFF(KCOUNT)=IERS
          CALL MPE_BROADCAST(IBUFF,KCOUNT,1,1,ITAG,0,0,0,IERR)
          IF(IERR.LT.0) THEN
            WRITE(IU06,*) 'OUTERS: MPE_BROADCAST RC=',IERR
            CALL MPEI_ABORT('MPE_BROADCAST ERROR in OUTERS')
          ENDIF
        ELSE
          KBUFRLENGTH=2*NERS+1
          CALL MPE_RECV(IBUFF,KBUFRLENGTH,1,1,ITAG,0,0,0,
     1                  KRCOUNT,KRFROM,KRTAG,IERR)
          IF(IERR.LT.0) CALL MPEI_ABORT('MPE_RECV ERROR in OUTERS')
          IF(KRTAG.NE.ITAG) CALL MPEI_ABORT
     1        ('MPE_RECV ERROR in OUTERS:  MISMATCHED TAGS (1)' )
          KCOUNT=0
          DO NGOU=1,NERS
             KCOUNT=KCOUNT+1
             IGERS(NGOU)=IBUFF(KCOUNT)
          ENDDO
          DO NGOU=1,NERS
             KCOUNT=KCOUNT+1
             IJERS(NGOU)=IBUFF(KCOUNT)
          ENDDO
          KCOUNT=KCOUNT+1
          IERS=IBUFF(KCOUNT)
          IF(KRCOUNT.NE.KCOUNT) CALL MPEI_ABORT
     1        ('MPE_RECV ERROR in OUTERS: WRONG MESSAGE LENGTH')
        ENDIF
       ENDIF
      ENDIF
C                                                                       
C*   1. LOOP OVER OUTPUT POINTS.                                        
C       ------------------------                                        
C                                                                       
      IF (CDTPRO.EQ.CDTERS.AND.IERS.GT.0) THEN                          

        ALLOCATE(U10(NIBLO),THW(NIBLO),US(NIBLO))

        IF(LMESSPASS) THEN
C       COLLECT NECESSARY FIELDS TO PROCESS 1
C
          IRECV=1
          ITAG=ITAG+1
          NSPFLD=1
          NSCFLD=6
          ALLOCATE(FLPTS(NSPFLD,IERS,NANG,NFRE))
          DO IJ = NINF,NSUP
            U10(IJ)=U10NEW(IJ)
            THW(IJ)=THWNEW(IJ)
            US(IJ)=USNEW(IJ)
          ENDDO
          CALL MPGATHERERSFILE(IRECV,ITAG,NSTART,NEND,NSPFLD,NSCFLD,
     1                         FL3,FLPTS,
     2                         EMEAN,FMEAN,THQ,U10,THW,US)
        ELSE
          DO IJ = NINF,NSUP
            U10(IJ)=U10NEW(IJ)
            THW(IJ)=THWNEW(IJ)
            US(IJ)=USNEW(IJ)
          ENDDO
        ENDIF

        IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
         XANG = REAL(NANG)                                              
         XFRE = REAL(NFRE)                                              
         DO 1002 NGOU=1,IERS                                            
            IF (IG.EQ.IGERS(NGOU)) THEN                                 
               IJ = IJERS(NGOU)                                         
               XLON = AMOWEP + REAL(IXLG(IJ,IG)-1)*ZDELLO(KXLT(IJ,IG))
               XLAT = AMOSOP + REAL(KXLT(IJ,IG)-1)*XDELLA               
C                                                                       
C*    1.1 WRITE INFORMATION TO FILE IU92.                               
C         -------------------------------                               
C                                                                       
               WRITE(IU92) XLON, XLAT, CDTPRO, XANG, XFRE,               
     1                     TH(1), FR(1), FREI                             
               WRITE(IU92) 4.*SQRT(EMEAN(IJ)), DEG*THQ(IJ),             
     1                     FMEAN(IJ), US(IJ), DEG*THW(IJ),        
     2                     U10(IJ)                                   
               IF(LMESSPASS) THEN
                 WRITE(IU92) ((FLPTS(1,NGOU,K,M),K=1,NANG),M=1,NFRE)
               ELSE
                 WRITE(IU92) ((FL3(IJ,K,M),K=1,NANG),M=1,NFRE)
               ENDIF
            ENDIF                                                       
 1002    CONTINUE                                                       
        ENDIF
      ENDIF                                                             
                                                                        
      IF (IG.EQ.IGL) THEN                                               
       IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
         CALL ERSFILE (IU91, IU92, -9, -9)                              
       ENDIF
      ENDIF                                                             

      IF( ALLOCATED(U10) ) DEALLOCATE(U10)
      IF( ALLOCATED(THW) ) DEALLOCATE(THW)
      IF( ALLOCATED(US) ) DEALLOCATE(US)
      IF( ALLOCATED(FLPTS) ) DEALLOCATE(FLPTS)

      RETURN                                                            
      END                                                               
