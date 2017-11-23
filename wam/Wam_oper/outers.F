      SUBROUTINE OUTERS (FL3, IJS, IJL, CDTPRO)
! -----------------------------------------------------------------     

!**** *OUTERS* -  OUTPUT OF SATELLITE COLOCATION SPECTRA.

!     H. GUNTHER     GKSS/ECMWF            JULY 1991                    
!     J. BIDLOT      ECMWF           FEBRARY 1996  MESSAGE PASSING
!     J. BIDLOT      ECMWF           JANUARY 1997  MODIFY DEFINITION
!                                                  OF TH0 

!*** INTERFACE.                                                         
!    ----------                                                         

!       *CALL  OUTERS (FL3, CDTPRO)
!          *FL3*    -   SPECTRUM.                                       
!          *CDTPRO* -   MODEL PROPAGATION TIME.                         

!    EXTERNALS.                                                         
!    ----------                                                         

!      ERSFILE   - FETCHES COLOCATION INFORMATION DONE BY PREUWA.       

!    METHOD.                                                            
!    -------                                                            

!      NONE.                                                            

!    REFERENCES.                                                        
!    -----------                                                        

!      NONE.                                                            

!---------------------------------------------------------------------- 

      USE YOWCOER  , ONLY : NERS     ,CDTERS   ,IERS     ,IJERS    ,
     &            IGERS
      USE YOWFRED  , ONLY : FR       ,TH       ,FRATIO
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,AMOWEP   ,AMOSOP   ,
     &            XDELLA   ,ZDELLO
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,THQ
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NIBLO
      USE YOWPCONS , ONLY : DEG      ,ZMISS
      USE YOWSPEC, ONLY   : NSTART   ,NEND     , 
     &            U10NEW   ,THWNEW   ,USNEW
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE MPL_MODULE

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "ersfile.intfb.h"
#include "mpgatherersfile.intfb.h"

      INTEGER, INTENT(IN) :: IJS, IJL
      REAL, DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: FL3
      CHARACTER(LEN=14), INTENT(IN) :: CDTPRO

      INTEGER :: IG
      INTEGER :: NGOU, KCOUNT, IJ, K, M 
      INTEGER :: IU91, IU92, ITAG
      INTEGER :: KBUFRLENGTH
      INTEGER :: IRECV, NSPFLD, NSCFLD
      INTEGER,ALLOCATABLE :: IBUFF(:)

      REAL :: XANG, XFRE, XLAT, XLON
      REAL,ALLOCATABLE :: FLPTS(:,:,:,:)
      REAL,ALLOCATABLE,DIMENSION(:) :: EMPTS,FMPTS,THQPTS,U10PTS,THWPTS,
     &                                 USPTS

! ----------------------------------------------------------------------

      IU91 = 91                                                         
      IU92 = 92                                                         

      IG = 1

        CALL ERSFILE (IU91, IU92, .false.)

        IF(NPROC.GT.1.AND.CDTPRO.EQ.CDTERS) THEN

!*  EXCHANGE THE VALUE OF NERS, IERS,IGERS, AND IJERS WITH THE OTHER PE

          ITAG=92
          ALLOCATE(IBUFF(2))

          IF(IRANK.EQ.1) THEN
            IBUFF(1) = NERS
            IBUFF(2) = IERS
          ENDIF
          CALL MPL_BROADCAST(IBUFF,KROOT=1,KTAG=ITAG,CDSTRING='OUTERS:')
          ITAG=ITAG+1
          IF(IRANK.NE.1) THEN
            NERS = IBUFF(1)
            IERS = IBUFF(2)
          ENDIF
          DEALLOCATE(IBUFF)

          IF(IERS.GT.0.) THEN
            KBUFRLENGTH=2*NERS
            ALLOCATE(IBUFF(KBUFRLENGTH))

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
            ENDIF

            CALL MPL_BROADCAST(IBUFF(1:KBUFRLENGTH),KROOT=1,KTAG=ITAG,
     &       CDSTRING='OUTERS 1:')

            IF(IRANK.NE.1) THEN
              ALLOCATE(IGERS(NERS))
              ALLOCATE(IJERS(NERS))

              KCOUNT=0
              DO NGOU=1,NERS
                KCOUNT=KCOUNT+1
                IGERS(NGOU)=IBUFF(KCOUNT)
              ENDDO
              DO NGOU=1,NERS
                KCOUNT=KCOUNT+1
                IJERS(NGOU)=IBUFF(KCOUNT)
              ENDDO
            ENDIF

            DEALLOCATE(IBUFF)
          ENDIF

        ENDIF

!*   1. LOOP OVER OUTPUT POINTS.                                        
!       ------------------------                                        

      IF (CDTPRO.EQ.CDTERS.AND.IERS.GT.0) THEN                          

!       COLLECT NECESSARY FIELDS TO PROCESS 1

          IRECV=1
          ITAG=ITAG+1
          NSPFLD=1
          NSCFLD=6

          ALLOCATE(EMPTS(IERS))
          ALLOCATE(FMPTS(IERS))
          ALLOCATE(THQPTS(IERS))
          ALLOCATE(U10PTS(IERS))
          ALLOCATE(THWPTS(IERS))
          ALLOCATE(USPTS(IERS))

          ALLOCATE(FLPTS(NSPFLD,IERS,NANG,NFRE))

          CALL MPGATHERERSFILE(IRECV,ITAG,NSTART,NEND,NSPFLD,NSCFLD,
     &                         IJS, IJL, FL3,FLPTS,
     &                         U10NEW(IJS),THWNEW(IJS),USNEW(IJS),
     &                         EMPTS,FMPTS,THQPTS,U10PTS,THWPTS,USPTS)

        IF(IRANK.EQ.1) THEN
          XANG = REAL(NANG)
          XFRE = REAL(NFRE)
          DO NGOU=1,IERS
              IJ = IJERS(NGOU)                                         
              XLON = AMOWEP + REAL(IXLG(IJ,IG)-1)*ZDELLO(KXLT(IJ,IG))
              XLAT = AMOSOP + REAL(KXLT(IJ,IG)-1)*XDELLA               

!*    1.1 WRITE INFORMATION TO FILE IU92.                               
!         -------------------------------                               

              IF (ITEST.GE.3) THEN
                WRITE(IU06,*)' SUB : outers  LOOP 1002 ngou=', ngou
                CALL FLUSH (IU06)
              ENDIF
              WRITE(IU92) XLON, XLAT, CDTPRO, XANG, XFRE,               
     &                    TH(1), FR(1), FRATIO
              IF (EMPTS(NGOU) .GT. 0.0 ) THEN
                WRITE(IU92) 4.*SQRT(EMPTS(NGOU)), DEG*THQPTS(NGOU),
     &                    FMPTS(NGOU), USPTS(NGOU), DEG*THWPTS(NGOU),
     &                    U10PTS(NGOU)
              ELSE
                WRITE(IU92) ZMISS, ZMISS, ZMISS, ZMISS, ZMISS, ZMISS
              ENDIF
              WRITE(IU92) ((FLPTS(1,NGOU,K,M),K=1,NANG),M=1,NFRE)
          ENDDO
        ENDIF
      ENDIF                                                             

      IF (CDTPRO.EQ.CDTERS) THEN
         CALL ERSFILE (IU91, IU92, .true.)
      ENDIF

      IF(ALLOCATED(IJERS))  DEALLOCATE (IJERS)
      IF(ALLOCATED(IGERS))  DEALLOCATE (IGERS)
      IF(ALLOCATED(EMPTS)) DEALLOCATE (EMPTS)
      IF(ALLOCATED(FMPTS)) DEALLOCATE (FMPTS)
      IF(ALLOCATED(THQPTS)) DEALLOCATE (THQPTS)
      IF(ALLOCATED(U10PTS)) DEALLOCATE (U10PTS)
      IF(ALLOCATED(THWPTS)) DEALLOCATE (THWPTS)
      IF(ALLOCATED(USPTS)) DEALLOCATE (USPTS)
      IF(ALLOCATED (FLPTS)) DEALLOCATE (FLPTS)

      END SUBROUTINE OUTERS
