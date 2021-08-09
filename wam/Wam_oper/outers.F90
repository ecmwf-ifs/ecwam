      SUBROUTINE OUTERS (FL1, IJS, IJL, CDTPRO)
! -----------------------------------------------------------------     

!**** *OUTERS* -  OUTPUT OF SATELLITE COLOCATION SPECTRA.

!     H. GUNTHER     GKSS/ECMWF            JULY 1991                    
!     J. BIDLOT      ECMWF           FEBRARY 1996  MESSAGE PASSING
!     J. BIDLOT      ECMWF           JANUARY 1997  MODIFY DEFINITION
!                                                  OF TH0 

!*** INTERFACE.                                                         
!    ----------                                                         

!       *CALL  OUTERS (FL1, CDTPRO)
!          *FL1*    -   SPECTRUM.                                       
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

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOER  , ONLY : NERS     ,CDTERS   ,IERS     ,IJERS    ,    &
     &            IGERS
      USE YOWFRED  , ONLY : FR       ,TH       ,FRATIO
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,AMOWEP   ,AMOSOP   ,    &
     &            XDELLA   ,ZDELLO
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NIBLO
      USE YOWPCONS , ONLY : DEG      ,ZMISS
      USE YOWSPEC, ONLY   : NSTART   ,NEND     ,                        &
     &            U10NEW   ,THWNEW   ,USNEW
      USE YOWTEST  , ONLY : IU06

      USE MPL_MODULE

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "ersfile.intfb.h"
#include "mpgatherersfile.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: FL1
      CHARACTER(LEN=14), INTENT(IN) :: CDTPRO

      INTEGER(KIND=JWIM) :: NGOU, KCOUNT, IJ, K, M 
      INTEGER(KIND=JWIM) :: IU91, IU92, ITAG
      INTEGER(KIND=JWIM) :: KBUFRLENGTH
      INTEGER(KIND=JWIM) :: IRECV, NSPFLD, NSCFLD
      INTEGER(KIND=JWIM),ALLOCATABLE :: IBUFF(:)

      REAL(KIND=JWRB) :: XANG, XFRE, XLAT, XLON
      REAL(KIND=JWRB),ALLOCATABLE :: FLPTS(:,:,:,:)
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: EM, FM, THQ
      REAL(KIND=JWRB),ALLOCATABLE,DIMENSION(:) :: EMPTS, FMPTS, THQPTS, &
     &                                U10PTS, THWPTS, USPTS

! ----------------------------------------------------------------------

      IU91 = 91                                                         
      IU92 = 92                                                         

        CALL ERSFILE (IU91, IU92, .false.)

        IF (NPROC > 1 .AND. CDTPRO == CDTERS) THEN

!*  EXCHANGE THE VALUE OF NERS, IERS,IGERS, AND IJERS WITH THE OTHER PE

          ITAG=92
          ALLOCATE(IBUFF(2))

          IF (IRANK == 1) THEN
            IBUFF(1) = NERS
            IBUFF(2) = IERS
          ENDIF
          CALL MPL_BROADCAST(IBUFF,KROOT=1,KTAG=ITAG,CDSTRING='OUTERS:')
          ITAG=ITAG+1
          IF (IRANK /= 1) THEN
            NERS = IBUFF(1)
            IERS = IBUFF(2)
          ENDIF
          DEALLOCATE(IBUFF)

          IF (IERS > 0.) THEN
            KBUFRLENGTH=2*NERS
            ALLOCATE(IBUFF(KBUFRLENGTH))

            IF (IRANK == 1) THEN
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

            CALL MPL_BROADCAST(IBUFF(1:KBUFRLENGTH),KROOT=1,KTAG=ITAG,  &
     &       CDSTRING='OUTERS 1:')

            IF (IRANK /= 1) THEN
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

      IF (CDTPRO == CDTERS .AND. IERS > 0) THEN

!       COMPUTE MEAN PARAMETERS
        CALL FEMEAN (IJS, IJL, FL1, EM, FM)
        CALL STHQ (IJS, IJL, FL1, THQ)

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

          CALL MPGATHERERSFILE(IRECV,ITAG,NSTART,NEND,NSPFLD,   &
     &                         NSCFLD,                          &
     &                         IJS, IJL, FL1, FLPTS,            &
     &                         EM, FM, THQ,                     &
     &                         U10NEW, THWNEW, USNEW,           &
     &                         EMPTS, FMPTS, THQPTS,            &
     &                         U10PTS, THWPTS, USPTS)

        IF (IRANK == 1) THEN
          XANG = REAL(NANG)
          XFRE = REAL(NFRE)
          DO NGOU=1,IERS
              IJ = IJERS(NGOU)                                         
              XLON = AMOWEP + REAL(IXLG(IJ)-1)*ZDELLO(KXLT(IJ))
              XLAT = AMOSOP + REAL(KXLT(IJ)-1)*XDELLA               

!*    1.1 WRITE INFORMATION TO FILE IU92.                               
!         -------------------------------                               

              WRITE(IU92) XLON, XLAT, CDTPRO, XANG, XFRE,               &
     &                    TH(1), FR(1), FRATIO
              IF (EMPTS(NGOU) > 0.0 ) THEN
                WRITE(IU92) 4.0_JWRB*SQRT(EMPTS(NGOU)), DEG*THQPTS(NGOU), &
     &                    FMPTS(NGOU), USPTS(NGOU), DEG*THWPTS(NGOU),     &
     &                    U10PTS(NGOU)
              ELSE
                WRITE(IU92) ZMISS, ZMISS, ZMISS, ZMISS, ZMISS, ZMISS
              ENDIF
              WRITE(IU92) ((FLPTS(1,NGOU,K,M),K=1,NANG),M=1,NFRE)
          ENDDO
        ENDIF
      ENDIF                                                             

      IF (CDTPRO == CDTERS) THEN
         CALL ERSFILE (IU91, IU92, .true.)
      ENDIF

      IF (ALLOCATED(IJERS))  DEALLOCATE (IJERS)
      IF (ALLOCATED(IGERS))  DEALLOCATE (IGERS)
      IF (ALLOCATED(EMPTS)) DEALLOCATE (EMPTS)
      IF (ALLOCATED(FMPTS)) DEALLOCATE (FMPTS)
      IF (ALLOCATED(THQPTS)) DEALLOCATE (THQPTS)
      IF (ALLOCATED(U10PTS)) DEALLOCATE (U10PTS)
      IF (ALLOCATED(THWPTS)) DEALLOCATE (THWPTS)
      IF (ALLOCATED(USPTS)) DEALLOCATE (USPTS)
      IF (ALLOCATED (FLPTS)) DEALLOCATE (FLPTS)

      END SUBROUTINE OUTERS
