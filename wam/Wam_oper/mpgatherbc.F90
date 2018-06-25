      SUBROUTINE MPGATHERBC(IRECV, IJS, IJL, NSCFLD,                    &
     &                      FL, EMEAN, THQ, FMEAN,                      &
     &                      FLPTS, EMPTS, TQPTS, FMPTS)

!***  *MPGATHERBC* - MESSAGE PASSING OF COARSE BOUNDARY CONDITIONS.


!     PURPOSE                                                       
!     -------                                                      

!     *MPGATHERBC* - GATHERS TO PE IRECV ALL COARSE BOUNDARY CONDITIONS

!      *IRECV*   - PE RECEIVING ALL CONTRIBUTIONS.
!      *IJS*     - INDEX OF FIRST GRIDPOINT.
!      *IJL*     - INDEX OF LAST GRIDPOINT.
!      *NSCFLD*  - NUMBER OF SCALAR FIELDS TO BE GATHERED
!      *FL*     - BLOCK OF SPECTRA.
!      *EMEAN*   - BLOCK TOTAL ENERGY
!      *THQ*     - BLOCK MEAN DIRECTION
!      *FMEAN*   - BLOCK MEAn FREQUENCY
!       THE FOLLOWING ARRAY WILL BE GATHERED TO COARSE BOUNDARY POINTS.
!                  FL   --> FLPTS
!                  EMEAN --> EMPTS
!                  THQ   --> TQPTS
!                  FMEAN --> FMPTS

!----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCPBO  , ONLY : NBOUNC   ,IJARC
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,LL1D
      USE YOWSPEC  , ONLY : IJ2NEWIJ ,NSTART   ,NEND
      USE YOWTEST  , ONLY : IU06
      USE MPL_MODULE

!----------------------------------------------------------------------
      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM) :: IRECV, ITAG, IJS, IJL, NSCFLD
      INTEGER(KIND=JWIM) :: MAXLENGTH, MESLENGTH
      INTEGER(KIND=JWIM) :: NGOU, IJ, M, K, IP, KCOUNT, IST, IND, IGOU
      INTEGER(KIND=JWIM) :: ISENDREQ
      INTEGER(KIND=JWIM) :: KRFROM, KRCOUNT, KRTAG
      INTEGER(KIND=JWIM), DIMENSION(NPROC) :: NPTS, IBS, IBL
      INTEGER(KIND=JWIM), DIMENSION(NBOUNC) :: IJBC, NGOUG

      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: EMEAN, FMEAN, THQ
      REAL(KIND=JWRB), DIMENSION(NBOUNC) :: EMPTS, TQPTS, FMPTS
      REAL(KIND=JWRB), DIMENSION(NBOUNC,NANG,NFRE) :: FLPTS
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NANG,NFRE) :: FL
      REAL(KIND=JWRB), ALLOCATABLE :: ZCOMBUFS(:), ZCOMBUFR(:)

!----------------------------------------------------------------------

      ITAG=1

!     1.0 DEFAULT ACTION IF NO FIELD GATHERING 
!         ------------------------------------
      IF (IRECV.EQ.0 .OR. NPROC.EQ.1) THEN
        DO NGOU=1,NBOUNC
          IJ=IJARC(NGOU)
          EMPTS(NGOU)=EMEAN(IJ)
          TQPTS(NGOU)=THQ(IJ)
          FMPTS(NGOU)=FMEAN(IJ)
          DO  M = 1, NFRE
            DO K = 1, NANG
              FLPTS(NGOU,K,M)=FL(IJ,K,M)
            ENDDO
          ENDDO
        ENDDO

      ELSE

!       DETERMINE WHICH PE'S NEED TO COMMUNICATE
!   should be computed only once and moved to a module
        KCOUNT=0
        DO IP=1,NPROC
          NPTS(IP)=0
          IBS(IP)=0
          IBL(IP)=0
          DO NGOU=1,NBOUNC
            IF (LL1D) THEN
              IJ = IJARC(NGOU)
            ELSE
              IJ = IJ2NEWIJ(IJARC(NGOU))
            ENDIF
            IF (IJ.GE.NSTART(IP) .AND. IJ.LE.NEND(IP)) THEN
              NPTS(IP)=NPTS(IP)+1
              KCOUNT=KCOUNT+1
              IJBC(KCOUNT)=IJ
              NGOUG(KCOUNT)=NGOU
            ENDIF
            IF (NPTS(IP).EQ.1) IBS(IP)=KCOUNT
          ENDDO
          IF (NPTS(IP).GT.0) THEN
            IBL(IP)=IBS(IP)+NPTS(IP)-1
          ENDIF
        ENDDO

        MESLENGTH=NSCFLD+NANG*NFRE
        MAXLENGTH=MESLENGTH*NBOUNC
        ALLOCATE(ZCOMBUFS(MAXLENGTH))

!     1.1 SEND TO THE PROCESS THAT GATHERS THE WHOLE FIELD
!         ------------------------------------------------

        IF (NPTS(IRANK).GT.0) THEN
          KCOUNT=MESLENGTH*(IBS(IRANK)-1)
          IST=MESLENGTH*(IBS(IRANK)-1)+1
          IND=IST+NPTS(IRANK)*MESLENGTH-1
          DO IGOU=IBS(IRANK),IBL(IRANK)
            IJ=IJBC(IGOU)
            KCOUNT=KCOUNT+1
            ZCOMBUFS(KCOUNT)=EMEAN(IJ)
            KCOUNT=KCOUNT+1
            ZCOMBUFS(KCOUNT)=THQ(IJ)
            KCOUNT=KCOUNT+1
            ZCOMBUFS(KCOUNT)=FMEAN(IJ)
            DO  M = 1, NFRE
              DO K = 1, NANG
                KCOUNT=KCOUNT+1
                ZCOMBUFS(KCOUNT)=FL(IJ,K,M)
              ENDDO
            ENDDO
          ENDDO
          IF (KCOUNT.NE.IND) THEN
            WRITE (IU06,*) ' ******************************************'
            WRITE (IU06,*) ' *                                        *'
            WRITE (IU06,*) ' *      FATAL ERROR SUB. MPGATHERBC       *'
            WRITE (IU06,*) ' *      ==========================        *'
            WRITE (IU06,*) ' *      KCOUNT SHOULD BE EQUAL TO IND     *'
            WRITE (IU06,*) ' *      KCOUNT = ',KCOUNT
            WRITE (IU06,*) ' *      IND = ',IND
            WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS       *'
            WRITE (IU06,*) ' *                                        *'
            WRITE (IU06,*) ' ******************************************'
            CALL ABORT1
          ENDIF

!         SEND BUFFER
          CALL GSTATS(673,0)
          CALL MPL_SEND(ZCOMBUFS(IST:IND),KDEST=IRECV,KTAG=ITAG,        &
     &                  KMP_TYPE=JP_NON_BLOCKING_STANDARD,              &
     &                  KREQUEST=ISENDREQ,                              &
     &                  CDSTRING='MPGATHERBC:')
          CALL GSTATS(673,1)
        ENDIF

!       RECEIVE CONTRIBUTIONS
!       ---------------------
        IF (IRANK.EQ.IRECV) THEN
          ALLOCATE(ZCOMBUFR(MAXLENGTH))
          CALL GSTATS(673,0)
          DO IP=1,NPROC
            IF (NPTS(IP).GT.0) THEN
              IST=MESLENGTH*(IBS(IP)-1)+1
              IND=IST+NPTS(IP)*MESLENGTH-1
              KCOUNT=IND-IST+1
              CALL MPL_RECV(ZCOMBUFR(IST:IND),                          &
     &                      KSOURCE=IP, KTAG=ITAG,                      &
     &                      KOUNT=KRCOUNT, KRECVTAG=KRTAG,              &
     &                      KMP_TYPE=JP_BLOCKING_STANDARD,              &
     &                      CDSTRING='MPGATHERBC:')
              IF (KRTAG.NE.ITAG) CALL MPL_ABORT                         &
     &            ('MPL_RECV ERROR in MPGATHERBC:  MISMATCHED TAGS' )
              IF (KRCOUNT.NE.KCOUNT) CALL MPL_ABORT                     &
     &            ('MPL_RECV ERROR in MPGATHERBC: WRONG MESSAGE LENGTH')
            ENDIF
          ENDDO
          CALL GSTATS(673,1)

!         DECODE BUFFER
          KCOUNT=0
          DO IGOU=1,NBOUNC
            NGOU=NGOUG(IGOU)
            KCOUNT=KCOUNT+1
            EMPTS(NGOU)=ZCOMBUFR(KCOUNT)
            KCOUNT=KCOUNT+1
            TQPTS(NGOU)=ZCOMBUFR(KCOUNT)
            KCOUNT=KCOUNT+1
            FMPTS(NGOU)=ZCOMBUFR(KCOUNT)
            DO  M = 1, NFRE
              DO K = 1, NANG
                KCOUNT=KCOUNT+1
                FLPTS(NGOU,K,M)=ZCOMBUFR(KCOUNT)
              ENDDO
            ENDDO
          ENDDO
          DEALLOCATE(ZCOMBUFR)
        ENDIF

        IF (NPTS(IRANK).GT.0) THEN
          CALL GSTATS(673,0)
          CALL MPL_WAIT(KREQUEST=ISENDREQ, CDSTRING='MPGATHERBC:')
          CALL GSTATS(673,1)
        ENDIF

        DEALLOCATE(ZCOMBUFS)
      ENDIF

      CALL MPL_BARRIER(CDSTRING='MPGATHERBC:')

      END SUBROUTINE MPGATHERBC
