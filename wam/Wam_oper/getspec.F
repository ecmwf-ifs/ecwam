      SUBROUTINE GETSPEC(FL,NSTART,NEND,IREAD)
C ----------------------------------------------------------------------
C     J. BIDLOT    ECMWF      SEPTEMBER 1997 
C
C*    PURPOSE.
C     --------
C     READ THE SPECTRA FROM DISK.
C
C**   INTERFACE.
C     ----------
C     *CALL* *GETSPEC(FL,NSTART,NEND,IREAD)
C     *FL*        ARRAY CONTAINING THE SPECTRA CONTRIBUTION ON EACH PE
C     *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
C     *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
C     *IREAD*     PROCESSOR WHICH WILL ACCESS THE FILE ON DISK 
C
C     METHOD.
C     -------
C     USE READFL TO READ IN THE SPECTRA DEPENDING ON THE USE OF THE PBIO
C     SOFTWARE OR NOT. PBIO WILL LIMIT THE SIZE OF THE ARRAY NECESSARY TO 
C     READ THE INPUT SPECTRA. THE READING IS ONLY DONE ON PE 1, THEREFORE
C     THE RELEVANT INFORMATION IS SENT TO THE OTHER PE'S USING MPDISTRIBFL
C
C     EXTERNALS.
C     ----------
C     GETENV
C     GRSTNAME
C     MPDISTRIBFL
C     MPE_BARRIER
C     READFL
C
C     REFERENCE.
C     ----------
C     NONE
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "comcout.h"
C
#include "comstat.h"
C
#include "comtest.h"
C
#include "comtext.h"
C
#include "txttext.h"
C
#include "commesspass.h"
C
#include "commpp.h"

      CHARACTER*9 CMBXSIZE
      CHARACTER*80 FILENAME
      LOGICAL LFRSTIME, LPBOPEN, LPBCLOSE
C
C ----------------------------------------------------------------------
      INTEGER,DIMENSION(NPROC) :: NSTART, NEND
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL
C ----------------------------------------------------------------------
C     define local allocatable arrays
      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: RFL


      DATA MDEL, LFRSTIME / 1, .TRUE./

      LPBOPEN = .TRUE.
      LPBCLOSE = .TRUE.

      IF(LMESSPASS) THEN

C
C DETERMINE MAXIMUM NUMBER OF FREQUENCIES WHICH CAN BE USED TO SEND 
C SPECTRA FROM ONE PE  TO THE OTHERS
C WITHOUT EXCEEDING THE MAXIMUM MAILBOX SIZE 
        IF(LFRSTIME) THEN
          CALL GETENV('VPP_MBX_SIZE',CMBXSIZE)
          READ(CMBXSIZE,'(I9)') MBXSIZE
          MINSIZE=NINT(FLOAT(NPROC-1)*NIBLO*NANG*NPRECI/NPROC)
          IF(MINSIZE.GT.MBXSIZE) THEN
            WRITE (IU06,*) '***************************************'
            WRITE (IU06,*) '*                                     *'
            WRITE (IU06,*) '* MESSAGE MAIL BOX SIZE IS TOO SMALL  *'
            WRITE (IU06,*) '* TO DISTRIBUTE THE RESTART SPECTRUM  *'
            WRITE (IU06,*) '* FILE                                *'
            WRITE (IU06,*) '* THE PRESCRIBED SIZE SIZE IS         *'
            WRITE (IU06,*) '* ', MBXSIZE,' BYTES                  *'
            WRITE (IU06,*) '* PLEASE INCREASE SIZE TO AT LEAST    *'
            WRITE (IU06,*) '* ', MINSIZE, 'BYTES                  *'
            WRITE (IU06,*) '* SEE ENVIRONMENT VARIABLE            *'
            WRITE (IU06,*) '* VPP_MBX_SIZE                        *' 
            WRITE (IU06,*) '*                                     *'
            WRITE (IU06,*) '***************************************'
            CALL ABORT1
          ENDIF

          MDEL=1
CCC          MDEL = NFRE 
CCC          DO WHILE (MINSIZE*MDEL.GE.MBXSIZE.AND.MDEL.GT.1)
CCC             MDEL=MDEL-1
CCC          ENDDO

          WRITE(IU06,*)''
          WRITE(IU06,*)' IN GETSPEC :'
          WRITE(IU06,*)'  MESSAGE MAIL BOX SIZE OF ',MBXSIZE,' BYTES'
          IF(MDEL.EQ.1) THEN
          WRITE(IU06,*)'  ',MDEL,' FREQUENCY IS DISTRIBUTED AT ONE TIME'
          ELSE
          WRITE(IU06,*)'  ',MDEL,
     &                ' FREQUENCIES ARE DISTRIBUTED EACH TIME'
          ENDIF
          WRITE(IU06,*)'  WHICH SHOULD RESULT IN A MAXIMUM USAGE OF'
          WRITE(IU06,*)'  ',MINSIZE*MDEL,
     &                ' BYTES FOR EACH MPDISTRIBFL CALL'
          WRITE(IU06,*)''
       
        ENDIF
       

        ISEND=IREAD

        IF(LPBIOIN) THEN


          IF (IRANK.EQ.ISEND)
     1        CALL GRSTNAME(IU06,CDTPRO,CDATEF,0,'BLS',FILENAME)

C use the ability of the PBIO software to read input by frequency band

          DO MLOOP=1,NFRE,MDEL

             MINF=MLOOP
             MSUP=MIN(MLOOP+MDEL-1,NFRE)


             ALLOCATE(RFL(0:NIBLO,NANG,MINF:MSUP))

C     READ RESTART SPECTRA FROM PE ISEND (IREAD) 
             IF (IRANK.EQ.ISEND) THEN

               LPBOPEN = .FALSE.
               LPBCLOSE = .FALSE.
               IF(MINF.EQ.1) LPBOPEN = .TRUE.
               IF(MSUP.EQ.NFRE) LPBCLOSE = .TRUE.

               CALL READFL(RFL,MINF,MSUP,FILENAME,IUNIT,
     1                     LPBOPEN,LPBCLOSE)
             ENDIF

C     PUT A BARRIER TO AVOID FLOODING THE CROSSBAR 
             CALL MPE_BARRIER(IERR)
             IF(IERR.LT.0) 
     1       CALL MPEI_ABORT('MPE_BARRIER ERROR in MPDISTRIB 1')

             CALL MPDISTRIBFL(ISEND,KTAG,NSTART,NEND,MINF,MSUP,
     1                        MINF,MSUP,RFL)
             KTAG=KTAG+1
             IF (ITEST.GE.2)
     1       WRITE(IU06,*)
     2       'SUB. GETSPEC: RESTART SPECTRUM COLLECTED, PART:',MLOOP

C     KEEP CORRESPONDING CONTRIBUTION TO FL
             DO M=MINF,MSUP
               DO K=1,NANG
                  DO IJ=NSTART(IRANK),NEND(IRANK)
                    FL(IJ,K,M) = RFL(IJ,K,M)
                  END DO
               END DO
             END DO

             DO M=MINF,MSUP
               DO K=1,NANG
                  FL(NINF-1,K,M) = 0. 
               END DO
             END DO

             DEALLOCATE(RFL)

          END DO


        ELSE

           ALLOCATE(RFL(0:NIBLO,NANG,NFRE))

C     READ RESTART SPECTRA FROM PE ISEND (IREAD)
           IF (IRANK.EQ.ISEND) THEN
             CALL READFL(RFL,1,NFRE,FILENAME,IUNIT,LPBOPEN,LPBCLOSE)
           ENDIF

           DO MLOOP=1,NFRE,MDEL
              MINF=MLOOP
              MSUP=MIN(MLOOP+MDEL-1,NFRE)

C     PUT A BARRIER TO AVOID FLOODING THE CROSSBAR 
             CALL MPE_BARRIER(IERR)
             IF(IERR.LT.0) 
     1       CALL MPEI_ABORT('MPE_BARRIER ERROR in MPDISTRIB 2')

             CALL MPDISTRIBFL(ISEND,KTAG,NSTART,NEND,MINF,MSUP,
     1                         1,NFRE,RFL)
             KTAG=KTAG+1
           ENDDO

C     KEEP CORRESPONDING CONTRIBUTION TO FL
             DO M=1,NFRE
               DO K=1,NANG
                  DO IJ=NSTART(IRANK),NEND(IRANK)
                    FL(IJ,K,M) = RFL(IJ,K,M)
                  END DO
               END DO
             END DO

             DO M=1,NFRE
               DO K=1,NANG
                  FL(NINF-1,K,M) = 0. 
               END DO
             END DO

           IF (ITEST.GE.2)
     1     WRITE(IU06,*)
     2     'SUB. GETSPEC: RESTART SPECTRUM COLLECTED'


           DEALLOCATE(RFL)

        ENDIF
C
C
      ELSE

        IF(LPBIOIN) THEN
           CALL GRSTNAME(IU06,CDTPRO,CDATEF,1,'BLS',FILENAME)
           LPBOPEN = .TRUE.
           LPBCLOSE = .TRUE.
        ENDIF

        CALL READFL(FL,1,NFRE,FILENAME,IUNIT,LPBOPEN,LPBCLOSE)

      END IF
C
        WRITE(IU06,*) ' SPECTRUM FILE READ IN........',
     1                ' CDTPRO  = ', CDTPRO
C
      RETURN
      END
