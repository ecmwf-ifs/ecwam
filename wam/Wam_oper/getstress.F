      SUBROUTINE GETSTRESS(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,NSTART,NEND,
     1                     IREAD)
C ----------------------------------------------------------------------
C     J. BIDLOT    ECMWF      SEPTEMBER 1997 
C
C*    PURPOSE.
C     --------
C     READS RESTART WIND AND STRESS FIELDS FROM DISK OR USE GRIB WIND
C     AND DRAG COEFFICIENT TO RECONSTRUCT ALL WIND AND STRESS FIELDS.
C
C**   INTERFACE.
C     ----------
C     *CALL**GETSTRESS(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,NSTART,NEND,IREAD)
C     *U10OLD*    WIND SPEED.
C     *THWOLD*    WIND DIRECTION (RADIANS).
C     *USOLD*     FRICTION VELOCITY.
C     *TAUW*      WAVE STRESS.
C     *Z0OLD*     ROUGHNESS LENGTH IN.
C     *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
C     *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
C     *IREAD*     PROCESSOR WHICH WILL ACCESS THE FILE ON DISK
C
C     METHOD.
C     -------
C     USE READSTRESS TO READ IN THE FILE DEPENDING ON THE PBIO
C     SOFTWARE OR NOT. 
C     THE READING IS ONLY DONE ON PE 1, THEREFORE
C     THE RELEVANT INFORMATION IS SENT TO THE OTHER PE'S USING 
C      MPDISTRIBSCFLD
C
C     EXTERNALS.
C     ----------
C     BUILDSTRESS
C     GETENV
C     GRSTNAME
C     MPDISTRIBSCFLD
C     MPE_BARRIER
C     READSTRESS
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
#include "comgrid.h"
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
C
#include "comwind.h"
C
      CHARACTER*9 CMBXSIZE
      CHARACTER*10 ZERO
      CHARACTER*80 FILENAME
C
C ----------------------------------------------------------------------
      INTEGER,DIMENSION(NPROC) :: NSTART, NEND
      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10OLD,THWOLD,USOLD,Z0OLD,TAUW
C ----------------------------------------------------------------------
C     define local allocatable arrays
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IBUF
      REAL,ALLOCATABLE,DIMENSION(:,:) :: RU10OLD,RTHWOLD,RUSOLD,
     1                                   RZ0OLD,RTAUW

      ZERO = ' '

      IF(LMESSPASS) THEN

        ISEND=IREAD
        IBUFLENGTH=40

        ALLOCATE(RU10OLD(NIBLO,NBLO),RTHWOLD(NIBLO,NBLO),
     1  RUSOLD(NIBLO,NBLO),RZ0OLD(NIBLO,NBLO),RTAUW(NIBLO,NBLO))

        ALLOCATE(IBUF(IBUFLENGTH))

        IF (IRANK.EQ.ISEND) THEN
           IF (LPBIOIN) THEN
               CALL GRSTNAME(IU06,CDATEA,CDATEF,0,'LAW',IRANK,NPROC,
     &                       FILENAME)
           ELSE 
               CALL GETREST (IU06, IU12, IU14, IU15, IGL,
     1                       CDATEA, CDATEF)
           ENDIF
        ENDIF

C       READ RESTART FILE(s) FROM PE ISEND (IREAD) 
        IF (IRANK.EQ.ISEND) THEN
           IF(LGRIBIN) THEN
             CALL BUILDSTRESS(RU10OLD,RTHWOLD,RUSOLD,RTAUW,RZ0OLD,
     &                        NEND(NPROC))
           ELSE
             CALL READSTRESS(RU10OLD,RTHWOLD,RUSOLD,RTAUW,RZ0OLD,
     1                       FILENAME)
           ENDIF
        ENDIF

C       BROADCAST THE 4 DATES FROM RESTART FILE TO THE OTHER PE'S
        IF((ISEND.EQ.0).OR.(NPROC.EQ.1)) THEN
      
        ELSE IF(IRANK.EQ.ISEND) THEN
           KCOUNT=0
           DO IC=1,10
              READ(CDTPRO(IC:IC),'(I1)') IBUF(IC)
           ENDDO
           DO IC=1,10
              READ(CDATEWO(IC:IC),'(I1)') IBUF(10+IC)
           ENDDO
           DO IC=1,10
              READ(CDAWIFL(IC:IC),'(I1)') IBUF(20+IC)
           ENDDO
           DO IC=1,10
              READ(CDATEFL(IC:IC),'(I1)') IBUF(30+IC)
           ENDDO

           CALL MPE_BROADCAST
     1              (IBUF,IBUFLENGTH,1,IRANK,KTAG,0,0,0,IERR)
            IF(IERR.LT.0) CALL MPEI_ABORT
     1                    ('MPE_BROADCAST ERROR in GETSTRESS' )

            KTAG=KTAG+1
        ELSE

           CALL MPE_RECV(IBUF,IBUFLENGTH,1,ISEND,KTAG,0,0,0,
     $                   KRCOUNT,KRFROM,KRTAG,IERR)
           IF(IERR.LT.0) CALL MPEI_ABORT
     1                      ('MPE_RECV ERROR in GETSTRESS' )
           IF(KRTAG.NE.KTAG) CALL MPEI_ABORT
     1       ('MPE_RECV ERROR in GETSTRESS MISMATCHED TAGS' )

            KTAG=KTAG+1

            DO IC=1,10
               WRITE(CDTPRO(IC:IC),'(I1)') IBUF(IC)
            ENDDO
            DO IC=1,10
               WRITE(CDATEWO(IC:IC),'(I1)') IBUF(10+IC)
            ENDDO
            DO IC=1,10
               WRITE(CDAWIFL(IC:IC),'(I1)') IBUF(20+IC)
            ENDDO
            DO IC=1,10
               WRITE(CDATEFL(IC:IC),'(I1)') IBUF(30+IC)
            ENDDO
        ENDIF

        IF(CDTPRO.EQ.'0000000000') CDTPRO = ZERO
        IF(CDATEWO.EQ.'0000000000')CDATEWO = ZERO
        IF(CDAWIFL.EQ.'0000000000')CDAWIFL = ZERO
        IF(CDATEFL.EQ.'0000000000') CDATEFL= ZERO


C       DISTRIBUTE THE DIFFERENT CONTRIBUTIONS TO THE OTHER PE'S

        CALL MPDISTRIBSCFLD(ISEND,KTAG,NSTART,NEND,RU10OLD)
        KTAG=KTAG+1
        CALL MPDISTRIBSCFLD(ISEND,KTAG,NSTART,NEND,RTHWOLD)
        KTAG=KTAG+1
        CALL MPDISTRIBSCFLD(ISEND,KTAG,NSTART,NEND,RUSOLD)
        KTAG=KTAG+1
        CALL MPDISTRIBSCFLD(ISEND,KTAG,NSTART,NEND,RTAUW)
        KTAG=KTAG+1
        CALL MPDISTRIBSCFLD(ISEND,KTAG,NSTART,NEND,RZ0OLD)
        KTAG=KTAG+1

        IF (ITEST.GE.2)
     1   WRITE(IU06,*)
     2   'SUB. GETSTRESS: RESTART WIND AND STRESS COLLECTED'

C       KEEP CORRESPONDING CONTRIBUTION 
        DO IG=1,NBLO
           DO IJ=NSTART(IRANK),NEND(IRANK)
              U10OLD(IJ,IG)=RU10OLD(IJ,IG)
              THWOLD(IJ,IG)=RTHWOLD(IJ,IG)
              USOLD(IJ,IG)=RUSOLD(IJ,IG)
              TAUW(IJ,IG)=RTAUW(IJ,IG)
              Z0OLD(IJ,IG)=RZ0OLD(IJ,IG)
           END DO
        END DO

C
       DEALLOCATE(RU10OLD,RTHWOLD,RUSOLD,RZ0OLD,RTAUW,IBUF)
C
      ELSE

        IF (LPBIOIN) THEN
          CALL GRSTNAME(IU06,CDATEA,CDATEF,0,'LAW',IRANK,NPROC,FILENAME)
        ELSE 
            CALL GETREST (IU06, IU12, IU14, IU15, IGL,
     1                    CDATEA, CDATEF)
        ENDIF

        IF(LGRIBIN) THEN
          CALL BUILDSTRESS(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,NIBLO)
        ELSE
          CALL READSTRESS(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,FILENAME)
        ENDIF


      END IF
C
      WRITE(IU06,*) ''
      WRITE(IU06,*) ' WIND AND STRESS FILE READ IN........',
     1              ' CDTPRO  = ', CDTPRO
C
      RETURN
      END
