      SUBROUTINE BUILDSTRESS(RU10OLD,RTHWOLD,RUSOLD,RTAUW,RZ0OLD, KEND)

C ----------------------------------------------------------------------
C     J. BIDLOT    ECMWF     APRIL 1998 
C
C*    PURPOSE.
C     --------
C     CREATES WIND AND STRESS FIELDS FROM GRIB WINDS AND CD.
C
C**   INTERFACE.
C     ----------
C     CALL *BUILDSTRESS*(RU10OLD,RTHWOLD,RUSOLD,RTAUW,RZ0OLD)
C     *RU10OLD*   WIND SPEED.
C     *RTHWOLD*   WIND DIRECTION (RADIANS).
C     *RUSOLD*    FRICTION VELOCITY.
C     *RTAUW*     WAVE STRESS.
C     *RZ0OLD*    ROUGHNESS LENGTH IN M.
C     *KEND*      USED UPPER LIMIT ON FIRST DIMENSION OF INPUT ARRAYS
C
C     METHOD.
C     -------
C
C     EXTERNALS.
C     ----------
C     *ABORT1*
C     *AIRSEA*
C     *GETWND*
C     *PBOPEN*
C     *PBREAD*
C     *PBCLOSE*
C     *READWGRIB*
C
C     REFERENCE.
C     ----------
C     NONE
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "parcons.h"
C
#include "parameter_wind.h"
C
#include "comcoup.h"
C
#include "comgrid.h"
C
#include "comstat.h"
C
#include "comtest.h"
C
#include "comunit.h"
C
#include "comwind.h"
C
#include "commesspass.h"
C
#include "commpp.h"
C
      LOGICAL LWNDFILE, LWAVEWIND, LCLOSEWND
      INTEGER IUGB(1),IVGB(1)
      REAL, ALLOCATABLE :: UWND(:,:),VWND(:,:)
      REAL, ALLOCATABLE :: CD(:,:)
      CHARACTER*24 FILNM
      REAL,DIMENSION(NIBLO,NBLO) :: RU10OLD,RTHWOLD,RUSOLD,RZ0OLD,RTAUW
C ----------------------------------------------------------------------
C
      CDATEWO = ' '
      CDAWIFL = ' '
      CDATEFL = ' '
      CDTPRO = CDATEA

C
C     GET ATMOSPHERIC MODEL 10m WINDS

      ILEN=1
      LWNDFILE=.TRUE.
      LCLOSEWND=.TRUE.
      RU10OLD = 0.
      RTHWOLD = 0.
      RUSOLD = 0.
      RTAUW = 0.
      RZ0OLD = 0.
      IF(.NOT.ALLOCATED(UWND)) ALLOCATE(UWND(NC,NR))
      IF(.NOT.ALLOCATED(VWND)) ALLOCATE(VWND(NC,NR))
      CALL GETWND (RU10OLD, RTHWOLD,1,NIBLO, KEND, CDTPRO, UWND, VWND,
     1              NC, NR, IUGB, IVGB, ILEN, LWNDFILE, LCLOSEWND)
      IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. GETWND DONE'
      IF(ALLOCATED(UWND)) DEALLOCATE(UWND)
      IF(ALLOCATED(VWND)) DEALLOCATE(VWND)
C
C     USE DATA FROM A FILE CONTAINING WIND SPEED MODIFIED BY
C     A PREVIOUS WAVE MODEL RUN ON WHICH THIS RESTART IS BASED
C     IF IT EXISTS, OTHERWISE PRESET WILL USE THE INFORMATION FROM
C     THE ATMOSPHERIC MODEL WIND SPEED ARCHIVE.
C
      LWAVEWIND=.FALSE.
C
      FILNM='uwavein'
      LIU = IECF_LEN(FILNM)
      FILNM=FILNM(1:LIU)
      INQUIRE(FILE=FILNM,EXIST=LWAVEWIND)
C
      IF(LWAVEWIND) THEN
C
        ILEV=2
        WRITE(IU06,*) ' '
        WRITE(IU06,*) 'A DATA FILE CONTAINING WIND SPEED INFORMATION'
        WRITE(IU06,*) 'AS PROVIDED BY A PREVIOUS WAVE MODEL RUN'
        WRITE(IU06,*) 'WAS USED TO UPDATE THE INPUT ATMOSPHERIC WINDS'
        WRITE(IU06,*) ' '
        WRITE(IU06,*) 'THE WINDS AND DRAG COEFFICIENT ARE ASSUMED TO'
        WRITE(IU06,*) 'HAVE BEEN DETERMINED FOR THE LOWEST LEVEL OF '
        WRITE(IU06,*) 'THE ATMOSPHERIC MODEL AT HEIGHT OF '
        WRITE(IU06,*) XNLEV(ILEV),' m'
        IF (ITEST.GT.0) CALL FLUSH(IU06) 
        IPARAM=245
        CALL READWGRIB(IU06,FILNM,IPARAM,CDTPRO,RU10OLD,1,NIBLO,KEND,
     &                 KZLEVUWAVE)
        IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. READWGRIB DONE FOR ',FILNM
C
      ELSE
C
        ILEV=1
        KZLEVUWAVE=10
        WRITE(IU06,*) ' '
        WRITE(IU06,*) '          !!!! NOTE !!!!'
        WRITE(IU06,*) ' '
        WRITE(IU06,*) 'NO INFORMATION ON WIND SPEEDS FROM THE'
        WRITE(IU06,*) 'WAVE MODEL WAS PROVIDED'
        WRITE(IU06,*) 'NO UPDATE OF THE INPUT ATMOSPHERIC WINDS'
        WRITE(IU06,*) 'WAS POSSIBLE'
        WRITE(IU06,*) ' '
        WRITE(IU06,*) 'THE WINDS AND DRAG COEFFICIENT ARE ASSUMED TO'
        WRITE(IU06,*) 'HAVE BEEN DETERMINED FOR HEIGHT AT '
        WRITE(IU06,*)  XNLEV(ILEV),' m'
        IF (ITEST.GT.0) CALL FLUSH(IU06) 
      ENDIF
C
      IF(.NOT.ALLOCATED(CD)) ALLOCATE(CD(NIBLO,NBLO))
C
C     INITIALISE CD USING THE FRICTION VELOCITY FOR TAUW=0.
C
      DO IG=1,IGL
         CALL AIRSEA (RU10OLD(1,IG), RTAUW(1,IG),RUSOLD(1,IG),
     1                RZ0OLD(1,IG), 1, KEND, ILEV)
         IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. AIRSEA DONE AT 1'
C
         DO IJ=1,KEND
            CDINV = (MIN(RU10OLD(IJ,IG)/RUSOLD(IJ,IG),100.0))**2
            CD(IJ,IG) = 1./CDINV
         ENDDO
      ENDDO
C
      IF(.NOT.LNOCDIN) THEN
        IPARAM=233
        FILNM='cdwavein'
        CALL READWGRIB(IU06,FILNM,IPARAM,CDTPRO,CD,1,NIBLO,KEND,KZLEVCD)
        IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. READWGRIB DONE FOR ',FILNM
C
C       TEST REFERENCE LEVEL FOR UWAVE AND CD
C
        IF(KZLEVUWAVE.NE.0.AND.KZLEVCD.NE.0) THEN
          IF(KZLEVUWAVE.NE.KZLEVCD) THEN
            WRITE(IU06,*)'************************************'
            WRITE(IU06,*)'*                                  *'
            WRITE(IU06,*)'* FATAL ERROR IN SUB BUILDSTRESS   *'
            WRITE(IU06,*)'* ==============================   *'
            WRITE(IU06,*)'* REFERENCE LEVEL FOR CD AND UWAVE *'
            WRITE(IU06,*)'* DIFFER. PROGRAM WILL ABORT       *'
            WRITE(IU06,*)'* REF. LEVEL CD = ',KZLEVCD
            WRITE(IU06,*)'* REF. LEVEL UWAVE = ',KZLEVUWAVE
            WRITE(IU06,*)'*                                  *'
            WRITE(IU06,*)'************************************'
            CALL ABORT1
          ELSE
            IF (KZLEVCD.EQ.INT(XNLEV(1))) THEN
              ILEV=1
            ELSE IF (KZLEVCD.EQ.INT(XNLEV(2))) THEN
              ILEV=2
            ELSE
              WRITE(IU06,*)'**************************************'
              WRITE(IU06,*)'*                                    *'
              WRITE(IU06,*)'* FATAL ERROR IN SUB BUILDSTRESS     *'
              WRITE(IU06,*)'* ==============================     *'
              WRITE(IU06,*)'* REFERENCE LEVEL FOR CD AND UWAVE   *'
              WRITE(IU06,*)'* ARE NOT RECOGNISED AS THOSE IN USE *'
              WRITE(IU06,*)'* PROGRAM WILL ABORT                 *'
              WRITE(IU06,*)'* REF. LEVEL CD = ',KZLEVCD
              WRITE(IU06,*)'* REF. LEVEL UWAVE = ',KZLEVUWAVE
              WRITE(IU06,*)'*                                    *'
              WRITE(IU06,*)'**************************************'
              CALL ABORT1
            ENDIF
          ENDIF
        ENDIF
C
C
C       COMPUTE TAUW,USOLD AND Z0OLD
C
        DO IG=1,IGL
           DO IJ=1,KEND
              CDSQRTINV = MIN(1./SQRT(CD(IJ,IG)),100.0)
              RZ0OLD(IJ,IG) = XNLEV(ILEV)*EXP(-XKAPPA*CDSQRTINV)
cc RUSOLD will first contain its square
CC the numerical relation between usold and u10old should always be as
CC in outgrid
              RUSOLD(IJ,IG) = CD(IJ,IG)*(RU10OLD(IJ,IG)**2+0.01)-0.0001
              CHARNOCK = G*RZ0OLD(IJ,IG)/RUSOLD(IJ,IG)
              CHARNOCK = MAX(CHARNOCK,ALPHA)
              RTAUW(IJ,IG) = RUSOLD(IJ,IG)*(1.-(ALPHA/CHARNOCK)**2)
              RUSOLD(IJ,IG) = SQRT(RUSOLD(IJ,IG))
           ENDDO
        ENDDO
      ENDIF
C
      IF(ALLOCATED(CD)) DEALLOCATE(CD)
C
      IF (ITEST.GE.2)
     1    WRITE(IU06,*) ' SUB. BUILDSTRESS: INPUT OF RESTART FILES DONE'
C
      RETURN
      END
