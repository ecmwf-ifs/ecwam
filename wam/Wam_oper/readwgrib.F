C------------------------------------------------------------------
C
      SUBROUTINE READWGRIB(IU06,FILNM,IPARAM,CDATE,FIELD,KINF,KSUP,
     &                     KEND, KZLEV)
C
C----------------------------------------------------------------------
C
C**** *READWGRIB*  READS FROM GRIB WAVE MODEL FIELD
C
C     J. BIDLOT   ECMWF   OCTOBER 1997 
C
C*    PURPOSE.
C     --------
C
C       INPUT FROM GRIB WAVE FIELD 
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *READWGRIB*(IU06,FILNM,IPARAM,CDATE,FIELD)
C
C*     VARIABLE.   TYPE.     PURPOSE.
C      ---------   -------   --------
C      *IU06*      INTEGER   OUTPUT UNIT FOR STANDARD OUTPUT.
C      *FILNM*     DATA INPUT FILENAME.
C      *IPARAM*    INTEGER   PARAMETER IDENTIFIER OF FIELD
C      *CDATE*     CHARACTER DATE OF THE REQUESTED FIELD 
C      *FIELD*     REAL      WAVE FIELD IN BLOCK FORMAT 
C      *KINF*      INTEGER   LOWER BOUND ON FIRST DIEMENSION OF FIELD
C      *KSUP*      INTEGER   UPPER BOUND ON FIRST DIEMENSION OF FIELD
C      *KEND*      INTEGER   USED UPPER LIMIT ON 1st DIMENSION OF FIELD
C      *KZLEV*     INTEGER   REFERENCE LEVEL IN full METER
C                           (SHOULD BE 0 EXCEPT FOR 233 AND 245 WHERE IT
C                           MIGHT BE DIFFERENT THAN 0 PROVIDED IT WAS
C                           CODED IN KSEC1(7) AND KSEC1(8) SEE *GRIBPAC*
C
C     METHOD.
C     -------
C      READS GRIB DATA USING INMARS, CHECK WHETHER GRID DEFINITION ARE 
C      COMPATIBLE. IN CASE THE INPUT DATA IS DEFINED ON A REGULAR LAT-
C      LON GRID, BUT THE MODEL USES AN IRREGULAR GRID, THE DATA WILL BE 
C      INTERPOLATED TO THE IRREGULAR GRID BY TAKING THE NEAREST GRID
C      POINTS.
C      THE GRID DATA ARE PUT INTO BLOCK FORMAT PROVIDED THEY ARE
C      DIFFERENT THAN ZMISS TO AVOID THE USE OF MISSING DATA INDICATOR
C      FOR LAND POINT FOLLOWING A POSSIBLE CHANGE IN THE LAND-SEA
C      MASK . FOR THAT REASON, 
C      THE BLOCK VALUE OF FIELD SHOUD BE INITIALISED PRIOR TO THE CALL
C      TO THIS ROUTINE.
C      
C     EXTERNALS.
C     ----------
C
C      *ABORT1*
C      *INMARS*
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C-------------------------------------------------------------------
C
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.
C
#include "param.h"
C
#include "parcons.h"
C
#include "comgrid.h"
C
#include "commap.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
C
C     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS
      REAL,DIMENSION(KINF:KSUP,NBLO) :: FIELD 
C
C-----------------------------------------------------------------------
      REAL,ALLOCATABLE ::  WORK(:,:)
      INTEGER,ALLOCATABLE :: KLONRGG(:)
C
      CHARACTER*12 CCDDATE,CDATE 
      CHARACTER*24 FILNM
      CHARACTER*40 MSG, TITL*100
C-----------------------------------------------------------------------
C
C
C*    1. INPUT OF GRIB DATA.
C     -----------------------
C
 1000 CONTINUE
C
      ALLOCATE(WORK(NGX,NGY),KLONRGG(NGY))

      CALL INMARS  ( FILNM, NGX, NGY, ZMISS, CCDDATE, IFORP, KPARAM,
     .               NGXCD, NGYCD, CDMOWEP, CDMOSOP, CDMOEAP, CDMONOP,
     .               KRGG, KLONRGG, KZLEV,WORK, IERR, MSG)
      IF ( IERR .EQ. 0 ) THEN
        WRITE(IU06,*)'            '
        WRITE(IU06,*)'SUB. READWGRIB: DATE=',CCDDATE,' IPARAM= ',IPARAM
      ELSE
        WRITE(IU06,*)'SUB. READWGRIB: IERR=',IERR,' IPARAM= ',IPARAM
        WRITE(IU06,*) MSG
        CALL ABORT1
      ENDIF
C
C
C*    1.1 ERROR MESSAGES
C     ------------------
C
 1100 CONTINUE
C
      IF (KPARAM.NE.IPARAM) THEN
         WRITE(IU06,*)'********************************'
         WRITE(IU06,*)'*                              *'
         WRITE(IU06,*)'* FATAL ERROR IN SUB READWGRIB *'
         WRITE(IU06,*)'* ===========================  *'
         WRITE(IU06,*)'*                              *'
         WRITE(IU06,*)'* GRIB PARAMETER  ',KPARAM
         WRITE(IU06,*)'* WAS READ INSTEAD OF ',IPARAM
         WRITE(IU06,*)'* IN FILE: ',FILNM 
         WRITE(IU06,*)'*                              *'
         WRITE(IU06,*)'********************************'
         CALL ABORT1
      ENDIF
      IF (CCDDATE.NE.CDATE) THEN
         WRITE(IU06,*)'**********************************'
         WRITE(IU06,*)'*                                *'
         WRITE(IU06,*)'* FATAL ERROR IN SUB READWGRIB   *'
         WRITE(IU06,*)'* ===========================    *'
         WRITE(IU06,*)'*                                *'
         WRITE(IU06,*)'* REQUESTED DATE IS NOT EQUAL TO *' 
         WRITE(IU06,*)'* RETRIEVED DATE.                *' 
         WRITE(IU06,*)'* IN FILE: ',FILNM 
         WRITE(IU06,*)'* CCDDATE = ',CCDDATE
         WRITE(IU06,*)'* CDATE = ',CDATE
         WRITE(IU06,*)'*                                *'
         WRITE(IU06,*)'**********************************'
         CALL ABORT1
      ENDIF
      CALL ADJUST (CDMOWEP, CDMOEAP)
      IF (AMOWEP.NE.CDMOWEP .OR. AMOEAP.NE.CDMOEAP .OR.
     +    AMONOP.NE.CDMONOP .OR. AMOSOP.NE.CDMOSOP) THEN
          WRITE(IU06,*) '***********************************'
          WRITE(IU06,*) '*                                 *'
          WRITE(IU06,*) '*  FATAL ERROR IN SUB. READWGRIB  *'
          WRITE(IU06,*) '*  ============================   *'
          WRITE(IU06,*) '*                                 *'
          WRITE(IU06,*) '*  THIS IS NOT THE SAME GRID AREA *'
          WRITE(IU06,*) '*  AS YOU HAD DELARED IN THE      *'
          WRITE(IU06,*) '*  COURSE MODEL GRID              *'
          WRITE(IU06,*) '* IN FILE: ',FILNM 
          WRITE(IU06,*) '*                                 *'
          WRITE(IU06,*) '*  AMOSOP: ', AMOSOP, 'CDMOSOP : ',CDMOSOP
          WRITE(IU06,*) '*  AMONOP: ', AMONOP, 'CDMONOP : ',CDMONOP
          WRITE(IU06,*) '*  AMOWEP: ', AMOWEP, 'CDMOWEP : ',CDMOWEP
          WRITE(IU06,*) '*  AMOEAP: ', AMOEAP, 'CDMOEAP : ',CDMOEAP
          WRITE(IU06,*) '*                                 *'
          WRITE(IU06,*) '*     THE PROGRAM ABORTS          *'
          WRITE(IU06,*) '***********************************'
          CALL ABORT1
      ENDIF
      IF (KRGG.EQ.0) THEN
        IF (NGXCD .NE. NGX) THEN
           WRITE(IU06,*)'********************************'
           WRITE(IU06,*)'*                              *'
           WRITE(IU06,*)'* FATAL ERROR IN SUB READWGRIB *'
           WRITE(IU06,*)'* ===========================  *'
           WRITE(IU06,*)'*                              *'
           WRITE(IU06,*)'* NGXCD SHOULD BE EQUAL TO NGX *'
           WRITE(IU06,*)'* IN FILE: ',FILNM
           WRITE(IU06,*)'* NGX  = ',NGX
           WRITE(IU06,*)'* NGXCD  = ',NGXCD
           WRITE(IU06,*)'*                              *'
           WRITE(IU06,*)'********************************'
           CALL ABORT1
        ENDIF
      ENDIF
      IF (NGYCD .NE. NGY) THEN
         WRITE(IU06,*)'********************************'
         WRITE(IU06,*)'*                              *'
         WRITE(IU06,*)'* FATAL ERROR IN SUB READWGRIB *'
         WRITE(IU06,*)'* ===========================  *'
         WRITE(IU06,*)'*                              *'
         WRITE(IU06,*)'* NGYCD SHOULD BE EQUAL TO NGY *'
         WRITE(IU06,*)'* IN FILE: ',FILNM
         WRITE(IU06,*)'* NGY  = ',NGY
         WRITE(IU06,*)'* NGYCD  = ',NGYCD
         WRITE(IU06,*)'*                              *'
         WRITE(IU06,*)'********************************'
         CALL ABORT1
      ENDIF
      IF (KRGG.NE.IRGG) THEN
         WRITE(IU06,*)'********************************'
         WRITE(IU06,*)'*                              *'
         WRITE(IU06,*)'* ERROR IN SUB READWGRIB       *'
         WRITE(IU06,*)'* =====================        *'
         WRITE(IU06,*)'*                              *'
         WRITE(IU06,*)'* KRGG IS NOT EQUAL TO IRGG    *'
         WRITE(IU06,*)'* KRGG  = ',KRGG
         WRITE(IU06,*)'* IRGG  = ',IRGG
         WRITE(IU06,*)'*                              *'
         IF(KRGG.EQ.0) THEN
         WRITE(IU06,*)'* WARNING INTERPOLATION WILL   *'
         WRITE(IU06,*)'* TAKE PLACE FROM REGULAR GRID *'
         WRITE(IU06,*)'* DATA TO IRREGULAR GRID       *'
         WRITE(IU06,*)'*                              *'
         WRITE(IU06,*)'********************************'
         ELSE
         WRITE(IU06,*)'* PROGRAM CANNOT INTERPOLATE   *'
         WRITE(IU06,*)'* BETWEEN THOSE 2 GRID         *'
         WRITE(IU06,*)'* IT WILL ABORT                *'
         WRITE(IU06,*)'*                              *'
         WRITE(IU06,*)'********************************'
         CALL ABORT1
         ENDIF
      ENDIF
      IF (KRGG.EQ.1) THEN
        DO J=1,NGY
          IF(KLONRGG(J).NE.NLONRGG(J)) THEN
             WRITE(IU06,*)'***********************************'
             WRITE(IU06,*)'*                                 *'
             WRITE(IU06,*)'* ERROR IN SUB READWGRIB          *'
             WRITE(IU06,*)'* =====================           *'
             WRITE(IU06,*)'*                                 *'
             WRITE(IU06,*)'* THE IRREGULAR GRID DEFINITION   *'
             WRITE(IU06,*)'* OF THE INPUT FILE DO NOT MATCH  *'
             WRITE(IU06,*)'* WITH THE MODEL ONE.             *' 
             WRITE(IU06,*)'* KLONRGG IS NOT EQUAL TO NLONRGG *'
             DO JJ=1,NGY
               WRITE(IU06,*)'* JJ = ',JJ
               WRITE(IU06,*)'* KLONRGG  = ',KLONRGG(JJ)
               WRITE(IU06,*)'* NLONRGG  = ',NLONRGG(JJ)
             ENDDO
             WRITE(IU06,*)'*                                 *'
             WRITE(IU06,*)'***********************************'
             CALL ABORT1
          ENDIF
        ENDDO
      ENDIF
C
C TRANSFORM GRID DATA TO BLOCK DATA
C
      IF(KRGG.EQ.0.AND.IRGG.EQ.1) THEN
C   INTERPOLATION IS REQUIRED : TAKE CLOSEST GRID POINT
        DELLO=XDELLO
        DELLA=XDELLA
        IF(KINF.EQ.1.AND.KSUP.EQ.NIBLO) THEN
          IG=1
             DO IJ = 1,KEND 
                IX = IXLG(IJ,IG)
                KX = KXLT(IJ,IG)
                XI = AMOWEP + REAL(IX-1)*ZDELLO(KX)-CDMOWEP
                XI = MOD(XI+720.,360.)
                XI = XI/DELLO+1.00001
                XJ = CDMONOP - AMOSOP - REAL(KX-1)*XDELLA
                XJ = XJ/DELLA+1.00001
                I  = NINT(XI)
                J  = NINT(XJ)
                IF(WORK(I,J).NE.ZMISS) FIELD(IJ,IG)=WORK(I,J)
             ENDDO
        ELSE
          DO IG=1,IGL
             DO IJ = IJS(IG),IJL(IG)
                IX = IXLG(IJ,IG)
                KX = KXLT(IJ,IG)
                XI = AMOWEP + REAL(IX-1)*ZDELLO(KX)-CDMOWEP
                XI = MOD(XI+720.,360.)
                XI = XI/DELLO+1.00001
                XJ = CDMONOP - AMOSOP - REAL(KX-1)*XDELLA
                XJ = XJ/DELLA+1.00001
                I  = NINT(XI)
                J  = NINT(XJ)
                IF(WORK(I,J).NE.ZMISS) FIELD(IJ,IG)=WORK(I,J)
             ENDDO
          ENDDO
        ENDIF
      ELSE
C   NO INTERPOLATION IS REQUIRED
        IF(KINF.EQ.1.AND.KSUP.EQ.NIBLO) THEN
          IG=1
             DO IJ = 1,KEND 
                IX = IXLG(IJ,IG)
                IY = NGY- KXLT(IJ,IG) +1
                IF(WORK(IX,IY).NE.ZMISS) FIELD(IJ,IG)=WORK(IX,IY)
             ENDDO
        ELSE
          DO IG=1,IGL
             DO IJ = IJS(IG),IJL(IG)
                IX = IXLG(IJ,IG)
                IY = NGY- KXLT(IJ,IG) +1
                IF(WORK(IX,IY).NE.ZMISS) FIELD(IJ,IG)=WORK(IX,IY)
             ENDDO
          ENDDO
        ENDIF
      ENDIF
C
      IF(ALLOCATED(WORK))DEALLOCATE(WORK)
      IF(ALLOCATED(KLONRGG))DEALLOCATE(KLONRGG)
C
      RETURN
      END
