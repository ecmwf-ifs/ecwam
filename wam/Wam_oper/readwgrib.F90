!-----------------------------------------------------------------------

      SUBROUTINE READWGRIB(IU06, FILNM, IPARAM, CDATE, MIJS, MIJL,      &
     &                     FIELD, KZLEV, LLONLYPOS, IREAD )

!-----------------------------------------------------------------------

!**** *READWGRIB*  READS FROM GRIB WAVE MODEL FIELD

!     J. BIDLOT   ECMWF   OCTOBER 1997 

!*    PURPOSE.
!     --------

!       SPECIAL INPUT FROM GRIB WAVE FIELD (SEE METHOD BELOW !!!) 

!**   INTERFACE.
!     ----------

!       *CALL* *READWGRIB*(IU06, FILNM, IPARAM, CDATE, MIJS, MIJL,
!    &                     FIELD, KZLEV, LLONLYPOS, IREAD )

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *IU06*      INTEGER   OUTPUT UNIT FOR STANDARD OUTPUT.
!      *FILNM*     DATA INPUT FILENAME.
!      *IPARAM*    INTEGER   PARAMETER IDENTIFIER OF FIELD
!      *CDATE*     CHARACTER DATE OF THE REQUESTED FIELD 
!      *MIJS*      INDEX OF FIRST GRIDPOINT
!      *MIJL*      INDEX OF LAST GRIDPOINT
!      *FIELD*     REAL      WAVE FIELD IN BLOCK FORMAT 
!      *KZLEV*     INTEGER   REFERENCE LEVEL IN full METER
!                           (SHOULD BE 0 EXCEPT FOR 233, 245 AND 249 WHERE IT
!                           MIGHT BE DIFFERENT THAN 0 PROVIDED IT WAS
!                           CODED, SEE *GRIBPAC*
!      *LLONLYPOS* LOGICAL   TAKE ONLY THE POSITIVE VALUES
!      *IREAD*     INTEGER  PROCESSOR WHICH WILL ACCESS THE FILE ON DISK


!     METHOD.
!     -------
!      READS GRIB DATA USING GRIB_API, CHECK WHETHER GRID DEFINITION ARE 
!      COMPATIBLE. IN CASE THE INPUT DATA IS DEFINED ON A REGULAR LAT-
!      LON GRID, BUT THE MODEL USES AN IRREGULAR GRID, THE DATA WILL BE 
!      INTERPOLATED TO THE IRREGULAR GRID BY TAKING THE NEAREST GRID
!      POINTS.
!      THE GRID DATA ARE PUT INTO BLOCK FORMAT PROVIDED THEY ARE
!      DIFFERENT THAN ZMISS TO AVOID THE USE OF MISSING DATA INDICATOR
!      FOR LAND POINT FOLLOWING A POSSIBLE CHANGE IN THE LAND-SEA
!      MASK . FOR THAT REASON, 
!      THE BLOCK VALUE OF FIELD SHOUD BE INITIALISED PRIOR TO THE CALL
!      TO THIS ROUTINE.

!     EXTERNALS.
!     ----------

!      *ABORT1*
!      *INWGRIB*

!     REFERENCE.
!     ----------

!       NONE.

!-------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWGRID  , ONLY : NLONRGG
      USE YOWMAP   , ONLY : IFROMIJ  ,JFROMIJ
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NBLO     ,NIBLO
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSTAT  , ONLY : NPROMA_WAM 

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!-----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "inwgrib.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IU06, IREAD, IPARAM
      INTEGER(KIND=JWIM), INTENT(IN) :: MIJS, MIJL
      INTEGER(KIND=JWIM), INTENT(INOUT) :: KZLEV

      REAL(KIND=JWRB),DIMENSION(MIJS:MIJL), INTENT(INOUT) :: FIELD 

      CHARACTER(LEN=14), INTENT(IN) :: CDATE
      CHARACTER(LEN=24), INTENT(IN) :: FILNM

      LOGICAL, INTENT(IN) :: LLONLYPOS


      INTEGER(KIND=JWIM) :: KPARAM
      INTEGER(KIND=JWIM) :: IJ, IX, JY
      INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA
      INTEGER(KIND=JWIM) :: KLONRGG(NGY)

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: WORK(NGX,NGY)

      CHARACTER(LEN=14) :: CCDDATE
      CHARACTER(LEN=40) MSG


!-----------------------------------------------------------------------

!*    1. INPUT OF GRIB DATA.
!     -----------------------

      IF (LHOOK) CALL DR_HOOK('READWGRIB',0,ZHOOK_HANDLE)

      CALL INWGRIB  (FILNM, IREAD, CCDDATE, KPARAM, KZLEV, WORK)

!*    SIMPLE CHECKS ON THE RETRIEVED DATA 
!     -----------------------------------

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


! TRANSFORM GRID DATA TO BLOCK DATA

          NPROMA=NPROMA_WAM

          IF(LLONLYPOS) THEN
            CALL GSTATS(1444,0)
!$OMP       PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ,IX,JY)
            DO JKGLO=MIJS,MIJL,NPROMA
              KIJS=JKGLO
              KIJL=MIN(KIJS+NPROMA-1,MIJL)
              DO IJ = KIJS, KIJL
                IX = IFROMIJ(IJ,1)
                JY = JFROMIJ(IJ,1)
                IF(WORK(IX,JY).NE.ZMISS .AND. WORK(IX,JY).GT.0.0_JWRB) FIELD(IJ)=WORK(IX,JY)
              ENDDO
            ENDDO
!$OMP       END PARALLEL DO
            CALL GSTATS(1444,1)

          ELSE
            CALL GSTATS(1444,0)
!$OMP       PARALLEL DO SCHEDULE(STATIC)  PRIVATE(JKGLO,KIJS,KIJL,IJ,IX,JY)
            DO JKGLO=MIJS,MIJL,NPROMA
              KIJS=JKGLO
              KIJL=MIN(KIJS+NPROMA-1,MIJL)
              DO IJ = KIJS, KIJL
                IX = IFROMIJ(IJ,1)
                JY = JFROMIJ(IJ,1)
                IF(WORK(IX,JY).NE.ZMISS) FIELD(IJ)=WORK(IX,JY)
              ENDDO
            ENDDO
!$OMP       END PARALLEL DO
            CALL GSTATS(1444,1)
          ENDIF

      IF (LHOOK) CALL DR_HOOK('READWGRIB',1,ZHOOK_HANDLE)

      END SUBROUTINE READWGRIB
