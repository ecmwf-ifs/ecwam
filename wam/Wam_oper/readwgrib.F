!-----------------------------------------------------------------------

      SUBROUTINE READWGRIB(IU06, FILNM, IPARAM, CDATE, FIELD,
     &                     KZLEV, LLONLYPOS, IREAD )

!-----------------------------------------------------------------------

!**** *READWGRIB*  READS FROM GRIB WAVE MODEL FIELD

!     J. BIDLOT   ECMWF   OCTOBER 1997 

!*    PURPOSE.
!     --------

!       INPUT FROM GRIB WAVE FIELD 

!**   INTERFACE.
!     ----------

!       *CALL* *READWGRIB*(IU06, FILNM, IPARAM, CDATE, FIELD,
!    &                     KZLEV, LLONLYPOS, IREAD )

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *IU06*      INTEGER   OUTPUT UNIT FOR STANDARD OUTPUT.
!      *FILNM*     DATA INPUT FILENAME.
!      *IPARAM*    INTEGER   PARAMETER IDENTIFIER OF FIELD
!      *CDATE*     CHARACTER DATE OF THE REQUESTED FIELD 
!      *FIELD*     REAL      WAVE FIELD IN BLOCK FORMAT 
!      *KZLEV*     INTEGER   REFERENCE LEVEL IN full METER
!                           (SHOULD BE 0 EXCEPT FOR 233, 245 AND 249 WHERE IT
!                           MIGHT BE DIFFERENT THAN 0 PROVIDED IT WAS
!                           CODED, SEE *GRIBPAC*
!      *LLONLYPOS* LOGICAL   TAKE ONLY THE POSITIVE VALUES
!      *IREAD*     INTEGER  PROCESSOR WHICH WILL ACCESS THE FILE ON DISK


!     METHOD.
!     -------
!      READS GRIB DATA USING INMARS, CHECK WHETHER GRID DEFINITION ARE 
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
!      *INMARS*

!     REFERENCE.
!     ----------

!       NONE.

!-------------------------------------------------------------------

      USE YOWGRID  , ONLY : NLONRGG  ,IGL      ,IJS      ,IJL
      USE YOWMAP   , ONLY : IFROMIJ  ,JFROMIJ
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NBLO     ,NIBLO
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSTAT  , ONLY : NWAM_BLKS

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

!-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IU06, IPARAM, KZLEV, IREAD
      INTEGER :: KPARAM 
      INTEGER :: IG 
      INTEGER :: IJ, IX, JY
      INTEGER :: JKGLO, KIJS, KIJL, NPROMA
      INTEGER :: KLONRGG(NGY)

      REAL :: ZHOOK_HANDLE
      REAL :: WORK(NGX,NGY)
      REAL,DIMENSION(NINF:NSUP,NBLO) :: FIELD 

      CHARACTER(LEN=14) :: CCDDATE,CDATE 
      CHARACTER(LEN=24) :: FILNM
      CHARACTER(LEN=40) MSG

      LOGICAL :: LLONLYPOS

!-----------------------------------------------------------------------

!*    1. INPUT OF GRIB DATA.
!     -----------------------

#ifdef ECMWF 
      IF (LHOOK) CALL DR_HOOK('READWGRIB',0,ZHOOK_HANDLE)
#endif

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

        DO IG=1,IGL
          NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1

          IF(LLONLYPOS) THEN
            CALL GSTATS(1444,0)
!$OMP       PARALLEL DO SCHEDULE(STATIC)
!$OMP+      PRIVATE(JKGLO,KIJS,KIJL,IJ,IX,JY)
            DO JKGLO=IJS(IG),IJL(IG),NPROMA
              KIJS=JKGLO
              KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
              DO IJ = KIJS, KIJL
                IX = IFROMIJ(IJ,IG)
                JY = JFROMIJ(IJ,IG)
                IF(WORK(IX,JY).NE.ZMISS .AND. WORK(IX,JY) .GT. 0.)
     &               FIELD(IJ,IG)=WORK(IX,JY)
              ENDDO
            ENDDO
!$OMP       END PARALLEL DO
            CALL GSTATS(1444,1)

          ELSE
            CALL GSTATS(1444,0)
!$OMP       PARALLEL DO SCHEDULE(STATIC) 
!$OMP+      PRIVATE(JKGLO,KIJS,KIJL,IJ,IX,JY)
            DO JKGLO=IJS(IG),IJL(IG),NPROMA
              KIJS=JKGLO
              KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
              DO IJ = KIJS, KIJL
                IX = IFROMIJ(IJ,IG)
                JY = JFROMIJ(IJ,IG)
                IF(WORK(IX,JY).NE.ZMISS) FIELD(IJ,IG)=WORK(IX,JY)
              ENDDO
            ENDDO
!$OMP       END PARALLEL DO
            CALL GSTATS(1444,1)
          ENDIF
        ENDDO

#ifdef ECMWF 
      IF (LHOOK) CALL DR_HOOK('READWGRIB',1,ZHOOK_HANDLE)
#endif

      RETURN
      END SUBROUTINE READWGRIB
