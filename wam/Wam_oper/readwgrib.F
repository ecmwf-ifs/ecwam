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
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,IRGG     ,AMOWEP   ,
     &            AMOSOP   ,AMOEAP   ,AMONOP   ,XDELLA   ,XDELLO   ,
     &            ZDELLO
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NBLO     ,NIBLO
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSTAT  , ONLY : NWAM_BLKS

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

!-----------------------------------------------------------------------

      INTEGER :: JKGLO,KIJS,KIJL,NPROMA
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

      IF (LHOOK) CALL DR_HOOK('READWGRIB',0,ZHOOK_HANDLE)

      CALL INWGRIB  (FILNM, IREAD, CCDDATE, KPARAM, KZLEV, WORK)
      KRGG=IRGG 

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

      IF(KRGG.EQ.0.AND.IRGG.EQ.1) THEN
!       INTERPOLATION IS REQUIRED : TAKE CLOSEST GRID POINT
!       (this part to be removed once we only use gribapi !)
        DELLO=XDELLO
        DELLA=XDELLA
        DO IG=1,IGL
          NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1

          IF(LLONLYPOS) THEN
            CALL GSTATS(1444,0)
!$OMP       PARALLEL DO SCHEDULE(STATIC)
!$OMP+      PRIVATE(JKGLO,KIJS,KIJL,IJ,IX,KX,XI,XJ,I,J)
            DO JKGLO=IJS(IG),IJL(IG),NPROMA
              KIJS=JKGLO
              KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
              DO IJ = KIJS, KIJL
                IX = IXLG(IJ,IG)
                KX = KXLT(IJ,IG)
                XI = AMOWEP + REAL(IX-1)*ZDELLO(KX)-CDMOWEP
                XI = MOD(XI+720.,360.)
                XI = XI/DELLO+1.00001
                XJ = CDMONOP - AMOSOP - REAL(KX-1)*XDELLA
                XJ = XJ/DELLA+1.00001
                I  = NINT(XI)
                J  = NINT(XJ)
                IF(WORK(I,J).NE.ZMISS .AND. WORK(I,J) .GT. 0.)
     &             FIELD(IJ,IG)=WORK(I,J)
              ENDDO
            ENDDO
!$OMP       END PARALLEL DO
            CALL GSTATS(1444,1)
          ELSE
            CALL GSTATS(1444,0)
!$OMP       PARALLEL DO SCHEDULE(STATIC)
!$OMP+      PRIVATE(JKGLO,KIJS,KIJL,IJ,IX,KX,XI,XJ,I,J)
            DO JKGLO=IJS(IG),IJL(IG),NPROMA
              KIJS=JKGLO
              KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
              DO IJ = KIJS, KIJL
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
!$OMP       END PARALLEL DO
            CALL GSTATS(1444,1)
          ENDIF
        ENDDO
      ELSE
!       NO INTERPOLATION IS REQUIRED
        DO IG=1,IGL
          NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1

          IF(LLONLYPOS) THEN
            CALL GSTATS(1444,0)
!$OMP       PARALLEL DO SCHEDULE(STATIC)
!$OMP+      PRIVATE(JKGLO,KIJS,KIJL,IJ,IX,IY)
            DO JKGLO=IJS(IG),IJL(IG),NPROMA
              KIJS=JKGLO
              KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
              DO IJ = KIJS, KIJL
                IX = IXLG(IJ,IG)
                IY = NGY- KXLT(IJ,IG) +1
                IF(WORK(IX,IY).NE.ZMISS .AND. WORK(IX,IY) .GT. 0.)
     &               FIELD(IJ,IG)=WORK(IX,IY)
              ENDDO
            ENDDO
!$OMP       END PARALLEL DO
            CALL GSTATS(1444,1)

          ELSE
            CALL GSTATS(1444,0)
!$OMP       PARALLEL DO SCHEDULE(STATIC) 
!$OMP+      PRIVATE(JKGLO,KIJS,KIJL,IJ,IX,IY)
            DO JKGLO=IJS(IG),IJL(IG),NPROMA
              KIJS=JKGLO
              KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
              DO IJ = KIJS, KIJL
                IX = IXLG(IJ,IG)
                IY = NGY- KXLT(IJ,IG) +1
                IF(WORK(IX,IY).NE.ZMISS) FIELD(IJ,IG)=WORK(IX,IY)
              ENDDO
            ENDDO
!$OMP       END PARALLEL DO
            CALL GSTATS(1444,1)
          ENDIF
        ENDDO
      ENDIF

      IF (LHOOK) CALL DR_HOOK('READWGRIB',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE READWGRIB
