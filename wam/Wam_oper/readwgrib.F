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
!                           CODED IN KSEC1(7) AND KSEC1(8) SEE *GRIBPAC*
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

      USE YOWCOUT  , ONLY : LGRIBAPI_I
      USE YOWGRIBHD, ONLY : NKSEK0   ,NKSEK1   ,NKSEK2   ,NKSEK3,
     &            NKSEK4   ,NPSEK2   ,NPSEK3   ,KSEK3    ,PSEK3
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

      IF(LGRIBAPI_I) THEN
        CALL INWGRIB  (FILNM, IREAD,
     &                 CCDDATE, KPARAM, KZLEV, WORK)
        KRGG=IRGG 
      ELSE
        CALL INMARS  (FILNM, NGX, NGY, ZMISS, NKSEK0, NKSEK1, NKSEK2,
     &                NKSEK3, NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
     &                CCDDATE, IFORP, KPARAM, NGXCD, NGYCD,
     &                CDMOWEP, CDMOSOP, CDMOEAP, CDMONOP,KRGG,KLONRGG,
     &                KZLEV,WORK, IERR, MSG, IU06,
     &                IREAD, IRANK, NPROC)
        IF ( IERR .EQ. 0 ) THEN
          WRITE(IU06,*)'            '
          WRITE(IU06,*)' SUB. READWGRIB RETURNING FROM INMARS: DATE=',
     &                CCDDATE,' IPARAM= ',IPARAM
        ELSE
          WRITE(IU06,*)' SUB. READWGRIB RETURNING FROM INAMRS: IERR=',
     &                IERR,' IPARAM= ',IPARAM
          WRITE(IU06,*) MSG
          CALL ABORT1
        ENDIF
        CALL ADJUST (CDMOWEP, CDMOEAP)
        KAMOWEP=NINT(AMOWEP*100.)
        KAMOEAP=NINT(AMOEAP*100.)
        KAMONOP=NINT(AMONOP*100.)
        KAMOSOP=NINT(AMOSOP*100.)
        KCDMOWEP=NINT(CDMOWEP*100.)
        KCDMOEAP=NINT(CDMOEAP*100.)
        KCDMONOP=NINT(CDMONOP*100.)
        KCDMOSOP=NINT(CDMOSOP*100.)
        IF (KAMOWEP.NE.KCDMOWEP .OR. KAMOEAP.NE.KCDMOEAP .OR.
     &      KAMONOP.NE.KCDMONOP .OR. KAMOSOP.NE.KCDMOSOP) THEN
          WRITE(IU06,*) '***********************************'
          WRITE(IU06,*) '*                                 *'
          WRITE(IU06,*) '*  FATAL ERROR IN SUB. READWGRIB  *'
          WRITE(IU06,*) '*  ============================   *'
          WRITE(IU06,*) '*                                 *'
          WRITE(IU06,*) '*  THIS IS NOT THE SAME GRID AREA *'
          WRITE(IU06,*) '*  AS YOU HAD DELARED IN THE      *'
          WRITE(IU06,*) '*  COARSE MODEL GRID              *'
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

      ENDIF


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
