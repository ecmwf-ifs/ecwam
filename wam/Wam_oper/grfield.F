      SUBROUTINE GRFIELD (WHME, USME, NSTART, NEND)

!--------------------------------------------------------------------

!**** *GRFIELD* - READS GRIDDED MEASUREMENTS AND DISTRIBUTES THEM.

!     B.HANSEN       ECMWF       JANUARY 1998

!     PURPOSE.
!     --------

!       PREPARE THE DATA FOR OPTIMAL INTERPOLATION.

!**   INTERFACE.
!     ----------

!       *CALL* *GRFIELD (WHME, USME, LDFAILD)*

!       *WHME*   REAL    OUTPUT FIELD OF MEASURED WAVE HEIGHTS.
!       *USME*   REAL    OUTPUT FIELD OF MEASURED SURFACE WINDS.
!       *NSTART* INTEGER INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!       *NEND*   INTEGER INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN


!     METHOD.
!     -------

!        READ PRE PREPARED SATELLITE DATA FROM FILE AND DISTRIBUTE TO
!        OTHER PES.

!      EXTERNALS.
!      ----------

!         MPBCASTGRDFLD  - BROADCAST 2-D DATA FIELD FROM 1 PE TO ALL
!                          OTHERS.

!      MODIFICATONS.
!      -------------
!         NONE.

!-----------------------------------------------------------------------


      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWICE   , ONLY : ICEMASK
      USE YOWMAP   , ONLY : IXLG     ,KXLT
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,KTAG
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWTEST  , ONLY : IU06
      USE YOWINTP  , ONLY : CDG



! ----------------------------------------------------------------------

      LOGICAL LLFAILED, LLOPENED
      INTEGER IREAD
      CHARACTER* 3 CLFID
      CHARACTER*12 CTIME
      REAL,   DIMENSION(NGX,NGY) :: WHME, USME
      INTEGER,DIMENSION(NPROC)   :: NSTART, NEND


!*     VARIABLE     TYPE      PURPOSE.
!      --------     ----      --------

!     *WHME*        REAL      GRID WAVE HEIGHT.
!     *USME*        REAL      GRID FRICTION VELOCITY.
!     *IREAD*       INTGER    PROCESSOR WHICH WILL READ THE INPUT FILES.
!     *CLFID        CHARACTER THREE CHARACTER FILE IDENTIFIER.

      DATA CLFID, IREAD /'RFL', 1/

! ----------------------------------------------------------------------

!*    1. OPEN FILE AND READ.
!        -------------------

      IF(LMESSPASS .AND. IRANK .EQ. IREAD .OR. .NOT.LMESSPASS) THEN
        IFIELD = 0
        LLFAILED = .FALSE.
        IFAIL = 1

!*       FIND SPARE UNIT TO ATTACH INPUT FILE.

        UNIT: DO IUME=33,99,1
          INQUIRE ( UNIT=IUME, OPENED=LLOPENED)
          IF ( .NOT. LLOPENED) THEN
            WRITE(IU06,'("  GRFIELD : TAKES UNIT ", I2,
     &       " AS ALTIMETER DATA INPUT UNIT")') IUME
            EXIT UNIT
          ENDIF
        ENDDO UNIT

!*       IF FOUND: SEEK INPUT FILE, READ HEADER AND FIELDS.
!*                 IFIELD MAY BE MISSING INDICATING THERE ARE NO FIELDS.

        IF (IUME.LT.100) THEN
          CALL CONFILE (IU06, IUME, CDTPRO, CLFID, IFAIL)
          IF(IFAIL.EQ.0) THEN
            READ (IUME,ERR=9000) CTIME, IFIELD
            READ (IUME,ERR=9000) WHME, USME
            CLOSE (IUME)
          ELSE
            LLFAILED = .TRUE.
            WRITE(IU06,'("WARNING IN GRFIELD:  W A R N I N G  !!!!!!")')
            WRITE(IU06,'("NO FILES FOUND"                            )')
            WRITE(IU06,'("WARNING IN GRFIELD:  W A R N I N G  !!!!!!")')
            CALL FLUSH(IU06)
          ENDIF
        ELSE

!*       IF NOT FOUND: ISSUE WARNING MESSAGE.

          LLFAILED = .TRUE.
          WRITE(IU06,'("WARNING IN GRFIELD:  W A R N I N G  !!!!!!!!")')
          WRITE(IU06,'("NO SPARE UNIT TO READ SATELLITE DATA FOUND"  )')
          WRITE(IU06,'("WARNING IN GRFIELD:  W A R N I N G  !!!!!!!!")')
          CALL FLUSH(IU06)
        ENDIF

!*       IF ALL ATTEMPTS OF THE NEW METHOD FAILED TRY THE OLD CALL TO
!*       GRDATA.

 9000   IF ( LLFAILED .OR. IFIELD .EQ. 0 ) THEN
          WRITE(IU06,'("WARNING IN GRFIELD:    W A R N I N G  !!!!!!")')
          WRITE(IU06,'("NO GRIDDED FIELDS FOUND IN FILE"             )')
          WRITE(IU06,'("TRY THE OLD CALL TO GRDATA                  ")')
          WRITE(IU06,'("WARNING IN GRFIELD:    W A R N I N G  !!!!!!")')
          CALL FLUSH(IU06)
          IF (IFAIL .EQ. 0) CLOSE (IUME)
          CALL GRDATA (WHME, USME, NSTART, NEND) 
        ENDIF
      ENDIF

! ----------------------------------------------------------------------

!*    2. SEND FIELDS TO ALL OTHER PES.
!        -----------------------------

      IF(LMESSPASS) THEN
        CALL MPBCASTGRDFLD(IREAD,KTAG,NGX,NGY,WHME)
        KTAG=KTAG+1
        CALL MPBCASTGRDFLD(IREAD,KTAG,NGX,NGY,USME)
        KTAG=KTAG+1
      ENDIF

! ----------------------------------------------------------------------

!     3. APPLY SCREENING OF DATA BASED ON MODEL FIELD
!        --------------------------------------------

!        EXCLUDE POINTS WHICH ARE ON MODEL ICE

      IF(LMESSPASS) THEN
        DO IJ = NSTART(1), NEND(NPROC)
          IF (ICEMASK(IJ,1).EQ.0) THEN
            IX = IXLG(IJ,1)
            IY = NGY- KXLT(IJ,1) +1
            USME(IX,IY)=-3.
            WHME(IX,IY)=-3.
          ENDIF
        ENDDO
      ELSE
        DO IG = 1, IGL
          DO IJ = IJS(IG),IJL(IG)
            IF (ICEMASK(IJ,IG).EQ.0) THEN
              IX = IXLG(IJ,IG)
              IY = NGY- KXLT(IJ,IG) +1
              USME(IX,IY)=-3.
              WHME(IX,IY)=-3.
            ENDIF
          ENDDO
        ENDDO
      ENDIF

!*       OBTAIN FRICTION VELOCITIES FROM OBSERVED WINDS AND THE DRAG.

      FV: DO J=1,NGY
        DO I=1,NGX
          IF (USME(I,J) .GT. 0) USME(I,J) = SQRT(CDG(I,J)) * USME(I,J)
          IF (CDG(I,J) .LT. 0.000) THEN
            USME(I,J)=-3.
            WHME(I,J)=-3.
          ENDIF
        ENDDO
      ENDDO FV

      RETURN
      END SUBROUTINE GRFIELD
