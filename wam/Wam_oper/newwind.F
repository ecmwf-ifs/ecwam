      SUBROUTINE NEWWIND (IJS,IJL,IG,IGL,CDATE,CDATEWH,NEWREAD,NEWFILE,
     &                    U10OLD,THWOLD,U10NEW,THWNEW,
     &                    ROAIRO, ROAIRN, ZIDLOLD,ZIDLNEW,
     &                    CICOVER, CITHICK, CIWA)
! ----------------------------------------------------------------------

!**** *NEWWIND* - HANDLING OF WINDS OUTSIDE MULTITASKED AREA    

!     P.A.E.M. JANSSEN  KNMI/ECMWF  SEPTEMBER 1994
!     J. BIDLOT         ECMWF       FEBRUARY 1996  MESSAGE PASSING 
!     S. ABDALLA        ECMWF       OCTOBER 2001   INCLUSION OF AIR
!                                                  DENSITY AND Zi/L
!     J. BIDLOT         ECMWF       AUGUST 2008: CLEAN-UP 

!*    PURPOSE.
!     --------

!       GETS FORCING FIELDS WHEN NEEDED.                                   

!**   INTERFACE.
!     ----------

!       *CALL* *NEWWIND (IJS, IJL, IG, IGL, CDATE, NEWREAD, NEWFILE,
!     1                  U10OLD,THWOLD,U10NEW,THWNEW, 
!     2                  ROAIRO, ROAIRN, ZIDLOLD,ZIDLNEW,
!     3                  CICOVER, CITHICK, CIWA)*
!      *IJS*     - INDEX OF FIRST GRIDPOINT
!      *IJL*     - INDEX OF LAST GRIDPOINT
!      *IG*      - BLOCK NUMBER
!      *IGL*     - NUMBER OF BLOCKS
!      *CDATE*   - START DATE OF SOURCE FUNCTION INTEGRATION
!      *NEWREAD* - TRUE IF NEW WINDS HAVE BEEN READ
!      *NEWFILE* - TRUE IF NEW WIND FILE HAS BEEN OPENED
!      *U10NEW*  - NEW WIND SPEED IN M/S.
!      *U10OLD*  - INTERMEDIATE STORAGE OF MODULUS OF WIND
!                  VELOCITY.
!      *THWNEW*  - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!      *THWOLD*  - INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
!                  WIND VELOCITY.
!      *ROAIRN*  - AIR DENSITY IN KG/M3.
!      *ROAIRO*  - INTERMEDIATE STORAGE OF AIR DENSITY.
!      *ZIDLNEW* - Zi/L (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!      *ZIDLOLD* - INTERMEDIATE STORAGE OF Zi/L.
!      *CICOVER* - SEA ICE COVER. 
!      *CITHICK* - SEA ICE THICKNESS. 
!      *CIWA*    - SEA ICE WAVE ATTENUATION FACTOR.

!     METHOD.
!     -------

!       READ NEW WINDS AND ASSOCIATED FIELDS WHEN NEEDED.

!     EXTERNALS.
!     ---------

!       *CREWFN*    - CREATE A WIND FILE NAME.
!       *INCDATE*   - UPDATE DATE TIME GROUP.

!     REFERENCE.
!     ----------

!       NONE

! ----------------------------------------------------------------------

!*    *PARAMETER*  FOR ARRAY DIMENSIONS.

      USE YOWPARAM , ONLY : NFRE     ,NBLO
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWSTAT  , ONLY : IDELWO   ,NWAM_BLKS
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWIND  , ONLY : CDA      ,CDAWIFL  ,CDATEFL  ,CDTNEXT  ,
     &            IIG      ,NSTORE   ,USNEXT   ,DSNEXT   ,ADSNEXT  ,
     &            ZIDLNEXT ,CICOVERNEXT, CITHICKNEXT

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK


! ----------------------------------------------------------------------
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 
      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10OLD,THWOLD,ROAIRO,ZIDLOLD
      REAL,DIMENSION(NINF:NSUP) :: U10NEW,THWNEW ,ROAIRN,ZIDLNEW
      REAL,DIMENSION(NINF:NSUP,NBLO) :: CICOVER, CITHICK 
      REAL,DIMENSION(NINF:NSUP,NFRE,NBLO) :: CIWA

! ----------------------------------------------------------------------

! Mod for OPENMP
      INTEGER JKGLO,KIJS,KIJL,NPROMA
! End Mod for OPENMP
      REAL ZHOOK_HANDLE

      LOGICAL NEWREAD, NEWFILE
      CHARACTER*14 CDATEWH, CDATE

      DATA ISTORE /0/

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('NEWWIND',0,ZHOOK_HANDLE)


!*    1. WINDS ARE TAKEN FROM INTERMEDIATE STORAGE.
!        ------------------------------------------
      IF (CDATE.LT.CDATEWH) THEN

        NPROMA=(IJL-IJS+1)/NWAM_BLKS+1
        CALL GSTATS(1492,0)
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ)
        DO JKGLO=IJS,IJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL)
          DO IJ = KIJS, KIJL
            U10NEW(IJ) = U10OLD(IJ,IG)
            THWNEW(IJ) = THWOLD(IJ,IG)
            ROAIRN(IJ) = ROAIRO(IJ,IG)
            ZIDLNEW(IJ) = ZIDLOLD(IJ,IG)
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1492,1)

      ELSE
!*    2. NEW WIND INPUT.
!        ---------------
        IF (CDATE.GE.CDATEFL) THEN
          IF (IG.EQ.1) THEN
            NEWFILE = .TRUE.
            ISTORE=1
            IF (ITEST.GE.2) THEN
              WRITE(IU06,*) '      SUB. NEWWIND: NEW WIND ',
     &         ' AT CDATE = ', CDATE
            ENDIF
          ENDIF
        ENDIF

!*    2.2 NEW WINDS ARE READ IN.
!         ----------------------
        IF(ISTORE.GT.NSTORE) THEN
          WRITE(IU06,*) ' ********************************************'
          WRITE(IU06,*) '      FATAL ERROR IN SUB. NEWWIND:           '
          WRITE(IU06,*) '      ISTORE > NSTORE !!!                    '
          WRITE(IU06,*) '      ISTORE = ',ISTORE 
          WRITE(IU06,*) '      NSTORE = ',NSTORE 
          WRITE(IU06,*) ' ********************************************'
          CALL ABORT1
        ENDIF

        CDA=CDTNEXT(ISTORE)
        IIG=IG

        NPROMA=(IJL-IJS+1)/NWAM_BLKS+1
        CALL GSTATS(1492,0)
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ)
        DO JKGLO=IJS,IJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL)
          DO IJ = KIJS, KIJL
            U10NEW(IJ)=USNEXT(IJ,ISTORE,IIG)
            THWNEW(IJ)=DSNEXT(IJ,ISTORE,IIG)
            ROAIRN(IJ)=ADSNEXT(IJ,ISTORE,IIG)
            ZIDLNEW(IJ)=ZIDLNEXT(IJ,ISTORE,IIG)
            CICOVER(IJ,IG)=CICOVERNEXT(IJ,ISTORE,IIG)
            CITHICK(IJ,IG)=CITHICKNEXT(IJ,ISTORE,IIG)
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1492,1)

        ISTORE=ISTORE+1

        CALL INCDATE(CDATEWH,IDELWO)
        NEWREAD = .TRUE.   

!       UPDATE THE SEA ICE REDUCTION FACTOR
        CALL CIREDUCE (CICOVER,CITHICK,CIWA)

      ENDIF   

      IF (LHOOK) CALL DR_HOOK('NEWWIND',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE NEWWIND
