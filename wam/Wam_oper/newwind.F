      SUBROUTINE NEWWIND (IJS,IJL,IG,IGL,CDATE,CDATEWH,NEWREAD,NEWFILE,
     &                    U10OLD,THWOLD,U10NEW,THWNEW,
     &                    USOLD, USNEW,
     &                    ROAIRO, ROAIRN, ZIDLOLD,ZIDLNEW,
     &                    CICOVER, CITHICK, CIWA,
     &                    TAUW)
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
!                        U10OLD,THWOLD,U10NEW,THWNEW, 
!                        USOLD, USNEW,
!                        ROAIRO, ROAIRN, ZIDLOLD,ZIDLNEW,
!                        CICOVER, CITHICK, CIWA,
!                        TAUW)
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
!      *USOLD*   - INTERMEDIATE STORAGE OF FRICTIOn VELOCITY
!      *USNEW*   - NEW FRICTIOn VELOCITY
!      *ROAIRN*  - AIR DENSITY IN KG/M3.
!      *ROAIRO*  - INTERMEDIATE STORAGE OF AIR DENSITY.
!      *ZIDLNEW* - Zi/L (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!      *ZIDLOLD* - INTERMEDIATE STORAGE OF Zi/L.
!      *CICOVER* - SEA ICE COVER. 
!      *CITHICK* - SEA ICE THICKNESS. 
!      *CIWA*    - SEA ICE WAVE ATTENUATION FACTOR.
!      *TAUW*    - WAVE STRESS

!     METHOD.
!     -------

!       READ NEW WINDS AND ASSOCIATED FIELDS WHEN NEEDED.

!     EXTERNALS.
!     ---------

!       *INCDATE*   - UPDATE DATE TIME GROUP.

!     REFERENCE.
!     ----------

!       NONE

! ----------------------------------------------------------------------

!*    *PARAMETER*  FOR ARRAY DIMENSIONS.

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWPARAM , ONLY : NFRE     ,NBLO
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWSTAT  , ONLY : IDELWO   ,NPROMA_WAM
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWIND  , ONLY : CDA      ,CDAWIFL  ,CDATEFL  ,CDTNEXT  ,
     &            NSTORE   ,FF_NEXT  ,WSPMIN_RESET_TAUW
      USE YOWWNDG  , ONLY : ICODE    ,ICODE_CPL
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL
      INTEGER, INTENT(IN) :: IG, IGL

      REAL,DIMENSION(NINF:NSUP,NBLO), INTENT(IN) :: U10OLD, THWOLD
      REAL,DIMENSION(NINF:NSUP,NBLO), INTENT(IN) :: USOLD
      REAL,DIMENSION(NINF:NSUP,NBLO), INTENT(IN) :: ROAIRO, ZIDLOLD
      REAL,DIMENSION(NINF:NSUP,NBLO), INTENT(INOUT) :: TAUW
      REAL,DIMENSION(NINF:NSUP), INTENT(INOUT) :: U10NEW, THWNEW
      REAL,DIMENSION(NINF:NSUP), INTENT(INOUT) :: USNEW
      REAL,DIMENSION(NINF:NSUP), INTENT(INOUT) :: ROAIRN, ZIDLNEW
      REAL,DIMENSION(NINF:NSUP,NBLO), INTENT(INOUT) :: CICOVER, CITHICK 
      REAL,DIMENSION(NINF:NSUP,NFRE,NBLO), INTENT(INOUT) :: CIWA

      CHARACTER(LEN=14), INTENT(IN) :: CDATE
      CHARACTER(LEN=14), INTENT(INOUT) :: CDATEWH

      LOGICAL, INTENT(INOUT) :: NEWREAD, NEWFILE


      INTEGER :: ICODE_WND
      INTEGER :: ISTORE
      INTEGER :: JKGLO, KIJS, KIJL, NPROMA, IJ

      REAL :: ZHOOK_HANDLE

      DATA ISTORE /0/

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('NEWWIND',0,ZHOOK_HANDLE)
#endif

      NPROMA=NPROMA_WAM

      IF(LWCOU) THEN
        ICODE_WND = ICODE
      ELSE
        ICODE_WND = ICODE_CPL
      ENDIF


!*    1. WINDS ARE TAKEN FROM INTERMEDIATE STORAGE.
!        ------------------------------------------
      IF (CDATE.LT.CDATEWH) THEN
        CALL GSTATS(1492,0)
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ)
        DO JKGLO=IJS,IJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL)
          DO IJ = KIJS, KIJL
            U10NEW(IJ) = U10OLD(IJ,IG)
            USNEW(IJ) = USOLD(IJ,IG)
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

        CALL GSTATS(1492,0)
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ)
        DO JKGLO=IJS,IJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL)
          IF (ICODE_WND == 3 ) THEN
            DO IJ = KIJS, KIJL
              U10NEW(IJ)=FF_NEXT(IJ,ISTORE)%WSWAVE
            ENDDO
          ELSE
            DO IJ = KIJS, KIJL
!              USNEW(IJ)=FF_NEXT(IJ,ISTORE)%USTAR
            ENDDO
          ENDIF
          DO IJ = KIJS, KIJL
            THWNEW(IJ)=FF_NEXT(IJ,ISTORE)%WDWAVE
            ROAIRN(IJ)=FF_NEXT(IJ,ISTORE)%AIRD
            ZIDLNEW(IJ)=FF_NEXT(IJ,ISTORE)%ZIDL
            CICOVER(IJ,IG)=FF_NEXT(IJ,ISTORE)%CIFR
            CITHICK(IJ,IG)=FF_NEXT(IJ,ISTORE)%CITH

! adapt first estimate of wave induced stress for low winds
            IF(U10NEW(IJ).LT.WSPMIN_RESET_TAUW) TAUW(IJ,IG)=0.0

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

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('NEWWIND',1,ZHOOK_HANDLE)
#endif
      RETURN
      END SUBROUTINE NEWWIND
