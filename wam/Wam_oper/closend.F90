      SUBROUTINE CLOSEND (IJS,IJL,CDATE,CDATEWH,NEWREAD,NEWFILE,        &
     &                    U10OLD,THWOLD,ROAIRO,ZIDLOLD,                 &
     &                    U10NEW,THWNEW,ROAIRN,ZIDLNEW) 

! ----------------------------------------------------------------------

!**** *CLOSEND* - HANDLING OF WINDS OUTSIDE MULTITASKED AREA    

!     P.A.E.M. JANSSEN  KNMI/ECMWF  SEPTEMBER 1994
!     J. BIDLOT         ECMWF       FEBRUARY  1996  MESSAGE PASSING
!     S. ABDALLA        ECMWF       OCTOBER   1996  AIR DENSITY AND Zi/L

!*    PURPOSE.
!     --------

!       READ WINDS WHEN NEEDED.                                   

!**   INTERFACE.
!     ----------

!     *CALL* *CLOSEND*(IJS, IJL, CDATE,CDATEWH,NEWREAD,NEWFILE,
!    1                 U10OLD,THWOLD,ROAIRO,ZIDLOLD,U10NEW,THWNEW,
!    2                 ROAIRN,ZIDLNEW) 
!      *IJS*    - INDEX OF FIRST GRIDPOINT
!      *IJL*    - INDEX OF LAST GRIDPOINT
!      *NEWREAD* - TRUE IF NEW WINDS HAVE BEEN READ
!      *NEWFILE* - TRUE IF NEW WIND FILE HAS BEEN OPENED
!      *U10NEW*    NEW WIND SPEED IN M/S.
!      *U10OLD*    INTERMEDIATE STORAGE OF MODULUS OF WIND
!                  VELOCITY.
!      *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!      *THWOLD*    INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
!                  WIND VELOCITY.
!      *ROAIRN*    AIR DENSITY IN KG/M3.
!      *ROAIRO*    INTERMEDIATE STORAGE OF AIR DENSITY.
!      *ZIDLNEW*   Zi/L (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!      *ZIDLOLD*   INTERMEDIATE STORAGE OF Zi/L.

!     METHOD.
!     -------

!       COPIES NEW TO OLD WINDS WHEN NEEDED. CLOSES FILES

!     EXTERNALS.
!     ---------

!       *INCDATE*   - UPDATE DATE TIME GROUP.

!     REFERENCE.
!     ----------

!       NONE

! ----------------------------------------------------------------------

!*    *PARAMETER*  FOR ARRAY DIMENSIONS.

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWSTAT  , ONLY : IDELPRO  ,IDELWI   ,NPROMA_WAM
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      CHARACTER(LEN=14), INTENT(IN) :: CDATEWH, CDATE
      LOGICAL, INTENT(INOUT) :: NEWREAD, NEWFILE

      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(INOUT) :: U10OLD, THWOLD, ROAIRO, ZIDLOLD
      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(IN) :: U10NEW, THWNEW, ROAIRN, ZIDLNEW 

      INTEGER(KIND=JWIM) :: IJ, JKGLO, KIJS, KIJL, NPROMA, IDELWH

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      LOGICAL :: LLEX

! ----------------------------------------------------------------------

!*    1. SAVE WIND INTO INTERMEDIATE STORAGE.
!     ----------------------------------------

      IF (LHOOK) CALL DR_HOOK('CLOSEND',0,ZHOOK_HANDLE)

      IF (NEWREAD) THEN
        NPROMA=NPROMA_WAM
        CALL GSTATS(1491,0)
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,IJ)
        DO JKGLO=IJS,IJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL)
          DO IJ=KIJS,KIJL
            U10OLD(IJ) = U10NEW(IJ)
            THWOLD(IJ) = THWNEW(IJ)
            ROAIRO(IJ) = ROAIRN(IJ)
            ZIDLOLD(IJ) = ZIDLNEW(IJ)
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1491,1)
        NEWREAD = .FALSE.
        IF (ITEST.GE.2) THEN
            WRITE(IU06,*) '       SUB. CLOSEND: STORE WIND ',           &
     &       ' AT CDATE = ',CDATE
        ENDIF
      ENDIF                                     

! ----------------------------------------------------------------------

!*    2. UPDATE WIND COUNTERS IF LAST BLOCK HAS BEEN DONE.
!        -------------------------------------------------

        IF (NEWFILE) THEN
!*        UPDATE WIND FILE TIME COUNTER AND UNITS.
          NEWFILE = .FALSE.
          IDELWH = MAX(IDELWI,IDELPRO)
          CALL INCDATE(CDAWIFL,IDELWH)
          CALL INCDATE(CDATEFL,IDELWH)
          IF (ITEST.GE.2) THEN
            WRITE(IU06,*) '      SUB. CLOSEND: NEW WINDFILE DATES'
            WRITE(IU06,*) '        DATE OF NEXT WIND FILE IS ',         &
     &       'CDAWIFL = ', CDAWIFL
            WRITE(IU06,*) '        FILE WILL BE ACCESSED  AT ',         &
     &       'CDATEFL = ', CDATEFL
          ENDIF
        ENDIF

!*      UPDATE WIND FIELD COUNTER.
        CDATEWO=CDATEWH

      IF (LHOOK) CALL DR_HOOK('CLOSEND',1,ZHOOK_HANDLE)

      END SUBROUTINE CLOSEND
