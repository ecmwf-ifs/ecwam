      SUBROUTINE WAMWND (IJS, IJL,                                      &
     &                   U10, US,                                       &
     &                   THW, ADS, ZIDL, CITH,                          &
     &                   LWCUR, ICODE_WND)

! ----------------------------------------------------------------------

!     MODIFIED  S. ABDALLA OCTOBER 2001  INCLUSION OF AIR DENSITY & Zi/L

!     MODIFIED J. BIDLOT 2008 :: INCLUDE SURFACE CURRENTS FOR
!                                RELATIVE WIND CALCULATION. 
!                                USE OF NEUTRAL WIND SPEED FROM WAM

!**** *WAMWND* - TRANSFORMS INPUT FORCING TO BLOCKED WAM POINTS.


!*    PURPOSE.
!     --------

!       CONVERTS THE INTERPOLATED INPUT FIELDS TO WAM BLOCKS FOR ALL
!       POINTS IN THE GRID ON A PE, EXCEPT FOR U and V CURRENTS !!!
!       SEE *GETCURR* (THEY ARE NEEDED OVER THE GRID POINT HALO).
!       U and V ARE HOWEVER NEEDED HERE !!!.

!**   INTERFACE.
!     ----------

!       *CALL WAMWND (IJS, IJL,
!    &                U10, US,
!    &                THW, ADS, ZIDL, CITH,
!    &                LWCUR, ICODE_WND)
!          *U10*  - INTERPOLATED WINDS AT ALL POINTS AND BLOCKS.
!          *US*   - INTERPOLATED FRICTION VELOCITY
!          *THW*  - INTERPOLATED WIND DIRECTION AT ALL POINTS.
!          *ADS*  - INTERPOLATED AIR DENSITY AT ALL POINTS.
!          *ZIDL* - INTERPOLATED Zi/L AT ALL POINTS
!                   (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!          *CITH* - SEA ICE THICKNESS. 
!          *IJS*    - INDEX OF FIRST GRIDPOINT
!          *IJL*    - INDEX OF LAST GRIDPOINT
!          *LWCUR*  - LOGICAL INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS
!          *ICODE_WND* - INTEGER INDICATES WHAT IS SAVED IN FIELDG%UWND, and
!                        FIELDG%VWND
!                        AND
!                        WHAT IS UPDATED:
!                        ICODE_WND=3 THEN U10, otherwise US is UPDATED !!
!                        ++++++++++++++++++++++++++++++++++++++++++++++++


!     METHOD.
!     -------

!       THE INTERPOLATED VALUES ARE TRANSFORMED TO
!       MAGNITUDE AND DIRECTION. INPUT MAY BE WIND IN 10M HEIGHT ,
!       SURFACE WINDS OR FRICTION VELOCITIES. THE INPUT GRID HAS TO BE
!       A LATITUDE/LONGITUDE GRID EITHER PERIODIC OR NON PERIODIC.

!     EXTERNALS.
!     ----------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWCURR  , ONLY : U        ,V
      USE YOWMAP   , ONLY : IFROMIJ  ,JFROMIJ
      USE YOWPCONS , ONLY : G        ,ZPI      ,ZMISS    ,EPSUS
      USE YOWPHYS  , ONLY : XKAPPA
      USE YOWSTAT  , ONLY : IREFRA   ,LRELWIND
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWIND  , ONLY : WSPMIN   ,FIELDG   ,LLWSWAVE ,LLWDWAVE ,    &
     &            RWFAC
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM), INTENT(IN) :: ICODE_WND
      REAL(KIND=JWRB), DIMENSION (IJS:IJL), INTENT(INOUT) :: U10, US
      REAL(KIND=JWRB), DIMENSION (IJS:IJL), INTENT(OUT) :: THW, ADS, ZIDL, CITH
      LOGICAL, INTENT(IN) :: LWCUR

      INTEGER(KIND=JWIM) :: IJ, IX, JY

      REAL(KIND=JWRB) :: RESCALE
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION (IJS:IJL) :: UU, VV, WSPEED

      LOGICAL :: LCORREL

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WAMWND',0,ZHOOK_HANDLE)

!     CORRECT FOR RELATIVE WINDS WITH RESPECT TO THE SURFACE CURRENTS.
      LCORREL=.FALSE.
      IF(LWCOU) THEN
!       In coupled experiments the winds are already relative to the currents
        IF(LWCUR .AND. (.NOT.LRELWIND) ) LCORREL=.TRUE.
      ELSE
        IF(LRELWIND .AND. (IREFRA.EQ.2 .OR. IREFRA.EQ.3)) LCORREL=.TRUE.
      ENDIF

!*    2. TRANSFORM GRIDDED WIND INPUT INTO BLOCK
!        ----------------------------------------

      DO IJ = IJS,IJL
        IX = IFROMIJ(IJ,1)
        JY = JFROMIJ(IJ,1)
        UU(IJ) = FIELDG(IX,JY)%UWND
        VV(IJ) = FIELDG(IX,JY)%VWND
        ADS(IJ) = FIELDG(IX,JY)%AIRD
        ZIDL(IJ)= FIELDG(IX,JY)%ZIDL
        CITH(IJ)= FIELDG(IX,JY)%CITH
      ENDDO


!*    3. PROCESS WINDS ACCORDING TO TYPE
!        ----------------------------------------

      SELECT CASE (ICODE_WND)
      CASE(3)
!       USE NEUTRAL WIND SPEED AND DIRECTION FROM A PREVIOUS WAM RUN.
        IF(LLWSWAVE .AND. LLWDWAVE) THEN
          DO IJ = IJS,IJL
            IX = IFROMIJ(IJ,1)
            JY = JFROMIJ(IJ,1)
            U10(IJ) = FIELDG(IX,JY)%WSWAVE
            THW(IJ) = FIELDG(IX,JY)%WDWAVE 
          ENDDO

!         THERE MIGHT BE POINTS FOR WHICH NO WAM VALUES ARE AVAILABLE
!         USE U and V FROM ATMOSPHERE INSTEAD
          DO IJ = IJS,IJL
            IF(U10(IJ).LE.0.0_JWRB) THEN
              WSPEED(IJ) = SQRT(UU(IJ)**2 + VV(IJ)**2)
              IF(WSPEED(IJ).GT.0.0_JWRB) THEN
                U10(IJ) = WSPEED(IJ)
                THW(IJ) = ATAN2(UU(IJ),VV(IJ))
              ELSE
                U10(IJ) = 0.0_JWRB
                THW(IJ) = 0.0_JWRB
              ENDIF
            ENDIF
          ENDDO

!         CORRECT FOR RELATIVE WINDS WITH RESPECT TO THE SURFACE CURRENTS.
          IF(LCORREL) THEN
            DO IJ = IJS,IJL
              UU(IJ) = U10(IJ)*SIN(THW(IJ))
              VV(IJ) = U10(IJ)*COS(THW(IJ))
              UU(IJ) = UU(IJ) - RWFAC*U(IJ,1)
              VV(IJ) = VV(IJ) - RWFAC*V(IJ,1)
              WSPEED(IJ) = SQRT(UU(IJ)**2 + VV(IJ)**2)
              IF(WSPEED(IJ).GT.0.0_JWRB) THEN
                U10(IJ) = WSPEED(IJ)
                THW(IJ) = ATAN2(UU(IJ),VV(IJ))
              ELSE
                U10(IJ) = 0.0_JWRB
                THW(IJ) = 0.0_JWRB
              ENDIF
            ENDDO
          ENDIF

!       USE NEUTRAL WIND SPEED FROM A PREVIOUS WAM RUN.
!       -----------------------------------------------
        ELSE IF(LLWSWAVE) THEN

          DO IJ = IJS,IJL
            IX = IFROMIJ(IJ,1)
            JY = JFROMIJ(IJ,1)

            IF(FIELDG(IX,JY)%WSWAVE.NE.ZMISS .AND.                      &
     &         FIELDG(IX,JY)%WSWAVE.GT.0.0_JWRB ) THEN
              WSPEED(IJ) = SQRT(UU(IJ)**2 + VV(IJ)**2)
              IF(WSPEED(IJ).GT.0.0_JWRB) THEN
                RESCALE=FIELDG(IX,JY)%WSWAVE/WSPEED(IJ)
                UU(IJ) = UU(IJ) * RESCALE 
                VV(IJ) = VV(IJ) * RESCALE
              ENDIF
            ENDIF
          ENDDO

          IF(LCORREL) THEN
            DO IJ = IJS,IJL
              UU(IJ) = UU(IJ) + RWFAC*U(IJ,1)
              VV(IJ) = VV(IJ) + RWFAC*V(IJ,1)
            ENDDO
          ENDIF

          DO IJ = IJS,IJL
            U10(IJ) = SQRT(UU(IJ)**2 + VV(IJ)**2)
            IF (U10(IJ).NE.0.0_JWRB) THEN 
              THW(IJ) = ATAN2(UU(IJ),VV(IJ))
            ELSE
              THW(IJ) = 0.0_JWRB
            ENDIF
          ENDDO


        ELSE

          IF(LCORREL) THEN
            DO IJ = IJS,IJL
              UU(IJ) = UU(IJ) + RWFAC*U(IJ,1)
              VV(IJ) = VV(IJ) + RWFAC*V(IJ,1)
            ENDDO
          ENDIF

          DO IJ = IJS,IJL
            U10(IJ) = SQRT(UU(IJ)**2 + VV(IJ)**2)
            IF (U10(IJ).NE.0.0_JWRB) THEN 
              THW(IJ) = ATAN2(UU(IJ),VV(IJ))
            ELSE
              THW(IJ) = 0.0_JWRB
            ENDIF
          ENDDO

        ENDIF

!       IMPOSE A MINIMUM WIND SPEED.
        DO IJ = IJS,IJL
          U10(IJ) = MAX(U10(IJ),WSPMIN)
        ENDDO

      CASE(1)

!*    3.2  INPUT IS FRICTION VELOCITY.
!          ---------------------------

        DO IJ = IJS,IJL
          US(IJ) = SQRT(UU(IJ)**2 + VV(IJ)**2)
          IF (US(IJ).NE.0.0_JWRB) THEN 
            THW(IJ) = ATAN2(UU(IJ),VV(IJ))
          ELSE
            THW(IJ) = 0.0_JWRB
          ENDIF
          US(IJ) = MAX(US(IJ),EPSUS)
        ENDDO

      CASE(2)

!*    3.3 INPUT WINDS ARE SURFACE STRESSES.
!         ---------------------------------

        DO IJ = IJS,IJL
          US(IJ) = SQRT(UU(IJ)**2 + VV(IJ)**2)
          IF (US(IJ).NE.0.0_JWRB) THEN 
            THW(IJ) = ATAN2(UU(IJ),VV(IJ))
          ELSE
            THW(IJ) = 0.0_JWRB
          ENDIF
          US(IJ) = SQRT(MAX(US(IJ),0.0_JWRB)/MAX(ADS(IJ),1.))
          US(IJ) = MAX(US(IJ),EPSUS)
        ENDDO

      END SELECT 

      DO IJ = IJS,IJL
        IF (THW(IJ).LT.0.0_JWRB) THW(IJ) = THW(IJ) + ZPI
      ENDDO

! ----------------------------------------------------------------------

!*    4. TEST OUTPUT OF WAVE MODEL BLOCKS
!        ---------------------------------

      IF (ITEST.GE.3) THEN
        WRITE (IU06,*) ' '
        WRITE (IU06,*) '      SUB. WAMWND:',                            &
     &   ' INPUT FORCING FIELDS CONVERTED TO BLOCKS'
        CALL FLUSH(IU06)
      ENDIF

      IF (LHOOK) CALL DR_HOOK('WAMWND',1,ZHOOK_HANDLE)

      END SUBROUTINE WAMWND
