      SUBROUTINE WAMWND (U10, THW, ADS, ZIDL, CITH, IJS, IJL, IG, LWCUR)

! ----------------------------------------------------------------------

!     MODIFIED  S. ABDALLA OCTOBER 2001  INCLUSION OF AIR DENSITY & Zi/L

!     MODIFIED J. BIDLOT 2008 :: INCLUDE SURFACE CURRENTS FOR
!                                RELATIVE WIND CALCULATION. 
!                                USE OF NEUTRAL WIND SPEED FROM WAM.

!**** *WAMWND* - TRANSFORMS INPUT FORCING TO BLOCKED WAM POINTS.


!*    PURPOSE.
!     --------

!       CONVERTS THE INTERPOLATED INPUT FIELDS TO WAM BLOCKS FOR ALL
!       POINTS IN THE GRID ON A PE, EXCEPT FOR U and V CURRENTS !!!
!       SEE *GETCURR* (THEY ARE NEEDED OVER THE GRID POINT HALO).
!       U and V ARE HOWEVER NEEDED HERE !!!.

!**   INTERFACE.
!     ----------

!       *CALL WAMWND (U10, THW, ADS, ZIDL, CITH, IJS, IJL, IG, LWCUR)*
!          *U10*  - INTERPOLATED WINDS AT ALL POINTS AND BLOCKS.
!          *THW*  - INTERPOLATED WIND DIRECTION AT ALL POINTS.
!          *ADS*  - INTERPOLATED AIR DENSITY AT ALL POINTS.
!          *ZIDL* - INTERPOLATED Zi/L AT ALL POINTS
!                   (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!          *CITH* - SEA ICE THICKNESS. 
!          *IJS*    - INDEX OF FIRST GRIDPOINT
!          *IJL*    - INDEX OF LAST GRIDPOINT
!          *IG*     - BLOCK NUMBER
!          *LWCUR*  - LOGICAL INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS


!     METHOD.
!     -------

!       THE INTERPOLATED VALUES ARE TRANSFORMED TO
!       MAGNITUDE AND DIRECTION. INPUT MAY BE WIND IN 10M HEIGHT ,
!       SURFACE WINDS OR FRICTION VELOCITIES. THE INPUT GRID HAS TO BE
!       A LATITUDE/LONGITUDE GRID EITHER PERIODIC OR NON PERIODIC.

!     EXTERNALS.
!     ----------

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWCOU    , XKAPPA
      USE YOWCURR  , ONLY : U        ,V
      USE YOWMAP   , ONLY : IFROMIJ  ,JFROMIJ
      USE YOWPCONS , ONLY : G        ,ZPI      ,ZMISS
      USE YOWSTAT  , ONLY : IREFRA   ,LRELWIND
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWIND  , ONLY : WSPMIN   ,FIELDG   ,LLWSWAVE ,LLWDWAVE ,
     &            RWFAC
      USE YOWWNDG  , ONLY : ICODE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: IJ, IX, JY, IJS, IJL, IG

      REAL, PARAMETER :: ALPHACH = 0.0185
      REAL :: UU, VV, WSPEED, RESCALE, USTAR, Z0, CD
      REAL :: ZHOOK_HANDLE
      REAL, DIMENSION (IJS:IJL) :: U10,THW,ADS,ZIDL,CITH

      LOGICAL :: LWCUR

! ----------------------------------------------------------------------

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('WAMWND',0,ZHOOK_HANDLE)
# endif

!*    2. TRANSFORM GRIDDED WIND INPUT INTO BLOCK
!*       AND TRANSFORM TO MAGNITUDE AND DIRECTION.
!*       AND SET LOWER WIND SPEED TO WSPMIN m/s
!        ----------------------------------------

      IF(LLWSWAVE .AND. LLWDWAVE) THEN
        DO IJ = IJS,IJL
          IX = IFROMIJ(IJ,IG)
          JY = JFROMIJ(IJ,IG)

!         USE NEUTRAL WIND SPEED AND DIRECTION FROM A PREVIOUS WAM RUN.
!         -------------------------------------------------------------
          U10(IJ) = FIELDG(IX,JY)%WSWAVE
          THW(IJ) = FIELDG(IX,JY)%WDWAVE 

!         THERE MIGHT BE POINTS FOR WHICH NO WAM VALUES ARE AVAILABLE
!         USE U and V FROM ATMOSPHERE INSTEAD
          IF(U10(IJ).LE.0.) THEN
            UU = FIELDG(IX,JY)%UWND
            VV = FIELDG(IX,JY)%VWND
            WSPEED = SQRT(UU**2 + VV**2)
            IF(WSPEED.GT.0.) THEN
              U10(IJ) = WSPEED
              THW(IJ) = ATAN2(UU,VV)
            ELSE
              U10(IJ) = 0.
              THW(IJ) = 0.
            ENDIF
          ENDIF
        ENDDO

        IF(LRELWIND .AND. (IREFRA.EQ.2 .OR. IREFRA.EQ.3) ) THEN
          DO IJ = IJS,IJL
            IX = IFROMIJ(IJ,IG)
            JY = JFROMIJ(IJ,IG)

            UU = U10(IJ)*SIN(THW(IJ))
            VV = U10(IJ)*COS(THW(IJ))
            UU = UU - RWFAC*U(IJ,IG)
            VV = VV - RWFAC*V(IJ,IG)
            WSPEED = SQRT(UU**2 + VV**2)
            IF(WSPEED.GT.0.) THEN
              U10(IJ) = WSPEED
              THW(IJ) = ATAN2(UU,VV)
            ELSE
              U10(IJ) = 0.
              THW(IJ) = 0.
            ENDIF
          ENDDO
        ENDIF

      ELSE

        DO IJ = IJS,IJL
          IX = IFROMIJ(IJ,IG)
          JY = JFROMIJ(IJ,IG)

          UU = FIELDG(IX,JY)%UWND
          VV = FIELDG(IX,JY)%VWND

!         USE NEUTRAL WIND SPEED FROM A PREVIOUS WAM RUN.
!         ----------------------------------------------
          IF(LLWSWAVE) THEN
            IF(FIELDG(IX,JY)%WSWAVE.NE.ZMISS .AND.
     &         FIELDG(IX,JY)%WSWAVE.GT.0. ) THEN
              WSPEED = SQRT(UU**2 + VV**2)
              IF(WSPEED.GT.0.) THEN
                RESCALE=FIELDG(IX,JY)%WSWAVE/WSPEED
                UU = UU * RESCALE 
                VV = VV * RESCALE
              ENDIF
            ENDIF
          ENDIF

!         CORRECT FOR RELATIVE WINDS WITH RESPECT TO THE SURFACE CURRENTS.
!         ----------------------------------------------------------------

          IF(LWCOU) THEN
            IF(LWCUR) THEN
!           In coupled experiments the winds are already relative to the currents
              IF(.NOT.LRELWIND) THEN
                UU = UU + RWFAC*U(IJ,IG)
                VV = VV + RWFAC*V(IJ,IG)
              ENDIF
            ENDIF
          ELSE
            IF(LRELWIND) THEN
              IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN
                UU = UU - RWFAC*U(IJ,IG)
                VV = VV - RWFAC*V(IJ,IG)
              ENDIF
            ENDIF
          ENDIF

          U10(IJ) = SQRT(UU**2 + VV**2)
          IF (U10(IJ).NE.0.) THEN 
            THW(IJ) = ATAN2(UU,VV)
          ELSE
            THW(IJ) = 0.
          ENDIF

        ENDDO
      ENDIF

      DO IJ = IJS,IJL
        IF (THW(IJ).LT.0.) THW(IJ) = THW(IJ) + ZPI
!       IMPOSE A MINIMUM WIND SPEED.
!       ----------------------------
        U10(IJ) = MAX(U10(IJ),WSPMIN)
      ENDDO

      DO IJ = IJS,IJL
        IX = IFROMIJ(IJ,IG)
        JY = JFROMIJ(IJ,IG)
        ADS(IJ) = FIELDG(IX,JY)%AIRD
        ZIDL(IJ)= FIELDG(IX,JY)%ZIDL
        CITH(IJ)= FIELDG(IX,JY)%CITH
      ENDDO

! ----------------------------------------------------------------------

!*    3. PROCESS WINDS ACCORDING TO TYPE
!        NOTHING TO DO FOR WIND SPEED U10 (ICODE = 3).
!        ---------------------------------------------

      IF (ICODE.EQ.1) THEN

!*    3.2  INPUT IS FRICTION VELOCITY.
!          ---------------------------

        DO IJ = IJS,IJL
          USTAR = MAX(0.01,U10(IJ))
          Z0  = ALPHACH/G*USTAR**2
          CD  = XKAPPA/LOG(10./Z0)
          U10(IJ) = USTAR/CD
          U10(IJ) = MAX(U10(IJ),WSPMIN)
        ENDDO

      ELSEIF (ICODE.EQ.2) THEN

!*    3.3 INPUT WINDS ARE SURFACE STRESSES.
!         ---------------------------------

        DO IJ = IJS,IJL
          USTAR = MAX (0.01, SQRT(U10(IJ)/ADS(IJ)))
          Z0  = ALPHACH/G*USTAR**2
          CD  = XKAPPA/LOG(10./Z0)
          U10(IJ) = USTAR/CD
          U10(IJ) = MAX(U10(IJ),WSPMIN)
        ENDDO
      ENDIF

!*    4. TEST OUTPUT OF WAVE MODEL BLOCKS
!        ---------------------------------

      IF (ITEST.GE.3) THEN
        WRITE (IU06,*) ' '
        WRITE (IU06,*) '      SUB. WAMWND:',
     &   ' INPUT FORCING FIELDS CONVERTED TO BLOCKS'
        CALL FLUSH(IU06)
      ENDIF
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('WAMWND',1,ZHOOK_HANDLE)
# endif
      RETURN
      END SUBROUTINE WAMWND
