      SUBROUTINE NEWWIND (IJS,IJL,IG,IGL,CDATE,CDATEWH,NEWREAD,NEWFILE,
     &                    U10OLD,THWOLD,U10NEW,THWNEW
     &                  , ROAIRO, ROAIRN, ZIDLOLD,ZIDLNEW)
! ----------------------------------------------------------------------

!**** *NEWWIND* - HANDLING OF WINDS OUTSIDE MULTITASKED AREA    

!     P.A.E.M. JANSSEN  KNMI/ECMWF  SEPTEMBER 1994
!     J. BIDLOT         ECMWF       FEBRUARY 1996  MESSAGE PASSING 
!     S. ABDALLA        ECMWF       OCTOBER 2001   INCLUSION OF AIR
!                                                  DENSITY AND Zi/L

!*    PURPOSE.
!     --------

!       READ WINDS WHEN NEEDED.                                   

!**   INTERFACE.
!     ----------

!       *CALL* *NEWWIND (IJS, IJL, IG, IGL, CDATE, NEWREAD, NEWFILE,
!     1                  U10OLD,THWOLD,U10NEW,THWNEW, 
!     2                  ROAIRO, ROAIRN, ZIDLOLD,ZIDLNEW)*
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

!     METHOD.
!     -------

!       READ NEW WINDS WHEN NEEDED.

!     EXTERNALS.
!     ---------

!       *CREWFN*    - CREATE A WIND FILE NAME.
!       *INCDATE*   - UPDATE DATE TIME GROUP.

!     REFERENCE.
!     ----------

!       NONE

! ----------------------------------------------------------------------

!*    *PARAMETER*  FOR ARRAY DIMENSIONS.

      USE YOWPARAM , ONLY : NBLO
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWSTAT  , ONLY : IDELWO
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWUNIT  , ONLY : IU17     ,IU18
      USE YOWWIND  , ONLY : CDA      ,CDAWIFL  ,CDATEFL  ,CDTNEXT  ,
     &            IIG      ,USNEXT   ,DSNEXT   ,ADSNEXT  ,ZIDLNEXT

! ----------------------------------------------------------------------
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 
      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10OLD,THWOLD, ROAIRO,ZIDLOLD
      REAL,DIMENSION(NINF:NSUP) :: U10NEW,THWNEW     , ROAIRN,ZIDLNEW

! ----------------------------------------------------------------------

      LOGICAL NEWREAD, NEWFILE

      CHARACTER*12 CDATEWH, CDATE

! ----------------------------------------------------------------------

!*    1. INITIALISATION.
!        ---------------


!*    1.1 INTERNAL TIME COUNTER AND WIND FILE UNIT.
!         -----------------------------------------

      IDNIN = IU17    

!*    1.2 WINDS ARE TAKEN FROM INTERMEDIATE STORAGE.
!         ------------------------------------------

      IF (CDATE.LT.CDATEWH) THEN
        DO IJ=IJS,IJL
          U10NEW(IJ) = U10OLD(IJ,IG)
          THWNEW(IJ) = THWOLD(IJ,IG)
          ROAIRN(IJ) = ROAIRO(IJ,IG)
          ZIDLNEW(IJ) = ZIDLOLD(IJ,IG)
        ENDDO
      ENDIF


! --------------------------------------------------------------------- 


!*    2. NEW WIND INPUT.
!        ---------------


      IF (CDATE.GE.CDATEWH) THEN

!*    2.1 WINDS HAVE TO BE READ FROM FILE.
!         --------------------------------

        IF (CDATE.GE.CDATEFL) THEN

!*    2.1.1 A NEW WIND FILE HAS TO BE USED.
!           -------------------------------

          IF (IG.EQ.1) THEN

!*    2.1.1.1 A NEW WIND FILE HAS TO BE ACCESSED.
!             -----------------------------------

!!!               CALL CREWFN (IU18, CDAWIFL)
            NEWFILE = .TRUE.
            IF (ITEST.GE.2) THEN
              WRITE(IU06,*) '      SUB. NEWWIND: NEW WIND FILE',
     &         ' AT CDATE = ', CDATE
              WRITE(IU06,*) '     WINDFILE DATE IS CDAWIFL ',
     &         CDAWIFL, ' UNIT IS IU18 = ',IU18
            ENDIF
          ENDIF
          IDNIN = IU18
        ENDIF

!*    2.2 NEW WINDS ARE READ IN.
!         ----------------------

        IF (ITEST.GE.2) THEN
          IF (ITESTB.GE.IG) THEN
            WRITE(IU06,*) '       SUB. NEWWIND: READ WIND ',
     &       'IDNIN = ', IDNIN, ' AT CDATE = ', CDATE
          ENDIF
        ENDIF
!!!!         READ(IDNIN) CDA, IIG, (U10NEW(L),L=IJS,IJL),
!!!!     1                         (THWNEW(L),L=IJS,IJL)
!!!!
!!!! read from file instead

        CDA=CDTNEXT
        IIG=IG
        DO L=IJS,IJL
          U10NEW(L)=USNEXT(L,IIG)
          THWNEW(L)=DSNEXT(L,IIG)
          ROAIRN(L)=ADSNEXT(L,IIG)
          ZIDLNEW(L)=ZIDLNEXT(L,IIG)
        ENDDO


        CALL INCDATE(CDATEWH,IDELWO)
        NEWREAD = .TRUE.   
      ENDIF   

      RETURN
      END SUBROUTINE NEWWIND
