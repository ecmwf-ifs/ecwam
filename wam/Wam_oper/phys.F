      SUBROUTINE PHYS (FL3,FL1,FL,IJS,IJL,IG,NSIZE,NWORK,CDATE,CDTPRA,
     &                 U10OLD,THWOLD,USOLD,TAUW,Z0OLD,ROAIRO,ZIDLOLD,
     &                 U10NEW,THWNEW,USNEW,Z0NEW,ROAIRN,ZIDLNEW,
     &                 SL,FCONST,L1STCALL)

! ----------------------------------------------------------------------

!**** *PHYS* - SET UP PROCESS TO RUN IMPLSCH AND PROPAGS

!     P.A.E.M. JANSSEN K.N.M.I./E.C.M.W.F. SEPTEMBER 1994
!     J. BIDLOT ECMWF FEBRUARY 1996-97 MESSAGE PASSING
!     S. ABDALLA  ECMWF  OCTOBER 2001  INCLUSION OF AIR DENSITY AND Zi/L.

!*    PURPOSE.
!     --------

!       BY SPLITTING UP INNER DO LOOP IN A NUMBER OF PARTS WE DISTRIBUTE
!       THE WORK OVER A NUMBER OF PROCESSORS.

!**   INTERFACE.
!     ----------

!       *CALL* *PHYS (FL3,FL,IJS,IJL,IG,NSIZE,NWORK,CDATE,CDTPRA,
!    1                U10OLD,THWOLD,USOLD,TAUW,Z0OLD,ROAIRO,ZIDLOLD,
!    2                U10NEW,THWNEW,USNEW,Z0NEW,ROAIRN,ZIDLNEW,
!    3                SL,FCONST,L1STCALL)*
!      *FL3*    - FREQUENCY SPECTRUM(INPUT AND OUTPUT).
!      *FL1*    - FREQUENCY SPECTRUM(INPUT).
!      *IJS*    - INDEX OF FIRST GRIDPOINT
!      *IJL*    - INDEX OF LAST GRIDPOINT
!      *IG*     - BLOCK NUMBER
!      *NSIZE*  - SIZE OF SUB TASK
!      *NWORK*  - TOTAL AMOUNT OF SUB TASKS
!      *CDATE*  - DATE OF RUN
!      *CDTPRA* - PROPAGATION DATE
!      *U10NEW*    NEW WIND SPEED IN M/S.
!      *U10OLD*    INTERMEDIATE STORAGE OF MODULUS OF WIND
!                  VELOCITY.
!      *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!      *THWOLD*    INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
!                  WIND VELOCITY.
!      *USNEW*     NEW FRICTION VELOCITY IN M/S.
!      *USOLD*     INTERMEDIATE STORAGE OF MODULUS OF FRICTION
!                  VELOCITY.
!      *Z0NEW*     ROUGHNESS LENGTH IN M.
!      *Z0OLD*     INTERMEDIATE STORAGE OF ROUGHNESS LENGTH IN
!                  M.
!      *TAUW*      WAVE STRESS IN (M/S)**2
!      *ROAIRN*    AIR DENSITY IN KG/M3.
!      *ROAIRO*    INTERMEDIATE STORAGE OF AIR DENSITY.
!      *ZIDLNEW*   Zi/L (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!      *ZIDLOLD*   INTERMEDIATE STORAGE OF Zi/L.
!      *SL*        REAL      TOTAL SOURCE FUNCTION ARRAY.
!      *FCONST*    REAL      = 1 FOR PROGNOSTIC FREQUENCY BANDS.
!                            = 0 FOR DIAGNOSTIC FREQUENCY BANDS.
!      *L1STCALL*  LOGICAL  SHOULD BE FA:SE AFTER FIRST CALL TO PROPAGS.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ---------

!       *IMPLSCH*    - INTEGRATES SOURCE TERMS.     

!     REFERENCE.
!     ----------

!       NONE.                                                     


! ----------------------------------------------------------------------

!*    *PARAMETER*  FOR ARRAY DIMENSIONS.

      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMLTSK , ONLY : NCHUNK
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO     ,NIBLO

! ----------------------------------------------------------------------

      CHARACTER*14 CDATE, CDTPRA
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL,FL3,FL1,SL
      REAL,DIMENSION(NINF:NSUP,NFRE) :: FCONST
      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10OLD,THWOLD,USOLD,Z0OLD,TAUW
     &                                 ,ROAIRO,ZIDLOLD
      REAL,DIMENSION(NINF:NSUP) :: U10NEW,THWNEW,USNEW,Z0NEW
     &                            ,ROAIRN,ZIDLNEW

      LOGICAL :: L1STCALL

! ----------------------------------------------------------------------


!*    1. INITIALISATION.
!        ---------------


! DETERMINE WHAT FORM OF PARALELLISM IS USED BASED ON THE TYPE
! OF MACHINE ARCHTECTURE

      IF(LMESSPASS) THEN

        IF (CDATE.EQ.CDTPRA) THEN
          IF(NIBLO.GT.1) THEN
            CALL PROPAGS(FL1, FL3, IJS, IJL, IG, SL, FCONST, L1STCALL)
          ELSE
            FL3 = FL1
          ENDIF
        ENDIF
        CALL IMPLSCH (FL3, FL, IJS, IJL, IG,
     &   THWOLD,USOLD,TAUW,Z0OLD,    ROAIRO,ZIDLOLD,
     &   U10NEW,THWNEW,USNEW,Z0NEW,  ROAIRN,ZIDLNEW,
     &   SL,FCONST)

      ELSE

 1000   CONTINUE  

        NN = INCINT(NCHUNK,1)
        IF (NN.LE.NWORK) THEN
          MIJS = IJS + (NN-1)*NSIZE
          MIJL = MIJS + NSIZE-1
          MIJS = MIN(MIJS,IJL)
          MIJL = MIN(MIJL,IJL)


          IF (CDATE.EQ.CDTPRA) THEN
            IF(NIBLO.GT.1) THEN
              CALL PROPAGS(FL1, FL3, MIJS, MIJL, IG, SL, FCONST)
            ELSE
              FL3 = FL1
            ENDIF
          ENDIF
          CALL IMPLSCH (FL3, FL, MIJS, MIJL, IG,
     &     THWOLD,USOLD,TAUW,Z0OLD,    ROAIRO,ZIDLOLD,
     &     U10NEW,THWNEW,USNEW,Z0NEW,  ROAIRN,ZIDLNEW,
     &     SL,FCONST)
          GO TO 1000
        ENDIF

      ENDIF

! ----------------------------------------------------------------------

      RETURN
      END SUBROUTINE PHYS
