      SUBROUTINE STRESS(IU06,ITEST)

! ----------------------------------------------------------------------

!**** *STRESS* - COMPUTATION OF TOTAL STRESS.

!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990
!     J. BIDLOT           ECMWF     SEPTEMBER 1996 : REMOVE Z0 DUE TO 
!                                   VISCOSITY AND ADD ZERO STRESS FOR
!                                   ZERO WIND.     
!     BJORN HANSEN        ECMWF     MAY 1997
!                                   STRESS FOR MORE THAN ONE LEVEL.
!     J. BIDLOT           ECMWF     OCTOBER 2004: USE QUADRATIC STEPPING
!                                                 FOR TAUW. COMPUTE THE
!                                                 SQRT OF TAUT HERE.

!*    PURPOSE.
!     ---------

!       TO GENERATE STRESS TABLE TAU(TAUW,U10).

!**   INTERFACE.
!     ----------

!       *CALL* *STRESS(IU06,ITEST)*
!          *IU06*  -  LOGICAL UNIT FOR PRINTER OUTPUT UNIT.
!          *ITEST* -  OUTPUT FLAG IF LE 0 NO EXTRA OUTPUT IS GENERATED

!     METHOD.
!     -------

!       A STEADY STATE WIND PROFILE IS ASSUMED.
!       THE WIND STRESS IS COMPUTED USING THE ROUGHNESSLENGTH

!                  Z1=Z0/SQRT(1-TAUW/TAU)

!       WHERE Z0 IS THE CHARNOCK RELATION , TAUW IS THE WAVE-
!       INDUCED STRESS AND TAU IS THE TOTAL STRESS.
!       WE SEARCH FOR STEADY-STATE SOLUTIONS FOR WHICH TAUW/TAU < 1.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : JPLEVC   ,ALPHA    ,XKAPPA   ,XNLEV
      USE YOWPCONS , ONLY : G
      USE YOWTABL  , ONLY : ITAUMAX  ,JUMAX    ,IUSTAR   ,IALPHA   ,
     &            USTARM   ,  JPLEVT   ,EPS1   ,UMAX     ,TAUT     ,
     &            DELTAUW  ,DELU

! ----------------------------------------------------------------------

      REAL, PARAMETER :: XM=0.50
      INTEGER, PARAMETER :: NITER=10

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *XM*        REAL      POWER OF TAUW/TAU IN ROUGHNESS LENGTH.
!      *NITER*     INTEGER   NUMBER OF ITERATIONS TO OBTAIN TOTAL STRESS

! ----------------------------------------------------------------------

!     0. ALLOCATE ARRAYS
!        ----------------

      IF(.NOT.ALLOCATED(TAUT)) ALLOCATE(TAUT(0:ITAUMAX,0:JUMAX,JPLEVT))

!*    1.DETERMINE TOTAL STRESS.
!       -----------------------

!*    1.1 INITIALISE CONSTANTS.
!         ---------------------

      TAUWMAX = USTARM 
      DELU    = UMAX/REAL(JUMAX)
      DELTAUW = TAUWMAX/REAL(ITAUMAX)

!*    1.2 DETERMINE STRESS.
!         -----------------

      DO JL=1,MIN(JPLEVT,JPLEVC)

        XL=XNLEV(JL)
        IF(ITEST.GE.1) THEN
          WRITE(IU06,*)' '
          WRITE(IU06,*)' STRESS FOR LEVEL HEIGHT ',XL,' m'
        ENDIF

        CDRAG = 0.0012875
        WCD = SQRT(CDRAG) 

        DO I=0,ITAUMAX

          ZTAUW   = (REAL(I)*DELTAUW)**2

          DO J=0,JUMAX
            UTOP    = REAL(J)*DELU
            USTOLD  = UTOP*WCD
            TAUOLD  = MAX(USTOLD**2, ZTAUW+EPS1)
            DO ITER=1,NITER
              X      = ZTAUW/TAUOLD
              UST    = SQRT(TAUOLD)
              Z0     = ALPHA*TAUOLD/(G)/(1.-X)**XM
              F      = UST-XKAPPA*UTOP/(LOG(XL/Z0))
              DELF   = 1.-XKAPPA*UTOP/(LOG(XL/Z0))**2*2./UST*
     &         (1.-(XM+1)*X)/(1.-X)
              UST    = UST-F/DELF
              TAUOLD =  MAX(UST**2., ZTAUW+EPS1)
            ENDDO
            TAUT(I,J,JL)  = SQRT(TAUOLD)
          ENDDO

        ENDDO

      ENDDO

!*    FORCE ZERO WIND TO HAVE ZERO STRESS

      DO JL=1,JPLEVT
        DO I=0,ITAUMAX
          TAUT(I,0,JL)=0.0
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE STRESS
