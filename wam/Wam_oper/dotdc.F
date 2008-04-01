      SUBROUTINE DOTDC (IJS, IJL, ISHALLO, FCONST)

! ----------------------------------------------------------------------

!**** *DOTDC* - SIGMA DOT FOR CURRENT REFRACTION.

!     H. GUNTHER   GKSS/ECMWF   17/02/91
!     J. BIDLOT    ECMWF FEBRUARY 1997 : MESSAGE PASSING
!     J. BIDLOT    ECMWF JANUARY 2004 : REMOVE READING FROM FILE.
!                                       THE SIGMA DOT ARE NOW IN SDOT.

!*    PURPOSE.
!     --------

!       SCATTER SIGMA/ SINH (2*K*D) TABLE.

!**   INTERFACE.
!     ----------

!       *CALL* *DOTDC (IJS, IJL, ISHALLO, FCONST)*
!          *IJS*     - INDEX OF FIRST GRID POINT.
!          *IJL*     - INDEX OF LAST  GRID POINT.
!          *ISHALLO* - SHALLOW WATER OPTION.
!          *FCONST * - SIGMA/SINH(2KD)

!     METHOD.
!     -------

!       SEQUENTIAL READ FROM UNIT.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NFRE
      USE YOWSHAL  , ONLY : TSIHKD   ,INDEP

! ----------------------------------------------------------------------

      REAL,DIMENSION(NINF:NSUP,NFRE):: FCONST

! ----------------------------------------------------------------------

!*    2. IF DEEP WATER RETURN.
!        ---------------------

      IF (ISHALLO.EQ.1) RETURN

!*    3. GATHER SIGMA /SINH(2KD) FROM TABLE.
!        -----------------------------------

      DO M=1,NFRE
        DO IJ=IJS,IJL
          FCONST(IJ,M) = TSIHKD(INDEP(IJ),M)
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE DOTDC
